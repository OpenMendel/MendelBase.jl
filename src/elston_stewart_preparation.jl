################################################################################
# This set of functions prepares OpenMendel to evaluate the likelihood
# of a collection of pedigrees by the Elston-Stewart algorithm.
################################################################################
#
# Required OpenMendel packages and modules.
#
# using DataStructures
# using GeneralUtilities
# using GeneticUtilities
# using ReadData

mutable struct Instruction
  instructions :: Int
  start :: Vector{Int}
  finish :: Vector{Int}
  operation :: Vector{Int}
  site :: Vector{Int}
  rank :: Vector{Int}
  space :: Vector{Int}
  extra :: Array{Vector{Int}, 1}
end # Instruction

export Instruction
export allele_consolidation, genotype_elimination!, orchestrate_likelihood

penetrance_and_prior_array = 1
transmission_array         = 2
pure_add                   = 3
multiply_and_add           = 4
pure_multiply              = 5
quit_processing_pedigree   = 6

export penetrance_and_prior_array, transmission_array, pure_add
export multiply_and_add, pure_multiply, quit_processing_pedigree
export model_recombination_fractions

"""
Orchestrate allele lumping, genotype consolidation,
and creation of the Elston-Stewart instructions.
"""
function orchestrate_likelihood(pedigree::Pedigree, person::Person,
  nuclear_family::NuclearFamily, locus::Locus,
  keyword::Dict{AbstractString, Any})

  pedigrees = pedigree.pedigrees
  populations = person.populations
  model_loci = locus.model_loci
  eliminate_genotypes = keyword["eliminate_genotypes"]
  lump_alleles = keyword["lump_alleles"]
  product_mode = keyword["product_mode"]
  complexity_threshold = keyword["complexity_threshold"]
  #
  # Define the frequency of the lumped allele at each model locus.
  #
  if lump_alleles
    lumped_frequency = zeros(pedigrees, populations, model_loci)
    for l = 1:model_loci
      loc = locus.model_locus[l]
      lumped_frequency[:, :, l] =
        allele_consolidation(pedigree, person, locus, loc)
    end
  else
    lumped_frequency = zeros(0, 0, 0)
  end
  #
  # Perform genotype elimination over the model loci.
  #
  if eliminate_genotypes
    for l = 1:model_loci
      loc = locus.model_locus[l]
      genotype_elimination!(pedigree, person, nuclear_family, locus, loc)
    end
  end
  #
  # Find an upper bound on the number of instructions.
  #
  bound = person.people
  if product_mode
    bound = bound * model_loci
  end
  bound = 6 * bound + pedigrees
  #
  # Allocate the components of the instruction vector and the first
  # and last instruction associated with each pedigree.
  #
  elston_stewart_count = zeros(pedigrees)
  instruction_start = zeros(Int, pedigrees)
  instruction_finish = zeros(Int, pedigrees)
  operation = zeros(Int, bound)
  site = zeros(Int, bound)
  rank = zeros(Int, bound)
  space = zeros(Int, bound)
  extra = Array{Array{Int, 1}}(undef, bound)
  for i = 1:bound
    extra[i] = zeros(Int, 0)
  end
  #
  # Find a good summation sequence.
  #
  instructions = prepare_elston_stewart(pedigree, person, locus,
    elston_stewart_count, instruction_start, instruction_finish, operation,
    site, rank, space, extra, bound, product_mode, complexity_threshold)
  #
  # Return the instructions for the Elston-Stewart algorithm.
  #
  instruction = Instruction(instructions, instruction_start, 
    instruction_finish, operation, site, rank, space, extra)
  locus.lumped_frequency = lumped_frequency
  return (instruction, elston_stewart_count)
end  # function orchestrate_likelihood

"""
Consolidate alleles in each pedigree at locus number loc.
The lumped allele is numbered alleles + 1,
where alleles counts the number of alleles at the locus.
The frequency of the lumped allele and the lumping status
of each pedigree is returned along with the altered
genotype set for each person.
Warning: Allele consolidation should precede genotype elimination.
"""
function allele_consolidation(pedigree::Pedigree, person::Person,
  locus::Locus, loc::Int)
  #
  # Initialize variables.
  #
  xlinked = locus.xlinked[loc]
  pedigrees = pedigree.pedigrees
  (populations, alleles) = size(locus.frequency[loc])
  maximum_genotypes = alleles^2
  lumped_allele = alleles + 1
  lumped = falses(pedigrees)
  lumped_frequency = ones(pedigrees, populations)
  #
  # Find the alleles observed or partially observed in
  # the current pedigree.
  #
  for ped = 1:pedigrees
    observed = BitSet()
    for i = pedigree.start[ped]:pedigree.twin_finish[ped]
      if xlinked && person.male[i]
        if length(person.genotype[i, loc]) < alleles
          for (gm, gp) in person.genotype[i, loc]
            push!(observed, gm)
          end
        end
      else
        if length(person.genotype[i, loc]) < maximum_genotypes
          for (gm, gp) in person.genotype[i, loc]
            push!(observed, gm)
            push!(observed, gp)
          end
        end
      end
    end
    #
    # Find the frequency of the lumped allele.
    #
    if length(observed) < alleles - 1
      lumped[ped] = true
      for a in observed
        f = locus.frequency[loc][:, a]
        lumped_frequency[ped, :] = lumped_frequency[ped, :] - f
      end
    else
      lumped_frequency[ped, :] .= 0.0
    end
    #
    # Construct a new genotype set for each person, replacing
    # each unobserved allele by the lumped allele.
    #
    for i = pedigree.start[ped]:pedigree.twin_finish[ped]
      new_genotype_set = Set{Tuple{Int, Int}}()
      for (gm, gp) in person.genotype[i, loc]
        if in(gm, observed)
          hm = gm
        else
          hm = lumped_allele
        end
        if in(gp, observed)
          hp = gp
        else
          hp = lumped_allele
        end
        push!(new_genotype_set, (hm, hp))
      end
      person.genotype[i, loc] = deepcopy(new_genotype_set)
    end
  end
  #
  # Recount the number of homozygous genotypes per person and return
  # the lumped genotypes and lumped allele frequencies.
  #
  count_homozygotes!(person, loc)
  return lumped_frequency
end # function allele_consolidation

"""
Perform genotype elimination at locus number loc.
Genotypes for person i are stored in the set person.genotype[i, loc].
"""
function genotype_elimination!(pedigree::Pedigree, person::Person,
  nuclear_family::NuclearFamily, locus::Locus, loc::Int)

  people = person.people
  pedigrees = pedigree.pedigrees
  families = nuclear_family.families
  xlinked = locus.xlinked[loc]
  stable_pedigree = falses(pedigrees)
  stable_family = falses(families)
  saved_genotype = Array{Set{Tuple{Int, Int}}}(undef, people)
  #
  # Begin elimination for locus loc.
  #
  progress = true
  while progress
  #
  # Loop over all nuclear families.
  #
    for fam = 1:families
      stable_family[fam] = true
      ped = nuclear_family.pedigree[fam]
      if stable_pedigree[ped]; continue; end
      #
      # Initialize the saved sets of genotypes to be empty.
      #
      mom = nuclear_family.mother[fam]
      dad = nuclear_family.father[fam]
      saved_genotype[mom] = Set{Tuple{Int, Int}}()
      saved_genotype[dad] = Set{Tuple{Int, Int}}()
      for sib in nuclear_family.sib[fam]
        saved_genotype[sib] = Set{Tuple{Int, Int}}()
      end
      #
      # Loop over all combinations of maternal and paternal genotypes.
      #
      for (m1, m2) in person.genotype[mom, loc]
        for (d1, d2) in person.genotype[dad, loc]
          family_match = true
          #
          # Loop over all siblings.
          #
          for sib in nuclear_family.sib[fam]
            #
            # Loop over all genotypes of the current sib.
            #
            sib_match = false
            for (s1, s2) in person.genotype[sib, loc]
              #
              # Check for compatibility between the sib genotype and
              # the parental genotypes.
              #
              maternal_match = (s1 == m1 || s1 == m2)
              if !xlinked || !person.male[sib]
                paternal_match = (s2 == d1 || s2 == d2)
                sib_match = maternal_match && paternal_match
              else
                sib_match = maternal_match
              end
              #
              # As soon as a compatible sib genotype is found, stop looking.
              #
              if sib_match; break; end
            end  # sib genotype loop
            #
            # If no genotype of the current sib matches, exit.
            #
            family_match = family_match && sib_match
            if !family_match; break; end
          end  # sib loop
          #
          # If all sibs provide a match, save all matching parental and
          # sibling genotypes.
          #
          if !family_match; continue; end
          push!(saved_genotype[mom], (m1, m2))
          push!(saved_genotype[dad], (d1, d2))
          for sib in nuclear_family.sib[fam]
            for (s1, s2) in person.genotype[sib, loc]
              maternal_match = (s1 == m1 || s1 == m2)
              if !xlinked || !person.male[sib]
                paternal_match = (s2 == d1 || s2 == d2)
              else
                paternal_match = true
              end
              if maternal_match && paternal_match
                push!(saved_genotype[sib], (s1, s2))
              end
            end
          end  # sib loop
        end  # dad loop
      end  # mom loop
      #
      # Replace the current genotype sets in the family by the saved
      # genotype sets.
      #
      if length(saved_genotype[mom]) < length(person.genotype[mom, loc])
        stable_family[fam] = false
        person.genotype[mom, loc] = intersect(saved_genotype[mom],
          person.genotype[mom, loc])
      end
      if length(saved_genotype[dad]) < length(person.genotype[dad, loc])
        stable_family[fam] = false
        person.genotype[dad, loc] = intersect(saved_genotype[dad],
          person.genotype[dad, loc])
      end
      for sib in nuclear_family.sib[fam]
        if length(saved_genotype[sib]) < length(person.genotype[sib, loc])
          stable_family[fam] = false
          person.genotype[sib, loc] = intersect(saved_genotype[sib],
            person.genotype[sib, loc])
        end
      end  # sib loop
    end  # family loop
    #
    # Genotype elimination continues until genotypes stabilize
    # across all pedigrees.
    #
    fill!(stable_pedigree, true)
    for fam = 1:families
      ped = nuclear_family.pedigree[fam]
      stable_pedigree[ped] = stable_pedigree[ped] && stable_family[fam]
    end
    progress = !all(stable_pedigree)
  end # progress loop
  #
  # Recount the number of homozygous genotypes per person.
  #
  count_homozygotes!(person, loc)
end # function genotype_elimination!

"""
Specify the steps of the Elston-Stewart algorithm for each pedigree.
"""
function prepare_elston_stewart(pedigree::Pedigree, person::Person,
  locus::Locus, elston_stewart_count::Vector{Float64},
  instruction_start::Vector{Int}, instruction_finish::Vector{Int},
  operation::Vector{Int}, site::Vector{Int}, rank::Vector{Int},
  space::Vector{Int}, extra::Vector{Vector{Int}}, bound::Int,
  product_mode::Bool, complexity_threshold::Float64)

  pedigrees = pedigree.pedigrees
  model_loci = locus.model_loci
  free_recombination = locus.free_recombination
  instructions = 0
  #
  # Record the homozygote counts at the model loci.
  #
  homozygotes = zeros(Int, person.people, model_loci)
  for l = 1:model_loci
    homozygotes[:, l] = person.homozygotes[:, locus.model_locus[l]]
  end
  #
  # Create the instructions pedigree by pedigree.
  #
  for ped = 1:pedigrees
    ped_start = pedigree.start[ped]
    ped_finish = pedigree.finish[ped]
    ped_size = ped_finish - ped_start + 1
    pedigree_product_mode = product_mode
    #
    # Set up the adjacency graph of the current pedigree. From the adjacency
    # graph, determine a summation sequence by the greedy algorithm.
    # Note that greedy changes the set array graph.
    #
    (nodes, graph, weight) = adjacency_graph(pedigree, person, locus,
      homozygotes, pedigree_product_mode, ped, model_loci, free_recombination)
    (lowest_cost_sum, permutation) = greedy(nodes, graph, weight, ped_size)
    #
    # When likelihood evaluation can proceed in product mode, check if
    # non-product mode will do better.
    #
    if product_mode && !free_recombination
      product_mode_cost = lowest_cost_sum
      pedigree_product_mode = false
      (nodes, graph, weight) = adjacency_graph(pedigree, person, locus,
        homozygotes, pedigree_product_mode, ped, model_loci, free_recombination)
      (lowest_cost_sum, permutation) = greedy(nodes, graph, weight, ped_size)
      #
      # When product mode is superior, revert to it.
      #
      if lowest_cost_sum > product_mode_cost
        pedigree_product_mode = true
        (nodes, graph, weight) = adjacency_graph(pedigree, person, locus,
          homozygotes, pedigree_product_mode, ped, model_loci,
          free_recombination)
        (lowest_cost_sum, permutation) = greedy(nodes, graph, weight, ped_size)
      end
    end
    #
    # Decide whether the current pedigree is too demanding for likelihood
    # computation. If not, set up the index sets corresponding to the
    # various arrays encountered in likelihood evaluation.
    #
    if lowest_cost_sum <= complexity_threshold
      (arrays, index_set) = create_index_set(pedigree, person, locus,
        nodes, ped, pedigree_product_mode)
      #
      # Prepare the instructions for the Elston-Stewart algorithm.
      #
      instruction_start[ped] = instructions + 1
      elston_stewart_count[ped] = 0.0
      instructions = prepare_instructions(pedigree, person, locus,
        elston_stewart_count, index_set, weight, permutation, extra,
        operation, rank, space, site, arrays, instructions, nodes,
        bound, ped, pedigree_product_mode)
      instruction_finish[ped] = instructions
    #
    # If the current pedigree is too complex, set the Elston-Stewart count
    # equal to the lowest cost sum.
    #
    else
      instructions = instructions + 1
      instruction_start[ped] = instructions
      elston_stewart_count[ped] = lowest_cost_sum
      instruction_finish[ped] = instructions
    end
  end
  return instructions
end # function prepare_elston_stewart

"""
Define the adjacency graph connecting the indices of summation.
"""
function adjacency_graph(pedigree::Pedigree, person::Person, locus::Locus,
  homozygotes::Matrix{Int}, pedigree_product_mode::Bool, ped::Int,
  model_loci::Int, free_recombination::Bool)

  ped_start = pedigree.start[ped]
  ped_finish = pedigree.finish[ped]
  ped_size = ped_finish - ped_start + 1
  offset = ped_start - 1
  #
  # In product mode, nodes of the adjacency graph are person-locus pairs.
  # In non-product mode, nodes are simply people.
  #
  if pedigree_product_mode
    nodes = ped_size * model_loci
  else
    nodes = ped_size
  end
  #
  # Each node has a set of neighbors graph[node]. The weight of a node equals
  # the number of its associated genotypes or multi-locus genotypes.
  #
  weight = ones(Int, nodes)
  graph = empties(nodes)
  #
  # In product mode, loop over all person-locus pairs. Number nodes by
  # putting them in a matrix with loci on rows and people on columns.
  # A node is numbered by the its column-major position in the matrix.
  #
  if pedigree_product_mode
    for i = ped_start:ped_finish
      mom = person.mother[i]
      dad = person.father[i]
      for j = 1:model_loci
        jloc = locus.model_locus[j]
        node = (i - ped_start) * model_loci + j
        weight[node] = length(person.genotype[i, jloc])
        for k = 1:model_loci
          kloc = locus.model_locus[k]
          self_connected = false
          if mom != 0
            if !loci_split(homozygotes[mom, :], j, k, free_recombination)           
              if j != k; self_connected = true; end
              push!(graph[node], (mom - ped_start) * model_loci + k)
            end
          end
          if dad != 0
            if !loci_split(homozygotes[dad, :], j, k, free_recombination)
              if j != k; self_connected = true; end
              push!(graph[node], (dad - ped_start) * model_loci + k)
            end
          end
          if length(person.children[i]) > 0 &&
            !loci_split(homozygotes[i, :], j, k, free_recombination)
            if j != k; self_connected = true; end
            for kid in person.children[i]
              push!(graph[node], (kid - ped_start) * model_loci + k)
            end
          end
          if self_connected
            push!(graph[node], (i - ped_start) * model_loci + k)
          end
        end
      end
    end
  else
  #
  # In non-product mode, two people are neighbors if and only if one
  # is the parent of the other. If this is true, then they occur as
  # indices in the same transmission array. To insure that nodes are
  # numbered starting with 1, subtract an offset from each person number.
  #
    for node = 1:nodes
      i = node + offset
      for kid in person.children[i]
        push!(graph[node], kid - offset)
      end
      if person.mother[i] != 0
        push!(graph[node], person.mother[i] - offset)
        push!(graph[node], person.father[i] - offset)
      end
      for l = 1:model_loci
        loc = locus.model_locus[l]
        weight[node] = weight[node] * length(person.genotype[i, loc])
      end
    end
  end
  return (nodes, graph, weight)
end # function adjacency_graph

"""
Decide if two model loci are split within person i
in the sense of being separated by an obligate heterozygous locus.
Under free recombination, the likelihood factors over different loci.
"""
function loci_split(homozygotes, locus_1::Int, locus_2::Int,
                    free_recombination::Bool)

  if free_recombination
    split_loci = locus_1 != locus_2
  else
    split_loci = false
    for loc = min(locus_1, locus_2) + 1:max(locus_1, locus_2) - 1
      if homozygotes[loc] == 0
        return true
      end
    end
  end
  return split_loci
end # function loci_split

"""
Select a good summation sequence by the greedy method
and store it in a permutation. Observe that the cost
of removing a node from the adjacency graph is a simple surrogate
for the number of arithmetic operations involved in summing over
the genotypes possible for the node.
"""
function greedy(nodes::Int, graph::Vector{BitSet},
                weight::Vector{Int}, ped_size::Int)

  log_cost = zeros(nodes)
  log_weight = zeros(nodes)
  #
  # Initialize the permutation and the log weight attached to each node.
  #
  permutation = collect(1:nodes)
  for node = 1:nodes
    log_weight[node] = log(weight[node])
  end
  lowest_cost_sum = 0.0
  #
  # Compute the log cost of eliminating a node.
  #
  for node = 1:nodes
    log_sum = log_weight[node]
    for element in graph[node]
      log_sum = log_sum + log_weight[element]
    end
    log_cost[node] = log_sum
  end
  #
  # Find the lowest cost node to eliminate. Note that the variable last
  # limits the search to nodes representing the current locus when we are
  # in product mode. Nodes are therefore considered in clumps starting
  # with the leftmost locus and ending with the rightmost locus.
  #
  for i = 1:nodes
    lowest_cost = log_cost[permutation[i]]
    lowest_cost_index = i
    last = (div(i - 1, ped_size) + 1) * ped_size
    for j = i + 1:last
      competing_cost = log_cost[permutation[j]]
      if competing_cost < lowest_cost
        (lowest_cost, lowest_cost_index) = (competing_cost, j)
      end
    end
    low = exp(lowest_cost)
    lowest_cost_sum = lowest_cost_sum + min(low, 1e20)
    #
    # Insert the lowest cost node into the permutation.
    #
    (permutation[i], permutation[lowest_cost_index]) =
      (permutation[lowest_cost_index], permutation[i])
    best = permutation[i]
    #
    # Expand the neighborhood of each neighbor of the lowest cost node and
    # Delete the lowest cost node and neighbor from each revised neighborhood.
    #
    for neighbor in graph[best]
      graph[neighbor] = union(graph[neighbor], graph[best])
      delete!(graph[neighbor], best)
      delete!(graph[neighbor], neighbor)
      #
      # Recompute the log cost of removing each neighbor.
      #
      log_sum = log_weight[neighbor]
      for element in graph[neighbor]
        log_sum = log_sum + log_weight[element]
      end
      log_cost[neighbor] = log_sum
    end
  end
  return (lowest_cost_sum, permutation)
end  # function greedy

"""
Create the index sets for all initial arrays in likelihood computation.
"""
function create_index_set(pedigree::Pedigree, person::Person,
  locus::Locus, nodes::Int, ped::Int, pedigree_product_mode::Bool)
  #
  # Allocate the index sets and a variable n tracking the current index set.
  #
  model_loci = locus.model_loci
  free_recombination = locus.free_recombination
  ped_start = pedigree.start[ped]
  ped_finish = pedigree.finish[ped]
  offset = ped_start - 1
  index_set = empties(6 * nodes)
  n = 0
  #
  # In product mode, loop over all person-locus pairs. Prepare for the
  # penetrance-prior arrays first.
  #
  if pedigree_product_mode
    for i = ped_start:ped_finish
      for l = 1:model_loci
        n = n + 1
        push!(index_set[n], (i - ped_start) * model_loci + l)
      end
    end
    #
    # Next prepare for the transmission arrays by looping over all children.
    # Transmission arrays will be considered parent by parent.
    #
    for i = ped_start:ped_finish
      if person.mother[i] != 0
        for parent = 1:2
          if parent == 1
            j = person.mother[i]
          else
            j = person.father[i]
          end
          #
          # Find a first and last locus in a block of loci not split by
          # an obligate heterozygous locus in the parent.
          #
          (first, last) = (1, 0)
          while last < model_loci
          #
          # Under free recombination, each transmission arrays completely 
          # factors.
          #
            if free_recombination
              last = first
            else
              last = model_loci
              for l = first + 1:model_loci
                loc = locus.model_locus[l]
                if person.homozygotes[j, loc] == 0
                  last = l
                  break
                end
              end
            end
            #
            # Insert the parent-locus pairs into the index set before the 
            # child-locus pairs.
            #
            n = n + 1
            spread = last - first + 1
            for l = first:last
              push!(index_set[n], (j - ped_start) * model_loci + l)
              push!(index_set[n], (i - ped_start) * model_loci + l)
            end
            #
            # Exit if the last locus has been hit. Otherwise, reset the 
            # first locus to the current last locus and find another block
            # of loci.
            #
            if last == model_loci; break; end
            if free_recombination; last = first + 1; end
            first = last
          end
        end
      end
    end
  else
    #
    # In nonproduct mode, prepare for the penetrance-prior arrays first.
    #
    for i = ped_start:ped_finish
      n = n + 1
      push!(index_set[n], i - offset)
    end
    #
    # Next prepare for the transmission arrays by looping over all children.
    #
    for i = ped_start:ped_finish
      if person.mother[i] != 0
        n = n + 1
        push!(index_set[n], person.mother[i] - offset)
        push!(index_set[n], i - offset)
        n = n + 1
        push!(index_set[n], person.father[i] - offset)
        push!(index_set[n], i - offset)
      end
    end
  end
  #
  # Return the initial index sets.
  #
  return (n, index_set)
end # function create_index_set

"""
Prepare the instructions for performing the Elston-Stewart algorithm.
"""
function prepare_instructions(pedigree::Pedigree, person::Person, locus::Locus,
  elston_stewart_count::Vector{Float64}, index_set::Vector{BitSet}, 
  weight::Vector{Int}, permutation::Vector{Int}, extra::Vector{Vector{Int}},
  operation::Vector{Int}, rank::Vector{Int}, space::Vector{Int},
  site::Vector{Int}, arrays::Int, instructions::Int, nodes::Int, bound::Int,
  ped::Int, pedigree_product_mode::Bool)

  model_loci = locus.model_loci
  ped_start = pedigree.start[ped]
  offset = ped_start - 1
  #
  # Let n be the index of the current instruction. Prepare a set
  # to hold which arrays are involved in the current operation.
  #
  n = instructions
  active_array = trues(6 * nodes)
  involved = BitSet()
  union_set = BitSet()
  #
  # Give the instructions for creating the initial arrays.
  #
  for i = 1:arrays
    n = n + 1
    if i <= nodes
      operation[n] = penetrance_and_prior_array
      extra[n] = zeros(Int, 3)
      #
      # Load the starting and finishing loci and the person involved
      # in a penetrance and prior array.
      #
      if pedigree_product_mode
        j = first(index_set[i]) - 1
        extra[n][1] = mod(j, model_loci) + 1
        extra[n][2] = extra[n][1]
        extra[n][3] = div(j, model_loci) + ped_start
      else
        extra[n][1] = 1
        extra[n][2] = model_loci
        extra[n][3] = first(index_set[i]) + offset
      end
    else
      operation[n] = transmission_array
      extra[n] = zeros(Int, 4)
      #
      # Load the starting and finishing loci and the parent and child
      # involved in a transmission array.
      #
      if pedigree_product_mode
        j = first(index_set[i]) - 1
        k = last(index_set[i]) - 1
        extra[n][1] = mod(j, model_loci) + 1
        extra[n][2] = mod(k, model_loci) + 1
        extra[n][3] = div(j, model_loci) + ped_start
        extra[n][4] = div(k, model_loci) + ped_start
      else
        extra[n][1] = 1
        extra[n][2] = model_loci
        extra[n][3] = first(index_set[i]) + offset
        extra[n][4] = last(index_set[i]) + offset
      end
    end
    #
    # Load the site, the rank, and the size of the current array.
    #
    site[n] = i
    rank[n] = length(index_set[i])
    space[n] = product_weights(index_set[i], weight)
  end
  #
  # Give the instructions for creating and destroying arrays in computing the
  # likelihood of the current pedigree.
  #
  for node = 1:nodes
    i = permutation[node]
    #
    # Search through the index sets and find those containing the node i
    # to be eliminated.
    #
    for iteration = 1:100000
      involved = BitSet()
      for k = 1:arrays
        if active_array[k] && i in index_set[k]
          push!(involved, k)
        end
      end
      #
      # Exit the loop if no index containing i remains.
      #
      if isempty(involved); break; end
      #
      # Update the instruction index.
      #
      n = n + 1
      #
      # For a pure add, record the array to be summed.
      #
      if length(involved) == 1
        operation[n] = pure_add
        j = first(involved)
        site[n] = j
        rank[n] = length(index_set[j])
        space[n] = product_weights(index_set[j], weight)
        elston_stewart_count[ped] = elston_stewart_count[ped] + space[n]
        active_array[j] = false
        break
      #
      # For a pure multiply or a multiply and add, choose a pair of arrays 
      # to combine.
      #
      else
        best_1 = 1
        best_2 = 2
        pivot = i
        #
        # For a pure multiply, find the best pair of arrays to multiply by a
        # greedy algorithm.
        #
        if length(involved) > 2
          pivot = 0
          best_size = Inf
          for j = 1:length(involved)
            k = select_set_element(involved, j)
            for l = j + 1:length(involved)
              m = select_set_element(involved, l)
              union_set = union(index_set[k], index_set[m])
              criterion = product_weights(union_set, weight)
              if criterion < best_size
                best_1 = j
                best_2 = l
                best_size = criterion
              end
            end
          end
        end
        #
        # Form the union of the index sets corresponding to the only or best
        # pair of arrays.
        #
        j = select_set_element(involved, best_1)
        k = select_set_element(involved, best_2)
        union_set = union(index_set[j], index_set[k])
        #
        # Combine nodes and insert the instructions for a pure multiply or
        # a multiply and add.
        #
        (rank[n], extra[n]) = combine_nodes(index_set[j], index_set[k],
          union_set, weight, j, k, n, pivot)
        #
        # Delete node i in a multiply and add. Also record the current 
        # operation.
        #
        space[n] = product_weights(union_set, weight)
        if length(involved) == 2
          elston_stewart_count[ped] = elston_stewart_count[ped] + 2 * space[n]
          delete!(union_set, i)
          operation[n] = multiply_and_add
          space[n] = product_weights(union_set, weight)
        else
          elston_stewart_count[ped] = elston_stewart_count[ped] + space[n]
          operation[n] = pure_multiply
        end
        #
        # Activate the array created and determine its location and space 
        # requirement.
        #
        arrays = arrays + 1
        index_set[arrays] = deepcopy(union_set)
        active_array[arrays] = true
        site[n] = arrays
        #
        # Deactivate the two arrays destroyed.
        #
        active_array[j] = false
        active_array[k] = false
      end
    end
  end
  #
  # Give the instruction for finishing the current pedigree and update the
  # number of instructions
  #
  n = n + 1
  operation[n] = quit_processing_pedigree
  return n
end # function prepare_instructions

"""
Consolidate the nodes for the multiplicand and product arrays into blocks.
Generate the instructions for a pure multiply or a multiply and add.
a and b are the index sets for the multiplicand arrays,
and a_union_b is the index set for the product array.
n is the current instruction index, and pivot is the node to
be eliminated in a multiply and add. a_location and b_location give
the locations of the multiplicand arrays in the evolving sequence of arrays.
"""
function combine_nodes(a::BitSet, b::BitSet, a_union_b::BitSet,
  weight::Vector{Int}, a_location::Int, b_location::Int, nodes::Int, pivot::Int)
  #
  # Initialize current positions in the original index sets and
  # the block index sets. Also initialize the position of the pivot
  # and the type indicator for block construction.
  #
  extent = zeros(Int, length(a_union_b))
  a_set_position = 1
  b_set_position = 1
  a_block_position = 0
  b_block_position = 0
  union_block_position = 0
  pivot_block_position = 0
  previous_type = -2
  a_stride = 0
  b_stride = 0
  #
  # Create and initialize arrays.
  #
  c = length(a_union_b)
  a_block = zeros(Int, c)
  b_block = zeros(Int, c)
  a_increment = zeros(Int, c)
  b_increment = zeros(Int, c)
  #
  # March through the union of the two original index sets, defining
  # blocks as we go. If the a set has been exhausted, flag its current
  # node as impossible; likewise for the b set.
  #
  while a_set_position <= length(a) || b_set_position <= length(b)
    if a_set_position > length(a)
      current_a_node = nodes + 1
    else
      current_a_node = select_set_element(a, a_set_position)
    end
    if b_set_position > length(b)
      current_b_node = nodes + 1
    else
      current_b_node = select_set_element(b, b_set_position)
    end
    #
    # A new block involving indices from the a set must be started.
    #
    if current_a_node < current_b_node
      new_type = 1
      current_node = current_a_node
      a_set_position = a_set_position + 1
    #
    # A new block involving indices from the b set must be started.
    #
    elseif current_a_node > current_b_node
      new_type = -1
      current_node = current_b_node
      b_set_position = b_set_position + 1
    #
    # The current block continues.
    #
    else
      new_type = 0
      current_node = current_a_node
      a_set_position = a_set_position + 1
      b_set_position = b_set_position + 1
    end
    #
    # A block boundary is hit whenever the pivot node for a multiply and
    # add is encountered.
    #
    if current_node == pivot
      previous_type = -2
      pivot_block_position = union_block_position + 1
    end
#
# Update the product of the weights associated with a continuing block.
#
    if new_type == previous_type
      extent[union_block_position] = extent[union_block_position] * 
        weight[current_node]
    #
    # For a new block, record its position, initial weight, and its membership 
    # in the a and/or b block index sets.
    #
    else
      union_block_position = union_block_position + 1
      extent[union_block_position] = weight[current_node]
      if new_type == 1 || new_type == 0
        a_block_position = a_block_position + 1
        a_block[a_block_position] = union_block_position
      end
      if new_type == -1 || new_type == 0
        b_block_position = b_block_position + 1
        b_block[b_block_position] = union_block_position
      end
      #
      # Record whether the current block involves a indices, b indices, or 
      # both. The pivot marks the beginning of a new block.
      #
      previous_type = new_type
      if current_node == pivot; previous_type = -2; end
    end
  end
  #
  # Insert the rank of the product array into the current instruction.
  #
  rank = union_block_position
  a_blocks = a_block_position
  b_blocks = b_block_position
  #
  # Allocate the extra array for the current instruction. Record the
  # locations of the multiplicand arrays. Also load the pivot range, which
  # determines whether a pure multiply or a multiply and add is performed.
  #
  extra = zeros(Int, 5 + 3 * rank)
  extra[1] = a_location
  extra[2] = b_location
  if pivot_block_position == 0
    pivot_range = 0
  else
    pivot_range = extent[pivot_block_position] - 1
  end
  extra[3] = pivot_range
  #
  # Calculate the increments for the multiplicand arrays. Also calculate
  # their strides for a multiply and sum.
  #
  a_block_position = 1
  b_block_position = 1
  a_prod = 1
  b_prod = 1
  for union_block_position = 1:rank
    if union_block_position == pivot_block_position
      a_stride = a_prod
      b_stride = b_prod
    end
    if a_block_position <= a_blocks && a_block[a_block_position] == union_block_position
      a_prod = a_prod * extent[union_block_position]
      a_increment[union_block_position] = 1
      a_block_position = a_block_position + 1
    else
      a_increment[union_block_position] = 1 - a_prod
    end
    if b_block_position <= b_blocks && b_block[b_block_position] == union_block_position
      b_prod = b_prod * extent[union_block_position]
      b_increment[union_block_position] = 1
      b_block_position = b_block_position + 1
    else
      b_increment[union_block_position] = 1 - b_prod
    end
  end
  #
  # Load the strides into the current instruction.
  #
  extra[4] = a_stride
  extra[5] = b_stride
  #
  # Adjust the increment of each index pivot to correct for the add
  # operation in a multiply and add.
  #
  for i = 1:pivot_block_position - 1
    a_increment[i] = a_increment[i] - a_stride * pivot_range
    b_increment[i] = b_increment[i] - b_stride * pivot_range
  end
  #
  # Load the extents and increments into the current instruction.
  #
  true_rank = rank
  if pivot_block_position != 0; true_rank = rank - 1; end
  j = 5
  for i = 1:rank
    if i != pivot_block_position
      j = j + 1
      extra[j] = extent[i]
      extra[j + true_rank] = a_increment[i]
      extra[j + 2 * true_rank] = b_increment[i]
    end
  end
  rank = true_rank
  return (rank, extra)
end # function combine_nodes

"""
Find the product of the weights of the elements of the integer set S.
"""
function product_weights(s::BitSet, weight::Vector{Int})

  p = 1
  for element in s
    p = p * weight[element]
  end
  return p
end # function product_weights

"""
Compute recombination fractions between adjacent model loci.
"""
function model_recombination_fractions(locus::Locus,
  keyword::Dict{AbstractString, Any})

  n = locus.model_loci - 1
  theta = zeros(2, n)
  choice = keyword["genetic_map_function"]
  for l = 1:n
    j = locus.model_locus[l]
    k = locus.model_locus[l + 1]
    d = abs(locus.morgans[1, k] - locus.morgans[1, j])
    theta[1, l] = map_function(d, choice)
    d = abs(locus.morgans[2, k] - locus.morgans[2, j])
    theta[2, l] = map_function(d, choice)
  end
  return theta
end # function model_recombination_fractions

