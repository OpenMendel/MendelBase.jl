################################################################################
# This set of functions implements the Elston-Stewart algorithm for
# evaluating the likelihood of a pedigree.
################################################################################
#
# Required OpenMendel packages and modules.
#
# using Search
# using SearchSetup
# using DataStructures
# using ReadData
# using ElstonStewartPreparation
# using ModelConstruction

export elston_stewart_loglikelihood

"""
Evaluate the loglikelihood of a collection of independent pedigrees
by the Elston-Stewart algorithm.
"""
function elston_stewart_loglikelihood(penetrance::Function, prior::Function, 
  transmission::Function, pedigree::Pedigree, person::Person,
  locus::Locus, parameter::Parameter, instruction::Instruction,
  keyword::Dict{AbstractString, Any})

  loglikelihood = 0.0
  fill!(pedigree.loglikelihood[:, 2], 0.0)
  par = parameter.par
  #
  # Given independence, loop over all pedigrees.
  #
  for ped = 1:pedigree.pedigrees
    #
    # Retrieve the population frequencies of the lumped alleles.
    #
    if length(locus.lumped_frequency) > 0
      for l = 1:locus.model_loci
        loc = locus.model_locus[l]
        locus.frequency[loc][:, end] = locus.lumped_frequency[ped, :, l]
      end
    end
    #
    # Compute the likelihood of the current pedigree.
    # Throw an error when the likelihood is 0.
    #
    (inconsistent, pedigree.loglikelihood[ped, 2]) =
      compute_likelihood(penetrance, prior, transmission, person, 
      locus, instruction, par, keyword, ped)
    if inconsistent
      name = pedigree.name[ped]
      throw(ArgumentError(
        "Error: in likelihood evaluation for pedigree $name.\n \n"))
    end
    loglikelihood = loglikelihood + pedigree.loglikelihood[ped, 2]
  end
  return loglikelihood
end # function elston_stewart_loglikelihood

"""
Compute the likelihood of a single pedigree
by the Elston-Stewart algorithm.
Various operations are performed on an array of arrays.
"""
function compute_likelihood(penetrance::Function, prior::Function, 
  transmission::Function, person::Person, locus::Locus,
  instruction::Instruction, par::Vector{Float64},
  keyword::Dict{AbstractString, Any}, ped::Int)
  #
  # Initialize the loglikelihood, the number of arrays, and the
  # array of arrays.
  #
  inconsistent = false
  loglikelihood = 0.0
  arrays = instruction.finish[ped] - instruction.start[ped]
  array = Array{Array{Float64, 1}}(undef, arrays)
  #
  # Now proceed instruction by instruction.
  #
  extent = 0
  increment1 = 0
  increment2 = 0
  for n = instruction.start[ped]:instruction.finish[ped]
    #
    # Fetch the site, rank, and size of the new array created by the
    # current instruction. Then allocate the array. If a multiply
    # and add operation is indicated, then rank pertains to the
    # product array before the pivot index is summed over.
    #
    operation = instruction.operation[n]
    if operation != quit_processing_pedigree
      k = instruction.site[n]
      if operation != pure_add
        rank = instruction.rank[n]
        space = instruction.space[n]
        array[k] = zeros(space)
      end
    end
    #
    # Each instruction involves one of six operations. Two of these,
    # "multiply_and_add" and "pure_multiply, " are treated simultaneously.
    #
    if operation == penetrance_and_prior_array
      #
      # Create a penetrance and prior array.
      #
      array[k] = construct_penetrance_prior!(penetrance, prior, person, 
        locus, instruction, par, keyword, array[k], n)
      (negative_entry, no_positive_entry) = check_array_entries(array[k])
      if negative_entry || no_positive_entry
        inconsistent = true
        println("Error: in penetrance or prior array.")
        break
      end
    #
    # Create a transmission array. 
    #
    elseif operation == transmission_array
      array[k] = construct_transmission!(transmission, person, locus, 
        instruction, par, keyword, array[k], n)
      (negative_entry, no_positive_entry) = check_array_entries(array[k])
      if negative_entry || no_positive_entry
        inconsistent = true
        println("Error: in transmission array.")
        break
      end
    #
    # Perform a pure add.
    #
    elseif operation == pure_add
      (loglikelihood, zero_likelihood) = pure_add_operation(loglikelihood, 
        array[k])
      if zero_likelihood
        inconsistent = true
        println("Error: in pure add operation.")
        break
      end
    #
    # Perform a multiply and add or a pure multiply.
    #
    elseif (operation == multiply_and_add || operation == pure_multiply)
      #
      # Fetch the positions of the multiplicand arrays. The number
      # pivot_range is 0 for a pure multiply and one less than the number
      # of entries to be added for a multiply and add. For a multiply
      # and add, record the summation strides for the muliplicand arrays.
      #
      i = instruction.extra[n][1]
      j = instruction.extra[n][2]
      pivot_range = instruction.extra[n][3]
      stride1 = instruction.extra[n][4]
      stride2 = instruction.extra[n][5]
      #
      # Insert into more easily accessed arrays, the extent of each dimension
      # of the product array and the increments to the positions within each
      # the multiplicand arrays.
      #
      if rank > 0
        extent = instruction.extra[n][6:5 + rank]
        increment1 = instruction.extra[n][6 + rank:5 + 2 * rank]
        increment2 = instruction.extra[n][6 + 2 * rank:5 + 3 * rank]
      end
      #
      # Perform the operation.
      #
      (loglikelihood, zero_likelihood, array[k]) = 
        multiply_add_operation!(array[i], array[j], array[k], extent, increment1,
        increment2, rank, pivot_range, stride1, stride2, loglikelihood)
      if zero_likelihood
        inconsistent = true
        println("Error: in a multiply and add multiply operation.")
        break
      end
    elseif operation == quit_processing_pedigree
    #
    # Terminate processing of the current pedigree.
    #
      break
    end
  end
  return (inconsistent, loglikelihood)
end # function compute_likelihood

"""
Construct a penetrance array or the product of a
penetrance and prior array.
"""
function construct_penetrance_prior!(penetrance::Function, prior::Function,
  person::Person, locus::Locus, instruction::Instruction, par::Vector{Float64},
  keyword::Dict{AbstractString, Any}, array::Vector{Float64}, n::Int)

  start = instruction.extra[n][1]
  finish = instruction.extra[n][2]
  i = instruction.extra[n][3]
  #
  # Construct i's multiple locus genotypes.
  #
  genotypes = genotype_count(person, locus, i, start, finish)
  multi_genotype = construct_multigenotypes(person, locus, start, finish,
                                            genotypes, i)
  #
  # Compute the penetrance of each genotype of i.
  #
  for j = 1:genotypes
    a = penetrance(person, locus, multi_genotype[:, :, j],
                   par, keyword, start, finish, i)
    if a < 0.0
      throw(ArgumentError(
        "Negative penetrance computed for $person.name[i].\n \n"))
    end
    #
    # Include the penetrance of each identical twin of i.
    #
    m = person.next_twin[i]
    while m > 0
      b = penetrance(person, locus, multi_genotype[:, :, j],
                     par, keyword, start, finish, m)
      if b < 0.0
        throw(ArgumentError(
          "Negative penetrance computed for $person.name[m].\n \n"))
      end
      a = a * b
      m = next_twin[m]
    end
    #
    # Include the prior when i is a founder.
    #
    if person.mother[i] == 0
      b = prior(person, locus, multi_genotype[:, :, j], par, keyword,
                start, finish, i)
      if b < 0.0
        throw(ArgumentError(
          "Negative prior computed for $person.name[i].\n \n"))
      end
      a = a * b
    end
    array[j] = a
  end
  return array
end # function construct_penetrance_prior!

"""
Construct a transmission array from parent i to child j.
"""
function construct_transmission!(transmission::Function, person::Person, 
  locus::Locus, instruction::Instruction, par::Vector{Float64},
  keyword::Dict{AbstractString, Any}, array::Vector{Float64}, n::Int)
  #
  # Fetch the starting & finishing loci and i & j.
  #
  start = instruction.extra[n][1]
  finish = instruction.extra[n][2]
  i = instruction.extra[n][3]
  j = instruction.extra[n][4]
  #
  # Construct the parent's multiple locus genotypes.
  #
  i_genotypes = genotype_count(person, locus, i, start, finish)
  multi_genotype = construct_multigenotypes(person, locus, start, finish,
                                            i_genotypes, i)
  #
  # Construct the child's multiple locus gametes.
  #
  maternal = !person.male[i]
  j_genotypes = genotype_count(person, locus, j, start, finish)
  gamete = construct_gametes(person, locus, start, finish, j_genotypes, j,
                             maternal)
  #
  # Insert the transmission probabilities into the transmission array in the
  # correct reverse dictionary order.
  #
  m = 0
  for l = 1:j_genotypes
    for k = 1:i_genotypes
      m = m + 1
      array[m] = transmission(person::Person, locus::Locus, gamete[:, l],
        multi_genotype[:, :, k], par, keyword, start, finish, i, j)
      if array[m] < 0.0
        throw(ArgumentError(
          "Negative transmission computed "*
          "from parent $person.name[i] to child $person.name[j].\n \n"))
      end
    end
  end
  return array
end # function construct_transmission!

"""
Count i's possible multilocus genotypes between loci start and finish.
"""
function genotype_count(person::Person, locus::Locus, i::Int, start::Int,
  finish::Int)

  count = 1
  for l = start:finish
    loc = locus.model_locus[l]
    count = count * length(person.genotype[i, loc])
  end
  return count
end # function genotype_count

"""
Check whether the array A has any negative entries and
any positive entries.
"""
function check_array_entries(array_a::Vector{Float64})

  negative_entry = minimum(array_a) < 0.0
  no_positive_entry = maximum(array_a) <= 0.0
  return (negative_entry, no_positive_entry)
end # function check_array_entries

"""
Carry out a pure add operation.
"""
function pure_add_operation(loglikelihood::Float64, array::Vector{Float64})

  array_sum = sum(array)
  if array_sum <= 0.0
    zero_likelihood = true
  else
    zero_likelihood = false
    loglikelihood = loglikelihood + log(array_sum)
  end
  return (loglikelihood, zero_likelihood)
end # function pure_add_operation

"""
Multiply array1 times array2 and sum on the pivot index if required.
Inserts the result in array3 and rescales to avoid underflows and overflows.
?? Can this be recoded to call a LAPACK routine??
"""
function multiply_add_operation!(array1::Vector{Float64},
  array2::Vector{Float64}, array3::Vector{Float64},
  extent::Vector{Int}, increment1::Vector{Int}, increment2::Vector{Int},
  rank::Int, pivot_range::Int, stride1::Int, stride2::Int,
  loglikelihood::Float64)

  address = ones(Int, rank)
  #
  # Initialize the positions in the multiplicand arrays and the reverse
  # dictionary address in the product array.
  #
  position1 = 1
  position2 = 1
  big = 0.0
  #
  # Compute the product array entry by entry.
  #
  for i = 1:length(array3)
    a = array1[position1] * array2[position2]
    #
    # Sum on the pivot index whenever pivot_range > 0.
    #
    for j = 1:pivot_range
      position1 = position1 + stride1
      position2 = position2 + stride2
      a = a + array1[position1] * array2[position2]
    end
    #
    # Save the result and keep a running maximum.
    #
    array3[i] = a
    big = max(big, a)
    #
    # Update the address in the product array and the positions in the
    # multiplicand arrays.
    #
    for j = 1:rank
      if address[j] < extent[j]
        address[j] = address[j] + 1
        position1 = position1 + increment1[j]
        position2 = position2 + increment2[j]
        break
      else
        address[j] = 1
      end
    end
  end
  #
  # Check for a zero-likelihood error. If none exists, rescale and
  # adjust the loglikelihood.
  #
  if big <= 0.0
    zero_likelihood = true
  else
    zero_likelihood = false
    loglikelihood = loglikelihood + log(big)
    array3 = array3 / big
  end
  return (loglikelihood, zero_likelihood, array3)
end # function multiply_add_operation!

"""
Create the multi-locus genotypes for person i from the single locus genotypes.
Only model loci between positions start and finish are considered.
"""
function construct_multigenotypes(person::Person, locus::Locus,
  start::Int, finish::Int, genotypes::Int, i::Int)

  multi_genotype = zeros(Int, 2, locus.model_loci, genotypes)
  #
  # Create the multiple locus genotypes by taking the appropriate
  # Cartesian product. Nested for loops here by computing the reverse
  # dictionary order location of each element of the product.
  #
  copies = 1
  increment = 1
  for l = start:finish
    loc = locus.model_locus[l]
    increment = increment * length(person.genotype[i, loc])
    initial_position = 1
    for (gm, gp) in person.genotype[i, loc]
      current = initial_position
      while current <= genotypes
        j = current
        for k = 1:copies
          multi_genotype[1, l, j] = gm
          multi_genotype[2, l, j] = gp
          j = j + 1
        end
        current = current + increment
      end
      initial_position = initial_position + copies
    end
    copies = increment
  end
  return multi_genotype
end # function construct_multigenotypes

"""
Create the maternal or paternal gametes for person i
from the single locus genotypes.
Only model loci between positions start and finish are considered.
"""
function construct_gametes(person::Person, locus::Locus,
  start::Int, finish::Int, genotypes::Int, i::Int, maternal::Bool)

  gamete = zeros(Int, locus.model_loci, genotypes)
  #
  # Create the multiple locus gametes by taking the appropriate Cartesian
  # product. Nested for loops are avoided here by computing the reverse
  # dictionary order location of each element of the product.
  #
  copies = 1
  increment = 1
  for l = start:finish
    loc = locus.model_locus[l]
    increment = increment * length(person.genotype[i, loc])
    initial_position = 1
    for (gm, gp) in person.genotype[i, loc]
      current = initial_position
      while current <= genotypes
        j = current
        for k = 1:copies
          if maternal
            gamete[l, j] = gm
          else
            gamete[l, j] = gp
          end
          j = j + 1
        end
        current = current + increment
      end
      initial_position = initial_position + copies
    end
    copies = increment
  end
  return gamete
end # function construct_gametes

