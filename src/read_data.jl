################################################################################
# This set of functions reads the pedigree and locus data from external files
# and then constructs the appropriate data frames and data structures.
################################################################################
#
# Other OpenMendel modules.
#
using SnpArrays
# using DataStructures
# using GeneralUtilities
# using GeneticUtilities
#
# External modules.
#
using DataFrames    # From package DataFrames.

export count_homozygotes!, read_external_data_files

"""
This function organizes reading in the data from external files.
All names of data files are stored in the relevant keyword values.
"""
function read_external_data_files(keyword::Dict{ASCIIString, Any})
  #
  # Set the field separator used in the data files.
  #
  fsep = keyword["field_separator"]
  #
  # Read the data from the appropriate dataframes.
  # Start by reading the locus frame.
  #
  if keyword["locus_file"] == ""
    locus_frame = DataFrame() # An empty dataframe
  else
    locus_frame = readtable(keyword["locus_file"], separator = fsep)
  end
  #
  # Read the phenotype frame.
  #
  if keyword["phenotype_file"] == ""
    phenotype_frame = DataFrame() # An empty dataframe
  else
    phenotype_frame = readtable(keyword["phenotype_file"], separator = fsep)
  end
  #
  # Read the mandatory pedigree frame. Plink .fam files are allowed.
  # Other pedigree files should contain a header line. The default
  # field separator is a comma but can be changed to another character. 
  #
  if contains(keyword["pedigree_file"], ".fam")
    pedigree_frame = read_plink_fam_file(keyword["pedigree_file"], keyword)
  else
    pedigree_frame = readtable(keyword["pedigree_file"], separator = fsep)
  end
  #
  # Add a column recording the order of entry of each person.
  #
  pedigree_frame[:EntryOrder] = @data(collect(1:size(pedigree_frame, 1)))
  #
  # Check that ancestral populations are present in both the
  # pedigree and locus frames.
  #
  check_populations(locus_frame, pedigree_frame, keyword)
  #
  # Assemble the locus, pedigree, and person data structures.
  #
  locus = locus_information(locus_frame, pedigree_frame, keyword)
  pedigree = pedigree_information(pedigree_frame)
  person = person_information(locus_frame, pedigree_frame, phenotype_frame,
    locus, pedigree, keyword)
  #
  # Complete the pedigree data structure, and assemble the nuclear data 
  # family structure.
  #
  pedigree_counts!(pedigree, person)
  nuclear_family = construct_nuclear_families(pedigree, person)
  #
  # Check the data structures for inconsistencies.
  #
  check_data_structures!(pedigree, person, locus, nuclear_family, keyword)
  #
  # Assemble the SNP Definition frame and the SNP data structure.
  # First check that the two SNP files are both present or both absent.
  #
  if (keyword["snpdefinition_file"] == "") != (keyword["snpdata_file"] == "")
    throw(ArgumentError(
      "Either the SNP definition or SNP data file was not specified."))
  end
  if keyword["snpdefinition_file"] == ""
    #
    # Construct empty SNP data structures.
    #
    snp_definition_frame = DataFrame() # An empty dataframe.
    snpmatrix = SnpArray(0,0)
    snpdata = snp_information(snp_definition_frame, person, snpmatrix, keyword)
  else
    if contains(keyword["snpdefinition_file"], ".bim")
      snp_definition_frame = read_plink_bim_file(keyword["snpdefinition_file"],
        keyword)
    else
      snp_definition_frame = readtable(keyword["snpdefinition_file"],
        separator = fsep)
    end
    #
    # Read the SNP bed file.
    #
    snpmatrix = SnpArray(keyword["snpdata_file"], people = person.people,
      snps = size(snp_definition_frame, 1))
    #
    # Assemble the SNP data structure.
    #
    snpdata = snp_information(snp_definition_frame, person, snpmatrix, keyword)
    (snpdata.maf, snpdata.minor_allele, snpdata.missings_per_snp,
      snpdata.missings_per_person) = summarize(snpmatrix)
  end
  return (pedigree, person, nuclear_family, locus, snpdata,
    locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame)
end # function read_external_data_files

"""
This function converts a Plink .fam file into a dataframe and 
substitutes blanks for missing values.
"""
function read_plink_fam_file(plink_fam_file::AbstractString,
  keyword::Dict{ASCIIString, Any})

  column_types = [UTF8String, UTF8String, UTF8String, UTF8String,
    UTF8String, UTF8String]
  column_names = [:Pedigree, :Person, :Father, :Mother, :Sex, :PlinkTrait]
  fsep = keyword["field_separator"]
  fam_dframe = readtable(plink_fam_file, header = false, separator = fsep, 
    eltypes = column_types, names = column_names)
  for i = 1:size(fam_dframe, 1)
    if fam_dframe[i, :Father] == "0"
      fam_dframe[i, :Father] = " "
    end
    if fam_dframe[i, :Mother] == "0"
     fam_dframe[i, :Mother] = " "
    end
    if fam_dframe[i, :PlinkTrait] == "-9" || fam_dframe[i, :PlinkTrait] == "0"
      fam_dframe[i, :PlinkTrait] = " "
    end
  end
  return fam_dframe
end # function read_plink_fam_file

"""
This function converts a Plink .bim file into a dataframe.
"""
function read_plink_bim_file(plink_bim_file::AbstractString,
  keyword::Dict{ASCIIString, Any})

  column_types = [UTF8String, UTF8String, Float64, Int, 
    UTF8String, UTF8String]
  column_names = [:Chromosome, :SNP, :CentiMorgans, :Basepairs, 
    :Allele1, :Allele2]
  fsep = keyword["field_separator"]
  bim_dframe = readtable(plink_bim_file, header = false, separator = fsep, 
    eltypes = column_types, names = column_names)
  return bim_dframe
end # function read_plink_bim_file

"""
This function copies information from a SNP definition dataframe into
the SNP data structure.
"""
function snp_information(snp_definition_frame::DataFrame, person::Person,
  snpmatrix, keyword::Dict{ASCIIString, Any})

  if keyword["snpdefinition_file"] == ""
    snpdata = SnpData(0, 0, Vector{AbstractString}(), Vector{AbstractString}(),
      Vector{AbstractString}(), Vector{Float64}(), Vector{Int}(),
      Vector{AbstractString}(), Vector{AbstractString}(), Vector{Float64}(),
      BitVector(), snpmatrix, Vector{Int}(), Vector{Int}())
    return snpdata
  end
  (snps, columns) = size(snp_definition_frame)
  column_names = names(snp_definition_frame)
  snp_name = blanks(snps)
  if :Locus in column_names
    for snp = 1:snps
      snp_name[snp] = string(snp_definition_frame[snp, :Locus])
    end
  end
  if :SNP in column_names
    for snp = 1:snps
      snp_name[snp] = string(snp_definition_frame[snp, :SNP])
    end
  end
  chromosome = blanks(snps)
  if :Chromosome in column_names
    for snp = 1:snps
      chromosome[snp] = string(snp_definition_frame[snp, :Chromosome])
    end
  end
  if :CentiMorgans in column_names
    centimorgans = snp_definition_frame[:CentiMorgans]
    centimorgans = convert(Vector{Float64}, centimorgans)
  else
    centimorgans = zeros(Float64, snps)
  end
  if :Basepairs in column_names
    basepairs = snp_definition_frame[:BasePairs]
  else
    basepairs = zeros(Int, snps)
  end
  allele1 = blanks(snps)
  if :Allele1 in column_names
    for snp = 1:snps
      allele1[snp] = string(snp_definition_frame[snp, :Allele1])
    end
  else
    fill!(allele1, "1")
  end
  allele2 = blanks(snps)
  if :Allele2 in column_names
    for snp = 1:snps
      allele2[snp] = string(snp_definition_frame[snp, :Allele2])
    end
  else
    fill!(allele1, "2")
  end
  #   
  # Return the SNP data structure
  #   
  people = person.people
  personid = person.name
  snpdata = SnpData(people, snps, personid, snp_name, chromosome, centimorgans,
    basepairs, allele1, allele2, Vector{Float64}(), BitVector(), snpmatrix, 
    Vector{Int}(), Vector{Int}())
  return snpdata 
end # function snp_information

"""
This function extracts locus information from the locus frame.
First extract relevant dimensions and fields.
"""
function locus_information(locus_frame::DataFrame, pedigree_frame::DataFrame,
                           keyword::Dict{ASCIIString, Any})
  #
  # Create a null locus structure.
  #
  if length(locus_frame) == 0
    a = Array(Array{ASCIIString, 1}, 1); a[1] = blanks(1)
    b = Array(Array{Float64, 2}, 1); b[1] = zeros(1, 1)
    c = zeros(1, 1, 1)
    locus = Locus(0, 0, 0, true, blanks(0), blanks(0), zeros(Int, 0),
                  zeros(2, 0), zeros(2, 0), trues(0), zeros(Int, 0),
                  a, b, c, zeros(Int, 0), zeros(Int, 0))
    return locus
  end
  #
  # Otherwise, create a regular locus structure.
  #
  pedigree_field = names(pedigree_frame)
  rows = length(locus_frame[:, :Locus])
  locus_field = names(locus_frame)
  columns = length(locus_field)
  locus_name = unique(locus_frame[:, :Locus])
  #
  # Check for errors and omissions.
  #
  if :Morgans in locus_field && !(:FemaleMorgans in locus_field)
    rename!(locus_frame, :Morgans, :FemaleMorgans)
    locus_field = names(locus_frame)
  end
  loci = 1
  for i = 2:rows
    if locus_frame[i, :Locus] != locus_frame[i - 1, :Locus]; loci = loci + 1; end
  end
  if loci != length(locus_name)
    throw(ArgumentError(
      "The alleles of each locus must be contiguous.\n \n"))
  end
  for i = 2:rows
    if locus_frame[i, :Locus] == locus_frame[i - 1, :Locus]
      a = string(locus_frame[i, :Chromosome])
      if locus_frame[i, :Chromosome] != locus_frame[i - 1, :Chromosome]
        throw(ArgumentError(
          "The chromosome of locus $a is inconsistent.\n \n"))
      end
      if :BasesPairs in locus_field
        if locus_frame[i, :BasesPairs] != locus_frame[i - 1, :BasesPairs]
          throw(ArgumentError(
            "The bp position of locus $a is inconsistent.\n \n"))
        end
      end
      if :FemaleMorgans in locus_field
        if locus_frame[i, :FemaleMorgans] != locus_frame[i - 1, :FemaleMorgans]
          throw(ArgumentError(
            "The map position of locus $a is inconsistent.\n \n"))
        end
      end
      if :MaleMorgans in locus_field
        if locus_frame[i, :MaleMorgans] != locus_frame[i - 1, :MaleMorgans]
          throw(ArgumentError(
            "The map position of locus $a is inconsistent.\n \n"))
        end
      end
    end
  end
  #
  # Sort the locus frame by chromosome and location. Warning:
  # Chromosome names should be 1, 2, ... 22, or X without further
  # adornment.
  #
  if :BasePairs in locus_field
    sort!(locus_frame::DataFrame, cols = [:Chromosome, :BasePairs])
  elseif :FemaleMorgans in locus_field
    sort!(locus_frame::DataFrame, cols = [:Chromosome, :FemaleMorgans])
  end
  locus_name = unique(locus_frame[:, :Locus])
  #
  # Initialize the marker location variables.
  #
  base_pairs = zeros(Int, loci)
  morgans = zeros(2, loci)
  #
  # Classify each locus by chromosome and number of alleles.
  #
  alleles = zeros(Int, loci)
  chromosome = blanks(loci)
  xlinked = falses(loci)
  loc = 0
  for i = 1:rows
    if i == 1 || locus_frame[i, :Locus] != locus_frame[i - 1, :Locus]
      loc = loc + 1
      chromosome[loc] = string(locus_frame[i, :Chromosome])
      #
      # When only base pair distances are available, equate 1e6
      # base pairs to a centiMorgan.
      #
      if :BasesPairs in locus_field
        bases_pairs[loc] = locus_frame[i, :BasePairs]
        if !(:FemaleMorgans in locus_field)
          morgans[:, loc] = bases_pairs[loc] / 1e8
        end
      end
      if :FemaleMorgans in locus_field
        morgans[1, loc] = locus_frame[i, :FemaleMorgans]
        if :MaleMorgans in locus_field
          morgans[2, loc] = locus_frame[i, :MaleMorgans]
        else
          morgans[2, loc] = morgans[1, loc]
        end
      end
      c = chromosome[loc][1]
      xlinked[loc] = c == 'X' || c == 'x'
      alleles[loc] = 1
    else
      alleles[loc] = alleles[loc] + 1
    end
  end
  #
  # Identify population names and number.
  #
  population_names = keyword["populations"]
  populations = length(keyword["populations"])
  #
  # If there are no designated populations, then find
  # the field in the Locus frame that contains the allele frequencies.
  #
  if populations == 0
    for i = 1:columns
      if typeof(locus_frame[i]) == DataArray{Float64, 1}
        s = symbol(locus_field[i])
        if s != :FemaleMorgans && s != :MaleMorgans
          populations = 1
          push!(population_names, string(locus_field[i]))
          break
        end
      end
    end
  end
  #
  # Find the observed loci. These are the ones common to the Pedigree
  # frame and the Locus field in the Locus frame.
  #
  B = convert(Vector{Symbol}, (unique(locus_frame[:, :Locus])))
  C = intersect(pedigree_field, B)
  #
  # Find the indices of the matching field names in the Pedigree frame.
  #
  observed_indices = findin(pedigree_field, C)
  #
  # Find the indices of the observed loci in the locus structure.
  #
  observed_loci = length(C)
  observed_locus = zeros(Int, observed_loci)
  locus_field_in_pedigree_frame = zeros(Int, observed_loci)
  i = 0
  for loc = 1:loci
    locus_symbol = convert(Symbol, locus_name[loc])
    if locus_symbol in C
      i = i + 1
      observed_locus[i] = loc
      for j in observed_indices
        if isequal(pedigree_field[j], locus_symbol)
          locus_field_in_pedigree_frame[i] = j
          break
        end
      end
    end
  end
  #
  # Eliminate unobserved loci.
  #
  loci = observed_loci
  locus_name = getindex(locus_name, observed_locus)
  chromosome = getindex(chromosome, observed_locus)
  base_pairs = getindex(base_pairs, observed_locus)
  temp = copy(morgans)
  morgans = zeros(2, length(observed_locus))
  morgans[1, :] = getindex(temp[1, :], observed_locus)
  morgans[2, :] = getindex(temp[2, :], observed_locus)
  xlinked = getindex(xlinked, observed_locus)
  alleles = getindex(alleles, observed_locus)
  #
  # Collect the allele names and population frequencies for each
  # observed locus.
  #
  allele_name = Array(Array{ASCIIString, 1}, loci)
  frequency = Array(Array{Float64, 2}, loci)
  loc = 0
  n = 0
  for i = 1:rows
    if i == 1 || locus_frame[i, :Locus] != locus_frame[i - 1, :Locus]
      loc = loc + 1
      n = 1
      l = findfirst(observed_locus, loc)
      if l != 0
        allele_name[l] = Array(ASCIIString, alleles[l])
        if populations == 0
          frequency[l] = zeros(1, alleles[l] + 1)
        else
          frequency[l] = zeros(populations, alleles[l] + 1)
        end
      end
    else
      n = n + 1
    end
    #
    # Fill the allele frequency array with input frequencies.
    #
    l = findfirst(observed_locus, loc)
    if l != 0
      allele_name[l][n] = string(locus_frame[i, :Allele])
      j = 0
      for pop in population_names
        j = j + 1
        if isna(locus_frame[i, symbol(pop)])
          frequency[l][j, n] = 0.0
        else
          frequency[l][j, n] = locus_frame[i, symbol(pop)]
        end
      end
    end
  end
  #
  # Check that allele frequencies are legal.
  #
  for loc = 1:loci
    a = locus_name[loc]
    j = 0
    for pop in population_names
      j = j + 1
      if any(frequency[loc][j, :] .< 0.0)
        if populations == 1
          println("Warning: Negative allele frequency at locus $a.")
        else
          println(
          "Warning: Negative allele frequency at locus $a and population $pop.")
        end
        frequency[loc][j, :] = NaN
      end
      if any(isnan, frequency[loc][j, :])
        frequency[loc][j, :] = NaN
      else
        total = sum(frequency[loc][j, :])
        if abs(total - 1.0) > 1e-5
          if populations == 1
            println("Warning: Allele frequencies at locus",
              " $a do not sum to 1.0.")
          else
            println("Warning: Allele frequencies at locus",
              " $a for population $pop do not sum to 1.0.")
          end
          frequency[loc][j, :] = frequency[loc][j, :] / total
        end
      end
    end
  end
  #
  # Compute the recombination fraction between each pair
  # of adjacent observed loci. Locations are measured in Morgans.
  # In the absence of distances, set all recombination fractions
  # equal to 0.5.
  #
  if loci > 1
    if :BasesPairs in locus_field || :FemaleMorgans in locus_field
      theta = zeros(2, loci - 1)
      for loc = 1:loci - 1
        d = abs(morgans[1, loc + 1] - morgans[1, loc])
        theta[1, loc] = map_function(d, "Haldane")
        d = abs(morgans[2, loc + 1] - morgans[2, loc])
        theta[2, loc] = map_function(d, "Haldane")
      end
    else
      theta = 0.5 * ones(2, loci - 1)
    end
  else
    theta = zeros(2, 0)
  end
  #
  # Record if observed loci freely recombine.
  #
  if loci == 1 || all(theta .>= 0.5)
    free_recombination = true
  else
    free_recombination = false
  end
  #
  # Allocate space for the frequency of the lumped allele.
  #
  pedigrees = length(unique(pedigree_frame[:, :Pedigree]))
  lumped_frequency = zeros(pedigrees, populations, loci)
  #
  # As a default, equate the set of model loci to the set of
  # observed loci and the trait locus to 0.
  #
  model_loci = loci
  model_locus = collect(1:loci)
  trait = 0
  #
  # Insert the various scalars and arrays in the locus data structure.
  #
  locus = Locus(loci, model_loci, trait, free_recombination, locus_name,
    chromosome, base_pairs, morgans, theta, xlinked, alleles, allele_name, 
    frequency, lumped_frequency, model_locus, locus_field_in_pedigree_frame)
  return locus
end # function locus_information

"""
This function extracts pedigree information from the pedigree frame.
Some fields must be supplied later.
"""
function pedigree_information(pedigree_frame::DataFrame)
  #
  # Initialize arrays.
  #
  people = length(pedigree_frame[:, :Person])
  pedigrees = length(unique(pedigree_frame[:, :Pedigree]))
  pedigree_name = blanks(pedigrees)
  start = zeros(Int, pedigrees)
  twin_finish = zeros(Int, pedigrees)
  finish = zeros(Int, pedigrees)
  individuals = zeros(Int, pedigrees)
  founders = zeros(Int, pedigrees)
  females = zeros(Int, pedigrees)
  males = zeros(Int, pedigrees)
  twins = zeros(Int, pedigrees)
  families = zeros(Int, pedigrees)
  #
  # Check whether each pedigree occupies a contiguous block.
  # If not, throw an error.
  #
  ped = 1
  for i = 2:people
    if pedigree_frame[i, :Pedigree] != pedigree_frame[i - 1, :Pedigree]
      ped = ped + 1
    end
    if ped > pedigrees
      throw(ArgumentError("Some pedigrees do not occur in a contiguous block." *
      " Please fix the data files.\n"))
    end
  end
  #
  # Find the name, start, and finish of each pedigree.
  #
  ped = 0
  for i = 1:people
    if i == 1 || pedigree_frame[i, :Pedigree] != pedigree_frame[i - 1, :Pedigree]
      ped = ped + 1
      start[ped] = i
      if i > 1; twin_finish[ped - 1] = i - 1; end
      pedigree_name[ped] = string(pedigree_frame[i, :Pedigree])
    end
  end
  twin_finish[pedigrees] = people
  #
  # Insert the gathered information into the Pedigree structure.
  #
  pedigree = Pedigree(pedigrees, pedigree_name, start, twin_finish, finish,
    individuals, founders, females, males, twins, families, zeros(pedigrees, 2))
  return pedigree
end # function pedigree_information

"""
This function extracts person information from the pedigree frame.
"""
function person_information(locus_frame::DataFrame, pedigree_frame::DataFrame,
  phenotype_frame::DataFrame, locus::Locus, pedigree::Pedigree,
  keyword::Dict{ASCIIString, Any})
  #
  # Initialize arrays by their default values.
  #
  people = length(pedigree_frame[:, :Person])
  pedigrees = pedigree.pedigrees
  allele_separator = keyword["allele_separator"]
  ordered_allele_separator = keyword["ordered_allele_separator"]
  mother = zeros(Int, people)
  father = zeros(Int, people)
  male = falses(people)
  #
  # Check to see what fields are present.
  #
  pedigree_field = names(pedigree_frame)
  parents_present = false
  if symbol(:Mother) in pedigree_field
    mother_string = pedigree_frame[:, :Mother]
    parents_present = true
  end
  if symbol(:Father) in pedigree_field
    father_string = pedigree_frame[:, :Father]
  else
    parents_present = false
  end
  #
  # Record who is an identical twin.
  #
  next_twin = zeros(Int, people)
  primary_twin = zeros(Int, people)
  twins_present = symbol(:Twin) in pedigree_field
  if twins_present
  #
  # Create pointers from one twin in a twin set to the next
  # twin. The last twin points to no-one. In each twin set
  # designate a primary twin.
  #
    twin_group = unique(pedigree_frame[:, :Twin])
    for t in twin_group
      if typeof(t) == Int && t == 0
        continue
      elseif isna(t)
        continue
      else
        twin_list = findin(pedigree_frame[:, :Twin], t)
        primary_twin[twin_list[1]] = twin_list[1]
        for i = 2:length(twin_list)
          next_twin[twin_list[i - 1]] = twin_list[i]
          primary_twin[twin_list[i]] = twin_list[1]
        end
      end
    end
  end
  #
  # Identify mothers and fathers.
  #
  pedigree_number = zeros(Int, people)
  person_name = blanks(people)
  if parents_present
    for ped = 1:pedigrees
      for i = pedigree.start[ped]:pedigree.twin_finish[ped]
        pedigree_number[i] = ped
        person_name[i] = string(pedigree_frame[i, :Person])
        mother_found = false
        father_found = false
        for j = pedigree.start[ped]:pedigree.twin_finish[ped]
          if j == i; continue; end
          if !isna(mother_string[i]) &&
             mother_string[i] == pedigree_frame[j, :Person]
            mother[i] = j
            mother_found = true
          end
          if !isna(father_string[i]) &&
             father_string[i] == pedigree_frame[j, :Person]
            father[i] = j
            father_found = true
          end
          if mother_found && father_found; break; end
        end
      end
    end
  end
  #
  # Redefine the parents of children of co-twins.
  #
  if twins_present
    for i = 1:people
      (j, k) = (mother[i], father[i])
      if j != 0 && primary_twin[j] != 0
        mother[i] = primary_twin[j]
      end
      if k != 0 && primary_twin[k] != 0
        father[i] = primary_twin[k]
      end
    end
  end
  #
  # Find a permutation arranging parents before their children
  # and putting co-twins at the end of each pedigree.
  #
  if parents_present || twins_present
    (perm, per) = loop(pedigree, father, mother, primary_twin)
    if per != 0
      throw(ArgumentError(
        "Person $person_name[per] is his/her own ancestor.\n \n"))
    end
    permute!(person_name, perm)
    #
    # Find the inverse permutation.
    #
    inverse_perm = collect(1:people)
    for i = 1:people
      inverse_perm[perm[i]] = i
    end
    #
    # Redefine parental indicators.
    #
    parent = copy(mother)
    for i = 1:people
      if parent[i] != 0
        mother[inverse_perm[i]] = inverse_perm[parent[i]]
      else
        mother[inverse_perm[i]] = 0
      end
    end
    parent = copy(father)
    for i = 1:people
      if parent[i] != 0
        father[inverse_perm[i]] = inverse_perm[parent[i]]
      else
        father[inverse_perm[i]] = 0
      end
    end
    #
    # Redine twin indicators.
    #
    if twins_present
      p_twin = copy(primary_twin)
      n_twin = copy(next_twin) 
      for i = 1:people
        if p_twin[i] != 0
          primary_twin[inverse_perm[i]] = inverse_perm[p_twin[i]]
        else
          primary_twin[inverse_perm[i]] = 0
        end
        if n_twin[i] != 0
          next_twin[inverse_perm[i]] = inverse_perm[n_twin[i]]
        else
          next_twin[inverse_perm[i]] = 0
        end
      end
    end
    #
    # Sort the pedigree dataframe according to the inverse permutation.
    #
    pedigree_frame = hcat(pedigree_frame, inverse_perm)
    sort!(pedigree_frame, cols = size(pedigree_frame, 2))
    delete!(pedigree_frame, :x1)
  end
  #
  # Record who is male.
  #
  male_symbols = keyword["male"]
  female_symbols = keyword["female"]
  if symbol(:Sex) in pedigree_field
    for i = 1:people
      if isna(pedigree_frame[i, :Sex]); continue; end
      s = pedigree_frame[i, :Sex]
      if !isa(parse(string(s), raise=false), Number); s = lowercase(s); end
      if s in male_symbols
        male[i] = true
      elseif s in female_symbols
        male[i] = false
      else
        per = person_name[i]
        throw(ArgumentError(
          "Person $per has sex indicator $s, " *
          "which matches neither female nor male.\n \n"))
      end
    end
  end
  #
  # Check for a switched mother and father.
  #
  for i = 1:people
    j = mother[i]
    k = father[i]
    if j != 0 && k != 0
      if male[j] && !male[k]
        (mother[i], father[i]) = (father[i], mother[i])
      end
    end
  end
  #
  # Search for quantitative variables in the pedigree frame that are not
  # admixture proportions.
  #
  variables = 0
  for i = 1:length(pedigree_field)
    if typeof(pedigree_frame[i]) == DataArray{Float64, 1}
      if !(pedigree_field[i] in keyword["populations"])
        variables = variables + 1
      end
    end
  end
  #
  # Fill in the variables values and names.
  #
  variable = zeros(people, variables)
  variable_name = blanks(variables)
  variables = 0
  for i = 1:length(pedigree_field)
    if typeof(pedigree_frame[i]) == DataArray{Float64, 1}
      if !(pedigree_field[i] in keyword["populations"])
        variables = variables + 1
        variable_name[variables] = string(pedigree_field[i])
        for per = 1:people
          if isna(pedigree_frame[per, i])
            variable[per, variables] = NaN
          else
            variable[per, variables] = pedigree_frame[per, i]
          end
        end
      end
    end
  end
  #
  # Record the disease status of each person.
  #
  disease_field = keyword["disease_status"]
  if disease_field != ""
    disease_field = symbol(disease_field)
    if disease_field in pedigree_field
      disease_status = blanks(people)
      status = pedigree_frame[:, disease_field]
      for i = 1:people
        if isna(status[i])
          disease_status[i] = ""
        else
          disease_status[i] = string(status[i])
        end
      end
    else
      disease_field = keyword["disease_status"]
      throw(ArgumentError(
        "Specified disease status field ($disease_field) " *
        "is not in the pedigree frame.\n \n"))
    end
  else
    disease_status = blanks(0)
  end
  #
  # Find the locus names and number of loci in the phenotype dataframe.
  #
  if length(phenotype_frame) > 0
    locus_name = unique(phenotype_frame[:, :Locus])
    loci = length(locus_name)
    rows = length(phenotype_frame[:, :Locus])
    phenotypes = zeros(Int, loci)
    phenotype = Array(Array{ASCIIString, 1}, loci)
    genotype_string = Array(Array{ASCIIString, 1}, loci)
    #
    # Find the number of phenotypes corresponding to each locus.
    #
    loc = 0
    for i = 1:rows
      if i == 1 || phenotype_frame[i, :Locus] != phenotype_frame[i - 1, :Locus]
        loc = loc + 1
      end
      phenotypes[loc] = phenotypes[loc] + 1
    end
    #
    # Collect the phenotypes and corresponding genotype strings for
    # each locus in the Phenotype frame.
    #
    loc = 0
    n = 0
    for i = 1:rows
      if i == 1 || phenotype_frame[i, :Locus] != phenotype_frame[i - 1, :Locus]
        loc = loc + 1
        n = 1
        phenotype[loc] = blanks(phenotypes[loc])
        genotype_string[loc] = blanks(phenotypes[loc])
      else
        n = n + 1
      end
      phenotype[loc][n] = string(phenotype_frame[i, :Phenotype])
      genotype_string[loc][n] = string(phenotype_frame[i, :Genotypes])
    end
    #
    # Record the correspondence between the loci in the Phenotype
    # frame and the loci in the Locus frame. The Locus frame may
    # contain codominant loci not in the Phenotype frame.
    #
    correspond = zeros(Int, loci)
    inverse_correspond = zeros(Int, locus.loci)
    for i = 1:loci
      for loc = 1:locus.loci
        if locus.name[loc] == locus_name[i]
          correspond[i] = loc
          inverse_correspond[loc] = i
          break
        end
      end
    end
    #
    # Reduce each phenotype to a set of integer pairs. Each integer
    # represents a numbered allele in the Locus frame.
    #
    genotype = Array(Array{Set{Tuple{Int, Int}}, 1}, loci)
    #
    # Loop over all loci in the Phenotype frame.
    #
    for loc = 1:loci
      cor_loc = correspond[loc]
      genotype[loc] = Array(Set{Tuple{Int, Int}}, phenotypes[loc])
      #
      # Loop over all phenotypes at the current locus.
      #
      for n = 1:phenotypes[loc]
        genotype[loc][n] = Set{Tuple{Int, Int}}()
        #
        # Split the genotype string into genotypes. Commas separate genotypes.
        #
        split_string = split(genotype_string[loc][n], ',')
        #
        # For each genotype consistent with the current phenotype, find
        # its constituent alleles and enter the allele pair into the
        # genotype set.
        #
        for g in split_string
          (double, a1, a2) = fetch_genotype(locus.allele_name[cor_loc], g,
            allele_separator, ordered_allele_separator)
          if a1 == 0
            locname = locus_name[loc]
            throw(ArgumentError(
              "Invalid genotype $g at locus $locname.\n \n"))
          end
          push!(genotype[loc][n], (a1, a2))
          if double; push!(genotype[loc][n], (a2, a1)); end
        end
      end
    end
  end
  #
  # Create the genotype sets for the observed loci.
  #
  observed_loci = locus.loci
  observed_genotype = Array(Set{Tuple{Int, Int}}, people, observed_loci)
  #
  # Loop over all observed loci.
  #
  errors = 0
  for loc = 1:observed_loci
    m = locus.alleles[loc]
    if length(phenotype_frame) > 0
      j = inverse_correspond[loc] # locus in the Phenotype structure.
    end
    #
    # Loop over all people.
    #
    for i = 1:people
      match = false
      if length(phenotype_frame) > 0
        #
        # First consider people with missing values.
        #
        observed_genotype[i, loc] = Set{Tuple{Int, Int}}()
        phen = pedigree_frame[i, locus.locus_field_in_pedigree_frame[loc]]
        if isna(phen) || phen == ""
          hemizygous = locus.xlinked[loc] && male[i]
          observed_genotype[i, loc] = full_genotype_set(m, hemizygous)
          match = true
        #
        # Next consider people with listed phenotypes.
        #
        else
          if j != 0
            phen = convert(ASCIIString, phen)
            for n = 1:phenotypes[j]
              p = phenotype[j][n]
              if isequal(phen, p)
                observed_genotype[i, loc] = genotype[j][n]
                match = true
                break
              end
            end
          end
        end
      end
      #
      # Finally attempt to decode a codominant genotype.
      #
      if !match
        g = pedigree_frame[i, locus.locus_field_in_pedigree_frame[loc]]
        if isna(g)
          g = ""
        else
          g = string(g)
        end
        if g == ""
          if locus.xlinked[loc] && male[i]
            observed_genotype[i, loc] = full_hemizygous_set(m)
          else
            observed_genotype[i, loc] = full_genotype_set(m, false)
          end
        else
          (double, a1, a2) = fetch_genotype(locus.allele_name[loc], g,
            allele_separator, ordered_allele_separator)
          observed_genotype[i, loc] = Set{Tuple{Int, Int}}()
          if a1 == 0
            locsym = locus.name[loc]
            persym = person_name[i]
            throw(ArgumentError(
              "Invalid genotype $g at locus $locsym for person $persym.\n \n"))
          end
          push!(observed_genotype[i, loc], (a1, a2))
          if double; push!(observed_genotype[i, loc], (a2, a1)); end
        end
      end
      if length(observed_genotype[i, loc]) == 0
        errors = errors + 1
        a = person_name[i]
        b = locus.name[loc]
        println("Error: no legal genotypes for person $a at locus $b.")
      end
    end # people loop
  end # observed_locus loop
  #
  # Shut down in the presence of invalid genotypes.
  #
  if errors > 0
    throw(ArgumentError(
      "Mendel terminated due to $errors invalid phenotypes.\n \n"))
  end
  #
  # Identify the ancestral populations and the corresponding
  # admixture fractions.
  #
  population_names = keyword["populations"]
  populations = max(1,length(keyword["populations"]))
  admixture = zeros(Float64, people, populations)
  if populations == 1
    fill!(admixture, 1.0)
  else
    for i = 1:people
      j = 1
      for pop in population_names
        if isna(pedigree_frame[i, symbol(pop)])
          admixture[i, j] = 0.0
        else
          admixture[i, j] = pedigree_frame[i, symbol(pop)]
        end
        j = j + 1
      end
    end
  end
  #
  # Create dummy arrays and return the person structure.
  #
  children = empties(people)
  spouse = empties(people)
  homozygotes = zeros(Int, people, observed_loci)
  person = Person(people, populations, person_name, pedigree_number, mother, 
    father, male, next_twin, primary_twin, admixture, children, spouse,
    observed_genotype, homozygotes, disease_status, variable,
    variable_name)
  return person
end # function person_information

"""
This function counts various features of a pedigree.
Mates is a set of parent pairs in a particular pedigree.
"""
function pedigree_counts!(pedigree::Pedigree, person::Person)

  for i = 1:pedigree.pedigrees
    founders = 0
    females = 0
    males = 0
    twins = 0
    co_twins = 0
    mates = Set{Tuple{Int, Int}}()
    for j = pedigree.start[i]:pedigree.twin_finish[i]
      if person.father[j] == 0
        founders = founders + 1
      else
        push!(mates, (person.mother[j], person.father[j]))
      end
      if person.male[j]
        males = males + 1
      else
        females = females + 1
      end
      if person.primary_twin[j] != 0
        twins = twins + 1
        if person.primary_twin[j] != j; co_twins = co_twins + 1; end
      end
    end
    pedigree.individuals[i] = males + females
    pedigree.founders[i] = founders
    pedigree.females[i] = females
    pedigree.males[i] = males
    pedigree.twins[i] = twins
    pedigree.families[i] = length(mates)
    pedigree.finish[i] = pedigree.twin_finish[i] - co_twins
  end
end # function pedigree_counts!

"""
This function identities all nuclear families, spouses, and children.
The variable k indexes the current nuclear family.
""" 
function construct_nuclear_families(pedigree::Pedigree, person::Person)

  fams = sum(pedigree.families)
  nuclear_family = NuclearFamily(fams, zeros(Int, fams), zeros(Int, fams),
    zeros(Int, fams), empties(fams))
  mother = 0
  father = 0
  k = 0
  #
  # Identify each nuclear family in a pedigree with a pair
  # of parents (mates). Collect the children and spouses of
  # each pedigree member in the process.
  #
  for i = 1:pedigree.pedigrees
    mates = Set{Tuple{Int, Int}}()
    for j = pedigree.start[i]:pedigree.twin_finish[i]
      if person.father[j] != 0
        push!(mates, (person.mother[j], person.father[j]))
        push!(person.children[person.mother[j]], j)
        push!(person.children[person.father[j]], j)
        push!(person.spouse[person.mother[j]], person.father[j])
        push!(person.spouse[person.father[j]], person.mother[j])
      end
    end
    #
    # Loop over the nuclear families and assign parents
    # and siblings.
    #
    for s in mates
      k = k + 1
      nuclear_family.pedigree[k] = i
      (mother, father) = s
      nuclear_family.mother[k] = mother
      nuclear_family.father[k] = father
      SetA = person.children[mother]
      SetB = person.children[father]
      nuclear_family.sib[k] = intersect(SetA, SetB)
    end
  end
  return nuclear_family
end # function construct_nuclear_families

"""
This function completes various data structures and performs
a few checks. Start by eliminating heterozygous genotypes
in a hemizygote, counting homozygotes, and filling in ethnic
admixture fractions.
"""
function check_data_structures!(pedigree::Pedigree, person::Person,
  locus::Locus, nuclear_family::NuclearFamily, keyword::Dict{ASCIIString, Any})
  #
  # Check for pedigree and person errors.
  #
  errors = preliminary_checks(pedigree, person)
  if errors != 0
    throw(ArgumentError(
      "Person or pedigree errors. OpenMendel terminated!\n \n"))
  end
  #
  # Compute the ethnic admixture proportions for non-founders.
  #
  ethnic_admixture!(person)
  #
  # Loop over all observed loci.
  #
  populations = person.populations
  for loc = 1:locus.loci
    #
    # Delete heterozygous genotypes for males at an x-linked locus.
    #
    check_hemizygous_genotypes!(person, locus, loc)
    #
    # Count the number of homozygous genotypes for each person.
    #
    count_homozygotes!(person, loc)
    #
    # Fill in missing allele frequencies if any frequencies are missing.
    #
    alleles = locus.alleles[loc]
    estimate_frequencies = false
    #
    # Check each population.
    #
    for pop = 1:populations
      for allele = 1:alleles
        if isnan(locus.frequency[loc][pop, allele])
          estimate_frequencies = true
          break
        end
      end
      if estimate_frequencies; break; end
    end
    #
    # Estimate missing allele frequencies by gene counting.
    #
    if estimate_frequencies
      pseudo_count = zeros(populations, alleles)
      fill!(pseudo_count, keyword["allele_pseudo_count"])
      gene_counting(person, locus, loc, pseudo_count)
    end
  end
end # function check_data_structures!

"""
This function returns a permutation putting founders at the
head of the pedigree, co-twins at the tail of the pedigree, and
everyone else arranged so that parents precede their children.
It also checks for directed cycles. For algorithmic details see:
Lawler E (1976) Combinatorial Optimization and Matroids. Holt,
Rinehart, and Winston.
"""
function loop(pedigree::Pedigree, father::Vector{Int},
  mother::Vector{Int}, primary_twin::Vector{Int})

  pedigrees = pedigree.pedigrees
  people = length(mother)
  perm = zeros(Int, people)
  #
  # Loop over all pedigrees.
  #
  for ped = 1:pedigrees
    start = pedigree.start[ped]
    twin_finish = pedigree.twin_finish[ped]
    #
    # Put all co-twins at the tail of the permutation.
    #
    n = twin_finish
    for i = start:twin_finish
      if primary_twin[i] != 0 && primary_twin[i] != i
        perm[n] = i
        n = n - 1
      end
    end
    #
    # Put all founders at the head of the permutation. For each
    # remaining person, perm encodes both a current permutation
    # location and how many of his parents have been eliminated as
    # possible candidates for a directed cycle.
    #
    m = start
    j = n
    for i = start:twin_finish
      if primary_twin[i] != 0 && primary_twin[i] != i; continue; end
      if mother[i] == 0
        perm[m] = i
        m = m + 1
      else
        perm[j] = i + 2 * people
        j = j - 1
      end
    end
    #
    # Check whether anyone has both parents eliminated. Such a person
    # cannot belong to a directed cycle. Eliminate this person and
    # compute his final permutation location. If no such person exists,
    # then the remaining people all belong to directed cycles.
    #
    m = start
    while m <= n
      k = 0
      for j = m:n
        if perm[j] <= people
          k = j
          break
        end
      end
      if k == 0
        error_person = mod(perm[m] - 1, people) + 1
        return (perm, error_person)
      end
      #
      # Swap the positions of the two people, and eliminate one of them.
      #
      temp = perm[k]
      perm[k] = perm[m]
      perm[m] = temp
      m = m + 1
      #
      # Find the children of the eliminated person and reduce their
      # Current parent counts by one.
      #
      for i = m:n
        j = mod(perm[i] - 1, people) + 1
        if father[j] == temp; perm[i] = perm[i] - people; end
        if mother[j] == temp; perm[i] = perm[i] - people; end
      end
    end
  end
  return (perm, 0)
end # function loop

"""
This function returns a full set of ordered genotypes.
"""
function full_genotype_set(n::Int, hemizygous::Bool)

  G = Set{Tuple{Int, Int}}()
  if hemizygous
    for i = 1:n
      push!(G, (i, i))
    end
  else
    for i = 1:n
      for j = 1:n
        push!(G, (i, j))
      end
    end
  end
  return G
end # function full_genotype_set

"""
This function recovers the constituent alleles from a genotype.
"""
function fetch_genotype(allele::Vector{ASCIIString}, genotype::AbstractString,
  allele_separator::ASCIIString, ordered_allele_separator::ASCIIString)
  #
  # Split the string genotype into two string alleles.
  #
  a = split(genotype, collect(allele_separator * ordered_allele_separator))
  if length(a) != 2
    return (false, 0, 0)
  end
  #
  # Identify the positions of the two alleles in the input list.
  #
  (a1, a2) = (0, 0)
  for i = 1:length(allele)
    if a[1] == allele[i]
      a1 = i
    end
    if a[2] == allele[i]
      a2 = i
    end
  end
  if a1 == 0 || a2 == 0
    return (false, 0, 0)
  end
  #
  # Decide whether one or two ordered genotypes are pertinent,
  # and return them.
  #
  if a1 == a2 || search(genotype, collect(ordered_allele_separator)) > 0
    return (false, a1, a2) # One genotype.
  else
    return (true, a1, a2) # Two genotypes.
  end
end # function fetch_genotype

"""
This function does a few pedigree and person checks.
"""
function preliminary_checks(pedigree::Pedigree, person::Person)

  errors = 0
  #
  # Check for repeated pedigree names.
  #
  (nonunique, repeat) = repeated_string(pedigree.name)
  if nonunique
    errors = errors + 1
    println("Error: Pedigree name $repeat is used for multiple pedigrees.")
  end
  #
  # Check for repeated person names.
  #
  for ped = 1:pedigree.pedigrees
    ped_start = pedigree.start[ped]
    ped_finish = pedigree.twin_finish[ped]
    (nonunique, repeat) = repeated_string(person.name[ped_start:ped_finish])
    if nonunique
      errors = errors + 1
      p = pedigree.name[ped]
      println("Error: Person name $repeat in pedigree $p is not unique.")
    end
  end
  #
  # Check for blank person ids.
  #
  for i = 1:person.people
    if person.name[i] == ""
      errors = errors + 1
      println("Error: Pedigree $person.pedigree[i] has a blank person name.")
    end
    #
    # Check for missing parents, parents of the wrong sex, or parents
    # the same as children.
    #
    mother_present = person.mother[i] != 0
    father_present = person.father[i] != 0
    if mother_present != father_present
      errors = errors + 1
      println("Error: Person $person.name[i] lacks one parent.")
    end
    if mother_present
      if person.male[person.mother[i]]
        errors = errors + 1
        mother_name = person.name[person.mother[i]]
        println("Error: The mother $mother_name is a male.")
      end
      if person.name[i] == person.name[person.mother[i]]
        errors = errors + 1
        println("Error: Person $person.name[i] is her own mother.")
      end
    end
    if father_present
      if !person.male[person.father[i]]
        errors = errors + 1
        father_name = person.name[person.father[i]]
        println("Error: The father $father_name is a female.")
      end
      if person.name[i] == person.name[person.father[i]]
        errors = errors + 1
        println("Error: Person $person.name[i] is his own father.")
      end
    end
    #
    # Check for twin inconsistencies.
    #
    if person.primary_twin[i] != 0
      j = person.next_twin[person.primary_twin[i]]
      if j == 0
        errors = errors + 1
        println("Error: The twin $person.name[i]] lacks identical siblings.")
      end
      while j != 0
        if person.male[i] != person.male[j]
          errors = errors + 1
          a = person.name[i]
          b = person.name[j]
          println("Error: Identical twins $a and $b are of opposite sex.")
        end
        j = person.next_twin[j]
      end
    end
  end
  return errors
end # function preliminary_checks

"""
This function deletes heterozygous genotypes ostensibly compatible
with x-linked male phenotypes.
"""
function check_hemizygous_genotypes!(person::Person, locus::Locus, loc::Int)

  if locus.xlinked[loc]
    for i = 1:person.people
      if !person.male[i]; continue; end
      for (gm, gp) in person.genotype[i, loc]
        if gm != gp
          delete!(person.genotype[i, loc], (gm, gp))
        end
      end
    end
  end
end # function check_hemizygous_genotypes!

"""
This function counts the number of homozygous genotypes at a locus.
"""
function count_homozygotes!(person::Person, loc::Int)

  fill!(person.homozygotes[:, loc], 0)
  for i = 1:person.people
    for (gm, gp) in person.genotype[i, loc]
      if gm == gp
        person.homozygotes[i, loc] = person.homozygotes[i, loc] + 1
      end
    end
  end
end # function count_homozygotes!

"""
This function computes the ethnic admixture fractions
for all non-founders. Warning: Parents must come before
children in each pedigree for averaging to work.
"""
function ethnic_admixture!(person::Person)

  for i = 1:person.people
    if person.mother[i] != 0
      a = person.admixture[person.mother[i], :]
      b = person.admixture[person.father[i], :]
      person.admixture[i, :] = 0.5 * (a + b)
    else
      for j = 1:length(person.admixture[i, :])
        if isna(person.admixture[i, j])
          person.admixture[i, j] = 0.0
        end
      end
      s = sum(person.admixture[i, :])
      person.admixture[i, :] = person.admixture[i, :] / s
    end
  end
end # function ethnic_admixture!

"""
This function estimates allele frequencies at locus loc
by gene counting (EM algorithm). People are treated as
unrelated. A pseudocount is required for each population
allele pair. Warning: Estimation should precede genotype
elimination and allele consolidation.
"""
function gene_counting(person::Person, locus::Locus, loc::Int,
                       pseudocount::Matrix{Float64})

  people = person.people
  xlinked = locus.xlinked[loc]
  populations = person.populations
  alleles = size(locus.frequency[loc], 2) - 1
  locus.frequency[loc][:, 1:alleles] = 1.0 / alleles
  maximum_genotypes = alleles^2
  #
  # In the EM loop, initialize allele counts with pseudocount.
  #
  for iteration = 1:1000
    allele_count = copy(pseudocount)
    #
    # Exclude people with no genotype information.
    #
    for i = 1:people
      if xlinked && person.male[i]
        if length(person.genotype[i, loc]) == alleles
          continue
        end
      else
        if length(person.genotype[i, loc]) == maximum_genotypes
          continue
        end
      end
      #
      # Compute the normalizing weight for the allele counts for
      # the current person. Fill mixed with composite allele
      # frequencies based on fractional ancestries and independence.
      #
      mixed = vec(person.admixture[i, :])' * locus.frequency[loc][:, 1:alleles]
      weight = 0.0
      if xlinked && person.male[i]
        for (gm, gp) in person.genotype[i, loc]
          weight = weight + mixed[gm]
        end
      else
        for (gm, gp) in person.genotype[i, loc]
          fm = mixed[gm]
          fp = mixed[gp]
          weight = weight + fm * fp
        end
      end
      #
      # Update the fractional allele counts.
      #
      if xlinked && person.male[i]
        for (gm, gp) in person.genotype[i, loc]
          for pop = 1:populations
            c = person.admixture[i, pop] * locus.frequency[loc][pop, gm]
            allele_count[pop, gm] = allele_count[pop, gm] + c / weight
          end
        end
      else
        for (gm, gp) in person.genotype[i, loc]
          fm = mixed[gm]
          fp = mixed[gp]
          for pop = 1:populations
            cm = person.admixture[i, pop] * locus.frequency[loc][pop, gm]
            allele_count[pop, gm] = allele_count[pop, gm] + cm * fp / weight
            cp = person.admixture[i, pop] * locus.frequency[loc][pop, gp]
            allele_count[pop, gp] = allele_count[pop, gp] + cp * fm / weight
          end
        end
      end
    end
    #
    # Update the EM estimates, and check for convergence.
    #
    frequency = copy(allele_count)
    for pop = 1:populations
      frequency[pop, :] = copy(frequency[pop, :] / sum(frequency[pop, :]))
    end
    if norm(frequency - locus.frequency[loc][:, 1:alleles]) < 1e-5
      locus.frequency[loc][:, 1:alleles] = frequency
      break
    else
      locus.frequency[loc][:, 1:alleles] = frequency
    end
  end # EM iteration loop
  return nothing
end # function gene_counting

"""
This function checks that ancestral populations are present in both the
pedigree and locus frames.
"""
function check_populations(locus_frame::DataFrame, pedigree_frame::DataFrame,
                           keyword::Dict{ASCIIString, Any})

  pedigree_field = names(pedigree_frame)
  for pop in keyword["populations"]
    if symbol(pop) in pedigree_field
      continue
    else
      throw(ArgumentError("Population $pop is not in the pedigree frame."))
    end
  end
  if length(locus_frame) != 0
    locus_field = names(locus_frame)
    for pop in keyword["populations"]
      if symbol(pop) in locus_field
        continue
      else
        throw(ArgumentError("Population $pop is not in the locus frame."))
      end
    end
  end
end # function check_populations

