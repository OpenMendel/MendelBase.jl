################################################################################
# Reads the pedigree and locus data from external files
# and then constructs the appropriate data frames and data structures.
################################################################################
#
# Required external modules.
#
# using CSV
# using DataFrames
# using Missings
#
# Required OpenMendel packges and modules.
#
# using DataStructures
# using GeneralUtilities
# using GeneticUtilities
# using SnpArrays

export count_homozygotes!, read_external_data_files

"""
This function organizes reading in the data from external files.
All names of data files are stored in the relevant keyword values.
"""
function read_external_data_files(keyword::Dict{AbstractString, Any})
  #
  # Recall the field separator (AKA delimiter) used in the data files.
  # The default is a comma but can be changed to another character. 
  #
  field_sep = keyword["field_separator"]
  #
  # Next, set the missing-value strings used in the data files.
  # The primary default is "NA" but this can be changed to
  # another string or set of strings using the keyword "missing_value".
  # In addition, when using a delimiter other than ' ',
  # an empty field is also labeled missing, e.g., field1,,field3
  # (note that empty means not even a space between the delimiters).
  # When the delimiter is ' ', we treat multiple spaces as a single delimiter.
  #
  null_strings = convert(Array{AbstractString, 1},
                         union([], keyword["missing_value"]))
  allow_padding = (field_sep == ' ')
  #
  # Read the data from the appropriate files. Put the data into dataframes.
  # All text data files must have one header row as the first line,
  # except Plink-format files which must have no header rows.
  # By default, all columns are typed as allowing missing values.
  #
  # Start by reading the locus file.
  #
  locus_file = string(keyword["locus_file"])
  if locus_file == ""
    locus_frame = DataFrame() # An empty dataframe
  else
    locus_frame = CSV.File(locus_file; header = 1,
      delim = field_sep, ignorerepeated = allow_padding,
      missingstrings = null_strings) |> DataFrame
  end
  #
  # Read the phenotype data into a frame.
  #
  phenotype_file = string(keyword["phenotype_file"])
  if phenotype_file == ""
    phenotype_frame = DataFrame() # An empty dataframe
  else
    phenotype_frame = CSV.File(phenotype_file; header = 1,
      delim = field_sep, ignorerepeated = allow_padding,
      missingstrings = null_strings) |> DataFrame
  end
  #
  # Read the data in the mandatory pedigree file into a frame.
  # Plink .fam files are allowed.
  # Other pedigree files should contain a header line.
  #
  pedigree_file = string(keyword["pedigree_file"])
  if occursin(".fam", pedigree_file)
    person_frame = read_plink_fam_file(pedigree_file, keyword)
  else
    person_frame = CSV.File(pedigree_file; header = 1,
      delim = field_sep, ignorerepeated = allow_padding,
      missingstrings = null_strings) |> DataFrame
  end
  #
  # Add a column recording the order of entry of each person.
  #
  person_frame.EntryOrder =
    collect(Union{Int64, Missings.Missing}, 1:size(person_frame, 1))
  #
  # Check that ancestral populations are present in both the
  # pedigree and locus frames.
  #
  check_populations(locus_frame, person_frame, keyword)
  #
  # Assemble the locus, pedigree, and person data structures.
  #
  locus = locus_information(locus_frame, person_frame, keyword)
  pedigree = pedigree_information(person_frame)
  person = person_information(locus_frame, person_frame, phenotype_frame,
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
      "Either the SNP definition or SNP data file was not specified.\n \n"))
  end
  snpdefinition_file = string(keyword["snpdefinition_file"])
  if snpdefinition_file == ""
    #
    # Construct empty SNP data structures.
    #
    snp_definition_frame = DataFrame() # An empty dataframe.
    snpmatrix = SnpArray(undef, 0, 0)
    snpdata = snp_information(snp_definition_frame, person, snpmatrix, keyword)
  else
    if occursin(".bim", snpdefinition_file)
      snp_definition_frame = read_plink_bim_file(snpdefinition_file, keyword)
    else
      snp_definition_frame = CSV.File(snpdefinition_file; header = 1,
      delim = field_sep, ignorerepeated = allow_padding,
      missingstrings = null_strings) |> DataFrame
    end
    #
    # Read the SNP bed file.
    #
    snpdata_file = string(keyword["snpdata_file"])
    snpmatrix = SnpArray(snpdata_file, person.people)
    #
    # Assemble the SNP data structure.
    #
    snpdata = snp_information(snp_definition_frame, person, snpmatrix, keyword)
    snpdata.maf = maf(snpmatrix)
    snpdata.minor_allele = minorallele(snpmatrix)
    snpdata.missings_per_snp = missingrate(snpmatrix, 1)
    snpdata.missings_per_person = missingrate(snpmatrix, 2)
  end
  return (pedigree, person, nuclear_family, locus, snpdata,
    locus_frame, phenotype_frame, person_frame, snp_definition_frame)
end # function read_external_data_files

"""
Read a Plink .fam file into a dataframe and 
substitutes blanks for missing values.
"""
function read_plink_fam_file(plink_fam_file::AbstractString,
  keyword::Dict{AbstractString, Any})
  #
  # Recall the field separator (AKA delimiter) specified for Plink data files.
  # The default delimiter is a space but can be changed to another character.
  #
  plink_field_sep = keyword["plink_field_separator"]
  #
  # Next, set the missing-value strings used in the data files.
  # The primary default is "NA" but this can be changed to
  # another string or set of strings using the keyword "plink_missing_value".
  # In addition, when using a delimiter other than ' ',
  # an empty field is also labeled missing, e.g., field1,,field3
  # (note that empty means not even a space between the delimiters).
  # When the delimiter is ' ', we treat multiple spaces as a single delimiter.
  # Note that some missing values are fixed by the Plink syntax
  # and these are accounted for in the code below.
  #
  plink_null_strings = convert(Array{AbstractString, 1},
                               union([], keyword["plink_missing_value"]))
  allow_padding = (plink_field_sep == ' ')
  #
  # By default, all columns are typed as allowing missing values.
  #
  column_types = [String, String, String, String, String, Float64]
  column_names = [:Pedigree, :Person, :Father, :Mother, :Sex, :Trait]
  fam_dframe = CSV.File(plink_fam_file;
    delim = plink_field_sep, ignorerepeated = allow_padding,
    types = column_types, header = column_names,
    missingstrings = plink_null_strings) |> DataFrame
  allowmissing!(fam_dframe)
###
### Read the fam data into a matrix of strings.
### The package DelimitedFiles contains readdlm().
###
##fam_matrix = readdlm(plink_fam_file, String)
##rows = size(fam_matrix, 1)
###
### Move the data from the fam matrix to a new fam dataframe.
###
##fam_dframe = DataFrame()
##fam_dframe[!, :Pedigree] = fam_matrix[:, 1]
##fam_dframe[!, :Person] = fam_matrix[:, 2]
##fam_dframe[!, :Father] = fam_matrix[:, 3]
##fam_dframe[!, :Mother] = fam_matrix[:, 4]
##fam_dframe[!, :Sex] = fam_matrix[:, 5]
##fam_dframe[!, :Trait] = map(x -> parse(Float64, x), fam_matrix[:, 6])
###
### Allow for missing data in some columns.
###
##allowmissing!(fam_dframe, (:Father, :Mother, :Trait))
###  fam_dframe[!, :Father] =
###    convert(Array{Union{String, Missings.Missing}, 1}, fam_dframe[:Father])
###  fam_dframe[!, :Mother] =
###    convert(Array{Union{String, Missings.Missing}, 1}, fam_dframe[:Mother])
###  fam_dframe[!, :Trait] =
###    convert(Array{Union{Float64, Missings.Missing}, 1}, fam_dframe[:Trait])
  #
  # Use the Julia-internal value "missing" for the data labeled with
  # the Plink-format missing value symbols.
  # Note: the trait value zero is not changed to missing!
  #
  for i = 1:size(fam_dframe, 1)
    if fam_dframe[i, :Father] == "0"
      fam_dframe[i, :Father] = missing
    end
    if fam_dframe[i, :Mother] == "0"
      fam_dframe[i, :Mother] = missing
    end
    if fam_dframe[i, :Trait] == -9 || fam_dframe[i, :Trait] == NaN
      fam_dframe[i, :Trait] = missing
    end
  end
  return fam_dframe
end # function read_plink_fam_file

"""
Read a Plink .bim file into a dataframe.
"""
function read_plink_bim_file(plink_bim_file::AbstractString,
  keyword::Dict{AbstractString, Any})
  #
  # Recall the field separator (AKA delimiter) specified for Plink data files.
  # The default delimiter is a space but can be changed to another character.
  #
  plink_field_sep = keyword["plink_field_separator"]
  #
  # Next, set the missing-value strings used in the data files.
  # The primary default is "NA" but this can be changed to
  # another string or set of strings using the keyword "plink_missing_value".
  # In addition, when using a delimiter other than ' ',
  # an empty field is also labeled missing, e.g., field1,,field3
  # (note that empty means not even a space between the delimiters).
  # When the delimiter is ' ', we treat multiple spaces as a single delimiter.
  #
  plink_null_strings = convert(Array{AbstractString, 1},
                               union([], keyword["plink_missing_value"]))
  allow_padding = (plink_field_sep == ' ')
  #
  # By default, all columns are typed as allowing missing values.
  #
  column_types = [String, String, Float64, Int64, String, String]
  column_names = [:Chromosome, :SNP, :CentiMorgans, :Basepairs, 
                  :Allele1, :Allele2]
  bim_dframe = CSV.File(plink_bim_file;
    delim = plink_field_sep, ignorerepeated = allow_padding,
    types = column_types, header = column_names,
    missingstrings = plink_null_strings) |> DataFrame
  allowmissing!(bim_dframe)
###
### Read the bim data into a matrix.
### The package DelimitedFiles contains readdlm().
###
##bim_matrix = readdlm(plink_bim_file, String)
##rows = size(bim_matrix, 1)
###
### Move the data from the bim matrix to a new bim dataframe.
###
##bim_dframe = DataFrame()
##bim_dframe[!, :Chromosome] = bim_matrix[:, 1]
##bim_dframe[!, :SNP] = bim_matrix[:, 2]
##bim_dframe[!, :CentiMorgans] = map(x -> parse(Float64, x), bim_matrix[:, 3])
##bim_dframe[!, :Basepairs] = map(x -> parse(Int64, x), bim_matrix[:, 4])
##bim_dframe[!, :Allele1] = bim_matrix[:, 5]
##bim_dframe[!, :Allele2] = bim_matrix[:, 6]
###
### Allow for missing data in some columns.
###
##allowmissing!(bim_dframe,
### (:Chromosome, :CentiMorgans, :Basepairs, :Allele1, :Allele2))
###  bim_dframe[!, :Chromosome] =
###  convert(Array{Union{String, Missings.Missing}, 1}, bim_dframe[:Chromosome])
###  bim_dframe[!, :CentiMorgans] =
###convert(Array{Union{Float64, Missings.Missing}, 1},bim_dframe[:CentiMorgans])
###  bim_dframe[!, :Basepairs] =
###   convert(Array{Union{Int64, Missings.Missing}, 1}, bim_dframe[:Basepairs])
###  bim_dframe[!, :Allele1] =
###    convert(Array{Union{String, Missings.Missing}, 1}, bim_dframe[:Allele1])
###  bim_dframe[!, :Allele2] =
###    convert(Array{Union{String, Missings.Missing}, 1}, bim_dframe[:Allele2])
  return bim_dframe
end # function read_plink_bim_file

"""
Copies information from a SNP definition dataframe into
the SNP data structure.
"""
function snp_information(snp_definition_frame::DataFrame, person::Person,
  snpmatrix::AbstractSnpArray, keyword::Dict{AbstractString, Any})
#
  if keyword["snpdefinition_file"] == ""
    snpdata = SnpDataStruct(0, 0, Vector{AbstractString}(),
      Vector{AbstractString}(), Vector{AbstractString}(), Vector{Float64}(),
      Vector{Int}(), Vector{AbstractString}(), Vector{AbstractString}(),
      Vector{Float64}(), BitVector(), snpmatrix, Vector{Int}(), Vector{Int}())
    return snpdata
  end
  (snps, columns) = size(snp_definition_frame)
  #
  # Strip any leading or trailing spaces from the field names.
  #
  for i in names(snp_definition_frame)
    rename!(snp_definition_frame, i => Symbol(strip(string(i))))
  end
  #
  # Fill-in the snp name array, which is not allowed to have missing values.
  #
  column_names = names(snp_definition_frame)
  if !(:SNP in column_names)
    if :Locus in column_names
      rename!(snp_definition_frame, :Locus => :SNP)
    elseif :Loci in column_names
      rename!(snp_definition_frame, :Loci => :SNP)
    elseif :locus in column_names
      rename!(snp_definition_frame, :locus => :SNP)
    elseif :loci in column_names
      rename!(snp_definition_frame, :loci => :SNP)
    elseif :Snp in column_names
      rename!(snp_definition_frame, :Snp => :SNP)
    elseif :snp in column_names
      rename!(snp_definition_frame, :snp => :SNP)
    end
    column_names = names(snp_definition_frame)
  end

  snp_name = blanks(snps)
  if :SNP in column_names
    for snp = 1:snps
      snp_string = string(snp)
      if ismissing(snp_definition_frame[snp, :SNP])
        throw(ArgumentError(
          "SNP $snp_string of the input data has no name.\n \n"))
      end
      if typeof(snp_definition_frame[snp, :SNP]) <: Number
        if snp_definition_frame[snp, :SNP] == NaN
          throw(ArgumentError(
            "SNP $snp_string of the input data has no recognizable name.\n \n"))
        end
      elseif typeof(snp_definition_frame[snp, :SNP]) <: AbstractString
        snp_definition_frame[snp, :SNP] = strip(snp_definition_frame[snp, :SNP])
        if snp_definition_frame[snp, :SNP] == ""
          throw(ArgumentError(
            "SNP $snp_string of the input data has no name.\n \n"))
        end
      else
        throw(ArgumentError(
          "SNP $snp_string of the input data has no recognizable name.\n \n"))
      end
      snp_name[snp] = string(snp_definition_frame[snp, :SNP])
    end
  else
    for snp = 1:snps
      snp_name[snp] = string(snp) # label SNPs 1, 2, 3, etc.
    end
  end
  #
  # Fill-in the chromosome array, which is allowed to have missing values.
  #
  if !(:Chromosome in column_names)
    if :chromosome in column_names
      rename!(snp_definition_frame, :chromosome => :Chromosome)
    elseif :Chr in column_names
      rename!(snp_definition_frame, :Chr => :Chromosome)
    elseif :chr in column_names
      rename!(snp_definition_frame, :chr => :SNP)
    end
    column_names = names(snp_definition_frame)
  end

  chromosome = blanks(snps)
  if :Chromosome in column_names
    for snp = 1:snps
      if ismissing(snp_definition_frame[snp, :Chromosome])
        chromosome[snp] = "autosome"
        continue
      end
      if typeof(snp_definition_frame[snp, :Chromosome]) <: Integer
        if snp_definition_frame[snp, :Chromosome] == NaN
          chromosome[snp] = "autosome"
          continue
        end
      elseif typeof(snp_definition_frame[snp, :Chromosome]) <: AbstractString
        snp_definition_frame[snp, :Chromosome] =
          strip(snp_definition_frame[snp, :Chromosome])
        if snp_definition_frame[snp, :Chromosome] == ""
          chromosome[snp] = "autosome"
          continue
        end
      else
        chromosome[snp] = "autosome"
        continue
      end
      chromosome[snp] = string(snp_definition_frame[snp, :Chromosome])
    end
  else
    fill!(chromosome, "autosome")
  end
  #
  # Fill-in the centimorgan array, which is allowed to have missing values.
  #
  if !(:CentiMorgans in column_names)
    if :CentiMorgan in column_names
      rename!(snp_definition_frame, :CentiMorgan => :CentiMorgans)
    elseif :Centimorgans in column_names
      rename!(snp_definition_frame, :Centimorgans => :CentiMorgans)
    elseif :Centimorgan in column_names
      rename!(snp_definition_frame, :Centimorgan => :CentiMorgans)
    elseif :centiMorgans in column_names
      rename!(snp_definition_frame, :centiMorgans => :CentiMorgans)
    elseif :centiMorgan in column_names
      rename!(snp_definition_frame, :centiMorgan => :CentiMorgans)
    elseif :centimorgans in column_names
      rename!(snp_definition_frame, :centimorgans => :CentiMorgans)
    elseif :centimorgan in column_names
      rename!(snp_definition_frame, :centimorgan => :CentiMorgans)
    elseif :CM in column_names
      rename!(snp_definition_frame, :CM => :CentiMorgans)
    elseif :cM in column_names
      rename!(snp_definition_frame, :cM => :CentiMorgans)
    elseif :cm in column_names
      rename!(snp_definition_frame, :cm => :CentiMorgans)
    end
    column_names = names(snp_definition_frame)
  end

  centimorgans = zeros(Float64, snps)
  if :CentiMorgans in column_names
    for snp = 1:snps
      if ismissing(snp_definition_frame[snp, :CentiMorgans])
        continue
      end
      if typeof(snp_definition_frame[snp, :CentiMorgans]) <: Real
        if snp_definition_frame[snp, :CentiMorgans] != NaN
          centimorgans[snp] =
            convert(Float64, snp_definition_frame[snp, :CentiMorgans])
        end
      elseif typeof(snp_definition_frame[snp, :CentiMorgans]) <: AbstractString
        snp_definition_frame[snp, :CentiMorgans] =
          strip(snp_definition_frame[snp, :CentiMorgans])
        cm_string = snp_definition_frame[snp, :CentiMorgans]
        if isa(Meta.parse(cm_string, raise = false), Real)
          centimorgans[snp] = convert(Float64, Meta.parse(cm_string))
        end
      end
    end
  end
  #
  # Fill-in the basepairs array, which is allowed to have missing values.
  #
  if !(:Basepairs in column_names)
    if :Basepair in column_names
      rename!(snp_definition_frame, :Basepair => :Basepairs)
    elseif :BasePairs in column_names
      rename!(snp_definition_frame, :BasePairs => :Basepairs)
    elseif :BasePair in column_names
      rename!(snp_definition_frame, :BasePair => :Basepairs)
    elseif :basepairs in column_names
      rename!(snp_definition_frame, :basepairs => :Basepairs)
    elseif :basepair in column_names
      rename!(snp_definition_frame, :basepair => :Basepairs)
    elseif :BP in column_names
      rename!(snp_definition_frame, :BP => :Basepairs)
    elseif :Bp in column_names
      rename!(snp_definition_frame, :Bp => :Basepairs)
    elseif :bp in column_names
      rename!(snp_definition_frame, :bp => :Basepairs)
    end
    column_names = names(snp_definition_frame)
  end

  basepairs = zeros(Int, snps)
  if :Basepairs in column_names
    for snp = 1:snps
      if ismissing(snp_definition_frame[snp, :Basepairs])
        continue
      end
      if typeof(snp_definition_frame[snp, :Basepairs]) <: Integer
        if snp_definition_frame[snp, :Basepairs] != NaN
          basepairs[snp] =
            convert(Int, snp_definition_frame[snp, :Basepairs])
        end
      elseif typeof(snp_definition_frame[snp, :Basepairs]) <: AbstractString
        snp_definition_frame[snp, :Basepairs] =
          strip(snp_definition_frame[snp, :Basepairs])
        bp_string = snp_definition_frame[snp, :Basepairs]
        if isa(Meta.parse(bp_string, raise = false), Integer)
          basepairs[snp] = convert(Int, Meta.parse(bp_string))
        end
      end
    end
  end
  #
  # Fill-in the allele arrays, which are allowed to have missing values.
  #
  if !(:Allele1 in column_names)
    if :allele1 in column_names
      rename!(snp_definition_frame, :allele1 => :Allele1)
    end
    column_names = names(snp_definition_frame)
  end
  if !(:Allele2 in column_names)
    if :allele2 in column_names
      rename!(snp_definition_frame, :allele2 => :Allele2)
    end
    column_names = names(snp_definition_frame)
  end

  allele1 = blanks(snps)
  if :Allele1 in column_names
    for snp = 1:snps
      if ismissing(snp_definition_frame[snp, :Allele1])
        allele1[snp] = "1"
        continue
      end
      if typeof(snp_definition_frame[snp, :Allele1]) <: Integer
        if snp_definition_frame[snp, :Allele1] == NaN ||
           snp_definition_frame[snp, :Allele1] == 0
          allele1[snp] = "1"
          continue
        end
      elseif typeof(snp_definition_frame[snp, :Allele1]) <: AbstractString
        snp_definition_frame[snp, :Allele1] =
          strip(snp_definition_frame[snp, :Allele1])
        if snp_definition_frame[snp, :Allele1] == ""
          allele1[snp] = "1"
          continue
        end
      else
        allele1[snp] = "1"
        continue
      end
      allele1[snp] = string(snp_definition_frame[snp, :Allele1])
    end
  else
    fill!(allele1, "1")
  end
  allele2 = blanks(snps)
  if :Allele2 in column_names
    for snp = 1:snps
      if ismissing(snp_definition_frame[snp, :Allele2])
        allele2[snp] = "2"
        continue
      end
      if typeof(snp_definition_frame[snp, :Allele2]) <: Integer
        if snp_definition_frame[snp, :Allele2] == NaN ||
           snp_definition_frame[snp, :Allele2] == 0
          allele2[snp] = "2"
          continue
        end
      elseif typeof(snp_definition_frame[snp, :Allele2]) <: AbstractString
        snp_definition_frame[snp, :Allele2] =
          strip(snp_definition_frame[snp, :Allele2])
        if snp_definition_frame[snp, :Allele2] == ""
          allele2[snp] = "2"
          continue
        end
      else
        allele2[snp] = "2"
        continue
      end
      allele2[snp] = string(snp_definition_frame[snp, :Allele2])
    end
  else
    fill!(allele2, "2")
  end
  #   
  # Return the SNP data structure.
  #   
  people = person.people
  personid = person.name
  snpdata = SnpDataStruct(people, snps, personid, snp_name, chromosome,
    centimorgans, basepairs, allele1, allele2, Vector{Float64}(), BitVector(), 
    snpmatrix, Vector{Int}(), Vector{Int}())
  return snpdata 
end # function snp_information

"""
Extracts locus information from the locus frame.
First extract relevant dimensions and fields.
"""
function locus_information(locus_frame::DataFrame, person_frame::DataFrame,
                           keyword::Dict{AbstractString, Any})
  #
  # If the locus_frame is empty, create a null locus structure.
  #
  if size(locus_frame, 2) == 0
    a = Array{Array{AbstractString, 1}}(undef, 1); a[1] = blanks(1)
    b = Array{Array{Float64, 2}}(undef, 1); b[1] = zeros(1, 1)
    c = zeros(1, 1, 1)
    locus = Locus(0, 0, 0, true, 0, 0, 0, zeros(2), blanks(0), blanks(0),
                  zeros(Int, 0), zeros(2, 0), zeros(2, 0), trues(0),
                  zeros(Int, 0), a, b, c, zeros(Int, 0), zeros(Int, 0))
    return locus
  end
  #
  # Strip any leading or trailing spaces from the field names.
  #
  for i in names(locus_frame)
    rename!(locus_frame, i => Symbol(strip(string(i))))
  end
  #
  # Fix some possible field naming issues.
  #
  locus_field = names(locus_frame)
  if !(:FemaleMorgans in locus_field)
    if :FemaleMorgan in locus_field
      rename!(locus_frame, :FemaleMorgan => :FemaleMorgans)
    elseif :Femalemorgans in locus_field
      rename!(locus_frame, :Femalemorgans => :FemaleMorgans)
    elseif :Femalemorgan in locus_field
      rename!(locus_frame, :Femalemorgan => :FemaleMorgans)
    elseif :femaleMorgans in locus_field
      rename!(locus_frame, :femaleMorgans => :FemaleMorgans)
    elseif :femaleMorgan in locus_field
      rename!(locus_frame, :femaleMorgan => :FemaleMorgans)
    elseif :femalemorgans in locus_field
      rename!(locus_frame, :femalemorgans => :FemaleMorgans)
    elseif :femalemorgan in locus_field
      rename!(locus_frame, :femalemorgan => :FemaleMorgans)
    elseif :Morgans in locus_field
      rename!(locus_frame, :Morgans => :FemaleMorgans)
    elseif :Morgan in locus_field
      rename!(locus_frame, :Morgan => :FemaleMorgans)
    elseif :morgans in locus_field
      rename!(locus_frame, :morgans => :FemaleMorgans)
    elseif :morgan in locus_field
      rename!(locus_frame, :morgan => :FemaleMorgans)
    end
  end
  if !(:MaleMorgans in locus_field)
    if :MaleMorgan in locus_field
      rename!(locus_frame, :MaleMorgan => :MaleMorgans)
    elseif :Malemorgans in locus_field
      rename!(locus_frame, :Malemorgans => :MaleMorgans)
    elseif :Malemorgan in locus_field
      rename!(locus_frame, :Malemorgan => :MaleMorgans)
    elseif :maleMorgans in locus_field
      rename!(locus_frame, :maleMorgans => :MaleMorgans)
    elseif :maleMorgan in locus_field
      rename!(locus_frame, :maleMorgan => :MaleMorgans)
    elseif :malemorgans in locus_field
      rename!(locus_frame, :malemorgans => :MaleMorgans)
    elseif :malemorgan in locus_field
      rename!(locus_frame, :malemorgan => :MaleMorgans)
    end
  end
  if !(:Basepairs in locus_field)
    if :Basepair in locus_field
      rename!(locus_frame, :Basepair => :Basepairs)
    elseif :BasePairs in locus_field
      rename!(locus_frame, :BasePairs => :Basepairs)
    elseif :BasePair in locus_field
      rename!(locus_frame, :BasePair => :Basepairs)
    elseif :basepairs in locus_field
      rename!(locus_frame, :basepairs => :Basepairs)
    elseif :basepair in locus_field
      rename!(locus_frame, :basepair => :Basepairs)
    elseif :BP in locus_field
      rename!(locus_frame, :BP => :Basepairs)
    elseif :Bp in locus_field
      rename!(locus_frame, :Bp => :Basepairs)
    elseif :bp in locus_field
      rename!(locus_frame, :bp => :Basepairs)
    end
  end
  if !(:Locus in locus_field)
    if :locus in locus_field
      rename!(locus_frame, :locus => :Locus)
    elseif :loci in locus_field
      rename!(locus_frame, :loci => :Locus)
    elseif :Loci in locus_field
      rename!(locus_frame, :Loci => :Locus)
    elseif :SNP in locus_field
      rename!(locus_frame, :SNP => :Locus)
    elseif :Snp in locus_field
      rename!(locus_frame, :Snp => :Locus)
    elseif :snp in locus_field
      rename!(locus_frame, :snp => :Locus)
    end
  end
  if !(:Chromosome in locus_field)
    if :chromosome in locus_field
      rename!(locus_frame, :chromosome => :Chromosome)
    elseif :Chr in locus_field
      rename!(locus_frame, :Chr => :Chromosome)
    elseif :chr in locus_field
      rename!(locus_frame, :chr => :Chromosome)
    end
  end
  if !(:Allele in locus_field)
    if :allele in locus_field
      rename!(locus_frame, :allele => :Allele)
    end
  end
  #
  # Check that certain fields are in the locus_frame
  # and of the proper type.
  #
  locus_field = names(locus_frame)
  if !(:Locus in locus_field)
    throw(ArgumentError(
      "The locus file must have a Locus or SNP field.\n \n"))
  end
  if !(:Allele in locus_field)
    throw(ArgumentError(
      "The locus file must have an Allele field.\n \n"))
  end
  if !(:Chromosome in locus_field)
    throw(ArgumentError(
      "The locus file must have a Chromosome field.\n \n"))
  end
  if :Basepairs in locus_field
    if typeof(locus_frame[:, :Basepairs]) !=
         Array{Union{Int64, Missings.Missing}, 1} &&
       typeof(locus_frame[:, :Basepairs]) != Array{Int64, 1}
      throw(ArgumentError(
      "The Basepairs column should only have integers or missing values.\n \n"))
    end
  end
  if :Morgans in locus_field
    if typeof(locus_frame[:, :Morgans]) !=
         Array{Union{Float64, Missings.Missing}, 1} &&
       typeof(locus_frame[:, :Morgans]) != Array{Float64, 1} &&
       typeof(locus_frame[:, :Morgans]) !=
         Array{Union{Int64, Missings.Missing}, 1} &&
       typeof(locus_frame[:, :Morgans]) != Array{Int64, 1}
      throw(ArgumentError(
        "All Morgans columns should only have numbers or missing values.\n \n"))
    end
  end
  if :FemaleMorgans in locus_field
    if typeof(locus_frame[:, :FemaleMorgans]) !=
         Array{Union{Float64, Missings.Missing}, 1} &&
       typeof(locus_frame[:, :FemaleMorgans]) != Array{Float64, 1} &&
       typeof(locus_frame[:, :FemaleMorgans]) !=
         Array{Union{Int64, Missings.Missing}, 1} &&
       typeof(locus_frame[:, :FemaleMorgans]) != Array{Int64, 1}
      throw(ArgumentError(
        "All Morgans columns should only have numbers or missing values.\n \n"))
    end
  end
  if :MaleMorgans in locus_field
    if typeof(locus_frame[:, :MaleMorgans]) !=
         Array{Union{Float64, Missings.Missing}, 1} &&
       typeof(locus_frame[:, :MaleMorgans]) != Array{Float64, 1} &&
       typeof(locus_frame[:, :MaleMorgans]) !=
         Array{Union{Int64, Missings.Missing}, 1} &&
       typeof(locus_frame[:, :MaleMorgans]) != Array{Int64, 1}
      throw(ArgumentError(
        "All Morgans columns should only have numbers or missing values.\n \n"))
    end
  end
  if :CentiMorgans in locus_field
    if typeof(locus_frame[:, :CentiMorgans]) !=
         Array{Union{Float64, Missings.Missing}, 1} &&
       typeof(locus_frame[:, :CentiMorgans]) != Array{Float64, 1} &&
       typeof(locus_frame[:, :CentiMorgans]) !=
         Array{Union{Int64, Missings.Missing}, 1} &&
       typeof(locus_frame[:, :CentiMorgans]) != Array{Int64, 1}
      throw(ArgumentError(
        "All Morgans columns should only have numbers or missing values.\n \n"))
    end
  end
  #
  # Determine some dimensions for the regular locus structure.
  #
  rows = length(locus_frame[:, :Locus])
  columns = length(locus_field)
  locus_name = unique(locus_frame[:, :Locus])
  #
  # Strip spaces from the loci, allele, and chromosome names.
  # Also check for missing values.
  #
  for i = 1:rows
    row_string = string(i)
    if ismissing(locus_frame[i, :Locus])
      throw(ArgumentError(
        "Row $row_string of the input loci data has no locus name.\n \n"))
    end
    s = string(locus_frame[i, :Locus])
    if !isa(Meta.parse(s, raise=false), Number)
      locus_frame[i, :Locus] = strip(s)
      if locus_frame[i, :Locus] == ""
        throw(ArgumentError(
          "Row $row_string of the input loci data has no locus name.\n \n"))
      end
    end
    if ismissing(locus_frame[i, :Allele])
      throw(ArgumentError(
        "Row $row_string of the input loci data has no allele name.\n \n"))
    end
    s = string(locus_frame[i, :Allele])
    if !isa(Meta.parse(s, raise=false), Number)
      locus_frame[i, :Allele] = strip(s)
      if locus_frame[i, :Allele] == ""
        throw(ArgumentError(
          "Row $row_string of the input loci data has no allele name.\n \n"))
      end
    end
    if !ismissing(locus_frame[i, :Chromosome])
      s = string(locus_frame[i, :Chromosome])
      if !isa(Meta.parse(s, raise=false), Number)
        locus_frame[i, :Chromosome] = strip(s)
##        if locus_frame[i, :Chromosome] == ""
##          locus_frame[i, :Chromosome] = missing
##        end
      end
    end
  end
  #
  # Check for errors and omissions.
  #
  loci = 1
  for i = 2:rows
    if locus_frame[i, :Locus] != locus_frame[i - 1, :Locus]
      loci = loci + 1
    end
  end
  if loci != length(locus_name)
    throw(ArgumentError(
      "The alleles of each locus must be contiguous.\n \n"))
  end
  for i = 2:rows
    if locus_frame[i, :Locus] == locus_frame[i - 1, :Locus]
      locus_string = string(locus_frame[i, :Locus])
      if locus_frame[i, :Chromosome] != locus_frame[i - 1, :Chromosome]
        throw(ArgumentError(
          "The chromosome label of locus $locus_string is inconsistent.\n \n"))
      end
      if :Basepairs in locus_field
        if locus_frame[i, :Basepairs] != locus_frame[i - 1, :Basepairs]
          throw(ArgumentError(
            "The basepairs position of locus $locus_string is inconsistent.\n \n"))
        end
      end
      if :FemaleMorgans in locus_field
        if locus_frame[i, :FemaleMorgans] != locus_frame[i - 1, :FemaleMorgans]
          throw(ArgumentError(
            "The map position of locus $locus_string is inconsistent.\n \n"))
        end
      end
      if :MaleMorgans in locus_field
        if locus_frame[i, :MaleMorgans] != locus_frame[i - 1, :MaleMorgans]
          throw(ArgumentError(
            "The map position of locus $locus_string is inconsistent.\n \n"))
        end
      end
    end
  end
  #
  # Sort the locus frame by chromosome and location.
  # Warning: Chromosome names should be 1, 2, ... 22, or X
  # without further adornment.
  #
  if :Basepairs in locus_field
    sort!(locus_frame, (:Chromosome, :Basepairs))
  elseif :FemaleMorgans in locus_field
    sort!(locus_frame, (:Chromosome, :FemaleMorgans))
  end
  locus_field = names(locus_frame)
  person_field = names(person_frame)
  #
  # Get the loci names in the new sorted order.
  #
##  locus_name = unique(locus_frame[:, :Locus])
  j = 1
  locus_name[1] = locus_frame[1, :Locus]
  for i = 2:rows
    if locus_frame[i, :Locus] != locus_name[j]
      j = j + 1
      locus_name[j] = locus_frame[i, :Locus]
    end
  end
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
      if !ismissing(locus_frame[i, :Chromosome])
        chromosome[loc] = string(locus_frame[i, :Chromosome])
      end
      #
      # When only basepair distances are available,
      # equate 1e6 base pairs to a centiMorgan.
      #
      if :Basepairs in locus_field
        if !ismissing(locus_frame[i, :Basepairs])
          base_pairs[loc] = locus_frame[i, :Basepairs]
        end
        if !(:FemaleMorgans in locus_field)
          morgans[:, loc] = base_pairs[loc] / 1e8
        end
      end
      if :FemaleMorgans in locus_field
        if !ismissing(locus_frame[i, :FemaleMorgans])
          morgans[1, loc] = locus_frame[i, :FemaleMorgans]
        end
        if :MaleMorgans in locus_field
          if !ismissing(locus_frame[i, :MaleMorgans])
            morgans[2, loc] = locus_frame[i, :MaleMorgans]
          end
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
      if typeof(locus_frame[:, i]) ==
            Array{Union{Float64, Missings.Missing}, 1} ||
         typeof(locus_frame[:, i]) == Array{Float64, 1}
        s = Symbol(locus_field[i])
        if s != :FemaleMorgans && s != :MaleMorgans
          populations = 1
          push!(population_names, string(locus_field[i]))
          break
        end
      end
    end
  end
  #
  # Find the observed loci. These are the loci common to
  # the Pedigree frame and the Locus field in the Locus frame.
  #
  A = unique(locus_frame[:, :Locus])
  B = Vector{Symbol}(undef, size(A, 1))
  i = 0
  for loc in A
    i = i+1
    B[i] = Symbol(loc)
  end
  C = intersect(person_field, B)
  #
  # Find the indices of the matching field names in the Pedigree frame.
  #
  observed_indices = findall((in)(C), person_field)
  #
  # Find the indices of the observed loci in the locus structure.
  #
  observed_loci = length(C)
  observed_locus = zeros(Int, observed_loci)
  locus_field_in_person_frame = zeros(Int, observed_loci)
  i = 0
  for loc = 1:loci
    locus_symbol = Symbol(locus_name[loc])
    if locus_symbol in C
      i = i + 1
      observed_locus[i] = loc
      for j in observed_indices
        if isequal(person_field[j], locus_symbol)
          locus_field_in_person_frame[i] = j
          break
        end
      end
    end
  end
  #
  # Eliminate unobserved loci in the following one-dimensional arrays.
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
  # Collect the allele names and population frequencies for each observed locus.
  #
  allele_name = Array{Array{AbstractString, 1}}(undef, loci)
  frequency = Array{Array{Float64, 2}}(undef, loci)
  loc = 0
  n = 0
  for i = 1:rows
    if i == 1 || locus_frame[i, :Locus] != locus_frame[i - 1, :Locus]
      loc = loc + 1
      n = 1
      l = findfirst(isequal(loc), observed_locus)
      if l != nothing
        allele_name[l] = Array{AbstractString}(undef, alleles[l])
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
    l = findfirst(isequal(loc), observed_locus)
    if l != nothing
      allele_name[l][n] = string(locus_frame[i, :Allele])
      j = 0
      for pop in population_names
        j = j + 1
        if ismissing(locus_frame[i, Symbol(pop)])
          frequency[l][j, n] = 0.0
        else
          frequency[l][j, n] = locus_frame[i, Symbol(pop)]
        end
      end
    end
  end
  #
  # Check that allele frequencies are legal.
  #
  for loc = 1:loci
    loc_name = locus_name[loc]
    j = 0
    for pop in population_names
      j = j + 1
      if any(frequency[loc][j, :] .< 0.0)
        if populations == 1
          println("Warning: Negative allele frequency at locus $loc_name.")
        else
          println("Warning: " *
          "Negative allele frequency at locus $loc_name and population $pop.")
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
              " $loc_name do not sum to 1.0.")
          else
            println("Warning: Allele frequencies at locus",
              " $loc_name for population $pop do not sum to 1.0.")
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
    if :Basepairs in locus_field || :FemaleMorgans in locus_field
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
  person_field = names(person_frame)
  if :Person in person_field
    people = length(person_frame[:, :Person])
  else
    people = length(person_frame[:, :Individual])
  end
  if :Pedigree in person_field
    pedigrees = length(unique(person_frame[:, :Pedigree]))
  else
    pedigrees = people
  end
  lumped_frequency = zeros(pedigrees, populations, loci)
  #
  # As a default, equate the set of model loci to the set of
  # observed loci and the trait locus to 0.
  #
  model_loci = loci
  model_locus = collect(1:loci)
  trait = 0
##
  mutable_locus = 0 # locus subject to mutation
  donor_allele = 0 # allele subject to mutation
  recipient_allele = 0 # allele mutated to
  mutation_rate = zeros(2)
  #
  # Insert the various scalars and arrays in the locus data structure.
  #
##
  locus = Locus(loci, model_loci, trait, free_recombination, 
    mutable_locus, donor_allele, recipient_allele, mutation_rate,
    locus_name, chromosome, base_pairs, morgans, theta, xlinked, 
    alleles, allele_name, frequency, lumped_frequency, model_locus,
    locus_field_in_person_frame)
  return locus
end # function locus_information

"""
Extracts pedigree information from the pedigree frame.
Some fields are supplied later.
"""
function pedigree_information(person_frame::DataFrame)
  #
  # Strip any leading or trailing spaces from the field names.
  #
  for i in names(person_frame)
    rename!(person_frame, i => Symbol(strip(string(i))))
  end
  #
  # Check dataframe field names.
  #
  person_field = names(person_frame)
  if !(:Pedigree in person_field)
    if :Pedigrees in person_field
      rename!(person_frame, :Pedigrees => :Pedigree)
    elseif :pedigree in person_field
      rename!(person_frame, :pedigree => :Pedigree)
    elseif :pedigrees in person_field
      rename!(person_frame, :pedigrees => :Pedigree)
    elseif :Ped in person_field
      rename!(person_frame, :Ped => :Pedigree)
    elseif :Peds in person_field
      rename!(person_frame, :Peds => :Pedigree)
    elseif :ped in person_field
      rename!(person_frame, :ped => :Pedigree)
    elseif :peds in person_field
      rename!(person_frame, :peds => :Pedigree)
    end
    person_field = names(person_frame)
  end
  if !(:Person in person_field)
    if :Persons in person_field
      rename!(person_frame, :Persons => :Person)
    elseif :person in person_field
      rename!(person_frame, :person => :Person)
    elseif :persons in person_field
      rename!(person_frame, :persons => :Person)
    elseif :Per in person_field
      rename!(person_frame, :Per => :Person)
    elseif :Pers in person_field
      rename!(person_frame, :Pers => :Person)
    elseif :per in person_field
      rename!(person_frame, :per => :Person)
    elseif :pers in person_field
      rename!(person_frame, :pers => :Person)
    end
    person_field = names(person_frame)
  end
  if !(:Individual in person_field)
    if :Individuals in person_field
      rename!(person_frame, :Individuals => :Individual)
    elseif :individual in person_field
      rename!(person_frame, :individual => :Individual)
    elseif :individuals in person_field
      rename!(person_frame, :individuals => :Individual)
    elseif :Ind in person_field
      rename!(person_frame, :Ind => :Individual)
    elseif :Inds in person_field
      rename!(person_frame, :Inds => :Individual)
    elseif :ind in person_field
      rename!(person_frame, :ind => :Individual)
    elseif :inds in person_field
      rename!(person_frame, :inds => :Individual)
    end
    person_field = names(person_frame)
  end
  if !(:Father in person_field)
    if :Fathers in person_field
      rename!(person_frame, :Fathers => :Father)
    elseif :father in person_field
      rename!(person_frame, :father => :Father)
    elseif :fathers in person_field
      rename!(person_frame, :fathers => :Father)
    elseif :Dad in person_field
      rename!(person_frame, :Dad => :Father)
    elseif :Dads in person_field
      rename!(person_frame, :Dads => :Father)
    elseif :dad in person_field
      rename!(person_frame, :dad => :Father)
    elseif :dads in person_field
      rename!(person_frame, :dads => :Father)
    end
    person_field = names(person_frame)
  end
  if !(:Mother in person_field)
    if :Mothers in person_field
      rename!(person_frame, :Mothers => :Mother)
    elseif :mother in person_field
      rename!(person_frame, :mother => :Mother)
    elseif :mothers in person_field
      rename!(person_frame, :mothers => :Mother)
    elseif :Mom in person_field
      rename!(person_frame, :Mom => :Mother)
    elseif :Moms in person_field
      rename!(person_frame, :Moms => :Mother)
    elseif :mom in person_field
      rename!(person_frame, :mom => :Mother)
    elseif :moms in person_field
      rename!(person_frame, :moms => :Mother)
    end
    person_field = names(person_frame)
  end
  if !(:Sex in person_field)
    if :Sexes in person_field
      rename!(person_frame, :Sexes => :Sex)
    elseif :Sexs in person_field
      rename!(person_frame, :Sexs => :Sex)
    elseif :sex in person_field
      rename!(person_frame, :sex => :Sex)
    elseif :sexes in person_field
      rename!(person_frame, :sexes => :Sex)
    elseif :sexs in person_field
      rename!(person_frame, :sexs => :Sex)
    end
    person_field = names(person_frame)
  end
  if !(:Twin in person_field)
    if :Twins in person_field
      rename!(person_frame, :Twins => :Twin)
    elseif :twin in person_field
      rename!(person_frame, :twin => :Twin)
    elseif :twins in person_field
      rename!(person_frame, :twins => :Twin)
    end
    person_field = names(person_frame)
  end
  #
  # Count the number of individuals in the pedigree data.
  #
  if !(:Person in person_field || :Individual in person_field)
    throw(ArgumentError("The pedigree data file does not contain a field\n" *
      "labeled Person or Individual. One such field is required.\n \n"))
  elseif :Person in person_field && :Individual in person_field
    throw(ArgumentError("The pedigree data file contains a field\n" *
      "labeled Person and a field labeled Individual.\n" *
      "It is required to have only one such field.\n \n"))
  end
  if :Person in person_field
    people = length(person_frame[:, :Person])
  else
    people = length(person_frame[:, :Individual])
  end
  #
  # Initialize arrays.
  #
  if :Pedigree in person_field
    pedigrees = length(unique(person_frame[:, :Pedigree]))
  else
    pedigrees = people
  end
  pedigree_name = blanks(pedigrees)
  pedstart = zeros(Int, pedigrees)
  twin_finish = zeros(Int, pedigrees)
  pedfinish = zeros(Int, pedigrees)
  individuals = zeros(Int, pedigrees)
  founders = zeros(Int, pedigrees)
  females = zeros(Int, pedigrees)
  males = zeros(Int, pedigrees)
  twins = zeros(Int, pedigrees)
  families = zeros(Int, pedigrees)
  #
  # If pedigree names are included in the input data,
  # check whether each individual has a pedigree listed.
  #
  if :Pedigree in person_field
    for i = 1:people
      a = string(i)
      if ismissing(person_frame[i, :Pedigree])
        throw(ArgumentError(
          "Individual $a of the pedigree data has no pedigree name.\n \n"))
      end
      if typeof(person_frame[i, :Pedigree]) <: Number
        if person_frame[i, :Pedigree] == NaN
          throw(ArgumentError(
            "Individual $a of the pedigree data has no pedigree name.\n \n"))
        end
      elseif typeof(person_frame[i, :Pedigree]) <: AbstractString
        person_frame[i, :Pedigree] = strip(person_frame[i, :Pedigree])
        if person_frame[i, :Pedigree] == ""
          throw(ArgumentError(
            "Individual $a of the pedigree data has no pedigree name.\n \n"))
        end
      else
        throw(ArgumentError(
          "Individual $a of the pedigree data has no pedigree name.\n \n"))
      end
    end
    #
    # Check whether each pedigree occupies a contiguous block.
    #
    ped = 1
    for i = 2:people
      if person_frame[i, :Pedigree] != person_frame[i - 1, :Pedigree]
        ped = ped + 1
      end
    end
    if ped > pedigrees
#      throw(ArgumentError("Some pedigrees are not in a contiguous block." *
#      " Please fix the data files.\n \n"))
      sort!(person_frame, :Pedigree)
    end
    #
    # Find the name, start, and finish of each pedigree.
    #
    ped = 0
    for i = 1:people
      if i == 1 || person_frame[i, :Pedigree] != person_frame[i - 1, :Pedigree]
        ped = ped + 1
        pedstart[ped] = i
        if i > 1; twin_finish[ped - 1] = i - 1; end
        pedigree_name[ped] = string(person_frame[i, :Pedigree])
      end
    end
    twin_finish[pedigrees] = people
  else
    for i = 1:people
      pedstart[i] = i
      pedfinish[i] = i
      twin_finish[i] = i
      pedigree_name[i] = "$i"
    end
  end
  #
  # Insert the gathered information into the Pedigree structure.
  #
  pedigree = Pedigree(pedigrees, pedigree_name, pedstart, twin_finish,
    pedfinish, individuals, founders, females, males, twins, families,
    zeros(pedigrees, 2))
  return pedigree
end # function pedigree_information

"""
Extracts person information from the pedigree frame.
"""
function person_information(locus_frame::DataFrame, person_frame::DataFrame,
  phenotype_frame::DataFrame, locus::Locus, pedigree::Pedigree,
  keyword::Dict{AbstractString, Any})
  #
  # Determine the number of individuals.
  # (It has previously been checked that exactly one of :Person or :Individual
  # is in the pedigree dataframe.)
  #
  person_field = names(person_frame)
  if :Person in person_field
    people = length(person_frame[:, :Person])
  else
    people = length(person_frame[:, :Individual])
  end
  #
  # Create array of person names, which is not allowed to have missing values.
  #
  person_name = blanks(people)
  if :Person in person_field
    for i = 1:people
      a = string(i)
      if ismissing(person_frame[i, :Person])
        throw(ArgumentError(
          "Individual # $a of the pedigree data has no name.\n \n"))
      end
      if typeof(person_frame[i, :Person]) <: Number
        if person_frame[i, :Person] == NaN
          throw(ArgumentError(
            "Individual # $a of the pedigree data has no name.\n \n"))
        end
      elseif typeof(person_frame[i, :Person]) <: AbstractString
        person_frame[i, :Person] = strip(person_frame[i, :Person])
        if person_frame[i, :Person] == ""
          throw(ArgumentError(
            "Individual # $a of the pedigree data has no name.\n \n"))
        end
      else
        throw(ArgumentError(
          "Individual # $a of the pedigree data has no name.\n \n"))
      end
      person_name[i] = string(person_frame[i, :Person])
    end
  else
    for i = 1:people
      a = string(i)
      if ismissing(person_frame[i, :Individual])
        throw(ArgumentError(
          "Individual # $a of the pedigree data has no name.\n \n"))
      end
      if typeof(person_frame[i, :Individual]) <: Number
        if person_frame[i, :Individual] == NaN
          throw(ArgumentError(
            "Individual # $a of the pedigree data has no name.\n \n"))
        end
      elseif typeof(person_frame[i, :Individual]) <: AbstractString
        person_frame[i, :Individual] = strip(person_frame[i, :Individual])
        if person_frame[i, :Individual] == ""
          throw(ArgumentError(
            "Individual # $a of the pedigree data has no name.\n \n"))
        end
      else
        throw(ArgumentError(
          "Individual # $a of the pedigree data has no name.\n \n"))
      end
      person_name[i] = string(person_frame[i, :Individual])
    end
  end
  #
  # Initialize arrays by their default values.
  #
  pedigrees = pedigree.pedigrees
  allele_separator = keyword["allele_separator"]
  ordered_allele_separator = keyword["ordered_allele_separator"]
  mother = zeros(Int, people)
  father = zeros(Int, people)
  male = falses(people)
  #
  # Create array of parental names, which is allowed to have missing values.
  #
  parents_present = :Mother in person_field && :Father in person_field
  mother_string = blanks(people)
  father_string = blanks(people)
  if parents_present
    for i = 1:people
      if ismissing(person_frame[i, :Mother]); continue; end
      if typeof(person_frame[i, :Mother]) <: Number
        if person_frame[i, :Mother] == NaN ||
           person_frame[i, :Mother] == 0
          continue
        end
      elseif typeof(person_frame[i, :Mother]) <: AbstractString
        person_frame[i, :Mother] = strip(person_frame[i, :Mother])
        if person_frame[i, :Mother] == ""; continue; end
      else
        continue
      end
      if ismissing(person_frame[i, :Father]); continue; end
      if typeof(person_frame[i, :Father]) <: Number
        if person_frame[i, :Father] == NaN ||
           person_frame[i, :Father] == 0
          continue
        end
      elseif typeof(person_frame[i, :Father]) <: AbstractString
        person_frame[i, :Father] = strip(person_frame[i, :Father])
        if person_frame[i, :Father] == ""; continue; end
      else
        continue
      end
      mother_string[i] = string(person_frame[i, :Mother])
      father_string[i] = string(person_frame[i, :Father])
    end
  end
  #
  # Initialize the data structure containing twin status.
  # Check if the twin field is present in the input data.
  # If so, strip spaces from their values.
  #
  next_twin = zeros(Int, people)
  primary_twin = zeros(Int, people)
  if :Twin in person_field
    twins_present = true
    for i = 1:people
      if ismissing(person_frame[i, :Twin]); continue; end
      if typeof(person_frame[i, :Twin]) <: AbstractString
        person_frame[i, :Twin] = strip(person_frame[i, :Twin])
      end
    end
  #
  # Create pointers from each twin in a twin set to the next twin.
  # The last twin points to no-one.
  # In each twin set designate a primary twin.
  #
    twin_group = unique(person_frame[:, :Twin])
    for t in twin_group
      if ismissing(t)
        continue
      elseif typeof(t) <: Number && (t == 0 || t == NaN)
        continue
      else
        twin_list = findall((in)(t), person_frame[:, :Twin])
        primary_twin[twin_list[1]] = twin_list[1]
        for i = 2:length(twin_list)
          next_twin[twin_list[i - 1]] = twin_list[i]
          primary_twin[twin_list[i]] = twin_list[1]
        end
      end
    end
  else
    twins_present = false
  end
  #
  # Identify mothers and fathers.
  #
  pedigree_number = zeros(Int, people)
  for ped = 1:pedigrees
    for i = pedigree.pedstart[ped]:pedigree.twin_finish[ped]
      pedigree_number[i] = ped
      if parents_present && mother_string[i] != ""
        mother_found = false
        father_found = false
        for j = pedigree.pedstart[ped]:pedigree.twin_finish[ped]
          if j == i; continue; end
          if :Person in person_field
            name_j = string(person_frame[j, :Person])
          else
            name_j = string(person_frame[j, :Individual])
          end
          if mother_string[i] == name_j
            mother[i] = j
            mother_found = true
          end
          if father_string[i] == name_j
            father[i] = j
            father_found = true
          end
          if mother_found && father_found; break; end
        end
        pednum_string = string(ped)
        pednam_string = pedigree.name[ped]
        pernam_string = person_name[i]
        if !mother_found
          mother_name = mother_string[i]
          throw(ArgumentError(
            "In pedigree number $pednum_string with name $pednam_string,\n" *
            "person named $pernam_string has mother named $mother_name,\n" *
            "who does not appear to be a member of the pedigree.\n \n"))
        end
        if !father_found
          father_name = father_string[i]
          throw(ArgumentError(
            "In pedigree number $pednum_string,\n" *
            "person named $pernam_string has father named $father_name,\n" *
            "who does not appear to be a member of the pedigree.\n \n"))
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
      per_string = person_name[per]
      throw(ArgumentError(
        "Person $per_string is his or her own ancestor.\n \n"))
    end
    permute!(person_name, perm)
    #
    # Find the inverse permutation.
    #
    inverse_perm = collect(Union{Int64, Missings.Missing}, 1:people)
    for i = 1:people
      inverse_perm[perm[i]] = i
    end
    #
    # Redefine parental indicators.
    #
    parentlabel = copy(mother)
    for i = 1:people
      if parentlabel[i] != 0
        mother[inverse_perm[i]] = inverse_perm[parentlabel[i]]
      else
        mother[inverse_perm[i]] = 0
      end
    end
    parentlabel = copy(father)
    for i = 1:people
      if parentlabel[i] != 0
        father[inverse_perm[i]] = inverse_perm[parentlabel[i]]
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
    # Sort the person dataframe according to the inverse permutation.
    #
    person_frame[!, :Inverse_Perm] = inverse_perm
    sort!(person_frame, :Inverse_Perm)
  end
  #
  # Record who is male.
  # NOTE: if a sex value is missing, that individual is recorded as female.
  #
  male_symbols = keyword["male"]
  female_symbols = keyword["female"]
  if :Sex in person_field
    for i = 1:people
      if ismissing(person_frame[i, :Sex]); continue; end
      if typeof(person_frame[i, :Sex]) <: Number
        if person_frame[i, :Sex] == NaN; continue; end
        sex_i = person_frame[i, :Sex]
      elseif typeof(person_frame[i, :Sex]) <: AbstractString
        person_frame[i, :Sex] = strip(person_frame[i, :Sex])
        if person_frame[i, :Sex] == ""; continue; end
        sex_i = lowercase(person_frame[i, :Sex])
      else
        continue
      end
      if sex_i in male_symbols
        male[i] = true
      elseif sex_i in female_symbols
        male[i] = false
      else
        per_string = person_name[i]
        sex_i_string = string(sex_i)
        throw(ArgumentError(
          "Person $per_string has sex indicator $sex_i_string, " *
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
      elseif male[j] && male[k]
        per = person_name[i]
        throw(ArgumentError(
          "Person $per has both biological parents listed as male.\n \n"))
      elseif !male[j] && !male[k]
        per = person_name[i]
        throw(ArgumentError(
          "Person $per has both biological parents listed as female.\n \n"))
      end
    end
  end
  #
  # Search for quantitative variables in the person frame that are not
  # admixture proportions.
  #
  variables = 0
  for i = 1:length(person_field)
    if typeof(person_frame[:, i]) ==
          Array{Union{Float64, Missings.Missing}, 1} ||
       typeof(person_frame[:, i]) == Array{Float64, 1}
      if !(string(person_field[i]) in keyword["populations"])
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
  for i = 1:length(person_field)
    if typeof(person_frame[:, i]) ==
          Array{Union{Float64, Missings.Missing}, 1} ||
       typeof(person_frame[:, i]) == Array{Float64, 1}
      if !(string(person_field[i]) in keyword["populations"])
        variables = variables + 1
        variable_name[variables] = string(person_field[i])
        for per = 1:people
          if ismissing(person_frame[per, i])
            variable[per, variables] = NaN
          else
            variable[per, variables] = person_frame[per, i]
          end
        end
      end
    end
  end
  #
  # Record the disease status of each person.
  #
  disease_field = keyword["trait"]
  if disease_field != ""
    disease_field = Symbol(disease_field)
    if disease_field in person_field
      disease_status = blanks(people)
      status = person_frame[:, disease_field]
      for i = 1:people
        if ismissing(status[i])
          disease_status[i] = ""
        else
          disease_status[i] = strip(string(status[i]))
        end
      end
    else
      disease_field = keyword["trait"]
      throw(ArgumentError(
        "Specified disease status field ($disease_field) " *
        "is not in the pedigree file.\n \n"))
    end
  else
    disease_status = blanks(0)
  end
  #
  # Find the locus names and number of loci in the phenotype dataframe.
  #
  if size(phenotype_frame, 2) > 0
    locus_name = unique(phenotype_frame[:, :Locus])
    loci = length(locus_name)
    rows = length(phenotype_frame[:, :Locus])
    phenotypes = zeros(Int, loci)
    phenotype = Array{Array{AbstractString, 1}}(undef, loci)
    genotype_string = Array{Array{AbstractString, 1}}(undef, loci)
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
    # Reduce each phenotype to a set of integer pairs.
    # Each integer represents a numbered allele in the Locus frame.
    #
    genotype = Array{Array{Set{Tuple{Int, Int}}, 1}}(undef, loci)
    #
    # Loop over all loci in the Phenotype frame.
    #
    for loc = 1:loci
      cor_loc = correspond[loc]
      genotype[loc] = Array{Set{Tuple{Int, Int}}}(undef, phenotypes[loc])
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
          g_string = convert(AbstractString, g)
          (double, a1, a2) = fetch_genotype(locus.allele_name[cor_loc],
              g_string, allele_separator, ordered_allele_separator)
          if a1 == 0
            locname = locus_name[loc]
            throw(ArgumentError(
              "Invalid genotype $g_string at locus $locname.\n \n"))
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
  observed_genotype = Array{Set{Tuple{Int, Int}}}(undef, people, observed_loci)
  #
  # Loop over all observed loci.
  #
  errors = 0
  for loc = 1:observed_loci
    m = locus.alleles[loc]
    if size(phenotype_frame, 2) > 0
      j = inverse_correspond[loc] # locus in the Phenotype structure.
    end
    #
    # Loop over all people.
    #
    for i = 1:people
      match = false
      if size(phenotype_frame, 2) > 0
        #
        # First consider people with missing values.
        #
        observed_genotype[i, loc] = Set{Tuple{Int, Int}}()
        phen = person_frame[i, locus.locus_field_in_person_frame[loc]]
        if ismissing(phen) || phen == ""
          hemizygous = locus.xlinked[loc] && male[i]
          observed_genotype[i, loc] = full_genotype_set(m, hemizygous)
          match = true
        #
        # Next consider people with listed phenotypes.
        #
        else
          if j != 0
            phen = convert(AbstractString, phen)
            for n = 1:phenotypes[j]
              p = phenotype[j][n]
              if isequal(phen, p)
                #
                # Use a copy() of the genotype, so that when one is changed,
                # all the other copies are not identically changed.
                #
                observed_genotype[i, loc] = copy(genotype[j][n])
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
        g = person_frame[i, locus.locus_field_in_person_frame[loc]]
        if ismissing(g)
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
          (double, a1, a2) = fetch_genotype(locus.allele_name[loc],
              g, allele_separator, ordered_allele_separator)
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
  populations = max(1, length(keyword["populations"]))
  admixture = zeros(Float64, people, populations)
  if populations == 1
    fill!(admixture, 1.0)
  else
    for i = 1:people
      j = 1
      for pop in population_names
        if ismissing(person_frame[i, Symbol(pop)])
          admixture[i, j] = 0.0
        else
          admixture[i, j] = person_frame[i, Symbol(pop)]
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
Counts various features of a pedigree.
Mates is a set of parent pairs in a particular pedigree.
"""
function pedigree_counts!(pedigree::Pedigree, person::Person)

  for i = 1:pedigree.pedigrees
    founders = 0
    females = 0
    males = 0
    twins = 0
    co_twins = 0
    #
    # Mates is a set of parent pairs in a particular pedigree.
    #
    mates = Set{Tuple{Int, Int}}()
    for j = pedigree.pedstart[i]:pedigree.twin_finish[i]
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
    pedigree.pedfinish[i] = pedigree.twin_finish[i] - co_twins
  end
end # function pedigree_counts!

"""
Identify all nuclear families, spouses, and children.
""" 
function construct_nuclear_families(pedigree::Pedigree, person::Person)

  fams = sum(pedigree.families)
  nuclear_family = NuclearFamily(fams, zeros(Int, fams), zeros(Int, fams),
    zeros(Int, fams), empties(fams))
  mother = 0
  father = 0
  k = 0 # Index of the current nuclear family.
  #
  # Identify each nuclear family in a pedigree with a pair
  # of parents (mates). Collect the children and spouses of
  # each pedigree member in the process.
  #
  for i = 1:pedigree.pedigrees
    mates = Set{Tuple{Int, Int}}()
    for j = pedigree.pedstart[i]:pedigree.twin_finish[i]
      if person.father[j] != 0
        push!(mates, (person.mother[j], person.father[j]))
        push!(person.children[person.mother[j]], j)
        push!(person.children[person.father[j]], j)
        push!(person.spouse[person.mother[j]], person.father[j])
        push!(person.spouse[person.father[j]], person.mother[j])
      end
    end
    #
    # Loop over the nuclear families and assign parents and siblings.
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
Complete various data structures and perform a few checks.
Start by eliminating heterozygous genotypes in a hemizygote,
counting homozygotes, and filling in ethnic admixture fractions.
"""
function check_data_structures!(pedigree::Pedigree, person::Person,
  locus::Locus, nuclear_family::NuclearFamily,
  keyword::Dict{AbstractString, Any})
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
Return a permutation putting founders at the head of the pedigree,
co-twins at the tail of the pedigree, and everyone else
arranged so that parents precede their children.
Also check for directed cycles. For algorithmic details see:
Lawler E (1976) Combinatorial Optimization and Matroids.
Holt, Rinehart, and Winston.
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
    pedstart = pedigree.pedstart[ped]
    twin_finish = pedigree.twin_finish[ped]
    #
    # Put all co-twins at the tail of the permutation.
    #
    n = twin_finish
    for i = pedstart:twin_finish
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
    m = pedstart
    j = n
    for i = pedstart:twin_finish
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
    m = pedstart
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
Return a full set of ordered genotypes.
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
Recover the constituent alleles from a genotype.
"""
function fetch_genotype(allele::Vector{AbstractString},
         genotype::AbstractString, allele_separator::AbstractString,
         ordered_allele_separator::AbstractString)
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
  # and return them. Only one is relevent if the two alleles are the same
  # or if the genotype included the ordered allele separator.
  #
  if a1 == a2 ||
      findfirst((in)(collect(ordered_allele_separator)), genotype) != nothing
    return (false, a1, a2) # One genotype.
  else
    return (true, a1, a2) # Two genotypes.
  end
end # function fetch_genotype

"""
Check for a few pedigree and person inconsistencies.
"""
function preliminary_checks(pedigree::Pedigree, person::Person)

  errors = 0
  #
  # Check for repeated pedigree names.
  #
  (nonunique, repeat) = repeated_string(convert(Array{String}, pedigree.name))
  if nonunique
    errors = errors + 1
    println("Error: Pedigree name $repeat is used for multiple pedigrees.")
  end
  #
  # Check for repeated person names.
  #
  for ped = 1:pedigree.pedigrees
    pedstart = pedigree.pedstart[ped]
    twinfinish = pedigree.twin_finish[ped]
    (nonunique, repeat) =
      repeated_string(convert(Array{String},person.name[pedstart:twinfinish]))
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
      person_ped = person.pedigree[i]
      println("Error: Pedigree $person_ped has a blank person name.")
    end
    #
    # Check for missing parents, parents of the wrong sex, or parents
    # the same as children.
    #
    mother_present = person.mother[i] != 0
    father_present = person.father[i] != 0
    if mother_present != father_present
      errors = errors + 1
      person_name = person.name[i]
      println("Error: Person $person_name lacks one parent.")
    end
    if mother_present
      if person.male[person.mother[i]]
        errors = errors + 1
        mother_name = person.name[person.mother[i]]
        println("Error: The mother $mother_name is a male.")
      end
      if person.name[i] == person.name[person.mother[i]]
        errors = errors + 1
        person_name = person.name[i]
        println("Error: Person $person_name is her own mother.")
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
        person_name = person.name[i]
        println("Error: Person $person_name is his own father.")
      end
    end
    #
    # Check for twin inconsistencies.
    #
    if person.primary_twin[i] != 0
      j = person.next_twin[person.primary_twin[i]]
      if j == 0
        errors = errors + 1
        person_name = person.name[i]
        println("Error: The twin $person_name lacks identical siblings.")
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
Delete heterozygous genotypes ostensibly compatible
with X-linked male phenotypes.
"""
function check_hemizygous_genotypes!(person::Person, locus::Locus, loc::Int)

  if locus.xlinked[loc]
    for i = 1:person.people
      if person.male[i]
        for (gm, gp) in person.genotype[i, loc]
          if gm != gp
            delete!(person.genotype[i, loc], (gm, gp))
          end
        end
      end
    end
  end
end # function check_hemizygous_genotypes!

"""
Count the number of homozygous genotypes at a locus.
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
Compute the ethnic admixture fractions for all non-founders.
Warning: Parents must come before children in each pedigree
for averaging to work.
"""
function ethnic_admixture!(person::Person)
#
  populations = person.populations
  for i = 1:person.people
    if person.mother[i] != 0
      a = person.admixture[person.mother[i], :]
      b = person.admixture[person.father[i], :]
      person.admixture[i, :] = (a + b) / 2
    else
##       for j = 1:length(person.admixture[i, :])
      for j = 1:populations
        if ismissing(person.admixture[i, j])
          person.admixture[i, j] = 0.0
        end
      end
      s = sum(person.admixture[i, :])
      person.admixture[i, :] = person.admixture[i, :] / s
    end
  end
end # function ethnic_admixture!

"""
Estimate allele frequencies at locus loc by gene counting (EM algorithm).
People are treated as unrelated. A pseudocount is required
for each population allele pair. Warning: Estimation should precede
genotype elimination and allele consolidation.
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
Check that ancestral populations are in both the person and locus frames.
"""
function check_populations(locus_frame::DataFrame, person_frame::DataFrame,
                           keyword::Dict{AbstractString, Any})

  person_field = names(person_frame)
  for pop in keyword["populations"]
    if Symbol(pop) in person_field
      continue
    else
      throw(ArgumentError("Population $pop is not in the pedigree frame.\n \n"))
    end
  end
  if size(locus_frame, 2) > 0
    locus_field = names(locus_frame)
    for pop in keyword["populations"]
      if Symbol(pop) in locus_field
        continue
      else
        throw(ArgumentError("Population $pop is not in the locus frame.\n \n"))
      end
    end
  end
end # function check_populations
