################################################################################
# This set of functions houses OpenMendel's keyword's default values
# and processes the user-defined keywords that specify
# the external files to read and the analyses to perform.
################################################################################
#
# Required external packages.
#
# using Distributions # From package Distributions.
# using GLM           # From package GLM.
#
# Required OpenMendel packages and modules.
#
# using Search
# using SearchSetup
# using DataStructures

export process_keywords!, set_keyword_defaults!

"""
This function processes the keywords from the control file or command line.
"""
function process_keywords!(keyword::Dict{AbstractString, Any},
                           control_file::AbstractString, args)
  #
  # If the user did not specify any control file or keywords as command line
  # arguments, then simply write out a logo and quit.
  #
  if control_file == "" && length(args) == 0
    println("OpenMendel is a package for the statistical analysis of")
    println("genetic data.")
    println(" \nThe data input files can be listed either on the command")
    println("line or in a named control file.")
    println(" \nFor more information, please see the documentation at")
    println("http://www.genetics.ucla.edu/software/Mendel_current_doc.pdf")
    println(" \n")
    throw(ArgumentError(
      "No genetic data was specified, so no analysis was performed."))
    return nothing
  end
  #
  # Start the set of modified_keywords at the empty set.
  #
  set_of_modified_keywords = Set{AbstractString}()
  #
  # Modify the keywords using the commands in the Control file. This will
  # change the working directory to the directory containing the Control file.
  #
  if control_file != ""
    keyword["control_file"] = control_file
    push!(set_of_modified_keywords, "control_file")
    open_control_file!(keyword)
    read_control_file!(keyword, set_of_modified_keywords)
  end
  #
  # Modify the keywords using the function arguments.
  #
  if length(args) != 0
    read_args!(keyword, set_of_modified_keywords, args)
  end
  #
  # Make any necessary revisions of the keywords.
  #
  revise_keywords!(keyword, set_of_modified_keywords)
  #
  # Open the output file and list the user-modified keywords.
  #
  open_output_file!(keyword)
  list_modified_keywords(keyword, set_of_modified_keywords)
  #
  # Initialize the random number generator.
  #
  srand(keyword["seed"])
  return nothing
end # process_keywords!

"""
This function constructs a dictionary of the basic keywords
and assignes them their default values. Keywords specific
to an analysis option should be first declared at the marked position
in the source code for that analysis option.
"""
function set_keyword_defaults!(keyword::Dict{AbstractString, Any})
  #
  # Set the user-modifiable keywords using the format:
  # keyword["some_keyword_name"] = default_value
  #
  keyword["affected_designator"] = "1"
  keyword["allele_pseudo_count"] = 0.1
  keyword["allele_separator"] = "/\\" # first character is used in output
  keyword["analysis_option"] = ""
  keyword["complexity_threshold"] = 5e7
  keyword["control_file"] = ""
  keyword["eliminate_genotypes"] = false
  keyword["female"] = Set{Any}(["female", "f", 'f', "2", '2', 2])
  keyword["field_separator"] = ','
  keyword["genetic_map_function"] = "Haldane" # "Haldane" or "Kosambi"
  keyword["locus_file"] = ""
  keyword["lump_alleles"] = false
  keyword["male"] = Set{Any}(["male", "m", 'm', "1", '1', 1])
  keyword["new_pedigree_file"] = ""
  keyword["ordered_allele_separator"] = "|" # first character is used in output
  keyword["output_file"] = "Mendel_Output.txt"
  keyword["pedigree_file"] = ""
  keyword["phenotype_file"] = ""
  keyword["plink_input_basename"] = ""
  keyword["plink_output_basename"] = ""
  keyword["populations"] = Set{AbstractString}()
  keyword["product_mode"] = true
  keyword["seed"] = 1234
  keyword["snpdata_file"] = ""
  keyword["snpdefinition_file"] = ""
  keyword["trait"] = ""
  keyword["manhattan_plot_file"] = ""
  keyword["repetitions"] = ""
## For imputation:
##  keyword["gradient_provided"] = false
##  keyword["reference_haplotypes"] = 0
##  keyword["window_width"] = 200 # width of imputation window
  #
  # Add default optimization keywords from package Search.
  #
  keyword = optimization_keywords!(keyword)
  #
  # Non-user-modifiable keywords.
  #
  keyword["keywords_naming_general_sets"] =
    Set{AbstractString}(["female", "male"])
  keyword["keywords_naming_string_sets"] =
    Set{AbstractString}(["populations"])

  return keyword
end # function set_keyword_defaults!

"""
This function reads the keywords from the command line arguments.
"""
function read_args!(keyword::Dict{AbstractString, Any},
                    set_of_modified_keywords::Set{AbstractString}, args)

  for arg in args
    (keyword_symbol, value) = arg
    keyword_string = "$keyword_symbol"
    if !isascii(keyword_string)
      throw(ArgumentError(
        "The specified keyword $keyword_string is not recognized.\n \n"))
    end
    keyword_string = lowercase(keyword_string)
    keyword_string = replace(keyword_string, '-', '_')
    if haskey(keyword, keyword_string)
      keyword[keyword_string] = value
      push!(set_of_modified_keywords, keyword_string)
    else
      throw(ArgumentError(
        "The specified keyword $keyword_string is not a known keyword.\n \n"))
    end
    #
    # If a Control file is listed at any point, immediately
    # modify the keywords using the commands in the Control file.
    #
    if keyword_string == "control_file"
      open_control_file!(keyword)
      read_control_file!(keyword, set_of_modified_keywords)
    end
  end
  return
end # function read_args!

"""
This function revises the keywords based on the set of modified keywords.
"""
function revise_keywords!(keyword::Dict{AbstractString, Any},
                          set_of_modified_keywords::Set{AbstractString})

  const DOUBLE_QUOTE_CHAR :: Char = Char(34) # to avoid coloring bugs in editors
  #
  # Turn strings containing the names of distribution or link functions
  # into the named functions themselves.
  #
  if "distribution" in set_of_modified_keywords
    distribution_name = keyword["distribution"]
    if isa(distribution_name, AbstractString)
      keyword["distribution"] = eval(parse(distribution_name, raise = false))
    end
  end
  if "link" in set_of_modified_keywords
    link_name = keyword["link"]
    if isa(link_name, AbstractString)
      keyword["link"] = eval(parse(link_name, raise = false))
    end
  end
  #
  # Turn the field separator into a character.
  #
  if "field_separator" in set_of_modified_keywords
    keyword["field_separator"] = keyword["field_separator"][1]
  end
  #
  # If a Plink input files basename has been provided,
  # set the field separator to a blank and set the appropriate filenames.
  #
  if "plink_input_basename" in set_of_modified_keywords
    plink_basename = keyword["plink_input_basename"]
    keyword["snpdefinition_file"] = string(plink_basename, ".bim")
    push!(set_of_modified_keywords, "snpdefinition_file")
    keyword["snpdata_file"] = string(plink_basename, ".bed")
    push!(set_of_modified_keywords, "snpdata_file")
    keyword["pedigree_file"] = string(plink_basename, ".fam")
    push!(set_of_modified_keywords, "pedigree_file")
  end
  #
  # If the affected designator has been modified, make sure it is a string.
  # If it's not been modified and the pedigree is a Plink .fam file,
  # then change the affected designator to "2".
  #
  if "affected_designator" in set_of_modified_keywords
    keyword["affected_designator"] = string(keyword["affected_designator"])
  elseif contains(keyword["pedigree_file"], ".fam")
    keyword["affected_designator"] = "2"
    push!(set_of_modified_keywords, "affected_designator")
  end
  #
  # If the user has modified any of the keywords that take a set of strings
  # as its value, then form the set from the list of comma separated strings
  # that should have been assigned as the keyword's value.
  # Note that each string is stripped of leading and trailing spaces.
  #
  for keyword_string in
    intersect(keyword["keywords_naming_string_sets"], set_of_modified_keywords)
    set_of_strings =
      Set{AbstractString}(split(keyword[keyword_string], ','; keep = false))
    new_set = Set{AbstractString}()
    for element_string in set_of_strings
      element_string = strip(element_string, ' ')
      element_string = strip(element_string, DOUBLE_QUOTE_CHAR)
      element_string = strip(element_string, ' ')
      push!(new_set, element_string)
    end
    keyword[keyword_string] = new_set
  end
  #
  # If the user has modified any of the keywords that take a set of items
  # as its value, then form the set from the list of comma separated items
  # that should have been assigned as the keyword's value.
  # Note that the list of items can include strings, chars, and numbers.
  #
  for keyword_string in
    intersect(keyword["keywords_naming_general_sets"], set_of_modified_keywords)
    set_of_strings =
      Set{AbstractString}(split(keyword[keyword_string], ','; keep = false))
    new_set = Set{Any}()
    for element_string in set_of_strings
      element_string = strip(element_string, ' ')
      if isa(parse(element_string, raise = false), Number)
        element_item = parse(element_string)
      elseif isa(parse(element_string, raise = false), Char)
        element_item = parse(element_string)
      else
        element_item = strip(element_string, DOUBLE_QUOTE_CHAR)
        element_item = strip(element_item, ' ')
      end
      push!(new_set, element_item)
    end
    keyword[keyword_string] = new_set
  end
  #
  # A named analysis option used to be mandatory.
  #
  # if !("analysis_option" in set_of_modified_keywords)
  #   throw(ArgumentError("No analysis option was specified.\n \n"))
  # end
  #
  # A pedigree file is mandatory.
  #
  if !("pedigree_file" in set_of_modified_keywords)
    throw(ArgumentError("No pedigree file was named.\n \n"))
  end
  return
end # function revise_keywords!

"""
This function reads the keywords from the Control file.
"""
function read_control_file!(keyword::Dict{AbstractString, Any},
                    set_of_modified_keywords::Set{AbstractString})

  const COMMENT_CHARS :: AbstractString = "#!"
  const DOUBLE_QUOTE_CHAR :: Char = Char(34) # to avoid coloring bugs in editors

  control_unit = keyword["control_unit"]
  #
  # Parse the Control file line by line.
  #
  for keyline in eachline(control_unit)
    keyline = chomp(keyline) # remove trailing "newline" characters
    keyline = rstrip(keyline, '\r') # remove trailing "return" characters
    keyline = strip(keyline) # remove leading and trailing whitespace
    if length(keyline) == 0; continue; end
    #
    # Delete from comment characters to the end of the line.
    #
    comment_position = search(keyline, collect(COMMENT_CHARS))
    if comment_position == 1
      keyline = ""
    elseif comment_position > 1
      keyline = keyline[1:comment_position-1]
    end
    keyline = strip(keyline) # remove leading and trailing whitespace
    if length(keyline) == 0; continue; end
    #
    # Split keyline into keyword_string and value_string.
    #
    equalsign_position = search(keyline, '=')
    if equalsign_position == 0
      throw(ArgumentError(
        """The specified keyword line "$keyline" has no equals sign.\n \n"""))
    end
    keyword_string = strip(keyline[1:equalsign_position-1])
    value_string = strip(keyline[equalsign_position+1:end])
    if length(keyword_string) == 0
      throw(ArgumentError(
        """The specified keyword line "$keyline" contains no keyword.\n \n"""))
    end
    #
    # Test whether the keyword string is valid.
    #
    if !isascii(keyword_string)
      throw(ArgumentError(
        "The specified keyword $keyword_string is not recognized.\n \n"))
    end
    keyword_string = lowercase(keyword_string)
    keyword_string = replace(keyword_string, '-', '_')
    if !haskey(keyword, keyword_string)
      throw(ArgumentError(
        "The specified keyword $keyword_string is not a known keyword.\n \n"))
    end
    #
    # Check whether the keyword value is valid,
    # allow a blank value at the keyword "field_separator".
    #
    if length(value_string) == 0
      if keyword_string == "field_separator"
        value_string = " "
      else
        throw(ArgumentError(
        """The keyword line "$keyline" contains no key-value.\n \n"""))
      end
    end
    #
    # To prevent possible infinite loops, ignore the keyword control_file.
    #
    if keyword_string == "control_file"; continue; end
    #
    # Change the value_string into a boolean or number, if appropriate.
    # Strip value_string of any leading or trailing double-quote characters.
    # Record that value in the keyword dictionary. If the value parses as
    # a character, change the assigned value to that character.
    #
    lc_value_string = lowercase(value_string)
    if lc_value_string == "true" || lc_value_string == "t"
      keyword[keyword_string] = true
    elseif lc_value_string == "false" || lc_value_string == "f"
      keyword[keyword_string] = false
    elseif isa(parse(value_string, raise = false), Number)
      keyword[keyword_string] = parse(value_string)
    elseif keyword_string in keyword["keywords_naming_string_sets"]
      keyword[keyword_string] = strip(value_string, ' ')
    elseif keyword_string in keyword["keywords_naming_general_sets"]
      keyword[keyword_string] = strip(value_string, ' ')
    else
      keyword[keyword_string] = strip(value_string, DOUBLE_QUOTE_CHAR)
      if isa(parse(keyword[keyword_string], raise = false), Char)
        keyword[keyword_string] = parse(keyword[keyword_string])
      end
    end
    #
    # Record this keyword as being modified by the user.
    #
    push!(set_of_modified_keywords, keyword_string)
  #
  # Procede to next line of the Control file.
  #
  end
  close(control_unit)
  return
end #function read_control_file!

"""
This function opens the Control file for reading.
It also sets the current working directory to
the directory that contains the Control file.
"""
function open_control_file!(keyword::Dict{AbstractString, Any})

  control_file = keyword["control_file"]
  #
  # First, check that the named control file exists and is readable.
  #
  if isfile(control_file)
    #
    # Second, close any previously opened control file.
    #
    if haskey(keyword, "control_unit")
      control_unit = keyword["control_unit"]
      if isopen(control_unit); close(control_unit); end
    end
    #
    # Finally, open the new control file and move to its directory.
    #
    control_unit = open(control_file, "r")
    keyword["control_unit"] = control_unit
    control_dir = dirname(control_file)
    if (isdir(control_dir)); cd(control_dir); end
  else
    throw(ArgumentError(
      """The named control file "$control_file" is not readable or is not present.\n \n"""))
  end
  return
end # function open_control_file!

"""
This function opens the output file.
"""
function open_output_file!(keyword::Dict{AbstractString, Any})
  #
  #  Check that the output file is writable. If so, open it.
  #
  output_file = keyword["output_file"]
  if output_file != ""
    output_unit = open(output_file, "w")
    keyword["output_unit"] = output_unit
  end
##  if isfile(output_file) && iswritable(output_file)
##    output_unit = open(output_file, "w")
##    keyword["output_unit"] = output_unit
##  elseif isfile(output_file) && !iswriteable(output_file)
##    throw(ArgumentError(
##      """The named output file "$output_file" is not writable.\n \n"""))
##  elseif iswritable(pwd())
##    output_unit = open(output_file, "w")
##    keyword["output_unit"] = output_unit
##  else
##    throw(ArgumentError(
##      """The named output file "$output_file" is not writable.\n \n"""))
##  end
  return
end # function open_output_file!

"""
This function outputs the list of user-modified keywords.
"""
function list_modified_keywords(keyword::Dict{AbstractString, Any},
                    set_of_modified_keywords::Set{AbstractString})

  output_unit = keyword["output_unit"]
  #
  # First output the current working directory.
  #
  working_dir = pwd()
  println("""The current working directory is "$working_dir".\n""")
  #
  # Output list of all the user-modified keywords.
  #
  if length(set_of_modified_keywords) != 0
    println("Keywords modified by the user:\n")
    for keyword_string in sort(collect(set_of_modified_keywords))
      value = keyword[keyword_string]
      println("  $keyword_string = $value")
    end
    print(" \n")
  else
    println(" OpenMendel is using all the keywords at their default values.\n")
  end
  return
end # function list_modified_keywords

