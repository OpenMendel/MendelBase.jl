using MendelBase, Base.Test

info("Unit tests for keywords.jl")

srand(123)

# The keywords function only exports 2 functions, and hence there will only be 
# two sets of unit tests. Private functions will be tested through these large
# testsets. 

@testset "set_default_keywords" begin
    keyword = set_keyword_defaults!(Dict{AbstractString, Any}())

    @test keyword["affected_designator"] == "1"
    @test keyword["allele_pseudo_count"] == 0.1
    @test keyword["allele_separator"] == "/\\" # first character is used in output
    @test keyword["analysis_option"] == ""
    @test keyword["complexity_threshold"] == 5e7
    @test keyword["control_file"] == ""
    @test keyword["eliminate_genotypes"] == false
    @test keyword["female"] == Set{Any}(["female", "f", 'f', "2", '2', 2])
    @test keyword["field_separator"] == ','
    @test keyword["genetic_map_function"] == "Haldane" # "Haldane" or "Kosambi"
    @test keyword["locus_file"] == ""
    @test keyword["lump_alleles"] == false
    @test keyword["male"] == Set{Any}(["male", "m", 'm', "1", '1', 1])
    @test keyword["new_pedigree_file"] == ""
    @test keyword["ordered_allele_separator"] == "|" # first character is used in output
    @test keyword["output_file"] == "Mendel_Output.txt"
    @test keyword["pedigree_file"] == ""
    @test keyword["phenotype_file"] == ""
    @test keyword["plink_input_basename"] == ""
    @test keyword["plink_output_basename"] == ""
    @test keyword["populations"] == Set{AbstractString}()
    @test keyword["product_mode"] == true
    @test keyword["seed"] == 1234
    @test keyword["snpdata_file"] == ""
    @test keyword["snpdefinition_file"] == ""
    @test keyword["trait"] == ""

    # optimization keywords from package Search.
    @test keyword["cases"] == 0 # number of cases in a least squares problem
    @test keyword["constraints"] == 0 # number of affine constraints
    @test keyword["goal"] == "maximize" # "minimize" or "maximize"
    @test keyword["output_unit"] == STDOUT # for output of Search iterations
    @test keyword["parameters"] == 1 # number of parameters
    @test keyword["points"] == 0 # number of points in a grid search
    @test keyword["standard_errors"] == false # true for parameter standard errors
    @test keyword["title"] == "" # problem title
    @test keyword["travel"] == "search" # "search" or "grid"

    # Non-user-modifiable keywords.
    @test keyword["keywords_naming_general_sets"] ==
        Set{AbstractString}(["female", "male"])
    @test keyword["keywords_naming_string_sets"] ==
        Set{AbstractString}(["populations"])
end

@testset "process_keywords" begin
    keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
    control_file = "control.txt" #sample control file from MendelBase's documentation
    process_keywords!(keyword, control_file, [])

    # test modify keywords
    @test keyword["locus_file"] == "gamete competition LocusFrame.txt"
    @test keyword["pedigree_file"] == "gamete competition PedigreeFrame.txt"
    @test keyword["output_file"] == "gamete competition Output.txt"
    @test keyword["trait"] == "ACE"
    @test keyword["affected_designator"] == "1"
    @test keyword["standard_errors"] == true 

    # test command line keyword modification
    # note: args is a list of whose element is of the form (keyword_symbol, value)
    args = [("proDUCT-mode", false), ("AlLeLE_pSeUdO_coUNt", 1.234)]
    process_keywords!(keyword, control_file, args)
    @test keyword["product_mode"] == false
    @test keyword["allele_pseudo_count"] == 1.234
end


@testset "are errors of process_keywords begin thrown" begin
    dict = Dict{AbstractString,Any}()
    control_file = "control.txt" #sample control file from MendelBase's documentation
    @test_throws(ArgumentError, process_keywords!(dict, "", [])) #empty set & args error 

    bad_arg = [("我們在懷念的演唱會", "禮貌的吻別～～～～～")] #non ascii error
    @test_throws(ArgumentError, process_keywords!(dict, "control_file", bad_arg))
    bad_arg2 = [("gg rekt ezpz", "ff15 reports for all")] # !haskey error
    @test_throws(ArgumentError, process_keywords!(dict, "control_file", bad_arg2))
    bad_arg3 = [("female", Set{Any}(["male", "m", 'm', "1", '1', 1]))] #changing un-modifiable keywords
    @test_throws(ArgumentError, process_keywords!(dict, "control_file", bad_arg3))
end


