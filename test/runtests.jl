module PkgTest

using MendelBase
using Base.Test

# write your own tests here

include("genetic_utilities_test.jl")
include("general_utilities_test.jl")
include("keywords_test.jl")
include("read_data_test.jl")
end # PkgTest module


# julia -e 'Pkg.test("MendelBase",coverage=true)'
# @show get_summary(process_file("src/read_data.jl"))
# @show get_summary(process_file("src/genetic_utilities.jl"))
# @show get_summary(process_file("src/general_utilities.jl"))
# @show get_summary(process_file("src/keywords.jl"))
# @show get_summary(process_file("src/data_structures.jl"))
# @show get_summary(process_file("src/elston_stewart_evaluation.jl"))
# @show get_summary(process_file("src/elston_stewart_preparation.jl"))