"""
This module includes all the base functions of OpenMendel.
"""
module MendelBase
#
# Define the data structures used by OpenMendel.
#
include("data_structures.jl")
#
# Include functions to process the keywords that specify
# the data files to use and the analysis to perform.
#
include("keywords.jl")
#
# Include useful general utilities.
#
include("general_utilities.jl")
#
# Include useful genetic utilities.
#
include("genetic_utilities.jl")
#
# Include functions to read the genetic data from external files.
#
include("read_data.jl")
#
# Include functions to construct the likelihood model to analyze.
#
include("model_construction.jl")
#
# Include functions to prepare for and carryout a pedigree likelihood evaluation
# via the Elston-Stewart algorithm.
#
include("elston_stewart_preparation.jl")
include("elston_stewart_evaluation.jl")
#
# Include specialized optimization routines.
#
include("optimize.jl")

end # module MendelBase
