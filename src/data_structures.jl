################################################################################
# The following data structures are the basis of OpenMendel.
# Complex data types allow related variables to be accessed as a whole.
# The plural of a variable indicates the number of instances of the variable.
# Consider removing those marked with <remove>.
################################################################################

export Person, Pedigree, NuclearFamily, Locus, SnpDataStruct

mutable struct Person
  people :: Int # number of people in sample
  populations :: Int # number of populations in sample
  name :: Vector{AbstractString} # <remove>
  pedigree :: Vector{Int} # <remove>
  mother :: Vector{Int} # index of mother
  father :: Vector{Int} # index of father
  male :: Vector{Bool} # true for male; false for female
  next_twin :: Vector{Int} # pointer to next twin
  primary_twin :: Vector{Int} # pointer to primary twin
  admixture :: Matrix{Float64} # person.admixture[person, population]
  children :: Vector{BitSet} # set of children
  spouse :: Vector{BitSet} # set of spouses
  genotype :: Matrix{Set{Tuple{Int, Int}}} # person.genotype[person, locus]
  homozygotes :: Matrix{Int} # person.homozygotes[person, locus]
  disease_status :: Vector{AbstractString} # <remove>
  variable :: Matrix{Float64} # person.variable[person, variable] <remove>
  variable_name :: Vector{AbstractString} # <remove>
end # Person

mutable struct Pedigree
  pedigrees :: Int # number of pedigrees in sample
  name :: Vector{AbstractString} # <remove>
  pedstart :: Vector{Int} # start of pedigree
  twin_finish :: Vector{Int} # end of pedigree including co-twins
  pedfinish :: Vector{Int} # end of pedigree excluding co-twins
  individuals :: Vector{Int} # count per pedigree including co-twins
  founders :: Vector{Int} # count per pedigree
  females :: Vector{Int} # count per pedigree
  males :: Vector{Int} # count per pedigree
  twins :: Vector{Int} # count per pedigree
  families :: Vector{Int} # count of nuclear families per pedigree
  loglikelihood :: Matrix{Float64} # pedigree.loglikelihood[pedigree, 1:2]
end # Pedigree

mutable struct NuclearFamily
  families :: Int # number of nuclear families across all pedigrees
  pedigree :: Vector{Int} # pedigree of nuclear family
  mother:: Vector{Int} # mother of nuclear family
  father :: Vector{Int} # father of nuclear family
  sib :: Vector{BitSet} # sib set of nuclear family
end # NuclearFamily

mutable struct Locus
  loci :: Int # number of loci in locus file
  model_loci :: Int # number of loci in current model
  trait :: Int # position of the trait locus among the model loci
  free_recombination :: Bool # true when model loci unlinked
##
  mutable_locus :: Int # locus subject to mutation
  donor_allele :: Int # allele subject to mutation
  recipient_allele :: Int # allele mutated to
  mutation_rate :: Vector{Float64} # female and male mutation rates
  name :: Vector{AbstractString} # <remove>
  chromosome :: Vector{AbstractString} # <remove>
  base_pairs :: Vector{Int} # <remove>
  morgans :: Matrix{Float64} # female and male genetic map distance <remove ?>
  theta :: Matrix{Float64} # model loci recombination fractions
  xlinked :: Vector{Bool} # true for xlinked; false for autosomal
  alleles :: Vector{Int} # number of alleles per locus
  allele_name :: Vector{Vector{AbstractString}} # locus.allele_name[loc][allele]
  frequency :: Vector{Matrix{Float64}} # locus.frequency[loc][population,allele]
  lumped_frequency :: Array{Float64, 3} # indices: pedigree, population, locus
  model_locus :: Vector{Int} # indices of loci currently modeled
  locus_field_in_person_frame :: Vector{Int}
end # Locus

#
# Old-style structure for SNP and person information.
# (Should be better coordinated with SnpArrays.)
#
mutable struct SnpDataStruct
  people::Int                           # number of rows (individuals)
  snps::Int                             # number of columns (snps)
  personid::Vector{AbstractString}      # names of individuals
  snpid::Vector{AbstractString}         # SNP ids
  chromosome::Vector{AbstractString}    # SNP chromosome
  genetic_distance::Vector{Float64}     # genetic distance
  basepairs::Vector{Int}                # SNP base pair position
  allele1::Vector{AbstractString}       # A1 code
  allele2::Vector{AbstractString}       # A2 code
  maf::Vector{Float64}                  # minor allele frequencies
  minor_allele::BitVector               # bitvector designating the minor allele
  snpmatrix::AbstractSnpArray           # matrix of genotypes or haplotypes
  missings_per_person::Vector{Float64}  # number of missing genotypes per person
  missings_per_snp::Vector{Float64}     # number of missing genotypes per snp
end # end SnpDataStruct

