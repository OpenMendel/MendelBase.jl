################################################################################
# The following data structures are the basis of OpenMendel.
# Complex data types allow related variables to be accessed as a whole.
# The plural of a variable indicates the number of instances of the variable.
################################################################################

export Person, Pedigree, NuclearFamily, Locus, SnpDataStruct

mutable struct Person
  people :: Int
  populations :: Int
  name :: Vector{AbstractString}
  pedigree :: Vector{Int}
  mother :: Vector{Int}
  father :: Vector{Int}
  male :: Vector{Bool}
  next_twin :: Vector{Int} # pointer to co-twins
  primary_twin :: Vector{Int} # pointer to primary twin
  admixture :: Matrix{Float64} # person.admixture[person, population]
  children :: Vector{BitSet}
  spouse :: Vector{BitSet}
  genotype :: Matrix{Set{Tuple{Int, Int}}} # person.genotype[person, locus]
  homozygotes :: Matrix{Int} # person.homozygotes[person, locus]
  disease_status :: Vector{AbstractString}
  variable :: Matrix{Float64} # person.variable[person, variable]
  variable_name :: Vector{AbstractString}
end # Person

mutable struct Pedigree
  pedigrees :: Int
  name :: Vector{AbstractString}
  start :: Vector{Int}
  twin_finish :: Vector{Int} # end of pedigree including co-twins
  finish :: Vector{Int} # end of pedigree excluding co-twins
  individuals :: Vector{Int} # count including co-twins
  founders :: Vector{Int}
  females :: Vector{Int}
  males :: Vector{Int}
  twins :: Vector{Int}
  families :: Vector{Int} # number of nuclear families per pedigree
  loglikelihood :: Matrix{Float64} # pedigree.loglikelihood[pedigree, 1:2]
end # Pedigree

mutable struct NuclearFamily
  families :: Int # number of nuclear families across all pedigrees
  pedigree :: Vector{Int}
  mother:: Vector{Int}
  father :: Vector{Int}
  sib :: Vector{BitSet}
end # NuclearFamily

mutable struct Locus
  loci :: Int
  model_loci :: Int
  trait :: Int # position of the trait locus among the model loci
  free_recombination :: Bool # true when model loci unlinked
  name :: Vector{AbstractString}
  chromosome :: Vector{AbstractString}
  base_pairs :: Vector{Int}
  morgans :: Matrix{Float64} # female and male genetic map distance
  theta :: Matrix{Float64} # model loci recombination fractions
  xlinked :: Vector{Bool}
  alleles :: Vector{Int}
  allele_name :: Vector{Vector{AbstractString}} # locus.allele_name[loc][allele]
  frequency :: Vector{Matrix{Float64}} # locus.frequency[loc][population,allele]
  lumped_frequency :: Array{Float64, 3} # indices: pedigree, population, locus
  model_locus :: Vector{Int} # indices of loci currently modeled
  locus_field_in_pedigree_frame :: Vector{Int}
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

