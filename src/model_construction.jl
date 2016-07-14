################################################################################
# This set of functions allows evaluation of Mendelian models.
# These functions can be modified by the user to create models
# appropriate to specific types of analyses.
################################################################################

# using DataStructures
# using GeneticUtilities

export penetrance, prior, transmission, initialize_optimization

"""
This function supplies a penetrance for individual i at a
multiple locus genotype starting and ending with particular
model loci. The current version is only appropriate for
multiplicative penetrances and fully penetrant loci.
"""
function penetrance(person::Person, locus::Locus,
  multi_genotype::Matrix{Int},
  par::Vector{Float64}, keyword::Dict{ASCIIString, Any},
  start::Int, finish::Int, i::Int)

  pen = 1.0
  for l = start:finish
    allele1 = multi_genotype[1, l]
    allele2 = multi_genotype[2, l]
    loc = locus.model_locus[l]
    p = 1.0 # for reduced penetrance let p depend on loc
    pen = p * pen
  end
  return pen
end # function penetrance

"""
This function supplies a prior probability for founder i.
His or her multiple locus genotype starts and ends with
particular model loci. The current version assumes
Hardy-Weinberg and linkage equilibrium.
"""
function prior(person::Person, locus::Locus,
  multi_genotype::Matrix{Int},
  par::Vector{Float64}, keyword::Dict{ASCIIString, Any},
  start::Int, finish::Int, i::Int)

  prior_prob = 1.0
  for l = start:finish
    loc = locus.model_locus[l]
    allele = multi_genotype[1, l]
    if keyword["analysis_option"] == "EstimateFrequencies"
      frequency = par[allele]
    else
      frequency = dot(vec(person.admixture[i, :]),
                      vec(locus.frequency[loc][:, allele]))
    end
    prior_prob = prior_prob * frequency
    if !locus.xlinked[loc] || !person.male[i]
      allele = multi_genotype[2, l]
      if keyword["analysis_option"] == "EstimateFrequencies"
        frequency = par[allele]
      else
        frequency = dot(vec(person.admixture[i, :]),
                        vec(locus.frequency[loc][:, allele]))
      end
      prior_prob = prior_prob * frequency
    end
  end
  return prior_prob
end # function prior

"""
This function supplies the probability that a parent i
with a particular multiple locus genotype transmits a
particular gamete to his or her child j. The multiple
locus genotype and gamete start and finish with particular
model loci. The current version assumes linked loci and
no chiasma interference.
"""
function transmission(person::Person, locus::Locus, gamete::Vector{Int},
  multi_genotype::Matrix{Int}, par::Vector{Float64}, 
  keyword::Dict{ASCIIString, Any}, start::Int, finish::Int, i::Int, j::Int)

  loc = locus.model_locus[start]
  xlinked = locus.xlinked[loc]
  #
  # For male to male inheritance at an x-linked locus,
  # set the transmission probability equal to 1.
  #
  if xlinked && person.male[i] && person.male[j]
    return 1.0
  end
  #
  # Consider the Gamete Competition model.
  #
  if keyword["analysis_option"] == "GameteCompetition"
    k = multi_genotype[1, start]
    l = multi_genotype[2, start]
    m = gamete[start]
    trans = 0.0
    if person.disease_status[j] == keyword["affected_designator"]
      if k == m; trans = trans + par[k] / (par[k] + par[l]); end
      if l == m; trans = trans + par[l] / (par[k] + par[l]); end
    else
      if k == m; trans = trans + 0.5; end
      if l == m; trans = trans + 0.5; end
    end
    return trans
  end
  #
  # Consider the linkage analysis models: Two Point and Location Scores.
  # First find the position information for Two Point analysis.
  #
  if keyword["analysis_option"] == "TwoPointLinkage"
    if length(par) == 2
      locus.theta[:, 1] = par
    elseif length(par) == 1
      locus.theta[:, 1] = par[1]
    end
  #
  # For the Location Scores model, find the map location of the trait locus.
  #
  elseif keyword["analysis_option"] == "LocationScores"
    trait_place = locus.model_locus[locus.trait]
    locus.morgans[1, trait_place] = par[1]
    locus.morgans[2, trait_place] = par[end]
    if locus.trait > 1
      left = locus.model_locus[locus.trait - 1]
      d = locus.morgans[1, trait_place] - locus.morgans[1, left]
      locus.theta[1, locus.trait - 1] = map_function(d, "Haldane")
      d = locus.morgans[2, trait_place] - locus.morgans[2, left]
      locus.theta[2, locus.trait - 1] = map_function(d, "Haldane")
    end
    if locus.trait < locus.model_loci
      right = locus.model_locus[locus.trait + 1]
      d = locus.morgans[1, right] - locus.morgans[1, trait_place]
      locus.theta[1, locus.trait] = map_function(d, "Haldane")
      d = locus.morgans[2, right] - locus.morgans[2, trait_place]
      locus.theta[2, locus.trait] = map_function(d, "Haldane")
    end
  end
  #
  # Store an indicator of the sex of the parent.
  #
  if person.male[i]
    i_sex = 2
  else
    i_sex = 1
  end
  #
  # Reduce the computations by considering only the heterozygous loci.
  # Use Trow's formula to express the recombination fraction
  # between two heterozygous loci in terms of the recombination
  # fractions between the adjacent loci that separate them.
  # Set the logical variable found to true when the first heterozygous
  # parental locus is found. Phase records the phase of the most
  # recent heterozygous parental locus.
  #
  trans = 1.0
  found = false
  phase = true
  r = 0.5
  for l = start:finish
    match1 = multi_genotype[1, l] == gamete[l]
    match2 = multi_genotype[2, l] == gamete[l]
    #
    # Check whether either the first or second parental gene at
    # the current locus matches the gamete gene at this locus.
    # If not, then return with 0 for the transmission probability.
    #
    if !match1 && !match2
      return 0.0
    end
    #
    # Check whether the current locus is heterozygous.
    #
    if match1 != match2
      if found
        if phase == match1
          trans = trans * (0.5 + r)
        else
          trans = trans * (0.5 - r)
        end
      else
        found = true
        if start == 1 || start == finish
          trans = 0.5
        else
          trans = 1.0
        end
      end
      phase = match1
      r = 0.5
    end
    if found && l < finish
      r = r * (1.0 - 2.0 * locus.theta[i_sex, l])
    end
  end
  if !found; trans = 1.0; end
  return trans
end # function transmission

"""
This function initializes a minimization problem.
"""
function initialize_optimization(locus::Locus, parameter::Parameter,
  keyword::Dict{ASCIIString, Any})
  #
  # Redefine the defaults as needed.
  # Set the values for Gamete Competition.
  #
  if keyword["analysis_option"] == "GameteCompetition"
    for i = 1:parameter.parameters
      parameter.par[i] = 1.0
      parameter.min[i] = 1e-5
      parameter.name[i] = "tau" * " $i"
      parameter.name[i] = rpad(parameter.name[i], 8, ' ')
    end
    #
    # Constrain the parameter for the most frequent allele.
    #
    loc = locus.model_locus[1]
    n = 0
    sum_freq = 0.0
    for allele = 1:locus.alleles[loc]
      if sum(locus.frequency[loc][:, allele]) > sum_freq
         n = allele
         sum_freq = sum(locus.frequency[loc][:, allele])
      end
    end
    parameter.constraint[1, n] = 1.0
    parameter.constraint_level[1] = 1.0
  #
  # Set the values to Estimate Frequencies.
  #
  elseif keyword["analysis_option"] == "EstimateFrequencies"
    loc = locus.model_locus[1]
    for i = 1:parameter.parameters
      parameter.par[i] = locus.frequency[loc][1, i]
      parameter.min[i] = 1e-5
      parameter.name[i] = "freq"*" $i"
      parameter.name[i] = rpad(parameter.name[i], 8, ' ')
    end
    #
    # Force the parameter frequencies to sum to 1.0.
    #
    parameter.constraint[1, :] = 1.0
    parameter.constraint_level[1] = 1.0
  #
  # Set the values for Two Point Linkage.
  #
  elseif keyword["analysis_option"] == "TwoPointLinkage"
    fill!(parameter.name, "theta")
    if parameter.travel == "grid"
      parameter.grid[:, 1] = [0.5, 0.4, 0.3, 0.2, 0.15, 0.1, 0.05, 0.01, 0.001]
    else
      fill!(parameter.par, 0.5 - 1e-5)
      fill!(parameter.min, 1e-5)
      fill!(parameter.max, 0.5)
      if parameter.parameters == 2
        parameter.name[1] = "xxtheta"
        parameter.name[2] = "xytheta"
      end
    end
  #
  # Initialize parameters for Location Scores.
  # The parameter vector specifies the map location of the trait locus.
  #
  elseif keyword["analysis_option"] == "LocationScores"
    if parameter.parameters == 1
      parameter.name[1] = "theta"
    else
      parameter.name[1] = "xxtheta"
      parameter.name[2] = "xytheta"
    end
    small = 1e-4
    if locus.trait > 1
      left = locus.model_locus[locus.trait - 1]
      parameter.min[1] = locus.morgans[1, left] + small
      parameter.min[end] = locus.morgans[end, left] + small
    end
    if locus.trait < locus.model_loci
      right = locus.model_locus[locus.trait + 1]
      parameter.max[1] = locus.morgans[1, right] - small
      parameter.max[end] = locus.morgans[end, right] - small
    end
    d = keyword["flanking_distance"]
    if locus.trait == 1
      right = locus.model_locus[2]
      parameter.min[1] = locus.morgans[1, right] - d
      parameter.min[end] = locus.morgans[end, right] - d
    elseif locus.trait == locus.model_loci
      left = locus.model_locus[locus.trait - 1]
      parameter.max[1] = locus.morgans[1, left] + d
      parameter.max[end] = locus.morgans[end, left] + d
    end
    if keyword["travel"] == "search"
      parameter.par = 0.5*(parameter.min + parameter.max)
      if keyword["gender_neutral"]
        a = (parameter.max[2] - parameter.min[2])
        b = (parameter.max[1] - parameter.min[1])
        c = a / b
        d = parameter.min[2] - c * parameter.min[1]
        parameter.constraint[1, 1] = c
        parameter.constraint[1, 2] = -1.0
        parameter.constraint_level[1] = c * parameter.min[1] - parameter.min[2]
      end
    else
      if parameter.points == 1
        parameter.grid[1, :] = 0.5*(parameter.min + parameter.max)
      else
        for j = 1:parameter.points
          a = (j - 1.0) / (parameter.points - 1.0)
          parameter.grid[j, :] = a * parameter.max + (1.0 - a) * parameter.min
        end
        parameter.min[1:end] = -Inf
        parameter.max[1:end] = Inf
      end
    end
  end
  return parameter
end # function initialize_optimization

