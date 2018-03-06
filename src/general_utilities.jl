################################################################################
# This set of functions provides various general utilities.
################################################################################

export empties, blanks, repeated_string, select_set_element, random_category
export normalize!, print_sample_stats, sample_mean_std, sample_stats
export simes_fdr, regress, glm_score_test

"""
Creates an array of blank strings.
"""
function blanks(n::Int)

  blank = Array{AbstractString}(n)
  for i = 1:n
    blank[i] = ""
  end
  return blank
end # function blanks

"""
Creates an array of empty integer sets.
"""
function empties(n::Int)

  empty = Array{IntSet}(n)
  for i = 1:n
    empty[i] = IntSet()
  end
  return empty
end # function empties

"""
Identifies a repeated string in an array of strings.
"""
function repeated_string(name::Vector{<:AbstractString})

  sorted = sort(name)
  for i = 2:length(name)
    if sorted[i - 1] == sorted[i]
      return (true, sorted[i])
    end
  end
  return (false, "")
end # function repeated_string

"""
Selects element number k of the integer set S.
"""
function select_set_element(set_s::IntSet, k::Int)

  j = 1
  for s in set_s
    if j == k
      return s
    else
      j = j + 1
    end
  end
  return nothing
end # function select_set_element

"""
Returns a random category (numbered 1, 2, and so forth)
according to a vector of category frequencies.
"""
function random_category(frequency::Vector{Float64})

  u = rand(1)
  for i = 1:length(frequency)
    if u[1] <= frequency[i]
      return i
    else
      u[1] = u[1] - frequency[i]
    end
  end
  return length(frequency)
end # function random_category

"""
Normalizes the vector x to have mean 0.0 and variance 1.0.
Missing values are ignored.
"""
function normalize!(x::Vector{Float64})

  (p, avg, ss) = zeros(3)
  for i = 1:length(x)
    if !isnan(x[i])
      p = p + 1.0
      d = x[i] - avg
      avg = avg + d / p
      ss = ss + (p - 1.0) * d^2 / p
    end
  end
  stddev = sqrt(ss / p)
  for i = 1:length(x)
    if !isnan(x[i])
      x[i] = (x[i] - avg) / stddev
    end
  end
  return x
end # function normalize!

"""
Computes a sample mean and standard deviation.
Missing values are ignored.
"""
function sample_mean_std(x::Vector{Float64})

  (p,avg,ss) = zeros(3)
  for i = 1:length(x)
    if !isnan(x[i])
      p = p + 1.0
      d = x[i] - avg
      avg = avg + d/p
      ss = ss + (p - 1.0)*(d^2)/p
    end
  end
  if p > 0.0
    stddev = sqrt(ss/p)
  else
    stddev = 0.0
  end
  return (avg, stddev)
end # function sample_mean_std

"""
Computes sample statistics for the data vector x.
Missing values equal NaN. The output consists of
the number of values present, the number of values missing,
the minimum, the 0.25 quantile, the median, the 0.75 quantile,
the maximum, the mean, the standard deviation, the skewness,
and the excess kurtosis.
"""
function sample_stats(x::Vector{Float64})

  y = copy(x)
  (p, n) = (0, length(x))
  for i = 1:n
    if !isnan(x[i])
      p = p + 1
      y[p] = x[i]
    end
  end
  return p, n - p, minimum(y[1:p]), quantile(y[1:p], 0.25), median(y[1:p]),
    quantile(y[1:p], 0.75), maximum(y[1:p]), mean(y[1:p]), std(y[1:p]),
    skewness(y[1:p]), kurtosis(y[1:p])
end # function sample_stats

"""
Prints the sample statistics output by the function sample_stats.
"""
function print_sample_stats(s, io::IO = STDOUT)

  println(io, "Summary Stats:")
  println(io, "Values Present:     ", s[1])
  println(io, "Values Absent:      ", s[2])
  println(io, "Minimum:            ", s[3])
  println(io, "1st Quartile:       ", s[4])
  println(io, "Median:             ", s[5])
  println(io, "3rd Quartile:       ", s[6])
  println(io, "Maximum:            ", s[7])
  println(io, "Mean:               ", s[8])
  println(io, "Standard Deviation: ", s[9])
  println(io, "Skewness:           ", s[10])
  println(io, "Excess Kurtosis:    ", s[11])
end # function print_sample_stats

"""
Performs the Simes false discovery rate (FDR) procedure discussed 
by Benjamini and Hochberg. All p-values at or below the threshold
are declared significant for the given FDR rate and number of tests.
"""
function simes_fdr(pvalue::Vector{Float64}, fdr::Vector{Float64}, tests::Int)

  n = length(fdr)
  threshold = zeros(n)
  number_passing = zeros(Int, n)
  perm = sortperm(pvalue)
  k = 1 # previous value of i
  for j = 1:n
    i = 0
    for i = k:length(pvalue)
      if tests * pvalue[perm[i]] > i * fdr[j]; break; end
    end
    if i == 1; continue; end
    k = i - 1
    threshold[j] = pvalue[perm[i - 1]]
    number_passing[j] = i - 1
  end
  return (number_passing, threshold)
end # function simes_fdr

"""
Performs either linear, logistic, or Poisson regression with a canonical 
link function. X is the design matrix, and y is the response vector.
The parameter estimates and loglikelihood are returned.
"""
function regress(X::Matrix{Float64}, y::Vector{Float64}, model::AbstractString)

  if model != "linear" && model != "logistic" && model != "Poisson"
    throw(ArgumentError(
      "The only model choices are linear, logistic, and Poisson.\n \n"))
  end
  #
  # Create the score vector, information matrix, estimate, a work
  # vector z, and the loglikelihood.
  #
  (n, p) = size(X)
  @assert n == length(y)
  score = zeros(p)
  information = zeros(p, p)
  estimate = zeros(p)
  z = zeros(n)
  #
  # Handle linear regression separately.
  #
  if model == "linear"
    BLAS.gemv!('N', 1.0, X, estimate, 0.0, z) # z = X * estimate
    BLAS.axpy!(-1.0, y, z) # z = z - y
    score = BLAS.gemv('T', -1.0, X, z) # score = - X' * (z - y)
    information = BLAS.gemm('T', 'N', X, X) # information = X' * X
    estimate = information \ score
    BLAS.gemv!('N', 1.0, X, estimate, 0.0, z) # z = X * estimate
    obj = - 0.5 * n * log(sum(abs2, y - z) / n) - 0.5 * n
    return (estimate, obj)
  end
  #
  # Prepare for logistic and Poisson regression by estimating the 
  # intercept.
  #
  if model == "logistic"
    estimate[1] = log(mean(y) / (1.0 - mean(y)))
  elseif model == "Poisson"
    estimate[1] = log(mean(y))
  else
    throw(ArgumentError(
      "The only model choices are linear, logistic, and Poisson.\n \n"))
  end
  #
  # Initialize the loglikelihood and the convergence criterion.
  #
  v = zeros(p)
  obj = 0.0
  old_obj = 0.0
  epsilon = 1e-6
  # 
  #  Estimate parameters by the scoring algorithm.
  #
  for iteration = 1:10
    #
    # Initialize the score and information.
    #
    fill!(score, 0.0)
    fill!(information, 0.0)
    #
    # Compute the score, information, and loglikelihood (obj).
    #
    BLAS.gemv!('N', 1.0, X, estimate, 0.0, z) # z = X * estimate
    clamp!(z, -20.0, 20.0) 
    if model == "logistic"
      z = exp.(-z)
      z = 1.0 ./ (1.0 + z)
      w = z .* (1.0 .- z)
      BLAS.axpy!(-1.0, y, z) # z = z - y
      score = BLAS.gemv('T', -1.0, X, z) # score = - X' * (z - y)
      w = sqrt.(w)
      scale!(w, X) # diag(w) * X
      information = BLAS.gemm('T', 'N', X, X) # information = X' * W * X
      w = 1.0 ./ w
      scale!(w, X)
    elseif model == "Poisson"
      z = exp.(z)
      w = copy(z)
      BLAS.axpy!(-1.0, y, z) # z = z - y
      score = BLAS.gemv('T', -1.0, X, z) # score = - X' * (z - y)
      w = sqrt.(w)
      scale!(w, X) # diag(w) * X
      information = BLAS.gemm('T', 'N', X, X) # information = X' * W * X
      w = 1.0 ./ w
      scale!(w, X)
    end
    #
    # Compute the scoring increment.
    #
    increment = information \ score
    #
    # Step halve to produce an increase in the loglikelihood.
    #
    steps = -1
    for step_halve = 0:3
      steps = steps + 1
      obj = 0.0
      estimate = estimate + increment
      BLAS.gemv!('N', 1.0, X, estimate, 0.0, z) # z = X * estimate
      clamp!(z, -20.0, 20.0)
      #
      # Compute the loglikelihood under the appropriate model.
      #
      if model == "logistic"
        z = exp.(-z)
        z = 1.0 ./ (1.0 + z)
        for i = 1:n
          if y[i] > 0.0
            obj = obj + log(z[i])
          else
            obj = obj + log(1.0 - z[i])
          end
        end
      elseif model == "Poisson"
        for i = 1:n
          q = exp(z[i])
          obj = obj + y[i] * z[i] - q
        end  
      end
      #
      # Check for an increase in the loglikelihood.
      #
      if old_obj < obj
        break
      else
        estimate = estimate - increment
        increment = 0.5 * increment
      end
    end
    #
    # Check for convergence.
    # 
    if iteration > 1 && abs(obj - old_obj) < epsilon * (abs(old_obj) + 1.0)
      return (estimate, obj)
    else
      old_obj = obj
    end
  end
  return (estimate, obj)
end # function regress

"""
Performs a score test for either linear, logistic, or Poisson regression 
with a canonical link function. X is the design matrix, y is the 
response vector, and estimate is the MLE under the null hypothesis.
"""
function glm_score_test(X::Matrix{Float64}, y::Vector{Float64}, 
  estimate::Vector{Float64}, model::AbstractString)

  if model != "linear" && model != "logistic" && model != "Poisson"
    throw(ArgumentError(
      "The only model choices are linear, logistic, and Poisson.\n \n"))
  end
  #
  # Initialize the score vector and information matrix.
  #
  (n, p) = size(X)
  @assert n == length(y)
  @assert p == length(estimate)
  score = zeros(p)
  information = zeros(p, p)
  z = zeros(n)
  #
  # Handle each model separately.
  #
  BLAS.gemv!('N', 1.0, X, estimate, 0.0, z) # z = X * estimate
  if model == "linear"
    BLAS.axpy!(-1.0, y, z) # z = z - y
    var_inv = n / norm(z)^2 # var = residual sum of squares / n
    score = BLAS.gemv('T', -var_inv, X, z) # score = - X' * (z - y) / var
    information = BLAS.gemm('T', 'N', var_inv, X, X) # informat = X' * X / var
  elseif model == "logistic"
    clamp!(z, -20.0, 20.0)
    z = exp.(-z)
    z = 1.0 ./ (1.0 + z)
    w = z .* (1.0 .- z)
    BLAS.axpy!(-1.0, y, z) # z = z - y
    score = BLAS.gemv('T', -1.0, X, z) # score = - X' * (z - y)
    w = sqrt.(w)
    scale!(w, X) # diag(w) * X
    information = BLAS.gemm('T', 'N', X, X) # information = X' * W * X 
    w = 1.0 ./ w
    scale!(w, X)
  elseif model == "Poisson"
    clamp!(z, -20.0, 20.0) 
    z = exp.(z)
    w = copy(z)
    BLAS.axpy!(-1.0, y, z) # z = z - y
    score = BLAS.gemv('T', -1.0, X, z) # score = - X' * (z - y)
    w = sqrt.(w)
    scale!(w, X) # diag(w) * X
    information = BLAS.gemm('T', 'N', X, X) # information = X' * W * X
    w = 1.0 ./ w
    scale!(w, X)
  end
  #
  # Compute the score statistic from the score and information.
  #
  z = information \ score
  score_test = dot(score, z)
  return score_test
end # function glm_score_test

# using GeneralUtilities
# n = 10
# a = blanks(n)
# println(a,"  ",typeof(a))
# b = empties(n)
# println(b,"  ",typeof(b))
# c = repeated_string(["av", "bc", "qrb", "bc"])
# println(c)
# integerset = IntSet([7, 5, 4, 10])
# d = select_set_element(integerset, 3)
# println(integerset,"  ",d)
# frequency = rand(n)
# frequency = frequency / sum(frequency)
# @time for i = 1:1000
#   s = random_category(frequency)
#   if s < 1 || s > n
#     println("category error")
#   end
# end
# n = 100
# x = randn(n)
# x[1] = NaN
# normalize!(x)
# println(sample_mean_std(x))
# s = sample_stats(x)
# print_sample_stats(s)
# pvalue = 0.1*rand(n)
# fdr = collect(0.1:0.1:0.5)
# tests = 1000
# (number_passing, threshold) = simes_fdr(pvalue, fdr, tests)
# X = [68., 49, 60, 68, 97, 82, 59, 50, 73, 39, 71, 95, 61, 72, 87, 
#   40, 66, 58, 58, 77]
# y = [75., 63, 57, 88, 88, 79, 82, 73, 90, 62, 70, 96, 76, 75, 85,
#   40, 74, 70, 75, 72]
# X = [ones(size(X,1)) X]
# (estimate, loglikelihood) = regress(X, y, "linear") # Jennrich section 1.4
# println(estimate)
# println(" ")
# X = [0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 1.75, 2.00, 2.25, 2.50, 2.75,
#   3.00, 3.25, 3.50, 4.00, 4.25, 4.50, 4.75, 5.00, 5.50]
# y = [0., 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1]
# X = [ones(size(X,1)) X]
# (estimate, loglikelihood) = regress(X, y, "logistic") # Wikipedia problem
# println(estimate)
# println(" ")
# X = [1.0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
# X = reshape(X, 14, 1)
# y = [0., 1 ,2 ,3, 1, 4, 9, 18, 23, 31, 20, 25, 37, 45]
# (estimate, loglikelihood) = regress(X, y, "Poisson") # Aids problem
# estimate = [0.0; estimate]
# println(estimate)
# X = [ones(size(X,1)) X]
# test = glm_score_test(X, y, estimate, "Poisson")
# println("score test = ",test)
# (estimate, loglikelihood) = regress(X, y, "Poisson") # Aids problem
# println(" ")
# println(estimate," ",loglikelihood) 