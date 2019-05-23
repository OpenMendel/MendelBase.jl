################################################################################
# This set of functions provides various general utilities.
################################################################################
#
# Required external packages.
#
# using Distributions
# using LinearAlgebra
# using Statistics

export blanks, empties, repeated_string, select_set_element, random_category
export standarize!, sample_mean_std, sample_stats, print_sample_stats
export simes_fdr, fast_regress, fast_score_test
export glm_regress, glm_score_statistic

"""
Creates an array of blank strings.
"""
function blanks(n::Int)
#
  blank = Array{String}(undef, n)
  for i = 1:n
    blank[i] = ""
  end
  return blank
end # function blanks

"""
Creates an array of empty ordered sets.
"""
function empties(n::Int)
#
  empty = Array{BitSet}(undef, n)
  for i = 1:n
    empty[i] = BitSet()
  end
  return empty
end # function empties

"""
Identifies a repeated string in an array of strings.
"""
function repeated_string(name::Array{String, 1})
#
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
function select_set_element(set_s::BitSet, k::Int)
#
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
function random_category(frequency::Vector{T} where T <: Real)
#
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
Standardizes the vector x to have mean 0.0 and variance 1.0.
Missing values are ignored.
"""
function standarize!(x::Vector) 
#
  avg = mean(skipmissing(x))
  stdev = std(skipmissing(x))
  if stdev > zero(eltype(stdev))
    x .= (x .- avg) ./ stdev
  else
    x .= x .- avg
  end
end # function normalize!

"""
Computes the sample mean and standard deviation.
Missing values are ignored.
"""
function sample_mean_std(x::Vector)
#
  return (mean(skipmissing(x)), std(skipmissing(x)))
end # function sample_mean_std

"""
Computes sample statistics for the data vector x.
Missing values are ignored. The output consists of
the number of values present, the number of values missing,
the minimum, the 0.25 quantile, the median, the 0.75 quantile,
the maximum, the mean, the standard deviation, the skewness,
and the excess kurtosis.
"""
function sample_stats(x::Vector)
#
  y = collect(skipmissing(x))
  (p, n) = (length(y), length(x))
  return p, n - p, minimum(y), quantile(y, 0.25), median(y),
    quantile(y, 0.75), maximum(y), mean(y), std(y),
    skewness(y), kurtosis(y)
end # function sample_stats

"""
Prints the sample statistics output by the function sample_stats.
"""
function print_sample_stats(s, io::IO = STDOUT)
#
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
function simes_fdr(pvalue::Vector{T} where T <: Real,
                   fdr::Vector{T} where T <: Real, tests::Int)
#
  n = length(fdr)
  threshold = zeros(typeof(pvalue[1]), n)
  number_passing = zeros(Int, n)
  perm = sortperm(pvalue)
  k = 1 # previous value of i
  for j = 1:n
    i = 0
    for m = k:length(pvalue)
      i = m
      if tests * pvalue[perm[m]] > m * fdr[j]; break; end
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
function fast_regress(X::Matrix{Float64}, y::Vector{Float64}, model::AbstractString)

  if model != "linear" && model != "logistic" && model != "poisson"
    throw(ArgumentError(
      "The only model choices are linear, logistic, and poisson.\n \n"))
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
  elseif model == "poisson"
    estimate[1] = log(mean(y))
  else
    throw(ArgumentError(
      "The only model choices are linear, logistic, and poisson.\n \n"))
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
      z = 1.0 ./ (1.0 .+ z)
      w = z .* (1.0 .- z)
      BLAS.axpy!(-1.0, y, z) # z = z - y
      score = BLAS.gemv('T', -1.0, X, z) # score = - X' * (z - y)
      w = sqrt.(w)
      lmul!(Diagonal(w), X) # X = diag(w) * X
      information = BLAS.gemm('T', 'N', X, X) # information = X' * W * X
      w = 1.0 ./ w
      lmul!(Diagonal(w), X) # X = diag(w) * X
    elseif model == "poisson"
      z = exp.(z)
      w = copy(z)
      BLAS.axpy!(-1.0, y, z) # z = z - y
      score = BLAS.gemv('T', -1.0, X, z) # score = - X' * (z - y)
      w = sqrt.(w)
      lmul!(Diagonal(w), X) # X = diag(w) * X
      information = BLAS.gemm('T', 'N', X, X) # information = X' * W * X
      w = 1.0 ./ w
      lmul!(Diagonal(w), X) # X = diag(w) * X
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
        z = 1.0 ./ (1.0 .+ z)
        for i = 1:n
          if y[i] > 0.0
            obj = obj + log(z[i])
          else
            obj = obj + log(1.0 - z[i])
          end
        end
      elseif model == "poisson"
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
end # function fast_regress

"""
Performs a score test for either linear, logistic, or Poisson regression 
with a canonical link function. X is the design matrix, y is the 
response vector, and estimate is the MLE under the null hypothesis.
"""
function fast_score_test(X::Matrix{Float64}, y::Vector{Float64}, 
  estimate::Vector{Float64}, model::AbstractString)

  if model != "linear" && model != "logistic" && model != "poisson"
    throw(ArgumentError(
      "The only model choices are linear, logistic, and poisson.\n \n"))
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
    z = 1.0 ./ (1.0 .+ z)
    w = z .* (1.0 .- z)
    BLAS.axpy!(-1.0, y, z) # z = z - y
    score = BLAS.gemv('T', -1.0, X, z) # score = - X' * (z - y)
    w = sqrt.(w)
    lmul!(Diagonal(w), X) # X = diag(w) * X
    information = BLAS.gemm('T', 'N', X, X) # information = X' * W * X 
    w = 1.0 ./ w
    lmul!(Diagonal(w), X) # X = diag(w) * X
  elseif model == "poisson"
    clamp!(z, -20.0, 20.0) 
    z = exp.(z)
    w = copy(z)
    BLAS.axpy!(-1.0, y, z) # z = z - y
    score = BLAS.gemv('T', -1.0, X, z) # score = - X' * (z - y)
    w = sqrt.(w)
    lmul!(Diagonal(w), X) # X = diag(w) * X
    information = BLAS.gemm('T', 'N', X, X) # information = X' * W * X
    w = 1.0 ./ w
    lmul!(Diagonal(w), X) # X = diag(w) * X
  end
  #
  # Compute the score statistic from the score and information.
  #
  z = information \ score
  score_test = dot(score, z)
  return score_test
end # function fast_score_test

"""
Performs generalized linear regression.
X is the design matrix, y is the response vector,
meanf is the value and derivative of the inverse link,
and varf is the variance function of the mean.
"""
function glm_regress(X::Matrix, y::Vector, meanf::Function, varf::Function,
  loglikelihood::Function)
#
  (n, p) = size(X)
  @assert n == length(y)
  (score, inform, beta) = (zeros(p), zeros(p, p), zeros(p))
  (x, z) = (zeros(p), zeros(n))
  ybar = mean(y)
  for iteration = 1:20 # find the intercept by Newton's method
    g = meanf(beta[1])
    beta[1] = beta[1] - clamp((g[1] - ybar) / g[2], -1.0, 1.0)
    if abs(g[1] - ybar) < 1e-10
      break
    end
  end
  (obj, old_obj, c, v) = (0.0, 0.0, 0.0, 0.0)
  epsilon = 1e-8
  for iteration = 1:100 # scoring algorithm
    fill!(score, 0.0)
    fill!(inform, 0.0)
    mul!(z, X, beta) # z = X * beta
    for i = 1:n
      f = meanf(z[i])
      v = varf(f[1])
      c = ((y[i] - f[1]) / v) * f[2]
      copyto!(x, X[i, :])
      BLAS.axpy!(c, x, score) # score = score + c * x
      c = f[2]^2 / v
      BLAS.ger!(c, x, x, inform) # inform = inform + c * x * x'
    end
    increment = inform \ score
    beta = beta + increment
    steps = -1
    fill!(score, 0.0)
    for step_halve = 0:3 # step halving
      obj = 0.0
      mul!(z, X, beta) # z = X * beta
      steps = steps + 1
      for i = 1:n
        f = meanf(z[i])
        v = varf(f[1])
        c = ((y[i] - f[1]) / v) * f[2]
        copyto!(x, X[i, :])
        BLAS.axpy!(c, x, score) # score = score + c * x
        obj = obj + loglikelihood(y[i], f[1]) 
      end
      if obj > old_obj
        break
      else
        beta = beta - increment
        increment = 0.5 * increment
      end
    end
    println(iteration," ",old_obj," ",obj," ",steps)
    if iteration > 1 && abs(obj - old_obj) < epsilon * (abs(old_obj) + 1.0)
      return (beta, obj)
    else
      old_obj = obj
    end
  end
  return (beta, obj)
end # function glm_regress

"""
Computes the score statistic for a generalized linear regression. 
X is the design matrix, y is the response vector,
beta is the vector of regression coefficients,
meanf is the value and derivative of the inverse link,
and varf is the variance function of the mean.
"""
function glm_score_statistic(X::Matrix, y::Vector, beta::Vector, 
  meanf::Function, varf::Function)
#
  (n, p) = size(X)
  @assert n == length(y)
  @assert p == length(beta)
  (score, inform) = (zeros(p), zeros(p, p))
  (x, z) = (zeros(p), zeros(n))
  mul!(z, X, beta) # z = X * beta
  for i = 1:n
    f = meanf(z[i])
    v = varf(f[1])
    c = ((y[i] - f[1]) / v) * f[2]
    copyto!(x, X[i, :])
    BLAS.axpy!(c, x, score) # score = score + c * x
    c = f[2]^2 / v
    BLAS.ger!(c, x, x, inform) # inform = inform + c * x * x'
  end
  x = inform \ score
  score_statistic = dot(score, x)
  return score_statistic
end # function glm_score_statistic


# Testing Routines
#
# n = 10
# a = blanks(n)
# println(a,"  ",typeof(a))
# b = empties(n)
# println(b,"  ",typeof(b))
# c = repeated_string(["av", "bc", "qrb", "bc"])
# println(c)
# integerset = BitSet([7, 5, 4, 10])
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
# x = Vector{Union{Missing, Float64}}(undef, n);
# x[1:n - 1] = randn(n - 1)
# println(typeof(x), size(x))
# standarize!(x)
# println(sample_mean_std(x))
# s = sample_stats(x)
# print_sample_stats(s)
# pvalue = 0.1*rand(n)
# fdr = collect(0.1:0.1:0.5)
# tests = 1000
# (number_passing, threshold) = simes_fdr(pvalue, fdr, tests)
# function GaussMean(u)
#   return [u, one(eltype(u))]
# end
# function GaussVar(mu)
#   return one(eltype(mu))
# end
# function GaussLoglikelihood(y, mu)
#   - (y - mu)^2 
# end
# X = [68., 49, 60, 68, 97, 82, 59, 50, 73, 39, 71, 95, 61, 72, 87, 
#   40, 66, 58, 58, 77];
# y = [75., 63, 57, 88, 88, 79, 82, 73, 90, 62, 70, 96, 76, 75, 85,
#   40, 74, 70, 75, 72];
# X = [ones(size(X, 1)) X];
# # Jennrich section 1.4 least squares problem
# (beta, obj) = glm(X, y, GaussMean, GaussVar, GaussLoglikelihood) 
# println(beta," ",obj)
# println(" ")
# function LogisticMean(u)
#   p = exp(u) 
#   p = p / (1 + p)
#   return [p, p * (1 - p)]
# end
# function LogisticVar(mu)
#   return mu * (1 - mu)
# end
# function LogisticLoglikelihood(y, mu)
#   y * log(mu) + (1 - y) * log(1 - mu) 
# end
# X = [0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 1.75, 2.00, 2.25, 2.50, 2.75,
#   3.00, 3.25, 3.50, 4.00, 4.25, 4.50, 4.75, 5.00, 5.50];
# y = [0., 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1];
# X = [ones(size(X,1)) X];
# # Wikipedia logistic regression problem
# (beta, obj) = glm(X, y, LogisticMean, LogisticVar, LogisticLoglikelihood)
# println(beta," ",obj)
# function PoissonMean(u)
#   p = exp(u)
#   return [p, p]
# end
# function PoissonVar(mu)
#   return mu
# end
# function PoissonLoglikelihood(y, mu)
#   y * log(mu) - mu 
# end
# println(" ")
# X = [1.0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14];
# X = reshape(X, 14, 1);
# y = [0., 1 ,2 ,3, 1, 4, 9, 18, 23, 31, 20, 25, 37, 45];
# # AIDs Poisson regression problem
# (beta, obj) = glm(X, y, PoissonMean, PoissonVar, PoissonLoglikelihood) 
# println(beta," ",obj)
# beta = [zero(1); beta];
# X = [ones(size(X,1)) X];
# score_statistic = glm_score_statistic(X, y, beta, PoissonMean, PoissonVar)
# println("score_statistic = ",score_statistic)
