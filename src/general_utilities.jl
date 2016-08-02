################################################################################
# This set of functions provides various general utilities.
################################################################################

export empties, blanks, repeated_string, select_set_element, random_category
export normalize!, print_sample_stats, sample_mean_std, sample_stats, simes_fdr

"""
Create an array of blank ASCII strings.
"""
function blanks(n::Int)

  blank = Array(ASCIIString, n)
  for i = 1:n
    blank[i] = ""
  end
  return blank
end # function blanks

"""
Create an array of empty integer sets.
"""
function empties(n::Int)

  empty = Array(IntSet, n)
  for i = 1:n
    empty[i] = IntSet()
  end
  return empty
end # function empties

"""
Identify a repeated string in an array of strings.
"""
function repeated_string(name::Vector{ASCIIString})

  sorted = sort(name)
  for i = 2:length(name)
    if sorted[i - 1] == sorted[i]
      return (true, sorted[i])
    end
  end
  return (false, "")
end # function repeated_string

"""
Select element number k of the integer set S.
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
Return a random category (numbered 1, 2, and so forth)
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
Normalize the vector x. Missing values are ignored.
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
Compute a sample mean and standard deviation.
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
Compute sample statistics for the data vector x.
Missing values equal NaN. The output consists
of the number of values present, the number of values
missing, the minimum, the 0.25 quantile, the median, the
0.75 quantile, the maximum, the mean, the standard
deviation, the skewness, and the excess kurtosis.
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
Print the sample statistics output by the function sample_stats.
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
Perform the Simes false discovery rate (FDR) procedure discussed 
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
