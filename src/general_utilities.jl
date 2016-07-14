################################################################################
# This set of functions provides various general utilities.
################################################################################

export empties, blanks, repeated_string, select_set_element, random_category
export normalize!, print_sample_stats, sample_mean_std, sample_stats

"""
This function creates an array of blank ASCII strings.
"""
function blanks(n::Int)

  blank = Array(ASCIIString, n)
  for i = 1:n
    blank[i] = ""
  end
  return blank
end # function blanks

"""
This function creates an array of empty integer sets.
"""
function empties(n::Int)

  empty = Array(IntSet, n)
  for i = 1:n
    empty[i] = IntSet()
  end
  return empty
end # function empties

"""
This function identifies a repeated string in an array of strings.
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
This function selects element number k of the integer set S.
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
end # function select_set_element

"""
This function returns a random category (numbered 1, 2, and so forth)
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
This function normalizes the vector x. Missing values are ignored.
"""
function normalize!(x::Vector{AbstractFloat})

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
This function computes a sample mean and standard deviation.
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
This function computes sample statistics for the data
vector x. Missing values equal NaN. The output consists
of the number of values present, the number of values
missing, the minimum, the 0.25 quantile, the median, the
0.75 quantile, the maximum, the mean, the standard
deviation, the skewness, and the excess kurtosis.
"""
function sample_stats(x::Vector{AbstractFloat})

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
This function prints the sample statistics output by the function sample_stats.
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

