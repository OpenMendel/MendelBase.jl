using MendelBase, Base.Test, StatsBase

info("Unit tests for genetic_utilities")

srand(123)
n, p = 1000, 100

@testset "blanks and empties constructors" begin
    @test_throws(ErrorException, blanks(-1))
    @test_throws(MethodError, blanks(rand()))
    @test_throws(MethodError, blanks(NaN))
    @test_throws(MethodError, blanks(Inf))
    blank0 = blanks(0)
    @test length(blank0) == 0
    x = rand(1:n)
    blankx = blanks(x)
    @test typeof(blankx) <: Array{AbstractString}
    @test eltype(blankx) == AbstractString
    @test length(blankx) == x

    @test_throws(ErrorException, empties(-1))
    @test_throws(MethodError, empties(rand()))
    @test_throws(MethodError, empties(NaN))
    @test_throws(MethodError, empties(Inf))
    empty0 = empties(0)
    @test length(empty0) == 0
    emptyx = empties(x)
    @test typeof(emptyx) <: Array{IntSet}
    @test eltype(emptyx) == IntSet
    @test length(emptyx) == x
end

@testset "repeated_string" begin
    a = repeated_string(["hi", "ho", "ha", "he"])
    @test a == (false, "")
    b = repeated_string(["a", "b", "c", "d", "d"])
    @test b == (true, "d")
    c = repeated_string(["*&^%#", "*&^%#", "f32<", "*jf*ijioej2"])
    @test c == (true, "*&^%#")
    d = repeated_string(["aunt", "uncle", "aunt", "uncle", "aunt"])
    @test d == (true, "aunt")
end

@testset "random_category" begin
    freq = [0.1, 0.2, 0.7] #category 1, 2, 3
    cat_1 = 0
    cat_2 = 0
    cat_3 = 0
    for i in 1:100000
        test = random_category(freq)
        if test == 1 cat_1 += 1 end
        if test == 2 cat_2 += 1 end
        if test == 3 cat_3 += 1 end
    end
    @test round(cat_1/10000, 1) == 1 #test if frequency of cat_1 approximately 10%
    @test round(cat_2/10000, 1) == 2
    @test round(cat_3/10000, 1) == 7
end

@testset "select_set_element" begin
    firstset = IntSet([1, 4, 555, 3, 23, 40, 23456])
    secondset = IntSet([2, 2, 2, 2, 2])
    @test select_set_element(firstset, 5) == 40
    @test select_set_element(secondset, 1) == 2
    @test typeof(select_set_element(secondset, 2)) == Void
    @test typeof(select_set_element(secondset, 20)) == Void
end

@testset "normalize!" begin
    x = [1.0, 2.0, 3.0] #note the importance of comma
    MendelBase.normalize!(x) # must use MendelBase.normalize! because Base package also have a normalize! function
    @test mean(x) ≈ 0.0
    @test var(x) * 2.0 / 3.0 ≈ 1.0 # var(x) returns sample variance, we compute population variance

    x = rand(n)
    MendelBase.normalize!(x)
    @test round(mean(x), 10) == 0.0 # need to round because julia apparently thinks 1e-16 is not approximately 0
    @test var(x) * 999.0 / 1000.0 ≈ 1.0 

    x = rand(n)
    x[1:120] = NaN
    x[801:1000] = NaN # insert some NaN
    MendelBase.normalize!(x)
    @test round(mean(x[121:800]), 10) == 0.0
    @test var(x[121:800]) * 679 / 680 ≈ 1.0 
end

@testset "sample_mean_std" begin
    x = [1.0, 2.0, 3.0]
    result = sample_mean_std(x)
    @test result[1] ≈ 2.0
    @test result[2] ≈ sqrt(var(x) * 2.0 / 3.0)

    x = randn(n)
    result = sample_mean_std(x)
    @test result[1] ≈ mean(x)
    @test result[2] ≈ sqrt(var(x) * 999.0 / 1000.0) #var(x) calculates sample variance, we need population variance

    x = rand(n)
    x[1:120] = NaN
    x[801:1000] = NaN # insert some NaN
    result = sample_mean_std(x)
    @test result[1] ≈ mean(x[121:800])
    @test result[2] ≈ sqrt(var(x[121:800]) * 679 / 680)
end

@testset "sample_stats" begin
    v = [1.0, 2.0, NaN, 4.0, NaN, -6.0, NaN]
    u = [1.0, 2.0, 4.0, -6.0]
    result = sample_stats(v)
    @test result[1] == 4
    @test result[2] == 3
    @test result[3] == -6.0
    @test result[4] == quantile(u, 0.25)
    @test result[5] == 1.5
    @test result[6] == quantile(u, 0.75)
    @test result[7] == 4.0
    @test result[8] == 1/4 # NaN values does not increase denominator
    @test result[9] == std(u)
    @test result[10] == skewness(u)
    @test result[11] == kurtosis(u)

    #test above again but with length 1000 random vectors including NaN's
    x = randn(n)
    y = x
    @test typeof(sample_stats(x)) == Tuple{Int64,Int64,Float64,
    Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64}
    for i in 1:10 push!(y, NaN) end
    append!(y, x)
    append!(x, x) # vector y is vector x with an extra NaN in the middle.
    @test sample_stats(y) == sample_stats(x)
end

@testset "simes false discovery rate" begin
    #only testing basic stuff
    y = rand(n)
    fdr = rand(n)
    tests = 10
    result = simes_fdr(y, fdr, tests)
    @test typeof(result) <: Tuple{Array{Int64}, Array{Float64}}
    @test length(result[1]) == n
    @test length(result[2]) == n
end

@testset "regress" begin
    X = rand(n, p)
    @test_throws(AssertionError, regress(X, rand(n - 1), "logistic")) # matrix dimension mismatch
    y = rand(n)
    @test_throws(ArgumentError, regress(X, y, "ahuehuehue"))

    test1 = regress(X, y, "logistic")
    @test typeof(test1) <: Tuple{Array, Float64}
    @test length(test1[1]) == p
    @test eltype(test1[1]) == Float64
    @test typeof(test1[2]) == Float64

    # alligator example for linear regression 
    # https://www.r-bloggers.com/simple-linear-regression-2/
    X = [3.87, 3.61, 4.33, 3.43, 3.81, 3.83, 3.46, 3.76, 3.50, 3.58, 4.19, 3.78, 3.71, 3.73, 3.78]
    y = [4.87, 3.93, 6.46, 3.33, 4.38, 4.70, 3.50, 4.50,
           3.58, 3.64, 5.90, 4.43, 4.38, 4.42, 4.25]
    X = [ones(size(X)) X]
    (estimate, loglikelihood) = regress(X, y, "linear")
    @test round(estimate[1], 4) == -8.4761
    @test round(estimate[2], 4) == 3.4311

    #wikipedia example for logistic regression
    X = [0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 1.75, 2.00, 2.25, 2.50, 2.75,
        3.00, 3.25, 3.50, 4.00, 4.25, 4.50, 4.75, 5.00, 5.50]
    y = [0., 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1]
    X = [ones(size(X,1)) X]
    (estimate, loglikelihood) = regress(X, y, "logistic")
    @test round(estimate[1], 4) == -4.0777
    @test round(estimate[2], 4) == 1.5046

    # random online example for poisson (page 10)
    # http://www.biostat.umn.edu/~dipankar/bmtry711.11/lecture_13.pdf
    X = [1.0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
    y = [0.0, 1, 2, 3, 1, 4, 9, 18, 23, 31, 20, 25, 37, 45]
    X = [ones(size(X)) X]
    (estimate, loglikelihood) = regress(X, y, "Poisson")
    @test round(estimate[1], 4) == 0.3396
    @test round(estimate[2], 4) == 0.2565
    @test round(loglikelihood, 4) == 472.0625
end

@testset "glm_score_test" begin
    # find some examples online... 
end