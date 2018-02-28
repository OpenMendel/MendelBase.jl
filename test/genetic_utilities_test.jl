using MendelBase, Base.Test, Distributions

info("Unit tests for genetic_utilities")

srand(123)
n = 1000

@testset "map_function" begin
    # input x = genetic distance, i.e. must be greater than 1
    # output = recombination fraction, which will always be < 0.5 
    x = rand() 
    h = map_function(x, "Haldane")
    k = map_function(x, "Kosambi")

    @test typeof(h) == Float64
    @test typeof(k) == Float64
    @test_throws(ArgumentError, map_function(x, "RandomString"))
    @test_throws(DomainError, map_function(-x, "Haldane"))
    @test_throws(DomainError, map_function(-Inf, "Haldane"))
    @test map_function(Inf, "Kosambi") == 0.5
    @test map_function(0.0, "Kosambi") == 0.0
    @test map_function(Inf, "Haldane") == 0.5
    @test map_function(0.0, "Haldane") == 0.0
    @test isnan(map_function(NaN, "Haldane")) 
    @test isnan(map_function(NaN, "Kosambi")) 

    for i in 1:10
        x = rand()
        h = map_function(x, "Haldane")
        k = map_function(x, "Kosambi")
        @test h ≈ 0.5 * (1.0 - exp(-2*x))
        @test k ≈ 0.5 * tanh(2*x)
    end
end

@testset "inverse_map_function" begin
    # input = recombination fraction (i.e. some probability) which must be 
    # greater than 0, and output should be infinite if input >= 0.5
    @test inverse_map_function(0.0, "Kosambi") == 0.0
    @test inverse_map_function(0.5, "Haldane") == Inf
    @test inverse_map_function(0.5, "Kosambi") == Inf
    @test_throws(DomainError, inverse_map_function(1.1, "Haldane"))
    @test_throws(DomainError, inverse_map_function(-rand(), "Kosambi"))
    @test_throws(DomainError, inverse_map_function(-Inf, "Kosambi"))
    @test isnan(inverse_map_function(NaN, "Kosambi"))
    @test isnan(inverse_map_function(NaN, "Haldane"))

    for i in 1:10
        x = rand(Uniform(0, 0.5))
        h = inverse_map_function(x, "Haldane")
        k = inverse_map_function(x, "Kosambi")
        @test h ≈ -0.5 * log(1.0 - 2.0 * x)
        @test k ≈ 0.25 * log((1.0 + 2.0 * x) / (1.0 - 2.0 * x))
    end
end

@testset "hardy_weinberg_test" begin
    x = fill(NaN, n)
    @test hardy_weinberg_test(x) == 1.0
    y = rand(n)
    @test typeof(hardy_weinberg_test(y)) == Float64

    # flower example http://www.radford.edu/rsheehy/Gen_flash/Tutorials/Chi-Square_tutorial/x2-tut.htm
    for i in 1:200 x[i] = 0.0 end
    for i in 201:600 x[i] = 1.0 end
    for i in 601:1000 x[i] = 2.0 end
    result = hardy_weinberg_test(x)
    ans = ccdf(Chisq(1), 27.77777777)
    @test round(result * 10e5, 6) == round(ans * 10e5, 6) #julia still lacks a significant digit rounding funciton... lol
end

@testset "xlinked_hardy_weinberg_test" begin
    # problem 4.9.2 in Ken Lange's statistical genetics book
    # 
    # Frequency of t is (2 * 63 + 55 + 74) / (2 * 130 + 112) ≈ 0.6855
    # Frequency of y is thus 1 - 0.6855 ≈ 0.3145
    # so expected female tt, for example, is 0.6855^2 * 130
    # and expected female ty = 2 * 0.6855 * 0.3145 * 130
    #
    #                Observed   expected       (obs - exp)      (obs - exp)^2 / exp
    #           tt      63       61.0855        1.9145          0.0600
    #females    ty      55       56.0533        1.0533          0.01979 
    #           yy      12       12.8583        0.8583          0.05729
    #
    #males      t       74       76.776         2.776           0.1004
    #           y       38       35.224         2.776           0.2188 
    #
    # so the sum of last column ~ chi square distribution with 2 degress of freedom.

    x = fill(NaN, n)
    males = trues(n)
    @test xlinked_hardy_weinberg_test(x, males) == 1.0
    y = rand(n)
    @test typeof(xlinked_hardy_weinberg_test(y, males)) == Float64

    x = zeros(242)
    males = trues(242)
    for i in 1:130 males[i] = false end # first 130 person are females

    counter = 1
    for i in 1:63 x[counter] = 0.0; counter += 1 end
    for i in 1:55 x[counter] = 1.0; counter += 1 end
    for i in 1:12 x[counter] = 2.0; counter += 1 end
    for i in 1:74 x[counter] = 0.0; counter += 1 end
    for i in 1:38 x[counter] = 2.0; counter += 1 end
    result = xlinked_hardy_weinberg_test(x, males)
    ans = ccdf(Chisq(2), 0.4568) # see above how to obtain 0.4568
    @test floor(result * 10e1) == floor(ans * 10e1) #manual calculation only agrees to 2 decimal places
end





