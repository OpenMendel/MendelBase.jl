using MendelBase, Base.Test, SnpArrays, DataFrames, DataStructures, CSV

info("Unit tests for read_data")

srand(123)

# read data from gamete example
keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
process_keywords!(keyword, "control.txt", []) 
computed = read_external_data_files(keyword)

reference = CSV.read("gamete competition PedigreeFrame.txt")
num = length(unique(reference[1]))
tot_num = length(reference[1])

@testset "pedigree constructor 1st test" begin
    computed_ped = computed[1]
    @test computed_ped.pedigrees == num
    x = zeros(num)
    for i in 1:length(unique(reference[1])) 
        x[i] = parse(Int64, computed_ped.name[i]) 
    end # puts unique pedigrees to vector x
    @test all(x .== unique(reference[1]))

    # from now on, will only test if the last entry is correct, since theres no point 
    # to write code to reproduce another code
    @test computed_ped.start[num] == 549
    @test computed_ped.twin_finish[num] == 553 # no twins in gametes example
    @test computed_ped.finish[num] == 553
    @test computed_ped.individuals[num] == 5
    @test computed_ped.founders[num] == 2 #in a pedigree, count # of people who don't know their father
    @test computed_ped.females[num] == 4 #in a pedigree, count # of females
    @test computed_ped.males[num] == 1
    @test all(computed_ped.twins .== 0) # since no twins
    @test computed_ped.families[num] == 1 # number of differernt pairs of parents
    @test all(computed_ped.loglikelihood .== 0.0) #loglikelihood never computed...?
end

@testset "person constructor 1st test" begin
    computed_person = computed[2]
    @test computed_person.people == 553
    @test computed_person.populations == 1 #population is a set with length 1
    @test computed_person.name[tot_num] == "5" # person's name is unique inside pedigrees
    @test computed_person.pedigree[tot_num] == 69 # 69th unique pedigree
    @test computed_person.mother[tot_num] == 550 # mother and father flipped in the first row of pedigree file?
    @test computed_person.father[tot_num] == 549  
    @test computed_person.male[tot_num] == false # person's name in each pedigree is sorted first
    @test all(computed_person.next_twin .== 0) 
    @test all(computed_person.primary_twin .== 0) 
    @test size(computed_person.admixture) == (tot_num, computed_person.populations)
    @test all(computed_person.admixture .== 1.0)
    @test computed_person.children[tot_num - 3] == IntSet([551, 552, 553])
    @test computed_person.spouse[tot_num - 3] == IntSet([549])
    @test size(computed_person.genotype) == (553, 10) # 553 lines with 9 SNPs + 1 CT
    @test size(computed_person.homozygotes) == (553, 10) 
    @test computed_person.homozygotes[tot_num,:] == [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    @test computed_person.disease_status[16] == "1" # first instance when ACE = 1
    @test size(computed_person.variable) == (553,0) #nothing has non-admixture variables
    @test computed_person.variable_name == [] # no variables
end

# nuclear family is defined as a family of only the spouse and their children, 
# so in this set of 69 pedigrees we have 115 nuclear families because people
# can be a part of more than 1 nuclear family.
@testset "NuclearFamily constructor 1st test" begin 
    computed_nuclear = computed[3]
    num = sum(computed[1].families)
    @test computed_nuclear.families == num 
    @test computed_nuclear.pedigree[num] == 69
    
    #first move all NA up, then sort everybody in the same pedigree based on "person" number
    #parents must come before children.
    @test computed_nuclear.mother[num] == 550
    @test computed_nuclear.father[num] == 549
    @test computed_nuclear.sib[num] == IntSet([551, 552, 553]) 
end

@testset "Locus constructor 1st test" begin
    computed_locus = computed[4]
    @test computed_locus.loci == 10
    @test computed_locus.model_loci == 10
    @test computed_locus.trait == 0 
    @test computed_locus.free_recombination == false #since all on the same chromosome
    @test length(computed_locus.name) == 10
    @test all(computed_locus.chromosome .== "AUTOSOME")
    @test all(computed_locus.base_pairs .== 0)
    @test size(computed_locus.morgans) == (2, computed_locus.loci)
    @test size(computed_locus.theta) == (2, computed_locus.loci - 1)
    @test all(computed_locus.theta .== map_function(1.0, "Haldane")) # d = 1.0 since 1.0 morgans...?
    @test all(computed_locus.xlinked .== false)
    @test all(computed_locus.alleles .== 2)
    @test typeof(computed_locus.allele_name) == Array{Array{AbstractString,1},1}
    @test eltype(computed_locus.frequency) == Array{Float64,2}
    @test eltype(computed_locus.lumped_frequency) == Float64
    @test length(computed_locus.model_locus) == 10
    @test length(computed_locus.locus_field_in_pedigree_frame) == 10
end





# read data from gwas example, inputted using plink files
keyword2 = set_keyword_defaults!(Dict{AbstractString, Any}())
process_keywords!(keyword2, "gwas 1 Control.txt", []) 
computed2 = read_external_data_files(keyword2)

@testset "pedigree constructor using plink" begin
    # Plink.fam files are the only file that does not require a header line!
    computed_ped = computed2[1]
    @test computed_ped.pedigrees == 2200 # everybody are from unique pedigrees
    @test computed_ped.name[2200] == "2200" # every pedigree name is unique, listed in order
    @test computed_ped.start == collect(1:1:2200)#since everybody unique pedigree, start = the position in list
    @test computed_ped.twin_finish == collect(1:1:2200)# no twins 
    @test computed_ped.finish == collect(1:1:2200)# since everybody unique pedigree, finish = position in list
    @test computed_ped.individuals == ones(2200) # 1 person in each pedigree
    @test computed_ped.founders == ones(2200) # of person who doesn't have parents = everyone
    @test computed_ped.females[2200] == 1 #last is a female
    @test computed_ped.males[2200] == 0 
    @test computed_ped.twins == zeros(2200) # again, no twins
    @test computed_ped.families == zeros(2200) #number of differernt pairs of parents, but nobody knows their parents
    @test all(computed_ped.loglikelihood .== 0.0) #loglikelihood never computed...?
end

@testset "person constructor using plink" begin
    computed_person = computed2[2]
    @test computed_person.people == 2200
    @test computed_person.populations == 1
    @test all(computed_person.name .== "1") #everyone is person 1 (2nd column of plink file)
    @test computed_person.pedigree == collect(1:1:2200)
    @test computed_person.mother == zeros(2200)
    @test computed_person.father == zeros(2200)
    @test computed_person.male[2200] == false
    @test all(computed_person.next_twin .== 0)
    @test all(computed_person.primary_twin .== 0)
    @test all(computed_person.admixture .== 1.0) #no race variation, so fill admixture with ones
    @test all(computed_person.children .== IntSet([])) #nobody has children
    @test all(computed_person.spouse .== IntSet([]))
    @test size(computed_person.genotype) == (2200, 0) # no SNPs
    @test size(computed_person.homozygotes) == (2200, 0)
    @test computed_person.disease_status == [] #no disease column for plink.fam
    @test size(computed_person.variable) == (2200, 1) # Trait is always one (and perhaps only) variable for plink.fam files
    @test computed_person.variable_name[1] == "Trait"
end

@testset "NuclearFamily constructor using plink" begin
    # there are no nuclear families in the gwas example because nobody
    # knows their parents. 
    computed_nuclear = computed2[3]
    @test computed_nuclear.families == 0 
    @test computed_nuclear.pedigree == []
    @test computed_nuclear.mother == []
    @test computed_nuclear.father == []
    @test computed_nuclear.sib == []
end

@testset "Locus constructor for plink files" begin
    # there really aren't too much to test since locus_frame is empty
    # and pedigree frames (read from plink.fam file) has no locus information, 
    # but this is still included for completion.
    computed_locus = computed2[4]
    @test computed_locus.loci == 0
    @test computed_locus.model_loci == 0
    @test computed_locus.trait == 0# this is always zero, as defined last line of locus_information?
    @test computed_locus.free_recombination == true # this is true but I couldn't figure out why
    @test computed_locus.name == []
    @test computed_locus.chromosome == []
    @test computed_locus.base_pairs == []
    @test size(computed_locus.morgans) == (2, computed_locus.loci)
    @test size(computed_locus.theta) == (2, computed_locus.loci)
    @test computed_locus.xlinked == []
    @test computed_locus.alleles == []
    @test computed_locus.allele_name[1][1] == ""
    @test computed_locus.frequency[1][1] == 0.0 
    @test computed_locus.lumped_frequency[1][1][1] == 0.0
    @test computed_locus.model_locus == []
    @test computed_locus.locus_field_in_pedigree_frame == []
end

# read data from gene dropping example
# used to test if phenotype frames are processed correctly
keyword3 = set_keyword_defaults!(Dict{AbstractString, Any}())
process_keywords!(keyword3, "genedropping Control.txt", []) 
computed3 = read_external_data_files(keyword3)

@testset "Some tests on phenotype frames" begin
    computed_locus = computed3[4]
    computed_locus.chromosome == ["1", "9", "X", "X"]
    computed_locus.name == ["Rh", "ABO", "Xg", "XSNP"]
    computed_locus.base_pairs == [0, 10, 0, 0]
    computed_locus.xlinked == [false, false, true, true]
    computed_locus.alleles == [2, 3, 2, 2]
    computed_locus.allele_name == [AbstractString["D","d"], 
        AbstractString["A","B","O"], AbstractString["+","-"],
        AbstractString["1","2"]]
    computed_locus.frequency[1] == [1.0 0.0 0.0; 0.81 0.19 0.0; 0.61 0.39 0.0]
    computed_locus.model_locus == [1, 2, 3, 4]
    computed_locus.locus_field_in_pedigree_frame == [13, 12, 14, 15]
end

@testset "Pedigree data structure" begin
    #
    # Uses the gamete competition PedigreeFrame as test case.
    # check founders at the head of the pedigree,
    # everyone else arranged so that parents precede their children.
    # (co-twins has to be at tail of the pedigree, but no co-twins in gametes example)
    # 
    computed_ped = computed[1]
    computed_person = computed[2]
    individual_count = computed_ped.individuals # used to count number of rows in data
    list_of_founders = computed_ped.founders # used to compute number of founder in each pedigree
    data = [computed_person.mother computed_person.father] #matrix of data after permuting individuals internally by open mendel
    
    founder_at_front = true
    parents_precede_child = true

    for i in 1:computed_ped.pedigrees 
        founder_count = list_of_founders[i]
        ped_start_row = computed_ped.start[i]

        for j in 1:founder_count
            if data[ped_start_row + (j - 1), 1] != 0 || 
                data[ped_start_row + (j - 1), 2] != 0
                founder_at_front = false # if founders has parent, then condition not met
            end 
        end

        num_of_non_founder = individual_count[i] - founder_count
        
        for j in 1:num_of_non_founder
            child_row = ped_start_row + founder_count + j - 1
            mother_row = data[child_row, 1]
            father_row = data[child_row, 2]
            if child_row <= mother_row || child_row <= father_row
                parents_precede_child = false # test if parent comes after children
            end
        end
    end

    @test founder_at_front == true
    @test parents_precede_child == true
end

@testset "gene counting" begin
    # this is EM alg
end

