# find tfpv candidates and slow manifolds
tfpvs, varieties = tfpvs_and_varieties(problem);

# use predetermined parameters
preset = (k₋₁ = 1) 
tfpvs, varieties = tfpvs_and_varieties(problem; preset=preset);

# s is too large
@test_throws AssertionError ReductionProblem(f, x, p, 2)

# dimensions match
@test length(tfpvs) == length(varieties) 

# type checks
@test isa(problem, ReductionProblem) 
@test isa(tfpvs, Vector{Vector{Bool}})
@test isa(varieties, Vector{Vector{Variety}})

i = findall([i == [true, true, true, false] for i in tfpvs])
i = i[1]
@test i == 7
