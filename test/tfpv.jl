# find tfpv candidates and slow manifolds
tfpvs, varieties = tfpvs_and_varieties(problem);

# s is too large
@test_throws AssertionError ReductionProblem(f, x, p, 2)

# dimensions match
@test length(tfpvs) == length(varieties) 

# type checks
@test isa(problem, ReductionProblem) 
@test isa(tfpvs, Vector{Vector{Bool}})
@test isa(varieties, Vector{Vector{Variety}})

@test problem.p_sf == problem.p[problem.idx_slow_fast] 
i = findall([i == [true, true, false] for i in tfpvs])
@test length(i) == 1 
i = i[1]
