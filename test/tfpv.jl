# find tfpv candidates and slow manifolds
tfpvs, manifolds = tfpvs_and_manifolds(prob);

# s is too large
@test_throws AssertionError ReductionProblem(f, x, p, 2)

# dimensions match
@test length(tfpvs) == length(manifolds) 

# type checks
@test isa(prob, ReductionProblem) 
@test isa(tfpvs, Vector{Vector{Bool}})
@test isa(manifolds, Vector{Vector{SlowManifold}})

@test prob.p_sf == prob.p[prob.idx_slow_fast] 
i = findall([i == [true, true, false] for i in tfpvs])
@test length(i) == 1 
i = i[1]
