# s is too large
@test_throws AssertionError ReductionProblem(f, x, p, 2)

# dimensions match
@test length(idx) == length(V) == length(dim_Y)

# type checks
@test isa(prob, ReductionProblem) 

@test prob.π == prob.θ[prob.idx_slow_fast] 
i = findall([i == [true, true, false] for i in idx])
@test length(i) == 1 
i = i[1]
