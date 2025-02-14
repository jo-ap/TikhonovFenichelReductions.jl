# s is too large
@test_throws AssertionError ReductionProblem(f, x, p, 2)

# dimensions match
@test length(sf_separations) == length(V) == length(dim_V)

# type checks
@test isa(prob, ReductionProblem) 
@test isa(sf_separations, Vector{Vector{Bool}})
@test isa(V, Vector{Vector{Vector{TikhonovFenichelReductions.FracFieldElem{TikhonovFenichelReductions.QQMPolyRingElem}}}})
@test isa(dim_V, Vector{Vector{Int}})

@test prob.p_sf == prob.p[prob.idx_slow_fast] 
i = findall([i == [true, true, false] for i in sf_separations])
@test length(i) == 1 
i = i[1]
