# test construction of problem
R, (e₀, k₁, k₋₁, k₂, S, C) = polynomial_ring(QQ, ["e₀", "k₁", "k₋₁", "k₂", "S", "C"])
Fp, _ = rational_function_field(QQ, ["e₀", "k₁", "k₋₁", "k₂"])
Rx, _ = polynomial_ring(Fp, ["S", "C"])
problem_manual = ReductionProblem(
  [-e₀*k₁*S + k₁*S*C + k₋₁*C, e₀*k₁*S - k₁*S*C - k₋₁*C - k₂*C],
  [S, C],
  [e₀, k₁, k₋₁, k₂],
  1,
  matrix([-e₀*k₁+k₁*C k₁*S+k₋₁; e₀*k₁-k₁*C -k₁*S-k₋₁-k₂]),
  f,
  parent(S//C),
  Fp,
  Rx,
  f(gens(Rx), gens(Fp)),
  gens(Rx)
)
@test all([getfield(problem,fn) == getfield(problem_manual,fn) for fn in fieldnames(ReductionProblem)])

# find tfpv candidates and slow manifolds
tfpvs, varieties = tfpvs_and_varieties(problem);

# use predetermined parameters
preset = (k₋₁ = 1,)
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
