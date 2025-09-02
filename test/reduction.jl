# Reduction instance 
red_i = Reduction(problem, tfpvs[i])
@test isa(red_i, Reduction)

# make variables available
S, C = system_components(problem)
e₀, k₁, k₋₁, k₂ = system_parameters(problem)

# test splitting of system
@test red_i.f0 == [-(k₁*e₀*S - (k₁*S + k₋₁)*C), k₁*e₀*S - (k₁*S + k₋₁)*C]
@test red_i.f1 == [0, -k₂*C]

# set everything that is requered for computation of reduced system
@test set_manifold!(red_i, [S, e₀*k₁*S//(k₁*S + k₋₁)])
@test set_point!(red_i, [S, e₀*k₁*S//(k₁*S + k₋₁)])
@test set_decomposition!(red_i, [1; -1], [k₁*S*C - e₀*k₁*S + k₋₁*C])
@test set_decomposition!(red_i, [k₁*S*C - e₀*k₁*S + k₋₁*C])
@test set_decomposition!(red_i, varieties[i][1])

# compute reduction
@test compute_reduction!(red_i)
@test red_i.g[1] == (-k₂*k₁*e₀*S*(k₁*S + k₋₁)//(k₁*k₋₁*e₀ + (k₁*S + k₋₁)^2))

# convenience function for reduction
red_i_conv = Reduction(problem, tfpvs[i], varieties[i][1], parent(S//C).([S, e₀*k₁*S//(k₁*S + k₋₁)]))
@test all([getfield(red_i_conv, name) == getfield(red_i, name) for name in fieldnames(Reduction)])

# find all possible slow varieties
all_varieties = unique_varieties(problem, varieties)

# explicit description of varieties
M = [
  problem._F.([S, 0]),
  problem._F.([-k₋₁//k₁, C]),
  problem._F.([k₋₁*C//(k₁*(e₀ - C)), C])
]
M_auto = [get_explicit_manifold(problem, v) for v in all_varieties]
@test all(all(M .== [m[1] for m in M_auto]))

# compute all reductions
R, idx_M = compute_all_reductions(problem, tfpvs, varieties, [m[1] for m in M_auto])

# Access the `Reduction` object 
reduction_4 = R[4][1]
@test reduction_4.M == [S, 0]
@test Matrix(reduction_4.Df0) == [-e₀*k₁+k₁*C k₁*S+k₋₁; e₀*k₁-k₁*C -k₁*S-k₋₁-k₂]
@test Matrix(reduction_4.Df0_at_x0) == [0 k₁*S+k₋₁; 0 -k₁*S-k₋₁-k₂]
@test all(Matrix(reduction_4.P) .== [S*k₁ + k₋₁; -S*k₁ - k₋₁ - k₂])
@test all(Matrix(reduction_4.Psi) .== [C])

# and the corresponding reduced system
@test reduction_4.g[1] == (-S*e₀*k₁*k₂)//(S*k₁ + k₋₁ + k₂)
