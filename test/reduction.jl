# Reduction instance 
red_i = Reduction(prob, tfpvs[i])
@test isa(red_i, Reduction)

# make variables available
S, C = system_components(prob)
e₀, k₁, k₋₁, k₂ = system_parameters(prob)

# test splitting of system
@test red_i.f0 == [-(k₁*e₀*S - (k₁*S + k₋₁)*C), k₁*e₀*S - (k₁*S + k₋₁)*C]
@test red_i.f1 == [0, -k₂*C]

# set everything that is requered for computation of reduced system
@test set_manifold!(red_i, [S, e₀*k₁*S//(k₁*S + k₋₁)])
@test set_point!(red_i, [S, e₀*k₁*S//(k₁*S + k₋₁)])
@test set_decomposition!(red_i, [1; -1], [k₁*S*C - e₀*k₁*S + k₋₁*C])
@test set_decomposition!(red_i, [k₁*S*C - e₀*k₁*S + k₋₁*C])
@test set_decomposition!(red_i, manifolds[i][1])

# compute reduction
@test compute_reduction!(red_i)
@test red_i.g[1] == (-k₂*k₁*e₀*S*(k₁*S + k₋₁)//(k₁*k₋₁*e₀ + (k₁*S + k₋₁)^2))

# convenience function for reduction
red_i_conv = Reduction(prob, tfpvs[i], manifolds[i][1], parent(S//C).([S, e₀*k₁*S//(k₁*S + k₋₁)]))
@test all([getfield(red_i_conv, name) == getfield(red_i, name) for name in fieldnames(Reduction)])

# find all possible slow manifolds
unique_manifolds = unique_slow_manifolds(prob, manifolds)

# explicit description of manifolds
F = parent(S//C)
M = [
  F.([S, 0]),
  F.([-k₋₁//k₁, C]),
  F.([S, e₀*k₁*S//(k₁*S + k₋₁)])
]

R, idx = compute_all_reductions(prob, tfpvs, manifolds, M);

# Access the `Reduction` object 
reduction_3 = R[3][1]
@test reduction_3.M == [S, 0]
@test Matrix(reduction_3.Df0) == [-e₀*k₁+k₁*C k₁*S+k₋₁; e₀*k₁-k₁*C -k₁*S-k₋₁-k₂]
@test Matrix(reduction_3.Df0_at_x0) == [0 k₁*S+k₋₁; 0 -k₁*S-k₋₁-k₂]
@test all(Matrix(reduction_3.P) .== [S*k₁ + k₋₁; -S*k₁ - k₋₁ - k₂])
@test all(Matrix(reduction_3.Psi) .== [C])

# and the corresponding reduced system
@test reduction_3.g[1] == (-S*e₀*k₁*k₂)//(S*k₁ + k₋₁ + k₂)
