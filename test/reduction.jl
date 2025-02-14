# Reduction instance 
red_i = Reduction(prob, sf_separations[i])
@test isa(red_i, Reduction)

# make variables available
S, C = system_components(prob)
e₀, k₁, k₋₁, k₂ = system_parameters(prob)

# test splitting of system
@test red_i.f0 == [-(k₁*e₀*S - (k₁*S + k₋₁)*C), k₁*e₀*S - (k₁*S + k₋₁)*C]
@test red_i.f1 == [0, -k₂*C]

# set everything that is requered for computation of reduced system
@test set_manifold!(red_i, [S, e₀*k₁*S//(k₁*S + k₋₁)])
@test set_point!(red_i, [0,0])
@test set_decomposition!(red_i, [k₁*S*C - e₀*k₁*S + k₋₁*C])
@test set_decomposition!(red_i, V[i][1])
@test set_decomposition!(red_i, [1; -1], [k₁*S*C - e₀*k₁*S + k₋₁*C])

# compute reduction
_, g_i = compute_reduction(red_i)
@test g_i[1] == (-k₂*k₁*e₀*S*(k₁*S + k₋₁)//(k₁*k₋₁*e₀ + (k₁*S + k₋₁)^2))

# bulk methods
V_unique = unique_slow_manifolds(prob, V, dim_V)
idx_similar = similar_reductions(V, V_unique[1])

R = compute_bulk_reductions(prob, sf_separations, idx_similar, [C], [S,0]; idx_components=[1]);

# Access the `Reduction` object with the indices as in `idx_similar`
reduction_3 = R[3][1]
@test reduction_3.M == [S, 0]
@test Matrix(reduction_3.Df0) == [-e₀*k₁+k₁*C k₁*S+k₋₁; e₀*k₁-k₁*C -k₁*S-k₋₁-k₂]
@test Matrix(reduction_3.Df0_at_x0) == [0 k₁*S+k₋₁; 0 -k₁*S-k₋₁-k₂]
@test all(Matrix(reduction_3.P) .== [S*k₁ + k₋₁; -S*k₁ - k₋₁ - k₂])
@test all(Matrix(reduction_3.Psi) .== [C])

# and the corresponding reduced system
g_3 = R[3][2]
@test g_3 == [(-S*e₀*k₁*k₂)//(S*k₁ + k₋₁ + k₂), 0]
