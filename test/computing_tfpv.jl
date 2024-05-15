# Reduction instance 
red_i = Reduction(prob, idx[i])
@test isa(red_i, Reduction)

# make variables available
S, C = prob.x
e₀, k₁, k₋₁, k₂ = prob.θ

# test splitting of system
@test red_i.f⁰ == [-(k₁*e₀*S - (k₁*S + k₋₁)*C), k₁*e₀*S - (k₁*S + k₋₁)*C]
@test red_i.f¹ == [0, -k₂*C]

# check computation of reduced system
@test set_manifold!(red_i, [S, e₀*k₁*S//(k₁*S + k₋₁)])
@test set_point!(red_i, [0,0])
@test set_decomposition!(red_i, [k₁*S*C - e₀*k₁*S + k₋₁*C])
@test set_decomposition!(red_i, [1; -1], [k₁*S*C - e₀*k₁*S + k₋₁*C])

_, g_i = compute_reduction(red_i)
@test g_i[1] == (-k₂*k₁*e₀*S*(k₁*S + k₋₁)//(k₁*k₋₁*e₀ + (k₁*S + k₋₁)^2))
