# Testset based on: 
# N. Kruff, C. Lax, V. Liebscher, and S. Walcher, ‘The Rosenzweig–MacArthur
# system via reduction of an individual based model’, J. Math. Biol., vol. 78,
# no. 1–2, pp. 413–439, Jan. 2019, doi: 10.1007/s00285-018-1278-y.

# Define system
x = ["B", "S", "H"]
p = ["α", "β", "γ", "δ", "η", "ρ"]
function f(x, p)
  B, S, H = x
  α, β, γ, δ, η, ρ = p
  return [
    ρ*B*(1-B) - α*B*H,
    -η*S + γ*B*H,
    β*S- δ*H + η*S - γ*B*H
  ]
end

# init problem
s = 2
problem = ReductionProblem(f, x, p, s)

# Restrict search via preset rates
preset = (α=1, β=0)
sf_separations, varieties = tfpvs_and_varieties(problem; preset=preset);
@test length(sf_separations) == length(varieties) == 6

# Find all reductions
sf_separations, varieties = tfpvs_and_varieties(problem);
@test length(sf_separations) == length(varieties) == 15

# Make variables and parameters available in Main namespace
B, S, H = system_components(problem)
α, β, γ, δ, η, ρ = system_parameters(problem)

# The Rosenzweig-MacArthur system corresponds to the TFPV candidate 15 
reduction = Reduction(problem, sf_separations[15])

# compute reduction
@test set_manifold!(reduction, [B, S, 0])
@test set_decomposition!(reduction, varieties[15][1])
@test set_decomposition!(reduction, [H])
@test compute_reduction!(reduction)

# check if g as the RHS of the reduced system is the same as in the paper 
dBdt = ρ*B*(1-B) - α*(η + β)*B*S//(δ + γ*B)
dSdt = -η*S + γ*(η + β)*B*S//(δ + γ*B)
@test all(iszero.(reduction.g[reduction.idx_components] .- [dBdt, dSdt])) 

# convenience function to compute reduction directly
@test Reduction(problem, sf_separations[15], varieties[15][1], parent(B).([B,S,0])) == reduction

## Compute multiple reductions at once

# Get all unique slow manifolds of dimension 2 onto which reductions exist
all_varieties = unique_varieties(problem, varieties)
V_ideals = [
 ideal([B]),
 ideal([B-1]),
 ideal([S]),
 ideal([H]),
 ideal([γ*B*H - η*S]),
 ideal([β*S - δ*H]),
 ideal([ρ*B + α*H - ρ])
]
@test all([ideal(V.gens_R) == I for (V,I) in zip(all_varieties, V_ideals)])

# test if all manifolds can be computed automatically
M_auto = [get_explicit_manifold(problem, V) for V in all_varieties]
@test all([m[2] for m in M_auto])

F = parent(B//S)
M = [
  F.([0,S,H]),
  F.([1,S,H]),
  F.([B,0,H]),
  F.([B,S,0]),
  F.([η*S//(γ*H),S,H]),
  F.([B,δ//β*H,H]),
  F.([1 - α//ρ*H,S,H])
]
@test all(M .== [m[1] for m in M_auto])

# compute multiple reductions
R, idx = compute_reductions(problem, sf_separations, varieties, all_varieties, [M[1] for M in M_auto])
@test R[(15,1)] == reduction

