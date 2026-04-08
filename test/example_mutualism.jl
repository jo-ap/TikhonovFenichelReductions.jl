# Testset based on: 
# J. Apelt and V. Liebscher, ‘Tikhonov-Fenichel reductions and their
# application to a novel modelling approach for mutualism’, Theoretical
# Population Biology, vol. 166, pp. 16–35, Sept. 2025, 
# doi: 10.1016/j.tpb.2025.08.004.


# Define ODE system
# here we substitute kᵢ := 1/Kᵢ
x = ["H","S","C"]
p = ["β₂","β₃","δ₁","δ₂","δ₃","μ₁","μ₂", "η", "k₁", "k₂", "k₃"]
function f(x, p)
  H, S, C = x
  β₂, β₃, δ₁, δ₂, δ₃, μ₁, μ₂, η, k₁, k₂, k₃ = p
  return [
    -δ₁*H - η*S*H + μ₁*C*(1-H*k₁),
    β₂*S*(1-S*k₂) - δ₂*S - η*S*H  + μ₂*C*(1-S*k₂),
    β₃*C*(1-C*k₃) - δ₃*C + η*S*H
  ]
end

# Find TFPVs
s = 2
problem = ReductionProblem(f, x, p, s)
tfpvs, varieties = tfpvs_and_varieties(problem; preset = (k₁ = 1, k₂ = 1, k₃ = 1));
@test length(tfpvs) == length(varieties) == 27

# get all unique varieties
uni_varieties = unique_varieties(problem, varieties)
H, S, C = system_components(problem)
β₂, β₃, δ₁, δ₂, δ₃, μ₁, μ₂, η, k₁, k₂, k₃ = system_parameters(problem)
V_ideals = [
  ideal([S]),
  ideal([H]),
  ideal([k₂*S - 1]),
  ideal([C]),
  ideal([k₁*H - 1]),
  ideal([δ₂*S + μ₂*k₂*S*C - μ₂*C]),
  ideal([δ₁*H + μ₁*k₁*H*C - μ₁*C]),
  ideal([k₃*C - 1]),
  ideal([β₃*k₃*C - β₃ + δ₃]),
  ideal([β₂*S + μ₂*C]),
  ideal([β₂*k₂*S - β₂ + δ₂]),
  ideal([β₂*k₂*S^2 - β₂*S + δ₂*S + μ₂*k₂*S*C - μ₂*C])
]
@test all([ideal(V.gens_R) == I for (V,I) in zip(uni_varieties, V_ideals)])

# test if manifolds can be computed automatically from varieties
M_auto = [get_explicit_manifold(problem, V) for V in uni_varieties]
@test all([m[2] for m in M_auto])

# manually set one manifold
manifolds = [m[1] for m in M_auto];
manifolds[7] = [H//QQ(1), S//QQ(1), δ₁*H//(μ₁*(1 - k₁*H))]

# compute all reductions 
reductions, idx = compute_reductions(problem, tfpvs, varieties, uni_varieties, manifolds);
@test length(reductions) == 35
@test length(idx) == length(uni_varieties) == length(manifolds) == 12
_idx = [
  [(1, 1), (8, 1), (9, 1), (22, 1), (23, 1), (25, 1), (26, 1)],
  [(1, 2), (11, 1), (12, 1)],
  [(2, 1), (22, 2), (24, 1)],
  [(2, 2), (3, 2), (4, 1), (5, 1), (6, 1), (7, 1), (14, 1), (15, 1), (16, 1), (17, 1), (18, 1), (19, 1), (20, 1), (21, 1)],
  [(3, 1)],
  [(10, 1)],
  [(13, 1)],
  [(14, 2)],
  [(18, 2)],
  [(24, 2)],
  [(25, 2)],
  [(27, 1)]
]
@test idx == _idx
g = [
	β₂*S*(1-S*k₂) - δ₂*S + μ₂*C*(1-S*k₂) - μ₁*C*(η*S)//(η*S + δ₁),
	β₃*C*(1-C*k₃) - δ₃*C + μ₁*C*(η*S)//(η*S + δ₁)
]
@test all(reductions[(12,1)].g[reductions[(12,1)].idx_components] .- g .== 0)
