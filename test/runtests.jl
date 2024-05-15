using Test
using TikhonovFenichelReductions

# Michaelis Menten kinetics as test case
x = ["S", "C"]
p = ["e₀", "k₁", "k₋₁", "k₂"]
function f(x, p)
  S, C = x
  e₀, k₁, k₋₁, k₂ = p
  return [
    -k₁*e₀*S + (k₁*S + k₋₁)*C,
    k₁*e₀*S - (k₁*S + k₋₁ + k₂)*C
  ]
end
prob = ReductionProblem(f, x, p, 1; idx_slow_fast = [true, true, false, true])
idx, G, (V, dim_Y) = tfpv_candidates(prob);

@testset "Finding TFPVs" begin
  include("finding_tfpv.jl") 
end

@testset "Computing TFPVs" begin
  include("computing_tfpv.jl") 
end
