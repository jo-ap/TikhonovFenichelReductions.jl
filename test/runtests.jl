using Aqua
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
sf_separations, V, dim_V = tfpv_candidates(prob);

@testset "TikhonovFenichelReductions Tests" begin

  @testset "Finding TFPVs" begin
    include("tfpv.jl") 
  end

  @testset "Computing Reductions" begin
    include("reduction.jl") 
  end

  @testset "CatalystExt" begin
    include("catalyst.jl")
  end

  @testset "Output" begin
    include("output.jl")
  end

  @testset "Aqua" begin
    Aqua.test_all(TikhonovFenichelReductions)
  end

end
