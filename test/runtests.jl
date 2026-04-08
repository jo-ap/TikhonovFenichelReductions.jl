using Aqua
using Oscar
using Test
using TikhonovFenichelReductions

# Test basic functionality with Michaelis-Menten kinetics example
x = ["S", "C"]
p = ["e₀", "k₁", "k₋₁", "k₂"]
function f(x, p)
  S, C = x
  e₀, k₁, k₋₁, k₂ = p
  return [-k₁*e₀*S + (k₁*S + k₋₁)*C, k₁*e₀*S - (k₁*S + k₋₁ + k₂)*C]
end
problem = ReductionProblem(f, x, p, 1)
@testset verbose=true "TikhonovFenichelReductions Basic Tests" begin
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
  @testset "Edge Cases" begin
    include("edge_cases.jl")
  end
end

# Complex examples as additional cases
@testset verbose=true "TikhonovFenichelReductions Modelling Examples" begin
  @testset "Rosenzweig-MacArthur Model" begin
    include("example_rosenzweig-macarthur.jl")
  end
  @testset "Mutualism Model" begin
    include("example_mutualism.jl")
  end
end


