# TikhonovFenichelReductions.jl is licensed under the GPL v3+ (see LICENSE).
#
# Johannes Apelt
# contact: johannes.apelt@posteo.de

module TikhonovFenichelReductions

export ReductionProblem, 
  system_components, system_parameters,
  tfpv_candidates, 
  tfpv_candidates_groebner,
  Reduction,
  print_tfpv, print_varieties, print_results, print_slow_fast, 
  print_system, print_reduced_system,
  rewrite_rational,
  slow_manifolds,
  jacobian_tfpv_on_manifold, jacobian_tfpv_at_x0,
  set_manifold!, set_point!, set_decomposition!, 
  compute_reduction!, 
  # get_reduced_system,
  # compute_directional_reduction,
  similar_reductions, unique_slow_manifolds, compute_bulk_reductions

## Dependencies
import Base 
import Pkg
using DocStringExtensions
using Oscar
using LaTeXStrings
using Latexify
using Oscar.AbstractAlgebra.Generic: FracField, FracFieldElem, MatSpaceElem, MPoly, PolyRing, Poly, RationalFunctionFieldElem

const PROJECT_TOML = Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))
const VERSION_NUMBER = VersionNumber(PROJECT_TOML["version"])

# greetings
function __init__()
  if isinteractive() && Base.JLOptions().banner != 0
    println("\n==== TikhonovFenichelReductions.jl v$(VERSION_NUMBER) ====\n")
    println("Compute Tikhonov-Fenichel reductions for polynomial ODE systems using Oscar.jl.\n")
  end
end

## load code
include("./tfpv.jl")
include("./reduction.jl")
include("./bulk_reductions.jl")
include("./display.jl")
include("./output.jl")


end # module TikhonovFenichelReductions
