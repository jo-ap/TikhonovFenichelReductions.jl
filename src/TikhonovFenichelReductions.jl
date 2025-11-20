# TikhonovFenichelReductions.jl is licensed under the GPL v3+ (see LICENSE).
#
# Johannes Apelt
# contact: johannes.apelt@posteo.de

module TikhonovFenichelReductions

export ReductionProblem, 
  system_components, system_parameters,
  tfpvs_and_varieties, 
  tfpvs_groebner,
  Reduction,
  Variety,
  get_varieties,
  print_tfpvs, print_varieties, print_results, print_slow_fast, 
  print_system, print_reduced_system,
  rewrite_rational,
  jacobian_tfpv_on_manifold, jacobian_tfpv_at_x0,
  set_manifold!, set_point!, set_decomposition!, 
  compute_reduction!, 
  # get_reduced_system,
  # compute_directional_reduction,
  unique_varieties, compute_reductions, find_varieties,
  get_explicit_manifold

## Dependencies
import Base 
import Pkg
using DocStringExtensions
using LaTeXStrings
using Latexify
using Logging
using Oscar
using Oscar.AbstractAlgebra.Generic: FracField, FracFieldElem, MatSpaceElem, MPoly, PolyRing, Poly, RationalFunctionField, RationalFunctionFieldElem

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
