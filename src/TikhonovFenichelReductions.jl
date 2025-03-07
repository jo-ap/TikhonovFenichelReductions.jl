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
  print_tfpv, print_varieties, print_results, print_slow_fast, print_system, print_reduced_system,
  slow_manifolds,
  jacobian_tfpv_on_manifold, jacobian_tfpv_at_x0,
  set_manifold!, set_point!, set_decomposition!, 
  compute_reduction!, 
  # get_reduced_system,
  # compute_directional_reduction,
  similar_reductions, unique_slow_manifolds, compute_bulk_reductions,
  print_reduced_system

## Dependencies
using DocStringExtensions
using OrderedCollections
using Oscar
using LaTeXStrings
using Latexify

using AbstractAlgebra.Generic: FracField, FracFieldElem, MatSpaceElem, MPoly, PolyRing, Poly, RationalFunctionFieldElem

import Pkg
TFR_VERSION = Pkg.project().version

# greetings
function __init__()
  if isinteractive() && Base.JLOptions().banner != 0
    println("\n==== TikhonovFenichelReductions.jl v$(TFR_VERSION) ====\n")
    println("Compute Tikhonov-Fenichel reductions for polynomial ODE systems.")
    println("This package relies on Oscar.jl.\n")
  end
end

## load code
include("./tfpv.jl")
include("./reduction.jl")
include("./bulk_reductions.jl")
include("./display.jl")
include("./output.jl")

end # module TikhonovFenichelReductions
