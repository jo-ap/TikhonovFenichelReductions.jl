# TikhonovFenichelReductions.jl is licensed under the GPL v3+ (see LICENSE).
#
# Johannes Apelt
# contact: johannes.apelt@posteo.de

module TikhonovFenichelReductions

export ReductionProblem, 
  system_components, system_parameters,
  tfpv_candidates, 
  tfpv_groebner_basis,
  print_tfpv, print_varieties, print_results,
  Reduction,
  show_slow_fast,
  slow_manifolds,
  jacobian_tfpv_on_manifold, jacobian_tfpv_at_x0,
  set_manifold!, set_point!, set_decomposition!, 
  compute_reduction, compute_directional_reduction,
  similar_reductions, unique_slow_manifolds, compute_bulk_reductions

# greetings
function __init__()
  if isinteractive() && Base.JLOptions().banner != 0
    println("\n==== TikhonovFenichelReductions.jl ====\n")
    println("Compute Tikhonov-Fenichel reductions for polynomial ODE systems.")
    println("This package relies on Oscar.jl.\n")
  end
end

## Dependencies
using DocStringExtensions
using Oscar
using LaTeXStrings
using Latexify

using AbstractAlgebra.Generic: FracField, FracFieldElem, MatSpaceElem, MPoly, PolyRing, Poly, RationalFunctionFieldElem

## load code
include("./tfpv.jl")
include("./reductions.jl")
include("./bulk_reductions.jl")
include("./output.jl")

end # module TikhonovFenichelReductions
