# TikhonovFenichelReductions.jl is licensed under the GPL v3+ (see LICENSE).
#
# Johannes Apelt
# contact: johannes.apelt@posteo.de

module TikhonovFenichelReductions

export ReductionProblem, 
  tfpv_candidates, 
  tfpv_groebner_basis,
  print_tfpv, print_varieties, print_results,
  Reduction,
  get_slow_manifolds,
  jacobian_tfpv, jacobian_tfpv_on_manifold,
  set_manifold!, set_point!, set_decomposition!, 
  compute_reduction, compute_directional_reduction

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

## load code
include("./Finding_TFPV.jl")
include("./Computing_TFPV.jl")
include("./Output.jl")

end # module TikhonovFenichelReductions
