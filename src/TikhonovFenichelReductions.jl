module TikhonovFenichelReductions

export ReductionProblem, 
  num2bin,
  filter_determinants, filter_dimension, tfpv_candidates, 
  print_candidates, print_varieties,
  Reduction,
  jacobian_tfpv, jacobian_tfpv_on_manifold,
  set_manifold!, set_point!, set_decomposition!, 
  compute_reduction

# greetings
function __init__()
  if isinteractive() && Base.JLOptions().banner != 0
    println("\n    ==== TikhonovFenichelReductions.jl ==== \n")
    println("Compute Tikhonov-Fenichel reductions for polynomial ODE systems.")
    println("This package relies on Oscar.jl.\n")
  end
end

## Dependencies
using Oscar
using LaTeXStrings
using Latexify

## load code
include("./Finding_TFPV.jl")
include("./Computing_TFPV.jl")
include("./Output.jl")

end # module TikhonovFenichelReductions
