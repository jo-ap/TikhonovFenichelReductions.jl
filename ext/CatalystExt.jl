module CatalystExt

using Catalyst
using DocStringExtensions
using TikhonovFenichelReductions 
using Oscar 

"""
    $(TYPEDSIGNATURES)

Constructor for `ReductionProblem` Type from a `ReactionSystem` defined with
`Catalyst.jl`.
"""
function TikhonovFenichelReductions.ReductionProblem(
  reaction_network::ReactionSystem,
  s::Int;
  idx_slow_fast::Vector{Bool}=[true for _ in parameters(reaction_network)])
  # generate input for ReductionProblem
  ode_system = convert(ODESystem, reaction_network)
  _f = Meta.eval(generate_function(ode_system)[1]);
  f(x,p) = Base.invokelatest(_f, x, p, 0)
  x = [replace(str, "(t)" => "") for str in string.(species(reaction_network))]
  p = string.(parameters(reaction_network))
  return ReductionProblem(f, x, p, s; idx_slow_fast=idx_slow_fast)
end

end
