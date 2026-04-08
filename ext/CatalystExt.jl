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
  s::Int)
  # generate input for ReductionProblem
  # Catalyst <16.0
  # ode_system = convert(ODESystem, reaction_network)
  # func = generate_function(ode_system)[1]
  # Catalyst >=16.0
  ode_system = ode_model(reaction_network)
  func = generate_rhs(ode_system)[1]
  _f = Meta.eval(func);
  f(x,p) = Base.invokelatest(_f, x, p, 0)
  x = [replace(str, "(t)" => "") for str in string.(species(reaction_network))]
  p = string.(parameters(reaction_network))
  return ReductionProblem(f, x, p, s)
end

end
