
using Oscar # optional (results in prettier printing and loads Oscar functions to Main namespace)
using Debugger
using Revise
using TikhonovFenichelReductions

## Define system

# dynamic variables
x = ["B", "S", "H"]

# parameters
p = ["α", "β", "γ", "δ", "η", "ρ"]

# RHS of ODE system ẋ = f(x, p)
function f(x, p)
  B, S, H = x
  α, β, γ, δ, η, ρ = p
  return [
    ρ*B*(1-B) - α*B*H,
    -η*S + γ*B*H,
    β*S- δ*H + η*S - γ*B*H
  ]
end

# dimension of the reduced system
s = 2

# create Problem
problem = ReductionProblem(f, x, p, s)

# find slow-fast separations that are TFPVs
@time sf_separations, V, dim_V = tfpv_candidates(problem);

print_results(problem, sf_separations, V, dim_V)

B, S, H = x = system_components(problem)
α, β, γ, δ, η, ρ = p = system_parameters(problem)

# The Rosenzweig-MaxArthur system corresponds to the TFPV candidate 15 (See section 3.3 in the paper).
# instantiate reduction 
reduction = Reduction(problem, sf_separations[15])
set_manifold!(reduction, [B, S, 0])
set_decomposition!(reduction, V[15][1])
g,_ = compute_reduction(reduction)





function rewrite_rational(term)
  p = numerator(term)
  q = denominator(term)
  h,r = divrem(p,q)
  println(string(h) * " + (" * string(r) * ")//(" * string(q) * ")")
  return h, r, q
end

rewrite_rational.(g);
