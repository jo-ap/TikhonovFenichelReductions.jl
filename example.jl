## Exemplary use of TFPV.jl 
#
# This example is based on:
# N. Kruff, C. Lax, V. Liebscher, and S. Walcher, ‘The Rosenzweig–MacArthur
# system via reduction of an individual based model’, J. Math. Biol., vol. 78,
# no. 1–2, pp. 413–439, Jan. 2019, doi: 10.1007/s00285-018-1278-y.

## Load package
using Oscar
using TikhonovFenichelReductions

## Define system
# dynamic variables
x = ["B", "S", "H"]

# parameters
θ = ["α", "β", "γ", "δ", "η", "ρ"]

# RHS of ODE system ẋ = f(x, θ)
function f(x, θ)
  B, S, H = x
  α, β, γ, δ, η, ρ = θ
  return [
    ρ*B*(1-B) - α*B*H,
    -η*S + γ*B*H,
    β*S- δ*H + η*S - γ*B*H
  ]
end

# dimension of the reduced system
s = 2

# create Problem
prob = ReductionProblem(f, x, θ, s)

# compute TFPV candidates
# if this takes too long, consider setting compute_primary_decomposition=false
idx, G, V = tfpv_candidates(prob);

# Output
print_candidates(idx, prob; latex=true) # Candidates for TFPVs
print_varieties(V, prob) # The irreducible components of V(f⁰)

# change order of filters 
idx_dim, _V = filter_dimension(prob)
idx_det, _G = filter_determinants(prob; idx=idx_dim)

# check if this is the same
all(idx_det .== idx)

## Make variables and parameters available in Main namespace

B, S, H = prob.x;
α, β, γ, δ, η, ρ = prob.θ;

## Compute a reduced system
# The Rosenzweig-MaxArthur system corresponds to the TFPV candidate π₁₆ (See section 3.3 in the paper).

# instantiate reduction 
reduction = Reduction(prob, idx[16]);

# look at the variety that contains the slow manifold
V[16] # ⟹ M₀ = {(B,S,0) | B,S ∈ ℝ}
set_manifold!(reduction, [B, S, 0])

# choose non-singular point
set_point!(reduction, [0,0,0])

# define product decomposion f⁰ = P⋅ψ (can be done via specifying ψ with V(ψ) = V(f⁰) in this case)
set_decomposition!(reduction, [H])

# compute the reduced system
_, g = compute_reduction(reduction)

g[1]
g[2]

# check if g as the RHS of the reduced system is the same as in the paper 
dBdt = ρ*B*(1-B) - α*(η + β)*B*S//(δ + γ*B)
dSdt = -η*S + γ*(η + β)*B*S//(δ + γ*B)
all(iszero.(g[1:2] .- [dBdt, dSdt])) # yes

