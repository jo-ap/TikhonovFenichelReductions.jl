
using Oscar # optional (results in prettier printing and loads Oscar functions to Main namespace)
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

using Debugger 
@enter Reduction(problem, sf_separations[15])

Ψ = V[15][1]
typeof(Ψ[1])

sf_separation = sf_separations[15]

R = parent(problem.f[1])
K = fraction_field(R)
_p = TikhonovFenichelReductions.get_tfpv(problem, sf_separation)
n = length(problem.x)
r = n - s
f0, f1 = TikhonovFenichelReductions.splitsystem(problem, sf_separation)
Df = jacobian(problem.f, problem.x)
T, _ = polynomial_ring(K, "λ")
M = K.(problem.x)
x0 = zeros(K, n)
P = zero_matrix(K,n,r)

P .= f0.//Ψ
DPsi = jacobian(Ψ, problem.x)

# compute reduced system
A = DPsi*P
Q = P*inv(A)*DPsi
Iₙ = diagonal_matrix(K(1), n)
f1 = matrix_space(K, n, 1)(f1)
f_red_raw = (Iₙ - Q)*f1
# reshape as vector
f_red = reshape(Matrix(f_red_raw), n)
a = K.(M)
  f_red_subs = [evaluate(f, a) for f in f_red]
  return f_red, f_red_subs
else
  @warn "Slow manifold has not been defined succesfully. Reduced system is only returned in raw form, i.e. the reduced components are not substituted according to the slow manfold."
  return f_red, nothing
end
 Df = jacobian(problem.f, problem.x)

x0 = K.([B, S, 0])
Df_x0 = TikhonovFenichelReductions.eval_mat(Df, x0)

  chi = charpoly(T, Df_x0)
Df_x0 = TikhonovFenichelReductions.eval_mat(Df, [x0; reduction.K.(reduction._p)])
  
# look at the variety that contains the slow manifold
V[15] # => M₀ = {(B,S,0) | B,S ∈ ℝ}
dim_V[15] # has dimension 2 
set_manifold!(reduction, [B, S, 0])

# define product decomposion f0 = P⋅Psi (can be done via specifying Psi with V(Psi) = V(f⁰) in this case)
set_decomposition!(reduction, [H])

# compute the reduced system
g_raw, g = compute_reduction(reduction);
display(g)

reduction = 
