# Examples

## Michaelis-Menten kinetics 
The Michaelis-Menten reaction scheme
```math
E + S \xrightleftharpoons[k_{-1}]{k_1} C \xrightleftharpoons[k_{-2}]{k_2} E + P
```
leads to the following ODE system (using the conserved quantities $E + C$ and
$S + C + P$)
```math
\begin{aligned}
    \frac{dS}{dt} &= -k_1 e_0 S + (k_1 S + k_{-1}) C \\ 
    \frac{dC}{dt} &= k_1 e_0 S - (k_1 S + k_{-1} + k_2) C + k_{-2} (e_0 - C) (s_0 - S - C)
\end{aligned}
```
where $e_0 = E(0), s_0 = S(0)$.

The classical one-dimensional equation for the Michaelis-Menten reaction
kinetics is typically derived using a quasi-steady state (QSS) assumption.
Following [goeke2013a](@cite), this approach can be formalized with
Tikhonov-Fenichel reduction theory.
Here, we see how this can be achieved using `TikhonovFenichelReductions.jl`.
In the following, we consider the irreversible and reversible case, i.e.
$k_{-2}=0$ and $k_{-2}>0$, respectively.

```@example 1
# Setup
using Oscar 
using TikhonovFenichelReductions

# Define the system
x = ["S", "C"]
p = ["e₀", "s₀", "k₁", "k₋₁", "k₂", "k₋₂"]
function f(x, p)
  S, C = x
  e₀, s₀, k₁, k₋₁, k₂, k₋₂ = p
  return [-k₁*e₀*S + (k₁*S + k₋₁)*C, k₁*e₀*S - (k₁*S + k₋₁ + k₂)*C + k₋₂*(e₀ - C)*(s₀ - S - C)]
end
``` 
### The irreversible case
We start by considering the irreversible case, i.e. $k_{-2} = 0$.
```@example 1
# RHS for the irreversible case
f1(x,p) = f(x, p .* [1,1,1,1,1,0])

# Initialize the problem
problem = ReductionProblem(f1, x, p, 1)
```
```@example 1
# Search for TFPVs (fix s₀ and k₋₂, since they are irrelevant)
tfpvs, varieties = tfpvs_and_varieties(problem; preset=(s₀ = 1, k₋₂ = 1));
print_tfpvs(problem, tfpvs)
```

```@example 1 
print_varieties(varieties)
```

We can find all unique manifolds and compute the corresponding reductions.
```@example 1 
unique_V = unique_varieties(problem, varieties)
M = [get_explicit_manifold(problem, V) for V in unique_V]
```

```@example 1 
@assert all([m[2] for m in M])
R, idx = compute_reductions(problem, tfpvs, varieties, unique_V, [m[1] for m in M])
```

The ad-hoc method to reduce this system found in most textbooks relies on the
assumption that (at least in the initial phase of the reaction)
$\frac{dC}{dt}\in\mathcal{O}(\varepsilon)$ for $\varepsilon>0$ small.
Thus, one approximates $\varepsilon=0$ and solves
```math 
0 = \frac{dC}{dt} = k_1 e_0 S - (k_1 S + k_{-1} + k_2) C,
``` 
for $C$, which yields 
```math 
C = \frac{k_1 e_0 S}{k_1 S + k_{-1} + k_2}
```
and thus 
```math 
\frac{dS}{dt} = -\frac{k_1 k_2 e_0 S}{k_1 S + k_{-1} + k_2}.
```
The same result can be obtained by assuming the initial enzyme concentration to
be small, i.e. $e_0\in\mathcal{O}(\varepsilon)$ while all other parameters are
in $\mathcal{O}(1)$. 
```@example 1 
print_reduced_system(R[(7,1)]; rewrite=false)
```
The difference is that with this approach the reduced system is derived as the
limit $\varepsilon \to 0$ instead of just setting $\varepsilon=0$.
Furthermore, with Tikhonov-Fenichel reduction theory we consider individual
processes independent from the components, whereas the ad-hoc approach implies
all processes constituting change in $C$ are slow even though most of these
processes also govern the development of $S$, which is not considered to change
slowly.

Thus, Tikhonov-Fenichel reduction theory is a mathematical sound approach to
reduce ODE systems, which already yields multiple timescale separations for the
given system.
However, for this example, the benefits of this approach become even more
apparent when we consider the reversible case.

### The reversible case
In the reversible case, i.e. $k_{-2}>0$, the QSS ad-hoc approach of assuming
$\frac{dC}{dt}=0$ and solving the resulting algebraic equation for $C$,
yields an expression involving a square root.
However, with Tikhonov-Fenichel theory, we can still obtain a reduction with
rational RHS. 
```@example 1
# Initialize the problem
problem = ReductionProblem(f, x, p, 1)
```
Assuming the small parameter is still $e_0$, we can obtain the following
reduction.
```@example 1 
tfpv = Bool[0,1,1,1,1,1]
V = get_varieties(problem, tfpv)
```
Only the first irreducible component of $\mathcal{V}(f^{(0)})$ has dimension 1
as desired.
```@example 1 
M, _ = get_explicit_manifold(problem, V[1])
reduction = Reduction(problem, tfpv, V[1], M)
print_reduced_system(reduction; rewrite=false)
```
As before, there are other possible choices for slow-fast separations of rates.
```@example 1 
tfpvs, varieties = tfpvs_and_varieties(problem)
print_tfpvs(problem, tfpvs)
```
```@example 1 
print_varieties(varieties)
```
And indeed the reversible case can be seen as a limit of the irreversible case
with $k_{-2}$ approaching zero.
```@example 1 
M, _ = get_explicit_manifold(problem, varieties[29][1])
reduction = Reduction(problem, tfpvs[29], varieties[29][1], M)
print_reduced_system(reduction; rewrite=false)
```
A thorough discussion of the differences between QSS and Tikhonov-Fenichel
theory can be found in [goeke2013a](@cite).

## A Mutualism model
A more involved example can be found in [apelt2025](@cite), where a model for
mutualism is investigated using Tikhonov-Fenichel reduction theory.
The main idea is to define the two-species system for the subpopulations of a
host $H$, a symbiont $S$ and their complex $C$, thereby accounting explicitly
for the state of individuals (interacting or living autarkically).
From this, one can then consider two-dimensional reduced systems.
```@example 2 
using Oscar 
using TikhonovFenichelReductions

# Define ODE system f
# here we substitute kᵢ := 1/Kᵢ
x = ["H","S","C"]
p = ["β₂","β₃","δ₁","δ₂","δ₃","μ₁","μ₂", "η", "k₁", "k₂", "k₃"]
function f(x, p)
  H, S, C = x
  β₂, β₃, δ₁, δ₂, δ₃, μ₁, μ₂, η, k₁, k₂, k₃ = p
  return [
    -δ₁*H - η*S*H + μ₁*C*(1-H*k₁),
    β₂*S*(1-S*k₂) - δ₂*S - η*S*H  + μ₂*C*(1-S*k₂),
    β₃*C*(1-C*k₃) - δ₃*C + η*S*H
  ]
end

# find TFPV candidates for dimension s=2 (fix carrying capacities)
problem = ReductionProblem(f, x, p, 2)
```
```@example 2
tfpvs, varieties = tfpvs_and_varieties(problem; preset = (k₁ = 1, k₂ = 1, k₃ = 1));

# make variables available in Main namespace
H, S, C = system_components(problem)
β₂, β₃, δ₁, δ₂, δ₃, μ₁, μ₂, η, k₁, k₂, k₃ = system_parameters(problem)

# compute all reductions 
unique_V = unique_varieties(problem, varieties)

M_auto = [get_explicit_manifold(problem, V) for V in unique_V]
@assert all([m[2] for m in M_auto])
manifolds = [m[1] for m in M_auto];
F = parent(H//S)
manifolds[7] = F.([H, S, δ₁*H//(μ₁*(1 - k₁*H))])
reductions, idx_M = compute_reductions(problem, tfpvs, varieties, unique_V, manifolds)
```
The system given in Eq. 12 in [apelt2025](@cite) is for instance
```@example 2 
reductions[(12,1)]
```
By reducing the original ODE system to the plane, we are able to investigate
the behaviour of the system analytically. 
In particular, this allows to get a better understanding of conditions for
which the mutualistic relation may collapse.
