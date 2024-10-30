# Getting Started

## Installation
Run
~~~
add https://github.com/jo-ap/tikhonovfenichelreductions.jl
~~~
in Julia package Mode.

## Example
Here we consider the derivation of the Rosenzweig-MacArthur model as a
reduction from a three dimensional system as demonstrated in
[kruff2019](@cite). 

Load the package and its dependency [Oscar.jl](https://www.oscar-system.org/). 
Note that loading Oscar is optional, but results in prettier printing of
types and imports useful functions.
```@example 1
using Oscar # optional
using TikhonovFenichelReductions
```
Define the components, parameters and the RHS of the system. 
```@example 1
# components
x = ["B", "S", "H"]

# parameters
p = ["α", "β", "γ", "δ", "η", "ρ"]

# RHS of ODE system ẋ = f(x, p), where f is polynomial in x and p
function f(x, p)
  B, S, H = x
  α, β, γ, δ, η, ρ = p
  return [
    ρ*B*(1-B) - α*B*H,
    -η*S + γ*B*H,
    β*S- δ*H + η*S - γ*B*H
  ]
end
```

Initialize the problem with desired dimension of the reduced system
```@example 1
# dimension of the reduced system
s = 2

# create problem
problem = ReductionProblem(f, x, p, s)
```
Compute all slow-fast separations that are TFPVs by using necessary conditions
for the existence of a reduction. 
```@example 1
idx, V, dim_V = tfpv_candidates(problem)
```

You can also get all general TFPVs by computing a Gröbner basis `G` that
reflects necessary conditions on the parameters of the system. 
Note that this is potentially a very computationally intensive task.
```@example 1
G = tfpv_groebner_basis(problem)
```

Show the results: All possible slow-fast separation of rates
```@example 1
print_tfpv(problem, idx)
```
the corresponding irreducible components containing the slow manifold and their
Krull dimension
```@example 1
print_varieties(V; dim_V=dim_V)
```

For further use, we need to make the variables (the system components and
parameters) available in the Main namespace.
Here we compute the reduction corresponding to TFPV 15, which is the 
Rosenzweig-MacArthur model.
```@example 1
B, S, H = problem.x
α, β, γ, δ, η, ρ = problem.p

# instantiate reduction 
reduction = Reduction(problem, idx[15])
```

To find the slow manifold, we can consider the affine variety.
```@example 1
V[15]
```
This means the slow manifold is 
``\mathcal{V}(f^{(0)}) = \mathcal{V}(H) = \{(B,S,0) \mid B,S \in \mathbb{R}\}``
with dimension ``s=2`` as desired:
```@example 1
dim_V[15] 
```
Therefore, we can set the slow manifold accordingly.
```@example 1
set_manifold!(reduction, [B, S, 0])
```

Note that in this case any generic point on the affine variety is non-singular
and can be chosen.  
Thus, we don't have to call `set_point!` to set a non-singular point explicitly
on whose neighbourhood the reduction exists. 

Lastly, we define a product decomposition ``f^{(0)} = P \cdot \psi`` with 
``\mathcal{V}(\psi) = \mathcal{V}(f^{(0)})``.
If ``n-s=1``, this can always be done by specifying only ``\psi``:
here:
```@example 1
set_decomposition!(reduction, [H])
```
Now we can compute the reduced system.
```@example 1
_, g = compute_reduction(reduction)
```
The first return value is the system before variables are substituted as
defined by the slow manifold.
The first two components ```g``` define the reduced system and can be rewritten as 
```math
\begin{align*}
\frac{dB}{dt} &= \rho B (1 - B) - \alpha(\eta + \beta) S \frac{B}{\delta + \gamma B} \\
\frac{dS}{dt} &= -\eta S + \gamma(\eta + \beta) S \frac{B}{\delta + \gamma B} \\
\end{align*}
```
which is exactly the Rosenzweig-MacArthur model.

The slow manifold is attractive if all non-zero eigenvalues of the Jacobian at
the non-singular `x0` point have negative real part. 
```@example 1 
# pretty-print Jacobian at x0
show(stdout, "text/plain", reduction.Df_x0)
```
Thus, the full system converges to the reduction as ``\varepsilon \to 0`` if 
``B \geq 0``. 

