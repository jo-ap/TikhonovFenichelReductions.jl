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
types and imports useful functions such as `factor`, that we might want to use
later in order to simplify the reduced systems.
```@example 1
using Oscar # optional
using TikhonovFenichelReductions
```
Define the components, parameters and the RHS of the system. 
```@example 1
# components
x = ["B", "S", "H"]

# parameters
θ = ["α", "β", "γ", "δ", "η", "ρ"]

# RHS of ODE system ẋ = f(x, θ), where f is polynomial in x and θ
function f(x, θ)
  B, S, H = x
  α, β, γ, δ, η, ρ = θ
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
prob = ReductionProblem(f, x, θ, s)
```
Compute all TFPV candidates by using necessary conditions for the existence of
a reduction (in this case the ideal defined by `G` can be generated by
monomials, which means these really are all TFPVs). 
```@example 1
idx, G, (V, dim_Y) = tfpv_candidates(prob)
```

Show the results: All possible slow-fast separation of rates
```@example 1
print_candidates(idx, prob)
```
the corresponding irreducible components containing the slow manifold
```@example 1
print_varieties(V, prob)
```
and their Krull dimension
```@example 1
dim_Y
```

For further use, we need to make the variables (the system components and
parameters) available in the Main namespace.
Here we compute the reduction corresponding to TFPV 16, which is the 
Rosenzweig-MacArthur model.
```@example 1
B, S, H = prob.x
α, β, γ, δ, η, ρ = prob.θ

# instantiate reduction 
reduction = Reduction(prob, idx[16])
```

To find the slow manifold, we can consider the affine variety.
```@example 1
V[16]
```
This means the slow manifold is 
``\mathcal{V}(f^{(0)}) = \mathcal{V}(H) = \{(B,S,0) \mid B,S \in \mathbb{R}\}``,
since there is only one component with dimension ``s=2`` as desired:
```@example 1
dim_Y[16] 
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
``\mathcal{V}(\psi) = \mathcal{V}(f^{(0)})`` in some neighbourhood. 
This can always be done by specifying only ``\psi`` if ``n-s=1``, as is the case
here:
```@example 1
set_decomposition!(reduction, [H])
```
Now we can compute the reduced system.
```@example 1
_, g = compute_reduction(reduction)
```
The first two components define the reduced system and can be rewritten as 
```math
\begin{align*}
\frac{dB}{dt} &= \rho B (1 - B) - \alpha(\eta + \beta) S \frac{B}{\delta + \gamma B} \\
\frac{dS}{dt} &= -\eta S + \gamma(\eta + \beta) S \frac{B}{\delta + \gamma B} \\
\end{align*}
```
which is exactly the Rosenzweig-MacArthur model.

The slow manifold is attractive if all non-zero eigenvalues of the Jacobian at
the non-singular `x₀` point have negative real part. 
```@example 1 
reduction.Df_x₀
```
Thus, the full system converges to the reduction as ``\varepsilon \to 0`` if 
``B > \frac{\delta}{\gamma}``. 

