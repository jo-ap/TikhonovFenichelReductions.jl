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
using TikhonovFenichelReductions
using TikhonovFenichelReductions.Oscar # optional
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
This returns a vector of boolean indices for each TFPV candidate, where ``0``
corresponds to a small and ``1`` to a large parameter.
In addition, we obtain the generators of the irreducible components of the
affine variety 
```math
\mathcal{V}_{\mathbb{C}}(f(\cdot,\pi^\star)) = \{x\in\mathbb{C}^n \mid \forall i : f_i(x,\pi^\star)=0\}
```
and their dimension, where ``\pi^\star`` is defined by the corresponding
slow-fast separation in `sf_separation`. 
Note that later have to check manually whether the variety taken in
``\mathbb{R}^n`` has the same dimension (i.e. if there exists a real
non-singular point).

```@example 1
sf_separation, V, dim_V = tfpv_candidates(problem)
```

You can also get all general TFPVs by computing a Gröbner basis `G` that
reflects necessary conditions on the parameters of the system. 
Note that this is potentially a very computationally intensive task.
```@example 1
G = tfpv_candidates_groebner(problem)
```

Show the results: All possible slow-fast separation of rates 
```@example 1
print_tfpv(problem, sf_separation)
```
and the corresponding irreducible components containing the slow manifold and their
Krull dimension
```@example 1
print_varieties(V, dim_V)
```

For further use, we need to make the variables (the system components and
parameters) available in the Main namespace.
Then, we can compute the reduction corresponding to TFPV 15, which is the
Rosenzweig-MacArthur model.
```@example 1
B, S, H = system_components(problem)
α, β, γ, δ, η, ρ = system_parameters(problem)

# instantiate reduction 
reduction = Reduction(problem, sf_separation[15])
```

To find the slow manifold, we can consider the affine variety defined by the vanishing of
```@example 1
V[15]
```
which has (complex) dimension
```@example 1
dim_V[15] 
```
as a variety. 
We can see, that there is only one irreducible component, so the slow manifold is 
``\mathcal{V}_{\mathbb{R}}(H) = \{(B,S,0) \mid B,S \in \mathbb{R}\}``
with dimension ``s=2`` as desired.

We need to define the slow manifold explicitly in order to check whether the
reduction exists. 
This also allows us to substitute the variables that got reduced according to
the slow manifold in the reduced system.
```@example 1
set_manifold!(reduction, [B, S, 0])
```

Note that in this case any generic point on the affine variety is non-singular
and can be chosen.  
Thus, we don't have to call `set_point!` to set a non-singular point explicitly
on whose neighbourhood the reduction exists. 
This also means ``\dim \mathcal{V}_{\mathbb{C}}(f(\cdot,\pi^\star)) = \dim
\mathcal{V}_{\mathbb{R}}(f(\cdot,\pi^\star))`` locally for the irreducible
component containing the non-singular point.

Lastly, we define a product decomposition ``f(\cdot,\pi^\star) = P \cdot \psi``
with ``\mathcal{V}(\psi) = \mathcal{V}(f(\cdot,\pi^\star))``.
This can be done by specifying only ``\psi`` using the generators of the
irreducible component of ``\mathcal{V}(f^{(0)})`` that corresponds to the slow
manifold.
Then, the methods automatically computes ``P``.
```@example 1
set_decomposition!(reduction, V[15][1])
```

Now we can compute and show the reduced system.
```@example 1
compute_reduction!(reduction)
print_reduced_system(reduction)
```
This updates the reduction object, which now contains the reduced system before
and after variables are substituted as defined by the slow manifold.
You can see how the reduction object is updated:
```@example 1 
reduction
```
The first two components of ```reduction.g``` define the reduced system and can
be rewritten as 
```math
\begin{align*}
\frac{dB}{dt} &= \rho B (1 - B) - \alpha(\eta + \beta) S \frac{B}{\delta + \gamma B} \\
\frac{dS}{dt} &= -\eta S + \gamma(\eta + \beta) S \frac{B}{\delta + \gamma B} \\
\end{align*}
```
which is exactly the Rosenzweig-MacArthur model.

The slow manifold is attractive if all non-zero eigenvalues of the Jacobian of
``f(\cdot, \pi^\star)`` at the non-singular `x0` point have negative real part.
```@example 1 
# pretty-print Jacobian at x0
J = jacobian_tfpv_at_x0(reduction)
show(stdout, "text/plain", J)
```
Thus, the full system converges to the reduction as ``\varepsilon \to 0`` if 
``B \geq 0``. 

