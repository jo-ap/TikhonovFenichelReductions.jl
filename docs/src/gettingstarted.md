# Getting Started

## Installation
Run
~~~
add https://github.com/jo-ap/tikhonovfenichelreductions.jl
~~~
in Julia package Mode.
!!! note "Installation"
    TikhonovFenichelReductions.jl relies on Oscar.jl, which requieres at least
    6GB of free memory for the installation (the Oscar-Team recommends at least
    16GB).
    On Windows, Oscar.jl needs to be installed using the Windows Subsystem for
    Linux (WSL).
    Instructions can be found in their
    [documentation](https://www.oscar-system.org/install/win/).
    


## Example
Here we consider the derivation of the Rosenzweig-MacArthur model as a
reduction from a three dimensional system as demonstrated in
[kruff2019](@cite). 

Load the package and its dependency [Oscar.jl](https://www.oscar-system.org/). 
Note that loading Oscar is optional, but results in prettier printing of
types and imports useful functions.
```@example 1
using TikhonovFenichelReductions
using Oscar # optional
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
### Slow-Fast Separations of Rates
We can find all slow-fast separations that are TFPVs by using necessary
conditions for the existence of a reduction. 
This returns a vector of boolean indices for each TFPV candidate, where ``0``
corresponds to a small parameters, i.e. those multiplied by ``\varepsilon``,
and ``1`` corresponds to a parameter in ``\mathcal{O}(1)``.
In addition, we obtain the generators of the irreducible components of the
affine variety 
```math
\mathcal{V}_{\mathbb{C}}(f(\cdot,\pi^\star)) = \{x\in\mathbb{C}^n \mid f(x,\pi^\star)=0\}
```
and their dimension, where ``\pi^\star`` is defined by the corresponding
slow-fast separation in `sf_separation`. 
This information is stored using the type `Variety`.
Note that later have to check manually whether the variety taken in
``\mathbb{R}^n`` has the same dimension (i.e. if there exists a real
non-singular point), which then at least locally renders this variety a
manifold.

```@example 1
tfpvs, varieties = tfpvs_and_varieties(problem)
```
Show the results: All possible slow-fast separation of rates 
```@example 1
print_tfpvs(problem, tfpvs)
```
and the corresponding irreducible components containing the slow manifold and their
Krull dimension
```@example 1
print_varieties(varieties)
```

For further use, we need to make the variables (the system components and
parameters) available in the Main namespace.
Then, we can compute the reduction corresponding to TFPV 15, which is the
Rosenzweig-MacArthur model.
```@example 1
B, S, H = system_components(problem)
α, β, γ, δ, η, ρ = system_parameters(problem)

# instantiate reduction 
reduction = Reduction(problem, tfpvs[15])
```

To find the slow manifold, we can consider the affine variety defined by the vanishing of
```@example 1
varieties[15][1].gens_R
```
which has (complex) dimension
```@example 1
varieties[15][1].dim
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
Then, the methods automatically computes ``P``, which relies on having exactly
``r`` generators of the variety corresponding to the slow manifold.
```@example 1
set_decomposition!(reduction, varieties[15][1])
```
In most cases, this method should work, but you can also manually set ``P`` and
``Psi``.

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

#### Bulk Computations 
`TikhonovFenichelReductions.jl` has methods that simplify the computation of
multiple reductions with the same manifold, since typically many different TFPVs
share the same slow manifold.
Additionally, the function `get_explicit_manifold` implements a heuristic 
to find an explicit parametric description of the slow manifold.
Thus, in simple enough cases, the computation of all reductions is fully
automatic. 

```@example 1
all_V = unique_varieties(problem, varieties)
all_M = [get_explicit_manifold(problem, V) for V in all_V] 

# get reduction and indices of varieties that correspond to unique slow manifolds
R, idx_M = compute_reductions(problem, tfpvs, varieties, all_V, [m[1] for m in all_M]);

R[(15,1)]
```

### General TFPVs 

One can also get all general TFPVs by computing a Gröbner basis `G` that
reflects necessary conditions on the parameters of the system. 
Note that depending on the input system and the drop in dimension this can be a
very computationally intensive task.
```@example 1
G = tfpvs_groebner(problem)
```
The TFPVs are points in the affine variety of `G`. 
Thus, in cases where `G` is simple enough, we can compute a minimal primary
decomposition to characterise the cases in which `G` vanishes. 
Note that here we use functions exported by `Oscar.jl`.
```@example 1 
I = ideal(G)
PD = primary_decomposition(I)
for Q in PD
    println(gens(Q[2]))
end
```
Thus, in this case, all TFPVs are slow-fast separations of rates.

In cases where `G` is more complicated, we can check whether the ideal generated
by `G` contains any monomials as follows. 
```@example 1
# polynomial ring ℚ[p,x]
R = parent(α)
# polynomial ring ℚ[p]
S, v = polynomial_ring(QQ, "_" .* p)
h = hom(S, R, system_parameters(problem))

# the ideal generated by G in the ring S
I = preimage(h, ideal(G));

# compute the saturation I:⟨π₁⋯πₘ⟩^∞
I_sat = saturation(I, ideal(prod(v)))
is_one(I_sat)
```









