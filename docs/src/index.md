# TikhonovFenichelReductions.jl

A Julia package for computing Tikhonov-Fenichel Parameter Values (TFPVs) for a
polynomial ODE system and the corresponding reductions (see 
[goeke2013a,goeke2014,goeke2015](@cite) for details).

## Overview
The general framework for this package is singular perturbation theory.
More precisely, we consider an ODE system of the form 
```math
\dot{x} = f(x,\pi, \varepsilon), \quad x(0)=x_0, x \in U\subseteq\mathbb{R}^n, \pi \in \Pi \subseteq \mathbb{R}^m,
```
where ``f \in \mathbb{R}[x,\pi]`` is polynomial and ``\varepsilon \geq 0`` is a
small parameter. The results from [goeke2013a,goeke2014,goeke2015](@cite) allow
us to compute a reduced system for ``\varepsilon \to 0`` in the sense of Tikhonov
[tikhonov1952](@cite) and Fenichel [fenichel1979](@cite) using methods from
commutative algebra and algebraic geometry. 

TikhonovFenichelReductions.jl implements methods for finding all possible TFPV
candidates, i.e. separations of parameters into small and large (which in turn
correspond to slow-fast separations of processes in the original system).
It also includes functions to simplify the computation of corresponding
reduced systems.
Note that this approach yields all possible time-scale separations of _rates_ and
not just _components_ as in the classical approach.

## Basic Usage 

### Finding TFPV candidates
The main functionality is:
~~~
# Setup
problem = ReductionProblem(f, x, θ, s, idx_slow_fast)

# Compute TFPV candidates
idx, G, (V, dim_Y) = tfpv_candidates(problem)
~~~
This sets up the problem. The input is a (Julia) function `f(x,θ)` defining
``f`` and a vector of strings with the names of the dynamic variables ``x`` and
parameters ``\theta``, respectively. 
``f`` needs to be a multivariate polynomial in ``x`` and ``\theta`` defining the
RHS of the original system. The dimension of the reduced system is ``s < n``.  
`π=θ[idx_slow_fast]` are the parameters that should be considered
to be either small or large, i.e. all other parameters are fixed. 

The output is a vector of boolean indices `idx`, where `0` corresponds to a
small and `1` to a large parameter. 
`G` is a Gröbner basis of an ideal in ``\mathbb{R}[\pi]`` (the elimination
ideal of ``f`` and all determinants of ``k\times k`` minors of the Jacobian for
any ``k>s``). 
All TFPVs are contained in the affine variety of `G`, i.e. all polynomials in
`G` need to vanish for a TFPV ``\pi^\star`` in order for a reduction to exist.
`V`contains the generators of primary ideals that correspond to the irreducible
components of the affine variety ``\mathcal{V}(f(\cdot,\pi))``. 
By affine variety we simply think of the zero set of some polynomials (or ideal,
for which it is enough to consider its generators) the in the corresponding
affine space.
In this case, we take the polynomials in ``f`` as elements of the ring
``\mathbb{R}[x]`` and ``\pi`` as coefficients (in fact, technically we work in
the polynomial ring ``\mathbb{R}(\pi)[x]``).

If finding TFPV candidates takes too long, you can disable the computation of
the irreducible components of `V` by setting
`compute_primary_decomposition=false`.
`dim_Y` contains the Krull dimension of the irreducible components. 
These are equal to the topological dimensions in ``\mathbb{C}[x]``, but can be
larger in ``\mathbb{R}[x]``. 
However, they are also equal over the reals if a component contains a
non-singular point --- which we will require later anyway. 
We can therefore use this information already to discard unwanted cases. 
The default behaviour is to compute a primary decomposition and to check if any
component's dimension is exactly `s`. 
If you want to disable this filter, you can set `exact_dimension=false` (this
only has an effect if `compute_primary_decomposition=true`).

The functions 
~~~
print_tfpv(idx, prob)
print_varieties(V, prob)
~~~
print out the results in a simple format directly to the REPL. 
You can use the additional argument `latex=true` to print output as LaTeX source
code. 

#### Limitations

The method `tfpv_candidates` can only find TFPV candidates for which a subset of
the original parameters is set to zero. 
For more complicated candidates, i.e. when some function of the parameters is small, 
one has to consider the Gröbner Basis `G` directly.
These 'alternative TFPVs' are rational functions whose vanishing results in the
vanishing of `G`.

### Computing the reduced system

To compute the reduced system, we first need to make the variables available in
the `Main` namespace (parsed to the appropriate types from 
[Oscar.jl](https://www.oscar-system.org/)).
These can be accessed as `problem.x` and `problem.θ`. 
Then, we initialise a reduction by constructing an instance of type `Reduction`. 
For the first slow-fast separation in `idx`, this can be done with 
~~~
reduction = Reduction(prob, idx[1])
~~~
Looking at `V[1]`, we can see what a possible slow manifold should be, 
since it is contained in one of the irreducible components of the vanishing set
of the ideal generated by the polynomials in `V[1]`
Each `V[i]` is a vector containing the generators of these ideals (as a vector
of polynomials).

Before we can compute the reduced system, we can check if the necessary
conditions for its existence are satisfied. 
These can be found in [1-3]. Essentially, we need to find a manifold of
dimension ``s`` contained in ``\mathcal{V}(f(\cdot, \theta^\star))``, a
non-singular point ``x_0`` on this manifold and a product decomposition 
``f(\cdot, \theta^\star) = P\cdot \psi`` 
that locally satisfies 
``\mathcal{V}(f(x, \theta^\star)) = \mathcal{V}(\psi)``.

~~~
# define the slow manifold M₀ (M₀ has to have the same format as problem.x)
set_manifold!(reduction, M₀)

# choose a non-singular point x₀ on M₀
set_point!(reduction, x₀)

# define a product decomposion f⁰ = P⋅ψ, where ψ is a r×1 matrix of polynomials
# locally satisfying V(ψ) = V(f⁰). 
set_decomposition!(reduction, P, ψ)
# or compute P automatically
set_decomposition!(reduction, ψ)
~~~

When setting the slow manifold, a check is performed that indicates whether a
generic point on the manifold can be non-singular.
If that is not the case, a warning is printed.
Otherwise the generic point is set for the reduction (the
method `set_point` is only useful when a specific point should be chosen, i.e.
one that simplifies the Jacobian at `x_0`). 

If you choose `ψ` as a vector of polynomials in ``\mathbb{R}[x,\pi]``, it should
be possible to automatically compute `P` (which is then a matrix of rational
functions). 
However, as of now this might fail in complicated cases, but is guaranteed to
work if `r=n-s=1`.

If all the conditions for the slow manifold, the non-singular point and the
product decomposition are satisfied (indicated by the return value `true` of
each function), we can compute the reduced system with
~~~
g_raw, g = compute_reduction(reduction)
~~~
where `g` corresponds to the RHS of the reduced system on `M₀`.
`g_raw` is the same function, but before substituting the reduced dynamic variables.

## Dependencies
This packages only works due to the great [Oscar.jl](https://www.oscar-system.org/) project.

## License
GNU GENERAL PUBLIC LICENSE, Version 3 or later (see LICENSE)