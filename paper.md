---
title: 'TikhonovFenichelReductions.jl: A systematic approach to geometric singular perturbation theory'
tags:
  - Julia
  - Tikhonov
  - Fenichel
  - singular perturbation theory
  - time scale separations 
  - dynamical systems

authors:
  - name: Johannes Apelt
    affiliation: 1
affiliations:
 - name: Institute of Mathematics and Computer Science, University of Greifswald, Germany
   index: 1
date: 29 November 2025
bibliography: paper.bib
---

# Summary
Singular perturbation theory is a mathematical toolbox that allows
dimensionality reduction of ODE systems whose components evolve on different
time scales emanating from the presence of a small parameter $\varepsilon$. 
More precisely, we obtain the *reduced system* as the limit
\begin{equation}
  \label{eq:tikhonov}
  \begin{aligned}
    \dot{u}             & = g(u,v) \\
    \varepsilon \dot{v} & = h(u,v) 
  \end{aligned}
\quad\xrightarrow[\varepsilon\to0]{}\quad
 \begin{aligned}
    \dot{u} & = g(u,v) \\
    0       & = h(u,v)
  \end{aligned}
\end{equation}
as stated in Tikhonov's theorem [-@tikhonov1952; Theorem 1.1 in @verhulst2007].
For autonomous systems, @fenichel1979 established a geometric singular
perturbation theory (GSPT) in coordinate-free settings [see @wechselberger2020].
Convergence properties then follow from the slow manifold 
$M_0=\{(u,v)\,\vert\,h(u,v)=0\}$, on which the reduced system is defined.

The algebraic approach to GSPT recently developed by Goeke and Walcher (and
colleagues) allows to systematically find all critical parameters admitting a
reduction for polynomial or rational ODE systems of the form
\begin{equation}
    \label{eq:slowfastsystem}
    \dot{x} = f(x,\pi) = f^{(0)}(x) + \varepsilon f^{(1)}(x) + \mathcal{O}(\varepsilon^2), 
    \quad x\in\mathbb{R}^n, \pi\in\mathbb{R}^m,
\end{equation}
i.e. we obtain all reductions for a system with a slow-fast separation of
processes instead of components as in the standard form (\ref{eq:tikhonov}),
which renders this a coordinate-free approach.
This can be achieved by evaluating necessary conditions for the existence of a
reduction for a system as in (\ref{eq:slowfastsystem}) [@goeke2013; @goeke2013a;
@goeke2014; @goeke2015].

`TikhonovFenichelReductions.jl` is a `Julia` [@bezanson2017] package
implementing this approach for polynomial ODE systems.
@apelt2025 provide a showcasing example and more detailed explanations.

# Statement of need
The ad-hoc approach to singular perturbation theory requires prior knowledge
about a suitable time scale separation and substantial mathematical effort to
compute the reduction.
The algebraic approach yields algorithmically accessible conditions for the
existence of a reduction, which allows us to find all reductions of a given
polynomial ODE system using methods from computational algebra, and simplifies
the computation of reduced systems [@goeke2015].
`TikhonovFenichelReductions.jl` makes the required computations easily
accessible (even for non-expert users) by utilizing `Oscar.jl` [@oscar2022;
@decker2025].

The author is not aware of any publicly available implementation of the theory
by Goeke and Walcher, but there exist an implementation for a related approach
for computing invariant manifolds by Roberts [-@roberts1997; -@roberts-webapp].

# Features
The main features provided by `TikhonovFenichelReductions.jl` are the search for
critical parameters admitting a reduction, so-called *Tikhonov-Fenichel
Parameter Values (TFPVs)*, and the computation of the corresponding reduced
systems.
Detailed explanations and examples are provided by @apelt2025 and in the
[documentation](https://jo-ap.github.io/TikhonovFenichelReductions.jl/stable/).

## Finding TFPVs
The package provides two methods for finding TFPVs admitting a reduction onto an
$s$-dimensional slow manifold.

Restricting the search for TFPVs to *slow-fast separations of rates*, i.e. TFPVs
with some parameters set to zero, is typically more efficient and yields the
TFPV candidates together with their slow manifolds given implicitly as the
irreducible components of the corresponding affine variety
$\mathcal{V}(f^{(0)})=\{x\in\mathbb{R}^n \mid f^{(0)}(x)=0\}$.

Computing a Gröbner basis $G$ with an elimination ordering for the state
variables of an ideal obtained from necessary conditions for the existence of a
reduction yields all TFPVs implicitly as $\mathcal{V}(G)$.
This may be infeasible to compute, but can reveal more complex expressions that
can be considered as small parameters.
However, introducing new parameters from these expressions typically allows us
to use the methods for slow-fast separations of rates, particularly the
computation of the reduced system.

## Computing reductions
Computing a reduction for a slow-fast separation of rates as in
Theorem 1 in @goeke2014 requires essentially two steps.
Let $Y$ be the irreducible component of $\mathcal{V}(f^{(0)})$ with dimension
$s$ corresponding to the slow manifold $M_0$ and $r=n-s$. 

First, we need to find a parametric representation of $Y$ as a manifold.
This may be done automatically by a provided heuristic. 
If this fails, $Y$ must be parameterized explicitly w.r.t. $s$ state variables
(the local coordinates) manually.

Next, we need to find a product decomposition $P\psi = f^{(0)}$, where $P$ is a
$n\times r$ matrix of rational functions and $\psi$ is a vector of polynomials
locally satisfying $\mathcal{V}(\psi)=\mathcal{V}(f^{(0)})$ and $\text{rank}\;P
= \text{rank}\;D\psi = r$.
If $Y$ is defined by $r$ generators or $r$ independent entries of $f^{(0)}$ can
be found by a heuristic, $P$ and $\psi$ can be computed automatically.

With this, the reduced system in the sense of Tikhonov is given by
$$
    \dot{x} = \left[1_n - P(x)(D\psi(x)P(x))^{-1} D\psi(x)\right] f^{(1)}(x).
$$ 

For complex systems, there may exist many TFPVs each admitting multiple
reductions, which makes their analysis intractable. 
Thus, `TikhonovFenichelReductions.jl` contains methods for finding all unique
slow manifolds (as varieties) and computing all reductions onto each of them.
This relies on the heuristics for finding a parametric representations of the
varieties and setting a product decomposition for $f^{(0)}$, but in most cases
this means we can compute all or at least most reductions fully automatically. 

## Integration with the `Julia` ecosystem
The input system can be given as a reaction network defined with `Catalyst.jl`
[@loman2023].
Because the reduced systems are represented using types from `Oscar.jl`, the
latter's functions can be used to aid the symbolic analysis. 
`Julia`'s support for metaprogramming  allows to perform further tasks such as a
numerical analysis without having to copy or parse code (see e.g.
[`TFRSimulations.jl`](https://github.com/jo-ap/TFRSimulations)).
For convenience, there are several methods for displaying the output, including
printing as \LaTeX{} source code via
[`Latexify.jl`](https://github.com/korsbo/Latexify.jl).

# Acknowledgements
This work was supported by a scholarship awarded by the University of Greifswald 
according to the "Landesgraduiertenförderungsgesetz (LGFG) MV".
I like to thank Volkmar Liebscher for his supervision and support,
Sebastian Walcher and Alexandra Goeke for developing the Tikhonov-Fenichel
reduction theory,
Leonard Schmitz for discussing aspects of the computational algebra, 
and the OSCAR team for their helpful replies to my questions.

# References
