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
date: 14 January 2026
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
existence of a reduction, which allows to find *all* reductions of a given
polynomial ODE system using methods from computational algebra, and simplifies
the computation of reduced systems [@goeke2015].
`TikhonovFenichelReductions.jl` makes the required computations easily
accessible (even for non-expert users) by utilizing `Oscar.jl` [@oscar2022;
@decker2025].

The author is not aware of any publicly available implementation of the theory
by Goeke and Walcher, but there exist an implementation for a related approach
for computing invariant manifolds by Roberts [-@roberts1997; -@roberts-webapp].

# Software design 
`TikhonovFenichelReductions.jl` is implemented in `Julia` due to its
flexibility, the use of multiple dispatch and the availability of the
feature-rich computer algebra system `Oscar.jl` [@oscar2022; @decker2025]. 

Core features are the search for critical parameters admitting a reduction,
so-called *Tikhonov-Fenichel Parameter Values (TFPVs)*, and the computation of
the corresponding reduced systems.
Crucially, this requires various computations with multiple symbolic
representations (i.e. different polynomial rings, rational function fields and
matrix spaces) and therefore parsing of data between the corresponding types in
`Oscar.jl`, which `TikhonovFenichelReductions.jl` performs hidden away from the
user. 
Thus, the user mostly works in an object-oriented manner with the types and
methods provided.

The package is essentially designed around two types:
`ReductionProblem`, which constructs all symbolic data types needed for the
search of TFPVs, 
and `Reduction`, which holds all relevant information for the reduced
system and the steps required to compute it.
The latter also contains the reduced system and other information parsed to the
appropriate types from `Oscar.jl`, that can be further used, e.g. for a symbolic
analysis.

# Features
Detailed explanations and examples are provided by @apelt2025 and in the
[documentation](https://jo-ap.github.io/TikhonovFenichelReductions.jl/stable/).

## Finding TFPVs
The package provides a method to obtain all possible TFPVs implicitly by
computing a Gröbner Basis and one to find *slow-fast separations of rates*,
i.e. TFPVs with some parameters set to zero.
Although the former method is an extensive search, the latter is usually better
suited in practice as it yields the TFPVs one is typically interested in
explicitly, is more efficient, and directly outputs the corresponding slow
manifolds (implicitly as affine varieties in phase space). 

## Computing reductions
Computing a reduction for a slow-fast separation of rates as in
Theorem 1 in @goeke2014 requires essentially two steps.
First, we need to provide a parametric representation of the slow manifold,
which is given as an irreducible component of the affine variety
$\mathcal{V}(f^{(0)})$.
Then, we need to find a product decomposition $f^{(0)}=P \psi$, where $P$ is
a $n\times r$ matrix of rational functions and $\psi$ is a vector of polynomials
locally satisfying $\mathcal{V}(\psi)=\mathcal{V}(f^{(0)})$ and $\text{rank}\;P
= \text{rank}\;D\psi = r$.

With this, the reduced system in the sense of Tikhonov is given by
$$
    \dot{x} = \left[1_n - P(x)(D\psi(x)P(x))^{-1} D\psi(x)\right] f^{(1)}(x).
$$ 

For convenience, the packages provides multiple heuristics, which automate the
steps required to compute a reduction and allows bulk computation of multiple
reductions at once.

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

# Research impact statement
Time scale separations are widely used in various areas of mathematical
modelling [@wechselberger2020]. 
However, as far as the author is aware, the systematic approach due to Goeke and
Walcher seems to be scarcely adopted, even though it comes with many advantages
compared to the ad-hoc approach. 
This package aims to make the theory more accessible and convenient to use. 

In the field of mathematical ecology (which the author is most familiar with),
the theory was successfully used by @kruff2019 and more recently by @apelt2025.
For the model introduced in the former, the application was relatively
straightforward, but more complex models as in the latter require some more
work.
In particular, the method for finding TFPVs that relies on the computation of a
Gröbner Basis may fail for complex systems. 
In this case, `TikhonovFenichelReductions.jl` enables and simplifies the search
for the most common TFPVs.
Given the ubiquity of time scale separation techniques in this field alone 
[see e.g. @abbott2020; @poggiale2004; @revilla2015], the package can be a
potentially useful tool for modellers.

# AI usage disclosure
No generative AI tools were used in the development of this software, the writing
of this manuscript, or the preparation of supporting materials.

# Acknowledgements
This work was supported by a scholarship awarded by the University of Greifswald 
according to the "Landesgraduiertenförderungsgesetz (LGFG) MV".
I like to thank Volkmar Liebscher for his supervision and support,
Sebastian Walcher and Alexandra Goeke for developing the Tikhonov-Fenichel
reduction theory,
Leonard Schmitz for discussing aspects of the computational algebra, 
and the OSCAR team for their great work and helpful replies to my questions.

# References
