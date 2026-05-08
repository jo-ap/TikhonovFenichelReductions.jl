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

Natural systems often evolve on different characteristic time scales, with
process rates differing substantially by orders of magnitude [@hek2010].
To study such systems, singular perturbation theory provides a mathematical
toolbox that translates this feature to ODE models and enables dimensionality
reduction by focusing on a single characteristic time scale, (typically the slow
one governing the long-term behaviour) [@wechselberger2020].

This separation of time scales is reflected by the presence of a small parameter
$\varepsilon>0$, which can be interpreted as the ratio between the two time
scales. Roughly speaking, the dynamics in slow time $\tau=\varepsilon t$ can be
approximated by replacing the fast part of the system with its steady state, so
that the fast dynamics effectively appear instantaneous on the slow time scale.

More precisely, we obtain the *reduced system* (or *reduction* in short) as the
limit

\begin{equation}
  \label{eq:tikhonov}
  \begin{aligned}
    \dot{u}             & = g(u,v) \\
    \varepsilon \dot{v} & = h(u,v)
  \end{aligned}
\quad\xrightarrow[\varepsilon\to0]{}\quad
 \begin{aligned}
    \dot{u} & = g(u,v), \quad u\in\mathbb{R}^s \\ 
    0       & = h(u,v), \quad v\in\mathbb{R}^r 
  \end{aligned}
\end{equation}

with $\dot{}=\mathrm{d}/\mathrm{d}\tau$ as stated in Tikhonov's theorem
[-@tikhonov1952; Theorem 1.1 in @verhulst2007]. The reduced system is defined
on the so-called slow manifold
$M_0=\{(u,v)\,\vert\,h(u,v)=0\}$.

For autonomous systems, the results of @fenichel1979 allow a clear geometric
interpretation --- hence the name geometric singular perturbation theory (GSPT)
--- and guarantee that for sufficiently small $\varepsilon>0$, the reduction
captures the behaviour of the full system.

In practice, GSPT can mitigate the realism-complexity trade-off, as modellers
can work with the reduction instead of the full system, which is often more
tractable.
Because the reduction remains embedded in the original system through the time
scale separation, the original interpretability is retained.
This allows one to follow the modelling paradigm:
Start with a complex but realistic model, reduce later.

The algebraic approach recently developed by Goeke \& Walcher (and colleagues)
allows one to work with GSPT systematically and is coordinate-free, i.e. we do
not rely on the standard form (\ref{eq:tikhonov}) with a time scale separation
of components.

Instead, we consider systems of the form 
\begin{equation}
    \label{eq:slowfastsystem}
    x' = f^{(0)}(x) + \varepsilon f^{(1)}(x) + \mathcal{O}(\varepsilon^2), 
    \quad x\in\mathbb{R}^n,
\end{equation}
where $f^{(0)}$ contains the fast processes and $f^{(1)}$ contains the slow
ones.
Note that system (\ref{eq:tikhonov}) can also be written in fast time $t$ as
\begin{equation}
    \label{eq:tikhonovfast}
    \begin{aligned}
        u' &= \varepsilon g(u,v) \\
        v' &= h(u,v)
    \end{aligned}
\end{equation}

where $'=\mathrm{d}/\mathrm{d}t$. 
Thus, we can translate (\ref{eq:slowfastsystem}) into
(\ref{eq:tikhonovfast}) by setting
$x = (u,v)^{\mathrm{T}}$, $f^{(0)} = (0, h)^{\mathrm{T}}$ and
$f^{(1)} = (g, 0)^{\mathrm{T}}$.
However, the reverse direction requires a suitable coordinate transformation,
which is not always feasible.
Fortunately, this transformation is not explicitly needed to obtain the reduced
system, which can be computed using (\ref{eq:reduction}).

The main benefit of this framework is however that it allows to find all
critical parameters admitting a reduction for polynomial or rational ODE systems
of the form (\ref{eq:slowfastsystem}).
Thus, modellers do not have to rely on prior knowledge to find suitable time
scale separations as they can systematically consider all of them.
This can be achieved by evaluating necessary conditions for the existence of a
reduction [@goeke2013; @goeke2013a; @goeke2014; @goeke2015].
`TikhonovFenichelReductions.jl` is a `Julia` [@bezanson2017] package
implementing this approach for polynomial ODE systems.
@apelt2025 provide a showcasing example and more detailed explanations.

# Statement of need

The ad-hoc approach to singular perturbation theory requires prior knowledge of
a suitable time scale separation and substantial mathematical effort to compute
the reduction.
In contrast, the algebraic framework developed by Goeke \& Walcher
systematically yields *all* possible reductions for a given polynomial ODE
system, enabling modellers to consider several scenarios defined by different
time scale separations [@goeke2015].
`TikhonovFenichelReductions.jl` makes the required computations easily
accessible, even for non-expert users.

To the author's knowledge, no publicly available implementation of this theory
currently exists.
However, there is an implementation of a related approach for computing
invariant manifolds by Roberts [-@roberts1997; -@roberts-webapp].

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

Detailed explanations and examples are provided by @apelt2025 and in the
[documentation](https://jo-ap.github.io/TikhonovFenichelReductions.jl/stable/).
The core features are summarized in the following.

## Finding TFPVs
The package provides a method to obtain all possible TFPVs implicitly by
computing a Gröbner Basis and one to find *slow-fast separations of rates*,
i.e. TFPVs with some parameters set to zero.
Although the former method is an exhaustive search, the latter is usually better
suited in practice as it yields the TFPVs one is typically interested in
explicitly, is more efficient, and directly outputs the corresponding slow
manifolds (implicitly as affine varieties in phase space). 
In both cases, one utilizes the correspondence between algebraic and geometric
properties of the slow manifold, which renders this an algebraic approach to
GSPT.

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
\begin{equation}
    \label{eq:reduction}
    x' = \left[1_n - P(x)(D\psi(x)P(x))^{-1} D\psi(x)\right] f^{(1)}(x).
\end{equation}

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
However, as far as the author is aware, the systematic approach due to Goeke \&
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
In this case, `TikhonovFenichelReductions.jl` may enable and simplify the search
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
