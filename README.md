# TikhonovFenichelReductions.jl

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://jo-ap.github.io/TikhonovFenichelReductions.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://jo-ap.github.io/TikhonovFenichelReductions.jl/dev)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14411664.svg)](https://doi.org/10.5281/zenodo.14411664)


A Julia package for computing Tikhonov-Fenichel Parameter Values (TFPVs) for
polynomial ODE systems and their corresponding reductions (see [1-3] for
details).
The package is described and showcased in [4].

## Outline
Consider an ODE system of the form 
```math
\dot{x} = f(x,\pi, \varepsilon), \quad x(0)=x_0, x \in U\subseteq\mathbb{R}^n, \pi \in \Pi \subseteq \mathbb{R}^m,
```
where $f \in \mathbb{R}[x,\pi]$ is polynomial and $\varepsilon \geq 0$ is a small
parameter. 
The results from [1-3] allow us to compute a reduced system for $`\varepsilon
\to 0`$ in the sense of Tikhonov [5] and Fenichel [6] using methods from
commutative algebra and algebraic geometry. 

TikhonovFenichelReductions.jl implements methods for finding all possible
TFPVs.
It also includes functions to simplify the computation of the corresponding
reduced systems.
Note that this approach yields all possible timescale separations of _rates_
instead of _components_ as in Tikhonov's theorem [7] (it is thus a coordinate
free approach to geometric singular perturbation theory).

## Installation
Run
~~~
add TikhonovFenichelReductions
~~~
in Julia package Mode.
Note that this package relies on [Oscar.jl](https://www.oscar-system.org/) and
therefore needs to be installed in the Windows Subsystem for Linux on Windows.

## References
[1] A. Goeke and S. Walcher, ‘Quasi-Steady State: Searching for and Utilizing Small Parameters’, in Recent Trends in Dynamical Systems, A. Johann, H.-P. Kruse, F. Rupp, and S. Schmitz, Eds., in Springer Proceedings in Mathematics & Statistics, vol. 35. Basel: Springer Basel, 2013, pp. 153–178. doi: [10.1007/978-3-0348-0451-6_8](http://link.springer.com/10.1007/978-3-0348-0451-6_8).

[2] A. Goeke and S. Walcher, ‘A constructive approach to quasi-steady state reductions’, J Math Chem, vol. 52, no. 10, pp. 2596–2626, Nov. 2014, doi: [10.1007/s10910-014-0402-5](http://link.springer.com/10.1007/s10910-014-0402-5).

[3] A. Goeke, S. Walcher, and E. Zerz, ‘Determining “small parameters” for quasi-steady state’, Journal of Differential Equations, vol. 259, no. 3, pp. 1149–1180, Aug. 2015, doi: [10.1016/j.jde.2015.02.038](https://linkinghub.elsevier.com/retrieve/pii/S0022039615001102).

[4] J. Apelt and V. Liebscher, ‘Tikhonov-Fenichel Reductions and Their Application to a Novel Modelling Approach for Mutualism’, Theoretical Population Biology, pp. 16–35, Dec. 2025, doi: [10.1016/j.tpb.2025.08.004u](https://doi.org/10.1016/j.tpb.2025.08.004).

[5] A. N. Tikhonov, ‘Systems of differential equations containing small parameters in the derivatives’, Mat. Sb. (N.S.), vol. 73, no. 3, pp. 575--586, 1952, <https://www.mathnet.ru/php/archive.phtml?wshow=paper&jrnid=sm&paperid=5548&option_lang=eng>.

[6] N. Fenichel, ‘Geometric singular perturbation theory for ordinary differential equations’, Journal of Differential Equations, vol. 31, no. 1, pp. 53–98, Jan. 1979, doi: [10.1016/0022-0396(79)90152-9](https://linkinghub.elsevier.com/retrieve/pii/0022039679901529)

[7] F. Verhulst, ‘Singular perturbation methods for slow–fast dynamics’, Nonlinear Dynamics, vol. 50, pp. 747–753, 2007, doi: [10.1007/s11071-007-9236-z](https://doi.org/10.1007/s11071-007-9236-z)

## License
GNU GENERAL PUBLIC LICENSE, Version 3 or later (see LICENSE)
