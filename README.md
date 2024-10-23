# TikhonovFenichelReductions.jl

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://jo-ap.github.io/TikhonovFenichelReductions.jl/dev)

A Julia package for computing Tikhonov-Fenichel Parameter Values (TFPVs) for a
polynomial ODE system and the corresponding reductions (see [1-3] for details).

*This package is still under development*

## Background
Consider an ODE system of the form 
```math
\dot{x} = f(x,\pi, \varepsilon), \quad x(0)=x_0, x \in U\subseteq\mathbb{R}^n, \pi \in \Pi \subseteq \mathbb{R}^m,
```
where $f \in \mathbb{R}[x,\pi]$ is polynomial and $\varepsilon \geq 0$ is a small
parameter. The results from [1-3] allow us to compute a reduced system for
$`\varepsilon \to 0`$ in the sense of Tikhonov [4] and Fenichel [5] using
methods from commutative algebra and algebraic geometry. 

TikhonovFenichelReductions.jl implements methods for finding all possible
TFPVs.
It also includes functions to simplify the computation of corresponding
reduced systems.
Note that this approach yields all possible timescale separations of _rates_
and not just _components_.

## Installation
Run
~~~
add https://gitlab.com/jo-ap/tikhonovfenichelreductions.jl
~~~
in Julia package Mode.

## Example
Have a look at the file `example.jl`, which uses this package to compute all
TFPVs as discussed in [6] and computes the Rosenzweig-MacArthur model as a
reduction from a three dimensional system.
 
## References
[1] A. Goeke and S. Walcher, ‘Quasi-Steady State: Searching for and Utilizing Small Parameters’, in Recent Trends in Dynamical Systems, A. Johann, H.-P. Kruse, F. Rupp, and S. Schmitz, Eds., in Springer Proceedings in Mathematics & Statistics, vol. 35. Basel: Springer Basel, 2013, pp. 153–178. doi: [10.1007/978-3-0348-0451-6_8](http://link.springer.com/10.1007/978-3-0348-0451-6_8).

[2] A. Goeke and S. Walcher, ‘A constructive approach to quasi-steady state reductions’, J Math Chem, vol. 52, no. 10, pp. 2596–2626, Nov. 2014, doi: [10.1007/s10910-014-0402-5](http://link.springer.com/10.1007/s10910-014-0402-5).

[3] A. Goeke, S. Walcher, and E. Zerz, ‘Determining “small parameters” for quasi-steady state’, Journal of Differential Equations, vol. 259, no. 3, pp. 1149–1180, Aug. 2015, doi: [10.1016/j.jde.2015.02.038](https://linkinghub.elsevier.com/retrieve/pii/S0022039615001102).

[4] A. N. Tikhonov, ‘Systems of differential equations containing small parameters in the derivatives’, Mat. Sb. (N.S.), vol. 73, no. 3, pp. 575--586, 1952, <https://www.mathnet.ru/php/archive.phtml?wshow=paper&jrnid=sm&paperid=5548&option_lang=eng>.

[5] N. Fenichel, ‘Geometric singular perturbation theory for ordinary differential equations’, Journal of Differential Equations, vol. 31, no. 1, pp. 53–98, Jan. 1979, doi: [10.1016/0022-0396(79)90152-9](https://linkinghub.elsevier.com/retrieve/pii/0022039679901529)

[6] N. Kruff, C. Lax, V. Liebscher, and S. Walcher, ‘The Rosenzweig–MacArthur system via reduction of an individual based model’, J. Math. Biol., vol. 78, no. 1–2, pp. 413–439, Jan. 2019, doi: [10.1007/s00285-018-1278-y](http://link.springer.com/10.1007/s00285-018-1278-y).

## Dependencies
This packages only works due to the great [Oscar.jl](https://www.oscar-system.org/) project.

## License
GNU GENERAL PUBLIC LICENSE, Version 3 or later (see LICENSE)
