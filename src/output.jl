## Output functions

assert_numeric_idx(idx::Vector{Bool}) = (1:length(idx))[idx]
assert_numeric_idx(idx::AbstractVector{Int}) = idx

"""
    $(TYPEDSIGNATURES)

Pad string `s` to predefined width or keep unchanged if `width ≤ length(str)`.
Setting `padfront=false` inserts padding after the string (left aligned)
instead of in front. 
"""
function padstring(s::String, width::Int; padfront::Bool=true)
  str_length = length(s)
  if width <= str_length
    str_out =  s
  else
    pad_str = repeat(" ", width - str_length)
    str_out = padfront ? pad_str * s : s * pad_str
  end
  return str_out
end

"""
    print_tfpvs(
      [io::IO,] 
      problem::ReductionProblem,
      sf_separations::Vector{Vector{Bool}}; 
      latex::Bool=false,
      idx::Union{AbstractVector{Int}, Vector{Bool}}=1:length(sf_separations)
    )

Print slow-fast separations to terminal or use `latex=true` to print LaTeX
instead. Optionally only print subset defined by (boolean or numeric) index set
`idx`.
"""
function print_tfpvs(
    io::IO,
    problem::ReductionProblem,
    sf_separations::Vector{Vector{Bool}}; 
    latex::Bool=false,
    idx::Union{AbstractVector{Int}, Vector{Bool}}=1:length(sf_separations)
  )
  _print_tfpvs(io, problem, sf_separations, latex, assert_numeric_idx(idx))
end
function print_tfpvs(
    problem::ReductionProblem,
    sf_separations::Vector{Vector{Bool}}; 
    latex::Bool=false,
    idx::Union{AbstractVector{Int}, Vector{Bool}}=1:length(sf_separations)
  )
  _print_tfpvs(stdout, problem, sf_separations, latex, assert_numeric_idx(idx))
end

function _print_tfpvs(
    io::IO, 
    problem::ReductionProblem,
    sf_separations::Vector{Vector{Bool}},
    latex::Bool,
    idx::AbstractVector{Int}
  )
  if latex
    parameters = [latexify(p; env=:raw) for p in string.(problem.p)]
    m = length(parameters)
    str = "\$\\begin{array}{r$(repeat('c',length(problem.p)))} i & $(join(parameters .* "^\\star", " & ")) \\\\ \n"
    for i in idx
      str *= "$i & " * join([sf_separations[i][k] ? parameters[k] : "\\cdot" for k = 1:m], " & ") * " \\\\ \n"
    end
    str *= "\\end{array}\$"
  else
    subscripts = ["₀","₁","₂","₃","₄","₅","₆","₇","₈","₉"]
    numbers = string.(0:9)
    max_width = ndigits(length(sf_separations)) 
    field_widths = [length(pᵢ) for pᵢ in string.(problem.p)]
    str = "π $(repeat(" ", max_width-1)) = ($(join(string.(problem.p), ", ")))\n"
    str *= repeat("_", length(str)-2) * "\n"
    for i in idx
      idx_fast = sf_separations[i]
      str *= "π" *
      join([subscripts[string(n) .== numbers][1] for n in string(i)], "") *
      "$(repeat(" ", max_width - ndigits(i))) = (" * 
      join(padstring.([idx_fast[k] ? string(problem.p[k]) : "⋅" for k = 1:length(idx_fast)], field_widths), ", ") *")\n"
    end
  end
  println(io, str)
end

"""
    print_varieties(
      [io::IO,]
      V::Vector{Vector{Variety}};
      latex::Bool=false,
      idx::Union{AbstractVector{Int}, Vector{Bool}}=1:length(V)
    ) 

Print generators of ideals corresponding to the irreducible components of
varieties `V(f0)` for TFPV candidates and their dimension as returned by
`tfpvs_and_varieties` (these are wrapped in the type `Variety`).
Use keyword argument `latex=true` to print LaTeX code instead.

See also: [`tfpvs_and_varieties`](@ref), [`print_tfpvs`](@ref), [`print_results`](@ref)
"""
function print_varieties(
    io::IO, 
    V::Vector{Vector{Variety}};
    latex::Bool=false,
    idx::Union{AbstractVector{Int}, Vector{Bool}}=1:length(V)
  )
  _print_varieties(io::IO, V, latex, assert_numeric_idx(idx))
end
function print_varieties(
    V::Vector{Vector{Variety}};
    latex::Bool=false,
    idx::Union{AbstractVector{Int}, Vector{Bool}}=1:length(V)
  )
  _print_varieties(stdout, V, latex, assert_numeric_idx(idx))
end

function _print_varieties(
    io::IO, 
    V::Vector{Vector{Variety}},
    latex::Bool,
    idx::AbstractVector{Int}
  )
  if latex
    str = "\\[\n\\begin{array}{rll} i & \\mathcal{V}(f(\\cdot, \\pi^\\star)) & \\dim_{\\textrm{Krull}} \\\\ \n"
    for i in idx
      str *= "$i"
      for j in eachindex(V[i])
        str *= " & \\langle " * join(latexify.(string.(V[i][j].gens_R); env=:raw, mult_symbol=""), ", ") * " \\rangle & $(V[i][j].dim) \\\\ \n"
      end
    end
    str *= "\\end{array}\n\\]"
  else
    subscripts = ["₀","₁","₂","₃","₄","₅","₆","₇","₈","₉"]
    numbers = string.(0:9)
    max_width = ndigits(length(V)) 
    pad = "$(repeat(" ", max_width + 4))"
    str = ""
    for i in idx
      V_str = ["[" * replace(join(string.(V[i][j].gens_R), ", ")) * "], $(V[i][j].dim)" for j in eachindex(V[i])]
      str *= "V" * 
        join([subscripts[string(n) .== numbers][1] for n in string(i)], "") * 
        "$(repeat(" ", max_width - ndigits(i))) : " * 
        join(V_str, "\n" * pad) * "\n"
    end
  end
  println(io::IO, str)
end

"""
    print_results(
      [io::IO,]
      problem::ReductionProblem,
      sf_separations::Vector{Vector{Bool}},
      V::Vector{Vector{Variety}};
      idx::Union{AbstractVector{Int}, Vector{Bool}}=1:length(V)
    )

Print slow-fast separations that are TFPVs and generators for the corresponding
irreducible components of `V(f0)` together with their Krull dimension as stored
in `Variety` object.

### Arguments 
- `problem`: `ReductionProblem` 
- `sf_separations`: Boolean indices defining all TFPVs `π⁺` that are slow-fast separations (0: slow, 1: fast).
- `V`: Generators for the irreducible component of the affine varietiy `V(f(⋅,π⁺))` for each slow-fast separation and their respective dimension.
- `idx`: (optional) index vector to include only certain TFPVs (boolean or numeric)

See also: [`tfpvs_and_varieties`](@ref), [`print_tfpvs`](@ref), [`print_varieties`](@ref)
"""
function print_results(
    io::IO,
    problem::ReductionProblem,
    sf_separations::Vector{Vector{Bool}},
    V::Vector{Vector{Variety}};
    idx::Union{AbstractVector{Int}, Vector{Bool}}=1:length(V)
  )
  _print_results(io, problem, sf_separations, V, assert_numeric_idx(idx))
end
function print_results(
    problem::ReductionProblem,
    sf_separations::Vector{Vector{Bool}},
    V::Vector{Vector{Variety}};
    idx::Union{AbstractVector{Int}, Vector{Bool}}=1:length(V)
  )
  _print_results(stdout, problem, sf_separations, V, assert_numeric_idx(idx))
end

function _print_results(
    io::IO,
    problem::ReductionProblem,
    sf_separations::Vector{Vector{Bool}},
    V::Vector{Vector{Variety}},
    idx::Union{AbstractVector{Int}, Vector{Bool}}
  )
  for i in idx
    println(io, "$i")
    idx_fast = sf_separations[i]
    println(io, " π̃: [" * join([idx_fast[k] ? string(problem.p[k]) : "0" for k = 1:length(idx_fast)], ", ") *"]")
    V_str = ["[" * replace(join(string.(V[i][k].gens_R), ", ")) * "], $(V[i][k].dim)" for k in eachindex(V[i])]
    println(io, " V: " * join(V_str, "\n    ") * "\n")
  end
end

"""
    print_reduced_system(
      [io::IO,] 
      reduction::Reduction; 
      rewrite::Bool=true,
      factor::Bool=false, 
      latex::Bool=false,
      local_coordinates::Bool=true
    )

Print the reduced system (after `compute_reduction!` has been called
successfully on the `Reduction` object). 
If the slow manifold was successfully specified, this returns the system on the
slow manifold. 
If  `local_coordinates=true`, the reduced system is only given in
`s`-dimensions (i.e. the local coordinates on the slow manifold). 
If the reduction could be computed but the slow manifold was not set
successfully, the reduced system is shown in the original phase space.

### Arguments 
- `reduction`: `Reduction` holding the reduced system 
- `rewrite`: Whether the RHS should be decomposed into polynomial and rational part (default is `true`), see `rewrite_rational`
- `factor`: Whether the polynomial parts should be factorized
- `latex`: Whether to print latex string 
- `local_coordinates`: Whether to present reduced system in local coordinates on slow manifold (i.e. as `s` dimensional ODE system)

See also: [`rewrite_rational`](@ref), [`print_tfpvs`](@ref), [`print_varieties`](@ref)
"""
function print_reduced_system(io::IO, reduction::Reduction; rewrite::Bool=true, factor::Bool=false, latex::Bool=false, local_coordinates::Bool=true)
  str = _get_reduced_system_str(reduction; rewrite=rewrite, factor=factor, padfront=0, local_coordinates=local_coordinates)
  if latex && str != ""
    str_latex = [replace(latexify(s; env=:raw, mult_symbol=""), "=" => "&=") for s in split(str, "\n")]
    str = "\\begin{aligned}" * join(str_latex, "\\\\ \n") * "\\end{aligned}"
  end
  println(io, str)
end
function print_reduced_system(reduction::Reduction; rewrite::Bool=true, factor::Bool=false, latex::Bool=false, local_coordinates::Bool=true)
  print_reduced_system(stdout, reduction; rewrite, factor, latex, local_coordinates)
end

function get_system_str(f, u; rewrite::Bool=true, factor::Bool=false, padfront::Int=0)
  dxdt = ["d$_u/dt" for _u in string.(u)]
  string_length = maximum(length.(dxdt))
  pad = repeat(" ", padfront)
  str_out = ["" for _ in f]
  if rewrite
    _f = rewrite_rational.(f; factor=factor)
    for i in eachindex(_f)
      h, r, q = _f[i]
      str = pad * padstring(dxdt[i], string_length; padfront=false) * " = "
      if r == 0 
        _h = h == 0 ? "0" : string(h) 
        str_out[i] = str *  replace(_h, " * " => "*", "1 * " => "") 
      else
        _h = h == 0 ? "" : string(h) * " + "
        str_out[i] = str * replace(_h * "(" * string(r) * ")//(" * string(q) * ")", " * " => "*", "1 * " => "") 
      end
    end
  else 
    for i in eachindex(f) 
      str_out[i] *= pad * padstring(dxdt[i], string_length; padfront=false) * " = " * string(f[i])
    end
  end
  return join(str_out, "\n")
end

function _get_reduced_system_str(reduction::Reduction; rewrite::Bool=true, factor::Bool=false, padfront::Int=0, local_coordinates::Bool=true)
  if !any(reduction.reduction_cached)
    @warn "Reduced system is not yet computed"
    return ""
  end
  g = reduction.reduction_cached[2] ? reduction.g : reduction.g_raw
  x = reduction.problem.x
  if local_coordinates
    g = g[reduction.idx_components]
    x = x[reduction.idx_components]
  end
  return get_system_str(g, x; rewrite=rewrite, factor=factor, padfront=padfront)
end

"""
    $(TYPEDSIGNATURES)

Decompose a rational function `f = p/q` into polynomial and rational part, i.e.
return `f = h + r/q`, where `p,q,h,r` are polynomials.

### Arguments 
- `factor`: Whether all polynomials should be factorized (default is `false`)

See also: [`print_reduced_system`](@ref)
"""
function rewrite_rational(term::FracFieldElem{QQMPolyRingElem}; factor::Bool=false)
  p = numerator(term)
  q = denominator(term)
  h,r = divrem(p,q)
  if factor 
    h = h == 0 ? h : Oscar.factor(h)
    r = r == 0 ? r : Oscar.factor(r)
    q = Oscar.factor(q)
  end
  return h, r, q
end

"""
    print_slow_fast([io::IO,] reduction::Reduction; latex::Bool=false)

Print slow and fast parameters.
"""
function print_slow_fast(io::IO, reduction::Reduction; latex::Bool=false)
  print(io, _get_slow_fast_str(reduction; padfront=0, latex=latex))
end
function print_slow_fast(reduction::Reduction; latex::Bool=false)
  print(_get_slow_fast_str(reduction; padfront=0, latex=latex))
end

function _get_slow_fast_str(reduction::Reduction; padfront::Int=0, latex::Bool=false)
  slow, fast = _get_slow_fast(reduction)
  pad = repeat(" ", padfront)
  if latex 
    p = [reduction.sf_separation[i] ? latexify(string(reduction.problem.p[i]); env=:raw) : "\\varepsilon " * latexify(string(reduction.problem.p[i]); env=:raw) for i in eachindex(reduction.sf_separation)]
    str = pad * "\$\\left(" * join(p, ", ")  * "\\right)\$"
  else
    str = pad * "slow: " * join(string.(slow), ", ") * "\n" * pad * "fast: " * join(string.(fast), ", ")
  end
  return str
end


"""
    print_system([io::IO,] problem::ReductionProblem; latex::Bool=false)

Print input ODE system.
"""
function print_system(io::IO, problem::ReductionProblem; latex::Bool=false)
  print(io, _get_system_str(string.(problem.x), string.(problem.f); latex=latex))
end
print_system(problem::ReductionProblem; latex::Bool=false) = print_system(stdout, problem; latex)

function _get_system_str(x::Vector{String}, f::Vector{String}; padfront::Int=0, latex::Bool=false)
  dxdt = ["d$u/dt" for u in x]
  if latex
    str_out = "\\[\n\\begin{aligned}\n"
    for i in eachindex(f) 
      str_out *= latexify(dxdt[i]; env=:raw) * " &= " * latexify(f[i]; mult_symbol="", env=:raw) * "\\\\\n"
    end
    str_out *= "\\end{aligned}\n\\]"
    return str_out
  else
    str_out = ["" for _ in f]
    string_length = maximum(length.(dxdt))
    pad = repeat(" ", padfront)
    for i in eachindex(f) 
      str_out[i] = pad * padstring(dxdt[i], string_length; padfront=false) * " = " * f[i]
    end
    return join(str_out, "\n")
  end
end

