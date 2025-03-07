## Output functions

boolean_to_numeric(idx_bool) = (1:length(idx_bool))[idx_bool]

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
    $(TYPEDSIGNATURES)

    Print slow-fast separations to terminal or use `latex=true` to print LaTeX
instead.
Optionally only print subset defined by (boolean or numeric) index set `idx`.
"""
function print_tfpv(
  io::IO,
  problem::ReductionProblem,
  sf_separations::Vector{Vector{Bool}}; 
  latex::Bool=false,
  idx::Union{AbstractVector{Int}, Vector{Bool}}=Int[])
  idx = idx == [] ? (1:length(sf_separations)) : idx
  _print_tfpv(io, problem, sf_separations, latex, idx)
end
function print_tfpv(
  problem::ReductionProblem,
  sf_separations::Vector{Vector{Bool}}; 
  latex::Bool=false,
  idx::Union{AbstractVector{Int}, Vector{Bool}}=Int[])
  print_tfpv(stdout, problem, sf_separations; latex=latex, idx=idx)
end

function _print_tfpv(
  io::IO, 
  problem::ReductionProblem,
  sf_separations::Vector{Vector{Bool}},
  latex::Bool,
  idx::AbstractVector{Int})
  if latex
    parameters = [latexify(p_sfᵢ; env=:raw) for p_sfᵢ in string.(problem.p_sf)]
    m = length(parameters)
    str = "\\begin{array}{r$(repeat('c',length(problem.p_sf)))} i & $(join(parameters .* "^\\star", " & ")) \\\\ \n"
    for i in idx
      str *= "$i & " * join([sf_separations[i][k] ? parameters[k] : "\\cdot" for k = 1:m], " & ") * " \\\\ \n"
    end
    str *= "\\end{array}"
  else
    subscripts = ["₀","₁","₂","₃","₄","₅","₆","₇","₈","₉"]
    numbers = string.(0:9)
    max_width = ndigits(length(sf_separations)) 
    field_widths = [length(p_sfᵢ) for p_sfᵢ in string.(problem.p_sf)]
    str = "π $(repeat(" ", max_width-1)) = ($(join(string.(problem.p_sf), ", ")))\n"
    str *= repeat("_", length(str)-2) * "\n"
    for i in idx
      idx_fast = sf_separations[i]
      str *= "π" *
      join([subscripts[string(n) .== numbers][1] for n in string(i)], "") *
      "$(repeat(" ", max_width - ndigits(i))) = (" * 
      join(padstring.([idx_fast[k] ? string(problem.p_sf[k]) : "0" for k = 1:length(idx_fast)], field_widths), ", ") *")\n"
    end
  end
  println(io, str)
end
function _print_tfpv(
  io::IO, 
  problem::ReductionProblem,
  sf_separations::Vector{Vector{Bool}},
  latex::Bool,
  idx::Vector{Bool})
  @assert length(sf_separations) == length(idx) "`sf_separations` and `idx` must have same length"
  _print_tfpv(io, problem, sf_separations, latex, boolean_to_numeric(idx))
end

"""
    $(TYPEDSIGNATURES)

Print generators of ideals corresponding to the irreducible components of
varieties `V(f0)` for TFPV candidates and their dimension (`V` and `dim_V` as
returned by or `tfpv_candidates`).
Use keyword argument `latex=true` to print LaTeX code instead.

See also: [`tfpv_candidates`](@ref), [`print_tfpv`](@ref), [`print_results`](@ref)
"""
function print_varieties(
  io::IO, 
  V::Vector{Vector{Vector{FracFieldElem{QQMPolyRingElem}}}},
  dim_V::Vector{Vector{Int}};
  latex::Bool=false,
  idx::Union{AbstractVector{Int}, Vector{Bool}}=Int[]) 
  idx = idx == [] ? (1:length(V)) : idx
  _print_varieties(io::IO, V, dim_V, latex, idx)
end
function print_varieties(
  V::Vector{Vector{Vector{FracFieldElem{QQMPolyRingElem}}}},
  dim_V::Vector{Vector{Int}};
  latex::Bool=false,
  idx::Union{AbstractVector{Int}, Vector{Bool}}=Int[]) 
  idx = idx == [] ? (1:length(V)) : idx
  _print_varieties(stdout, V, dim_V, latex, idx)
end

function _print_varieties(
  io::IO, 
  V::Vector{Vector{Vector{FracFieldElem{QQMPolyRingElem}}}},
  dim_V::Vector{Vector{Int}},
  latex::Bool,
  idx::AbstractVector{Int})
  if latex
    str = "\\begin{array}{rll} i & \\mathcal{V}(f(\\cdot, \\pi^\\star)) & \\dim_{\\textrm{Krull}} \\\\ \n"
    for i in idx
      str *= "$i"
      for j in eachindex(V[i])
        str *= " & \\langle " * join(latexify.(string.(V[i][j]); env=:raw, cdot=false), ", ") * " \\rangle & $(dim_V[i][j]) \\\\ \n"
      end
    end
    str *= "\\end{array}"
  else
    subscripts = ["₀","₁","₂","₃","₄","₅","₆","₇","₈","₉"]
    numbers = string.(0:9)
    max_width = ndigits(length(V)) 
    pad = "$(repeat(" ", max_width + 4))"
    for i in idx
      V_str = ["[" * replace(join(string.(V[i][j]), ", ")) * "], $(dim_V[i][j])" for j in eachindex(V[i])]
      str = "V" * 
        join([subscripts[string(n) .== numbers][1] for n in string(i)], "") * 
        "$(repeat(" ", max_width - ndigits(i))) : " * 
        join(V_str, "\n" * pad)
    end
  end
  println(io::IO, str)
end
function _print_varieties(
  io::IO,
  V::Vector{Vector{Vector{FracFieldElem{QQMPolyRingElem}}}},
  dim_V::Vector{Vector{Int}},
  latex::Bool,
  idx::Vector{Bool})
  idx = idx == [] ? (1:length(V)) : boolean_to_numeric(idx)
  println(_print_varieties(io, V, dim_V, latex, idx))
end

"""
    $(TYPEDSIGNATURES)

Print slow-fast separations that are TFPVs and generators for the corresponding
irreducible components of `V(f0)` together with their Krull dimension.

### Arguments 
- `problem`: `ReductionProblem` 
- `sf_separations`: Boolean indices defining all TFPVs `π⁺` that are slow-fast separations (0: slow, 1: fast).
- `V`: Generators for the irreducible component of the affine varietiy `V(f(⋅,π⁺))` for each slow-fast separation.
- `dim_V`: Krull dimensions of the ideals generated by the elements in `V` (this equals the dimension of `V(f(⋅,π⁺))` in the affine space `ℂⁿ` or its real part, given that it contains a real non-singular point).
- `idx`: (optional) index vector to include only certain TFPVs (boolean or numeric)

See also: [`tfpv_candidates`](@ref), [`print_tfpv`](@ref), [`print_varieties`](@ref)
"""
function print_results(
  io::IO,
  problem::ReductionProblem,
  sf_separations::Vector{Vector{Bool}},
  V,
  dim_V::Vector{Vector{Int}};
  idx::Union{AbstractVector{Int}, Vector{Bool}} = Int[])
  _print_results(io, problem, sf_separations, V, dim_V, idx)
end
function print_results(
  problem::ReductionProblem,
  sf_separations::Vector{Vector{Bool}},
  V,
  dim_V::Vector{Vector{Int}};
  idx::Union{AbstractVector{Int}, Vector{Bool}} = Int[])
  _print_results(stdout, problem, sf_separations, V, dim_V, idx)
end

function _print_results(
  io::IO,
  problem::ReductionProblem,
  sf_separations::Vector{Vector{Bool}},
  V,
  dim_V::Vector{Vector{Int}},
  idx::Union{AbstractVector{Int}, Vector{Bool}} = Int[])
  idx = idx == [] ? (1:length(sf_separations)) : idx
  for i in idx
    println(io, "$i")
    idx_fast = sf_separations[i]
    println(io, " π̃: [" * join([idx_fast[k] ? string(problem.p_sf[k]) : "0" for k = 1:length(idx_fast)], ", ") *"]")
    V_str = ["[" * replace(join(string.(V[i][k]), ", ")) * "], $(dim_V[i][k])" for k in eachindex(V[i])]
    println(io, " V: " * join(V_str, "\n    ") * "\n")
  end
end
function _print_results(
  io::IO,
  problem::ReductionProblem,
  sf_separations::Vector{Vector{Bool}},
  V,
  dim_V::Vector{Vector{Int}},
  idx::Vector{Bool})
  @assert length(sf_separations) == length(V) == length(dim_V) == length(idx) "`sf_separations`, `dim_V` and `V` all need to have the same length as `idx`"
  _print_results(io, problem, sf_separations, V, dim_V, boolean_to_numeric(idx))
end

"""
    $(TYPEDSIGNATURES)

Print the reduced system (after `compute_reduction!` has been called
successfully on the `Reduction` object). 
If the slow manifold was successfully specified, this returns the system in
coordinates on this manifold (i.e. only the `s`-dimensional system without the
reduced components is shown), if the reduction could be computed but the slow
manifold was not set successfully, the reduced system is shown in the original
phase space.

### Arguments 
- `reduction`: `Reduction` holding the reduced system 
- `rewrite`: Whether the RHS should be decomposed into polynomial and rational part (default is `true`), see `rewrite_rational`
- `factor`: Whether the polynomial parts should be factorized

See also: [`rewrite_rational`](@ref), [`print_tfpv`](@ref), [`print_varieties`](@ref)
"""
function print_reduced_system(io::IO, reduction::Reduction; rewrite::Bool=true, factor::Bool=false)
  println(io, _get_reduced_system_str(reduction; rewrite=rewrite, factor=factor, padfront=0))
end
function print_reduced_system(reduction::Reduction; rewrite::Bool=true, factor::Bool=false)
  println(stdout, _get_reduced_system_str(reduction; rewrite=rewrite, factor=factor, padfront=0))
end

function _get_reduced_system_str(reduction::Reduction; rewrite::Bool=true, factor::Bool=false, padfront::Int=0)
  if !any(reduction.reduction_cached)
    @warn "Reduced system is not yet computed"
    return ""
  end
  g = reduction.reduction_cached[2] ? reduction.g : reduction.g_raw
  x = reduction.reduction_cached[2] ? reduction.x[reduction.idx_components] : reduction.x
  dxdt = ["d$u/dt" for u in string.(x)]
  string_length = maximum(length.(dxdt))
  pad = repeat(" ", padfront)
  str_out = ["" for _ in g]
  if rewrite
    _g = rewrite_rational.(g; factor=factor)
    for i in eachindex(_g)
      h, r, q = _g[i]
      dudt = dxdt[i]
      str = pad * padstring(dudt, string_length; padfront=false) * " = "
      if r == 0 
        _h = h == 0 ? "0" : string(h)
        str_out[i] = str *  replace(_h, " * " => "*", "1 * " => "") 
      else
        _h = h == 0 ? "" : string(h)
        str_out[i] = str * replace(_h * " + (" * string(r) * ")//(" * string(q) * ")", " * " => "*", "1 * " => "") 
      end
    end
  else 
    for i in eachindex(g) 
      str_out[i] *= pad * padstring(dxdt[i], string_length; padfront=false) * " = " * string(g[i])
    end
  end
  return join(str_out, "\n")
end

"""
    $(TYPEDSIGNATURES)

Decompose a rational function `f = p/q` into polynomial and rational part, i.e.
return `f = h + r/q`, where `p,q,h,r` are polynomials.

### Arguments 
- `factor`: Indicate whether the RHS should be decomposed into polynomial and
rational part (default is `true`)

See also: [`print_reduced_system`](@ref)
"""
function rewrite_rational(term::FracFieldElem{QQMPolyRingElem}; factor::Bool=false)
  p = numerator(term)
  q = denominator(term)
  h,r = divrem(p,q)
  if factor 
    h = h == 0 ? h : Oscar.factor(h)
    r = r == 0 ? r : Oscar.factor(r)
  end
  return h, r, q
end

"""
    $(TYPEDSIGNATURES)

Print slow and fast parameters.
"""
function print_slow_fast(io::IO, reduction::Reduction)
  print(io, _get_slow_fast_str(reduction; padfront=0))
end
function print_slow_fast(reduction::Reduction)
  print(_get_slow_fast_str(reduction; padfront=0))
end

function _get_slow_fast_str(reduction::Reduction; padfront::Int=0)
  slow, fast = _get_slow_fast(reduction)
  pad = repeat(" ", padfront)
  str = pad * "slow: " * join(string.(slow), ", ") * "\n" * pad * "fast: " * join(string.(fast), ", ")
  return str
end


"""
    $(TYPEDSIGNATURES)

Print ODE system.
"""
function print_system(io::IO, problem::ReductionProblem; latex::Bool=false)
  print(io, _get_system_str(string.(problem.x), string.(problem.f); latex=latex))
end
print_system(problem::ReductionProblem; latex::Bool=false) = print_system(stdout, problem; latex)

function _get_system_str(x::Vector{String}, f::Vector{String}; padfront::Int=0, latex::Bool=false)
  dxdt = ["d$u/dt" for u in x]
  if latex
    str_out = "\\begin{align}\n"
    for i in eachindex(f) 
      str_out *= latexify(dxdt[i]; env=:raw) * " &= " * latexify(f[i]; env=:raw) * "\\\\\n"
    end
    str_out *= "\\end{align}"
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

