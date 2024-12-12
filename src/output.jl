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
function print_tfpv(problem::ReductionProblem,
                    sf_separations::Vector{Vector{Bool}}; 
                    latex::Bool=false,
                    idx::Union{AbstractVector{Int}, Vector{Bool}}=Int[])
  idx = idx == [] ? (1:length(sf_separations)) : idx
  _print_tfpv(problem, sf_separations, latex, idx)
end

function _print_tfpv(problem::ReductionProblem,
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
    println(str)
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
    print(str)
  end
  return
end
function _print_tfpv(problem::ReductionProblem,
                     sf_separations::Vector{Vector{Bool}},
                     latex::Bool,
                     idx::Vector{Bool})
  @assert length(sf_separations) == length(idx) "`sf_separations` and `idx` must have same length"
  _print_tfpv(problem, sf_separations, latex, boolean_to_numeric(idx))
end

"""
    $(TYPEDSIGNATURES)

Print generators of ideals corresponding to the irreducible components of
varieties `V(f0)` for TFPV candidates and their dimension (`V` and `dim_V` as
returned by or `tfpv_candidates`).
Use keyword argument `latex=true` to print LaTeX code instead.

See also: [`tfpv_candidates`](@ref), [`print_tfpv`](@ref), [`print_results`](@ref)
"""
function print_varieties(V::Vector{Vector{Vector{MPoly{RationalFunctionFieldElem{QQFieldElem, QQMPolyRingElem}}}}},
                         dim_V::Vector{Vector{Int}};
                         latex::Bool=false,
                         idx::Union{AbstractVector{Int}, Vector{Bool}}=Int[]) 
  idx = idx == [] ? (1:length(V)) : idx
  _print_varieties(V, dim_V, latex, idx)
end

function _print_varieties(V::Vector{Vector{Vector{MPoly{RationalFunctionFieldElem{QQFieldElem, QQMPolyRingElem}}}}},
                         dim_V::Vector{Vector{Int}},
                         latex::Bool,
                         idx::AbstractVector{Int})
  if latex
    str_latex = "\\begin{array}{rll} i & \\mathcal{V}(f(\\cdot, \\pi^\\star)) & \\dim_{\\textrm{Krull}} \\\\ \n"
    for i in idx
      str_latex *= "$i"
      for j in eachindex(V[i])
        str_latex *= " & \\langle " * join(latexify.(string.(V[i][j]); env=:raw, cdot=false), ", ") * " \\rangle & $(dim_V[i][j]) \\\\ \n"
      end
    end
    str_latex *= "\\end{array}"
    println(str_latex)
    return 
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
      println(str)
    end
    return
  end
end
function _print_varieties(V::Vector{Vector{Vector{MPoly{RationalFunctionFieldElem{QQFieldElem, QQMPolyRingElem}}}}},
                         dim_V::Vector{Vector{Int}},
                         latex::Bool,
                         idx::Vector{Bool})
  idx = idx == [] ? (1:length(V)) : boolean_to_numeric(idx)
  _print_varieties(V, dim_V, latex, idx)
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
  problem::ReductionProblem,
  sf_separations::Vector{Vector{Bool}},
  V,
  dim_V::Vector{Vector{Int}};
  idx::Union{AbstractVector{Int}, Vector{Bool}} = Int[])
  _print_results(problem, sf_separations, V, dim_V, idx)
end

function _print_results(
  problem::ReductionProblem,
  sf_separations::Vector{Vector{Bool}},
  V,
  dim_V::Vector{Vector{Int}},
  idx::Union{AbstractVector{Int}, Vector{Bool}} = Int[])
  idx = idx == [] ? (1:length(sf_separations)) : idx
  for i in idx
    println("$i")
    idx_fast = sf_separations[i]
    println(" π̃: [" * join([idx_fast[k] ? string(problem.p_sf[k]) : "0" for k = 1:length(idx_fast)], ", ") *"]")
    V_str = ["[" * replace(join(string.(V[i][k]), ", ")) * "], $(dim_V[i][k])" for k in eachindex(V[i])]
    println(" V: " * join(V_str, "\n    ") * "\n")
  end
end
function _print_results(
  problem::ReductionProblem,
  sf_separations::Vector{Vector{Bool}},
  V,
  dim_V::Vector{Vector{Int}},
  idx::Vector{Bool})
  @assert length(sf_separations) == length(V) == length(dim_V) == length(idx) "`sf_separations`, `dim_V` and `V` all need to have the same length as `idx`"
  _print_results(problem, sf_separations, V, dim_V, boolean_to_numeric(idx))
end


function print_system(reduction::Reduction; latex=false)
  if latex
    latexstring("f(x, \\tilde{\\pi}) = $(latexify(string.(reduction.f0); env=:raw)) + \\varepsilon $(latexify(string.(reduction.f1); env=:raw))")
  else
    "f(x,π̃) = (" * join(string.(reduction.f0), ", ") * ") + ε(" * join(string.(reduction.f1), ", ") * ")"
  end
end 

