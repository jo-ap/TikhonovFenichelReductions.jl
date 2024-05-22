## Output functions

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

Print TFPV candidates to terminal or return LaTeXString that can be included in
tex file.
"""
function print_candidates(idx::Vector{Vector{Bool}}, prob::ReductionProblem; latex=false)
  if latex
    parameters = [latexify(πᵢ; env=:raw) for πᵢ in string.(prob.π)]
    m = length(parameters)
    str = "\\begin{array}{r$(repeat('c',length(string(prob.π))))} i & $(join(parameters, "^* & ")) \\\\ \n"
    for i in eachindex(idx)
      str *= "$i & " * join([idx[i][k] ? parameters[k] : "\\cdot" for k = 1:m], " & ") * " \\\\ \n"
    end
    str *= "\\end{array}"
    return latexstring(str)
  else
    subscripts = ["₀","₁","₂","₃","₄","₅","₆","₇","₈","₉"]
    numbers = string.(0:9)
    max_width = ndigits(length(idx)) 
    field_widths = [length(πᵢ) for πᵢ in string.(prob.π)]
    str = "π̃ $(repeat(" ", max_width-1)) = ($(join(string.(prob.π), ", ")))\n"
    str *= repeat("_", length(str)-2) * "\n"
    for i = 1:length(idx)
      idx_fast = idx[i]
      str *= "π" *
        join([subscripts[string(n) .== numbers][1] for n in string(i)], "") *
        "$(repeat(" ", max_width - ndigits(i))) = (" * 
        join(padstring.([idx_fast[k] ? string(prob.π[k]) : "0" for k = 1:length(idx_fast)], field_widths), ", ") *")\n"
    end
    print(str)
    return nothing
  end
end

"""
    $(TYPEDSIGNATURES)

Print generators of ideals in corresponding to irreducible components of varieties
for TFPV candidates to terminal (`V` as returned by `filter_dimension` or
`tfpv_candidates`).

See also: [`filter_dimension`](@ref), [`tfpv_candidates`](@ref)
"""
function print_varieties(V, prob::ReductionProblem; latex=false)
  if latex
    str_latex = "\\begin{array}{rl} i & \\mathcal{V}(f(\\cdot, \\tilde{\\pi_i})) \\\\ \n"
    for i = 1:length(V)
      Q = V[i]
      str = ["\\mathcal{V}(" * join(latexify.(string.(Q[i]); env=:raw, cdot=false), ", ") * ")" for i in eachindex(Q)]
      str_out = join(str, " \\cup ")
      str_latex *= "$i & " * str_out * "\\\\ \n"
    end
    str_latex *= "\\end{array}"
    return latexstring(str_latex)
  else
    subscripts = ["₀","₁","₂","₃","₄","₅","₆","₇","₈","₉"]
    numbers = string.(0:9)
    max_width = ndigits(length(V)) 
    for i = 1:length(V)
      V_str = ["𝑉(" * replace(join(string.(Y), ", "),  "=" => " = " ) * 
        ")" for Y∈V[i]]
      str = "V" * 
        join([subscripts[string(n) .== numbers][1] for n in string(i)], "") * 
        "$(repeat(" ", max_width - ndigits(i))) = " * 
        join(V_str, " ∪ ")  
      println(str)
    end
    return nothing
  end
end

function print_system(reduction::Reduction; latex=false)
  if latex
    latexstring("f(x, \\tilde{\\pi}) = $(latexify(string.(reduction.f⁰); env=:raw)) + \\varepsilon $(latexify(string.(reduction.f¹); env=:raw))")
  else
    "f(x,π̃) = (" * join(string.(reduction.f⁰), ", ") * ") + ε(" * join(string.(reduction.f¹), ", ") * ")"
  end
end 



