import Base 

function _get_slow_fast(reduction::Reduction)
  p = reduction.p_sf
  sf_separation = reduction.sf_separation
  slow = p[.!sf_separation]
  fast = p[sf_separation]
  return slow, fast
end

"""
    $(TYPEDSIGNATURES)

Print slow and fast parameters.
"""
function show_slow_fast(reduction::Reduction; padfront=true)
  slow, fast = _get_slow_fast(reduction)
  pad = padfront ? " " : ""
  println(pad * "slow: " * join(string.(slow), ", "))
  println(pad * "fast: " * join(string.(fast), ", "))
end

function Base.show(io::IO, ::MIME"text/plain", problem::ReductionProblem)
  println(io, "ReductionProblemlem for dimension s = $(problem.s)")
  if !all(problem.idx_slow_fast)
    println(io, " x    = [" * join(string.(problem.x), ", ") * "]");
    println(io, " p    = [" * join(string.(problem.p), ", ") * "]");
    println(io, " p_sf = [" * join(string.(problem.p_sf), ", ") * "]");
  else
    println(io, " x = [" * join(string.(problem.x), ", ") * "]");
    println(io, " p = [" * join(string.(problem.p), ", ") * "]");
  end
  print(io, " f: Function $(problem._f) from $(typeof(problem._f).name.module) defining ODE system") 
  max_width = maximum(length.(string.(problem.x)))
  for (x,dxdt) in zip(problem.x, problem.f)
    print(io, "\n   " * padstring("d$x/dt", max_width + 4; padfront=false) * " = " * string(dxdt))
  end
end

function Base.show(io::IO, ::MIME"text/plain", reduction::Reduction)
  println(io, "Reduction for dimension s = $(reduction.s) with")
  slow, fast = _get_slow_fast(reduction)
  println(io, "  slow: " * join(string.(slow), ", "))
  print(io, "  fast: " * join(string.(fast), ", "))
  if all(reduction.success)
    print(io, "\n M = [" * join(string.(reduction.M), ", ") * "]")
    print(io, "\n P = $(reduction.P)")
    print(io, "\n Ψ = [" * join([string(ψ) for ψ in reduction.Psi], ", ") * "]") # bug in Oscar? string.(reduction.Psi) throws error
    print(io, "\n Df(x_0) = $(reduction.Df0_at_x0)")
  end
end

