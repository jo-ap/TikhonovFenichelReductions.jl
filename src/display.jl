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
function show_slow_fast(reduction::Reduction)
  slow, fast = _get_slow_fast(reduction)
  println(" slow: " * join(string.(slow), ", "))
  println(" fast: " * join(string.(fast), ", "))
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
  println(io, " f: Function $(problem._f) from $(typeof(problem._f).name.module) defining ODE system") 
  for (x,dxdt) in zip(problem.x, problem.f)
    println(io, "   d$(string(x))/dt = " * string(dxdt))
  end
end

function Base.show(io::IO, ::MIME"text/plain", reduction::Reduction)
  println(io, "Reduction for dimension s = $(reduction.s) with")
  slow, fast = _get_slow_fast(reduction)
  println(io, "  slow: " * join(string.(slow), ", "))
  println(io, "  fast: " * join(string.(fast), ", "))
  if all(reduction.success)
    println(io, " M = [" * join(string.(reduction.M), ", ") * "]")
    println(io, " P = $(reduction.P)")
    println(io, " Ψ = [" * join([string(ψ) for ψ in reduction.Psi], ", ") * "]") # bug in Oscar? string.(reduction.Psi) throws error
    println(io, " Df(x_0) = $(reduction.Df_x0)")
  end
end

