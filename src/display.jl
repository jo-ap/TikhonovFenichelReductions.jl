
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
  # print(io, " f: Function $(problem._f) from $(typeof(problem._f).name.module) defining ODE system\n") 
  print(io, " ODE system:\n") 
  print(io, _get_system_str(string.(problem.x), string.(problem.f); padfront=2))
end

function Base.show(io::IO, ::MIME"text/plain", reduction::Reduction)
  println(io, "Reduction for dimension s = $(reduction.s) with")
  print(io, _get_slow_fast_str(reduction; padfront=1))
  if all(reduction.success)
    print(io, "\n M      = [" * join(string.(reduction.M), ", ") * "]")
    print(io, "\n P      = $(reduction.P)")
    print(io, "\n Ψ      = [" * join([string(ψ) for ψ in reduction.Psi], ", ") * "]") 
    print(io, "\n x₀     = [" * join(string.(reduction.x0), ", ") * "]")
    print(io, "\n Df(x₀) = $(reduction.Df0_at_x0)")
  end
  if any(reduction.reduction_cached)
    str = _get_reduced_system_str(reduction; padfront=2)
    str_reduced_system = reduction.reduction_cached[2] ?  "\n Reduced system:" :  "\n Formal reduced system:"
    print(io, str_reduced_system * "\n" * str)
  end
end

