import Base 

"""
    $(TYPEDSIGNATURES)

Print slow and fast parameters.
"""
function show_slow_fast(reduction::Reduction)
  p = reduction.p_sf
  sf_separation = reduction.sf_separation
  println(" slow: " * join(string.(p[.!sf_separation]), ", "))
  println(" fast: " * join(string.(p[sf_separation]), ", "))
end

function Base.show(io::IO, ::MIME"text/plain", problem::ReductionProblem)
  println("ReductionProblemlem for dimension s = $(problem.s)")
  if !all(problem.idx_slow_fast)
    println(" x    = [" * join(string.(problem.x), ", ") * "]");
    println(" p    = [" * join(string.(problem.p), ", ") * "]");
    println(" p_sf = [" * join(string.(problem.p_sf), ", ") * "]");
  else
    println(" x = [" * join(string.(problem.x), ", ") * "]");
    println(" p = [" * join(string.(problem.p), ", ") * "]");
  end
  println(" f: Function $(problem._f) from $(typeof(problem._f).name.module) defining ODE system") 
  for (x,dxdt) in zip(problem.x, problem.f)
    println("   d$(string(x))/dt = " * string(dxdt))
  end
end

function Base.show(io::IO, ::MIME"text/plain", reduction::Reduction)
  println("Reduction for dimension s = $(reduction.s) with")
  show_slow_fast(reduction)
  # p = reduction.p_sf
  # sf_separation = reduction.sf_separation
  # println(" rates:")
  # println("  slow: " * join(string.(p[.!sf_separation]), ", "))
  # println("  fast: " * join(string.(p[sf_separation]), ", "))
  if all(reduction.success)
    println(" M = [" * join(string.(reduction.M), ", ") * "]")
    println(" P = $(reduction.P)")
    println(" Ψ = [" * join([string(ψ) for ψ in reduction.Psi], ", ") * "]") # bug in Oscar? string.(reduction.Psi) throws error
    println(" Df(x_0) = $(reduction.Df_x0)")
  end
end

