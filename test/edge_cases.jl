# System with constant RHS
x = ["u", "v"]
p = ["a", "b"]
function f(x,p)
  u,v = x 
  a,b = p
  return [a*u + b*v, b]
end
problem = ReductionProblem(f, x, p, 1)

@test isa(problem.f, Vector{QQMPolyRingElem})

tfpvs, varieties = tfpvs_and_varieties(problem)
@test tfpvs == [[true, false]]
@test isa(varieties, Vector{Vector{Variety}})

@test unique_varieties(problem, varieties) == varieties[1]
M = get_explicit_manifold(problem, varieties[1][1])
@test M[2]

reduction = Reduction(problem, tfpvs[1], varieties[1][1], M[1])
u,v = system_components(problem)
a,b = system_parameters(problem)
@test reduction.g == parent(u//v).([0, b])

function f(x,p)
  u,v = x 
  a,b = p
  return [a, b]
end
problem = ReductionProblem(f, x, p, 1)
tfpvs, varieties = tfpvs_and_varieties(problem)
@test tfpvs == Bool[]
@test varieties == Variety[]


