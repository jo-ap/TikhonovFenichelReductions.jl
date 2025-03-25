## Convenience methods for cases with a large number of reductions 

# computing multiple reductions at once and 
"""
    $(TYPEDSIGNATURES)

Find reductions onto the same slow manifold and return their numeric index.

### Arguments
- `V`: `SlowManifold` object with Generators for the irreducible components of `V(f0)` for each slow-fast separation as returned by `tfpvs_and_manifolds`
- `Y`: Generators for one irreducible component that corresponds to the slow manifold

See also: [`tfpvs_and_manifolds`](@ref)
"""
function similar_reductions(V::Vector{Vector{SlowManifold}}, Y::SlowManifold)
  idx_component = Tuple[]
  for k in eachindex(V) 
    i = findfirst([Y.ideal == Vₖ.ideal for Vₖ in V[k]])
    if !isnothing(i)
      push!(idx_component, (k,i))
    end
  end
  return idx_component
end
  
"""
    $(TYPEDSIGNATURES)

Find all implicitly given manifolds (i.e. varieties) as returned by
`tfpvs_and_manifolds` for a given explicit (parametric) description `M` (same
form as the components of the system `problem.x`).

See also: [`compute_all_reductions`](@ref)
"""
function find_implicit_manifolds(
    problem::ReductionProblem, 
    manifolds::Vector{Vector{SlowManifold}},
    M::Union{Vector{QQMPolyRingElem}, Vector{FracFieldElem{QQMPolyRingElem}}}
  )
  @assert length(M) == length(problem.x) "`M` must have the same size as `problem.x`"
  F = parent(problem.x[1]//problem.x[2])
  v = F.([problem.p; M])
  idx = Tuple{Int, Int}[]
  for i in eachindex(manifolds)
    for j in eachindex(manifolds[i])
      if all([evaluate(f, v) == 0 for f in manifolds[i][j].gens_R])
        push!(idx, (i,j))
      end
    end
  end
  return idx
end

"""
    $(TYPEDSIGNATURES)

Compute reductions for TFPVs `sf_separations` onto each of the explicitly given
slow manifolds in `M`.
Returns the reductions `R` in the same format as `manifolds` and indices `idx`
such that `M[k]` corresponds to `manifolds[i][j]` and `R[i][j]` for each
`(i,j)` in `idx[k]`.

### Description
All possible choices of slow manifolds can be obtained with `unique_slow_manifolds`.
For each of these (that one is interested in), a parameterized form, i.e. the
slow manifold in phase space defined by the remaining components of the system,
must be provided.
The backwards step, i.e. finding all `manifolds` for which this is an explicit
description, is then done automatically with `find_implicit_manifolds`.

See also: [`unique_slow_manifolds`](@ref), [`compute_all_reductions`](@ref), [`compute_reduction!`](@ref), [`Reduction`](@ref), [`set_manifold!`](@ref), [`set_decomposition!`](@ref)
"""
function compute_all_reductions(
  problem::ReductionProblem,
  sf_separations::Vector{Vector{Bool}},
  manifolds::Vector{Vector{SlowManifold}},
  M::Union{Vector{Vector{QQMPolyRingElem}}, Vector{Vector{FracFieldElem{QQMPolyRingElem}}}}
)
  idx_M = [find_implicit_manifolds(problem, manifolds, _M) for _M in M]
  R = [[Reduction(problem, sf_separations[i]) for _ in 1:length(manifolds[i])] for i in eachindex(sf_separations)]
  for k in eachindex(M)
    println("M[$k] = (" * join(string.(M[k]), ", ") * ")\n")
    for (i,j) in idx_M[k]
      println("Reduction $i\n" * _get_slow_fast_str(R[i][j]; padfront=1))
      if set_manifold!(R[i][j], M[k])
        set_decomposition!(R[i][j], manifolds[i][j])
        if all(R[i][j].success)
          compute_reduction!(R[i][j]);
          println(_get_reduced_system_str(R[i][j]; padfront=1))
        end
      end
      println("")
    end
  end
  return R, idx_M
end

"""
    $(TYPEDSIGNATURES)

Obtain all possible choices of slow manifolds with dimension `s` given `V` and
`dim_V` as returned by `tfpvs_and_manifolds`.

See also: [`tfpvs_and_manifolds`](@ref)
"""
function unique_slow_manifolds(
    problem::ReductionProblem,
    V::Vector{Vector{SlowManifold}}
  )
  V_flat = vcat(V...)
  V_flat = V_flat[[v.dim == problem.s for v in V_flat]]
  return unique(v -> v.ideal, V_flat)
end


