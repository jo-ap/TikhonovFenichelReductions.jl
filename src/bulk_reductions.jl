## Convenience methods for cases with a large number of reductions 

# computing multiple reductions at once and 
"""
    $(TYPEDSIGNATURES)

Find reductions onto the same slow manifold and return their numeric index.

### Arguments
- `V`: `Variety` object with Generators for the irreducible components of `V(f0)` for each slow-fast separation as returned by `tfpvs_and_varieties`
- `Y`: Generators for one irreducible component that corresponds to the slow manifold

See also: [`tfpvs_and_varieties`](@ref)
"""
function similar_reductions(V::Vector{Vector{Variety}}, Y::Variety)
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

Find all varieties in `tfpvs_and_varieties` that are equal to `V`.  

See also: [`compute_reductions`](@ref)
"""
function find_varieties(varieties::Vector{Vector{Variety}}, V::Variety)
  idx = Tuple{Int, Int}[]
  for i in eachindex(varieties)
    for j in eachindex(varieties[i])
      if V.ideal == varieties[i][j].ideal
        push!(idx, (i,j))
      end
    end
  end
  return idx
end

"""
    $(TYPEDSIGNATURES)

Compute all reductions onto each of the explicitly given slow manifolds in `M`
with corresponding varieties `V`.
Returns the reductions `R` as a dictionary with keys `(i,j)`, such that
`R[(i,j)]` is the reduction for `tfpv[i]` and the variety `V[i][j]`.

### Description
All possible choices of slow manifolds can be obtained with
[`unique_varieties`](@ref).
For each of these, a parameterized form, i.e. the slow manifold in phase space
defined by the remaining components of the system, must be provided.
In most cases this can be obtained with [`get_explicit_manifold`](@ref).
The backwards step, i.e. finding all `manifolds` for which this is an explicit
description, is then with `find_varieties`.

!!! warning
    If `set_decomposition!` fails, the reduction cannot be computed. This triggers a warning, but the `Reduction` instance will be returned (without containing the reduced system).

See also: [`unique_varieties`](@ref), [`compute_reduction!`](@ref), [`Reduction`](@ref), [`set_manifold!`](@ref), [`set_decomposition!`](@ref)
"""
function compute_reductions(
    problem::ReductionProblem,
    tfpvs::Vector{Vector{Bool}},
    varieties::Vector{Vector{Variety}},
    V::Vector{Variety},
    M::Union{Vector{Vector{QQMPolyRingElem}}, Vector{Vector{FracFieldElem{QQMPolyRingElem}}}}
  )
  # disable all info or debug messages
  logger = ConsoleLogger(stderr, Logging.Warn)
  d = Dict{Tuple{Int,Int}, Reduction}()
  idx_V = [find_varieties(varieties, v) for v in V]
  with_logger(logger) do 
    for k in eachindex(idx_V)
      for (i,j) in idx_V[k]
        d[(i,j)] = Reduction(problem, tfpvs[i], varieties[i][j], M[k]) 
        if any(.!d[(i,j)].reduction_cached) && !d[(i,j)].no_reduction
          @warn "Reduction $i.$j could not be computed"
        end
      end
    end
  end
  return d, idx_V
end


"""
    $(TYPEDSIGNATURES)

Obtain all possible choices of slow manifolds with dimension `s` given the
manifolds as returned by `tfpvs_and_varieties`.

See also: [`tfpvs_and_varieties`](@ref)
"""
function unique_varieties(
    problem::ReductionProblem,
    V::Vector{Vector{Variety}}
  )
  V_flat = vcat(V...)
  V_flat = V_flat[[v.dim == problem.s for v in V_flat]]
  return unique(v -> v.ideal, V_flat)
end


