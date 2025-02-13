## Convenience methods for cases with a large number of reductions 

# computing multiple reductions at once and 
"""
    $(TYPEDSIGNATURES)

Find reductions onto the same slow manifold and return their numeric index.

### Arguments
- `V`: Generators for the irreducible components of `V(f0)` for each slow-fast separation as returned by `tfpv_candidates`
- `Y`: Generators for one irreducible component that corresponds to the slow manifold

See also: [`tfpv_candidates`](@ref)
"""
function similar_reductions(V::Vector{Vector{Vector{FracFieldElem{QQMPolyRingElem}}}},
                            Y::Vector{FracFieldElem{QQMPolyRingElem}})
  idx_bool = [false for _ in V]
  for k in eachindex(V) 
    idx_bool[k] = any([Y == Vₖ for Vₖ in V[k]])
  end
  idx_num = (1:length(V))[idx_bool]
  return idx_num
end

"""
    $(TYPEDSIGNATURES)

Compute reductions for TFPVs `sf_separations[idx]` onto slow manifold `M`.

### Description
`Psi` for the product decomposition `f0=P⋅Psi` and the slow manifold `M` depend
on the irreducible component `Y` of `V(f0)`. 
Thus, for a particular choice of `Y`, we can compute all reduced systems for
TFPVs for which `Y` is an irreducible component of `V(f0)` at once (the same
choice for `Psi` and `M` works in all those cases).
This function is essentially a convenience wrapper calling `set_manifold!`,
`set_decomposition!` and `compute_reduction` to speed up the process of
computing reduced systems onto the same slow manifold.

All possible choices of slow manifolds can be obtained with `unique_slow_manifolds`.
To find all TFPVs that have slow manifolds in common use `similar_reductions`.

See also: [`unique_slow_manifolds`](@ref), [`similar_reductions`](@ref), [`compute_reduction`](@ref), [`Reduction`](@ref), [`set_manifold!`](@ref), [`set_decomposition!`](@ref)
"""
function compute_bulk_reductions(problem::ReductionProblem,
                                 sf_separation::Vector{Vector{Bool}},
                                 idx::Union{UnitRange{Int}, Vector{Int}},
                                 Psi,
                                 M::VecOrMat;
                                 idx_components::Vector{Int} = Int[])
  G = Tuple{Int,Tuple{Reduction, Vector{FracFieldElem{QQMPolyRingElem}}}}[]
  idx_components = (idx_components == []) ? (1:length(problem.x)) : idx_components
  cnt = 1
  for i in idx
    red_i = Reduction(problem, sf_separation[i]);
    println("Reduction $i")
    show_slow_fast(red_i)
    if set_manifold!(red_i, M)
      set_decomposition!(red_i, Psi)
      if all(red_i.success)
        _, g_i = compute_reduction(red_i);
        display(g_i[idx_components])
        push!(G, (i, (red_i, g_i)))
      end
    end
    cnt += 1
    println("")
  end
  G = OrderedDict([i => rg_i for (i,rg_i) in G])
  return G
end

"""
    $(TYPEDSIGNATURES)

Obtain all possible choices of slow manifolds with dimension `s` given `V` and
`dim_V` as returned by `tfpv_candidates`.

See also: [`tfpv_candidates`](@ref)
"""
function unique_slow_manifolds(problem::ReductionProblem,
                               V::Vector{Vector{Vector{FracFieldElem{QQMPolyRingElem}}}},
                               dim_V::Vector{Vector{Int}})
  V_flat = vcat(V...)
  dim_flat = vcat(dim_V...)
  V_flat_all_dim = V_flat[dim_flat .== problem.s]
  return unique(V_flat_all_dim)
end

# function components_not_zero(V::Vector{Vector{Vector{MPoly{RationalFunctionFieldElem{QQFieldElem, QQMPolyRingElem}}}}}, sf_separations::Vector{Int})
#   I = ideal.(V)
#   R = parent(V[1][1])
#   u = gens(R)
#   keep = [false for _ in V]
#   for i in eachindex(I)
#     keep[i] = !any([ideal_membership(x, I[i]) for x in u[sf_separations]])
#   end
#   return any.(keep)
# end

