
"""
Type that defines a Tikhonov-Fenichel reduction problem, i.e. a polynomial ODE
system, for which all slow-fast separations of rates yielding a reduction onto
dimension `s` are considered. 

The type `QQMPolyRingElem` is used in [Oscar.jl](https://www.oscar-system.org/)
to represent elements of a polynomial ring; here we this is `ℚ[p,x]`.

### Fields
    $(TYPEDFIELDS)

"""
struct ReductionProblem 
  "RHS of ODE system as a vector of polynomials"
  f::Vector{QQMPolyRingElem}
  "Vector of dynamic variables"
  x::Vector{QQMPolyRingElem}
  "Vector of all parameters"
  p::Vector{QQMPolyRingElem}
  "Dimension of reduced system"
  s::Int
  "Vector of parameters that may be small (all others are considered fixed)"
  p_sf::Vector{QQMPolyRingElem}
  "Boolean index, such that `p_sf=p[idx_slow_fast]`"
  idx_slow_fast::Vector{Bool}
  "RHS of ODE system as a Julia function with signature `_f(x,p)`"
  _f::Function
end

"""
Type that holds information on the slow manifold as an irreducible component of
the variety `V(f0)`.

    $(TYPEDFIELDS)
"""
struct SlowManifold
  "associated ideal in Rx"
  ideal::MPolyIdeal
  "generators of associated ideal parsed to R"
  gens_R::Vector{QQMPolyRingElem}
  "groebner basis of associated ideal in Rx"
  groebner_basis::Vector{MPoly{RationalFunctionFieldElem{QQFieldElem, QQMPolyRingElem}}}
  "groebner basis of associated ideal parsed to R"
  groebner_basis_R::Vector{QQMPolyRingElem}
  "Matrix T, such that `gens(I)*T=groebner_basis`"
  T::MatSpaceElem{MPoly{RationalFunctionFieldElem{QQFieldElem, QQMPolyRingElem}}}
  "Krull dimension of associated ideal"
  dim::Int64
end

function SlowManifold(I::MPolyIdeal, dim::Int64, problem::ReductionProblem) 
  _I = ideal(rewrite_variety_generator.(gens(I)))
  gens_R = [parse_variety_generator(problem, rewrite_variety_generator(g)) for g in gens(_I)]
  @assert I == _I
  G, _T = groebner_basis_with_transformation_matrix(_I; complete_reduction=true)
  G = gens(G)
  # cancel_monomial_leading_coefficients!(G, _T)
  T = transpose(_T)
  gb_R = [parse_variety_generator(problem, rewrite_variety_generator(g)) for g in G]
  return SlowManifold(_I, gens_R, G, gb_R, T, dim)
end

# function cancel_monomial_leading_coefficients!(G::Vector{<:MPolyRingElem}, T::MatSpaceElem)
#   lc = leading_coefficient.(G)
#   for j in eachindex(lc)
#     if !isone(lc[j]) && is_monomial(G[j])
#       for i in 1:size(T,1)
#         T[i,j] *= 1//lc[j]
#       end
#       G[j] *= 1//lc[j]
#     end
#   end
#   return nothing
# end

"""
    $(TYPEDSIGNATURES)

Constructor for `ReductionProblem` Type.

### Arguments 
- `f(x,p)::Function`: Julia function defining the RHS of the ODE system 
- `x::Vector{String}`: Vector of dynamic variables 
- `p::Vector{String}`: Vector of all parameters
- `s::Int`: Dimension of reduced system 
- `idx_slow_fast::Vector{Bool}`: Boolean index for all rates that are either small or large (all others are considered fixed)

### Description
This function is used to set up the problem of finding Tikhonov-Fenichel
Parameter Values for dimension `s`.
The names of all variables and parameters in the system are parsed to
appropriate types in
[Oscar.jl](https://www.oscar-system.org/),
so that the necessary conditions for the existence of a reduction onto an
`s`-dimensional slow manifold can be evaluated.

See also: [`tfpvs_and_manifolds`](@ref), [`tfpvs_groebner`](@ref)
"""
function ReductionProblem(
  f::Function, 
  x::Vector{T}, 
  p::Vector{T}, 
  s::Int; 
  idx_slow_fast::Vector{Bool}=[true for _ in p]) where T<:Union{String,Symbol}
  @assert s < length(x) "the dimension of the reduced system must be smaller than that of the full system, i.e. s < n"
  _, _x, _p, _f = parse_system(f, x, p)
  _p_sf = _p[idx_slow_fast]
  ReductionProblem(_f, _x, _p, s, _p_sf, idx_slow_fast, f)
end

"""
    $(TYPEDSIGNATURES)

Convenience function to get components of ODE system from instance of `ReductionProblem`. 

See also: [`ReductionProblem`](@ref), [`system_parameters`](@ref)
"""
system_components(problem::ReductionProblem) = problem.x

"""
    $(TYPEDSIGNATURES)

Convenience function to get parameters of ODE system from instance of `ReductionProblem`. 

See also: [`ReductionProblem`](@ref), [`system_components`](@ref)
"""
system_parameters(problem::ReductionProblem) = problem.p

## Helper Functions
"""
    $(TYPEDSIGNATURES) 

Parse dynamic variables `x` and parameters `p` of polynomial ODE system so that
they can be used with Oscar.jl. Return the polynomial Ring `R = ℚ[x,p]`
together with `x`, `p` and `f` parsed to the appropriate OSCAR types.
"""
function parse_system(f::Function, x::Vector{String}, p::Vector{String})
  R, v = polynomial_ring(QQ, [p..., x...])
  _p = v[1:length(p)]
  _x = v[length(p)+1:end]
  _f = f(_x, _p)
  return R, _x, _p, _f
end

"""
    $(TYPEDSIGNATURES) 

Parse a polynomial in the ring `Rx = QQ(p)[x]` into `R = QQ[x,p]`.

See also: [`rewrite_variety_generator`](@ref)
"""
function parse_variety_generator(problem::ReductionProblem, Y::MPoly{RationalFunctionFieldElem{QQFieldElem,QQMPolyRingElem}})
  ce = coefficients_and_exponents(Y)
  s = parent(problem.f[1])(0)
  for (c,α) in ce 
    p = evaluate(numerator(c), problem.p)
    q = denominator(c)
    @assert q == 1 "Generators for slow manifold are not preprocessed correctly"
    s += p*prod(problem.x.^α)
  end
  return s
end

"""
    $(TYPEDSIGNATURES) 

Rewrite a generating polynomial of a variety as a subset in the phase space,
i.e. a polynomial in `Rx = QQ(p)[x]`, by multiplying with parameters occuring
in the denominator.
Thus, the resulting polynomial can be parsed to the ring `R = QQ[p,x]`.

See also: [`parse_variety_generator`](@ref)
"""
function rewrite_variety_generator(p::MPoly{RationalFunctionFieldElem{QQFieldElem,QQMPolyRingElem}})
  while true
    c = coefficients(p)
    for _c in c
      denom = denominator(_c)
      if denom != 1 
        p = denom * p 
        break
      end
    end
    break
  end
  return p
end

"""
    $(TYPEDSIGNATURES)

Convert integer `i` into a boolean vector of length `n` representing `i` in
binary.
"""
function num2bin(i::Int, n::Int)
  sf_separation = bitstring(i)[(end-n+1):end]
  sf_separation = [i == '1' for i in sf_separation]
end

"""
    $(TYPEDSIGNATURES)

Compute determinants of all possible k×k minors of quadratic matrix `M` for
k > `r`.
"""
function get_determinants(M::MatSpaceElem{QQMPolyRingElem}, r::Int)
  n = size(M)[1]
  @assert n == size(M)[2] "M must be a quadratic matrix."
  # get all valid combinations of rows and columns
  combinations = num2bin.(1:2^n, n)
  combinations = combinations[sum.(combinations) .> r]
  l_s = sum([binomial(n, k)^2 for k = (r+1):(n-1)])
  d = Vector{QQMPolyRingElem}(undef, l_s+1)
  idx = collect(1:n)
  i = 1
  for row in combinations 
    for col in combinations[sum.(combinations) .== sum(row)]
      d[i] = det(M[idx[row], idx[col]])
      i += 1
    end
  end
  return d
end

"""
    $(TYPEDSIGNATURES)

Check if all polynomials in `F` vanish for all parameters in `p_sf` set to zero.
"""
function allvanish(F::Vector{QQMPolyRingElem}, p_sf::Vector{QQMPolyRingElem})
  z = zeros(parent(p_sf[1]), length(p_sf))
  for f∈F
    if !iszero(evaluate(f, p_sf, z))
      return false
    end
  end
  return true
end

"""
    $(TYPEDSIGNATURES)

Use the necessary conditions imposed by the RHS of the polynomial ODE System
and its Jacobian to filter out TFPV candidates.
All TFPVs must result in the vanishing of the returned set of polynomials.

### Description
If `π⁺` is a TFPV yielding a reduction onto an `s`-dimensional slow manifold,
there exist a point `x₀`, such that 
- `f(x₀,π⁺)=0`
- for any `k>s` the determinants of all `k×k` minors of `D₁f(x₀,π⁺)` vanish
These properties can be used to filter out possible TFPV candidates.

We are interested in partial solutions of the system of polynomials defined by
the conditions above.
In particular, we only consider conditions on the parameters, since there might
be multiple slow manifolds and reductions.
Thus, we eliminate the dynamic variables from the ideal generated by the
polynomial conditions above.
This function computes a generating set for this elimination ideal and all
TFPVs lie in its vanishing set.

See also: [`tfpvs_and_manifolds`](@ref) 
"""
function tfpvs_groebner(problem::ReductionProblem)
  # number of dimensions to reduce the system by
  r = length(problem.x) - problem.s
  # Jacobian of the full system
  J = jacobian(problem.f, problem.x)
  # determinants of all k×k minors of J for k>r
  d = get_determinants(J, r)
  # All polynomials that generate the ideal used to determine TFPVs
  poly_gens = [problem.f; d]
  # Build ideal from polynomial expressions 
  I = ideal(parent(problem.f[1]), poly_gens)
  # eliminate dynamic variables
  Iₓ = eliminate(I, problem.x)
  # Generating set for Iₓ
  G = groebner_basis(Iₓ; complete_reduction=true)
  G = gens(G)
  return G
end

## Filter functions
# """
#     $(TYPEDSIGNATURES)

# Use the Krull dimension of the ideal generated by the unperturbed part of the
# system to filter out TFPV candidates form all possible slow-fast separations. 

# ### Description
# The slow manifold on which the reduced system is defined is contained in
# `V(f0)`, the affine variety of `f0`, i.e. the zero set of the fast /
# unperturbed part of the system (when we consider the entries as polynomials in
# the dynamic varieties `x`, so that the variety is a subset of the phase space).
# Thus, `V(f0)` needs to have an irreducible component with dimension `s`, which
# we can check using the Krull dimension of the corresponding ideal.

# By default, all 2ᵐ-2 possible slow-fast separations of the m parameters are
# checked, but if `idx::Vector{Vector{Bool}}` is defined, only those candidates
# are considered (0: slow variable, 1: fast variable).

# If `compute_primary_decomposition=true` (default behaviour) is set, the
# function attempts to compute a primary decomposition of the ideal corresponding
# to the unperturbed part of the system `f0`.

# If in addition, `exact_dimension=true` (default), only those candidates are
# kept for which the affine variety has exactly dimension `problem.s`. Keep in
# mind that this is the Krull dimension, so the topological dimension over the
# reals can be smaller. However, the existence of a non-singular point in the
# irreducible component implies that topological and Krull dimension are equal. 
# Since we require such a point later, we may use the exact dimension already. 

# See also: [`tfpvs_groebner`](@ref)
# """
# function dimension_criterion(
#   problem::ReductionProblem; 
#   idx::Vector{Vector{Bool}}=[Bool[]],
#   compute_primary_decomposition::Bool=true, 
#   exact_dimension::Bool=true)
#   # redefine RHS of ode system: Interpret parameters as coefficients and only
#   # use dynamic variables as variables for the polynomial ring
#   x_str = string.(problem.x)
#   p_str = string.(problem.p)
#   K, p = rational_function_field(QQ, p_str)
#   R, x = polynomial_ring(K, x_str)
#   p_sf = p[problem.idx_slow_fast]
#   # filter TFPV candidates 
#   idx = length(idx[1]) == 0 ? num2bin.(1:(2^length(problem.p_sf)-2), length(problem.p_sf)) : idx
#   idx_candidates = zeros(Bool, length(idx))
#   V = compute_primary_decomposition ? [] : nothing
#   dim_components = compute_primary_decomposition ? Vector{Vector{Int}}() : nothing
#   cnt = 1
#   for i in idx    
#     # Get unperturbed part of system (fast part)
#     _p = p
#     _p[problem.idx_slow_fast] = p_sf .* i 
#     f0 = problem._f(x, _p)
#     # Check if Krull dimension is at least s
#     I = ideal(R, f0)
#     if dim(I) >= problem.s
#       # compute the irreducible components of V(f0)
#       if compute_primary_decomposition 
#         PD = primary_decomposition(I)
#         Y = [gens(groebner_basis(Q[2]; complete_reduction=true)) for Q in PD]
#         Y_dim = dim.([ideal(R, Yᵢ) for Yᵢ in Y])
#         # check if there is an irreducible component of the affine variety
#         # V(f0) with dimension exactly s
#         if exact_dimension && any(Y_dim .== problem.s)
#           push!(V, Y)
#           push!(dim_components, Y_dim)
#           idx_candidates[cnt] = true
#         end
#       # do not compute the irreducible components of V(f0)
#       else
#         idx_candidates[cnt] = true
#       end
#     end
#     cnt += 1
#   end
#   # return list with slow-fast separation candidates as defined by boolean indices
#   _idx = idx[idx_candidates]
#   return _idx, (V, dim_components)
# end


"""
    $(TYPEDSIGNATURES)

Return parameter vector where all slow rates are set to zero.
"""
function get_tfpv(p, idx_slow_fast, sf_separation)
  _p = copy(p)
  _p[idx_slow_fast] = _p[idx_slow_fast].*sf_separation
  return _p
end

function update_cofficients(f, p)
  if f == 0
    return f 
  else
    coeffs = [c for c in coefficients(f)]
    mons = [m for m in monomials(f)]
    @assert sum(coeffs .* mons) == f
    coeffs_new = [evaluate(c, p) for c in coeffs]
    return sum(coeffs_new .* mons)
  end
end

"""
    $(TYPEDSIGNATURES)

Find all slow-fast separations `π⁺` that are TFPVs by using the necessary conditions
- the affine variety `V(f0)` contains an irreducible component `Y` of dimension `s` 
- the `s`-th coefficient of the characteristic polynomial of `D₁f(x,π⁺)` is non-zero for `x∈Y`

### Description
The irreducible components are obtained by computing a minimal primary decomposition. 
The Jacobian at a point in an irreducible component `Y` is constructed
symbolically by computing normal forms with respect to a Gröbner basis `G`, s.t.
`V(G)=Y`. 

To obtain all general TFPVs and not just slow-fast separations, one can use the
function `tfpvs_groebner`.

See also: [`tfpvs_groebner`](@ref), [`print_results`](@ref), [`print_tfpvs`](@ref), [`print_slow_manifolds`](@ref)
"""
function tfpvs_and_manifolds(problem::ReductionProblem)
  # check all possible slow-fast separations for sufficient conditions to be a TFPV for dimension s
  slow_fast = num2bin.(1:(2^length(problem.p_sf)-2), length(problem.p_sf)) 
  # define all polynomials and the Jacobian of f in ℝ(p_sf)[x]
  F, _p = rational_function_field(QQ, string.(problem.p))
  _, _x = polynomial_ring(F, string.(problem.x))
  f = problem._f(_x, _p)
  J = matrix(parent(f[1]), [[derivative(fᵢ, xᵢ) for xᵢ in _x] for fᵢ in f])
  # keep track of which slow-fast separation satisfies the conditions and save
  # irreducible components of V(f0) and their dimensions
  idx_keep = zeros(Bool, length(slow_fast))
  slow_manifolds = Vector{Vector{SlowManifold}}(undef, length(slow_fast))
  for i in eachindex(slow_fast)
    keep_i, Y, dim_Y = get_slow_manifolds(problem, _x, _p, slow_fast[i], J)
    if any(keep_i)
      idx_keep[i] = true
      slow_manifolds[i] = [SlowManifold(Y[k], dim_Y[k], problem) for k in eachindex(Y)]
    end
  end
  return slow_fast[idx_keep], slow_manifolds[idx_keep]
end

function get_slow_manifolds(problem, _x, _p, slow_fast, J)
  tfpv_candidate = get_tfpv(_p, problem.idx_slow_fast, slow_fast)
  f0 = problem._f(_x, tfpv_candidate)
  I = ideal(f0)
  if dim(I) < problem.s
    return [false], nothing, nothing 
  else
    PD = primary_decomposition(I)
    Y = [Q[2] for Q in PD]
    dim_Y = dim.(Y)
    # set the slow parameters to zero in Jacobian
    J_p_sf = map(f -> update_cofficients(f, tfpv_candidate), J)
    keep_i = [false for _ in Y]
    for j in eachindex(Y) 
      # check if dimension of irreducible component Yⱼ is as desired
      if dim_Y[j] == problem.s
        # substitute x = x0 ∈ Yⱼ
        # we need that Y[j] is the radical of Qᵢ in order for normal forms work !
        M = map(f -> normal_form(f, Y[j]), J_p_sf)
        # Let Χ(τ) = τⁿ+ σₙ₋₁(x,p_sf)τ⁽ⁿ⁻¹⁾ + … + σ₁(x,p_sf)τ + σ₀(x,p_sf) be the characteristic polynomial 
        # for x0 ∈ Yⱼ and p_sf⁺ a TFPV for dimension s we have σₛ(x0,p_sf⁺) ≠ 0 
        keep_i[j] = coeff(charpoly(M), problem.s) != 0
      end 
    end
  end
  return keep_i, Y, dim_Y
end
