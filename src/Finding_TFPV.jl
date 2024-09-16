"""
    $(TYPEDEF)

Type that defines a Tikhonov-Fenichel reduction problem, i.e. a polynomial ODE
system, for which all slow-fast separations of rates yielding a reduction onto
dimension `s` are considered. 

### Fields 
- `f::Vector{QQMPolyRingElem}`: RHS of ODE system as a vector of polynomials
- `x::Vector{QQMPolyRingElem}`: Vector of dynamic variables
- `θ::Vector{QQMPolyRingElem}`: Vector of all parameters
- `s::Int`: Dimension of reduced system
- `π::Vector{QQMPolyRingElem}`: Vector of parameters to be considered slow or fast
- `idx_slow_fast::Vector{Bool}`: Boolean index, such that `π=ϑ[idx_slow_fast]`
- `J::AbstractAlgebra.Generic.MatSpaceElem{QQMPolyRingElem}`: Jacobian of `f`
- `_f::Function`: RHS of ODE system as a Julia function with arguments `x` and `ϑ`

The type `QQMPolyRingElem` is used in [Oscar.jl](https://www.oscar-system.org/)
to represent elements of a polynomial ring; here this is `ℝ[x,θ]`.
"""
mutable struct ReductionProblem 
  f::Vector{QQMPolyRingElem}
  x::Vector{QQMPolyRingElem}
  θ::Vector{QQMPolyRingElem}
  s::Int
  π::Vector{QQMPolyRingElem}
  idx_slow_fast::Vector{Bool}
  J::AbstractAlgebra.Generic.MatSpaceElem{QQMPolyRingElem}
  _f::Function
end

"""
    $(TYPEDSIGNATURES)

Constructor for `ReductionProblem` Type.

### Arguments 
- `f(x,θ)::Function`: Julia function defining the RHS of system 
- `x::Vector{String}`: Vector of dynamic variables 
- `θ::Vector{String}`: Vector of all parameters
- `s::Int`: Dimension of reduced system 
- `idx_slow_fast::Vector{Bool}`: Boolean index for all rates that are either small or large (all others are considered fixed)

### Description
This function is used to setup the problem of finding Tikhonov-Fenichel
Parameter Values, i.e. slow-fast separations of rates that yield a reduction
onto an `s`-dimensional system. 
The names of all variables and parameters in the system are parsed to
appropriate types in
[Oscar.jl](https://www.oscar-system.org/),
so that the necessary conditions for the existence of such a reduction can be
evaluated.

See also: [`tfpv_candidates`](@ref)
"""
function ReductionProblem(
  f::Function, 
  x::Vector{String}, 
  θ::Vector{String}, 
  s::Int; 
  idx_slow_fast::Vector{Bool}=Bool[])
  @assert s < length(x) "the dimension of the reduced system must be smaller than that of the full system, i.e. s < n"
  _, _x, _θ, _f = parse_system(f, x, θ)
  idx_slow_fast = length(idx_slow_fast) > 0 ? idx_slow_fast : [true for i in 1:length(θ)]   
  _π = _θ[idx_slow_fast]
  J = jacobian(_f, _x)
  ReductionProblem(_f, _x, _θ, s, _π, idx_slow_fast, J, f)
end

## Helper Functions
"""
    $(TYPEDSIGNATURES) 

Parse dynamic variables `x` and parameters `θ` of polynomial ODE system so that
they can be used with Oscar.jl. Return the polynomial Ring `R = ℚ[x,θ]`
together with `x`, `θ` and `f` parsed to the appropriate OSCAR types.
"""
function parse_system(f::Function, x::Vector{String}, θ::Vector{String})
  R, v = polynomial_ring(QQ, [x..., θ...])
  _x = v[1:length(x)]
  _θ = v[length(x)+1:end]
  _f = f(_x, _θ)
  return R, _x, _θ, _f
end

"""
    $(TYPEDSIGNATURES)

Convert integer `i` into a boolean vector of length `n` representing `i` in
binary.
"""
function num2bin(i::Int, n::Int)
  idx = bitstring(i)[(end-n+1):end]
  idx = [i == '1' for i in idx]
end

"""
    $(TYPEDSIGNATURES)

Compute determinants of all possible k×k minors of quadratic matrix `M` for
k > `r`.
"""
function get_determinants(M::AbstractAlgebra.Generic.MatSpaceElem{QQMPolyRingElem}, r::Int)
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

Check if all polynomials in `F` vanish for all parameters in `π` set to zero.
"""
function allvanish(F::Vector{QQMPolyRingElem}, π::Vector{QQMPolyRingElem})
  z = zeros(parent(π[1]), length(π))
  for f∈F
    if !iszero(evaluate(f, π, z))
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
- `τ` divides the characteristic polynomial `Χ(τ)` of `D₁f(x₀,π⁺)` exactly with power `s`
These properties can be used to filter out possible TFPV candidates.

We are interested in partial solutions of the system of polynomials defined by
the conditions above.
In particular, we only consider conditions on the parameters, since there might
be multiple slow manifolds and reductions.
Thus, we eliminate the dynamic variables from the ideal and obtain a set of
polynomials whose variables are the parameters of the ODE system.
This function computes a generating set for this elimination ideal and all TFPVs
lie in its vanishing set.

See also: [`tfpv_candidates`](@ref) 
"""
function determinants_criterion(problem::ReductionProblem)
  # number of dimensions to reduce the system by
  r = length(problem.x) - problem.s
  # determinants of all k×k minors of J for k>r
  d = get_determinants(problem.J, r)
  # All polynomials that generate the ideal used to determine TFPVs
  poly_gens = [problem.f; d]
  # Build ideal from polynomial expressions 
  I = ideal(parent(problem.f[1]), poly_gens)
  # eliminate dynamic variables
  Iₓ = eliminate(I, problem.x)
  # Generating set for Iₓ
  G = gens(Iₓ)
  return G
end

## Filter functions
"""
    $(TYPEDSIGNATURES)

Use the Krull dimension of the ideal generated by the unperturbed part of the
system to filter out TFPV candidates form all possible slow-fast separations. 

### Description
The slow manifold on which the reduced system is defined is contained in
`𝑉(f⁰)`, the affine variety of `f⁰`, i.e. the zero set of the fast /
unperturbed part of the system (when we consider the entries as polynomials in
the dynamic varieties `x`, so that the variety is a subset of the phase space).
Thus, `𝑉(f⁰)` needs to have an irreducible component with dimension `s`, which
we can check using the Krull dimension of the corresponding ideal.

By default, all 2ᵐ-2 possible slow-fast separations of the m parameters are
checked, but if `idx::Vector{Vector{Bool}}` is defined, only those candidates
are considered (0: slow variable, 1: fast variable).

If `compute_primary_decomposition=true` (default behaviour) is set, the
function attempts to compute a primary decomposition of the ideal corresponding
to the unperturbed part of the system `f⁰`.

If in addition, `exact_dimension=true` (default), only those candidates are
kept for which the affine variety has exactly dimension `problem.s`. Keep in
mind that this is the Krull dimension, so the topological dimension over the
reals can be smaller. However, the existence of a non-singular point in the
irreducible component implies that topological and Krull dimension are equal. 
Since we require such a point later, we may use the exact dimension already. 

See also: [`determinants_criterion`](@ref)
"""
function dimension_criterion(
  problem::ReductionProblem; 
  idx::Vector{Vector{Bool}}=[Bool[]],
  compute_primary_decomposition::Bool=true, 
  exact_dimension::Bool=true)
  # redefine RHS of ode system: Interpret parameters as coefficients and only
  # use dynamic variables as variables for the polynomial ring
  x_str = string.(problem.x)
  θ_str = string.(problem.θ)
  K, θ = rational_function_field(QQ, θ_str)
  R, x = polynomial_ring(K, x_str)
  π = θ[problem.idx_slow_fast]
  # filter TFPV candidates 
  idx = length(idx[1]) == 0 ? num2bin.(1:(2^length(problem.π)-2), length(problem.π)) : idx
  idx_candidates = zeros(Bool, length(idx))
  V = compute_primary_decomposition ? [] : nothing
  dim_components = compute_primary_decomposition ? Vector{Vector{Int}}() : nothing
  cnt = 1
  for i in idx    
    # Get unperturbed part of system (fast part)
    _θ = θ
    _θ[problem.idx_slow_fast] = π .* i 
    f⁰ = problem._f(x, _θ)
    # Check if Krull dimension is at least s
    I = ideal(R, f⁰)
    if dim(I) >= problem.s
      # compute the irreducible components of V(f⁰)
      if compute_primary_decomposition 
        PD = primary_decomposition(I)
        Y = [gens(groebner_basis(Q[2]; complete_reduction=true)) for Q in PD]
        Y_dim = dim.([ideal(R, Yᵢ) for Yᵢ in Y])
        # check if there is an irreducible component of the affine variety
        # V(f⁰) with dimension exactly s
        if exact_dimension && any(Y_dim .== problem.s)
          push!(V, Y)
          push!(dim_components, Y_dim)
          idx_candidates[cnt] = true
        end
      # do not compute the irreducible components of V(f⁰)
      else
        idx_candidates[cnt] = true
      end
    end
    cnt += 1
  end
  # return list with slow-fast separation candidates as defined by boolean indices
  _idx = idx[idx_candidates]
  return _idx, (V, dim_components)
end


"""
    $(TYPEDSIGNATURES)

Return parameter vector where all slow rates are set to zero.
"""
function get_tfpv(p, idx_slow_fast, idx)
  _p = copy(p)
  _p[idx_slow_fast] = _p[idx_slow_fast].*idx
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

Find all simple TFPVs by checking whether the affine variety `V(f⁰)` contains
an irreducible component `Y` of dimension `s` and the characteristic polynomial
of `D₁f(x,π⁺)`, where `x∈Y` and `π⁺` is a simple TFPV candidate (i.e. a
slow-fast separation of rates).
"""
function tfpv_candidates(problem)
  # check all possible slow-fast separations for sufficient conditions to be a TFPV for dimension s
  slow_fast = num2bin.(1:(2^length(problem.π)-2), length(problem.π)) 
  # define all polynomials and the Jacobian of f in ℝ(π)[x]
  F, _p = rational_function_field(QQ, string.(problem.θ))
  R, _x = polynomial_ring(F, string.(problem.x))
  f = problem._f(_x, _p)
  J = matrix(parent(f[1]), [[derivative(fᵢ, xᵢ) for xᵢ in _x] for fᵢ in f])
  # keep track of which slow-fast separation satisfies the conditions and save
  # irreducible components of V(f⁰) and their dimensions
  idx_keep = Vector{Vector{Bool}}(undef, length(slow_fast))
  components = Vector{Vector{MPolyIdeal}}()
  dim_components = Vector{Vector{Int64}}()
  for i in eachindex(slow_fast)
    sf = slow_fast[i]
    tfpv_candidate = get_tfpv(_p, problem.idx_slow_fast, sf)
    f⁰ = problem._f(_x, tfpv_candidate)
    I = ideal(f⁰)
    dim_I = dim(I)
    if dim_I < problem.s
      idx_keep[i] = [false]
    else
      PD = primary_decomposition(ideal(f⁰))
      Y = [Q[2] for Q in PD]
      dim_Y = dim.(Y)
      # set the slow parameters to zero in Jacobian
      J_π = map(f -> update_cofficients(f, tfpv_candidate), J)
      keep_i = [false for _ in Y]
      for j in eachindex(Y) 
        # check if dimension of irreducible component Yⱼ is as desired
        if dim_Y[j] == problem.s
          # substitute x = x₀ ∈ Yⱼ
          M = map(f -> normal_form(f, Y[j]), J_π)
          # Let Χ(τ) = τⁿ+ σₙ(x,π)τ⁽ⁿ⁻¹⁾ + … + σ₁(x,π)τ + σ₀(x,π) be the characteristic polynomial 
          # for x₀ ∈ Yⱼ and π⁺ a TFPV for dimension s we have σₛ(x₀,π⁺) ≠ 0 
          keep_i[j] = coeff(charpoly(M), problem.s) != 0
        end 
      end
      idx_keep[i] = keep_i
      if any(keep_i)
        push!(components, Y)
        push!(dim_components, dim_Y)
      end
    end
  end
  idx_tfpv = any.(idx_keep)
  gens_components = [gens.(groebner_basis.(v; complete_reduction=true)) for v in components]
  return slow_fast[idx_tfpv], gens_components, dim_components, idx_keep[idx_tfpv]
end

# filter based on elimination ideal
function filter_determinants(problem)
  G = determinants_criterion(problem)
  slowfast = num2bin.(1:(2^length(problem.π)-2), length(problem.π))
  idx_keep = [allvanish(G, problem.π[k .== false]) for k in slowfast]
  return slowfast[idx_keep], G
end

function tfpv_candidates_elimination_based(problem; 
                                compute_primary_decomposition::Bool=true,
                                exact_dimension::Bool=true)
  # Use determinants as first filter
  idx_det, G = filter_determinants(problem)
  # Use Krull dimension as second filter
  idx_TFPV, V = dimension_criterion(problem; idx=idx_det, compute_primary_decomposition=compute_primary_decomposition, exact_dimension=exact_dimension) 
  return idx_TFPV, G, V
end
