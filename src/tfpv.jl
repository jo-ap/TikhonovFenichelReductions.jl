
"""
Type that defines a Tikhonov-Fenichel reduction problem, i.e. a polynomial ODE
system, for which all slow-fast separations of rates yielding a reduction onto
dimension `s` are considered. 

The type `QQMPolyRingElem` is used in [Oscar.jl](https://www.oscar-system.org/)
to represent elements of a polynomial ring with coefficients in `ℚ`.

### Fields
    $(TYPEDFIELDS)
"""
struct ReductionProblem 
  "RHS of ODE system as a vector of polynomials"
  f::Vector{QQMPolyRingElem}
  "Vector of state variables"
  x::Vector{QQMPolyRingElem}
  "Vector of all parameters"
  p::Vector{QQMPolyRingElem}
  "Dimension of reduced system"
  s::Int
  "Jacobian of `f` w.r.t. `x`"
  J::MatSpaceElem{QQMPolyRingElem}
  "RHS of ODE system as a Julia function with signature `_f(x,p)`"
  _f::Function
  "ℚ(p,x): Fraction field over ℚ[p,x]"
  _F::FracField{QQMPolyRingElem}
  "ℚ(p): Field of rational functions in the parameters"
  _Fp::RationalFunctionField{QQFieldElem, QQMPolyRingElem}
  "ℚ(p)[x]: Polynomial ring in x over the field `Fp`"
  _Rx::MPolyRing{RationalFunctionFieldElem{QQFieldElem, QQMPolyRingElem}}
  "RHS of ODE system in polynomial ring Rx=ℚ(p)[x]"
  _f_Rx::Vector{MPoly{RationalFunctionFieldElem{QQFieldElem, QQMPolyRingElem}}}
  "Vector of state variables in polynomial ring Rx=ℚ(p)[x]"
  _x_Rx::Vector{MPoly{RationalFunctionFieldElem{QQFieldElem, QQMPolyRingElem}}}
end

"""
    $(TYPEDSIGNATURES)

Constructor for `ReductionProblem` Type.

### Arguments 
- `f(x,p)::Function`: Julia function defining the RHS of the ODE system 
- `x::Vector{String}`: Vector of state variables 
- `p::Vector{String}`: Vector of all parameters
- `s::Int`: Dimension of reduced system 

### Description
This function is used to set up the problem of finding Tikhonov-Fenichel
Parameter Values for dimension `s`.
The names of all variables and parameters in the system are parsed to
appropriate types in
[`Oscar.jl`](https://www.oscar-system.org/),
so that the necessary conditions for the existence of a reduction onto an
`s`-dimensional slow manifold can be evaluated.

See also: [`tfpvs_and_varieties`](@ref), [`tfpvs_groebner`](@ref)
"""
function ReductionProblem(
    f::Function,
    x::Vector{T},
    p::Vector{T},
    s::Int;
  ) where {T<:Union{String,Symbol}}
  @assert s < length(x) "the dimension of the reduced system must be smaller than that of the full system, i.e. s < n"
  R, Fp, Rx, f_R, x_R, p_R, f_Rx = parse_system(f, x, p)
  F = fraction_field(R)
  J = jacobian(f_R, x_R)
  x_Rx = gens(Rx)
  ReductionProblem(f_R, x_R, p_R, s, J, f, F, Fp, Rx, f_Rx, x_Rx)
end

"""
Type that holds information on the slow manifold as an irreducible component of
the variety `V(f0)`.

    $(TYPEDFIELDS)
"""
struct Variety
  "associated ideal in Rx=ℚ(p)[x]"
  ideal::MPolyIdeal
  "generators of associated ideal parsed to R=ℚ[p,x]"
  gens_R::Vector{QQMPolyRingElem}
  "groebner basis of associated ideal in Rx=ℚ(p)[x]"
  groebner_basis::Vector{MPoly{RationalFunctionFieldElem{QQFieldElem, QQMPolyRingElem}}}
  "groebner basis of associated ideal parsed to R=ℚ[p,x]"
  groebner_basis_R::Vector{QQMPolyRingElem}
  "Matrix `T`, such that `gens(I)*T=groebner_basis`"
  T::MatSpaceElem{MPoly{RationalFunctionFieldElem{QQFieldElem, QQMPolyRingElem}}}
  "Krull dimension of associated ideal"
  dim::Int64
end

function Variety(I::MPolyIdeal, dim::Int64, problem::ReductionProblem) 
  _I = ideal(rewrite_variety_generator.(gens(I)))
  gens_R = [parse_variety_generator(problem, rewrite_variety_generator(g)) for g in gens(_I)]
  @assert I == _I
  G, _T = groebner_basis_with_transformation_matrix(_I; complete_reduction=true)
  G = gens(G)
  # cancel_monomial_leading_coefficients!(G, _T)
  T = transpose(_T)
  gb_R = [parse_variety_generator(problem, rewrite_variety_generator(g)) for g in G]
  return Variety(_I, gens_R, G, gb_R, T, dim)
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

Parse state variables `x` and parameters `p` of polynomial ODE system so that
they can be used with Oscar.jl. 
"""
function parse_system(f::Function, x::Vector{T}, p::Vector{T}) where T<:Union{String, Symbol}
  # polynomial ring in p and x
  R,_  = polynomial_ring(QQ, [p..., x...])
  _p = gens(R)[1:length(p)]
  _x = gens(R)[length(p)+1:end]
  _f = f(_x, _p)
  # rational function field in p and polynomial ring in x 
  Fp,_ = rational_function_field(QQ, p)
  Rx,_ = polynomial_ring(Fp, x)
  _f_Rx = f(gens(Rx), gens(Fp))
  return R, Fp, Rx, _f, _x, _p, _f_Rx
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

Check if all polynomials in `F` vanish for all parameters in `p_slow` set to zero.
"""
function allvanish(F::Vector{QQMPolyRingElem}, p_slow::Vector{QQMPolyRingElem})
  z = zeros(parent(p_slow[1]), length(p_slow))
  for f∈F
    if !iszero(evaluate(f, p_slow, z))
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
Thus, we eliminate the state variables from the ideal generated by the
polynomial conditions above.
This function computes a generating set for this elimination ideal and all
TFPVs lie in its vanishing set.

See also: [`tfpvs_and_varieties`](@ref) 
"""
function tfpvs_groebner(problem::ReductionProblem)
  # number of dimensions to reduce the system by
  r = length(problem.x) - problem.s
  # determinants of all k×k minors of J for k>r
  d = get_determinants(problem.J, r)
  # All polynomials that generate the ideal used to determine TFPVs
  poly_gens = [problem.f; d]
  # Build ideal from polynomial expressions 
  I = ideal(parent(problem.f[1]), poly_gens)
  # eliminate state variables
  Iₓ = eliminate(I, problem.x)
  # Generating set for Iₓ
  G = groebner_basis(Iₓ; complete_reduction=true)
  G = gens(G)
  return G
end

"""
    $(TYPEDSIGNATURES)

Return parameter vector where all slow rates are set to zero.
"""
get_tfpv(p, sf_separation) = p.*sf_separation

function get_f0_Rx(
    problem::ReductionProblem,
    sf_separation::Vector{RationalFunctionFieldElem{QQFieldElem, QQMPolyRingElem}}
  )::Vector{MPoly{RationalFunctionFieldElem{QQFieldElem, QQMPolyRingElem}}}
  x = gens(problem._Rx)
  return problem._f(x, sf_separation)
end

function update_cofficients(f::MPoly{RationalFunctionFieldElem{QQFieldElem, QQMPolyRingElem}}, p::Vector{RationalFunctionFieldElem{QQFieldElem, QQMPolyRingElem}})
  if f == 0
    return f 
  else
    coeffs = [c for c in coefficients(f)]
    mons = [m for m in monomials(f)]
    coeffs_new = [evaluate(c, p) for c in coeffs]
    return sum(coeffs_new .* mons)
  end
end

"""
    $(TYPEDSIGNATURES)

Find all slow-fast separations `π⁺` that are TFPVs by using the necessary conditions
- the affine variety `V(f0)` contains an irreducible component `Y` of dimension `s` 
- the `s`-th coefficient of the characteristic polynomial of `D₁f(x,π⁺)` is non-zero for `x∈Y`

`preset` can be used to restrict the search for TFPVs to slow-fast separations
of rates with predefined slow/fast rates, e.g. for `preset=(α=0, β=1)` only TFPVs
with α ∈ O(ε) and β ∈ O(1) are considered.

### Description
The irreducible components are obtained by computing a minimal primary decomposition. 
The Jacobian at a point in an irreducible component `Y` is constructed
symbolically by computing normal forms with respect to a Gröbner basis `G`, s.t.
`G` generates the vanishing ideal of `Y`.
To obtain all general TFPVs and not just slow-fast separations, one can use the
function `tfpvs_groebner`.

See also: [`tfpvs_groebner`](@ref), [`print_results`](@ref), [`print_tfpvs`](@ref), [`print_varieties`](@ref)
"""
function tfpvs_and_varieties(problem::ReductionProblem; preset::NamedTuple=NamedTuple())
  # check possible slow-fast separations for sufficient conditions to be a TFPV for dimension s
  # keep track of which slow-fast separation satisfies the conditions and save
  # irreducible components of V(f0) and their dimensions
  tfpv_candidates = get_tfpv_candidates(problem, preset)
  N = length(tfpv_candidates)
  idx_keep = zeros(Bool, N)
  varieties = Vector{Vector{Variety}}(undef, N)
  for i in 1:N
    tfpv_candidate = get_tfpv(gens(problem._Fp), tfpv_candidates[i])
    # define all polynomials and the Jacobian of f in ℝ(p)[x]
    f0 = get_f0_Rx(problem, tfpv_candidate)
    J_p_slow = jacobian(f0, problem._x_Rx)
    I = ideal(f0)
    if dim(I) >=  problem.s 
      keep_i, Y, dim_Y = _get_varieties(problem, I, J_p_slow)
      if keep_i
        idx_keep[i] = true
        varieties[i] = [Variety(Y[k], dim_Y[k], problem) for k in eachindex(Y)]
      end
    end
  end
  return tfpv_candidates[idx_keep], varieties[idx_keep]
end

function get_tfpv_candidates(problem, preset)
  if isempty(preset)
    # check all slow-fast separations
    return num2bin.(1:(2^length(problem.p)-2), length(problem.p))
  else 
    # only consider slow-fast separations that match presets
    slow = zeros(Bool, length(problem.p))
    fast = zeros(Bool, length(problem.p))
    p = Symbol.(problem.p)
    for k in keys(preset)
      idx = findfirst(k .== p)
      @assert !isnothing(idx) "$(string(k)) is not a parameter"
      if preset[k] == 0
        slow[idx] = 1 
      else
        fast[idx] = 1 
      end
    end
    # free parameters
    idx_inner = .!(slow .|| fast)
    m = sum(idx_inner)
    tfpv_candidates_inner = num2bin.(0:(2^m-1), m)
    # set predefined and free parameters
    tfpv_candidates = [zeros(Bool, length(problem.p)) for _ in tfpv_candidates_inner]
    for i in eachindex(tfpv_candidates_inner)
      tfpv_candidates[i][fast] .= true 
      tfpv_candidates[i][idx_inner] = tfpv_candidates_inner[i]
    end
    return tfpv_candidates
  end
end

function _get_varieties(
    problem::ReductionProblem,
    I::MPolyIdeal{MPoly{RationalFunctionFieldElem{QQFieldElem, QQMPolyRingElem}}},
    J::MatSpaceElem{MPoly{RationalFunctionFieldElem{QQFieldElem, QQMPolyRingElem}}}
  )
  PD = primary_decomposition(I)
  Y = [Q[2] for Q in PD]
  dim_Y = dim.(Y)
  # set the slow parameters to zero in Jacobian
  keep_i = false
  for j in eachindex(Y) 
    # check if dimension of irreducible component Yⱼ is as desired
    if dim_Y[j] == problem.s 
      # substitute x = x0 ∈ Yⱼ
      # we need that Y[j] is the radical of Qᵢ in order for normal forms work !
      M = map(f -> normal_form(f, Y[j]), J)
      # Let Χ(τ) = τⁿ+ σₙ₋₁(x,p)τ⁽ⁿ⁻¹⁾ + … + σ₁(x,p)τ + σ₀(x,p) be the characteristic polynomial 
      # for x0 ∈ Yⱼ and p⁺ a TFPV for dimension s we have σₛ(x0,p⁺) ≠ 0 
      keep_i = keep_i || coeff(charpoly(M), problem.s) != 0
      if keep_i 
        break
      end
    end 
  end
  return keep_i, Y, dim_Y
end





