## Computing a reduction
"""
Type that holds all information to compute a Tikhonov-Fenichel reduction for a
given slow-fast separation of rates.

### Fields
    $(TYPEDFIELDS)

"""
mutable struct Reduction
  "slow-fast separation (0: slow, 1: fast)"
  sf_separation::Vector{Bool}
  "Dimension of reduced system (= dimension of slow manifold)"
  s::Int
  "Polynomial ring `ℚ[p,x]`"
  R::QQMPolyRing
  "Polynomial ring `ℚ(p)[x]`"
  Rx::MPolyRing{RationalFunctionFieldElem{QQFieldElem, QQMPolyRingElem}}
  "Rational function field `ℚ(p,x)`"
  F::FracField{QQMPolyRingElem}
  "Dynamic variables of system "
  x::Vector{QQMPolyRingElem}
  "Parameters of the system"
  p::Vector{QQMPolyRingElem}
  "Parameters of the system where small ones are set to 0"
  _p::Vector{QQMPolyRingElem}
  "Parameters that can be small (all others are considered fixed)"
  p_sf::Vector{QQMPolyRingElem}
  "Boolean indices, s.t. `p_sf=p[idx_slow_fast]`"
  idx_slow_fast::Vector{Bool}
  "RHS of system as vector with elements of ring `R`"
  f::Vector{QQMPolyRingElem}
  "RHS of system as julia function with signature `_f(x,p)`"
  _f::Function
  "Fast part of the system"
  f0::Vector{QQMPolyRingElem}
  "Slow part of the system"
  f1::Vector{QQMPolyRingElem}
  "Jacobian of `f0`"
  Df0::MatSpaceElem{QQMPolyRingElem}
  "Jacobian of `f0` at the non-singular point `x0`"
  Df0_at_x0::MatSpaceElem{FracFieldElem{QQMPolyRingElem}}
  "Polynomial ring `Q(p,x)[λ]` for characteristic polynomial of `Df0`"
  T::PolyRing{FracFieldElem{QQMPolyRingElem}}
  "Characteristic polynomial of `Df0`"
  chi::Poly{FracFieldElem{QQMPolyRingElem}}
  "Components of the system on slow manifold"
  M::Vector{FracFieldElem{QQMPolyRingElem}}
  "Non-singular point in the irreducible component of `V(f0)` containing the slow manifold"
  x0::Vector{FracFieldElem{QQMPolyRingElem}}
  "Matrix with rational functions, such that `f0=P⋅Psi`"
  P::MatSpaceElem{FracFieldElem{QQMPolyRingElem}}
  "Vector with polynomials, such that `f0=P⋅Psi`"
  Psi::MatSpaceElem{FracFieldElem{QQMPolyRingElem}}
  "Jacobian of `Psi`"
  DPsi::MatSpaceElem{FracFieldElem{QQMPolyRingElem}}
  "Indicates whether slow manifold `M`, non-singular point `x0` and product decomposition `f0=P⋅Psi` have been set successfully to allow the computation of the reduced system"
  success::Vector{Bool}
  "Boolean indices of components that determine the flow on the slow manifold"
  idx_components::Vector{Bool}
  "Reduced system in general form (before substituting variables according to the slow manifold)"
  g_raw::Vector{FracFieldElem{QQMPolyRingElem}}
  "Reduced system on the slow manifold (`s`-dimensional)"
  g::Vector{FracFieldElem{QQMPolyRingElem}}
  "Whether `g` and `g_raw` are already computed"
  reduction_cached::Vector{Bool}
end

# conversion functions
function R_to_Rx(f::QQMPolyRingElem, reduction::Reduction)
  _p = gens(base_ring(reduction.Rx))
  _x = gens(reduction.Rx)
  f_Rx = zero(reduction.Rx)
  for (c,a) in coefficients_and_exponents(f)
    f_Rx += c*prod([_p; _x].^a)
  end
  return f_Rx
end
function Rx_to_F(f::MPoly{RationalFunctionFieldElem{QQFieldElem, QQMPolyRingElem}}, reduction::Reduction)
  _p = gens(reduction.F)[1:length(reduction.p)]
  _x = gens(reduction.F)[length(reduction.p)+1:end]
  f_F = zero(reduction.F)
  for (c,a) in coefficients_and_exponents(f)
    p = evaluate(numerator(c), _p)
    q = evaluate(denominator(c), _p)
    f_F += p*prod(_x.^a)//q
  end
  return f_F
end

"""
    $(TYPEDSIGNATURES)

Split RHS of system into fast/unperturbed and slow/perturbed part for a given
slow-fast separation of rates.
"""
function splitsystem(f::Vector{QQMPolyRingElem}, p_sf::Vector{QQMPolyRingElem}, sf_separation::Vector{Bool})
  R = parent(f[1])
  f0 = [evaluate(fᵢ, p_sf[.!sf_separation], zero.(p_sf[.!sf_separation])) for fᵢ in f]
  f1 = f .- f0 
  return f0, f1
end

"""
    $(TYPEDSIGNATURES)

Compute Jacobian of `f` with respect to `x`.
"""
function jacobian(f, x)
  matrix(parent(f[1]), [[derivative(fᵢ, xᵢ) for xᵢ in x] for fᵢ in f])
end

"""
    $(TYPEDSIGNATURES)

Compute the irreducible components of `V(f0)` and their dimensions for a given
TFPV candidate. 
If there exist a reduction, the corresponding slow manifold must be contained
in one of these components.

### Arguments
- `problem`: Reduction problem type holding information on system and dimension of reduction.
- `sf_separation`: Boolean index indicating slow-fast separation of rates (0: small, 1: large).

### Description
This function can be used if one wants to check whether a particular slow-fast
separation of rates yields a reduction for any dimension.
If the dimension of an irreducible component of `V(f0)` differs from what was
defined with `ReductionProblem`, the constructor `Reduction` can be called with
the additional argument `s` specifying the dimension.

See also: [`Reduction`](@ref)

"""
function slow_manifolds(problem::ReductionProblem, sf_separation::Vector{Bool})
  F, p = rational_function_field(QQ, string.(problem.p))
  _, x = polynomial_ring(F, string.(problem.x))
  p_sf = p
  p_sf[sf_separation .== 0] .= F(0)
  f = problem._f(x,p_sf)
  PD = primary_decomposition(ideal(f))
  Q = [q[2] for q in PD]
  dim_Q = dim.(Q)
  return Q, dim_Q
end

"""
    $(TYPEDSIGNATURES)
Constructor for `Reduction` Type.

### Arguments
- `problem`: Reduction problem type holding information on system and dimension of reduction.
- `sf_separation`: Boolean index indicating slow-fast separation of rates (0: small, 1: large).
- `s::Int`: (optional) Dimension of slow manifold. Can be specified if a reduction corresponding to a TFPV for dimension different from `problem.s` should be considered (e.g. for manually computing a reduction for a given slow-fast separation that is not necessarily obtained via `tfpv_candidates`).

See also: [`set_manifold!`](@ref) [`set_decomposition!`](@ref)
"""
function Reduction(problem::ReductionProblem, sf_separation::Vector{Bool}; s::Union{Nothing,Int}=nothing)
  s = isnothing(s) ? problem.s : s
  R = parent(problem.f[1])
  F = fraction_field(R)
  Fp, _ = rational_function_field(QQ, string.(problem.p))
  Rx, _ = polynomial_ring(Fp, string.(problem.x))
  _p = copy(problem.p)
  _p[problem.idx_slow_fast] = problem.p_sf.*sf_separation
  n = length(problem.x)
  r = n - s
  f0, f1 = splitsystem(problem.f, problem.p_sf, sf_separation)
  Df0 = jacobian(problem.f, problem.x)
  T, _ = polynomial_ring(F, "λ")
  M = F.(problem.x)
  x0 = zeros(F, n)
  P = zero_matrix(F,n,r)
  Psi = zero_matrix(F,r,1)
  DPsi = zero_matrix(F,r,n)
  Df0_at_x0 = matrix(F, Matrix(Df0))
  return Reduction(sf_separation,
                   s,
                   R,
                   Rx, 
                   F,
                   problem.x,
                   problem.p,
                   _p,
                   problem.p_sf,
                   problem.idx_slow_fast,
                   problem.f,
                   problem._f,
                   f0,
                   f1,
                   Df0,
                   Df0_at_x0,
                   T,
                   T(0),
                   M,
                   x0,
                   P,
                   Psi,
                   DPsi,
                   zeros(Bool, 3),
                   zeros(Bool, n),
                   zeros(F, n),
                   zeros(F, s),
                   zeros(Bool,2)
                   )
end

"""
    $(TYPEDSIGNATURES)
Convenience function that constructs an object of type `Reduction` and calls
`set_manifold!` and `set_decomposition!`.

### Arguments
- `problem`: Reduction problem type holding information on system and dimension of reduction.
- `sf_separation`: Boolean index indicating slow-fast separation of rates (0: small, 1: large).
- `V`: Generators of affine variety corresponding to the the slow manifold 
- `M`: Slow manifold in explicit form
- `s::Int`: (optional) Dimension of slow manifold. Can be specified if a reduction corresponding to a TFPV for dimension different from `problem.s` should be considered (e.g. for manually computing a reduction for a given slow-fast separation that is not necessarily obtained via `tfpv_candidates`).

See also: [`set_manifold!`](@ref) [`set_decomposition!`](@ref) [`tfpv_candidates`](@ref)
"""
function Reduction(
  problem::ReductionProblem,
  sf_separation::Vector{Bool},
  V::Vector{QQMPolyRingElem},
  M::Vector{<:RingElem}; 
  s::Union{Nothing,Int}=nothing)
  reduction = Reduction(problem, sf_separation; s=s)
  sm = set_manifold!(reduction, M)
  sd = set_decomposition!(reduction, V)
  if sm & sd 
    compute_reduction!(reduction)
    return reduction
  end
end 

function parse_ring(R, x)
  try x = R.(x)
  catch
    println("Cannot parse $x into $R")
  end
  return x
end

"""
    $(TYPEDSIGNATURES)

Set the slow manifold by defining the values of the components of the system.
Note that `M` must be defined as a vector with the same length as the system's
components, i.e. `reduction.x`.

See also: [`Reduction`](@ref), [`set_decomposition!`](@ref), [`set_point!`](@ref)
"""
function set_manifold!(reduction::Reduction, M::AbstractVector)::Bool
  M = parse_ring(reduction.F, M)
  n = length(reduction.x)
  @assert length(M) == n "The slow manifold M must be defined in $n components."
  _f0 = [evaluate(fᵢ, [reduction.F.(reduction.p); M]) for fᵢ in reduction.f0]
  f_vanishes = all(iszero.(_f0))
  if !f_vanishes 
    @warn "f0 does no vanish on the slow manifold"
  else
    reduction.M = M
    reduction.success[1] = true
    reduction.idx_components = reduction.M .== reduction.x
    # check if there exists a reduction
    # an eigenvalue λ has to factor the characteristic polynomial of the
    # jacobian at a non-singular point exactly with power s, i.e. there is no
    # reduction if the polynomial is given by λ^(s+1)•r(λ)
    coeffs = collect(coefficients(charpoly(jacobian_tfpv_on_manifold(reduction))))
    if all(coeffs[1:reduction.s + 1] .== 0)
      @info "There exists no reduction onto this manifold"
      # check if a generic point on the slow manifold is non-singular
    elseif !_set_point!(reduction, M)
      @warn "Could not set generic non-singular point on slow manifold"
    end
  end
  return f_vanishes
end

function eval_mat(M, x, v)
  m,n = size(M)
  _M = copy(M)
  for i = 1:m
    for j = 1:n
      _M[i,j] = evaluate(M[i,j], x, v)
    end
  end
  return _M
end
function eval_mat(M, v)
  m,n = size(M)
  F = parent(v[1])
  _M = matrix(F, Matrix(M))
  for i = 1:m
    for j = 1:n
      _M[i,j] = evaluate(M[i,j], v)
    end
  end
  return _M
end


"""
    $(TYPEDSIGNATURES)

Return the Jacobian `D₁f(x,π⁺)` at a generic point `x` on the slow manifold for
a TFPV `π⁺`.

See also: [`Reduction`](@ref), [`set_manifold!`](@ref)
"""
function jacobian_tfpv_on_manifold(reduction::Reduction)
  eval_mat(reduction.Df0, [reduction.F.(reduction._p); reduction.M])
end

"""
    $(TYPEDSIGNATURES)

Return the Jacobian `D₁f(x₀,π⁺)` at the point `x₀` on the slow manifold for a
TFPV `π⁺`.

See also: [`Reduction`](@ref), [`set_point!`](@ref), [`set_manifold!`](@ref)
"""
function jacobian_tfpv_at_x0(reduction::Reduction)
  eval_mat(reduction.Df0, [reduction.F.(reduction._p); reduction.x0])
end


function _set_point!(reduction::Reduction, x0::AbstractVector)::Bool
  x0 = parse_ring(reduction.F, x0)
  n = length(reduction.x)
  @assert length(x0) == n "The point x0 must have $n components."
  # compute characteristic polynomial
  Df0_at_x0 = eval_mat(reduction.Df0, [reduction.F.(reduction._p); x0])
  chi = charpoly(reduction.T, Df0_at_x0)
  # check condition for coefficients
  c = collect(coefficients(chi))
  check_chi = all(iszero.(c[1:reduction.s])) && !iszero(c[reduction.s+1])
  if check_chi
    reduction.x0 = x0
    reduction.Df0_at_x0 = Df0_at_x0
    reduction.chi = chi
    reduction.success[2] = true
  end
  return check_chi 
end

"""
    $(TYPEDSIGNATURES)

Set non-singular point on irreducible component of `V(f0)` corresponding to the slow manifold. 
Typically, this can be done automatically by setting the slow manifold.

See also: [`set_manifold!`](@ref), [`set_decomposition!`](@ref), [`Reduction`](@ref)
"""
function set_point!(reduction::Reduction, x0::AbstractVector)::Bool
  retval = _set_point!(reduction, x0)
  if !retval
    @warn "The eigenvalue λ does not factor the characteristic polynomial of D₁f(x0,p_sf) with power s=$(reduction.s)"
  end
  return retval 
end
  
"""
    $(TYPEDSIGNATURES)

Set product decomposition `f0=P⋅Psi` locally satisfying `V(f0)=V(Psi)`, where `P`
is a matrix of rational functions  and `Psi` is a vector of polynomials.

See also: [`set_manifold!`](@ref) [`set_point!`](@ref) [`Reduction`](@ref)
"""
function set_decomposition!(reduction::Reduction, P::Union{MatSpaceElem,VecOrMat}, Psi)
  _set_decomposition!(reduction, P, Psi)
end

function _set_decomposition!(reduction::Reduction, P::MatSpaceElem, Psi)
  n = length(reduction.x)
  r = n - reduction.s
  try Psi = reshape(Psi, r, 1)
  catch
    println("Psi must be of size $r or $r×1")
  end
  DPsi = jacobian(reshape(Psi, r), reduction.x)
  DPsi = parent(reduction.DPsi)(reduction.F.(DPsi))
  Psi = reduction.F.(Psi)
  Psi = parent(reduction.Psi)(Psi)
  P = reduction.F.(P)
  # check if product decomposition is correct
  # parse f0 as matrix, so that they can be compared
  M = P*Psi
  is_equal = all(iszero.(M .- parent(M)(reduction.f0)))
  if is_equal
    reduction.P = P
    reduction.Psi = Psi
    reduction.DPsi = DPsi
    reduction.success[3] = true
  else
    @warn "P⋅Psi ≠ f0: Product decomposition cannot be set"
  end
  return is_equal
end
function _set_decomposition!(reduction::Reduction, P::VecOrMat, Psi)
  n = length(reduction.x)
  r = n - reduction.s
  try P = reshape(P, n, r)
  catch
    println("P must be of size $n×$r")
  end
  P = reduction.F.(P)
  P = parent(reduction.P)(P)
  _set_decomposition!(reduction, P, Psi)
end
  
# try computing matrix of rational functions P from Psi
function get_decomposition(reduction::Reduction, Psi::Vector{QQMPolyRingElem})
  if size(Psi, 1) == 1
    return reduction.f0.//Psi
  else 
    p = gens(base_ring(reduction.Rx))
    x = gens(reduction.Rx)
    _f0 = reduction._f(x, p .* reduction.sf_separation)
    _Psi = [R_to_Rx(p, reduction) for p in Psi]
    U, Q, H = reduce_with_quotients_and_unit(_f0, _Psi)
    if all(H .== 0)
      P = matrix(reduction.F, [Rx_to_F(f, reduction) for f in U*Q])
      return P, Psi
    end
    _Psi = find_independent_polys(reduction.f0)
    U, Q, H = reduce_with_quotients_and_unit(reduction.f0, _Psi)
    if all(H .== 0)
      @info "Automatically updated Psi"
      return U*Q, _Psi
    end
    @warn "Could not set P automatically."
    return nothing, Psi
  end
end

function find_independent_polys(f::Vector)
  idx = [false for _ in f]
  for i in eachindex(f)
    if is_algebraically_independent([f[idx]..., f[i]])
      idx[i] = true
    end 
  end 
  return f[idx]
end

"""
    $(TYPEDSIGNATURES)

Try to automatically compute matrix of rational functions `P` from given vector
of polynomials `Psi`, such that `f0=P⋅Psi` and `V(f0)=V(Psi)` holds locally.

NOTE: This always works if the drop in dimension `r=n-s=1`, but is experimental
for `r>1`

### Description
Typically, `Psi` can be chosen from `s` independent entries of `f0`. 
Practically one can consider the generators of the ideals defining the
irreducible components of `V(f0)` as entries for `Psi` (possibly rewriting the
rational equations as polynomials by multiplying appropriately with
parameters occurring in a denominator).
"""
function set_decomposition!(reduction::Reduction, Psi)
  P, Psi = get_decomposition(reduction, Psi)
  isnothing(P) ? false : set_decomposition!(reduction, P, Psi)
end

"""
    $(TYPEDSIGNATURES)

Compute the reduced system after the slow manifold, non-singular point and
product decomposition have been set successfully.

The function returns true if the reduced system was computed successfully.
The reduction in raw form, i.e. before substituting the variables `x` according
to the slow manifold is set `reduction.g_raw` while the `s`-dimensional
reduction on the slow manifold is given by `reduction.g`.
A safe getter function for this is `get_reduced_system(reduction::Reduction)`.

See also: [`set_manifold!`](@ref), [`set_decomposition!`](@ref), [`set_point!`](@ref), [`compute_bulk_reductions`](@ref), [`Reduction`](@ref), [`print_reduced_system`](@ref)
"""
function compute_reduction!(reduction::Reduction)
  # Check if P-Psi-composition is defined 
  if reduction.success[3]
    # check if non-singular point is defined
    if !reduction.success[2]
      @info "The non-singular point on the slow manifold has not been set successfully. Trying to compute the reduction anyway."
    end
    # dimensions
    n = length(reduction.x)
    # compute reduced system
    A = reduction.DPsi*reduction.P
    Q = reduction.P*inv(A)*reduction.DPsi
    Iₙ = diagonal_matrix(reduction.F(1), n)
    f1 = matrix_space(reduction.F, n, 1)(reduction.f1)
    f_red_raw = (Iₙ - Q)*f1
    # reshape as vector
    f_red = reshape(Matrix(f_red_raw), n)
    reduction.g_raw = f_red
    reduction.reduction_cached[1] = true
  else
    # reduction cannot be computed
    @error "Reduced system cannot be computed. You need to set a valid product decomposition, i.e. P and Psi such that f0 = P⋅Psi"
    return false
  end
  # Check if slow manifold is set 
  if reduction.success[1]
    # substitute components according to slow manifold
    a = reduction.F.([reduction.p; reduction.M])
    f_red_subs = [evaluate(f, a) for f in f_red]
    reduction.g = f_red_subs[reduction.idx_components]
    reduction.reduction_cached[2] = true
  else
    @warn "Slow manifold has not been defined succesfully. Reduced system is only returned in raw form, i.e. the reduced components are not substituted according to the slow manfold."
  end
  return true
end

# function get_reduced_system(reduction::Reduction; on_manifold=true)
#   # compute the reduced system, if this was not done before 
#   warn = true
#   if !any(reduction.reduction_cached)
#     ret = compute_reduction!(reduction)
#     warn = false
#     if !ret 
#       return 
#     end
#   end
#   # return reduced system
#   if on_manifold 
#     if reduction.reduction_cached[2] 
#       return reduction.g
#     elseif warn
#       @warn "Slow manifold has not been defined succesfully. Reduced system is only returned in raw form, i.e. the reduced components are not substituted according to the slow manfold."
#     end 
#     return 
#   end
#   return reduction.g_raw
# end


# function compute_directional_reduction(reduction::Reduction, γ::Function)
#   # check if curve satisfies γ(0) = π⁺
#   @assert all(γ(reduction.p, 0) .== reduction._p) "The curve γ must satisfy γ(0) = π⁺"
#   # compute γ'(0)
#   Rδ, δ = polynomial_ring(reduction.R, :δ)
#   dγ = derivative.(γ(reduction.p, δ))
#   dγ_0 = matrix(reduction.F ,length(reduction.p), 1, evaluate.(dγ, 0))
#   # compute D₂f(x, p_sfˣ)
#   D₂f = reduction.F.(eval_mat(jacobian(reduction.f, reduction.p), reduction.p, reduction._p))
#   # Check if P-Psi-composition is defined 
#   if reduction.success[3]
#     # check if non-singular point is defined
#     if !reduction.success[2]
#       @info "The non-singular point on the slow manifold has not been set successfully. Trying to compute the reduction anyway."
#     end
#     # dimensions
#     n = length(reduction.x)
#     # compute reduced system
#     A = reduction.DPsi*reduction.P
#     Q = reduction.P*inv(A)*reduction.DPsi
#     Iₙ = diagonal_matrix(reduction.F(1), n)
#     f_red_raw = (Iₙ - Q)*D₂f*dγ_0
#     # reshape as vector
#     f_red = reshape(Matrix(f_red_raw), n)
#   else
#     # reduction cannot be computed
#     @error "Reduced system cannot be computed. You need to set a valid product decomposition, i.e. P and Psi such that f0 = P⋅Psi"
#     return, nothing
#   end
#   # Check if slow manifold is set 
#   if reduction.success[1]
#     # substitute components according to slow manifold
#     a = reduction.F.([reduction.p; reduction.M])
#     f_red_subs = [evaluate(f, a) for f in f_red]
#     return f_red, f_red_subs
#   else
#     @warn "Slow manifold has not been defined succesfully. Reduced system is only returned in raw form, i.e. the reduced components are not substituted according to the slow manfold."
#     return f_red, nothing
#   end
# end

function _get_slow_fast(reduction::Reduction)
  p = reduction.p_sf
  sf_separation = reduction.sf_separation
  slow = p[.!sf_separation]
  fast = p[sf_separation]
  return slow, fast
end
