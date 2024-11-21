## Computing a reduction
"""
    $(TYPEDEF)

Type that holds all information to compute a Tikhonov-Fenichel reduction for a
given slow-fast separation of rates.

### Fields 
- `sf_separation::Vector{Bool}`: slow-fast separation (0: slow, 1: fast)
- `s::Int`: Dimension of reduced system (= dimension of slow manifold)
- `R::QQMPolyRing`: Ring over rationals in `x` and `p`
- `x::Vector{QQMPolyRingElem}`: Dynamic variables of system 
- `p::Vector{QQMPolyRingElem}`: All parameters
- `_p::Vector{QQMPolyRingElem}`: All parameters, where slow parameters are set to 0
- `p_sf::Vector{QQMPolyRingElem}`: Parameters, that are considered to be either small or large (all others are considered fixed)
- `idx_slow_fast::Vector{Bool}`: Boolean indices, s.t. `p_sf=p[idx_slow_fast]`
- `f::Vector{QQMPolyRingElem}`: RHS of system as vector with elements of ring `R`
- `f0::Vector{QQMPolyRingElem}`: Fast / unperturbed part of system as vector with elements of ring `R`
- `f1::Vector{QQMPolyRingElem}`: Slow / perturbed part of system as vector with elements of ring `R`
- `Df::MatSpaceElem{QQMPolyRingElem}`: Jacobian of `f`
- `Df_x0::MatSpaceElem{FracFieldElem{QQMPolyRingElem}}`: Jacobian of `f` at non-singular point `x0`
- `T::PolyRing{FracFieldElem{QQMPolyRingElem}}`: Ring in `x` over Fraction field `K`
- `chi::Poly{FracFieldElem{QQMPolyRingElem}}`: Characteristic polynomial of `Df_x0`
- `M::Vector{FracFieldElem{QQMPolyRingElem}}`: Slow manifold defined in all components of system
- `x0::Vector{FracFieldElem{QQMPolyRingElem}}`: Non-singular point in the irreducible component of V(f0) containing the slow manifold
- `K::FracField{QQMPolyRingElem}`: Fraction field in `p`
- `P::MatSpaceElem{FracFieldElem{QQMPolyRingElem}}`: Matrix with rational functions, such that `f0=P⋅Psi`
- `Psi::MatSpaceElem{FracFieldElem{QQMPolyRingElem}}`: Vector with polynomials, such that `f0=P⋅Psi`
- `DPsi::MatSpaceElem{FracFieldElem{QQMPolyRingElem}}`: Jacobian of `Psi`
- `success::Vector{Bool}`: Indicates whether slow manifold `M`, non-singular point `x0` and product decomposition `f0=P⋅Psi` have been set successfully
"""
mutable struct Reduction
  sf_separation::Vector{Bool}
  s::Int
  R::QQMPolyRing
  x::Vector{QQMPolyRingElem}
  p::Vector{QQMPolyRingElem}
  _p::Vector{QQMPolyRingElem}
  p_sf::Vector{QQMPolyRingElem}
  idx_slow_fast::Vector{Bool}
  f::Vector{QQMPolyRingElem}
  f0::Vector{QQMPolyRingElem}
  f1::Vector{QQMPolyRingElem}
  Df::MatSpaceElem{QQMPolyRingElem}
  Df_x0::MatSpaceElem{FracFieldElem{QQMPolyRingElem}}
  T::PolyRing{FracFieldElem{QQMPolyRingElem}}
  chi::Poly{FracFieldElem{QQMPolyRingElem}}
  M::Vector{FracFieldElem{QQMPolyRingElem}}
  x0::Vector{FracFieldElem{QQMPolyRingElem}}
  K::FracField{QQMPolyRingElem}
  P::MatSpaceElem{FracFieldElem{QQMPolyRingElem}}
  Psi::MatSpaceElem{FracFieldElem{QQMPolyRingElem}}
  DPsi::MatSpaceElem{FracFieldElem{QQMPolyRingElem}}
  success::Vector{Bool}
end


"""
    $(TYPEDSIGNATURES)

Print slow and fast parameters.
"""
function show_slow_fast(reduction::Reduction)
  p = reduction.p
  sf_separation = reduction.sf_separation
  println("slow: " * join(string.(p[.!sf_separation]), ", "))
  println("fast: " * join(string.(p[sf_separation]), ", "))
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
function jacobian(f::Vector{QQMPolyRingElem}, x::Vector{QQMPolyRingElem})
  matrix(parent(f[1]), [[derivative(fᵢ, xᵢ) for xᵢ in x] for fᵢ in f])
end

"""
    $(TYPEDSIGNATURES)

Compute the irreducible components of V(f⁰) and their dimensions for a given
TFPV candidate. 
If there exist a reduction, the corresponding slow manifold must be contained
in one of these components.

### Arguments
- `problem`: Reduction problem type holding information on system and dimension of reduction.
- `sf_separation`: Boolean index indicating slow-fast separation of rates (0: small, 1: large).

### Description
This function can be used if one wants to check whether a particular slow-fast
separation of rates yields a reduction for any dimension.
If the dimension of an irreducible component of V(f⁰) differs from what was
defined with `ReductionProblem`, the constructor `Reduction` can be called with
the additional argument `s` specifying the dimension.

See also: [`Reduction`](@ref)

"""
function slow_manifolds(problem::ReductionProblem, sf_separation::Vector{Bool})
  F, p = rational_function_field(QQ, string.(problem.p))
  R, x = polynomial_ring(F, string.(problem.x))
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
  K = fraction_field(R)
  _p = copy(problem.p)
  _p[problem.idx_slow_fast] = problem.p_sf.*sf_separation
  n = length(problem.x)
  r = n - s
  f0, f1 = splitsystem(problem.f, problem.p_sf, sf_separation)
  Df = jacobian(problem.f, problem.x)
  T, _ = polynomial_ring(K, "λ")
  M = K.(problem.x)
  x0 = zeros(K, n)
  P = zero_matrix(K,n,r)
  Psi = zero_matrix(K,r,1)
  DPsi = zero_matrix(K,r,n)
  Df_x0 = matrix(K, Matrix(Df))
  return Reduction(sf_separation, s, R, problem.x, problem.p, _p, problem.p_sf, problem.idx_slow_fast, problem.f, f0, f1, Df, Df_x0, T, T(0), M, x0, K, P, Psi, DPsi, zeros(Bool, 3))
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
  M = parse_ring(reduction.K, M)
  n = length(reduction.x)
  @assert length(M) == n "The slow manifold M must be defined in $n components."
  _f0 = [evaluate(fᵢ, [M; reduction.K.(reduction.p)]) for fᵢ in reduction.f0]
  f_vanishes = all(iszero.(_f0))
  if !f_vanishes 
    @warn "f0 does no vanish on the slow manifold"
  else
    reduction.M = M
    reduction.success[1] = true
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
  K = parent(v[1])
  _M = matrix(K, Matrix(M))
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
  eval_mat(reduction.Df, [reduction.M; reduction.K.(reduction._p)])
end

"""
    $(TYPEDSIGNATURES)

Return the Jacobian `D₁f(x₀,π⁺)` at the point `x₀` on the slow manifold for a
TFPV `π⁺`.

See also: [`Reduction`](@ref), [`set_point!`](@ref), [`set_manifold!`](@ref)
"""
function jacobian_tfpv_at_x0(reduction::Reduction)
  eval_mat(reduction.Df, [reduction.x0; reduction.K.(reduction._p)])
end


function _set_point!(reduction::Reduction, x0::AbstractVector)::Bool
  x0 = parse_ring(reduction.K, x0)
  n = length(reduction.x)
  @assert length(x0) == n "The point x0 must have $n components."
  # compute characteristic polynomial
  Df_x0 = eval_mat(reduction.Df, [x0; reduction.K.(reduction._p)])
  chi = charpoly(reduction.T, Df_x0)
  # check condition for coefficients
  c = collect(coefficients(chi))
  check_chi = all(iszero.(c[1:reduction.s])) && !iszero(c[reduction.s+1])
  if check_chi
    reduction.x0 = x0
    reduction.Df_x0 = Df_x0
    reduction.chi = chi
    reduction.success[2] = true
  end
  return check_chi 
end

"""
    $(TYPEDSIGNATURES)

Set non-singular point on irreducible component of V(f⁰) corresponding to the slow manifold. 
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

Set product decomposition `f0=P⋅Psi` locally satisfying `𝑉(f0)=𝑉(Psi)`, where `P`
is a matrix of rational functions  and `Psi` is a vector of polynomials.

See also: [`set_manifold!`](@ref) [`set_point!`](@ref) [`Reduction`](@ref)
"""
function set_decomposition!(reduction::Reduction, P::MatSpaceElem, Psi)
  n = length(reduction.x)
  r = n - reduction.s
  try Psi = reshape(Psi, r, 1)
  catch
    println("Psi must be of size $r or $r×1")
  end
  DPsi = jacobian(reshape(Psi, r), reduction.x)
  DPsi = parent(reduction.DPsi)(reduction.K.(DPsi))
  Psi = reduction.K.(Psi)
  Psi = parent(reduction.Psi)(Psi)
  P = reduction.K.(P)
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
end,
function set_decomposition!(reduction::Reduction, P::VecOrMat, Psi)
  n = length(reduction.x)
  r = n - reduction.s
  try P = reshape(P, n, r)
  catch
    println("P must be of size $n×$r")
  end
  P = reduction.K.(P)
  P = parent(reduction.P)(P)
  set_decomposition!(reduction, P, Psi)
end
  
# try computing matrix of rational functions P from Psi
function get_P(reduction::Reduction, Psi) 
  if size(Psi, 1) == 1
    return reduction.f0.//Psi
  else 
    U, Q, H = reduce_with_quotients_and_unit(reduction.f0, Psi)
    if any(H .!= 0)
      @warn "Could not automatically compute P"
      return nothing
    else
      return U*Q
    end
  end
end

"""
    $(TYPEDSIGNATURES)

Try to automatically compute matrix of rational functions `P` from given vector
of polynomials `Psi`, such that `f0=P⋅Psi` and `𝑉(f0)=𝑉(Psi)` holds locally.

NOTE: This works always if the drop in dimension `r=n-s=1`, but is experimental
for `r>1`

### Description
Typically, `Psi` can be chosen from `s` independent entries of `f0`. 
Practically one can consider the generators of the ideals defining the
irreducible components of `𝑉(f0)` as entries for `Psi` (possibly rewriting the
rational equations as polynomials by multiplying appropriately with
parameters occurring in a denominator).
"""
function set_decomposition!(reduction::Reduction, Psi)
  P = get_P(reduction, Psi)
  isnothing(P) ? false : set_decomposition!(reduction, P, Psi)
end


"""
    $(TYPEDSIGNATURES)

Compute the reduced system after the slow manifold, non-singular point and
product decomposition have been set successfully.

The function returns a tuple containing the reduced system in raw form and with
variables substituted according to the slow manifold.

See also: [`set_manifold!`](@ref), [`set_decomposition!`](@ref), [`set_point!`](@ref), [`compute_bulk_reductions`](@ref), [`Reduction`](@ref)
"""
function compute_reduction(reduction::Reduction)
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
    Iₙ = diagonal_matrix(reduction.K(1), n)
    f1 = matrix_space(reduction.K, n, 1)(reduction.f1)
    f_red_raw = (Iₙ - Q)*f1
    # reshape as vector
    f_red = reshape(Matrix(f_red_raw), n)
  else
    # reduction cannot be computed
    @error "Reduced system cannot be computed. You need to set a valid product decomposition, i.e. P and Psi such that f0 = P⋅Psi"
    return nothing, nothing
  end
  # Check if slow manifold is set 
  if reduction.success[1]
    # substitute components according to slow manifold
    a = reduction.K.([reduction.M; reduction.p])
    f_red_subs = [evaluate(f, a) for f in f_red]
    return f_red, f_red_subs
  else
    @warn "Slow manifold has not been defined succesfully. Reduced system is only returned in raw form, i.e. the reduced components are not substituted according to the slow manfold."
    return f_red, nothing
  end
end


function compute_directional_reduction(reduction::Reduction, γ::Function)
  # check if curve satisfies γ(0) = π⁺
  @assert all(γ(reduction.p, 0) .== reduction._p) "The curve γ must satisfy γ(0) = π⁺"
  # compute γ'(0)
  Rδ, δ = polynomial_ring(reduction.R, :δ)
  dγ = derivative.(γ(reduction.p, δ))
  dγ_0 = matrix(reduction.K ,length(reduction.p), 1, evaluate.(dγ, 0))
  # compute D₂f(x, p_sfˣ)
  D₂f = reduction.K.(eval_mat(jacobian(reduction.f, reduction.p), reduction.p, reduction._p))
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
    Iₙ = diagonal_matrix(reduction.K(1), n)
    f_red_raw = (Iₙ - Q)*D₂f*dγ_0
    # reshape as vector
    f_red = reshape(Matrix(f_red_raw), n)
  else
    # reduction cannot be computed
    @error "Reduced system cannot be computed. You need to set a valid product decomposition, i.e. P and Psi such that f0 = P⋅Psi"
    return nothing, nothing
  end
  # Check if slow manifold is set 
  if reduction.success[1]
    # substitute components according to slow manifold
    a = reduction.K.([reduction.M; reduction.p])
    f_red_subs = [evaluate(f, a) for f in f_red]
    return f_red, f_red_subs
  else
    @warn "Slow manifold has not been defined succesfully. Reduced system is only returned in raw form, i.e. the reduced components are not substituted according to the slow manfold."
    return f_red, nothing
  end
end
