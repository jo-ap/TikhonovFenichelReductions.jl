## Computing a reduction
"""
    $(TYPEDEF)

Type that holds all information to compute a Tikhonov-Fenichel reduction for a
given slow-fast separation of rates.

### Fields 
- `idx::Vector{Bool}`: slow-fast separation (0: slow, 1: fast)
- `s::Int`: Dimension of reduced system (= dimension of slow manifold)
- `R::QQMPolyRing`: Ring over rationals in `x` and `θ`
- `x::Vector{QQMPolyRingElem}`: Dynamic variables of system 
- `θ::Vector{QQMPolyRingElem}`: All parameters
- `_θ::Vector{QQMPolyRingElem}`: All parameters, where slow parameters are set to 0
- `π::Vector{QQMPolyRingElem}`: Parameters, that are considered to be either small or large
- `idx_slow_fast::Vector{Bool}`: Boolean indices, s.t. `π=θ[idx_slow_fast]`
- `f::Vector{QQMPolyRingElem}`: RHS of system as vector with elements of ring `R`
- `f⁰::Vector{QQMPolyRingElem}`: Fast / unperturbed part of system as vector with elements of ring `R`
- `f¹::Vector{QQMPolyRingElem}`: Slow / perturbed part of system as vector with elements of ring `R`
- `Df::AbstractAlgebra.Generic.MatSpaceElem{QQMPolyRingElem}`: Jacobian of `f`
- `Df_x₀::AbstractAlgebra.Generic.MatSpaceElem{AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}}`: Jacobian of `f` at non-singular point `x₀`
- `T::AbstractAlgebra.Generic.PolyRing{AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}}`: Ring in `x` over Fraction field `K`
- `chi::AbstractAlgebra.Generic.Poly{AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}}`: Characteristic polynomial of `Df_x₀`
- `M::Vector{AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}}`: Slow manifold defined in all components of system
- `x₀::Vector{AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}}`: Non-singular point in the irreducible component of V(f⁰) containing the slow manifold
- `K::AbstractAlgebra.Generic.FracField{QQMPolyRingElem}`: Fraction field in `θ`
- `P::AbstractAlgebra.Generic.MatSpaceElem{AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}}`: Matrix with rational functions, such that `f⁰=P⋅ψ`
- `ψ::AbstractAlgebra.Generic.MatSpaceElem{AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}}`: Vector with polynomials, such that `f⁰=P⋅ψ`
- `Dψ::AbstractAlgebra.Generic.MatSpaceElem{AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}}`: Jacobian of `ψ`
- `success::Vector{Bool}`: Indicates whether slow manifold `M`, non-singular point `x₀` and product decomposition `f⁰=P⋅ψ` have been set successfully
"""
mutable struct Reduction
  idx::Vector{Bool}
  s::Int
  R::QQMPolyRing
  x::Vector{QQMPolyRingElem}
  θ::Vector{QQMPolyRingElem}
  _θ::Vector{QQMPolyRingElem}
  π::Vector{QQMPolyRingElem}
  idx_slow_fast::Vector{Bool}
  f::Vector{QQMPolyRingElem}
  f⁰::Vector{QQMPolyRingElem}
  f¹::Vector{QQMPolyRingElem}
  Df::AbstractAlgebra.Generic.MatSpaceElem{QQMPolyRingElem}
  Df_x₀::AbstractAlgebra.Generic.MatSpaceElem{AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}}
  T::AbstractAlgebra.Generic.PolyRing{AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}}
  chi::AbstractAlgebra.Generic.Poly{AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}}
  M::Vector{AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}}
  x₀::Vector{AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}}
  K::AbstractAlgebra.Generic.FracField{QQMPolyRingElem}
  P::AbstractAlgebra.Generic.MatSpaceElem{AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}}
  ψ::AbstractAlgebra.Generic.MatSpaceElem{AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}}
  Dψ::AbstractAlgebra.Generic.MatSpaceElem{AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}}
  success::Vector{Bool}
end

"""
    $(TYPEDSIGNATURES)

Split RHS of system into fast/unperturbed and slow/perturbed part for a given
slow-fast separation of rates.
"""
function splitsystem(f::Vector{QQMPolyRingElem}, π::Vector{QQMPolyRingElem}, idx::Vector{Bool})
  R = parent(f[1])
  f⁰ = [evaluate(fᵢ, π[.!idx], zero.(π[.!idx])) for fᵢ in f]
  f¹ = f .- f⁰ 
  return f⁰, f¹
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
- `idx`: Boolean index indicating slow-fast separation of rates (0: small, 1: large).

### Description
This function can be used if one wants to check whether a particular slow-fast
separation of rates yields a reduction for any dimension.
If the dimension of an irreducible component of V(f⁰) differs from what was
defined with `ReductionProblem`, the constructor `Reduction` can be called with
the additional argument `s` specifying the dimension.

See also: [Reduction](@ref)

"""
function slow_manifold_candidates(problem::ReductionProblem, idx::Vector{Bool})
  F, p = rational_function_field(QQ, string.(problem.θ))
  R, x = polynomial_ring(F, string.(problem.x))
  π = p
  π[idx .== 0] .= F(0)
  f = problem._f(x,π)
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
- `idx`: Boolean index indicating slow-fast separation of rates (0: small, 1: large).
- `s::Int`: (optional) Dimension of slow manifold. Can be specified if a reduction corresponding to a TFPV for dimension different from `problem.s` should be considered (e.g. for manually computing a reduction for a given slow-fast separation that is not necessarily obtained via `tfpv_candidates`).

See also: [set_manifold!](@ref), [set_decomposition!](@ref)
"""
function Reduction(problem::ReductionProblem, idx::Vector{Bool}; s::Union{Nothing,Int}=nothing)
  s = isnothing(s) ? problem.s : s
  R = parent(problem.f[1])
  K = fraction_field(R)
  _θ = copy(problem.θ)
  _θ[problem.idx_slow_fast] = problem.π.*idx
  n = length(problem.x)
  r = n - s
  f⁰, f¹ = splitsystem(problem.f, problem.π, idx)
  Df = jacobian(problem.f, problem.x)
  T, _ = polynomial_ring(K, "λ")
  M = K.(problem.x)
  x₀ = zeros(K, n)
  P = zero_matrix(K,n,r)
  ψ = zero_matrix(K,r,1)
  Dψ = zero_matrix(K,r,n)
  Df_x₀ = matrix(K, Matrix(Df))
  return Reduction(idx, s, R, problem.x, problem.θ, _θ, problem.π, problem.idx_slow_fast, problem.f, f⁰, f¹, Df, Df_x₀, T, T(0), M, x₀, K, P, ψ, Dψ, zeros(Bool, 3))
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

See also: [Reduction](@ref), [set_decomposition!](@ref), [set_point!](@ref)
"""
function set_manifold!(reduction::Reduction, M::AbstractVector)::Bool
  M = parse_ring(reduction.K, M)
  n = length(reduction.x)
  @assert length(M) == n "The slow manifold M must be defined in $n components."
  _f⁰ = [evaluate(fᵢ, [M; reduction.K.(reduction.θ)]) for fᵢ in reduction.f⁰]
  f_vanishes = all(iszero.(_f⁰))
  if !f_vanishes 
    @warn "f⁰ does no vanish on the slow manifold"
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
    elseif _set_point!(reduction, M)
      println("Set generic non-singular point on slow manifold")
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

function jacobian_tfpv_on_manifold(reduction::Reduction)
  eval_mat(reduction.Df, [reduction.M; reduction.K.(reduction._θ)])
end

function jacobian_tfpv(reduction::Reduction)
  eval_mat(reduction.Df, reduction.K.([reduction.x; reduction._θ]))
end

function _set_point!(reduction::Reduction, x₀::AbstractVector)::Bool
  x₀ = parse_ring(reduction.K, x₀)
  n = length(reduction.x)
  @assert length(x₀) == n "The point x₀ must have $n components."
  # compute characteristic polynomial
  Df_x₀ = eval_mat(reduction.Df, [x₀; reduction.K.(reduction._θ)])
  chi = charpoly(reduction.T, Df_x₀)
  # check condition for coefficients
  c = collect(coefficients(chi))
  check_chi = all(iszero.(c[1:reduction.s])) && !iszero(c[reduction.s+1])
  if check_chi
    reduction.x₀ = x₀
    reduction.Df_x₀ = Df_x₀
    reduction.chi = chi
    reduction.success[2] = true
  end
  return check_chi 
end

"""
    $(TYPEDSIGNATURES)

Set non-singular point on irreducible component of V(f⁰) corresponding to the slow manifold. 
Typically, this can be done automatically by setting the slow manifold.

See also: [set_manifold!](@ref), [set_decomposition!](@ref), [Reduction](@ref)
"""
function set_point!(reduction::Reduction, x₀::AbstractVector)::Bool
  retval = _set_point!(reduction, x₀)
  if !retval
    @warn "The eigenvalue λ does not factor the characteristic polynomial of D₁f(x₀,π) with power s=$(reduction.s)"
  end
  return retval 
end
  
"""
    $(TYPEDSIGNATURES)

Set product decomposition `f⁰=P⋅ψ` locally satisfying `𝑉(f⁰)=𝑉(ψ)`, where `P`
is a matrix of rational functions  and `ψ` is a vector of polynomials.

### Variants:
- `set_decomposition!(reduction::Reduction, P::AbstractAlgebra.Generic.MatSpaceElem, ψ)`: Manually specify `P` and `ψ`
- `set_decomposition!(reduction::Reduction,  ψ)`: Try to compute `P` automatically. This works always for `s=1`, but may fail if `s>1`.

### Description
Typically, `ψ` can be chosen from `s` independent entries of `f⁰`. 
Practically one can consider the generators of the ideals defining the
irreducible components of `𝑉(f⁰)` as entries for `ψ` (possibly rewriting the
rational equations as polynomials by multiplying appropriately with
parameters occurring in a denominator).

See also: [set_manifold!](@ref), [set_point!](@ref), [Reduction](@ref)
"""
function set_decomposition!(reduction::Reduction, P::AbstractAlgebra.Generic.MatSpaceElem, ψ)
  n = length(reduction.x)
  r = n - reduction.s
  try ψ = reshape(ψ, r, 1)
  catch
    println("ψ must be of size $r or $r×1")
  end
  Dψ = jacobian(reshape(ψ, r), reduction.x)
  Dψ = parent(reduction.Dψ)(reduction.K.(Dψ))
  ψ = reduction.K.(ψ)
  ψ = parent(reduction.ψ)(ψ)
  P = reduction.K.(P)
  # check if product decomposition is correct
  # parse f⁰ as matrix, so that they can be compared
  M = P*ψ
  is_equal = all(iszero.(M .- parent(M)(reduction.f⁰)))
  if is_equal
    reduction.P = P
    reduction.ψ = ψ
    reduction.Dψ = Dψ
    reduction.success[3] = true
  else
    @warn "P⋅ψ ≠ f⁰: Product decomposition cannot be set"
  end
  return is_equal
end
function set_decomposition!(reduction::Reduction, P::VecOrMat, ψ)
  n = length(reduction.x)
  r = n - reduction.s
  try P = reshape(P, n, r)
  catch
    println("P must be of size $n×$r")
  end
  P = reduction.K.(P)
  P = parent(reduction.P)(P)
  set_decomposition!(reduction, P, ψ)
end
# try computing matrix of rational functions P from ψ
function get_P(reduction::Reduction, ψ::VecOrMat) 
  if size(ψ, 1) == 1
    P = reduction.f⁰.//ψ
  else 
    U, Q, H = reduce_with_quotients_and_unit(reduction.f⁰, ψ)
    P = U*Q
    if any(H .!= 0)
      @warn "Could not automatically compute P"
    end
  end
  return P
end

function set_decomposition!(reduction::Reduction, ψ::VecOrMat)
  P = get_P(reduction, ψ)
  set_decomposition!(reduction, P, ψ)
end

# Experimental: Try guessing P and ψ automatically
# function get_decomposition(reduction)
#   G = groebner_basis(ideal(reduction.f⁰); complete_reduction=true)
#   ψ = gens(G)
#   while true 
#     U, Q, H = reduce_with_quotients_and_unit(reduction.f⁰, ψ)
#     idx = [any(Q[:,i] .!= 0) for i in 1:size(Q, 2)]
#     if all(idx) 
#       return U*Q, ψ
#     else
#       ψ = ψ[idx] 
#     end
#   end
# end


"""
    $(TYPEDSIGNATURES)

Compute the reduced system after the slow manifold, non-singular point and
product decomposition have been set successfully.

The function returns a tuple containing the reduced system in raw form and with
variables substituted according to the slow manifold.

See also: [set_manifold!](@ref), [set_decomposition!](@ref), [set_point!](@ref), [Reduction](@ref)
"""
function compute_reduction(reduction::Reduction)
  # Check if P-ψ-composition is defined 
  if reduction.success[3]
    # check if non-singular point is defined
    if !reduction.success[2]
      @info "The non-singular point on the slow manifold has not been set successfully. Trying to compute the reduction anyway."
    end
    # dimensions
    n = length(reduction.x)
    # compute reduced system
    A = reduction.Dψ*reduction.P
    Q = reduction.P*inv(A)*reduction.Dψ
    Iₙ = diagonal_matrix(reduction.K(1), n)
    F¹ = matrix_space(reduction.K, n, 1)(reduction.f¹)
    f_red_raw = (Iₙ - Q)*F¹
    # reshape as vector
    f_red = reshape(Matrix(f_red_raw), n)
  else
    # reduction cannot be computed
    @error "Reduced system cannot be computed. You need to set a valid product decomposition, i.e. P and ψ such that f⁰ = P⋅ψ"
    return nothing, nothing
  end
  # Check if slow manifold is set 
  if reduction.success[1]
    # substitute components according to slow manifold
    a = reduction.K.([reduction.M; reduction.θ])
    f_red_subs = [evaluate(f, a) for f in f_red]
    return f_red, f_red_subs
  else
    @warn "Slow manifold has not been defined succesfully. Reduced system is only returned in raw form, i.e. the reduced components are not substituted according to the slow manfold."
    return f_red, nothing
  end
end


function compute_directional_reduction(reduction::Reduction, γ::Function)
  # check if curve satisfies γ(0) = π⁺
  @assert all(γ(reduction.θ, 0) .== reduction._θ) "The curve γ must satisfy γ(0) = π⁺"
  # compute γ'(0)
  Rδ, δ = polynomial_ring(reduction.R, :δ)
  dγ = derivative.(γ(reduction.θ, δ))
  dγ_0 = matrix(reduction.K ,length(reduction.θ), 1, evaluate.(dγ, 0))
  # compute D₂f(x, πˣ)
  D₂f = reduction.K.(eval_mat(jacobian(reduction.f, reduction.θ), reduction.θ, reduction._θ))
  # Check if P-ψ-composition is defined 
  if reduction.success[3]
    # check if non-singular point is defined
    if !reduction.success[2]
      @info "The non-singular point on the slow manifold has not been set successfully. Trying to compute the reduction anyway."
    end
    # dimensions
    n = length(reduction.x)
    # compute reduced system
    A = reduction.Dψ*reduction.P
    Q = reduction.P*inv(A)*reduction.Dψ
    Iₙ = diagonal_matrix(reduction.K(1), n)
    f_red_raw = (Iₙ - Q)*D₂f*dγ_0
    # reshape as vector
    f_red = reshape(Matrix(f_red_raw), n)
  else
    # reduction cannot be computed
    @error "Reduced system cannot be computed. You need to set a valid product decomposition, i.e. P and ψ such that f⁰ = P⋅ψ"
    return nothing, nothing
  end
  # Check if slow manifold is set 
  if reduction.success[1]
    # substitute components according to slow manifold
    a = reduction.K.([reduction.M; reduction.θ])
    f_red_subs = [evaluate(f, a) for f in f_red]
    return f_red, f_red_subs
  else
    @warn "Slow manifold has not been defined succesfully. Reduced system is only returned in raw form, i.e. the reduced components are not substituted according to the slow manfold."
    return f_red, nothing
  end
end
