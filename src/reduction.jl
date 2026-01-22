## Computing a reduction

"""
Type that holds all information to compute a Tikhonov-Fenichel reduction for a
given slow-fast separation of rates.

### Fields
    $(TYPEDFIELDS)

"""
mutable struct Reduction
  "information on input system"
  problem::ReductionProblem
  "slow-fast separation (0: slow, 1: fast)"
  sf_separation::Vector{Bool}
  "Parameters of the system where small ones are set to 0"
  _p::Vector{QQMPolyRingElem}
  "RHS of system as vector with elements of ring `R`"
  f0::Vector{QQMPolyRingElem}
  "Slow part of the system"
  f1::Vector{QQMPolyRingElem}
  "Second order terms in ε"
  higher_order_terms::Vector{QQMPolyRingElem}
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
  "Whether existence of reduction onto the given manifold can be excluded"
  no_reduction::Bool
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
  "Reduced system on the slow manifold "
  g::Vector{FracFieldElem{QQMPolyRingElem}}
  "Whether `g` and `g_raw` are already computed"
  reduction_cached::Vector{Bool}
end

# conversion functions
function R_to_Rx(problem::ReductionProblem, f::QQMPolyRingElem)
  _p = gens(base_ring(problem._Rx))
  _x = gens(problem._Rx)
  f_Rx = zero(problem._Rx)
  for (c,a) in coefficients_and_exponents(f)
    f_Rx += c*prod([_p; _x].^a)
  end
  return f_Rx
end
function Rx_to_F(problem::ReductionProblem, f::MPoly{RationalFunctionFieldElem{QQFieldElem, QQMPolyRingElem}})
  _p = gens(problem._F)[1:length(problem.p)]
  _x = gens(problem._F)[length(problem.p)+1:end]
  f_F = zero(problem._F)
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
Terms in O(ε²) are discarded.
"""
function split_system(f::Vector{QQMPolyRingElem}, p::Vector{QQMPolyRingElem}, sf_separation::Vector{Bool})
  R = parent(f[1])
  f0 = [evaluate(fᵢ, p[.!sf_separation], zero.(p[.!sf_separation])) for fᵢ in f]
  f1 = f .- f0 
  # remove all terms from f1 that are in O(ε²)
  r = [R(0) for _ in f1]
  # get indices of parameters considered for slow-fast separation
  idx_p = [any(v .== p) for v in gens(R)]
  idx_p_num = (1:ngens(R))[idx_p]
  idx_slow_all = idx_p_num[.!sf_separation]
  # store all terms in O(ε²) in 
  for i in eachindex(f1)
    ts = terms(f[i])
    for t in ts
      e = [exponent(t, 1, j) for j in idx_slow_all]
      if sum(e) > 1 
        r[i] += t
      end
    end
  end
  f1 = f1 .- r
  # make sure splitting worked
  @assert all(f0 .+ f1 .+ r .== f) "Something went wrong with splitting the system, please file a bug report at https://github.com/jo-ap/TikhonovFenichelReductions.jl"
  return f0, f1, r 
end

"""
    $(TYPEDSIGNATURES)

Compute Jacobian of `f` with respect to `x`.
"""
function jacobian(f::Vector{T}, x::Vector{T}) where T<:MPolyRingElem
  return matrix(parent(f[1]), [[derivative(fᵢ, xⱼ) for xⱼ in x] for fᵢ in f])
end
function jacobian(f::MatSpaceElem{T}, x::Vector{T}) where T<:MPolyRingElem
  @assert size(f, 2) == 1 "f must be n×1 Matrix"
  f = reshape(Matrix(f), size(f,1))
  return jacobian(f, x)
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
function get_varieties(problem::ReductionProblem, sf_separation::Vector{Bool})
  tfpv_candidate = get_tfpv(gens(problem._Fp), sf_separation)
  f0 = get_f0_Rx(problem, tfpv_candidate)
  J = jacobian(f0, problem._x_Rx)
  I = ideal(f0)
  _get_varieties(problem, I, J)
end

"""
    $(TYPEDSIGNATURES)
Constructor for `Reduction` Type.

### Arguments
- `problem`: Reduction problem type holding information on system and dimension of reduction.
- `sf_separation`: Boolean index indicating slow-fast separation of rates (0: small, 1: large).
- `s::Int`: (optional) Dimension of slow manifold. Can be specified if a reduction corresponding to a TFPV for dimension different from `problem.s` should be considered (e.g. for manually computing a reduction for a given slow-fast separation that is not necessarily obtained via `tfpvs_and_varieties`).

See also: [`set_manifold!`](@ref) [`set_decomposition!`](@ref)
"""
function Reduction(problem::ReductionProblem, sf_separation::Vector{Bool}; s::Int=problem.s)
  _p = problem.p .* sf_separation
  n = length(problem.x)
  r = n - s
  f0, f1, higher_order_terms = split_system(problem.f, problem.p, sf_separation)
  Df0 = jacobian(problem.f, problem.x)
  T, _ = polynomial_ring(problem._F, "λ")
  M = problem._F.(problem.x)
  x0 = zeros(problem._F, n)
  P = zero_matrix(problem._F,n,r)
  Psi = zero_matrix(problem._F,r,1)
  DPsi = zero_matrix(problem._F,r,n)
  Df0_at_x0 = matrix(problem._F, Matrix(Df0))
  return Reduction(problem,
                   sf_separation,
                   _p,
                   f0,
                   f1,
                   higher_order_terms,
                   Df0,
                   Df0_at_x0,
                   T,
                   T(0),
                   M,
                   false,
                   x0,
                   P,
                   Psi,
                   DPsi,
                   zeros(Bool, 3),
                   zeros(Bool, n),
                   zeros(problem._F, n),
                   zeros(problem._F, n),
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
- `s::Int`: (optional) Dimension of slow manifold. Can be specified if a reduction corresponding to a TFPV for dimension different from `problem.s` should be considered (e.g. for manually computing a reduction for a given slow-fast separation that is not necessarily obtained via `tfpvs_and_varieties`).

See also: [`set_manifold!`](@ref) [`set_decomposition!`](@ref) [`tfpvs_and_varieties`](@ref)
"""
function Reduction(
  problem::ReductionProblem,
  sf_separation::Vector{Bool},
  V::Union{Variety, Vector{QQMPolyRingElem}},
  M::Vector{<:RingElem}; 
  s::Int=problem.s)
  reduction = Reduction(problem, sf_separation; s=s)
  set_manifold!(reduction, M)
  set_decomposition!(reduction, V)
  if all(reduction.success)
    compute_reduction!(reduction)
  end
  return reduction
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
components, i.e. `reduction.problem.x`.

See also: [`Reduction`](@ref), [`set_decomposition!`](@ref), [`set_point!`](@ref)
"""
function set_manifold!(reduction::Reduction, M::AbstractVector)::Bool
  M = parse_ring(reduction.problem._F, M)
  n = length(reduction.problem.x)
  @assert length(M) == n "The slow manifold M must be defined in $n components."
  _f0 = [evaluate(fᵢ, [reduction.problem._F.(reduction.problem.p); M]) for fᵢ in reduction.f0]
  f_vanishes = all(iszero.(_f0))
  if !f_vanishes 
    @warn "f0 does no vanish on the slow manifold"
  else
    reduction.M = M
    reduction.success[1] = true
    reduction.idx_components = reduction.M .== reduction.problem.x
    # check if there exists a reduction
    # an eigenvalue λ has to factor the characteristic polynomial of the
    # jacobian at a non-singular point exactly with power s, i.e. there is no
    # reduction if the polynomial is given by λ^(s+1)•r(λ)
    coeffs = collect(coefficients(charpoly(jacobian_tfpv_on_manifold(reduction))))
    if all(coeffs[1:reduction.problem.s + 1] .== 0)
      @info "There exists no reduction onto this manifold"
      reduction.no_reduction = true
      return false
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
  eval_mat(reduction.Df0, [reduction.problem._F.(reduction._p); reduction.M])
end

"""
    $(TYPEDSIGNATURES)

Return the Jacobian `D₁f(x₀,π⁺)` at the point `x₀` on the slow manifold for a
TFPV `π⁺`.

See also: [`Reduction`](@ref), [`set_point!`](@ref), [`set_manifold!`](@ref)
"""
function jacobian_tfpv_at_x0(reduction::Reduction)
  eval_mat(reduction.Df0, [reduction.problem._F.(reduction._p); reduction.x0])
end

function _set_point!(reduction::Reduction, x0::AbstractVector)::Bool
  x0 = parse_ring(reduction.problem._F, x0)
  n = length(reduction.problem.x)
  @assert length(x0) == n "The point x0 must have $n components."
  # compute characteristic polynomial
  Df0_at_x0 = eval_mat(reduction.Df0, [reduction.problem._F.(reduction._p); x0])
  chi = charpoly(reduction.T, Df0_at_x0)
  # check condition for coefficients
  c = collect(coefficients(chi))
  check_chi = all(iszero.(c[1:reduction.problem.s])) && !iszero(c[reduction.problem.s+1])
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
    @warn "The eigenvalue λ does not factor the characteristic polynomial of D₁f(x0,p_sf) with power s=$(reduction.problem.s)"
  end
  return retval 
end
  
"""
    $(TYPEDSIGNATURES)

Set product decomposition `f0=P⋅Psi` locally satisfying `V(f0)=V(Psi)`, where `P`
is a matrix of rational functions  and `Psi` is a vector of polynomials.

See also: [`set_manifold!`](@ref) [`set_point!`](@ref) [`Reduction`](@ref)
"""
function set_decomposition!(
    reduction::Reduction,
    P::Union{MatSpaceElem,VecOrMat},
    Psi::Union{MatSpaceElem,Vector{QQMPolyRingElem}}
  )
  _set_decomposition!(reduction, P, Psi)
end

"""
    $(TYPEDSIGNATURES)

Try to automatically compute matrix of rational functions `P` from given vector
of polynomials `Psi`, such that `f0=P⋅Psi` and `V(f0)=V(Psi)` holds locally.

!!! warning 
    This always works if the drop in dimension `r=n-s=1`, but may fail for `r>1` (if the number of generators for the irreducible component of `V(f0)` is greater than `r`).

### Description
`Psi` can be chosen from `r` algebraically independent entries of `f0`. 
Practically, one can use the generators of the ideals defining the
irreducible components of `V(f0)` as entries for `Psi` (possibly rewriting the
rational equations as polynomials by multiplying appropriately with
parameters occurring in a denominator).
"""
function set_decomposition!(reduction::Reduction, Psi::Union{Variety,Vector{QQMPolyRingElem}})
  P, Psi = get_decomposition(reduction, Psi)
  isnothing(P) ? false : set_decomposition!(reduction, P, Psi)
end

function _set_decomposition!(reduction::Reduction, P::MatSpaceElem, Psi)
  n = length(reduction.problem.x)
  r = n - reduction.problem.s
  @assert size(Psi, 1) == r && size(Psi, 2) == 1 "Psi must be of size $r or $r×1"
  DPsi = jacobian(Psi, reduction.problem.x)
  DPsi = parent(reduction.DPsi)(reduction.problem._F.(DPsi))
  Psi = reduction.problem._F.(Psi)
  Psi = parent(reduction.Psi)(Psi)
  P = reduction.problem._F.(P)
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
  n = length(reduction.problem.x)
  r = n - reduction.problem.s
  @assert size(P,1) == n && size(P,2) == r "P must be of size $n×$r"
  P = reshape(P, n, r)
  P = reduction.problem._F.(P)
  P = parent(reduction.P)(P)
  _set_decomposition!(reduction, P, Psi)
end
  
# try computing matrix of rational functions P from Psi
function get_decomposition(reduction::Reduction, variety::Variety)
  R = parent(reduction.problem.x[1])
  r = length(reduction.problem.x) - reduction.problem.s
  p = gens(base_ring(reduction.problem._Rx))
  x = gens(reduction.problem._Rx)
  _f0 = reduction.problem._f(x, get_tfpv(p, reduction.sf_separation))
  # use generators for irreducible component as entries for Psi
  if length(variety.gens_R) == r 
    Psi = matrix(R, reshape(variety.gens_R, r, 1))
    U, Q, H = reduce_with_quotients_and_unit(_f0, variety.groebner_basis)
    if all(H .== 0)
      _P = U*Q*variety.T
      P = matrix([Rx_to_F(reduction.problem, f) for f in _P])
      return P, Psi
    end
  else
    # try to find r independent entries in f0 and use these as entries for Psi
    idx_independent = find_independent_polys(_f0)
    if sum(idx_independent) == r
      Psi = matrix(R, reshape(reduction.f0[idx_independent], r, 1))
      G, T = groebner_basis_with_transformation_matrix(ideal(_f0[idx_independent]); complete_reduction=true)
      _T = transpose(T)
      U, Q, H = reduce_with_quotients_and_unit(_f0, G)
      if all(H .== 0)
        _P = U*Q*_T
        P = matrix([Rx_to_F(reduction.problem, f) for f in _P])
        return P, Psi
      end
    end
  end
  @warn "Could not set product decomposition automatically."
  return nothing, nothing
end
function get_decomposition(reduction::Reduction, Psi::Vector{QQMPolyRingElem})
  r = length(reduction.problem.x) - reduction.problem.s
  @assert size(Psi, 1) == r "Psi must have length r=$r"
  if size(Psi, 1) == 1
    return reduction.f0.//Psi, Psi
  else 
    p = gens(base_ring(reduction.problem._Rx))
    x = gens(reduction.problem._Rx)
    _f0 = reduction.problem._f(x, p .* reduction.sf_separation)
    _Psi = [R_to_Rx(reduction.problem, p) for p in Psi]
    U, Q, H = reduce_with_quotients_and_unit(_f0, _Psi)
    if all(H .== 0)
      P = matrix(reduction.problem._F, [Rx_to_F(reduction.problem, f) for f in U*Q])
      return P, Psi
    end
    @warn "Could not set P automatically."
    return nothing, Psi
  end
end

function find_independent_polys(f::Vector{<:MPolyRingElem})
  idx = [false for _ in f]
  for i in eachindex(f)
    if is_algebraically_independent([f[idx]..., f[i]])
      idx[i] = true
    end 
  end 
  return idx
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

See also: [`set_manifold!`](@ref), [`set_decomposition!`](@ref), [`set_point!`](@ref), [`compute_reductions`](@ref), [`Reduction`](@ref), [`print_reduced_system`](@ref)
"""
function compute_reduction!(reduction::Reduction)
  # Check if P-Psi-composition is defined 
  if reduction.success[3]
    # check if non-singular point is defined
    if !reduction.success[2]
      @info "The non-singular point on the slow manifold has not been set successfully. Trying to compute the reduction anyway."
    end
    # dimensions
    n = length(reduction.problem.x)
    # compute reduced system
    A = reduction.DPsi*reduction.P
    Q = reduction.P*inv(A)*reduction.DPsi
    Iₙ = diagonal_matrix(reduction.problem._F(1), n)
    f1 = matrix_space(reduction.problem._F, n, 1)(reduction.f1)
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
    a = reduction.problem._F.([reduction.problem.p; reduction.M])
    f_red_subs = [evaluate(f, a) for f in f_red]
    reduction.g = f_red_subs
    reduction.reduction_cached[2] = true
  else
    @warn "Slow manifold has not been defined succesfully. Reduced system is only returned in raw form, i.e. the reduced components are not substituted according to the slow manfold."
  end
  return true
end

function _get_slow_fast(reduction::Reduction)
  sf_separation = reduction.sf_separation
  slow = reduction.problem.p[.!sf_separation]
  fast = reduction.problem.p[sf_separation]
  return slow, fast
end

function is_linear(p::MPolyRingElem, i::Int)
  e = [exponent(p, k, i) for k=1:length(p)]
  return any(e .> 0) && all(e .<= 1)
end

# function solve_linear(p::Union{MPolyRingElem, AbstractAlgebra.Generic.FracFieldElem}, i::Int)
function solve_linear(p::MPolyRingElem, i::Int)
  R = parent(p)
  if is_linear(p, i)
    # p = s + t = 0
    t = evaluate(p, [i], R.([0]))
    s = p - t
    # s = c*v = -t => v = -t//c
    c = evaluate(s, [i], R.([1]))
    return -t//c
  end
  return nothing
end

function occurs_in_poly(p::T, v::Vector{T}) where T<:MPolyRingElem 
  _vars = vars(p)
  return [x in _vars for x in v]
end
function occurs_in_poly(p::AbstractAlgebra.Generic.FracFieldElem{T}, v::Vector{T}) where T<:MPolyRingElem
  x_occurs_num = occurs_in_poly(numerator(p), v) 
  x_occurs_den = occurs_in_poly(denominator(p), v)
  return x_occurs_num .|| x_occurs_den
end

function is_correct_manifold(problem, M, variety)
  # check if f^{(0)}(x)=0 for x ∈ M 
  for p in variety.groebner_basis_R 
    if evaluate(p, problem.x, M) != 0
      return false
    end
  end
  # check if the number of local coordinates is correct
  v = zeros(Int, length(M))
  for m in M 
    v .+= occurs_in_poly(m, problem.x)
  end
  return sum(v .> 0) == problem.s
end

"""
    $(TYPEDSIGNATURES)

Heuristic approach to get an explicit (i.e. parameterized) representation of
the slow manifold from a variety.
If succesfull, the function returns the (attempted) explicit manifold together
with a boolean value indicating whether the manifold could be computed
automatically.
"""
function get_explicit_manifold(problem::ReductionProblem, variety::Variety)
  R = parent(problem.x[1])
  # index of state variables in ring R
  idx_x = (length(problem.p)+1):(length(problem.p) + length(problem.x))
  M = fraction_field(R).(problem.x)
  G = sort(variety.groebner_basis_R; by=length)
  G_old = zeros(R, length(G))
  # try solving polynomials in G that are linear in one component until s local
  # coordinates are found
  while !all(G .== G_old)
    # halt if G does not change anymore
    G_old = G 
    for k in eachindex(G)
      # get all variables that occur
      v = [occurs_in_poly(g, problem.x) for g in G]
      if any(v[k]) 
        # check if variable i can be used to solve a polynomial
        for i in idx_x[v[k]]
          if is_linear(G[k], i)
            val = solve_linear(G[k], i)
            M = [evaluate(m, [i], [val]) for m in M]
            # keep only numerators of resulting rational functions
            G = [numerator(evaluate(g, [i], [val])) for g in G]
            break
          end
        end
      end
    end
    if all(G .== 0)
      break
    end
  end
  return M, is_correct_manifold(problem, M, variety)
end

