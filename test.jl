using Oscar
import Oscar: IdealGens, AbstractAlgebra.Generic.MatSpaceElem

F, (α, β) = rational_function_field(QQ, [:α, :β]);
R, (x, y) = polynomial_ring(F, [:x, :y]);
I = ideal([2*x, α*x + β*y]);
G1 = Oscar.groebner_basis(I; complete_reduction=true)
G2, T2 = groebner_basis_with_transformation_matrix(I; complete_reduction=true);

function reduce_leading_coefficients(G::IdealGens, T::MatSpaceElem)
  G_gens = gens(G)
  lc = leading_coefficient.(G_gens)
  for j in eachindex(lc)
    if !isone(lc[j])
      for i in 1:size(T,2)
        T[i,j] *= 1//lc[j]
      end
      G_gens[j] *= 1//lc[j]
    end
  end
  _G = IdealGens(parent(G_gens[1]), G_gens, G.ord)
  _G.isGB = G.isGB
  _G.isReduced = G.isReduced
  return _G, T
end
function reduce_leading_coefficients(G::IdealGens)
  G_gens = gens(G)
  lc = leading_coefficient.(G_gens)
  for j in eachindex(lc)
    if !isone(lc[j])
      G_gens[j] *= 1//lc[j]
    end
  end
  _G = IdealGens(parent(G_gens[1]), G_gens, G.ord)
  _G.isGB = G.isGB
  _G.isReduced = G.isReduced
  return _G
end

_G1 = reduce_leading_coefficients(G1)
@time _G2, _T2 = reduce_leading_coefficients((G2, T2))
gens(I)*_T2 

GB = groebner_basis_with_transformation_matrix(I; complete_reduction=true)
GB2 = reduce_leading_coefficient(GB)

GB = groebner_basis_with_transformation_matrix(I; complete_reduction=true)

G, T = groebner_basis_with_transformation_matrix(I; complete_reduction=true, reduce_leading_coefficient=true)

G2
U,Q,H = reduce_with_quotients_and_unit(gens(I), gens(G1));
T1 = transpose(inv(U*Q))
gens(I)*T1 == gens(G1)




I = ideal([x, 2*x + 3*y])
G, T = groebner_basis_with_transformation_matrix(I; complete_reduction=true)
G

_, (x, y) = polynomial_ring(QQ, [:x, :y])
groebner_basis_with_transformation_matrix(ideal([x, 2*x + 3*y]); complete_reduction=true)[1]



