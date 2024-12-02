#ifndef A2D_ISOTROPIC_SHELL_STRESS_H
#define A2D_ISOTROPIC_SHELL_STRESS_H

#include <type_traits>

#include "../../a2ddefs.h"
#include "../a2dmat.h"
#include "../a2dstack.h"
#include "../a2dtest.h"
#include "../core/a2dsymmatveccore.h"

namespace A2D {

template <typename T, class Data>
A2D_FUNCTION void IsotropicShellStressCore(
  const T& E, const T& nu, const T& thick, const T& tOffset,
  T strain[], T stress[]
) {

  // stress and strain are length 9 vectors in pointer form now
  // need to compute A,B,D,As matrices etc
  T C[6];
  Data::evalTangentStiffness2D(E, nu, C);

  // use scale operations to not have to store A,B,D in memory (more memory efficient)
  // much cheaper than storing full 9x9 ABD matrix

  // Nij = A * Eij; A = C * thick
  A2D::SymMatVecCoreScale3x3<T,false>(thick, C, strain, stress);
  // Nij += B * kij; B = C * -thick * tOffset
  A2D::SymMatVecCoreScale3x3<T,true>(-thick * tOffset, C, &strain[3], stress);
  
  // Mij = B * Eij; B = C * -thick tOffset
  A2D::SymMatVecCoreScale3x3<T,false>(-thick * tOffset, C, strain, &stress[3]);
  // Mij += D * kij; D = t^3/12 * A + tOffset^2 * t^2 * A (the units on second term don't seem correct; error in TACS?)
  T I = thick*thick*thick/12.0;
  A2D::SymMatVecCoreScale3x3<T,true>(I, C, &strain[3], &stress[3]);
  A2D::SymMatVecCoreScale3x3<T,true>(thick*thick*tOffset*tOffset, C, &strain[3], &stress[3]);

  // transverse shear moments M13, M23
  T As = Data::transverseShearCorrectionFactor * thick * C[5];
  stress[6] = As * strain[6];
  stress[7] = As * strain[7];

  // drill stress M12
  T drill = Data::drillingRegularization * As;
  stress[8] = drill * strain[8]; 
}

template <typename T, class Data>
A2D_FUNCTION void IsotropicShellStressReverseCore(
  const T& E, const T& nu, const T& thick, const T& tOffset,
  T stressb[], T strainb[]
) {

  // stress and strain are length 9 vectors in pointer form now
  // need to compute A,B,D,As matrices etc
  T C[6];
  Data::evalTangentStiffness2D(E, nu, C);

  // use scale operations to not have to store A,B,D in memory (more memory efficient)
  // much cheaper than storing full 9x9 ABD matrix
  
  // reverse of y = A * x (each 3x3 submatrix part) is xb = A^T * yb (except A is sym so A^T = A)
  // and xb = A * yb (so we can use same A2D::SymMat operations in reverse no transpose one)

  // Nij = A * Eij; A = C * thick
  A2D::SymMatVecCoreScale3x3<T,false>(thick, C, stressb, strainb);
  // Nij += B * kij; B = C * -thick * tOffset
  A2D::SymMatVecCoreScale3x3<T,true>(-thick * tOffset, C, stressb, &strainb[3]);
  
  // Mij = B * Eij; B = C * -thick tOffset
  A2D::SymMatVecCoreScale3x3<T,false>(-thick * tOffset, C, &stressb[3], strainb);
  // Mij += D * kij; D = t^3/12 * A + tOffset^2 * t^2 * A (the units on second term don't seem correct; error in TACS?)
  T I = thick*thick*thick/12.0;
  A2D::SymMatVecCoreScale3x3<T,true>(I, C, &stressb[3], &strainb[3]);
  A2D::SymMatVecCoreScale3x3<T,true>(thick*thick*tOffset*tOffset, C, &stressb[3], &strainb[3]);

  // transverse shear moments M13, M23
  T As = Data::transverseShearCorrectionFactor * thick * C[5];
  strainb[6] = As * stressb[6];
  strainb[7] = As * stressb[7];

  // drill stress M12
  T drill = Data::drillingRegularization * As;
  strainb[8] = drill * stressb[8]; 
}

// TODO : hessian reverse core? prob the same as this one since linear

template <typename T, class Data>
A2D_FUNCTION void IsotropicShellStress(
  const T& E, const T &nu, const T& thick, const T& tOffset,
  const Vec<T,9>& strain, Vec<T,9>& stress) {
    
  IsotropicShellStressCore(E, nu, thick, tOffset, get_data(strain), get_data(stress));
}

template <class Data, class strainType, class stressType>
class IsotropicShellStressExpr {
 public:
  // Extract the numeric type to use
  typedef typename get_object_numeric_type<strainType>::type T;

  // Extract the dimensions of the underlying matrix
  static constexpr int nstrain = get_vec_size<strainType>::size;
  static constexpr int nstress = get_vec_size<strainType>::size;

  // make sure the correct sizes
  static_assert((nstrain == nstress),
                "Shell Stress Expression does not have right size..");

  static constexpr ADiffType adstrain = get_diff_type<strainType>::diff_type;
  static constexpr ADiffType adstress = get_diff_type<stressType>::diff_type;

  // Get the differentiation order from the output
  static constexpr ADorder order = get_diff_order<stressType>::order;

  // Make sure that the order matches
  static_assert(get_diff_order<strainType>::order == order,
                "ADorder does not match");

  // constructor
  A2D_FUNCTION IsotropicShellStressExpr(const T& E, const T& nu, const T& thick, const T& tOffset, 
    strainType& strain, stressType& stress) 
      : E(E), nu(nu), thick(thick), tOffset(tOffset), strain(strain), stress(stress) {}

  A2D_FUNCTION void eval() {
    IsotropicShellStressCore<T,Data>(E, nu, thick, tOffset, get_data(strain), get_data(stress));
  }

  A2D_FUNCTION void bzero() { stress.bzero(); }

  template <ADorder forder>
  A2D_FUNCTION void forward() {
    static_assert(
        !(order == ADorder::FIRST and forder == ADorder::SECOND),
        "Can't perform second order forward with first order objects");
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;

    IsotropicShellStressCore<T,Data>(
      E, nu, thick, tOffset, 
      GetSeed<seed>::get_data(strain), GetSeed<seed>::get_data(stress));
  }

  A2D_FUNCTION void reverse() {
    constexpr ADseed seed = ADseed::b;
    IsotropicShellStressReverseCore<T,Data>(
      E, nu, thick, tOffset, 
      GetSeed<seed>::get_data(stress), GetSeed<seed>::get_data(strain));
  }

  A2D_FUNCTION void hzero() { stress.hzero(); }

  A2D_FUNCTION void hreverse() {
    constexpr ADseed seed = ADseed::h;
    IsotropicShellStressReverseCore<T,Data>(
      E, nu, thick, tOffset, 
      GetSeed<seed>::get_data(stress), GetSeed<seed>::get_data(strain));
  }

  const T& E;
  const T& nu;
  const T& thick;
  const T& tOffset;
  strainType& strain;
  stressType& stress;
};

// fix this later
template <typename T, class Data>
A2D_FUNCTION auto IsotropicShellStress(
  const T& E, const T& nu, const T& thick, const T& tOffset,
  A2D::ADObj<A2D::Vec<T,9>>& strain, A2D::ADObj<A2D::Vec<T,9>>& stress) {

  return IsotropicShellStressExpr<Data, A2D::ADObj<A2D::Vec<T,9>>, A2D::ADObj<A2D::Vec<T,9>>>(E, nu, thick, tOffset, strain, stress);
}

// template <typename T, class Data, class strainType, class stressType>
// A2D_FUNCTION auto IsotropicShellStress(
//   T& E, T& nu, T& thick, T& tOffset,
//   A2DObj<strainType>& strain, A2DObj<stressType>& stress) {

//   return IsotropicShellStressExpr<Data, strainType, stressType>(E, nu, thick, tOffset, strain, stress);
// }

// TODO : make unittests for this
// namespace Test {
// }  // namespace Test

}  // namespace A2D

#endif  // A2D_ISOTROPIC_SHELL_STRESS_H
