#ifndef ELASTICITY_3D_H
#define ELASTICITY_3D_H

#include "a2dtmp.h"
#include <cstddef>

/*
  Scatter the variables stored at the nodes to the a data structure
  that stores the element variables at each element node.
*/
template <class ConnArray,
          class NodeArray,
          class ElementArray>
void element_scatter( ConnArray& conn,
                      NodeArray& X,
                      ElementArray& Xe ){
  // Loop over the elements
  for ( std::size_t i = 0; i < conn.extent(0); i++ ){

    // Loop over each element nodes
    for ( std::size_t j = 0; j < conn.extent(1); j++ ){
      const std::size_t index = conn(i, j);

      // Loop over the variables
      for ( std::size_t k = 0; k < X.extent(1); k++ ){
        Xe(i, j, k) = X(index, k);
      }
    }
  }
}

class NonlinearElasticity3D {
public:
  static const int NUM_VARS = 3;

  template<typename T, class IdxType, class QuadPointData>
  static T compute_energy( IdxType i, IdxType j, QuadPointData& data,
                           A2D::Mat<T, 3, 3>& Jinv,
                           A2D::Mat<T, 3, 3>& Uxi ){
    typedef A2D::SymmMat<T, 3> SymmMat3x3;
    typedef A2D::Mat<T, 3, 3> Mat3x3;

    T mu(data(i, j, 0)), lambda(data(i, j, 1));
    Mat3x3 Ux;
    SymmMat3x3 E, S;
    T output;

    A2D::Mat3x3MatMult(Uxi, Jinv, Ux);
    A2D::Mat3x3GreenStrain(Ux, E);
    A2D::Symm3x3IsotropicConstitutive(mu, lambda, E, S);
    A2D::Symm3x3SymmMultTrace(S, E, output);

    return output;
  }

  template<typename T, class IdxType, class QuadPointData>
  static T compute_residual( IdxType i, IdxType j, QuadPointData& data,
                             T wdetJ,
                             A2D::Mat<T, 3, 3>& Jinv,
                             A2D::Mat<T, 3, 3>& Uxi0,
                             A2D::Mat<T, 3, 3>& Uxib ){
    typedef A2D::SymmMat<T, 3> SymmMat3x3;
    typedef A2D::Mat<T, 3, 3> Mat3x3;

    T mu(data(i, j, 0)), lambda(data(i, j, 1));
    Mat3x3 Ux0, Uxb;
    SymmMat3x3 E0, Eb;
    SymmMat3x3 S0, Sb;

    A2D::ADMat<Mat3x3> Uxi(Uxi0, Uxib);
    A2D::ADMat<Mat3x3> Ux(Ux0, Uxb);
    A2D::ADMat<SymmMat3x3> S(S0, Sb);
    A2D::ADMat<SymmMat3x3> E(E0, Eb);
    A2D::ADScalar<T> output;

    auto mult = A2D::Mat3x3MatMult(Uxi, Jinv, Ux);
    auto strain = A2D::Mat3x3GreenStrain(Ux, E);
    auto constitutive = A2D::Symm3x3IsotropicConstitutive(mu, lambda, E, S);
    auto trace = A2D::Symm3x3SymmMultTrace(S, E, output);

    output.bvalue = wdetJ;

    trace.reverse();
    constitutive.reverse();
    strain.reverse();
    mult.reverse();

    return output.value;
  }

  template<typename T, class IdxType, class QuadPointData>
  static T compute_jacobian( IdxType i, IdxType j, QuadPointData& data,
                             T wdetJ,
                             A2D::Mat<T, 3, 3>& Jinv,
                             A2D::Mat<T, 3, 3>& Uxi0,
                             A2D::Mat<T, 3, 3>& Uxib,
                             A2D::SymmTensor<T, 3, 3>& jac ){
    typedef A2D::SymmMat<T, 3> SymmMat3x3;
    typedef A2D::Mat<T, 3, 3> Mat3x3;

    T mu(data(i, j, 0)), lambda(data(i, j, 1));
    Mat3x3 Ux0, Uxb, Uxp, Uxh;
    Mat3x3 Uxip, Uxih;
    SymmMat3x3 E0, Eb, Ep, Eh;
    SymmMat3x3 S0, Sb, Sp, Sh;

    A2D::A2DMat<Mat3x3> Uxi(Uxi0, Uxib, Uxip, Uxih);
    A2D::A2DMat<Mat3x3> Ux(Ux0, Uxb, Uxp, Uxh);
    A2D::A2DMat<SymmMat3x3> S(S0, Sb, Sp, Sh);
    A2D::A2DMat<SymmMat3x3> E(E0, Eb, Ep, Eh);
    A2D::A2DScalar<T> output;

    auto mult = A2D::Mat3x3MatMult(Uxi, Jinv, Ux);
    auto strain = A2D::Mat3x3GreenStrain(Ux, E);
    auto constitutive = A2D::Symm3x3IsotropicConstitutive(mu, lambda, E, S);
    auto trace = A2D::Symm3x3SymmMultTrace(S, E, output);

    output.bvalue = wdetJ;

    trace.reverse();
    constitutive.reverse();
    strain.reverse();
    mult.reverse();

    for ( int i = 0; i < 3; i++ ){
      for ( int j = 0; j < 3; j++ ){
        // Zero the derivatives
        Uxih.zero();
        Uxh.zero();
        Sh.zero();
        Eh.zero();

        Uxip.zero();
        Uxip(i, j) = 1.0;

        mult.hforward();
        strain.hforward();
        constitutive.hforward();
        trace.hreverse();
        constitutive.hreverse();
        strain.hreverse();
        mult.hreverse();

        for ( int k = 0; k < 3; k++ ){
          for ( int l = 0; l < 3; l++ ){
            jac(i, j, k, l) = Uxih(k, l);
          }
        }
      }
    }

    return output.value;
  }
};


class LinearElasticity3D {
public:
  static const int NUM_VARS = 3;

  template<typename T, class IdxType, class QuadPointData>
  static T compute_energy( IdxType i, IdxType j, QuadPointData& data,
                           A2D::Mat<T, 3, 3>& Jinv,
                           A2D::Mat<T, 3, 3>& Uxi ){
    typedef A2D::SymmMat<T, 3> SymmMat3x3;
    typedef A2D::Mat<T, 3, 3> Mat3x3;

    T mu(data(i, j, 0)), lambda(data(i, j, 1));
    Mat3x3 Ux;
    SymmMat3x3 E, S;
    T output;

    A2D::Mat3x3MatMult(Uxi, Jinv, Ux);
    A2D::Mat3x3LinearGreenStrain(Ux, E);
    A2D::Symm3x3IsotropicConstitutive(mu, lambda, E, S);
    A2D::Symm3x3SymmMultTrace(S, E, output);

    return output;
  }

  template<typename T, class IdxType, class QuadPointData>
  static T compute_residual( IdxType i, IdxType j, QuadPointData& data,
                             T wdetJ,
                             A2D::Mat<T, 3, 3>& Jinv,
                             A2D::Mat<T, 3, 3>& Uxi0,
                             A2D::Mat<T, 3, 3>& Uxib ){
    typedef A2D::SymmMat<T, 3> SymmMat3x3;
    typedef A2D::Mat<T, 3, 3> Mat3x3;

    T mu(data(i, j, 0)), lambda(data(i, j, 1));
    Mat3x3 Ux0, Uxb;
    SymmMat3x3 E0, Eb;
    SymmMat3x3 S0, Sb;

    A2D::ADMat<Mat3x3> Uxi(Uxi0, Uxib);
    A2D::ADMat<Mat3x3> Ux(Ux0, Uxb);
    A2D::ADMat<SymmMat3x3> S(S0, Sb);
    A2D::ADMat<SymmMat3x3> E(E0, Eb);
    A2D::ADScalar<T> output;

    auto mult = A2D::Mat3x3MatMult(Uxi, Jinv, Ux);
    auto strain = A2D::Mat3x3LinearGreenStrain(Ux, E);
    auto constitutive = A2D::Symm3x3IsotropicConstitutive(mu, lambda, E, S);
    auto trace = A2D::Symm3x3SymmMultTrace(S, E, output);

    output.bvalue = wdetJ;

    trace.reverse();
    constitutive.reverse();
    strain.reverse();
    mult.reverse();

    return output.value;
  }

  template<typename T, class IdxType, class QuadPointData>
  static T compute_jacobian( IdxType i, IdxType j, QuadPointData& data,
                             T wdetJ,
                             A2D::Mat<T, 3, 3>& Jinv,
                             A2D::Mat<T, 3, 3>& Uxi0,
                             A2D::Mat<T, 3, 3>& Uxib,
                             A2D::SymmTensor<T, 3, 3>& jac ){
    typedef A2D::SymmMat<T, 3> SymmMat3x3;
    typedef A2D::Mat<T, 3, 3> Mat3x3;

    T mu(data(i, j, 0)), lambda(data(i, j, 1));
    Mat3x3 Ux0, Uxb, Uxp, Uxh;
    Mat3x3 Uxip, Uxih;
    SymmMat3x3 E0, Eb, Ep, Eh;
    SymmMat3x3 S0, Sb, Sp, Sh;

    A2D::A2DMat<Mat3x3> Uxi(Uxi0, Uxib, Uxip, Uxih);
    A2D::A2DMat<Mat3x3> Ux(Ux0, Uxb, Uxp, Uxh);
    A2D::A2DMat<SymmMat3x3> S(S0, Sb, Sp, Sh);
    A2D::A2DMat<SymmMat3x3> E(E0, Eb, Ep, Eh);
    A2D::A2DScalar<T> output;

    auto mult = A2D::Mat3x3MatMult(Uxi, Jinv, Ux);
    auto strain = A2D::Mat3x3LinearGreenStrain(Ux, E);
    auto constitutive = A2D::Symm3x3IsotropicConstitutive(mu, lambda, E, S);
    auto trace = A2D::Symm3x3SymmMultTrace(S, E, output);

    output.bvalue = wdetJ;

    trace.reverse();
    constitutive.reverse();
    strain.reverse();
    mult.reverse();

    for ( int i = 0; i < 3; i++ ){
      for ( int j = 0; j < 3; j++ ){
        // Zero the derivatives
        Uxih.zero();
        Uxh.zero();
        Sh.zero();
        Eh.zero();

        Uxip.zero();
        Uxip(i, j) = 1.0;

        mult.hforward();
        strain.hforward();
        constitutive.hforward();
        trace.hreverse();
        constitutive.hreverse();
        strain.hreverse();
        mult.hreverse();

        for ( int k = 0; k < 3; k++ ){
          for ( int l = 0; l < 3; l++ ){
            jac(i, j, k, l) = Uxih(k, l);
          }
        }
      }
    }

    return output.value;
  }
};


const double GaussQuadPts2[] = { -0.577350269189626, 0.577350269189626 };
const double GaussQuadWts2[] = {  1.0,               1.0 };

class HexQuadrature {
public:
  static const int NUM_QUAD_PTS = 8;

  static int getNumQuadPoints(){ return NUM_QUAD_PTS; }
  static void getQuadPoint( const int index, double pt[] ){
    pt[0] = GaussQuadPts2[index % 2];
    pt[1] = GaussQuadPts2[(index % 4)/2];
    pt[2] = GaussQuadPts2[index/4];
  }
  static double getQuadWeight( const int index ){
    return (GaussQuadWts2[index % 2] *
            GaussQuadWts2[(index % 4)/2] *
            GaussQuadWts2[index/4]);
  }
};

template <class Quadrature>
class HexBasis {
public:
  static const int NUM_NODES = 8;

  /*
    Interpolate from element-oriented data to quadrature point data
  */
  template<const int num_vars,
           class ElementArray,
           class QuadPointArray>
  static void interp( ElementArray& input,
                      QuadPointArray& output ){
    for ( std::size_t j = 0; j < Quadrature::getNumQuadPoints(); j++ ){
      double pt[3];
      Quadrature::getQuadPoint(j, pt);

      double n1[2], n2[2], n3[2];
      n1[0] = 0.5*(1.0 - pt[0]);
      n1[1] = 0.5*(1.0 + pt[0]);
      n2[0] = 0.5*(1.0 - pt[1]);
      n2[1] = 0.5*(1.0 + pt[1]);
      n3[0] = 0.5*(1.0 - pt[2]);
      n3[1] = 0.5*(1.0 + pt[2]);

      for ( std::size_t i = 0; i < input.extent(0); i++ ){
        for ( int ii = 0; ii < num_vars; ii++ ){
          output(i, j, ii) =
            n3[0] * (n2[0] * (n1[0] * input(i, 0, ii) + n1[1] * input(i, 1, ii)) +
                     n2[1] * (n1[0] * input(i, 2, ii) + n1[1] * input(i, 3, ii))) +
            n3[1] * (n2[0] * (n1[0] * input(i, 4, ii) + n1[1] * input(i, 5, ii)) +
                     n2[1] * (n1[0] * input(i, 6, ii) + n1[1] * input(i, 7, ii)));
        }
      }
    }
  }

  /*
    Compute the Jacobian transformation at each quadrature point
  */
  template<typename T,
           class ElementNodeArray,
           class QuadPointDetJArray,
           class QuadPointJacobianArray>
  static void compute_jtrans( ElementNodeArray& X,
                              QuadPointDetJArray& detJ,
                              QuadPointJacobianArray& Jinv ){
    for ( std::size_t j = 0; j < Quadrature::getNumQuadPoints(); j++ ){
      double pt[3];
      Quadrature::getQuadPoint(j, pt);

      double n1[2], n2[2], n3[2];
      n1[0] = 0.5*(1.0 - pt[0]);
      n1[1] = 0.5*(1.0 + pt[0]);
      n2[0] = 0.5*(1.0 - pt[1]);
      n2[1] = 0.5*(1.0 + pt[1]);
      n3[0] = 0.5*(1.0 - pt[2]);
      n3[1] = 0.5*(1.0 + pt[2]);

      double n1x[2], n2x[2], n3x[2];
      n1x[0] = -0.5;
      n1x[1] = 0.5;
      n2x[0] = -0.5;
      n2x[1] = 0.5;
      n3x[0] = -0.5;
      n3x[1] = 0.5;

      for ( std::size_t i = 0; i < X.extent(0); i++ ){
        // Compute the Jacobian transformation
        A2D::Mat<T, 3, 3> J;
        for ( int ii = 0; ii < 3; ii++ ){
          J(ii, 0) =
            n3[0] * (n2[0] * (n1x[0] * X(i, 0, ii) + n1x[1] * X(i, 1, ii)) +
                     n2[1] * (n1x[0] * X(i, 2, ii) + n1x[1] * X(i, 3, ii))) +
            n3[1] * (n2[0] * (n1x[0] * X(i, 4, ii) + n1x[1] * X(i, 5, ii)) +
                     n2[1] * (n1x[0] * X(i, 6, ii) + n1x[1] * X(i, 7, ii)));

          J(ii, 1) =
            n3[0] * (n2x[0] * (n1[0] * X(i, 0, ii) + n1[1] * X(i, 1, ii)) +
                     n2x[1] * (n1[0] * X(i, 2, ii) + n1[1] * X(i, 3, ii))) +
            n3[1] * (n2x[0] * (n1[0] * X(i, 4, ii) + n1[1] * X(i, 5, ii)) +
                     n2x[1] * (n1[0] * X(i, 6, ii) + n1[1] * X(i, 7, ii)));

          J(ii, 2) =
            n3x[0] * (n2[0] * (n1[0] * X(i, 0, ii) + n1[1] * X(i, 1, ii)) +
                      n2[1] * (n1[0] * X(i, 2, ii) + n1[1] * X(i, 3, ii))) +
            n3x[1] * (n2[0] * (n1[0] * X(i, 4, ii) + n1[1] * X(i, 5, ii)) +
                      n2[1] * (n1[0] * X(i, 6, ii) + n1[1] * X(i, 7, ii)));
        }

        // Compute the 3x3 matrix inverse
        A2D::Mat<T, 3, 3> jinv;
        A2D::Mat3x3Inverse(J, jinv);
        A2D::Mat3x3Det(J, detJ(i, j));

        // Copy values of the inverse of the Jacobian
        for ( int ii = 0; ii < 3; ii++ ){
          for ( int jj = 0; jj < 3; jj++ ){
            Jinv(i, j, ii, jj) = jinv(ii, jj);
          }
        }
      }
    }
  }

  template<typename T,
           const int num_vars,
           class ElementSolutionArray,
           class QuadPointGradientArray>
  static void gradient( ElementSolutionArray& U,
                        QuadPointGradientArray& Uxi ){
    for ( std::size_t j = 0; j < Quadrature::getNumQuadPoints(); j++ ){
      double pt[3];
      Quadrature::getQuadPoint(j, pt);

      double n1[2], n2[2], n3[2];
      n1[0] = 0.5*(1.0 - pt[0]);
      n1[1] = 0.5*(1.0 + pt[0]);
      n2[0] = 0.5*(1.0 - pt[1]);
      n2[1] = 0.5*(1.0 + pt[1]);
      n3[0] = 0.5*(1.0 - pt[2]);
      n3[1] = 0.5*(1.0 + pt[2]);

      double n1x[2], n2x[2], n3x[2];
      n1x[0] = -0.5;
      n1x[1] = 0.5;
      n2x[0] = -0.5;
      n2x[1] = 0.5;
      n3x[0] = -0.5;
      n3x[1] = 0.5;

      for ( std::size_t i = 0; i < U.extent(0); i++ ){
        for ( int ii = 0; ii < num_vars; ii++ ){
          Uxi(i, j, ii, 0) =
            n3[0] * (n2[0] * (n1x[0] * U(i, 0, ii) + n1x[1] * U(i, 1, ii)) +
                     n2[1] * (n1x[0] * U(i, 2, ii) + n1x[1] * U(i, 3, ii))) +
            n3[1] * (n2[0] * (n1x[0] * U(i, 4, ii) + n1x[1] * U(i, 5, ii)) +
                     n2[1] * (n1x[0] * U(i, 6, ii) + n1x[1] * U(i, 7, ii)));

          Uxi(i, j, ii, 1) =
            n3[0] * (n2x[0] * (n1[0] * U(i, 0, ii) + n1[1] * U(i, 1, ii)) +
                     n2x[1] * (n1[0] * U(i, 2, ii) + n1[1] * U(i, 3, ii))) +
            n3[1] * (n2x[0] * (n1[0] * U(i, 4, ii) + n1[1] * U(i, 5, ii)) +
                     n2x[1] * (n1[0] * U(i, 6, ii) + n1[1] * U(i, 7, ii)));

          Uxi(i, j, ii, 2) =
            n3x[0] * (n2[0] * (n1[0] * U(i, 0, ii) + n1[1] * U(i, 1, ii)) +
                      n2[1] * (n1[0] * U(i, 2, ii) + n1[1] * U(i, 3, ii))) +
            n3x[1] * (n2[0] * (n1[0] * U(i, 4, ii) + n1[1] * U(i, 5, ii)) +
                      n2[1] * (n1[0] * U(i, 6, ii) + n1[1] * U(i, 7, ii)));
        }
      }
    }
  }

  template<typename T,
           class Model,
           class QuadPointModelDataArray,
           class QuadPointDetJArray,
           class QuadPointJacobianArray,
           class QuadPointGradientArray>
  static void energy( QuadPointModelDataArray& Edata,
                      QuadPointDetJArray& detJ,
                      QuadPointJacobianArray& Jinv,
                      QuadPointGradientArray& Uxi,
                      T& energy ){
    for ( std::size_t j = 0; j < Quadrature::getNumQuadPoints(); j++ ){
      double weight = Quadrature::getQuadWeight(j);

      for ( std::size_t i = 0; i < Uxi.extent(0); i++ ){
        // Extract Jinv
        A2D::Mat<T, 3, 3> Jinv0;
        for ( int ii = 0; ii < 3; ii++ ){
          for ( int jj = 0; jj < 3; jj++ ){
            Jinv0(ii, jj) = Jinv(i, j, ii, jj);
          }
        }

        // Extract Uxi0
        A2D::Mat<T, Model::NUM_VARS, 3> Uxi0;
        for ( int ii = 0; ii < Model::NUM_VARS; ii++ ){
          for ( int jj = 0; jj < 3; jj++ ){
            Uxi0(ii, jj) = Uxi(i, j, ii, jj);
          }
        }

        energy += weight * detJ(i, j) * Model::compute_energy(i, j, Edata, Jinv0, Uxi0);
      }
    }
  }

  template<typename T,
           class Model,
           class QuadPointModelDataArray,
           class QuadPointDetJArray,
           class QuadPointJacobianArray,
           class QuadPointGradientArray,
           class ElementResidualArray>
  static void residuals( QuadPointModelDataArray& Edata,
                         QuadPointDetJArray& detJ,
                         QuadPointJacobianArray& Jinv,
                         QuadPointGradientArray& Uxi,
                         ElementResidualArray& res ){
    for ( std::size_t j = 0; j < Quadrature::getNumQuadPoints(); j++ ){
      double pt[3];
      Quadrature::getQuadPoint(j, pt);
      double weight = Quadrature::getQuadWeight(j);

      double n1[2], n2[2], n3[2];
      n1[0] = 0.5*(1.0 - pt[0]);
      n1[1] = 0.5*(1.0 + pt[0]);
      n2[0] = 0.5*(1.0 - pt[1]);
      n2[1] = 0.5*(1.0 + pt[1]);
      n3[0] = 0.5*(1.0 - pt[2]);
      n3[1] = 0.5*(1.0 + pt[2]);

      double n1x[2], n2x[2], n3x[2];
      n1x[0] = -0.5;
      n1x[1] = 0.5;
      n2x[0] = -0.5;
      n2x[1] = 0.5;
      n3x[0] = -0.5;
      n3x[1] = 0.5;

      double N1x[8], N2x[8], N3x[8];
      N1x[0] = n3[0] * n2[0] * n1x[0];
      N1x[1] = n3[0] * n2[0] * n1x[1];
      N1x[2] = n3[0] * n2[1] * n1x[0];
      N1x[3] = n3[0] * n2[1] * n1x[1];
      N1x[4] = n3[1] * n2[0] * n1x[0];
      N1x[5] = n3[1] * n2[0] * n1x[1];
      N1x[6] = n3[1] * n2[1] * n1x[0];
      N1x[7] = n3[1] * n2[1] * n1x[1];

      N2x[0] = n3[0] * n2x[0] * n1[0];
      N2x[1] = n3[0] * n2x[0] * n1[1];
      N2x[2] = n3[0] * n2x[1] * n1[0];
      N2x[3] = n3[0] * n2x[1] * n1[1];
      N2x[4] = n3[1] * n2x[0] * n1[0];
      N2x[5] = n3[1] * n2x[0] * n1[1];
      N2x[6] = n3[1] * n2x[1] * n1[0];
      N2x[7] = n3[1] * n2x[1] * n1[1];

      N3x[0] = n3x[0] * n2[0] * n1[0];
      N3x[1] = n3x[0] * n2[0] * n1[1];
      N3x[2] = n3x[0] * n2[1] * n1[0];
      N3x[3] = n3x[0] * n2[1] * n1[1];
      N3x[4] = n3x[1] * n2[0] * n1[0];
      N3x[5] = n3x[1] * n2[0] * n1[1];
      N3x[6] = n3x[1] * n2[1] * n1[0];
      N3x[7] = n3x[1] * n2[1] * n1[1];

      for ( std::size_t i = 0; i < Uxi.extent(0); i++ ){
        A2D::Mat<T, 3, 3> Jinv0;
        A2D::Mat<T, Model::NUM_VARS, 3> Uxi0, Uxib;

        // Extract Jinv
        for ( int ii = 0; ii < 3; ii++ ){
          for ( int jj = 0; jj < 3; jj++ ){
            Jinv0(ii, jj) = Jinv(i, j, ii, jj);
          }
        }

        // Extract Uxi0
        for ( int ii = 0; ii < Model::NUM_VARS; ii++ ){
          for ( int jj = 0; jj < 3; jj++ ){
            Uxi0(ii, jj) = Uxi(i, j, ii, jj);
          }
        }

        Model::compute_residual(i, j, Edata, weight * detJ(i, j), Jinv0, Uxi0, Uxib);

        for ( int ii = 0; ii < Model::NUM_VARS; ii++ ){
          for ( int k = 0; k < 8; k++ ){
            res(i, k, ii) +=
              N1x[k] * Uxib(ii, 0) + N2x[k] * Uxib(ii, 1) + N3x[k] * Uxib(ii, 2);
          }
        }
      }
    }
  }

  template<typename T,
           class Model,
           class QuadPointModelDataArray,
           class QuadPointDetJArray,
           class QuadPointJacobianArray,
           class QuadPointGradientArray,
           class ElementResidualArray>
  static void jacobians( QuadPointModelDataArray& Edata,
                         QuadPointDetJArray& detJ,
                         QuadPointJacobianArray& Jinv,
                         QuadPointGradientArray& Uxi,
                         ElementResidualArray& jac ){
    for ( std::size_t j = 0; j < Quadrature::getNumQuadPoints(); j++ ){
      double pt[3];
      Quadrature::getQuadPoint(j, pt);
      double weight = Quadrature::getQuadWeight(j);

      double n1[2], n2[2], n3[2];
      n1[0] = 0.5*(1.0 - pt[0]);
      n1[1] = 0.5*(1.0 + pt[0]);
      n2[0] = 0.5*(1.0 - pt[1]);
      n2[1] = 0.5*(1.0 + pt[1]);
      n3[0] = 0.5*(1.0 - pt[2]);
      n3[1] = 0.5*(1.0 + pt[2]);

      double n1x[2], n2x[2], n3x[2];
      n1x[0] = -0.5;
      n1x[1] = 0.5;
      n2x[0] = -0.5;
      n2x[1] = 0.5;
      n3x[0] = -0.5;
      n3x[1] = 0.5;

      double N1x[8], N2x[8], N3x[8];
      N1x[0] = n3[0] * n2[0] * n1x[0];
      N1x[1] = n3[0] * n2[0] * n1x[1];
      N1x[2] = n3[0] * n2[1] * n1x[0];
      N1x[3] = n3[0] * n2[1] * n1x[1];
      N1x[4] = n3[1] * n2[0] * n1x[0];
      N1x[5] = n3[1] * n2[0] * n1x[1];
      N1x[6] = n3[1] * n2[1] * n1x[0];
      N1x[7] = n3[1] * n2[1] * n1x[1];

      N2x[0] = n3[0] * n2x[0] * n1[0];
      N2x[1] = n3[0] * n2x[0] * n1[1];
      N2x[2] = n3[0] * n2x[1] * n1[0];
      N2x[3] = n3[0] * n2x[1] * n1[1];
      N2x[4] = n3[1] * n2x[0] * n1[0];
      N2x[5] = n3[1] * n2x[0] * n1[1];
      N2x[6] = n3[1] * n2x[1] * n1[0];
      N2x[7] = n3[1] * n2x[1] * n1[1];

      N3x[0] = n3x[0] * n2[0] * n1[0];
      N3x[1] = n3x[0] * n2[0] * n1[1];
      N3x[2] = n3x[0] * n2[1] * n1[0];
      N3x[3] = n3x[0] * n2[1] * n1[1];
      N3x[4] = n3x[1] * n2[0] * n1[0];
      N3x[5] = n3x[1] * n2[0] * n1[1];
      N3x[6] = n3x[1] * n2[1] * n1[0];
      N3x[7] = n3x[1] * n2[1] * n1[1];

      for ( std::size_t i = 0; i < Uxi.extent(0); i++ ){
        A2D::Mat<T, 3, 3> Jinv0;
        A2D::Mat<T, Model::NUM_VARS, 3> Uxi0, Uxib;

        // Extract Jinv
        for ( int ii = 0; ii < 3; ii++ ){
          for ( int jj = 0; jj < 3; jj++ ){
            Jinv0(ii, jj) = Jinv(i, j, ii, jj);
          }
        }

        // Extract Uxi0
        for ( int ii = 0; ii < Model::NUM_VARS; ii++ ){
          for ( int jj = 0; jj < 3; jj++ ){
            Uxi0(ii, jj) = Uxi(i, j, ii, jj);
          }
        }

        // The Jacobian of the energy
        A2D::SymmTensor<T, Model::NUM_VARS, 3> ja;

        Model::compute_jacobian(i, j, Edata, weight * detJ(i, j), Jinv0, Uxi0, Uxib, ja);

        for ( int ky = 0; ky < 8; ky++ ){
          for ( int kx = 0; kx < 8; kx++ ){
            for ( int iy = 0; iy < Model::NUM_VARS; iy++ ){
              for ( int ix = 0; ix < Model::NUM_VARS; ix++ ){
                jac(i, ky, kx, iy, ix) +=
                  N1x[ky] * ( N1x[kx] * ja(iy, 0, ix, 0) +
                              N2x[kx] * ja(iy, 0, ix, 1) +
                              N3x[kx] * ja(iy, 0, ix, 2) ) +
                  N2x[ky] * ( N1x[kx] * ja(iy, 1, ix, 0) +
                              N2x[kx] * ja(iy, 1, ix, 1) +
                              N3x[kx] * ja(iy, 1, ix, 2) ) +
                  N3x[ky] * ( N1x[kx] * ja(iy, 2, ix, 0) +
                              N2x[kx] * ja(iy, 2, ix, 1) +
                              N3x[kx] * ja(iy, 2, ix, 2) );
              }
            }
          }
        }
      }
    }
  }
};

#endif // ELASTICITY_3D_H
