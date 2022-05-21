#ifndef ELASTICITY_3D_H
#define ELASTICITY_3D_H

#include "a2dtmp.h"

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

  template<typename T, class IdxType, class QuadPointData, class MatType>
  static T compute_energy( IdxType i, IdxType j, QuadPointData& data, MatType& Ux ){
    typedef SymmMat<T, 3> SymmMat3x3;

    T mu(data(i, j, 0)), lambda(data(i, j, 1));
    T output, outputb(1.0);

    Mat3x3GreenStrain<MatType, SymmMat3x3> strain(Ux, E);
    Symm3x3IsotropicConstitutive<T, SymmMat3x3, SymmMat3x3> stress(mu, lambda, E, S);
    Symm3x3SymmMultTrace<SymmMat3x3, SymmMat3x3, T> trace(E, S, output);

    return outpu;
  }

  template<typename T, class IdxType, class QuadPointData, class MatType>
  static T compute_residual( IdxType i, IdxType j, QuadPointData& data,
                             MatType& Uxi, MatType& Uxd ){
    typedef SymmMat<T, 3> SymmMat3x3;

    T mu(data(i, j, 0)), lambda(data(i, j, 1));
    T output, outputb(1.0);

    Mat3x3GreenStrain<MatType, SymmMat3x3> strain(Ux, E);
    Symm3x3IsotropicConstitutive<T, SymmMat3x3, SymmMat3x3> stress(mu, lambda, E, S);
    Symm3x3SymmMultTrace<SymmMat3x3, SymmMat3x3, T> trace(E, S, output);

    trace.reverse(outputb, Eb, Sb);
    stress.reverse(Sb, Eb);
    strain.reverse(Eb, Uxd);

    return output;
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
           class QuadPointJacobianArray,
           class QuadPointDetJArray>
  static void compute_jtrans( ElementNodeArray& X,
                              QuadPointJacobianArray& Jinv,
                              QuadPointDetJArray& detJ ){
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
        A2D::Mat<T, 3, 3> JinvObj;
        detJ(i, j) = A2D::Mat3x3InverseCore<T>(J, JinvObj);

        // Copy values of the determinant of the inverse of the Jacobian
        for ( int ii = 0; ii < 3; ii++ ){
          for ( int jj = 0; jj < 3; jj++ ){
            Jinv(i, j, ii, jj) = JinvObj(ii, jj);
          }
        }
      }
    }
  }

  template<typename T,
           const int num_vars,
           class ElementSolutionArray,
           class QuadPointJacobianArray,
           class QuadPointGradientArray>
  static void gradient( ElementSolutionArray& U,
                        QuadPointJacobianArray& Jinv,
                        QuadPointGradientArray& Ux ){
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
        A2D::Mat<T, 3, 3> JinvObj(Jinv, i, j);

        for ( int ii = 0; ii < num_vars; ii++ ){
          A2D::Vec<T, 3> UxiObj;
          UxiObj(0) =
            n3[0] * (n2[0] * (n1x[0] * U(i, 0, ii) + n1x[1] * U(i, 1, ii)) +
                     n2[1] * (n1x[0] * U(i, 2, ii) + n1x[1] * U(i, 3, ii))) +
            n3[1] * (n2[0] * (n1x[0] * U(i, 4, ii) + n1x[1] * U(i, 5, ii)) +
                     n2[1] * (n1x[0] * U(i, 6, ii) + n1x[1] * U(i, 7, ii)));

          UxiObj(1) =
            n3[0] * (n2x[0] * (n1[0] * U(i, 0, ii) + n1[1] * U(i, 1, ii)) +
                     n2x[1] * (n1[0] * U(i, 2, ii) + n1[1] * U(i, 3, ii))) +
            n3[1] * (n2x[0] * (n1[0] * U(i, 4, ii) + n1[1] * U(i, 5, ii)) +
                     n2x[1] * (n1[0] * U(i, 6, ii) + n1[1] * U(i, 7, ii)));

          UxiObj(2) =
            n3x[0] * (n2[0] * (n1[0] * U(i, 0, ii) + n1[1] * U(i, 1, ii)) +
                      n2[1] * (n1[0] * U(i, 2, ii) + n1[1] * U(i, 3, ii))) +
            n3x[1] * (n2[0] * (n1[0] * U(i, 4, ii) + n1[1] * U(i, 5, ii)) +
                      n2[1] * (n1[0] * U(i, 6, ii) + n1[1] * U(i, 7, ii)));

          A2D::Vec<T, 3> UxObj;
          A2D::MatTrans3x3VecMult(JinvObj, UxiObj, UxObj);

          for ( int jj = 0; jj < 3; jj++ ){
            Ux(i, j, ii, jj) = UxObj(jj);
          }
        }
      }
    }
  }

  template<typename T,
           class Model,
           class QuadPointGradientArray,
           class QuadPointModelDataArray,
           class QuadPointJacobianArray,
           class QuadPointDetJArray>
  static void energy( QuadPointGradientArray& Ux,
                      QuadPointModelDataArray& Edata,
                      QuadPointJacobianArray& Jinv,
                      QuadPointDetJArray& detJ,
                      T& energy ){
    for ( std::size_t j = 0; j < Quadrature::getNumQuadPoints(); j++ ){
      double weight = Quadrature::getQuadWeight(j);

      for ( std::size_t i = 0; i < Ux.extent(0); i++ ){
        // Extract U,x to UxObj
        Mat<T, Model::NUM_VARS, 3> UxObj(Ux, i, j);
        for ( int ii = 0; ii < Model::NUM_VARS; ii++ ){
          for ( int jj = 0; jj < 3; jj++ ){
            UxObj(ii, jj) = Ux(i, j, ii, jj);
          }
        }

        energy += weight * detJ(i, j) * Model::getEnergy(i, j, Edata, UxObj);
      }
    }
  }

  template<typename T,
           class Model,
           class QuadPointGradientArray,
           class QuadPointModelDataArray,
           class QuadPointJacobianArray,
           class QuadPointDetJArray,
           class ElementResidualArray>
  static void residuals( QuadPointGradientArray& Ux,
                         QuadPointModelDataArray& Edata,
                         QuadPointJacobianArray& Jinv,
                         QuadPointDetJArray& detJ,
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

      for ( std::size_t i = 0; i < Ux.extent(0); i++ ){
        // Extract U,x
        Mat<T, Model::NUM_VARS, 3> UxObj(Ux, i, j);
        Mat<T, Model::NUM_VARS, 3> dUxObj;

        // Extract the energy and residual
        dUxObj.zero();
        Model::getEnergyAndResidual(i, j, Edata, UxObj, dUxObj);

        for ( int ii = 0; ii < Model::NUM_VARS; ii++ ){
          A2D::Vec3<T> UxiObj;
          UxiObj(0) =
            n3[0] * (n2[0] * (n1x[0] * U(i, 0, ii) + n1x[1] * U(i, 1, ii)) +
                     n2[1] * (n1x[0] * U(i, 2, ii) + n1x[1] * U(i, 3, ii))) +
            n3[1] * (n2[0] * (n1x[0] * U(i, 4, ii) + n1x[1] * U(i, 5, ii)) +
                     n2[1] * (n1x[0] * U(i, 6, ii) + n1x[1] * U(i, 7, ii)));

          UxiObj(1) =
            n3[0] * (n2x[0] * (n1[0] * U(i, 0, ii) + n1[1] * U(i, 1, ii)) +
                     n2x[1] * (n1[0] * U(i, 2, ii) + n1[1] * U(i, 3, ii))) +
            n3[1] * (n2x[0] * (n1[0] * U(i, 4, ii) + n1[1] * U(i, 5, ii)) +
                     n2x[1] * (n1[0] * U(i, 6, ii) + n1[1] * U(i, 7, ii)));

          UxiObj(2) =
            n3x[0] * (n2[0] * (n1[0] * U(i, 0, ii) + n1[1] * U(i, 1, ii)) +
                      n2[1] * (n1[0] * U(i, 2, ii) + n1[1] * U(i, 3, ii))) +
            n3x[1] * (n2[0] * (n1[0] * U(i, 4, ii) + n1[1] * U(i, 5, ii)) +
                      n2[1] * (n1[0] * U(i, 6, ii) + n1[1] * U(i, 7, ii)));

          A2D::Vec3<T> UxObj;
          A2D::MatTrans3x3VecMult<T> multObj(JinvObj, UxiObj, UxObj);

          for ( int jj = 0; jj < 3; jj++ ){
            Ux(i, j, ii, jj) = UxObj(jj);
          }
        }
      }
    }
  }
};

#endif // ELASTICITY_3D_H
