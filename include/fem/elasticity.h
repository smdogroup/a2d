#ifndef A2D_ELASTICITY_H
#define A2D_ELASTICITY_H

#include "array.h"
#include "basis.h"
#include "constitutive.h"
#include "element.h"
#include "functional.h"
#include "utils/a2dprofiler.h"

namespace A2D {

/**
 * @brief Store all the type information about the PDE in one place
 *
 * @tparam spatial_dim 2 or 3
 * @tparam I index type
 * @tparam T data type
 */
template <index_t spatial_dim, typename I, typename T>
class ElasticityPDEInfo {
 public:
  static_assert(spatial_dim == 3 or spatial_dim == 2,
                "spatial_dim must be 2 or 3");
  static const index_t SPATIAL_DIM = spatial_dim;
  static const index_t vars_per_node = SPATIAL_DIM;
  static const index_t null_space_dim = SPATIAL_DIM * (SPATIAL_DIM + 1) / 2;
  static const index_t data_per_point = 2;
  static const index_t dvs_per_point = 1;

  // Set compile-time dimensions for arrays
  using BCsArray = A2D::MultiArrayNew<I* [2]>;
  using NodeArray = A2D::MultiArrayNew<T* [SPATIAL_DIM]>;
  using SolutionArray = A2D::MultiArrayNew<T* [vars_per_node]>;
  using DesignArray = A2D::MultiArrayNew<T* [dvs_per_point]>;
  using NullSpaceArray = A2D::MultiArrayNew<T* [vars_per_node][null_space_dim]>;

  // Jacobian matrix
  typedef A2D::BSRMat<I, T, vars_per_node, vars_per_node> SparseMat;

  // Sparse matrix multigrid type
  typedef A2D::BSRMatAmg<I, T, vars_per_node, null_space_dim> SparseAmg;

  static void compute_null_space(NodeArray& X, NullSpaceArray& B) {
    A2D::BLAS::zero(B);
    if (SPATIAL_DIM == 3) {
      for (I i = 0; i < B.extent(0); i++) {
        B(i, 0, 0) = 1.0;
        B(i, 1, 1) = 1.0;
        B(i, 2, 2) = 1.0;

        // Rotation about the x-axis
        B(i, 1, 3) = X(i, 2);
        B(i, 2, 3) = -X(i, 1);

        // Rotation about the y-axis
        B(i, 0, 4) = X(i, 2);
        B(i, 2, 4) = -X(i, 0);

        // Rotation about the z-axis
        B(i, 0, 5) = X(i, 1);
        B(i, 1, 5) = -X(i, 0);
      }
    } else {
      for (I i = 0; i < B.extent(0); i++) {
        B(i, 0, 0) = 1.0;
        B(i, 1, 1) = 1.0;

        // Rotation about the z-axis
        B(i, 0, 2) = X(i, 1);
        B(i, 1, 2) = -X(i, 0);
      }
    }
  }
};

template <typename I, typename T, class BasisOps>
class NonlinElasticityElement
    : public ElementBasis<I, T, ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I, T>,
                          BasisOps> {
 public:
  // Finite-element basis class
  static const index_t SPATIAL_DIM = BasisOps::SPATIAL_DIM;
  static const index_t NUM_VARS = SPATIAL_DIM;  // Number of variables per node
  static const index_t NUM_DATA = 2;  // Data points per quadrature point

  // Short cut for base class name
  typedef ElementBasis<I, T, ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I, T>,
                       BasisOps>
      base;

  NonlinElasticityElement(const index_t nelems)
      : ElementBasis<I, T, ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I, T>,
                     BasisOps>(nelems) {}

  template <typename IdxType>
  NonlinElasticityElement(const index_t nelems, const IdxType conn_[])
      : ElementBasis<I, T, ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I, T>,
                     BasisOps>(nelems, conn_) {}

  T energy() {
    auto data = this->get_quad_data();
    auto detJ = this->get_detJ();
    auto Jinv = this->get_Jinv();
    auto Uxi = this->get_quad_gradient();

    T engry = BasisOps::template integrate<T, NUM_VARS>(
        detJ, Jinv, Uxi,
        A2D_LAMBDA(index_t i, index_t j, T wdetJ,
                   A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> & Jinv0,
                   A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> & Uxi0)
            ->T {
              T mu(data(i, j, 0)), lambda(data(i, j, 1));
              A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> Ux;
              A2D::SymmMat<T, SPATIAL_DIM> E;
              T output;

              A2D::MatMatMult(Uxi0, Jinv0, Ux);
              A2D::MatGreenStrain(Ux, E);
              A2D::SymmIsotropicEnergy(mu, lambda, E, output);

              return wdetJ * output;
            });

    return engry;
  }

  void add_residual(typename ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I,
                                               T>::SolutionArray& res) {
    // Allocate the element residual
    typename base::ElemResArray elem_res("elem_res", this->nelems);

    // Retrieve the element data
    auto data = this->get_quad_data();
    auto detJ = this->get_detJ();
    auto Jinv = this->get_Jinv();
    auto Uxi = this->get_quad_gradient();

    BasisOps::template residuals<T, NUM_VARS>(
        detJ, Jinv, Uxi,
        A2D_LAMBDA(index_t i, index_t j, T wdetJ,
                   A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> & Jinv0,
                   A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> & Uxi0,
                   A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> & Uxib)
            ->void {
              T mu(data(i, j, 0)), lambda(data(i, j, 1));
              A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> Ux0, Uxb;
              A2D::SymmMat<T, SPATIAL_DIM> E0, Eb;

              A2D::ADMat<A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM>> Uxi(Uxi0, Uxib);
              A2D::ADMat<A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM>> Ux(Ux0, Uxb);
              A2D::ADMat<A2D::SymmMat<T, SPATIAL_DIM>> E(E0, Eb);
              A2D::ADScalar<T> output;

              auto mult = A2D::MatMatMult(Uxi, Jinv0, Ux);
              auto strain = A2D::MatGreenStrain(Ux, E);
              auto energy = A2D::SymmIsotropicEnergy(mu, lambda, E, output);

              output.bvalue = wdetJ;

              energy.reverse();
              strain.reverse();
              mult.reverse();
            },
        elem_res);

    VecElementGatherAdd(this->get_conn(), elem_res, res);
  }

  void add_jacobian(
      typename ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I, T>::SparseMat& J) {
    // Time this function
    Timer timer("NonlinElasticityElement::add_jacobian()");
    typename base::ElemJacArray elem_jac("elem_jac", this->nelems);

    // Retrieve the element data
    auto data = this->get_quad_data();
    auto detJ = this->get_detJ();
    auto Jinv = this->get_Jinv();
    auto Uxi = this->get_quad_gradient();

    BasisOps::template jacobians<T, NUM_VARS>(
        detJ, Jinv, Uxi,
        A2D_LAMBDA(index_t i, index_t j, T wdetJ,
                   A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> & Jinv0,
                   A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> & Uxi0,
                   A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> & Uxib,
                   A2D::SymmTensor<T, SPATIAL_DIM, SPATIAL_DIM> & jac)
            ->void {
              T mu(data(i, j, 0)), lambda(data(i, j, 1));
              A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> Ux0, Uxb;
              A2D::SymmMat<T, SPATIAL_DIM> E0, Eb;

              const int N = SPATIAL_DIM * SPATIAL_DIM;
              A2D::A2DMat<N, A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM>> Uxi(Uxi0,
                                                                        Uxib);
              A2D::A2DMat<N, A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM>> Ux(Ux0,
                                                                       Uxb);
              A2D::A2DMat<N, A2D::SymmMat<T, SPATIAL_DIM>> E(E0, Eb);
              A2D::A2DScalar<N, T> output;

              // Set up the seed values
              for (int k = 0; k < N; k++) {
                A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM>& Up = Uxi.pvalue(k);
                Up(k / SPATIAL_DIM, k % SPATIAL_DIM) = 1.0;
              }

              auto mult = A2D::MatMatMult(Uxi, Jinv0, Ux);
              auto strain = A2D::MatGreenStrain(Ux, E);
              auto energy = A2D::SymmIsotropicEnergy(mu, lambda, E, output);

              output.bvalue = wdetJ;

              energy.reverse();
              strain.reverse();
              mult.reverse();

              mult.hforward();
              strain.hforward();
              energy.hreverse();
              strain.hreverse();
              mult.hreverse();

              for (int k = 0; k < N; k++) {
                A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM>& Uxih = Uxi.hvalue(k);
                for (int i = 0; i < SPATIAL_DIM; i++) {
                  for (int j = 0; j < SPATIAL_DIM; j++) {
                    jac(i, j, k / (int)SPATIAL_DIM, k % (int)SPATIAL_DIM) =
                        Uxih(i, j);
                  }
                }
              }
            },
        elem_jac);

    A2D::BSRMatAddElementMatrices(this->get_conn(), elem_jac, J);
  }

  void add_adjoint_dfddata(typename ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I,
                                                      T>::SolutionArray& psi,
                           typename base::QuadDataArray& dfdx) {
    // Retrieve the element data
    auto conn = this->get_conn();
    auto detJ = this->get_detJ();
    auto Jinv = this->get_Jinv();
    auto Uxi = this->get_quad_gradient();
    auto data = this->get_quad_data();

    // Compute the element adjoint data
    typename base::ElemSolnArray pe("pe", this->nelems);
    typename base::QuadGradArray pxi("pxi", this->nelems);

    VecElementScatter(conn, psi, pe);
    BasisOps::template gradient<T, NUM_VARS>(pe, pxi);

    // Compute the product
    BasisOps::template adjoint_product<T, NUM_VARS>(
        detJ, Jinv, Uxi, pxi,
        A2D_LAMBDA(index_t i, index_t j, T wdetJ,
                   A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> & Jinv0,
                   A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> & Uxi0,
                   A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> & Pxi0)
            ->void {
              A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> Uxib, Ux0, Uxb;
              A2D::SymmMat<T, SPATIAL_DIM> E0, Eb;

              const int N = 1;
              A2D::A2DMat<N, A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM>> Uxi(Uxi0,
                                                                        Uxib);
              A2D::A2DMat<N, A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM>> Ux(Ux0,
                                                                       Uxb);
              A2D::A2DMat<N, A2D::SymmMat<T, SPATIAL_DIM>> E(E0, Eb);
              A2D::A2DScalar<N, T> output;
              A2D::A2DScalar<N, T> mu(data(i, j, 0)), lambda(data(i, j, 1));

              // Set the seed values
              A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM>& Psi = Uxi.pvalue(0);
              for (int k2 = 0; k2 < SPATIAL_DIM; k2++) {
                for (int k1 = 0; k1 < SPATIAL_DIM; k1++) {
                  Psi(k1, k2) = Pxi0(k1, k2);
                }
              }

              auto mult = A2D::MatMatMult(Uxi, Jinv0, Ux);
              auto strain = A2D::MatGreenStrain(Ux, E);
              auto energy = A2D::SymmIsotropicEnergy(mu, lambda, E, output);

              output.bvalue = wdetJ;

              energy.reverse();
              strain.reverse();
              mult.reverse();

              mult.hforward();
              strain.hforward();
              energy.hreverse();

              dfdx(i, j, 0) = mu.hvalue[0];
              dfdx(i, j, 1) = lambda.hvalue[0];
            });
  }
};

template <typename I, typename T, class BasisOps>
class LinElasticityElement
    : public ElementBasis<I, T, ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I, T>,
                          BasisOps> {
 public:
  // Finite-element basis class
  static const index_t SPATIAL_DIM = BasisOps::SPATIAL_DIM;
  static const index_t NUM_VARS = SPATIAL_DIM;  // Number of variables per node
  static const index_t NUM_DATA = 2;  // Data points per quadrature point

  // Short cut for base class name
  typedef ElementBasis<I, T, ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I, T>,
                       BasisOps>
      base;

  LinElasticityElement(const index_t nelems)
      : ElementBasis<I, T, ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I, T>,
                     BasisOps>(nelems) {}

  template <typename IdxType>
  LinElasticityElement(const index_t nelems, const IdxType conn_[])
      : ElementBasis<I, T, ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I, T>,
                     BasisOps>(nelems, conn_) {}

  T energy() {
    auto data = this->get_quad_data();
    auto detJ = this->get_detJ();
    auto Jinv = this->get_Jinv();
    auto Uxi = this->get_quad_gradient();

    T elem_energy = BasisOps::template integrate<T, NUM_VARS>(
        detJ, Jinv, Uxi,
        A2D_LAMBDA(index_t i, index_t j, T wdetJ,
                   A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> & Jinv0,
                   A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> Uxi0)
            ->T {
              T mu(data(i, j, 0)), lambda(data(i, j, 1));
              A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> Ux;
              A2D::SymmMat<T, SPATIAL_DIM> E;
              T output;

              A2D::MatMatMult(Uxi0, Jinv0, Ux);
              A2D::MatLinearGreenStrain(Ux, E);
              A2D::SymmIsotropicEnergy(mu, lambda, E, output);

              return wdetJ * output;
            });

    return elem_energy;
  }

  void add_residual(typename ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I,
                                               T>::SolutionArray& res) {
    // Allocate the element residual
    typename base::ElemResArray elem_res("elem_res", this->nelems);

    // Retrieve the element data
    auto data = this->get_quad_data();
    auto detJ = this->get_detJ();
    auto Jinv = this->get_Jinv();
    auto Uxi = this->get_quad_gradient();

    BasisOps::template residuals<T, NUM_VARS>(
        detJ, Jinv, Uxi,
        A2D_LAMBDA(index_t i, index_t j, T wdetJ,
                   A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> & Jinv0,
                   A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> & Uxi0,
                   A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> & Uxib)
            ->void {
              T mu(data(i, j, 0)), lambda(data(i, j, 1));
              A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> Ux0, Uxb;
              A2D::SymmMat<T, SPATIAL_DIM> E0, Eb;

              A2D::ADMat<A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM>> Uxi(Uxi0, Uxib);
              A2D::ADMat<A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM>> Ux(Ux0, Uxb);
              A2D::ADMat<A2D::SymmMat<T, SPATIAL_DIM>> E(E0, Eb);
              A2D::ADScalar<T> output;

              auto mult = A2D::MatMatMult(Uxi, Jinv0, Ux);
              auto strain = A2D::MatLinearGreenStrain(Ux, E);
              auto energy = A2D::SymmIsotropicEnergy(mu, lambda, E, output);

              output.bvalue = wdetJ;

              energy.reverse();
              strain.reverse();
              mult.reverse();
            },
        elem_res);

    VecElementGatherAdd(this->get_conn(), elem_res, res);
  }

  void add_jacobian(
      typename ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I, T>::SparseMat& J) {
    Timer timer("LinElasticityElement::add_jacobian()");
    typename base::ElemJacArray elem_jac("elem_jac", this->nelems);

    // Retrieve the element data
    auto data = this->get_quad_data();
    auto detJ = this->get_detJ();
    auto Jinv = this->get_Jinv();
    auto Uxi = this->get_quad_gradient();

    BasisOps::template jacobians<T, NUM_VARS>(
        detJ, Jinv, Uxi,
        A2D_LAMBDA(index_t i, index_t j, T wdetJ,
                   A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> & Jinv0,
                   A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> & Uxi0,
                   A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> & Uxib,
                   A2D::SymmTensor<T, SPATIAL_DIM, SPATIAL_DIM> & jac)
            ->void {
              T mu(data(i, j, 0)), lambda(data(i, j, 1));
              A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> Ux0, Uxb;
              A2D::SymmMat<T, SPATIAL_DIM> E0, Eb;

              const int N = SPATIAL_DIM * SPATIAL_DIM;
              A2D::A2DMat<N, A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM>> Uxi(Uxi0,
                                                                        Uxib);
              A2D::A2DMat<N, A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM>> Ux(Ux0,
                                                                       Uxb);
              A2D::A2DMat<N, A2D::SymmMat<T, SPATIAL_DIM>> E(E0, Eb);
              A2D::A2DScalar<N, T> output;

              // Set up the seed values
              for (int k = 0; k < N; k++) {
                A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM>& Up = Uxi.pvalue(k);
                Up(k / SPATIAL_DIM, k % SPATIAL_DIM) = 1.0;
              }

              auto mult = A2D::MatMatMult(Uxi, Jinv0, Ux);
              auto strain = A2D::MatLinearGreenStrain(Ux, E);
              auto energy = A2D::SymmIsotropicEnergy(mu, lambda, E, output);

              output.bvalue = wdetJ;

              energy.reverse();
              strain.reverse();
              mult.reverse();

              mult.hforward();
              strain.hforward();
              energy.hreverse();
              strain.hreverse();
              mult.hreverse();

              for (int k = 0; k < N; k++) {
                A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM>& Uxih = Uxi.hvalue(k);
                for (int i = 0; i < SPATIAL_DIM; i++) {
                  for (int j = 0; j < SPATIAL_DIM; j++) {
                    jac(i, j, k / (int)SPATIAL_DIM, k % (int)SPATIAL_DIM) =
                        Uxih(i, j);
                  }
                }
              }
            },
        elem_jac);

    A2D::BSRMatAddElementMatrices(this->get_conn(), elem_jac, J);
  }

  void add_adjoint_dfddata(typename ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I,
                                                      T>::SolutionArray& psi,
                           typename base::QuadDataArray& dfdx) {
    // Retrieve the element data
    auto conn = this->get_conn();
    auto detJ = this->get_detJ();
    auto Jinv = this->get_Jinv();
    auto Uxi = this->get_quad_gradient();
    auto data = this->get_quad_data();

    // Compute the element adjoint data
    typename base::ElemSolnArray pe("pe", this->nelems);
    typename base::QuadGradArray pxi("pxi", this->nelems);

    VecElementScatter(conn, psi, pe);
    BasisOps::template gradient<T, NUM_VARS>(pe, pxi);

    // Compute the product
    BasisOps::template adjoint_product<T, NUM_VARS>(
        detJ, Jinv, Uxi, pxi,
        A2D_LAMBDA(index_t i, index_t j, T wdetJ,
                   A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> & Jinv0,
                   A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> & Uxi0,
                   A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> & Pxi0)
            ->void {
              A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> Uxib, Ux0, Uxb;
              A2D::SymmMat<T, SPATIAL_DIM> E0, Eb;

              const int N = 1;
              A2D::A2DMat<N, A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM>> Uxi(Uxi0,
                                                                        Uxib);
              A2D::A2DMat<N, A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM>> Ux(Ux0,
                                                                       Uxb);
              A2D::A2DMat<N, A2D::SymmMat<T, SPATIAL_DIM>> E(E0, Eb);
              A2D::A2DScalar<N, T> output;
              A2D::A2DScalar<N, T> mu(data(i, j, 0)), lambda(data(i, j, 1));

              // Set the seed values
              A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM>& Psi = Uxi.pvalue(0);
              for (int k2 = 0; k2 < SPATIAL_DIM; k2++) {
                for (int k1 = 0; k1 < SPATIAL_DIM; k1++) {
                  Psi(k1, k2) = Pxi0(k1, k2);
                }
              }

              auto mult = A2D::MatMatMult(Uxi, Jinv0, Ux);
              auto strain = A2D::MatLinearGreenStrain(Ux, E);
              auto energy = A2D::SymmIsotropicEnergy(mu, lambda, E, output);

              output.bvalue = wdetJ;

              energy.reverse();
              strain.reverse();
              mult.reverse();

              mult.hforward();
              strain.hforward();
              energy.hreverse();

              dfdx(i, j, 0) = mu.hvalue[0];
              dfdx(i, j, 1) = lambda.hvalue[0];
            });
  }
};

template <typename I, typename T, class BasisOps>
class TopoIsoConstitutive final
    : public ConstitutiveBase<I, T,
                              ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I, T>> {
 public:
  static const index_t SPATIAL_DIM = BasisOps::SPATIAL_DIM;
  static const index_t dvs_per_point =
      ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I, T>::dvs_per_point;
  static const index_t nodes_per_elem = BasisOps::NUM_NODES;
  static const index_t quad_pts_per_elem = BasisOps::quadrature::NUM_QUAD_PTS;

  using ElemDesignArray =
      A2D::MultiArrayNew<T* [nodes_per_elem][dvs_per_point]>;
  using QuadDesignArray =
      A2D::MultiArrayNew<T* [quad_pts_per_elem][dvs_per_point]>;

  TopoIsoConstitutive(
      std::shared_ptr<ElementBasis<
          I, T, ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I, T>, BasisOps>>
          element,
      T q, T E, T nu, T density, T design_stress, T beta = 20.0,
      T xoffset = 0.5)
      : q(q),
        xoffset(xoffset),
        beta(beta),
        E(E),
        nu(nu),
        density(density),
        design_stress(design_stress),
        element(element) {
    Timer t("TopoIsoConstitutive::TopoIsoConstitutive");
    xe = ElemDesignArray("xe", element->nelems);
    xq = QuadDesignArray("xq", element->nelems);
    mu = 0.5 * E / (1.0 + nu);
    lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
  }

  // Penalization value
  const T q;

  // Heaviside filter approximation
  const T xoffset;
  const T beta;

  // Constitutive data
  const T E;
  const T nu;
  const T density;
  const T design_stress;

  // Get the Lame parameters
  void get_lame_parameters(T& mu_, T& lambda_) {
    mu_ = mu;
    lambda_ = lambda;
  }

  /*
    Set the design variables values into the element object
  */
  void set_design_vars(
      typename ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I, T>::DesignArray& x) {
    // Set the design variable values
    // x -> xe -> xq
    auto conn = element->get_conn();
    VecElementScatter(conn, x, xe);
    BasisOps::template interp<dvs_per_point>(xe, xq);

    auto data = element->get_quad_data();
    for (I i = 0; i < data.extent(0); i++) {
      for (I j = 0; j < data.extent(1); j++) {
        T rho_exp = std::exp(-beta * (xq(i, j, 0) - xoffset));
        T rho = 1.0 / (1.0 + rho_exp);
        T penalty = rho / (1.0 + q * (1.0 - rho));

        data(i, j, 0) = mu * penalty;
        data(i, j, 1) = lambda * penalty;
      }
    }
  }

  /*
    Compute the derivative of the adjoint-residual product data w.r.t. x
  */
  void add_adjoint_dfdx(typename ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I,
                                                   T>::SolutionArray& psi,
                        typename ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I,
                                                   T>::DesignArray& dfdx) {
    typename ElementBasis<I, T, ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I, T>,
                          BasisOps>::QuadDataArray dfddata("dfddata",
                                                           element->nelems);

    // Compute the product of the adjoint with the derivatives of
    // the residuals w.r.t. the element data
    element->add_adjoint_dfddata(psi, dfddata);

    // Set the result into the dfdxq array
    QuadDesignArray dfdxq("dfdxq", element->nelems);
    for (I i = 0; i < dfddata.extent(0); i++) {
      for (I j = 0; j < dfddata.extent(1); j++) {
        T rho_exp = std::exp(-beta * (xq(i, j, 0) - xoffset));
        T rho = 1.0 / (1.0 + rho_exp);
        T denom = (1.0 + q * (1.0 - rho));
        T dpenalty = (q + 1.0) / (denom * denom);
        dpenalty *= beta * rho_exp * rho * rho;

        dfdxq(i, j, 0) =
            dpenalty * (mu * dfddata(i, j, 0) + lambda * dfddata(i, j, 1));
      }
    }

    ElemDesignArray dfdxe("dfdxe", element->nelems);
    BasisOps::template interpReverseAdd<dvs_per_point>(dfdxq, dfdxe);

    auto conn = element->get_conn();
    VecElementGatherAdd(conn, dfdxe, dfdx);
  }

  std::shared_ptr<ElementBasis<
      I, T, ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I, T>, BasisOps>>
  get_element() {
    return element;
  }

  ElemDesignArray& get_elem_design() { return xe; }
  QuadDesignArray& get_quad_design() { return xq; }

 private:
  // Reference to the element class
  std::shared_ptr<ElementBasis<
      I, T, ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I, T>, BasisOps>>
      element;

  // Design variable views
  ElemDesignArray xe;
  QuadDesignArray xq;

  // Parameter value
  T mu, lambda;
};

/*
  Evaluate the volume of the structure, given the constitutive class
*/
template <typename I, typename T, class BasisOps>
class TopoVolume
    : public ElementFunctional<I, T,
                               ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I, T>> {
 public:
  static const index_t SPATIAL_DIM = BasisOps::SPATIAL_DIM;
  static const index_t vars_per_node =
      ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I, T>::vars_per_node;
  static const index_t dvs_per_point =
      ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I, T>::dvs_per_point;

  TopoVolume(std::shared_ptr<TopoIsoConstitutive<I, T, BasisOps>> con)
      : con(con) {}

  T eval_functional() {
    auto element = con->get_element();
    auto detJ = element->get_detJ();
    auto xq = con->get_quad_design();

    T integral = BasisOps::template integrate<T>(
        detJ, A2D_LAMBDA(index_t i, index_t j, T wdetJ)->T {
          return xq(i, j, 0) * wdetJ;
        });

    return integral;
  }
  void add_dfdx(typename ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I,
                                           T>::DesignArray& dfdx) {
    auto element = con->get_element();
    auto detJ = element->get_detJ();
    auto xq = con->get_quad_design();

    typename TopoIsoConstitutive<I, T, BasisOps>::QuadDesignArray dfdxq(
        "dfdxq", element->nelems);

    T integral = BasisOps::template integrate<T>(
        detJ, A2D_LAMBDA(index_t i, index_t j, T wdetJ)->T {
          dfdxq(i, j, 0) = wdetJ;
          return xq(i, j, 0) * wdetJ;
        });

    typename TopoIsoConstitutive<I, T, BasisOps>::ElemDesignArray dfdxe(
        "dfdxe", element->nelems);
    BasisOps::template interpReverseAdd<dvs_per_point>(dfdxq, dfdxe);

    auto conn = element->get_conn();
    VecElementGatherAdd(conn, dfdxe, dfdx);
  }

 private:
  std::shared_ptr<TopoIsoConstitutive<I, T, BasisOps>> con;
};

/*
  Evalute the KS functional of the stress, given the constitutive class
*/
template <typename I, typename T, class BasisOps>
class TopoVonMisesAggregation
    : public ElementFunctional<I, T,
                               ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I, T>> {
 public:
  static const index_t SPATIAL_DIM = BasisOps::SPATIAL_DIM;
  static const int vars_per_node =
      ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I, T>::vars_per_node;
  static const index_t dvs_per_point =
      ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I, T>::dvs_per_point;

  TopoVonMisesAggregation(
      std::shared_ptr<TopoIsoConstitutive<I, T, BasisOps>> con,
      T weight = 100.0)
      : weight(weight), con(con) {
    offset = 0.0;
    integral = 1.0;
  }

  // The KS aggregation weight
  const T weight;

  /*
    Reset the maximum value
  */
  void compute_offset() {
    auto element = con->get_element();
    auto data = element->get_quad_data();
    auto detJ = element->get_detJ();
    auto Jinv = element->get_Jinv();
    auto Uxi = element->get_quad_gradient();
    auto xq = con->get_quad_design();

    T ys = con->design_stress;
    T mu, lambda;
    con->get_lame_parameters(mu, lambda);
    T qval = con->q;

    // Compute the maximum value over all quadrature points
    offset = BasisOps::template maximum<T, vars_per_node>(
        detJ, Jinv, Uxi,
        A2D_LAMBDA(index_t i, index_t j, T wdetJ,
                   A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> & Jinv0,
                   A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> & Uxi0)
            ->T {
              A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> Ux;
              A2D::SymmMat<T, SPATIAL_DIM> E, S;
              T vm, trS, trSS;

              A2D::MatMatMult(Uxi0, Jinv0, Ux);
              A2D::MatGreenStrain(Ux, E);
              A2D::SymmIsotropicConstitutive(mu, lambda, E, S);
              A2D::SymmTrace(S, trS);
              A2D::SymmSymmMultTrace(S, S, trSS);

              // Compute the penalty = (q + 1) * x/(q * x + 1)
              T penalty =
                  (qval + 1.0) * xq(i, j, 0) / (qval * xq(i, j, 0) + 1.0);

              // von Mises = 1.5 * tr(S * S) - 0.5 * tr(S)**2;
              vm = penalty * (1.5 * trSS - 0.5 * trS * trS) / ys;

              return vm;
            });
  }

  T eval_functional() {
    auto element = con->get_element();
    auto data = element->get_quad_data();
    auto detJ = element->get_detJ();
    auto Jinv = element->get_Jinv();
    auto Uxi = element->get_quad_gradient();
    auto xq = con->get_quad_design();

    T ys = con->design_stress;
    T mu, lambda;
    con->get_lame_parameters(mu, lambda);
    T qval = con->q;
    T off = offset;
    T wgt = weight;

    // Compute the maximum value over all quadrature points
    integral = BasisOps::template maximum<T, vars_per_node>(
        detJ, Jinv, Uxi,
        A2D_LAMBDA(index_t i, index_t j, T wdetJ,
                   A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> & Jinv0,
                   A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> & Uxi0)
            ->T {
              A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> Ux;
              A2D::SymmMat<T, SPATIAL_DIM> E, S;
              T vm, trS, trSS;

              A2D::MatMatMult(Uxi0, Jinv0, Ux);
              A2D::MatGreenStrain(Ux, E);
              A2D::SymmIsotropicConstitutive(mu, lambda, E, S);
              A2D::SymmTrace(S, trS);
              A2D::SymmSymmMultTrace(S, S, trSS);

              // Compute the penalty = (q + 1) * x/(q * x + 1)
              T penalty =
                  (qval + 1.0) * xq(i, j, 0) / (qval * xq(i, j, 0) + 1.0);

              // von Mises = 1.5 * tr(S * S) - 0.5 * tr(S)**2;
              vm = penalty * (1.5 * trSS - 0.5 * trS * trS) / ys;

              return wdetJ * std::exp(wgt * (vm - off));
            });

    return offset + std::log(integral) / weight;
  }

  void add_dfdu(typename ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I,
                                           T>::SolutionArray& dfdu) {
    auto element = con->get_element();
    typename ElementBasis<I, T, ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I, T>,
                          BasisOps>::ElemResArray elem_dfdu("elem_dfdu",
                                                            element->nelems);

    auto data = element->get_quad_data();
    auto detJ = element->get_detJ();
    auto Jinv = element->get_Jinv();
    auto Uxi = element->get_quad_gradient();
    auto xq = con->get_quad_design();

    T ys = con->design_stress;
    T mu, lambda;
    con->get_lame_parameters(mu, lambda);
    T qval = con->q;
    T off = offset;
    T wgt = weight;
    T intgrl = integral;

    BasisOps::template residuals<T, vars_per_node>(
        detJ, Jinv, Uxi,
        A2D_LAMBDA(index_t i, index_t j, T wdetJ,
                   A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> & Jinv0,
                   A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> & Uxi0,
                   A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> & Uxib)
            ->void {
              A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> Ux0, Uxb;
              A2D::SymmMat<T, SPATIAL_DIM> E0, Eb;
              A2D::SymmMat<T, SPATIAL_DIM> S0, Sb;

              A2D::ADMat<A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM>> Uxi(Uxi0, Uxib);
              A2D::ADMat<A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM>> Ux(Ux0, Uxb);
              A2D::ADMat<A2D::SymmMat<T, SPATIAL_DIM>> E(E0, Eb);
              A2D::ADMat<A2D::SymmMat<T, SPATIAL_DIM>> S(S0, Sb);
              A2D::ADScalar<T> trS, trSS;

              auto mult = A2D::MatMatMult(Uxi, Jinv0, Ux);
              auto strain = A2D::MatGreenStrain(Ux, E);
              auto cons = A2D::SymmIsotropicConstitutive(mu, lambda, E, S);
              auto trace1 = A2D::SymmTrace(S, trS);
              auto trace2 = A2D::SymmSymmMultTrace(S, S, trSS);

              // Compute the penalty = (q + 1) * x/(q * x + 1)
              T penalty =
                  (qval + 1.0) * xq(i, j, 0) / (qval * xq(i, j, 0) + 1.0);

              // von Mises = 1.5 * tr(S * S) - 0.5 * tr(S)**2;
              T vm = penalty *
                     (1.5 * trSS.value - 0.5 * trS.value * trS.value) / ys;

              T scale =
                  wdetJ * penalty * std::exp(wgt * (vm - off)) / (ys * intgrl);

              trSS.bvalue = 1.5 * scale;
              trS.bvalue = -trS.value * scale;

              trace2.reverse();
              trace1.reverse();
              cons.reverse();
              strain.reverse();
              mult.reverse();
            },
        elem_dfdu);

    VecElementGatherAdd(element->get_conn(), elem_dfdu, dfdu);
  }

  void add_dfdx(typename ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I,
                                           T>::DesignArray& dfdx) {
    auto element = con->get_element();
    auto data = element->get_quad_data();
    auto detJ = element->get_detJ();
    auto Jinv = element->get_Jinv();
    auto Uxi = element->get_quad_gradient();
    auto xq = con->get_quad_design();

    T ys = con->design_stress;
    T mu, lambda;
    con->get_lame_parameters(mu, lambda);
    T qval = con->q;
    T off = offset;
    T wgt = weight;
    T intgrl = integral;

    typename TopoIsoConstitutive<I, T, BasisOps>::QuadDesignArray dfdxq(
        "dfdxq", element->nelems);

    BasisOps::template maximum<T, vars_per_node>(
        detJ, Jinv, Uxi,
        A2D_LAMBDA(index_t i, index_t j, T wdetJ,
                   A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> & Jinv0,
                   A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> & Uxi0)
            ->T {
              A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> Ux;
              A2D::SymmMat<T, SPATIAL_DIM> E, S;
              T trS, trSS;

              A2D::MatMatMult(Uxi0, Jinv0, Ux);
              A2D::MatGreenStrain(Ux, E);
              A2D::SymmIsotropicConstitutive(mu, lambda, E, S);
              A2D::SymmTrace(S, trS);
              A2D::SymmSymmMultTrace(S, S, trSS);

              // Compute the penalty = (q + 1) * x/(q * x + 1)
              T denom = (qval * xq(i, j, 0) + 1.0) * (qval * xq(i, j, 0) + 1.0);
              T dpenalty = (qval + 1.0) / denom;

              // von Mises = 1.5 * tr(S * S) - 0.5 * tr(S)**2;
              T vm = (1.5 * trSS - 0.5 * trS * trS) / ys;

              T scale = wdetJ * vm * std::exp(wgt * (vm - off)) / intgrl;

              dfdxq(i, j, 0) += scale * dpenalty;

              return wdetJ * std::exp(wgt * (vm - off));
            });

    typename TopoIsoConstitutive<I, T, BasisOps>::ElemDesignArray dfdxe(
        "dfdxe", element->nelems);
    BasisOps::template interpReverseAdd<dvs_per_point>(dfdxq, dfdxe);

    auto conn = element->get_conn();
    VecElementGatherAdd(conn, dfdxe, dfdx);
  }

 private:
  std::shared_ptr<TopoIsoConstitutive<I, T, BasisOps>> con;
  T offset;    // Offset value for computing the KS function value
  T integral;  // Integral: int_{Omega} e^{weight*(vm - offset)} dOmega
};

}  // namespace A2D

#endif  // A2D_ELASTICITY_H
