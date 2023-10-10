#ifndef A2D_FE_MAPPING_H
#define A2D_FE_MAPPING_H

#include "a2dcore.h"
#include "multiphysics/febasis.h"
#include "multiphysics/fespace.h"

namespace A2D {

template <
    class Geometry, class Space, typename T,
    std::enable_if_t<get_diff_type<Geometry>::diff_type == ADiffType::PASSIVE,
                     bool> = true,
    std::enable_if_t<get_diff_type<Space>::diff_type == ADiffType::PASSIVE,
                     bool> = true,
    std::enable_if_t<get_diff_type<T>::diff_type == ADiffType::PASSIVE, bool> =
        true>
KOKKOS_FUNCTION void RefElementTransform(const Geometry& geo, const Space& in,
                                         T& detJ, Space& out) {
  static_assert(Geometry::dim == Space::dim,
                "Spatial and finite-element space dimensions must agree");

  const Mat<T, Geometry::dim, Geometry::dim>& J = get_grad<0>(geo);
  Mat<T, Geometry::dim, Geometry::dim> Jinv;
  MatInv(J, Jinv);
  MatDet(J, detJ);
  in.transform(detJ, J, Jinv, out);
}

template <class Geometry, class Space, typename T>
class RefElementTransformConstGeoExpr {
 public:
  // In this case, the geometry must be passive
  static_assert(get_diff_type<Geometry>::diff_type == ADiffType::PASSIVE,
                "Must use passive geometry type");

  // Set the type of matrix to use - here use a constant matrix
  using JinvType = Mat<T, Geometry::dim, Geometry::dim>;

  KOKKOS_FUNCTION RefElementTransformConstGeoExpr(const Geometry& geo,
                                                  Space& in, T& detJ,
                                                  Space& out)
      : geo(geo), in(in), detJ(detJ), out(out) {}

  KOKKOS_FUNCTION void eval() {
    const Mat<T, Geometry::dim, Geometry::dim>& J = get_grad<0>(geo);
    MatInv(J, Jinv);
    MatDet(J, detJ);
    in.value().transform(detJ, J, Jinv, out.value());
  }

  KOKKOS_FUNCTION void bzero() { out.bzero(); }

  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;
    const Mat<T, Geometry::dim, Geometry::dim>& J = get_grad<0>(geo);
    GetSeed<seed>::get_obj(in).transform(detJ, J, Jinv,
                                         GetSeed<seed>::get_obj(out));
  }

  KOKKOS_FUNCTION void reverse() {
    const Mat<T, Geometry::dim, Geometry::dim>& J = get_grad<0>(geo);
    out.bvalue().btransform(detJ, J, Jinv, in.bvalue());
  }

  KOKKOS_FUNCTION void hzero() { out.hzero(); }

  KOKKOS_FUNCTION void hreverse() {
    const Mat<T, Geometry::dim, Geometry::dim>& J = get_grad<0>(geo);
    out.hvalue().btransform(detJ, J, Jinv, in.hvalue());
  }

 private:
  // Internal data
  Mat<T, Geometry::dim, Geometry::dim> Jinv;

  // Input/output data
  const Geometry& geo;
  Space& in;
  T detJ;
  Space& out;
};

template <
    class Geometry, class Space, typename T,
    std::enable_if_t<get_diff_type<Geometry>::diff_type == ADiffType::PASSIVE,
                     bool> = true,
    std::enable_if_t<get_diff_type<T>::diff_type == ADiffType::PASSIVE, bool> =
        true>
KOKKOS_FUNCTION auto RefElementTransform(const Geometry& geo, ADObj<Space>& in,
                                         T& detJ, ADObj<Space>& out) {
  return RefElementTransformConstGeoExpr<Geometry, ADObj<Space>, T>(geo, in,
                                                                    detJ, out);
}

template <
    class Geometry, class Space, typename T,
    std::enable_if_t<get_diff_type<Geometry>::diff_type == ADiffType::PASSIVE,
                     bool> = true,
    std::enable_if_t<get_diff_type<T>::diff_type == ADiffType::PASSIVE, bool> =
        true>
KOKKOS_FUNCTION auto RefElementTransform(const Geometry& geo, A2DObj<Space>& in,
                                         T& detJ, A2DObj<Space>& out) {
  return RefElementTransformConstGeoExpr<Geometry, A2DObj<Space>, T>(geo, in,
                                                                     detJ, out);
}

template <class Geometry, class Space, class dtype>
class RefElementTransformExpr {
 public:
  // Set the numeric type to use
  typedef typename get_object_numeric_type<Geometry>::type T;

  // In this case, the geometry and space objects must be the same AD type
  static constexpr ADorder order = get_diff_order<Geometry>::order;
  static_assert((order == get_diff_order<Space>::order &&
                 order == get_diff_order<dtype>::order),
                "All objects must be the same AD order");

  static const int dim = remove_a2dobj<Geometry>::type::dim;

  // Set the type of matrix to use - here use a constant matrix
  using JType = ADObjSelect<ADiffType::ACTIVE, order, Mat<T, dim, dim>&>;
  using JinvType = ADObjSelect<ADiffType::ACTIVE, order, Mat<T, dim, dim>>;

  KOKKOS_FUNCTION RefElementTransformExpr(Geometry& geo_, Space& in_,
                                          dtype& detJ_, Space& out_)
      : geo(geo_),
        in(in_),
        detJ(detJ_),
        out(out_),
        J(get_grad<0>(geo)),
        det(J, detJ),
        inv(J, Jinv) {}

  KOKKOS_FUNCTION RefElementTransformExpr(const RefElementTransformExpr& src)
      : geo(src.geo),
        in(src.in),
        detJ(src.detJ),
        out(src.out),
        J(src.J),
        det(J, src.detJ),
        inv(J, Jinv) {}

  KOKKOS_FUNCTION void eval() {
    det.eval();
    inv.eval();
    in.value().transform(detJ.value(), J.value(), Jinv.value(), out.value());
  }

  KOKKOS_FUNCTION void bzero() {
    out.bzero();
    detJ.bzero();
  }

  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;

    det.template forward<forder>();
    inv.template forward<forder>();
    if constexpr (forder == ADorder::FIRST) {
      in.value().template forward_transform<forder>(in.bvalue(), detJ, J, Jinv,
                                                    out.bvalue());
    } else {
      in.value().template forward_transform<forder>(in.pvalue(), detJ, J, Jinv,
                                                    out.pvalue());
    }
  }

  KOKKOS_FUNCTION void reverse() {
    in.value().reverse_transform(out.bvalue(), detJ, J, Jinv, in.bvalue());
    det.reverse();
    inv.reverse();
  }

  KOKKOS_FUNCTION void hzero() {
    out.hzero();
    detJ.hzero();
  }

  KOKKOS_FUNCTION void hreverse() {
    in.value().hreverse_transform(out.hvalue(), out.bvalue(), in.pvalue(), detJ,
                                  J, Jinv, in.hvalue());
    det.hreverse();
    inv.hreverse();
  }

 private:
  // Input/output data
  Geometry& geo;
  Space& in;
  dtype& detJ;
  Space& out;

  // Internal data
  JType J;
  JinvType Jinv;

  // The expression objects for determinant and inverse operations
  MatDetExpr<JType, dtype> det;
  MatInvExpr<JType, JinvType> inv;
};

template <class Geometry, class Space, typename T>
KOKKOS_FUNCTION auto RefElementTransform(ADObj<Geometry>& geo, ADObj<Space>& in,
                                         ADObj<T>& detJ, ADObj<Space>& out) {
  return RefElementTransformExpr<ADObj<Geometry>, ADObj<Space>, ADObj<T>>(
      geo, in, detJ, out);
}

template <class Geometry, class Space, typename T>
KOKKOS_FUNCTION auto RefElementTransform(A2DObj<Geometry>& geo,
                                         A2DObj<Space>& in, A2DObj<T>& detJ,
                                         A2DObj<Space>& out) {
  return RefElementTransformExpr<A2DObj<Geometry>, A2DObj<Space>, A2DObj<T>>(
      geo, in, detJ, out);
}

template <
    class Geometry, class Space, typename T,
    std::enable_if_t<get_diff_type<Geometry>::diff_type == ADiffType::PASSIVE,
                     bool> = true,
    std::enable_if_t<get_diff_type<Space>::diff_type == ADiffType::PASSIVE,
                     bool> = true,
    std::enable_if_t<get_diff_type<T>::diff_type == ADiffType::PASSIVE, bool> =
        true>
KOKKOS_FUNCTION void BoundaryElementTransform(const Geometry& geo,
                                              const Space& in, T& detJ,
                                              Space& out) {
  const Mat<T, Geometry::dim, Geometry::dim - 1>& J = get_grad<0>(geo);
  Mat<T, Geometry::dim, Geometry::dim> Jinv;

  if constexpr (Geometry::dim == 2) {
    Vec<T, 2> x0, n;
    MatColumnToVec(0, J, x0);
    VecNorm(n, detJ);
  } else if constexpr (Geometry::dim == 3) {
    Vec<T, 3> x0, x1, n;
    MatColumnToVec(0, J, x0);
    MatColumnToVec(1, J, x1);
    VecCross(x0, x1, n);
    VecNorm(n, detJ);
  }

  in.transform(detJ, J, Jinv, out);
}

template <class Geometry, class Space, typename T>
class BoundaryElementTransformConstGeoExpr {
 public:
  // In this case, the geometry must be passive
  static_assert(get_diff_type<Geometry>::diff_type == ADiffType::PASSIVE,
                "Must use passive geometry type");

  KOKKOS_FUNCTION BoundaryElementTransformConstGeoExpr(const Geometry& geo,
                                                       Space& in, T& detJ,
                                                       Space& out)
      : geo(geo), in(in), detJ(detJ), out(out) {}

  KOKKOS_FUNCTION void eval() {
    const Mat<T, Geometry::dim, Geometry::dim - 1>& J = get_grad<0>(geo);

    if constexpr (Geometry::dim == 2) {
      Vec<T, 2> x0;
      MatColumnToVec(0, J, x0);
      VecNorm(x0, detJ);
    } else if constexpr (Geometry::dim == 3) {
      Vec<T, 3> x0, x1, n;
      MatColumnToVec(0, J, x0);
      MatColumnToVec(1, J, x1);
      VecCross(x0, x1, n);
      VecNorm(n, detJ);
    }

    in.value().transform(detJ, J, Jinv, out.value());
  }

  KOKKOS_FUNCTION void bzero() { out.bzero(); }

  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {}

  KOKKOS_FUNCTION void reverse() {
    const Mat<T, Geometry::dim, Geometry::dim - 1>& J = get_grad<0>(geo);
    out.bvalue().btransform(detJ, J, Jinv, in.bvalue());
  }

  KOKKOS_FUNCTION void hzero() { out.hzero(); }

  KOKKOS_FUNCTION void hreverse() {}

 private:
  // Internal data
  Vec<T, Geometry::dim> x0, x1, n;
  Mat<T, Geometry::dim, Geometry::dim> Jinv;

  // Input/output data
  const Geometry& geo;
  Space& in;
  T& detJ;
  Space& out;
};

template <class Geometry, class Space, typename T>
KOKKOS_FUNCTION auto BoundaryElementTransform(const Geometry& geo,
                                              ADObj<Space>& in, T& detJ,
                                              ADObj<Space>& out) {
  return BoundaryElementTransformConstGeoExpr<Geometry, ADObj<Space>, T>(
      geo, in, detJ, out);
}

template <class Geometry, class Space, typename T>
KOKKOS_FUNCTION auto BoundaryElementTransform(const Geometry& geo,
                                              A2DObj<Space>& in, T& detJ,
                                              A2DObj<Space>& out) {
  return BoundaryElementTransformConstGeoExpr<Geometry, A2DObj<Space>, T>(
      geo, in, detJ, out);
}

namespace Test {

template <typename T, class Geometry, class Space>
class RefElementTransformConstGeoTest
    : public A2D::Test::A2DTest<T, Space, Space> {
 public:
  using Input = VarTuple<T, Space>;
  using Output = VarTuple<T, Space>;

  RefElementTransformConstGeoTest() { this->set_rand(Geometry::ncomp, geo); }

  // the constant geometry object
  Geometry geo;

  // Assemble a string to describe the test
  std::string name() { return "RefElementTransformConstGeo"; }

  // Evaluate the function
  Output eval(const Input& x) {
    Space in, out;
    T detJ;
    x.get_values(in);
    RefElementTransform(geo, in, detJ, out);
    return MakeVarTuple<T>(out);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& x, Input& g) {
    ADObj<Space> in, out;
    T detJ;
    x.get_values(in.value());
    auto stack = MakeStack(RefElementTransform(geo, in, detJ, out));
    seed.get_values(out.bvalue());
    stack.reverse();
    g.set_values(in.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) {
    A2DObj<Space> in, out;
    T detJ;
    x.get_values(in.value());
    p.get_values(in.pvalue());
    auto stack = MakeStack(RefElementTransform(geo, in, detJ, out));
    seed.get_values(out.bvalue());
    hval.get_values(out.hvalue());
    stack.reverse();
    stack.hforward();
    stack.hreverse();
    h.set_values(in.hvalue());
  }
};

template <typename T, class Geometry, class Space>
class RefElementTransformTest
    : public A2D::Test::A2DTest<T, Space, Geometry, Space> {
 public:
  using Input = VarTuple<T, Geometry, Space>;
  using Output = VarTuple<T, Space>;

  // Assemble a string to describe the test
  std::string name() { return "RefElementTransform"; }

  // Evaluate the function
  Output eval(const Input& x) {
    Geometry geo;
    Space in, out;
    T detJ;
    x.get_values(geo, in);
    RefElementTransform(geo, in, detJ, out);
    return MakeVarTuple<T>(out);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& x, Input& g) {
    ADObj<Geometry> geo;
    ADObj<Space> in, out;
    ADObj<T> detJ;
    x.get_values(geo.value(), in.value());
    auto stack = MakeStack(RefElementTransform(geo, in, detJ, out));
    seed.get_values(out.bvalue());
    stack.reverse();
    g.set_values(geo.bvalue(), in.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) {
    A2DObj<Geometry> geo;
    A2DObj<Space> in, out;
    A2DObj<T> detJ;
    x.get_values(geo.value(), in.value());
    p.get_values(geo.pvalue(), in.pvalue());
    auto stack = MakeStack(RefElementTransform(geo, in, detJ, out));
    seed.get_values(out.bvalue());
    hval.get_values(out.hvalue());
    stack.reverse();
    stack.hforward();
    stack.hreverse();
    h.set_values(geo.hvalue(), in.hvalue());
  }
};

template <typename T, class Geometry, class Space>
class RefElementTransformDetTest
    : public A2D::Test::A2DTest<T, T, Geometry, Space> {
 public:
  using Input = VarTuple<T, Geometry, Space>;
  using Output = VarTuple<T, T>;

  // Assemble a string to describe the test
  std::string name() { return "RefElementTransformDetTest"; }

  // Evaluate the function
  Output eval(const Input& x) {
    Geometry geo;
    Space in, out;
    T detJ;
    x.get_values(geo, in);
    RefElementTransform(geo, in, detJ, out);
    return MakeVarTuple<T>(detJ);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& x, Input& g) {
    ADObj<Geometry> geo;
    ADObj<Space> in, out;
    ADObj<T> detJ;
    x.get_values(geo.value(), in.value());
    auto stack = MakeStack(RefElementTransform(geo, in, detJ, out));
    seed.get_values(detJ.bvalue());
    stack.reverse();
    g.set_values(geo.bvalue(), in.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) {
    A2DObj<Geometry> geo;
    A2DObj<Space> in, out;
    A2DObj<T> detJ;
    x.get_values(geo.value(), in.value());
    p.get_values(geo.pvalue(), in.pvalue());
    auto stack = MakeStack(RefElementTransform(geo, in, detJ, out));
    seed.get_values(detJ.bvalue());
    hval.get_values(detJ.hvalue());
    stack.reverse();
    stack.hforward();
    stack.hreverse();
    h.set_values(geo.hvalue(), in.hvalue());
  }
};

template <typename T, class Geometry, class Space>
bool RefElementTransformTestHelper(bool component, bool write_output) {
  bool passed = true;

  RefElementTransformConstGeoTest<T, Geometry, Space> test1;
  passed = passed && A2D::Test::Run(test1, component, write_output);

  RefElementTransformTest<T, Geometry, Space> test2;
  passed = passed && A2D::Test::Run(test2, component, write_output);

  RefElementTransformDetTest<T, Geometry, Space> test3;
  passed = passed && A2D::Test::Run(test3, component, write_output);

  return passed;
}

bool RefElementTransformTestAll(bool component, bool write_output) {
  bool passed = true;

  using T = std::complex<double>;
  const int dim = 3;
  using Geometry = FESpace<T, dim, H1Space<T, dim, dim>>;
  // using Space = FESpace<T, dim, L2Space<T, 4, dim>, L2Space<T, 1, dim>,
  //                       H1Space<T, 2, dim>, H1Space<T, 1, dim>>;
  // using Space = FESpace<T, dim, HdivSpace<T, dim>>;
  using Space = FESpace<T, dim, H1Space<T, 2, dim>>;

  passed = RefElementTransformTestHelper<T, Geometry, Space>(component,
                                                             write_output);

  return passed;
}

}  // namespace Test

/**
 * @brief 3D or 2D volume transform transform
 *
 */
template <typename T, index_t dim>
class InteriorMapping {
 public:
  template <class FiniteElementGeometry>
  KOKKOS_FUNCTION InteriorMapping(const FiniteElementGeometry& geo, T& detJ)
      : J(get_grad<0>(geo)), detJ(detJ) {
    // Compute the inverse of the transformation
    MatInv(J, Jinv);

    // Compute the determinant of the Jacobian matrix
    MatDet(J, detJ);
  }

  template <class FiniteElementSpace>
  KOKKOS_FUNCTION void transform(const FiniteElementSpace& in,
                                 FiniteElementSpace& out) {
    in.transform(detJ, J, Jinv, out);
  }

  template <class FiniteElementSpace>
  KOKKOS_FUNCTION void btransform(const FiniteElementSpace& in,
                                  FiniteElementSpace& out) {
    in.btransform(detJ, J, Jinv, out);
  }

  template <class FiniteElementSpace, class QMatType>
  KOKKOS_FUNCTION void jtransform(const QMatType& mat_in, QMatType& mat_out) {
    constexpr index_t ncomp = QMatType::nrows;

    FiniteElementSpace pref, p, Jp;

    for (index_t k = 0; k < ncomp; k++) {
      pref.zero();
      pref[k] = T(1.0);
      transform(pref, p);

      // Compute mat-vec multiplication
      Jp.zero();

      // TODO: use MatVec operation instead
      for (index_t j = 0; j < ncomp; j++) {
        for (index_t i = 0; i < ncomp; i++) {
          Jp[i] += mat_in(i, j) * p[j];
        }
      }
      btransform(Jp, pref);
      for (index_t i = 0; i < ncomp; i++) {
        mat_out(i, k) = pref[i];
      }
    }
  }

 private:
  const Mat<T, dim, dim>& J;
  Mat<T, dim, dim> Jinv;
  T& detJ;
};

/**
 * @brief surface transformation - TODO: this needs to be fixed
 */
template <typename T, index_t dim>
class SurfaceMapping {
 public:
  static const index_t dim_surf = dim - 1;

  template <class FiniteElementGeometry>
  KOKKOS_FUNCTION SurfaceMapping(const FiniteElementGeometry& geo, T& detJ)
      : detJ(detJ) {
    const Mat<T, dim, dim>& Jxi = geo.template get<0>().get_grad();
    Vec<T, dim> x, y, nA;
    if constexpr (dim == 2) {
      // Find the nA = vector of the distorted bound
      nA(0) = Jxi(0, 0);
      nA(1) = Jxi(1, 0);
      detJ = std::sqrt(nA(0) * nA(0) + nA(1) * nA(1));
    } else if constexpr (dim == 3) {
      // Find the nA = (Area) * normal direction
      x(0) = Jxi(0, 0);  // d(x, y, z)/dxi
      x(1) = Jxi(1, 0);
      x(2) = Jxi(2, 0);

      y(0) = Jxi(0, 1);  // d(x, y, z)/deta
      y(1) = Jxi(1, 1);
      y(2) = Jxi(2, 1);

      nA(0) = x(1) * y(2) - x(2) * y(1);
      nA(1) = x(2) * y(0) - x(0) * y(2);
      nA(2) = x(0) * y(1) - x(1) * y(0);

      detJ = std::sqrt(nA(0) * nA(0) + nA(1) * nA(1) + nA(2) * nA(2));
    }
  }

  template <class FiniteElementSpace>
  KOKKOS_FUNCTION void transform(const FiniteElementSpace& in,
                                 FiniteElementSpace& out) {
    in.transform(detJ, J, Jinv, out);
  }

  template <class FiniteElementSpace>
  KOKKOS_FUNCTION void btransform(const FiniteElementSpace& in,
                                  FiniteElementSpace& out) {
    in.btransform(detJ, J, Jinv, out);
  }

 private:
  // Determinant of the Jacobian transformation
  T& detJ;

  // J with the normal direction added
  Mat<T, dim_surf, dim_surf> J, Jinv;  // TODO
};

}  // namespace A2D

#endif  //  A2D_FE_MAPPING_H
