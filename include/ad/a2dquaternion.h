#ifndef A2D_QUATERNION_H
#define A2D_QUATERNION_H

#include <iostream>

#include "../a2ddefs.h"
#include "a2dgemm.h"
#include "a2dmat.h"
#include "a2dvec.h"

namespace A2D {

template <typename T>
void QuaternionMatrixCore(const T q[], T C[]) {
  C[0] = 1.0 - 2.0 * (q[2] * q[2] + q[3] * q[3]);
  C[1] = 2.0 * (q[1] * q[2] + q[3] * q[0]);
  C[2] = 2.0 * (q[1] * q[3] - q[2] * q[0]);

  C[3] = 2.0 * (q[2] * q[1] - q[3] * q[0]);
  C[4] = 1.0 - 2.0 * (q[1] * q[1] + q[3] * q[3]);
  C[5] = 2.0 * (q[2] * q[3] + q[1] * q[0]);

  C[6] = 2.0 * (q[3] * q[1] + q[2] * q[0]);
  C[7] = 2.0 * (q[3] * q[2] - q[1] * q[0]);
  C[8] = 1.0 - 2.0 * (q[1] * q[1] + q[2] * q[2]);
}

template <typename T>
void QuaternionMatrixForwardCore(const T q[], const T qd[], T Cd[]) {
  Cd[0] = -4.0 * (q[2] * qd[2] + q[3] * qd[3]);
  Cd[1] = 2.0 * (q[1] * qd[2] + q[3] * qd[0] + qd[1] * q[2] + qd[3] * q[0]);
  Cd[2] = 2.0 * (q[1] * qd[3] - q[2] * qd[0] + qd[1] * q[3] - qd[2] * q[0]);

  Cd[3] = 2.0 * (q[2] * qd[1] - q[3] * qd[0] + qd[2] * q[1] - qd[3] * q[0]);
  Cd[4] = -4.0 * (q[1] * qd[1] + q[3] * qd[3]);
  Cd[5] = 2.0 * (q[2] * qd[3] + q[1] * qd[0] + qd[2] * q[3] + qd[1] * q[0]);

  Cd[6] = 2.0 * (q[3] * qd[1] + q[2] * qd[0] + qd[3] * q[1] + qd[2] * q[0]);
  Cd[7] = 2.0 * (q[3] * qd[2] - q[1] * qd[0] + qd[3] * q[2] - qd[1] * q[0]);
  Cd[8] = -4.0 * (q[1] * qd[1] + q[2] * qd[2]);
}

template <typename T>
void QuaternionMatrixReverseCore(const T q[], const T dC[], T r[]) {
  r[0] += 2.0 * (q[3] * (dC[1] - dC[3]) + q[2] * (dC[6] - dC[2]) +
                 q[1] * (dC[5] - dC[7]));
  r[1] += 2.0 * (q[0] * (dC[5] - dC[7]) - 2.0 * q[1] * (dC[4] + dC[8]) +
                 q[2] * (dC[1] + dC[3]) + q[3] * (dC[2] + dC[6]));
  r[2] += 2.0 * (q[0] * (dC[6] - dC[2]) + q[1] * (dC[1] + dC[3]) -
                 2.0 * q[2] * (dC[0] + dC[8]) + q[3] * (dC[7] + dC[5]));
  r[3] += 2.0 * (q[0] * (dC[1] - dC[3]) + q[1] * (dC[2] + dC[6]) +
                 q[2] * (dC[5] + dC[7]) - 2.0 * q[3] * (dC[0] + dC[4]));
}

template <typename T>
void QuaternionMatrix(const Vec<T, 4>& q, Mat<T, 3, 3>& C) {
  QuaternionMatrixCore<T>(get_data(q), get_data(C));
}

template <class qtype, class Ctype>
class QuaternionMatrixExpr {
 public:
  // Extract the numeric type to use
  using T = object_scalar_t<Ctype>;

  // Extract the underlying sizes of the matrix
  static constexpr int M = get_matrix_rows<Ctype>::size;
  static constexpr int N = get_matrix_columns<Ctype>::size;

  // Extract the dimensions of the underlying vectors
  static constexpr int L = get_vec_size<qtype>::size;

  static_assert(M == N && N == 3, "Rotation matrix dimension must be 3");
  static_assert(L == 4, "Quaternion dimension must be 4");

  A2D_FUNCTION QuaternionMatrixExpr(qtype& q, Ctype& C) : q(q), C(C) {}

  A2D_FUNCTION void eval() {
    QuaternionMatrixCore<T>(get_data(q), get_data(C));
  }

  A2D_FUNCTION void bzero() { C.bzero(); }

  template <ADorder forder>
  A2D_FUNCTION void forward() {
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;
    QuaternionMatrixForwardCore<T>(get_data(q), GetSeed<seed>::get_data(q),
                                   GetSeed<seed>::get_data(C));
  }

  A2D_FUNCTION void reverse() {
    constexpr ADseed seed = ADseed::b;
    QuaternionMatrixReverseCore<T>(get_data(q), GetSeed<seed>::get_data(C),
                                   GetSeed<seed>::get_data(q));
  }

  A2D_FUNCTION void hzero() { C.hzero(); }

  A2D_FUNCTION void hreverse() {
    QuaternionMatrixReverseCore<T>(get_data(q), GetSeed<ADseed::h>::get_data(C),
                                   GetSeed<ADseed::h>::get_data(q));
    QuaternionMatrixReverseCore<T>(GetSeed<ADseed::p>::get_data(q),
                                   GetSeed<ADseed::b>::get_data(C),
                                   GetSeed<ADseed::h>::get_data(q));
  }

 private:
  qtype& q;
  Ctype& C;
};

template <class qtype, class Ctype>
A2D_FUNCTION auto QuaternionMatrix(ADObj<qtype>& q, ADObj<Ctype>& C) {
  return QuaternionMatrixExpr<ADObj<qtype>, ADObj<Ctype>>(q, C);
}

template <class qtype, class Ctype>
A2D_FUNCTION auto QuaternionMatrix(A2DObj<qtype>& q, A2DObj<Ctype>& C) {
  return QuaternionMatrixExpr<A2DObj<qtype>, A2DObj<Ctype>>(q, C);
}

template <typename T>
void QuaternionAngularVelocityCore(const T q[], const T qdot[], T omega[]) {
  omega[0] =
      2.0 * (q[0] * qdot[1] - qdot[0] * q[1] + q[3] * qdot[2] - q[2] * qdot[3]);
  omega[1] =
      2.0 * (q[0] * qdot[2] - qdot[0] * q[2] + q[1] * qdot[3] - q[3] * qdot[1]);
  omega[2] =
      2.0 * (q[0] * qdot[3] - qdot[0] * q[3] + q[2] * qdot[1] - q[1] * qdot[2]);
}

template <typename T>
void QuaternionAngularVelocityForwardCore(const T q[], const T qdot[],
                                          const T qd[], const T qdotd[],
                                          T omega[]) {
  omega[0] = 2.0 * (qd[0] * qdot[1] - qdotd[0] * q[1] + qd[3] * qdot[2] -
                    qd[2] * qdot[3] + q[0] * qdotd[1] - qdot[0] * qd[1] +
                    q[3] * qdotd[2] - q[2] * qdotd[3]);
  omega[1] = 2.0 * (qd[0] * qdot[2] - qdotd[0] * q[2] + qd[1] * qdot[3] -
                    qd[3] * qdot[1] + q[0] * qdotd[2] - qdot[0] * qd[2] +
                    q[1] * qdotd[3] - q[3] * qdotd[1]);
  omega[2] = 2.0 * (qd[0] * qdot[3] - qdotd[0] * q[3] + qd[2] * qdot[1] -
                    qd[1] * qdot[2] + q[0] * qdotd[3] - qdot[0] * qd[3] +
                    q[2] * qdotd[1] - q[1] * qdotd[2]);
}

template <typename T>
void QuaternionAngularVelocityReverseCore(const T q[], const T qdot[],
                                          const T omegab[], T qb[], T qdotb[]) {
  qb[0] +=
      2.0 * (qdot[1] * omegab[0] + qdot[2] * omegab[1] + qdot[3] * omegab[2]);
  qb[1] +=
      2.0 * (-qdot[0] * omegab[0] + qdot[3] * omegab[1] - qdot[2] * omegab[2]);
  qb[2] +=
      2.0 * (-qdot[3] * omegab[0] - qdot[0] * omegab[1] + qdot[1] * omegab[2]);
  qb[3] +=
      2.0 * (qdot[2] * omegab[0] - qdot[1] * omegab[1] - qdot[0] * omegab[2]);

  qdotb[0] -= 2.0 * (q[1] * omegab[0] + q[2] * omegab[1] + q[3] * omegab[2]);
  qdotb[1] += 2.0 * (q[0] * omegab[0] - q[3] * omegab[1] + q[2] * omegab[2]);
  qdotb[2] += 2.0 * (q[3] * omegab[0] + q[0] * omegab[1] - q[1] * omegab[2]);
  qdotb[3] += 2.0 * (-q[2] * omegab[0] + q[1] * omegab[1] + q[0] * omegab[2]);
}

template <typename T>
void QuaternionAngularVelocity(const Vec<T, 4>& q, const Vec<T, 4>& qdot,
                               Vec<T, 3>& omega) {
  QuaternionAngularVelocityCore<T>(get_data(q), get_data(qdot),
                                   get_data(omega));
}

template <class qtype, class wtype>
class QuaternionAngularVelocityExpr {
 public:
  // Extract the numeric type to use
  using T = object_scalar_t<wtype>;

  // Extract the underlying sizes of the matrix
  static constexpr int M = get_vec_size<wtype>::size;

  // Extract the dimensions of the underlying vectors
  static constexpr int L = get_vec_size<qtype>::size;

  static_assert(M == 3, "Rotation matrix dimension must be 3");
  static_assert(L == 4, "Quaternion dimension must be 4");

  A2D_FUNCTION QuaternionAngularVelocityExpr(qtype& q, qtype& qdot,
                                             wtype& omega)
      : q(q), qdot(qdot), omega(omega) {}

  A2D_FUNCTION void eval() {
    QuaternionAngularVelocityCore<T>(get_data(q), get_data(qdot),
                                     get_data(omega));
  }

  A2D_FUNCTION void bzero() { omega.bzero(); }

  template <ADorder forder>
  A2D_FUNCTION void forward() {
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;
    QuaternionAngularVelocityForwardCore<T>(
        get_data(q), get_data(qdot), GetSeed<seed>::get_data(q),
        GetSeed<seed>::get_data(qdot), GetSeed<seed>::get_data(omega));
  }

  A2D_FUNCTION void reverse() {
    constexpr ADseed seed = ADseed::b;
    QuaternionAngularVelocityReverseCore<T>(
        get_data(q), get_data(qdot), GetSeed<seed>::get_data(omega),
        GetSeed<seed>::get_data(q), GetSeed<seed>::get_data(qdot));
  }

  A2D_FUNCTION void hzero() { omega.hzero(); }

  A2D_FUNCTION void hreverse() {
    QuaternionAngularVelocityReverseCore<T>(
        get_data(q), get_data(qdot), GetSeed<ADseed::h>::get_data(omega),
        GetSeed<ADseed::h>::get_data(q), GetSeed<ADseed::h>::get_data(qdot));

    QuaternionAngularVelocityReverseCore<T>(
        GetSeed<ADseed::p>::get_data(q), GetSeed<ADseed::p>::get_data(qdot),
        GetSeed<ADseed::b>::get_data(omega), GetSeed<ADseed::h>::get_data(q),
        GetSeed<ADseed::h>::get_data(qdot));
  }

 private:
  qtype& q;
  qtype& qdot;
  wtype& omega;
};

template <class qtype, class wtype>
A2D_FUNCTION auto QuaternionAngularVelocity(ADObj<qtype>& q, ADObj<qtype>& qdot,
                                            ADObj<wtype>& omega) {
  return QuaternionAngularVelocityExpr<ADObj<qtype>, ADObj<wtype>>(q, qdot,
                                                                   omega);
}

template <class qtype, class wtype>
A2D_FUNCTION auto QuaternionAngularVelocity(A2DObj<qtype>& q,
                                            A2DObj<qtype>& qdot,
                                            A2DObj<wtype>& omega) {
  return QuaternionAngularVelocityExpr<A2DObj<qtype>, A2DObj<wtype>>(q, qdot,
                                                                     omega);
}

}  // namespace A2D

#endif  // A2D_QUATERNION_H
