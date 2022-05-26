#ifndef A2D_OBJS_H
#define A2D_OBJS_H

namespace A2D {

template <typename T, int N>
class Vec {
 public:
  typedef T type;

  Vec() {
    for (int i = 0; i < N; i++) {
      x[i] = 0.0;
    }
  }
  Vec(const T* vals) {
    for (int i = 0; i < N; i++) {
      x[i] = vals[i];
    }
  }
  template <class VecType>
  Vec(const VecType& vec) {
    for (int i = 0; i < N; i++) {
      x[i] = vec(i);
    }
  }
  void zero() {
    for (int i = 0; i < N; i++) {
      x[i] = 0.0;
    }
  }
  template <class IdxType>
  T& operator()(const IdxType i) {
    return x[i];
  }
  template <class IdxType>
  const T& operator()(const IdxType i) const {
    return x[i];
  }

  T x[N];
};

template <typename T, int M, int N>
class Mat {
 public:
  typedef T type;

  Mat() {
    for (int i = 0; i < M * N; i++) {
      A[i] = 0.0;
    }
  }
  Mat(const T* vals) {
    for (int i = 0; i < M * N; i++) {
      A[i] = vals[i];
    }
  }
  template <class MatType>
  Mat(const MatType& mat) {
    for (int i = 0; i < M; i++) {
      for (int j = 0; j < N; j++) {
        A[N * i + j] = mat(i, j);
      }
    }
  }
  void zero() {
    for (int i = 0; i < M * N; i++) {
      A[i] = 0.0;
    }
  }
  template <class IdxType>
  T& operator()(const IdxType i, const IdxType j) {
    return A[N * i + j];
  }
  template <class IdxType>
  const T& operator()(const IdxType i, const IdxType j) const {
    return A[N * i + j];
  }

  T A[M * N];
};

template <typename T, int N>
class SymmMat {
 public:
  typedef T type;
  static const int MAT_SIZE = (N * (N + 1)) / 2;

  SymmMat() {
    for (int i = 0; i < MAT_SIZE; i++) {
      A[i] = 0.0;
    }
  }
  SymmMat(const T* vals) {
    for (int i = 0; i < MAT_SIZE; i++) {
      A[i] = vals[i];
    }
  }
  void zero() {
    for (int i = 0; i < MAT_SIZE; i++) {
      A[i] = 0.0;
    }
  }
  template <class IdxType>
  T& operator()(const IdxType i, const IdxType j) {
    if (i >= j) {
      return A[j + i * (i + 1) / 2];
    } else {
      return A[i + j * (j + 1) / 2];
    }
  }
  template <class IdxType>
  const T& operator()(const IdxType i, const IdxType j) const {
    if (i >= j) {
      return A[j + i * (i + 1) / 2];
    } else {
      return A[i + j * (j + 1) / 2];
    }
  }

  T A[MAT_SIZE];
};

/*
  Full 4-th order tensor without any symmetry.

  This tensor class is required for mixed second order derivatives of
  two matrices X and Y such that

  A(i, j, k, l) = d^2 f/dX(i, j) dY(k, l)
*/
template <typename T, int M, int N, int P, int Q>
class Tensor {
 public:
  typedef T type;
  static const int TENSOR_SIZE = M * N * P * Q;
  Tensor() {
    for (int i = 0; i < TENSOR_SIZE; i++) {
      A[i] = 0.0;
    }
  }
  template <class IdxType>
  T& operator()(const IdxType i, const IdxType j, const IdxType k,
                const IdxType l) {
    return A[l + Q * (k + P * (j + N * i))];
  }
  template <class IdxType>
  const T& operator()(const IdxType i, const IdxType j, const IdxType k,
                      const IdxType l) const {
    return A[l + Q * (k + P * (j + N * i))];
  }

  T A[TENSOR_SIZE];
};

/*
  Basic 4-th order symmetric tensor.

  This class stores a symmetric tensor found by taking the second order
  derivatives of a non-symmetric matrix

  A(i, j, k, l) = d^2 f/dX(i, j) dX(k, l)

  As a result:

  A(i, j, k, l) = A(k, l, i, j).
*/
template <typename T, int M, int N>
class SymmTensor {
 public:
  typedef T type;
  static const int TENSOR_SIZE = (M * N * (M * N + 1)) / 2;
  SymmTensor() {
    for (int i = 0; i < TENSOR_SIZE; i++) {
      A[i] = 0.0;
    }
  }
  template <class IdxType>
  T& operator()(const IdxType i, const IdxType j, const IdxType k,
                const IdxType l) {
    const int ii = N * i + j;
    const int jj = N * k + l;

    if (ii >= jj) {
      return A[jj + ii * (ii + 1) / 2];
    } else {
      return A[ii + jj * (jj + 1) / 2];
    }
  }
  template <class IdxType>
  const T& operator()(const IdxType i, const IdxType j, const IdxType k,
                      const IdxType l) const {
    const int ii = N * i + j;
    const int jj = N * k + l;

    if (ii >= jj) {
      return A[jj + ii * (ii + 1) / 2];
    } else {
      return A[ii + jj * (jj + 1) / 2];
    }
  }

  T A[TENSOR_SIZE];
};

/*
  Mixed symmetric

  A(i, j, k, l) = d^2 f/d X[i, j] d Y[k, l]
*/
/*
template <typename T, int M, int N, P>
class Symm {
  typedef T type;
  static const int TENSOR_SIZE = ;
};
*/
/*


  A(i, j, k, l) = d^2 f/dX[i, j] dX[k, l]

  where X[i, j] = X[j, i]

  As a result:

  A(i, j, k, l) = A(k, l, i, j)
  A(i, j, k, l) = A(j, i, k, l)
  A(i, j, k, l) = A(i, j, l, k)
*/
/*
template <typename T, int M>
class SymmSymmTensor {
  typedef T type;
  static const int TENSOR_SIZE = ;
}
*/
}  // namespace A2D

#endif  // A2D_OBJS_H
