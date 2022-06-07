#ifndef BLOCK_NUMERIC_H
#define BLOCK_NUMERIC_H

namespace A2D {

double fabs(std::complex<double> a) {
  if (a.real() >= 0.0) {
    return a.real();
  } else {
    return -a.real();
  }
}

double fabs(double a) {
  if (a >= 0.0) {
    return a;
  } else {
    return -a;
  }
}

/*
  Compute y = A * x
*/
template <typename T, int M, int N, class AType, class xType, class yType>
inline void blockGemv(const AType& A, const xType& x, yType& y) {
  for (int i = 0; i < M; i++) {
    T prod = 0.0;
    for (int j = 0; j < N; j++) {
      prod += A(i, j) * x(j);
    }
    y(i) = prod;
  }
}

/*
  Compute y += A * x
*/
template <typename T, int M, int N, class AType, class xType, class yType>
inline void blockGemvAdd(const AType& A, const xType& x, yType& y) {
  for (int i = 0; i < M; i++) {
    T prod = 0.0;
    for (int j = 0; j < N; j++) {
      prod += A(i, j) * x(j);
    }
    y(i) += prod;
  }
}

/*
  Compute y -= A * x
*/
template <typename T, int M, int N, class AType, class xType, class yType>
inline void blockGemvSub(const AType& A, const xType& x, yType& y) {
  for (int i = 0; i < M; i++) {
    T prod = 0.0;
    for (int j = 0; j < N; j++) {
      prod += A(i, j) * x(j);
    }
    y(i) -= prod;
  }
}

/*
  Compute y += scale * A * x
*/
template <typename T, int M, int N, class AType, class xType, class yType>
inline void blockGemvScale(T scale, const AType& A, const xType& x, yType& y) {
  for (int i = 0; i < M; i++) {
    T prod = 0.0;
    for (int j = 0; j < N; j++) {
      prod += A(i, j) * x(j);
    }
    y(i) += scale * prod;
  }
}

/*
  Compute: C = A * B

  A in M x N
  B in N x P
  C in M x P
*/
template <typename T, int M, int N, int P, class AType, class BType,
          class CType>
inline void blockGemm(const AType& A, const BType& B, CType& C) {
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < P; j++) {
      T prod = 0.0;
      for (int k = 0; k < N; k++) {
        prod += A(i, k) * B(k, j);
      }
      C(i, j) = prod;
    }
  }
}

/*
  Compute: C += A * B

  A in M x N
  B in N x P
  C in M x P
*/
template <typename T, int M, int N, int P, class AType, class BType,
          class CType>
inline void blockGemmAdd(const AType& A, const BType& B, CType& C) {
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < P; j++) {
      T prod = 0.0;
      for (int k = 0; k < N; k++) {
        prod += A(i, k) * B(k, j);
      }
      C(i, j) += prod;
    }
  }
}

/*
  Compute: C -= A * B

  A in M x N
  B in N x P
  C in M x P
*/
template <typename T, int M, int N, int P, class AType, class BType,
          class CType>
inline void blockGemmSub(const AType& A, const BType& B, CType& C) {
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < P; j++) {
      T prod = 0.0;
      for (int k = 0; k < N; k++) {
        prod += A(i, k) * B(k, j);
      }
      C(i, j) -= prod;
    }
  }
}

/*
  Compute: C += scale * A * B

  A in M x N
  B in N x P
  C in M x P
*/
template <typename T, int M, int N, int P, class AType, class BType,
          class CType>
inline void blockGemmScale(T scale, const AType& A, const BType& B, CType& C) {
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < P; j++) {
      T prod = 0.0;
      for (int k = 0; k < N; k++) {
        prod += A(i, k) * B(k, j);
      }
      C(i, j) += scale * prod;
    }
  }
}

/*
  Compute: Ainv = A^{-1} with pivoting
*/
template <typename T, int N, class AType, class AinvType, class IType>
int blockInverse(AType& A, AinvType& Ainv, IType& ipiv) {
  int fail = 0;

  for (int k = 0; k < N - 1; k++) {
    // Find the maximum value and use it as the pivot
    int r = k;
    T maxv = A(k, k);
    for (int j = k + 1; j < N; j++) {
      T t = A(j, k);
      if (fabs(t) > fabs(maxv)) {
        maxv = t;
        r = j;
      }
    }

    ipiv(k) = r;

    // If a swap is required, swap the rows
    if (r != k) {
      for (int j = 0; j < N; j++) {
        T t = A(k, j);
        A(k, j) = A(r, j);
        A(r, j) = t;
      }
    }

    if (fabs(A(k, k)) == 0.0) {
      fail = k + 1;
      return fail;
    }

    for (int i = k + 1; i < N; i++) {
      A(i, k) = A(i, k) / A(k, k);
    }

    for (int i = k + 1; i < N; i++) {
      for (int j = k + 1; j < N; j++) {
        A(i, j) -= A(i, k) * A(k, j);
      }
    }
  }

  // Now, compute the matrix-inverse
  for (int k = 0; k < N; k++) {
    int ip = k;
    for (int i = 0; i < N - 1; i++) {
      if (ip == ipiv(i)) {
        ip = i;
      } else if (ip == i) {
        ip = ipiv(i);
      }
    }

    for (int i = 0; i < ip; i++) {
      Ainv(i, k) = 0.0;
    }

    Ainv(ip, k) = 1.0;

    for (int i = ip + 1; i < N; i++) {
      Ainv(i, k) = 0.0;
      for (int j = ip; j < i; j++) {
        Ainv(i, k) -= A(i, j) * Ainv(j, k);
      }
    }

    for (int i = N - 1; i >= 0; i--) {
      for (int j = i + 1; j < N; j++) {
        Ainv(i, k) -= A(i, j) * Ainv(j, k);
      }
      Ainv(i, k) = Ainv(i, k) / A(i, i);
    }
  }

  return fail;
}

}  // namespace A2D

#endif  // BLOCK_NUMERIC_H
