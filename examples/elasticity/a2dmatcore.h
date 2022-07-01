#ifndef A2D_MAT_CORE_H
#define A2D_MAT_CORE_H

namespace A2D {

/**
 * @brief mat-mat multiplication C = A * B
 * @param A, M-by-N matrix
 * @param B, N-by-P matrix
 * @param C, M-by-P matrix, output
 */
template <int M, int N, int P, class AMatType, class BMatType, class CMatType>
void MatMatMultCore(const AMatType& A, const BMatType& B, CMatType& C) {
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < P; j++) {
      C(i, j) = 0.0;
      for (int k = 0; k < N; k++) {
        C(i, j) += A(i, k) * B(k, j);
      }
    }
  }
}

/**
 * @brief mat-mat multiplication C = AT * B
 * @param A, N-by-M matrix
 * @param B, N-by-P matrix
 * @param C, M-by-P matrix, output
 */
template <int M, int N, int P, class AMatType, class BMatType, class CMatType>
void MatTransMatMultCore(const AMatType& A, const BMatType& B, CMatType& C) {
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < P; j++) {
      C(i, j) = 0.0;
      for (int k = 0; k < N; k++) {
        C(i, j) += A(k, i) * B(k, j);
      }
    }
  }
}

/**
 * @brief mat-mat multiplication C = A * BT
 * @param A, M-by-N matrix
 * @param B, P-by-N matrix
 * @param C, M-by-P matrix, output
 */
template <int M, int N, int P, class AMatType, class BMatType, class CMatType>
void MatMatTransMultCore(const AMatType& A, const BMatType& B, CMatType& C) {
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < P; j++) {
      C(i, j) = 0.0;
      for (int k = 0; k < N; k++) {
        C(i, j) += A(i, k) * B(j, k);
      }
    }
  }
}

/**
 * @brief mat-mat multiplication C = AT * BT
 * @param A, N-by-M matrix
 * @param B, P-by-N matrix
 * @param C, M-by-P matrix, output
 */
template <int M, int N, int P, class AMatType, class BMatType, class CMatType>
void MatTransMatTransMultCore(const AMatType& A, const BMatType& B,
                              CMatType& C) {
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < P; j++) {
      C(i, j) = 0.0;
      for (int k = 0; k < N; k++) {
        C(i, j) += A(k, i) * B(j, k);
      }
    }
  }
}

template <class AMatType, class BMatType, class CMatType>
inline void Symm3x3SymmMultCore(const AMatType& A, const BMatType& B,
                                CMatType& C) {
  C(0, 0) = (A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0) + A(0, 2) * B(2, 0));
  C(0, 1) = (A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1) + A(0, 2) * B(2, 1));
  C(0, 2) = (A(0, 0) * B(0, 2) + A(0, 1) * B(1, 2) + A(0, 2) * B(2, 2));
  C(1, 0) = (A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0) + A(1, 2) * B(2, 0));
  C(1, 1) = (A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1) + A(1, 2) * B(2, 1));
  C(1, 2) = (A(1, 0) * B(0, 2) + A(1, 1) * B(1, 2) + A(1, 2) * B(2, 2));
  C(2, 0) = (A(2, 0) * B(0, 0) + A(2, 1) * B(1, 0) + A(2, 2) * B(2, 0));
  C(2, 1) = (A(2, 0) * B(0, 1) + A(2, 1) * B(1, 1) + A(2, 2) * B(2, 1));
  C(2, 2) = (A(2, 0) * B(0, 2) + A(2, 1) * B(1, 2) + A(2, 2) * B(2, 2));
}

template <typename T, class AMatType, class BMatType, class CMatType>
inline void Symm3x3SymmMultScaleCore(T scale, const AMatType& A,
                                     const BMatType& B, CMatType& C) {
  C(0, 0) = scale * (A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0) + A(0, 2) * B(2, 0));
  C(0, 1) = scale * (A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1) + A(0, 2) * B(2, 1));
  C(0, 2) = scale * (A(0, 0) * B(0, 2) + A(0, 1) * B(1, 2) + A(0, 2) * B(2, 2));
  C(1, 0) = scale * (A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0) + A(1, 2) * B(2, 0));
  C(1, 1) = scale * (A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1) + A(1, 2) * B(2, 1));
  C(1, 2) = scale * (A(1, 0) * B(0, 2) + A(1, 1) * B(1, 2) + A(1, 2) * B(2, 2));
  C(2, 0) = scale * (A(2, 0) * B(0, 0) + A(2, 1) * B(1, 0) + A(2, 2) * B(2, 0));
  C(2, 1) = scale * (A(2, 0) * B(0, 1) + A(2, 1) * B(1, 1) + A(2, 2) * B(2, 1));
  C(2, 2) = scale * (A(2, 0) * B(0, 2) + A(2, 1) * B(1, 2) + A(2, 2) * B(2, 2));
}

template <class AMatType, class BMatType, class CMatType>
inline void Symm3x3SymmMultAddCore(const AMatType& A, const BMatType& B,
                                   CMatType& C) {
  C(0, 0) += (A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0) + A(0, 2) * B(2, 0));
  C(0, 1) += (A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1) + A(0, 2) * B(2, 1));
  C(0, 2) += (A(0, 0) * B(0, 2) + A(0, 1) * B(1, 2) + A(0, 2) * B(2, 2));
  C(1, 0) += (A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0) + A(1, 2) * B(2, 0));
  C(1, 1) += (A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1) + A(1, 2) * B(2, 1));
  C(1, 2) += (A(1, 0) * B(0, 2) + A(1, 1) * B(1, 2) + A(1, 2) * B(2, 2));
  C(2, 0) += (A(2, 0) * B(0, 0) + A(2, 1) * B(1, 0) + A(2, 2) * B(2, 0));
  C(2, 1) += (A(2, 0) * B(0, 1) + A(2, 1) * B(1, 1) + A(2, 2) * B(2, 1));
  C(2, 2) += (A(2, 0) * B(0, 2) + A(2, 1) * B(1, 2) + A(2, 2) * B(2, 2));
}

template <class AMatType, class BMatType, class CMatType>
inline void Symm3x3SymmMultSubCore(const AMatType& A, const BMatType& B,
                                   CMatType& C) {
  C(0, 0) -= (A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0) + A(0, 2) * B(2, 0));
  C(0, 1) -= (A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1) + A(0, 2) * B(2, 1));
  C(0, 2) -= (A(0, 0) * B(0, 2) + A(0, 1) * B(1, 2) + A(0, 2) * B(2, 2));
  C(1, 0) -= (A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0) + A(1, 2) * B(2, 0));
  C(1, 1) -= (A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1) + A(1, 2) * B(2, 1));
  C(1, 2) -= (A(1, 0) * B(0, 2) + A(1, 1) * B(1, 2) + A(1, 2) * B(2, 2));
  C(2, 0) -= (A(2, 0) * B(0, 0) + A(2, 1) * B(1, 0) + A(2, 2) * B(2, 0));
  C(2, 1) -= (A(2, 0) * B(0, 1) + A(2, 1) * B(1, 1) + A(2, 2) * B(2, 1));
  C(2, 2) -= (A(2, 0) * B(0, 2) + A(2, 1) * B(1, 2) + A(2, 2) * B(2, 2));
}

template <typename T, class AMatType, class BMatType, class CMatType>
inline void Symm3x3SymmMultAddScaleCore(T scale, const AMatType& A,
                                        const BMatType& B, CMatType& C) {
  C(0, 0) +=
      scale * (A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0) + A(0, 2) * B(2, 0));
  C(0, 1) +=
      scale * (A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1) + A(0, 2) * B(2, 1));
  C(0, 2) +=
      scale * (A(0, 0) * B(0, 2) + A(0, 1) * B(1, 2) + A(0, 2) * B(2, 2));
  C(1, 0) +=
      scale * (A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0) + A(1, 2) * B(2, 0));
  C(1, 1) +=
      scale * (A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1) + A(1, 2) * B(2, 1));
  C(1, 2) +=
      scale * (A(1, 0) * B(0, 2) + A(1, 1) * B(1, 2) + A(1, 2) * B(2, 2));
  C(2, 0) +=
      scale * (A(2, 0) * B(0, 0) + A(2, 1) * B(1, 0) + A(2, 2) * B(2, 0));
  C(2, 1) +=
      scale * (A(2, 0) * B(0, 1) + A(2, 1) * B(1, 1) + A(2, 2) * B(2, 1));
  C(2, 2) +=
      scale * (A(2, 0) * B(0, 2) + A(2, 1) * B(1, 2) + A(2, 2) * B(2, 2));
}

template <class AMatType, class BMatType, class CMatType>
inline void Symm3x3MatMultCore(const AMatType& A, const BMatType& B,
                               CMatType& C) {
  C(0, 0) = (A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0) + A(0, 2) * B(2, 0));
  C(0, 1) = (A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1) + A(0, 2) * B(2, 1));
  C(0, 2) = (A(0, 0) * B(0, 2) + A(0, 1) * B(1, 2) + A(0, 2) * B(2, 2));
  C(1, 0) = (A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0) + A(1, 2) * B(2, 0));
  C(1, 1) = (A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1) + A(1, 2) * B(2, 1));
  C(1, 2) = (A(1, 0) * B(0, 2) + A(1, 1) * B(1, 2) + A(1, 2) * B(2, 2));
  C(2, 0) = (A(2, 0) * B(0, 0) + A(2, 1) * B(1, 0) + A(2, 2) * B(2, 0));
  C(2, 1) = (A(2, 0) * B(0, 1) + A(2, 1) * B(1, 1) + A(2, 2) * B(2, 1));
  C(2, 2) = (A(2, 0) * B(0, 2) + A(2, 1) * B(1, 2) + A(2, 2) * B(2, 2));
}

template <typename T, class AMatType, class BMatType, class CMatType>
inline void Symm3x3MatMultScaleCore(T scale, const AMatType& A,
                                    const BMatType& B, CMatType& C) {
  C(0, 0) = scale * (A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0) + A(0, 2) * B(2, 0));
  C(0, 1) = scale * (A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1) + A(0, 2) * B(2, 1));
  C(0, 2) = scale * (A(0, 0) * B(0, 2) + A(0, 1) * B(1, 2) + A(0, 2) * B(2, 2));
  C(1, 0) = scale * (A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0) + A(1, 2) * B(2, 0));
  C(1, 1) = scale * (A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1) + A(1, 2) * B(2, 1));
  C(1, 2) = scale * (A(1, 0) * B(0, 2) + A(1, 1) * B(1, 2) + A(1, 2) * B(2, 2));
  C(2, 0) = scale * (A(2, 0) * B(0, 0) + A(2, 1) * B(1, 0) + A(2, 2) * B(2, 0));
  C(2, 1) = scale * (A(2, 0) * B(0, 1) + A(2, 1) * B(1, 1) + A(2, 2) * B(2, 1));
  C(2, 2) = scale * (A(2, 0) * B(0, 2) + A(2, 1) * B(1, 2) + A(2, 2) * B(2, 2));
}

template <class AMatType, class BMatType, class CMatType>
inline void Symm3x3MatMultAddCore(const AMatType& A, const BMatType& B,
                                  CMatType& C) {
  C(0, 0) += (A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0) + A(0, 2) * B(2, 0));
  C(0, 1) += (A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1) + A(0, 2) * B(2, 1));
  C(0, 2) += (A(0, 0) * B(0, 2) + A(0, 1) * B(1, 2) + A(0, 2) * B(2, 2));
  C(1, 0) += (A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0) + A(1, 2) * B(2, 0));
  C(1, 1) += (A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1) + A(1, 2) * B(2, 1));
  C(1, 2) += (A(1, 0) * B(0, 2) + A(1, 1) * B(1, 2) + A(1, 2) * B(2, 2));
  C(2, 0) += (A(2, 0) * B(0, 0) + A(2, 1) * B(1, 0) + A(2, 2) * B(2, 0));
  C(2, 1) += (A(2, 0) * B(0, 1) + A(2, 1) * B(1, 1) + A(2, 2) * B(2, 1));
  C(2, 2) += (A(2, 0) * B(0, 2) + A(2, 1) * B(1, 2) + A(2, 2) * B(2, 2));
}

template <class AMatType, class BMatType, class CMatType>
inline void Symm3x3MatMultSubCore(const AMatType& A, const BMatType& B,
                                  CMatType& C) {
  C(0, 0) -= (A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0) + A(0, 2) * B(2, 0));
  C(0, 1) -= (A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1) + A(0, 2) * B(2, 1));
  C(0, 2) -= (A(0, 0) * B(0, 2) + A(0, 1) * B(1, 2) + A(0, 2) * B(2, 2));
  C(1, 0) -= (A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0) + A(1, 2) * B(2, 0));
  C(1, 1) -= (A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1) + A(1, 2) * B(2, 1));
  C(1, 2) -= (A(1, 0) * B(0, 2) + A(1, 1) * B(1, 2) + A(1, 2) * B(2, 2));
  C(2, 0) -= (A(2, 0) * B(0, 0) + A(2, 1) * B(1, 0) + A(2, 2) * B(2, 0));
  C(2, 1) -= (A(2, 0) * B(0, 1) + A(2, 1) * B(1, 1) + A(2, 2) * B(2, 1));
  C(2, 2) -= (A(2, 0) * B(0, 2) + A(2, 1) * B(1, 2) + A(2, 2) * B(2, 2));
}

template <typename T, class AMatType, class BMatType, class CMatType>
inline void Symm3x3MatMultAddScaleCore(T scale, const AMatType& A,
                                       const BMatType& B, CMatType& C) {
  C(0, 0) +=
      scale * (A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0) + A(0, 2) * B(2, 0));
  C(0, 1) +=
      scale * (A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1) + A(0, 2) * B(2, 1));
  C(0, 2) +=
      scale * (A(0, 0) * B(0, 2) + A(0, 1) * B(1, 2) + A(0, 2) * B(2, 2));
  C(1, 0) +=
      scale * (A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0) + A(1, 2) * B(2, 0));
  C(1, 1) +=
      scale * (A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1) + A(1, 2) * B(2, 1));
  C(1, 2) +=
      scale * (A(1, 0) * B(0, 2) + A(1, 1) * B(1, 2) + A(1, 2) * B(2, 2));
  C(2, 0) +=
      scale * (A(2, 0) * B(0, 0) + A(2, 1) * B(1, 0) + A(2, 2) * B(2, 0));
  C(2, 1) +=
      scale * (A(2, 0) * B(0, 1) + A(2, 1) * B(1, 1) + A(2, 2) * B(2, 1));
  C(2, 2) +=
      scale * (A(2, 0) * B(0, 2) + A(2, 1) * B(1, 2) + A(2, 2) * B(2, 2));
}

template <class AMatType, class BMatType, class CMatType>
inline void Symm3x3MatTransMultCore(const AMatType& A, const BMatType& B,
                                    CMatType& C) {
  C(0, 0) = (A(0, 0) * B(0, 0) + A(0, 1) * B(0, 1) + A(0, 2) * B(0, 2));
  C(0, 1) = (A(0, 0) * B(1, 0) + A(0, 1) * B(1, 1) + A(0, 2) * B(1, 2));
  C(0, 2) = (A(0, 0) * B(2, 0) + A(0, 1) * B(2, 1) + A(0, 2) * B(2, 2));
  C(1, 0) = (A(1, 0) * B(0, 0) + A(1, 1) * B(0, 1) + A(1, 2) * B(0, 2));
  C(1, 1) = (A(1, 0) * B(1, 0) + A(1, 1) * B(1, 1) + A(1, 2) * B(1, 2));
  C(1, 2) = (A(1, 0) * B(2, 0) + A(1, 1) * B(2, 1) + A(1, 2) * B(2, 2));
  C(2, 0) = (A(2, 0) * B(0, 0) + A(2, 1) * B(0, 1) + A(2, 2) * B(0, 2));
  C(2, 1) = (A(2, 0) * B(1, 0) + A(2, 1) * B(1, 1) + A(2, 2) * B(1, 2));
  C(2, 2) = (A(2, 0) * B(2, 0) + A(2, 1) * B(2, 1) + A(2, 2) * B(2, 2));
}

template <typename T, class AMatType, class BMatType, class CMatType>
inline void Symm3x3MatTransMultScaleCore(T scale, const AMatType& A,
                                         const BMatType& B, CMatType& C) {
  C(0, 0) = scale * (A(0, 0) * B(0, 0) + A(0, 1) * B(0, 1) + A(0, 2) * B(0, 2));
  C(0, 1) = scale * (A(0, 0) * B(1, 0) + A(0, 1) * B(1, 1) + A(0, 2) * B(1, 2));
  C(0, 2) = scale * (A(0, 0) * B(2, 0) + A(0, 1) * B(2, 1) + A(0, 2) * B(2, 2));
  C(1, 0) = scale * (A(1, 0) * B(0, 0) + A(1, 1) * B(0, 1) + A(1, 2) * B(0, 2));
  C(1, 1) = scale * (A(1, 0) * B(1, 0) + A(1, 1) * B(1, 1) + A(1, 2) * B(1, 2));
  C(1, 2) = scale * (A(1, 0) * B(2, 0) + A(1, 1) * B(2, 1) + A(1, 2) * B(2, 2));
  C(2, 0) = scale * (A(2, 0) * B(0, 0) + A(2, 1) * B(0, 1) + A(2, 2) * B(0, 2));
  C(2, 1) = scale * (A(2, 0) * B(1, 0) + A(2, 1) * B(1, 1) + A(2, 2) * B(1, 2));
  C(2, 2) = scale * (A(2, 0) * B(2, 0) + A(2, 1) * B(2, 1) + A(2, 2) * B(2, 2));
}

template <class AMatType, class BMatType, class CMatType>
inline void Symm3x3MatTransMultAddCore(const AMatType& A, const BMatType& B,
                                       CMatType& C) {
  C(0, 0) += (A(0, 0) * B(0, 0) + A(0, 1) * B(0, 1) + A(0, 2) * B(0, 2));
  C(0, 1) += (A(0, 0) * B(1, 0) + A(0, 1) * B(1, 1) + A(0, 2) * B(1, 2));
  C(0, 2) += (A(0, 0) * B(2, 0) + A(0, 1) * B(2, 1) + A(0, 2) * B(2, 2));
  C(1, 0) += (A(1, 0) * B(0, 0) + A(1, 1) * B(0, 1) + A(1, 2) * B(0, 2));
  C(1, 1) += (A(1, 0) * B(1, 0) + A(1, 1) * B(1, 1) + A(1, 2) * B(1, 2));
  C(1, 2) += (A(1, 0) * B(2, 0) + A(1, 1) * B(2, 1) + A(1, 2) * B(2, 2));
  C(2, 0) += (A(2, 0) * B(0, 0) + A(2, 1) * B(0, 1) + A(2, 2) * B(0, 2));
  C(2, 1) += (A(2, 0) * B(1, 0) + A(2, 1) * B(1, 1) + A(2, 2) * B(1, 2));
  C(2, 2) += (A(2, 0) * B(2, 0) + A(2, 1) * B(2, 1) + A(2, 2) * B(2, 2));
}

template <class AMatType, class BMatType, class CMatType>
inline void Symm3x3MatTransMultSubCore(const AMatType& A, const BMatType& B,
                                       CMatType& C) {
  C(0, 0) -= (A(0, 0) * B(0, 0) + A(0, 1) * B(0, 1) + A(0, 2) * B(0, 2));
  C(0, 1) -= (A(0, 0) * B(1, 0) + A(0, 1) * B(1, 1) + A(0, 2) * B(1, 2));
  C(0, 2) -= (A(0, 0) * B(2, 0) + A(0, 1) * B(2, 1) + A(0, 2) * B(2, 2));
  C(1, 0) -= (A(1, 0) * B(0, 0) + A(1, 1) * B(0, 1) + A(1, 2) * B(0, 2));
  C(1, 1) -= (A(1, 0) * B(1, 0) + A(1, 1) * B(1, 1) + A(1, 2) * B(1, 2));
  C(1, 2) -= (A(1, 0) * B(2, 0) + A(1, 1) * B(2, 1) + A(1, 2) * B(2, 2));
  C(2, 0) -= (A(2, 0) * B(0, 0) + A(2, 1) * B(0, 1) + A(2, 2) * B(0, 2));
  C(2, 1) -= (A(2, 0) * B(1, 0) + A(2, 1) * B(1, 1) + A(2, 2) * B(1, 2));
  C(2, 2) -= (A(2, 0) * B(2, 0) + A(2, 1) * B(2, 1) + A(2, 2) * B(2, 2));
}

template <typename T, class AMatType, class BMatType, class CMatType>
inline void Symm3x3MatTransMultAddScaleCore(T scale, const AMatType& A,
                                            const BMatType& B, CMatType& C) {
  C(0, 0) +=
      scale * (A(0, 0) * B(0, 0) + A(0, 1) * B(0, 1) + A(0, 2) * B(0, 2));
  C(0, 1) +=
      scale * (A(0, 0) * B(1, 0) + A(0, 1) * B(1, 1) + A(0, 2) * B(1, 2));
  C(0, 2) +=
      scale * (A(0, 0) * B(2, 0) + A(0, 1) * B(2, 1) + A(0, 2) * B(2, 2));
  C(1, 0) +=
      scale * (A(1, 0) * B(0, 0) + A(1, 1) * B(0, 1) + A(1, 2) * B(0, 2));
  C(1, 1) +=
      scale * (A(1, 0) * B(1, 0) + A(1, 1) * B(1, 1) + A(1, 2) * B(1, 2));
  C(1, 2) +=
      scale * (A(1, 0) * B(2, 0) + A(1, 1) * B(2, 1) + A(1, 2) * B(2, 2));
  C(2, 0) +=
      scale * (A(2, 0) * B(0, 0) + A(2, 1) * B(0, 1) + A(2, 2) * B(0, 2));
  C(2, 1) +=
      scale * (A(2, 0) * B(1, 0) + A(2, 1) * B(1, 1) + A(2, 2) * B(1, 2));
  C(2, 2) +=
      scale * (A(2, 0) * B(2, 0) + A(2, 1) * B(2, 1) + A(2, 2) * B(2, 2));
}

template <class AMatType, class BMatType, class CMatType>
inline void Mat3x3SymmMultCore(const AMatType& A, const BMatType& B,
                               CMatType& C) {
  C(0, 0) = (A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0) + A(0, 2) * B(2, 0));
  C(0, 1) = (A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1) + A(0, 2) * B(2, 1));
  C(0, 2) = (A(0, 0) * B(0, 2) + A(0, 1) * B(1, 2) + A(0, 2) * B(2, 2));
  C(1, 0) = (A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0) + A(1, 2) * B(2, 0));
  C(1, 1) = (A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1) + A(1, 2) * B(2, 1));
  C(1, 2) = (A(1, 0) * B(0, 2) + A(1, 1) * B(1, 2) + A(1, 2) * B(2, 2));
  C(2, 0) = (A(2, 0) * B(0, 0) + A(2, 1) * B(1, 0) + A(2, 2) * B(2, 0));
  C(2, 1) = (A(2, 0) * B(0, 1) + A(2, 1) * B(1, 1) + A(2, 2) * B(2, 1));
  C(2, 2) = (A(2, 0) * B(0, 2) + A(2, 1) * B(1, 2) + A(2, 2) * B(2, 2));
}

template <typename T, class AMatType, class BMatType, class CMatType>
inline void Mat3x3SymmMultScaleCore(T scale, const AMatType& A,
                                    const BMatType& B, CMatType& C) {
  C(0, 0) = scale * (A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0) + A(0, 2) * B(2, 0));
  C(0, 1) = scale * (A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1) + A(0, 2) * B(2, 1));
  C(0, 2) = scale * (A(0, 0) * B(0, 2) + A(0, 1) * B(1, 2) + A(0, 2) * B(2, 2));
  C(1, 0) = scale * (A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0) + A(1, 2) * B(2, 0));
  C(1, 1) = scale * (A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1) + A(1, 2) * B(2, 1));
  C(1, 2) = scale * (A(1, 0) * B(0, 2) + A(1, 1) * B(1, 2) + A(1, 2) * B(2, 2));
  C(2, 0) = scale * (A(2, 0) * B(0, 0) + A(2, 1) * B(1, 0) + A(2, 2) * B(2, 0));
  C(2, 1) = scale * (A(2, 0) * B(0, 1) + A(2, 1) * B(1, 1) + A(2, 2) * B(2, 1));
  C(2, 2) = scale * (A(2, 0) * B(0, 2) + A(2, 1) * B(1, 2) + A(2, 2) * B(2, 2));
}

template <class AMatType, class BMatType, class CMatType>
inline void Mat3x3SymmMultAddCore(const AMatType& A, const BMatType& B,
                                  CMatType& C) {
  C(0, 0) += (A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0) + A(0, 2) * B(2, 0));
  C(0, 1) += (A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1) + A(0, 2) * B(2, 1));
  C(0, 2) += (A(0, 0) * B(0, 2) + A(0, 1) * B(1, 2) + A(0, 2) * B(2, 2));
  C(1, 0) += (A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0) + A(1, 2) * B(2, 0));
  C(1, 1) += (A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1) + A(1, 2) * B(2, 1));
  C(1, 2) += (A(1, 0) * B(0, 2) + A(1, 1) * B(1, 2) + A(1, 2) * B(2, 2));
  C(2, 0) += (A(2, 0) * B(0, 0) + A(2, 1) * B(1, 0) + A(2, 2) * B(2, 0));
  C(2, 1) += (A(2, 0) * B(0, 1) + A(2, 1) * B(1, 1) + A(2, 2) * B(2, 1));
  C(2, 2) += (A(2, 0) * B(0, 2) + A(2, 1) * B(1, 2) + A(2, 2) * B(2, 2));
}

template <class AMatType, class BMatType, class CMatType>
inline void Mat3x3SymmMultSubCore(const AMatType& A, const BMatType& B,
                                  CMatType& C) {
  C(0, 0) -= (A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0) + A(0, 2) * B(2, 0));
  C(0, 1) -= (A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1) + A(0, 2) * B(2, 1));
  C(0, 2) -= (A(0, 0) * B(0, 2) + A(0, 1) * B(1, 2) + A(0, 2) * B(2, 2));
  C(1, 0) -= (A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0) + A(1, 2) * B(2, 0));
  C(1, 1) -= (A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1) + A(1, 2) * B(2, 1));
  C(1, 2) -= (A(1, 0) * B(0, 2) + A(1, 1) * B(1, 2) + A(1, 2) * B(2, 2));
  C(2, 0) -= (A(2, 0) * B(0, 0) + A(2, 1) * B(1, 0) + A(2, 2) * B(2, 0));
  C(2, 1) -= (A(2, 0) * B(0, 1) + A(2, 1) * B(1, 1) + A(2, 2) * B(2, 1));
  C(2, 2) -= (A(2, 0) * B(0, 2) + A(2, 1) * B(1, 2) + A(2, 2) * B(2, 2));
}

template <typename T, class AMatType, class BMatType, class CMatType>
inline void Mat3x3SymmMultAddScaleCore(T scale, const AMatType& A,
                                       const BMatType& B, CMatType& C) {
  C(0, 0) +=
      scale * (A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0) + A(0, 2) * B(2, 0));
  C(0, 1) +=
      scale * (A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1) + A(0, 2) * B(2, 1));
  C(0, 2) +=
      scale * (A(0, 0) * B(0, 2) + A(0, 1) * B(1, 2) + A(0, 2) * B(2, 2));
  C(1, 0) +=
      scale * (A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0) + A(1, 2) * B(2, 0));
  C(1, 1) +=
      scale * (A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1) + A(1, 2) * B(2, 1));
  C(1, 2) +=
      scale * (A(1, 0) * B(0, 2) + A(1, 1) * B(1, 2) + A(1, 2) * B(2, 2));
  C(2, 0) +=
      scale * (A(2, 0) * B(0, 0) + A(2, 1) * B(1, 0) + A(2, 2) * B(2, 0));
  C(2, 1) +=
      scale * (A(2, 0) * B(0, 1) + A(2, 1) * B(1, 1) + A(2, 2) * B(2, 1));
  C(2, 2) +=
      scale * (A(2, 0) * B(0, 2) + A(2, 1) * B(1, 2) + A(2, 2) * B(2, 2));
}

template <class AMatType, class BMatType, class CMatType>
inline void MatTrans3x3SymmMultCore(const AMatType& A, const BMatType& B,
                                    CMatType& C) {
  C(0, 0) = (A(0, 0) * B(0, 0) + A(1, 0) * B(1, 0) + A(2, 0) * B(2, 0));
  C(0, 1) = (A(0, 0) * B(0, 1) + A(1, 0) * B(1, 1) + A(2, 0) * B(2, 1));
  C(0, 2) = (A(0, 0) * B(0, 2) + A(1, 0) * B(1, 2) + A(2, 0) * B(2, 2));
  C(1, 0) = (A(0, 1) * B(0, 0) + A(1, 1) * B(1, 0) + A(2, 1) * B(2, 0));
  C(1, 1) = (A(0, 1) * B(0, 1) + A(1, 1) * B(1, 1) + A(2, 1) * B(2, 1));
  C(1, 2) = (A(0, 1) * B(0, 2) + A(1, 1) * B(1, 2) + A(2, 1) * B(2, 2));
  C(2, 0) = (A(0, 2) * B(0, 0) + A(1, 2) * B(1, 0) + A(2, 2) * B(2, 0));
  C(2, 1) = (A(0, 2) * B(0, 1) + A(1, 2) * B(1, 1) + A(2, 2) * B(2, 1));
  C(2, 2) = (A(0, 2) * B(0, 2) + A(1, 2) * B(1, 2) + A(2, 2) * B(2, 2));
}

template <typename T, class AMatType, class BMatType, class CMatType>
inline void MatTrans3x3SymmMultScaleCore(T scale, const AMatType& A,
                                         const BMatType& B, CMatType& C) {
  C(0, 0) = scale * (A(0, 0) * B(0, 0) + A(1, 0) * B(1, 0) + A(2, 0) * B(2, 0));
  C(0, 1) = scale * (A(0, 0) * B(0, 1) + A(1, 0) * B(1, 1) + A(2, 0) * B(2, 1));
  C(0, 2) = scale * (A(0, 0) * B(0, 2) + A(1, 0) * B(1, 2) + A(2, 0) * B(2, 2));
  C(1, 0) = scale * (A(0, 1) * B(0, 0) + A(1, 1) * B(1, 0) + A(2, 1) * B(2, 0));
  C(1, 1) = scale * (A(0, 1) * B(0, 1) + A(1, 1) * B(1, 1) + A(2, 1) * B(2, 1));
  C(1, 2) = scale * (A(0, 1) * B(0, 2) + A(1, 1) * B(1, 2) + A(2, 1) * B(2, 2));
  C(2, 0) = scale * (A(0, 2) * B(0, 0) + A(1, 2) * B(1, 0) + A(2, 2) * B(2, 0));
  C(2, 1) = scale * (A(0, 2) * B(0, 1) + A(1, 2) * B(1, 1) + A(2, 2) * B(2, 1));
  C(2, 2) = scale * (A(0, 2) * B(0, 2) + A(1, 2) * B(1, 2) + A(2, 2) * B(2, 2));
}

template <class AMatType, class BMatType, class CMatType>
inline void MatTrans3x3SymmMultAddCore(const AMatType& A, const BMatType& B,
                                       CMatType& C) {
  C(0, 0) += (A(0, 0) * B(0, 0) + A(1, 0) * B(1, 0) + A(2, 0) * B(2, 0));
  C(0, 1) += (A(0, 0) * B(0, 1) + A(1, 0) * B(1, 1) + A(2, 0) * B(2, 1));
  C(0, 2) += (A(0, 0) * B(0, 2) + A(1, 0) * B(1, 2) + A(2, 0) * B(2, 2));
  C(1, 0) += (A(0, 1) * B(0, 0) + A(1, 1) * B(1, 0) + A(2, 1) * B(2, 0));
  C(1, 1) += (A(0, 1) * B(0, 1) + A(1, 1) * B(1, 1) + A(2, 1) * B(2, 1));
  C(1, 2) += (A(0, 1) * B(0, 2) + A(1, 1) * B(1, 2) + A(2, 1) * B(2, 2));
  C(2, 0) += (A(0, 2) * B(0, 0) + A(1, 2) * B(1, 0) + A(2, 2) * B(2, 0));
  C(2, 1) += (A(0, 2) * B(0, 1) + A(1, 2) * B(1, 1) + A(2, 2) * B(2, 1));
  C(2, 2) += (A(0, 2) * B(0, 2) + A(1, 2) * B(1, 2) + A(2, 2) * B(2, 2));
}

template <class AMatType, class BMatType, class CMatType>
inline void MatTrans3x3SymmMultSubCore(const AMatType& A, const BMatType& B,
                                       CMatType& C) {
  C(0, 0) -= (A(0, 0) * B(0, 0) + A(1, 0) * B(1, 0) + A(2, 0) * B(2, 0));
  C(0, 1) -= (A(0, 0) * B(0, 1) + A(1, 0) * B(1, 1) + A(2, 0) * B(2, 1));
  C(0, 2) -= (A(0, 0) * B(0, 2) + A(1, 0) * B(1, 2) + A(2, 0) * B(2, 2));
  C(1, 0) -= (A(0, 1) * B(0, 0) + A(1, 1) * B(1, 0) + A(2, 1) * B(2, 0));
  C(1, 1) -= (A(0, 1) * B(0, 1) + A(1, 1) * B(1, 1) + A(2, 1) * B(2, 1));
  C(1, 2) -= (A(0, 1) * B(0, 2) + A(1, 1) * B(1, 2) + A(2, 1) * B(2, 2));
  C(2, 0) -= (A(0, 2) * B(0, 0) + A(1, 2) * B(1, 0) + A(2, 2) * B(2, 0));
  C(2, 1) -= (A(0, 2) * B(0, 1) + A(1, 2) * B(1, 1) + A(2, 2) * B(2, 1));
  C(2, 2) -= (A(0, 2) * B(0, 2) + A(1, 2) * B(1, 2) + A(2, 2) * B(2, 2));
}

template <typename T, class AMatType, class BMatType, class CMatType>
inline void MatTrans3x3SymmMultAddScaleCore(T scale, const AMatType& A,
                                            const BMatType& B, CMatType& C) {
  C(0, 0) +=
      scale * (A(0, 0) * B(0, 0) + A(1, 0) * B(1, 0) + A(2, 0) * B(2, 0));
  C(0, 1) +=
      scale * (A(0, 0) * B(0, 1) + A(1, 0) * B(1, 1) + A(2, 0) * B(2, 1));
  C(0, 2) +=
      scale * (A(0, 0) * B(0, 2) + A(1, 0) * B(1, 2) + A(2, 0) * B(2, 2));
  C(1, 0) +=
      scale * (A(0, 1) * B(0, 0) + A(1, 1) * B(1, 0) + A(2, 1) * B(2, 0));
  C(1, 1) +=
      scale * (A(0, 1) * B(0, 1) + A(1, 1) * B(1, 1) + A(2, 1) * B(2, 1));
  C(1, 2) +=
      scale * (A(0, 1) * B(0, 2) + A(1, 1) * B(1, 2) + A(2, 1) * B(2, 2));
  C(2, 0) +=
      scale * (A(0, 2) * B(0, 0) + A(1, 2) * B(1, 0) + A(2, 2) * B(2, 0));
  C(2, 1) +=
      scale * (A(0, 2) * B(0, 1) + A(1, 2) * B(1, 1) + A(2, 2) * B(2, 1));
  C(2, 2) +=
      scale * (A(0, 2) * B(0, 2) + A(1, 2) * B(1, 2) + A(2, 2) * B(2, 2));
}

template <class AMatType, class BMatType, class CMatType>
inline void Mat3x3MatMultCore(const AMatType& A, const BMatType& B,
                              CMatType& C) {
  C(0, 0) = (A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0) + A(0, 2) * B(2, 0));
  C(0, 1) = (A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1) + A(0, 2) * B(2, 1));
  C(0, 2) = (A(0, 0) * B(0, 2) + A(0, 1) * B(1, 2) + A(0, 2) * B(2, 2));
  C(1, 0) = (A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0) + A(1, 2) * B(2, 0));
  C(1, 1) = (A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1) + A(1, 2) * B(2, 1));
  C(1, 2) = (A(1, 0) * B(0, 2) + A(1, 1) * B(1, 2) + A(1, 2) * B(2, 2));
  C(2, 0) = (A(2, 0) * B(0, 0) + A(2, 1) * B(1, 0) + A(2, 2) * B(2, 0));
  C(2, 1) = (A(2, 0) * B(0, 1) + A(2, 1) * B(1, 1) + A(2, 2) * B(2, 1));
  C(2, 2) = (A(2, 0) * B(0, 2) + A(2, 1) * B(1, 2) + A(2, 2) * B(2, 2));
}

template <typename T, class AMatType, class BMatType, class CMatType>
inline void Mat3x3MatMultScaleCore(T scale, const AMatType& A,
                                   const BMatType& B, CMatType& C) {
  C(0, 0) = scale * (A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0) + A(0, 2) * B(2, 0));
  C(0, 1) = scale * (A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1) + A(0, 2) * B(2, 1));
  C(0, 2) = scale * (A(0, 0) * B(0, 2) + A(0, 1) * B(1, 2) + A(0, 2) * B(2, 2));
  C(1, 0) = scale * (A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0) + A(1, 2) * B(2, 0));
  C(1, 1) = scale * (A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1) + A(1, 2) * B(2, 1));
  C(1, 2) = scale * (A(1, 0) * B(0, 2) + A(1, 1) * B(1, 2) + A(1, 2) * B(2, 2));
  C(2, 0) = scale * (A(2, 0) * B(0, 0) + A(2, 1) * B(1, 0) + A(2, 2) * B(2, 0));
  C(2, 1) = scale * (A(2, 0) * B(0, 1) + A(2, 1) * B(1, 1) + A(2, 2) * B(2, 1));
  C(2, 2) = scale * (A(2, 0) * B(0, 2) + A(2, 1) * B(1, 2) + A(2, 2) * B(2, 2));
}

template <class AMatType, class BMatType, class CMatType>
inline void Mat3x3MatMultAddCore(const AMatType& A, const BMatType& B,
                                 CMatType& C) {
  C(0, 0) += (A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0) + A(0, 2) * B(2, 0));
  C(0, 1) += (A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1) + A(0, 2) * B(2, 1));
  C(0, 2) += (A(0, 0) * B(0, 2) + A(0, 1) * B(1, 2) + A(0, 2) * B(2, 2));
  C(1, 0) += (A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0) + A(1, 2) * B(2, 0));
  C(1, 1) += (A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1) + A(1, 2) * B(2, 1));
  C(1, 2) += (A(1, 0) * B(0, 2) + A(1, 1) * B(1, 2) + A(1, 2) * B(2, 2));
  C(2, 0) += (A(2, 0) * B(0, 0) + A(2, 1) * B(1, 0) + A(2, 2) * B(2, 0));
  C(2, 1) += (A(2, 0) * B(0, 1) + A(2, 1) * B(1, 1) + A(2, 2) * B(2, 1));
  C(2, 2) += (A(2, 0) * B(0, 2) + A(2, 1) * B(1, 2) + A(2, 2) * B(2, 2));
}

template <class AMatType, class BMatType, class CMatType>
inline void Mat3x3MatMultSubCore(const AMatType& A, const BMatType& B,
                                 CMatType& C) {
  C(0, 0) -= (A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0) + A(0, 2) * B(2, 0));
  C(0, 1) -= (A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1) + A(0, 2) * B(2, 1));
  C(0, 2) -= (A(0, 0) * B(0, 2) + A(0, 1) * B(1, 2) + A(0, 2) * B(2, 2));
  C(1, 0) -= (A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0) + A(1, 2) * B(2, 0));
  C(1, 1) -= (A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1) + A(1, 2) * B(2, 1));
  C(1, 2) -= (A(1, 0) * B(0, 2) + A(1, 1) * B(1, 2) + A(1, 2) * B(2, 2));
  C(2, 0) -= (A(2, 0) * B(0, 0) + A(2, 1) * B(1, 0) + A(2, 2) * B(2, 0));
  C(2, 1) -= (A(2, 0) * B(0, 1) + A(2, 1) * B(1, 1) + A(2, 2) * B(2, 1));
  C(2, 2) -= (A(2, 0) * B(0, 2) + A(2, 1) * B(1, 2) + A(2, 2) * B(2, 2));
}

template <typename T, class AMatType, class BMatType, class CMatType>
inline void Mat3x3MatMultAddScaleCore(T scale, const AMatType& A,
                                      const BMatType& B, CMatType& C) {
  C(0, 0) +=
      scale * (A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0) + A(0, 2) * B(2, 0));
  C(0, 1) +=
      scale * (A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1) + A(0, 2) * B(2, 1));
  C(0, 2) +=
      scale * (A(0, 0) * B(0, 2) + A(0, 1) * B(1, 2) + A(0, 2) * B(2, 2));
  C(1, 0) +=
      scale * (A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0) + A(1, 2) * B(2, 0));
  C(1, 1) +=
      scale * (A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1) + A(1, 2) * B(2, 1));
  C(1, 2) +=
      scale * (A(1, 0) * B(0, 2) + A(1, 1) * B(1, 2) + A(1, 2) * B(2, 2));
  C(2, 0) +=
      scale * (A(2, 0) * B(0, 0) + A(2, 1) * B(1, 0) + A(2, 2) * B(2, 0));
  C(2, 1) +=
      scale * (A(2, 0) * B(0, 1) + A(2, 1) * B(1, 1) + A(2, 2) * B(2, 1));
  C(2, 2) +=
      scale * (A(2, 0) * B(0, 2) + A(2, 1) * B(1, 2) + A(2, 2) * B(2, 2));
}

template <class AMatType, class BMatType, class CMatType>
inline void Mat3x3MatTransMultCore(const AMatType& A, const BMatType& B,
                                   CMatType& C) {
  C(0, 0) = (A(0, 0) * B(0, 0) + A(0, 1) * B(0, 1) + A(0, 2) * B(0, 2));
  C(0, 1) = (A(0, 0) * B(1, 0) + A(0, 1) * B(1, 1) + A(0, 2) * B(1, 2));
  C(0, 2) = (A(0, 0) * B(2, 0) + A(0, 1) * B(2, 1) + A(0, 2) * B(2, 2));
  C(1, 0) = (A(1, 0) * B(0, 0) + A(1, 1) * B(0, 1) + A(1, 2) * B(0, 2));
  C(1, 1) = (A(1, 0) * B(1, 0) + A(1, 1) * B(1, 1) + A(1, 2) * B(1, 2));
  C(1, 2) = (A(1, 0) * B(2, 0) + A(1, 1) * B(2, 1) + A(1, 2) * B(2, 2));
  C(2, 0) = (A(2, 0) * B(0, 0) + A(2, 1) * B(0, 1) + A(2, 2) * B(0, 2));
  C(2, 1) = (A(2, 0) * B(1, 0) + A(2, 1) * B(1, 1) + A(2, 2) * B(1, 2));
  C(2, 2) = (A(2, 0) * B(2, 0) + A(2, 1) * B(2, 1) + A(2, 2) * B(2, 2));
}

template <typename T, class AMatType, class BMatType, class CMatType>
inline void Mat3x3MatTransMultScaleCore(T scale, const AMatType& A,
                                        const BMatType& B, CMatType& C) {
  C(0, 0) = scale * (A(0, 0) * B(0, 0) + A(0, 1) * B(0, 1) + A(0, 2) * B(0, 2));
  C(0, 1) = scale * (A(0, 0) * B(1, 0) + A(0, 1) * B(1, 1) + A(0, 2) * B(1, 2));
  C(0, 2) = scale * (A(0, 0) * B(2, 0) + A(0, 1) * B(2, 1) + A(0, 2) * B(2, 2));
  C(1, 0) = scale * (A(1, 0) * B(0, 0) + A(1, 1) * B(0, 1) + A(1, 2) * B(0, 2));
  C(1, 1) = scale * (A(1, 0) * B(1, 0) + A(1, 1) * B(1, 1) + A(1, 2) * B(1, 2));
  C(1, 2) = scale * (A(1, 0) * B(2, 0) + A(1, 1) * B(2, 1) + A(1, 2) * B(2, 2));
  C(2, 0) = scale * (A(2, 0) * B(0, 0) + A(2, 1) * B(0, 1) + A(2, 2) * B(0, 2));
  C(2, 1) = scale * (A(2, 0) * B(1, 0) + A(2, 1) * B(1, 1) + A(2, 2) * B(1, 2));
  C(2, 2) = scale * (A(2, 0) * B(2, 0) + A(2, 1) * B(2, 1) + A(2, 2) * B(2, 2));
}

template <class AMatType, class BMatType, class CMatType>
inline void Mat3x3MatTransMultAddCore(const AMatType& A, const BMatType& B,
                                      CMatType& C) {
  C(0, 0) += (A(0, 0) * B(0, 0) + A(0, 1) * B(0, 1) + A(0, 2) * B(0, 2));
  C(0, 1) += (A(0, 0) * B(1, 0) + A(0, 1) * B(1, 1) + A(0, 2) * B(1, 2));
  C(0, 2) += (A(0, 0) * B(2, 0) + A(0, 1) * B(2, 1) + A(0, 2) * B(2, 2));
  C(1, 0) += (A(1, 0) * B(0, 0) + A(1, 1) * B(0, 1) + A(1, 2) * B(0, 2));
  C(1, 1) += (A(1, 0) * B(1, 0) + A(1, 1) * B(1, 1) + A(1, 2) * B(1, 2));
  C(1, 2) += (A(1, 0) * B(2, 0) + A(1, 1) * B(2, 1) + A(1, 2) * B(2, 2));
  C(2, 0) += (A(2, 0) * B(0, 0) + A(2, 1) * B(0, 1) + A(2, 2) * B(0, 2));
  C(2, 1) += (A(2, 0) * B(1, 0) + A(2, 1) * B(1, 1) + A(2, 2) * B(1, 2));
  C(2, 2) += (A(2, 0) * B(2, 0) + A(2, 1) * B(2, 1) + A(2, 2) * B(2, 2));
}

template <class AMatType, class BMatType, class CMatType>
inline void Mat3x3MatTransMultSubCore(const AMatType& A, const BMatType& B,
                                      CMatType& C) {
  C(0, 0) -= (A(0, 0) * B(0, 0) + A(0, 1) * B(0, 1) + A(0, 2) * B(0, 2));
  C(0, 1) -= (A(0, 0) * B(1, 0) + A(0, 1) * B(1, 1) + A(0, 2) * B(1, 2));
  C(0, 2) -= (A(0, 0) * B(2, 0) + A(0, 1) * B(2, 1) + A(0, 2) * B(2, 2));
  C(1, 0) -= (A(1, 0) * B(0, 0) + A(1, 1) * B(0, 1) + A(1, 2) * B(0, 2));
  C(1, 1) -= (A(1, 0) * B(1, 0) + A(1, 1) * B(1, 1) + A(1, 2) * B(1, 2));
  C(1, 2) -= (A(1, 0) * B(2, 0) + A(1, 1) * B(2, 1) + A(1, 2) * B(2, 2));
  C(2, 0) -= (A(2, 0) * B(0, 0) + A(2, 1) * B(0, 1) + A(2, 2) * B(0, 2));
  C(2, 1) -= (A(2, 0) * B(1, 0) + A(2, 1) * B(1, 1) + A(2, 2) * B(1, 2));
  C(2, 2) -= (A(2, 0) * B(2, 0) + A(2, 1) * B(2, 1) + A(2, 2) * B(2, 2));
}

template <typename T, class AMatType, class BMatType, class CMatType>
inline void Mat3x3MatTransMultAddScaleCore(T scale, const AMatType& A,
                                           const BMatType& B, CMatType& C) {
  C(0, 0) +=
      scale * (A(0, 0) * B(0, 0) + A(0, 1) * B(0, 1) + A(0, 2) * B(0, 2));
  C(0, 1) +=
      scale * (A(0, 0) * B(1, 0) + A(0, 1) * B(1, 1) + A(0, 2) * B(1, 2));
  C(0, 2) +=
      scale * (A(0, 0) * B(2, 0) + A(0, 1) * B(2, 1) + A(0, 2) * B(2, 2));
  C(1, 0) +=
      scale * (A(1, 0) * B(0, 0) + A(1, 1) * B(0, 1) + A(1, 2) * B(0, 2));
  C(1, 1) +=
      scale * (A(1, 0) * B(1, 0) + A(1, 1) * B(1, 1) + A(1, 2) * B(1, 2));
  C(1, 2) +=
      scale * (A(1, 0) * B(2, 0) + A(1, 1) * B(2, 1) + A(1, 2) * B(2, 2));
  C(2, 0) +=
      scale * (A(2, 0) * B(0, 0) + A(2, 1) * B(0, 1) + A(2, 2) * B(0, 2));
  C(2, 1) +=
      scale * (A(2, 0) * B(1, 0) + A(2, 1) * B(1, 1) + A(2, 2) * B(1, 2));
  C(2, 2) +=
      scale * (A(2, 0) * B(2, 0) + A(2, 1) * B(2, 1) + A(2, 2) * B(2, 2));
}

template <class AMatType, class BMatType, class CMatType>
inline void MatTrans3x3MatMultCore(const AMatType& A, const BMatType& B,
                                   CMatType& C) {
  C(0, 0) = (A(0, 0) * B(0, 0) + A(1, 0) * B(1, 0) + A(2, 0) * B(2, 0));
  C(0, 1) = (A(0, 0) * B(0, 1) + A(1, 0) * B(1, 1) + A(2, 0) * B(2, 1));
  C(0, 2) = (A(0, 0) * B(0, 2) + A(1, 0) * B(1, 2) + A(2, 0) * B(2, 2));
  C(1, 0) = (A(0, 1) * B(0, 0) + A(1, 1) * B(1, 0) + A(2, 1) * B(2, 0));
  C(1, 1) = (A(0, 1) * B(0, 1) + A(1, 1) * B(1, 1) + A(2, 1) * B(2, 1));
  C(1, 2) = (A(0, 1) * B(0, 2) + A(1, 1) * B(1, 2) + A(2, 1) * B(2, 2));
  C(2, 0) = (A(0, 2) * B(0, 0) + A(1, 2) * B(1, 0) + A(2, 2) * B(2, 0));
  C(2, 1) = (A(0, 2) * B(0, 1) + A(1, 2) * B(1, 1) + A(2, 2) * B(2, 1));
  C(2, 2) = (A(0, 2) * B(0, 2) + A(1, 2) * B(1, 2) + A(2, 2) * B(2, 2));
}

template <typename T, class AMatType, class BMatType, class CMatType>
inline void MatTrans3x3MatMultScaleCore(T scale, const AMatType& A,
                                        const BMatType& B, CMatType& C) {
  C(0, 0) = scale * (A(0, 0) * B(0, 0) + A(1, 0) * B(1, 0) + A(2, 0) * B(2, 0));
  C(0, 1) = scale * (A(0, 0) * B(0, 1) + A(1, 0) * B(1, 1) + A(2, 0) * B(2, 1));
  C(0, 2) = scale * (A(0, 0) * B(0, 2) + A(1, 0) * B(1, 2) + A(2, 0) * B(2, 2));
  C(1, 0) = scale * (A(0, 1) * B(0, 0) + A(1, 1) * B(1, 0) + A(2, 1) * B(2, 0));
  C(1, 1) = scale * (A(0, 1) * B(0, 1) + A(1, 1) * B(1, 1) + A(2, 1) * B(2, 1));
  C(1, 2) = scale * (A(0, 1) * B(0, 2) + A(1, 1) * B(1, 2) + A(2, 1) * B(2, 2));
  C(2, 0) = scale * (A(0, 2) * B(0, 0) + A(1, 2) * B(1, 0) + A(2, 2) * B(2, 0));
  C(2, 1) = scale * (A(0, 2) * B(0, 1) + A(1, 2) * B(1, 1) + A(2, 2) * B(2, 1));
  C(2, 2) = scale * (A(0, 2) * B(0, 2) + A(1, 2) * B(1, 2) + A(2, 2) * B(2, 2));
}

template <class AMatType, class BMatType, class CMatType>
inline void MatTrans3x3MatMultAddCore(const AMatType& A, const BMatType& B,
                                      CMatType& C) {
  C(0, 0) += (A(0, 0) * B(0, 0) + A(1, 0) * B(1, 0) + A(2, 0) * B(2, 0));
  C(0, 1) += (A(0, 0) * B(0, 1) + A(1, 0) * B(1, 1) + A(2, 0) * B(2, 1));
  C(0, 2) += (A(0, 0) * B(0, 2) + A(1, 0) * B(1, 2) + A(2, 0) * B(2, 2));
  C(1, 0) += (A(0, 1) * B(0, 0) + A(1, 1) * B(1, 0) + A(2, 1) * B(2, 0));
  C(1, 1) += (A(0, 1) * B(0, 1) + A(1, 1) * B(1, 1) + A(2, 1) * B(2, 1));
  C(1, 2) += (A(0, 1) * B(0, 2) + A(1, 1) * B(1, 2) + A(2, 1) * B(2, 2));
  C(2, 0) += (A(0, 2) * B(0, 0) + A(1, 2) * B(1, 0) + A(2, 2) * B(2, 0));
  C(2, 1) += (A(0, 2) * B(0, 1) + A(1, 2) * B(1, 1) + A(2, 2) * B(2, 1));
  C(2, 2) += (A(0, 2) * B(0, 2) + A(1, 2) * B(1, 2) + A(2, 2) * B(2, 2));
}

template <class AMatType, class BMatType, class CMatType>
inline void MatTrans3x3MatMultSubCore(const AMatType& A, const BMatType& B,
                                      CMatType& C) {
  C(0, 0) -= (A(0, 0) * B(0, 0) + A(1, 0) * B(1, 0) + A(2, 0) * B(2, 0));
  C(0, 1) -= (A(0, 0) * B(0, 1) + A(1, 0) * B(1, 1) + A(2, 0) * B(2, 1));
  C(0, 2) -= (A(0, 0) * B(0, 2) + A(1, 0) * B(1, 2) + A(2, 0) * B(2, 2));
  C(1, 0) -= (A(0, 1) * B(0, 0) + A(1, 1) * B(1, 0) + A(2, 1) * B(2, 0));
  C(1, 1) -= (A(0, 1) * B(0, 1) + A(1, 1) * B(1, 1) + A(2, 1) * B(2, 1));
  C(1, 2) -= (A(0, 1) * B(0, 2) + A(1, 1) * B(1, 2) + A(2, 1) * B(2, 2));
  C(2, 0) -= (A(0, 2) * B(0, 0) + A(1, 2) * B(1, 0) + A(2, 2) * B(2, 0));
  C(2, 1) -= (A(0, 2) * B(0, 1) + A(1, 2) * B(1, 1) + A(2, 2) * B(2, 1));
  C(2, 2) -= (A(0, 2) * B(0, 2) + A(1, 2) * B(1, 2) + A(2, 2) * B(2, 2));
}

template <typename T, class AMatType, class BMatType, class CMatType>
inline void MatTrans3x3MatMultAddScaleCore(T scale, const AMatType& A,
                                           const BMatType& B, CMatType& C) {
  C(0, 0) +=
      scale * (A(0, 0) * B(0, 0) + A(1, 0) * B(1, 0) + A(2, 0) * B(2, 0));
  C(0, 1) +=
      scale * (A(0, 0) * B(0, 1) + A(1, 0) * B(1, 1) + A(2, 0) * B(2, 1));
  C(0, 2) +=
      scale * (A(0, 0) * B(0, 2) + A(1, 0) * B(1, 2) + A(2, 0) * B(2, 2));
  C(1, 0) +=
      scale * (A(0, 1) * B(0, 0) + A(1, 1) * B(1, 0) + A(2, 1) * B(2, 0));
  C(1, 1) +=
      scale * (A(0, 1) * B(0, 1) + A(1, 1) * B(1, 1) + A(2, 1) * B(2, 1));
  C(1, 2) +=
      scale * (A(0, 1) * B(0, 2) + A(1, 1) * B(1, 2) + A(2, 1) * B(2, 2));
  C(2, 0) +=
      scale * (A(0, 2) * B(0, 0) + A(1, 2) * B(1, 0) + A(2, 2) * B(2, 0));
  C(2, 1) +=
      scale * (A(0, 2) * B(0, 1) + A(1, 2) * B(1, 1) + A(2, 2) * B(2, 1));
  C(2, 2) +=
      scale * (A(0, 2) * B(0, 2) + A(1, 2) * B(1, 2) + A(2, 2) * B(2, 2));
}

template <class AMatType, class BMatType, class CMatType>
inline void MatTrans3x3MatTransMultCore(const AMatType& A, const BMatType& B,
                                        CMatType& C) {
  C(0, 0) = (A(0, 0) * B(0, 0) + A(1, 0) * B(0, 1) + A(2, 0) * B(0, 2));
  C(0, 1) = (A(0, 0) * B(1, 0) + A(1, 0) * B(1, 1) + A(2, 0) * B(1, 2));
  C(0, 2) = (A(0, 0) * B(2, 0) + A(1, 0) * B(2, 1) + A(2, 0) * B(2, 2));
  C(1, 0) = (A(0, 1) * B(0, 0) + A(1, 1) * B(0, 1) + A(2, 1) * B(0, 2));
  C(1, 1) = (A(0, 1) * B(1, 0) + A(1, 1) * B(1, 1) + A(2, 1) * B(1, 2));
  C(1, 2) = (A(0, 1) * B(2, 0) + A(1, 1) * B(2, 1) + A(2, 1) * B(2, 2));
  C(2, 0) = (A(0, 2) * B(0, 0) + A(1, 2) * B(0, 1) + A(2, 2) * B(0, 2));
  C(2, 1) = (A(0, 2) * B(1, 0) + A(1, 2) * B(1, 1) + A(2, 2) * B(1, 2));
  C(2, 2) = (A(0, 2) * B(2, 0) + A(1, 2) * B(2, 1) + A(2, 2) * B(2, 2));
}

template <typename T, class AMatType, class BMatType, class CMatType>
inline void MatTrans3x3MatTransMultScaleCore(T scale, const AMatType& A,
                                             const BMatType& B, CMatType& C) {
  C(0, 0) = scale * (A(0, 0) * B(0, 0) + A(1, 0) * B(0, 1) + A(2, 0) * B(0, 2));
  C(0, 1) = scale * (A(0, 0) * B(1, 0) + A(1, 0) * B(1, 1) + A(2, 0) * B(1, 2));
  C(0, 2) = scale * (A(0, 0) * B(2, 0) + A(1, 0) * B(2, 1) + A(2, 0) * B(2, 2));
  C(1, 0) = scale * (A(0, 1) * B(0, 0) + A(1, 1) * B(0, 1) + A(2, 1) * B(0, 2));
  C(1, 1) = scale * (A(0, 1) * B(1, 0) + A(1, 1) * B(1, 1) + A(2, 1) * B(1, 2));
  C(1, 2) = scale * (A(0, 1) * B(2, 0) + A(1, 1) * B(2, 1) + A(2, 1) * B(2, 2));
  C(2, 0) = scale * (A(0, 2) * B(0, 0) + A(1, 2) * B(0, 1) + A(2, 2) * B(0, 2));
  C(2, 1) = scale * (A(0, 2) * B(1, 0) + A(1, 2) * B(1, 1) + A(2, 2) * B(1, 2));
  C(2, 2) = scale * (A(0, 2) * B(2, 0) + A(1, 2) * B(2, 1) + A(2, 2) * B(2, 2));
}

template <class AMatType, class BMatType, class CMatType>
inline void MatTrans3x3MatTransMultAddCore(const AMatType& A, const BMatType& B,
                                           CMatType& C) {
  C(0, 0) += (A(0, 0) * B(0, 0) + A(1, 0) * B(0, 1) + A(2, 0) * B(0, 2));
  C(0, 1) += (A(0, 0) * B(1, 0) + A(1, 0) * B(1, 1) + A(2, 0) * B(1, 2));
  C(0, 2) += (A(0, 0) * B(2, 0) + A(1, 0) * B(2, 1) + A(2, 0) * B(2, 2));
  C(1, 0) += (A(0, 1) * B(0, 0) + A(1, 1) * B(0, 1) + A(2, 1) * B(0, 2));
  C(1, 1) += (A(0, 1) * B(1, 0) + A(1, 1) * B(1, 1) + A(2, 1) * B(1, 2));
  C(1, 2) += (A(0, 1) * B(2, 0) + A(1, 1) * B(2, 1) + A(2, 1) * B(2, 2));
  C(2, 0) += (A(0, 2) * B(0, 0) + A(1, 2) * B(0, 1) + A(2, 2) * B(0, 2));
  C(2, 1) += (A(0, 2) * B(1, 0) + A(1, 2) * B(1, 1) + A(2, 2) * B(1, 2));
  C(2, 2) += (A(0, 2) * B(2, 0) + A(1, 2) * B(2, 1) + A(2, 2) * B(2, 2));
}

template <class AMatType, class BMatType, class CMatType>
inline void MatTrans3x3MatTransMultSubCore(const AMatType& A, const BMatType& B,
                                           CMatType& C) {
  C(0, 0) -= (A(0, 0) * B(0, 0) + A(1, 0) * B(0, 1) + A(2, 0) * B(0, 2));
  C(0, 1) -= (A(0, 0) * B(1, 0) + A(1, 0) * B(1, 1) + A(2, 0) * B(1, 2));
  C(0, 2) -= (A(0, 0) * B(2, 0) + A(1, 0) * B(2, 1) + A(2, 0) * B(2, 2));
  C(1, 0) -= (A(0, 1) * B(0, 0) + A(1, 1) * B(0, 1) + A(2, 1) * B(0, 2));
  C(1, 1) -= (A(0, 1) * B(1, 0) + A(1, 1) * B(1, 1) + A(2, 1) * B(1, 2));
  C(1, 2) -= (A(0, 1) * B(2, 0) + A(1, 1) * B(2, 1) + A(2, 1) * B(2, 2));
  C(2, 0) -= (A(0, 2) * B(0, 0) + A(1, 2) * B(0, 1) + A(2, 2) * B(0, 2));
  C(2, 1) -= (A(0, 2) * B(1, 0) + A(1, 2) * B(1, 1) + A(2, 2) * B(1, 2));
  C(2, 2) -= (A(0, 2) * B(2, 0) + A(1, 2) * B(2, 1) + A(2, 2) * B(2, 2));
}

template <typename T, class AMatType, class BMatType, class CMatType>
inline void MatTrans3x3MatTransMultAddScaleCore(T scale, const AMatType& A,
                                                const BMatType& B,
                                                CMatType& C) {
  C(0, 0) +=
      scale * (A(0, 0) * B(0, 0) + A(1, 0) * B(0, 1) + A(2, 0) * B(0, 2));
  C(0, 1) +=
      scale * (A(0, 0) * B(1, 0) + A(1, 0) * B(1, 1) + A(2, 0) * B(1, 2));
  C(0, 2) +=
      scale * (A(0, 0) * B(2, 0) + A(1, 0) * B(2, 1) + A(2, 0) * B(2, 2));
  C(1, 0) +=
      scale * (A(0, 1) * B(0, 0) + A(1, 1) * B(0, 1) + A(2, 1) * B(0, 2));
  C(1, 1) +=
      scale * (A(0, 1) * B(1, 0) + A(1, 1) * B(1, 1) + A(2, 1) * B(1, 2));
  C(1, 2) +=
      scale * (A(0, 1) * B(2, 0) + A(1, 1) * B(2, 1) + A(2, 1) * B(2, 2));
  C(2, 0) +=
      scale * (A(0, 2) * B(0, 0) + A(1, 2) * B(0, 1) + A(2, 2) * B(0, 2));
  C(2, 1) +=
      scale * (A(0, 2) * B(1, 0) + A(1, 2) * B(1, 1) + A(2, 2) * B(1, 2));
  C(2, 2) +=
      scale * (A(0, 2) * B(2, 0) + A(1, 2) * B(2, 1) + A(2, 2) * B(2, 2));
}

}  // namespace A2D

#endif  // A2D_MAT_CORE_H
