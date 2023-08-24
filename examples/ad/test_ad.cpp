#include <complex>
#include <cstdlib>
#include <iostream>
#include <string>

#include "ad/a2dgemm.h"
#include "ad/a2dstack.h"
#include "ad/core/a2dgemmcore.h"

template <typename I, typename T>
void print_row_major_matrix(const std::string name, I M, I N, T mat[]) {
  std::cout << name << ":\n";
  for (I i = 0; i < M; i++) {
    for (I j = 0; j < N; j++) {
      std::cout << std::setw(15) << mat[i * N + j];
    }
    std::cout << "\n";
  }
}

// A set of reference implementations for the purpose of unit testing
namespace RefImpl {
// Transpose A.
// On entry, A is an M-by-N matrix, on exit, A is an N-by-M matrix
template <typename T, int M, int N>
void Transpose(T A[]) {
  T tmp[M * N];

  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) {
      tmp[M * j + i] = A[N * i + j];
    }
  }

  for (int i = 0; i < M * N; i++) {
    A[i] = tmp[i];
  }
}

// Computes C {=,+=} alpha * A * B, where A is M-by-P, B is P-by-N, C is
// M-by-N
template <typename T, int M, int P, int N, bool additive = false>
void MatMatMult(const T A[], const T B[], T C[], T alpha = T(1.0)) {
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) {
      if constexpr (!additive) {
        C[N * i + j] = T(0.0);
      }
      for (int p = 0; p < P; p++) {
        // C(i, j) = A(i, p) * B(p, j)
        C[N * i + j] += alpha * A[P * i + p] * B[N * p + j];
      }
    }
  }
}

/**
 * @brief Given F = F(A) and Adot, compute Fdot = dF/dA * Adot using complex
 * step
 *
 * Note: F and A are all matrices, eval has the following signature:
 *
 *   eval(const std::complex<T> Ac[], std::complex<T> Fc[])
 *
 * Note: On exit, Fdot is incremented, hence this function can be called
 * multiple times to evaluate forward derivatives using complex step of a
 * function that has multiple inputs: F(A, B, C, ...)
 */
template <typename T, int Am, int An, int Fm, int Fn, class Functor>
void complex_step(const T A[], const T Adot[], T Fdot[], Functor& eval,
                  T h = T(1e-30)) {
  std::complex<T> Ac[Am * An];    // A in complex
  std::complex<T> Ap[Am * An];    // A perturbed
  std::complex<T> dFdA[Fm * Fn];  // df/dA(k, l)

  for (int n = 0; n < Am * An; n++) {
    Ac[n] = A[n];
  }

  for (int k = 0; k < Am; k++) {
    for (int l = 0; l < An; l++) {
      // Reset A perturbed
      for (int n = 0; n < Am * An; n++) {
        Ap[n] = A[n];
      }

      // Perturb A(k, l) only
      Ap[k * An + l] += std::complex<T>(0.0, h);

      // Evaluate dF(i, j)/dA(k, l) for all i, j
      eval(Ap, dFdA);

      // Compute Fdot(i, j) += dF(i, j)/dA(k, l) * Adot(k, l) for all i, j
      for (int i = 0; i < Fm; i++) {
        for (int j = 0; j < Fn; j++) {
          Fdot[i * Fn + j] += dFdA[i * Fn + j].imag() / h * Adot[k * An + l];
        }
      }
    }
  }
}

}  // namespace RefImpl

template <int M, int N, int P>
void test_matmatmult_forward() {
  using T = double;

  std::srand(30067);  // Set rng seed

  int constexpr Am = M, An = P;
  int constexpr Bm = P, Bn = N;
  int constexpr Cm = M, Cn = N;

  A2D::Mat<T, Am, An> A, Ab;
  A2D::Mat<T, Bm, Bn> B, Bb;
  A2D::Mat<T, Cm, Cn> C, Cb;

  for (int i = 0; i < Am * An; i++) {
    A.data()[i] = T(std::rand()) / RAND_MAX;
    Ab.data()[i] = T(std::rand()) / RAND_MAX;
  }

  for (int i = 0; i < Bm * Bn; i++) {
    B.data()[i] = T(std::rand()) / RAND_MAX;
    Bb.data()[i] = T(std::rand()) / RAND_MAX;
  }

  A2D::Mat<T, Cm, Cn> Cref, Cbref;

  RefImpl::MatMatMult<T, Am, An, Bn>(A.data(), B.data(), Cref.data());

  A2D::ADMat<A2D::Mat<T, Am, An>> Aobj(A, Ab);
  A2D::ADMat<A2D::Mat<T, Bm, Bn>> Bobj(B, Bb);
  A2D::ADMat<A2D::Mat<T, Cm, Cn>> Cobj(C, Cb);

  auto expr = A2D::MatMatMult(Aobj, Bobj, Cobj);

  expr.template forward<A2D::ADorder::FIRST>();
  print_row_major_matrix("Cb", Cm, Cn, Cb.data());

  auto evalA = [=](const std::complex<T> Ac[], std::complex<T> Cc[]) mutable {
    std::complex<T> Bc[Bm * Bn];
    for (int n = 0; n < Bm * Bn; n++) {
      Bc[n] = B.data()[n];
    }
    A2D::MatMatMultCore<std::complex<T>, Am, An, Bm, Bn, Cm, Cn>(Ac, Bc, Cc);
  };

  auto evalB = [=](const std::complex<T> Bc[], std::complex<T> Cc[]) mutable {
    std::complex<T> Ac[Am * An];
    for (int n = 0; n < Am * An; n++) {
      Ac[n] = A.data()[n];
    }
    A2D::MatMatMultCore<std::complex<T>, Am, An, Bm, Bn, Cm, Cn>(Ac, Bc, Cc);
  };

  RefImpl::complex_step<T, Am, An, Cm, Cn>(A.data(), Ab.data(), Cbref.data(),
                                           evalA);
  RefImpl::complex_step<T, Bm, Bn, Cm, Cn>(B.data(), Bb.data(), Cbref.data(),
                                           evalB);

  print_row_major_matrix("Cb_ref", Cm, Cn, Cbref.data());
}

void test_matmatmult_op() {
  using T = double;

  A2D::Mat<T, 3, 3> A, B, C, D, E, F, G, Ab, Bb, Cb, Db, Eb, Fb, Gb, Ap, Bp, Cp,
      Dp, Ep, Fp, Gp, Ah, Bh, Ch, Dh, Eh, Fh, Gh;
  A2D::ADMat<A2D::Mat<T, 3, 3>> Amat(A, Ab), Bmat(B, Bb), Cmat(C, Cb),
      Dmat(D, Db), Emat(E, Eb), Fmat(F, Fb), Gmat(G, Gb);

  A2D::A2DMat<A2D::Mat<T, 3, 3>> A2Dmat(A, Ab, Ap, Ah), B2Dmat(B, Bb, Bp, Bh),
      C2Dmat(C, Cb, Cp, Ch);

  // Test all active/passive combinations
  A2D::MatMatMult(Amat, Bmat, Cmat);
  A2D::MatMatMult(A, Bmat, Cmat);
  A2D::MatMatMult(A, B, C);
  A2D::MatMatMult(Amat, B, Cmat);

  auto expr1 = A2D::MatMatMult(Amat, Bmat, Cmat);  // C = A * B
  auto expr2 = A2D::MatMatMult(Cmat, Dmat, Emat);  // E = C * D
  auto expr3 = A2D::MatMatMult(Emat, Fmat, Gmat);  // G = E * F
  auto stack = A2D::MakeStack(expr1, expr2, expr3);

  stack.reverse();
  stack.forward();

  auto exprH = A2D::MatMatMult(A2Dmat, B2Dmat, C2Dmat);
  exprH.forward<A2D::ADorder::FIRST>();
  exprH.reverse();
  exprH.forward<A2D::ADorder::SECOND>();
  exprH.hreverse();
}

int main() {
  test_matmatmult_forward<1, 1, 1>();
  test_matmatmult_forward<2, 2, 2>();
  test_matmatmult_forward<3, 3, 3>();
  test_matmatmult_forward<3, 4, 5>();
  test_matmatmult_op();
  return 0;
}