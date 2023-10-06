#ifndef A2D_SYM_MAT_EIGS_H
#define A2D_SYM_MAT_EIGS_H

#include "a2ddefs.h"
#include "utils/a2dlapack.h"

namespace A2D {

template <typename T>
KOKKOS_FUNCTION void SymEigs2x2(const T* A, T* eigs, T* Q = nullptr) {
  T tr = A[0] + A[2];
  T diff = A[0] - A[2];
  T discrm = sqrt(diff * diff + 4.0 * A[1] * A[1]);
  T det = A[0] * A[2] - A[1] * A[1];

  // Compute the eigenvalues such that eigs[0] <= eigs[1]
  if (std::real(tr) > 0.0) {
    eigs[1] = 0.5 * (tr + discrm);
    eigs[0] = det / eigs[1];
  } else if (std::real(tr) < 0.0) {
    eigs[0] = 0.5 * (tr - discrm);
    eigs[1] = det / eigs[0];
  } else {
    eigs[0] = -0.5 * discrm;
    eigs[1] = 0.5 * discrm;
  }

  if (Q != nullptr) {
    // Compute the eigenvector components
    T u = 1.0, v = 0.0;
    if (std::real(A[1]) != 0.0) {
      if (std::real(diff) > 0.0) {
        T a = 0.5 * (diff + discrm);
        T inv = 1.0 / sqrt(a * a + A[1] * A[1]);
        u = inv * A[1];
        v = -inv * a;
      } else {
        T a = 0.5 * (diff - discrm);
        T inv = 1.0 / sqrt(a * a + A[1] * A[1]);
        u = inv * a;
        v = inv * A[1];
      }
    }

    // Set the eigenvector components
    Q[0] = u;
    Q[1] = -v;
    Q[2] = v;
    Q[3] = u;
  }
}

template <typename T, int N>
KOKKOS_FUNCTION void SymEigsForward(const T* eigs, const T* Q, const T* Ad,
                                    T* eigsd) {
  // eigsd[k] = Q[j, k]^{T} * Ad[j, i] * Q[i, k]
  for (int k = 0; k < N; k++) {
    T value = 0.0;
    for (int j = 0; j < N; j++) {
      const T* ad = &Ad[j * (j + 1) / 2];

      int i = 0;
      for (; i < j; i++) {
        value += ad[0] * Q[k + i * N] * Q[k + j * N];
        ad++;
      }
      for (; i < N; i++) {
        value += ad[0] * Q[k + i * N] * Q[k + j * N];
        ad += i + 1;
      }
    }
    eigsd[k] += value;
  }
}

template <typename T, int N>
KOKKOS_FUNCTION void SymEigsReverse(const T* eigs, const T* Q, const T* beigs,
                                    T* bA) {
  // bA[i, j] = Q[i, k] * beigs[k] * Q[j, k]
  for (int j = 0; j < N; j++) {
    for (int i = 0; i <= j; i++) {
      T value = 0.0;
      for (int k = 0; k < N; k++) {
        value += beigs[k] * Q[k + i * N] * Q[k + j * N];
      }

      if (i == j) {
        bA[0] += value;
      } else {
        bA[0] += 2.0 * value;
      }
      bA++;
    }
  }
}

// Compute the reverse mode for the second derivative of the eigenvalues
template <typename T, int N>
KOKKOS_FUNCTION void SymEigsHReverse(const T* eigs, const T* Q, const T* beigs,
                                     const T* Ap, T* Ah) {
  // Compute Bp = Q^{T} * Ap * Q
  // Bp[i, j] = Q[k, i] * Ap[k, l] * Q[l, j]
  T Bp[N * N];
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      T value = 0.0;
      for (int l = 0; l < N; l++) {
        const T* ap = &Ap[l * (l + 1) / 2];

        int k = 0;
        for (; k < l; k++) {
          value += ap[0] * Q[i + k * N] * Q[j + l * N];
          ap++;
        }
        for (; k < N; k++) {
          value += ap[0] * Q[i + k * N] * Q[j + l * N];
          ap += k + 1;
        }
      }

      Bp[j + i * N] = value;
    }
  }

  // Multiply by the F-matrix
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      // Bp[i, j] = 2.0 * beigs[i] / (eigs[i] - eigs[j])
      if (i == j || eigs[i] == eigs[j]) {
        Bp[j + i * N] = T(0.0);
      } else {
        Bp[j + i * N] *= 2.0 * beigs[i] / (eigs[i] - eigs[j]);
      }
    }
  }

  // Compute Ah += Q * Bp * Q^{T}
  // Ah[i, j] += Q[i, k] * Bp[k, l] * Q[j, l]
  for (int j = 0; j < N; j++) {
    for (int i = 0; i <= j; i++) {
      T value = 0.0;
      for (int k = 0; k < N; k++) {
        for (int l = 0; l < N; l++) {
          value += Bp[l + k * N] * Q[k + i * N] * Q[l + j * N];
        }
      }
      if (i != j) {
        for (int k = 0; k < N; k++) {
          for (int l = 0; l < N; l++) {
            value += Bp[l + k * N] * Q[k + j * N] * Q[l + i * N];
          }
        }
      }

      Ah[0] += value;
      Ah++;
    }
  }
}

/**
 * Compute the eigenvalues and eigenvectors of a symmetric eigenvalue problem
 */
template <typename T, int N>
KOKKOS_FUNCTION void SymEigs(const SymMat<T, N>& S, Vec<T, N>& eigs) {
  if constexpr (N == 2) {
    SymEigs2x2(get_data(S), get_data(eigs));
  } else {
    int n = N, ldz = N, info;
    T Ap[N * (N + 1) / 2], W[N];
    double work[3 * N];

    // Copy the data over since LAPACK over-writes the matrix
    for (int i = 0; i < N * (N + 1) / 2; i++) {
      Ap[i] = S[i];
    }
    LAPACKdspev("N", "U", &n, Ap, W, nullptr, &ldz, work, &info);

    // Set the eigenvalues
    for (int i = 0; i < N; i++) {
      eigs(i) = W[i];
    }
  }
}

template <class Stype, class etype>
class SymEigsExpr {
 public:
  // Extract the numeric type to use
  typedef typename get_object_numeric_type<etype>::type T;

  // Extract the dimensions of the underlying matrices
  static constexpr int N = get_symmatrix_size<Stype>::size;

  // Make sure we have the correct size
  static constexpr int K = get_vec_size<etype>::size;
  static_assert(K == N, "Vector of eigenvalues must be correct size");

  KOKKOS_FUNCTION SymEigsExpr(Stype& S, etype& eigs) : S(S), eigs(eigs) {}

  KOKKOS_FUNCTION void eval() {
    if constexpr (N == 2) {
      SymEigs2x2(get_data(S), get_data(eigs), get_data(Q));
    } else {
      auto S0 = S.value();
      auto eigs0 = eigs.value();

      int n = N, ldz = N, info;
      double Ap[N * (N + 1) / 2], Z[N * N], W[N];
      double work[3 * N];

      // Copy the data over since LAPACK over-writes the matrix
      for (int i = 0; i < N * (N + 1) / 2; i++) {
        Ap[i] = std::real(S0[i]);
      }
      LAPACKdspev("V", "U", &n, Ap, W, Z, &ldz, work, &info);

      // Copy the data back to the matrix
      for (int j = 0; j < N; j++) {
        for (int i = 0; i < N; i++) {
          Q(i, j) = T(Z[i + j * N]);
        }
      }

      // Set the eigenvalues
      for (int i = 0; i < N; i++) {
        eigs0(i) = T(W[i]);
      }
    }
  }

  KOKKOS_FUNCTION void bzero() { eigs.bzero(); }

  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;
    SymEigsForward<T, N>(get_data(eigs), get_data(Q),
                         GetSeed<seed>::get_data(S),
                         GetSeed<seed>::get_data(eigs));
  }

  KOKKOS_FUNCTION void reverse() {
    constexpr ADseed seed = ADseed::b;
    SymEigsReverse<T, N>(get_data(eigs), get_data(Q),
                         GetSeed<seed>::get_data(eigs),
                         GetSeed<seed>::get_data(S));
  }

  KOKKOS_FUNCTION void hzero() { eigs.hzero(); }

  KOKKOS_FUNCTION void hreverse() {
    SymEigsReverse<T, N>(get_data(eigs), get_data(Q),
                         GetSeed<ADseed::h>::get_data(eigs),
                         GetSeed<ADseed::h>::get_data(S));
    SymEigsHReverse<T, N>(
        get_data(eigs), get_data(Q), GetSeed<ADseed::b>::get_data(eigs),
        GetSeed<ADseed::p>::get_data(S), GetSeed<ADseed::h>::get_data(S));
  }

 private:
  Stype& S;
  etype& eigs;
  Mat<T, N, N> Q;
};

template <class Stype, class etype>
KOKKOS_FUNCTION auto SymEigs(ADObj<Stype>& S, ADObj<etype>& eigs) {
  return SymEigsExpr<ADObj<Stype>, ADObj<etype>>(S, eigs);
}

template <class Stype, class etype>
KOKKOS_FUNCTION auto SymEigs(A2DObj<Stype>& S, A2DObj<etype>& eigs) {
  return SymEigsExpr<A2DObj<Stype>, A2DObj<etype>>(S, eigs);
}

namespace Test {
template <typename T, int N>
class SymEigsTest : public A2DTest<T, Vec<T, N>, SymMat<T, N>> {
 public:
  using Input = VarTuple<T, SymMat<T, N>>;
  using Output = VarTuple<T, Vec<T, N>>;

  // Assemble a string to describe the test
  std::string name() {
    std::stringstream s;
    s << "SymEigsTest<" << N << ">";
    return s.str();
  }

  // Evaluate the matrix-matrix product
  Output eval(const Input& x) {
    SymMat<T, N> S;
    Vec<T, N> eigs;
    x.get_values(S);
    SymEigs(S, eigs);
    return MakeVarTuple<T>(eigs);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& x, Input& g) {
    ADObj<SymMat<T, N>> S;
    ADObj<Vec<T, N>> eigs;
    x.get_values(S.value());
    auto stack = MakeStack(SymEigs(S, eigs));
    seed.get_values(eigs.bvalue());
    stack.reverse();
    g.set_values(S.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) {
    A2DObj<SymMat<T, N>> S;
    A2DObj<Vec<T, N>> eigs;
    x.get_values(S.value());
    p.get_values(S.pvalue());
    auto stack = MakeStack(SymEigs(S, eigs));
    seed.get_values(eigs.bvalue());
    hval.get_values(eigs.hvalue());
    stack.hproduct();
    h.set_values(S.hvalue());
  }
};

bool SymEigsTestAll(bool component = false, bool write_output = true) {
  using Tc = std::complex<double>;

  bool passed = true;
  for (int i = 0; i < 10; i++) {
    SymEigsTest<Tc, 2> test1;
    passed = passed && Run(test1, component, write_output);
  }

  // for (int i = 0; i < 10; i++) {
  //   SymEigsTest<Tc, 3> test1;
  //   passed = passed && Run(test1, component, write_output);
  // }

  return passed;
}

}  // namespace Test

}  // namespace A2D

#endif  // A2D_SYM_MAT_EIGS_H
