#ifndef A2D_LAPACK_REAL_H
#define A2D_LAPACK_REAL_H
#include <complex>

namespace A2D {

// Real BLAS/LAPACK functions
extern "C" {
extern double ddot_(int *n, double *x, int *incx, double *y, int *incy);
extern double dnrm2_(int *n, double *x, int *incx);
extern void daxpy_(int *n, double *a, double *x, int *incx, double *y,
                   int *incy);
extern void dscal_(int *n, double *a, double *x, int *incx);
extern void dsyrk_(const char *uplo, const char *trans, int *n, int *k,
                   double *alpha, double *a, int *lda, double *beta, double *c,
                   int *ldc);
extern void dtpsv_(const char *uplo, const char *transa, const char *diag,
                   int *n, double *a, double *x, int *incx);
extern void dgemv_(const char *c, int *m, int *n, double *alpha, double *a,
                   int *lda, double *x, int *incx, double *beta, double *y,
                   int *incy);
extern void dgemm_(const char *ta, const char *tb, int *m, int *n, int *k,
                   double *alpha, double *a, int *lda, double *b, int *ldb,
                   double *beta, double *c, int *ldc);
extern void dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);
extern void dgetrs_(const char *c, int *n, int *nrhs, double *a, int *lda,
                    int *ipiv, double *b, int *ldb, int *info);
extern void dpptrf_(const char *c, int *n, double *ap, int *info);
extern void dpptrs_(const char *c, int *n, int *nrhs, double *ap, double *rhs,
                    int *ldrhs, int *info);
}

// Complex BLAS/LAPACK functions
extern "C" {
extern std::complex<double> zdotu_(int *n, std::complex<double> *x, int *incx,
                                   std::complex<double> *y, int *incy);
extern double dznrm2_(int *n, std::complex<double> *x, int *incx);
extern void zaxpy_(int *n, std::complex<double> *a, std::complex<double> *x,
                   int *incx, std::complex<double> *y, int *incy);
extern void zscal_(int *n, std::complex<double> *a, std::complex<double> *x,
                   int *incx);
extern void zsyrk_(const char *uplo, const char *trans, int *n, int *k,
                   std::complex<double> *alpha, std::complex<double> *a,
                   int *lda, std::complex<double> *beta,
                   std::complex<double> *c, int *ldc);
extern void ztpsv_(const char *uplo, const char *transa, const char *diag,
                   int *n, std::complex<double> *a, std::complex<double> *x,
                   int *incx);
extern void zgemv_(const char *c, int *m, int *n, std::complex<double> *alpha,
                   std::complex<double> *a, int *lda, std::complex<double> *x,
                   int *incx, std::complex<double> *beta,
                   std::complex<double> *y, int *incy);
extern void zgemm_(const char *ta, const char *tb, int *m, int *n, int *k,
                   std::complex<double> *alpha, std::complex<double> *a,
                   int *lda, std::complex<double> *b, int *ldb,
                   std::complex<double> *beta, std::complex<double> *c,
                   int *ldc);
extern void zgetrf_(int *m, int *n, std::complex<double> *a, int *lda,
                    int *ipiv, int *info);
extern void zgetrs_(const char *c, int *n, int *nrhs, std::complex<double> *a,
                    int *lda, int *ipiv, std::complex<double> *b, int *ldb,
                    int *info);
extern void zpptrf_(const char *c, int *n, std::complex<double> *ap, int *info);
extern void zpptrs_(const char *c, int *n, int *nrhs, std::complex<double> *ap,
                    std::complex<double> *rhs, int *ldrhs, int *info);
}

/* BLAS/LAPACK functions for double precision real datatype */

// Level 1 BLAS routines

inline double BLASdot(int *n, double *x, int *incx, double *y, int *incy) {
  return ddot_(n, x, incx, y, incy);
}
inline double BLASnrm2(int *n, double *x, int *incx) {
  return dnrm2_(n, x, incx);
}
inline void BLASaxpy(int *n, double *a, double *x, int *incx, double *y,
                     int *incy) {
  return daxpy_(n, a, x, incx, y, incy);
}
inline void BLASscal(int *n, double *a, double *x, int *incx) {
  return dscal_(n, a, x, incx);
}

// Compute C := alpha*A*A**double + beta*C or C := alpha*A**double*A + beta*C
inline void BLASsyrk(const char *uplo, const char *trans, int *n, int *k,
                     double *alpha, double *a, int *lda, double *beta,
                     double *c, int *ldc) {
  return dsyrk_(uplo, trans, n, k, alpha, a, lda, beta, c, ldc);
}

// Solve A*x = b or A^double*x = b where A is in packed format
inline void BLAStpsv(const char *uplo, const char *transa, const char *diag,
                     int *n, double *a, double *x, int *incx) {
  return dtpsv_(uplo, transa, diag, n, a, x, incx);
}

// Level 2 BLAS routines
// y = alpha * A * x + beta * y, for a general matrix
inline void BLASgemv(const char *c, int *m, int *n, double *alpha, double *a,
                     int *lda, double *x, int *incx, double *beta, double *y,
                     int *incy) {
  return dgemv_(c, m, n, alpha, a, lda, x, incx, beta, y, incy);
}

// Level 3 BLAS routines
// C := alpha*op( A )*op( B ) + beta*C,
inline void BLASgemm(const char *ta, const char *tb, int *m, int *n, int *k,
                     double *alpha, double *a, int *lda, double *b, int *ldb,
                     double *beta, double *c, int *ldc) {
  return dgemm_(ta, tb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

// General factorization routines
inline void LAPACKgetrf(int *m, int *n, double *a, int *lda, int *ipiv,
                        int *info) {
  return dgetrf_(m, n, a, lda, ipiv, info);
}
inline void LAPACKgetrs(const char *c, int *n, int *nrhs, double *a, int *lda,
                        int *ipiv, double *b, int *ldb, int *info) {
  return dgetrs_(c, n, nrhs, a, lda, ipiv, b, ldb, info);
}

// Factorization of packed-storage matrices
inline void LAPACKpptrf(const char *c, int *n, double *ap, int *info) {
  return dpptrf_(c, n, ap, info);
}
inline void LAPACKpptrs(const char *c, int *n, int *nrhs, double *ap,
                        double *rhs, int *ldrhs, int *info) {
  return dpptrs_(c, n, nrhs, ap, rhs, ldrhs, info);
}

/* BLAS/LAPACK functions for double precision complex datatype */

// Level 1 BLAS routines
inline std::complex<double> BLASdot(int *n, std::complex<double> *x, int *incx,
                                    std::complex<double> *y, int *incy) {
  return zdotu_(n, x, incx, y, incy);
}
inline double BLASnrm2(int *n, std::complex<double> *x, int *incx) {
  return dznrm2_(n, x, incx);
}
inline void BLASaxpy(int *n, std::complex<double> *a, std::complex<double> *x,
                     int *incx, std::complex<double> *y, int *incy) {
  return zaxpy_(n, a, x, incx, y, incy);
}
inline void BLASscal(int *n, std::complex<double> *a, std::complex<double> *x,
                     int *incx) {
  return zscal_(n, a, x, incx);
}

// Compute C := alpha*A*A**std::complex<double> + beta*C or C :=
// alpha*A**std::complex<double>*A + beta*C
inline void BLASsyrk(const char *uplo, const char *trans, int *n, int *k,
                     std::complex<double> *alpha, std::complex<double> *a,
                     int *lda, std::complex<double> *beta,
                     std::complex<double> *c, int *ldc) {
  return zsyrk_(uplo, trans, n, k, alpha, a, lda, beta, c, ldc);
}

// Solve A*x = b or A^std::complex<double>*x = b where A is in packed format
inline void BLAStpsv(const char *uplo, const char *transa, const char *diag,
                     int *n, std::complex<double> *a, std::complex<double> *x,
                     int *incx) {
  return ztpsv_(uplo, transa, diag, n, a, x, incx);
}

// Level 2 BLAS routines
// y = alpha * A * x + beta * y, for a general matrix
inline void BLASgemv(const char *c, int *m, int *n, std::complex<double> *alpha,
                     std::complex<double> *a, int *lda, std::complex<double> *x,
                     int *incx, std::complex<double> *beta,
                     std::complex<double> *y, int *incy) {
  return zgemv_(c, m, n, alpha, a, lda, x, incx, beta, y, incy);
}

// Level 3 BLAS routines
// C := alpha*op( A )*op( B ) + beta*C,
inline void BLASgemm(const char *ta, const char *tb, int *m, int *n, int *k,
                     std::complex<double> *alpha, std::complex<double> *a,
                     int *lda, std::complex<double> *b, int *ldb,
                     std::complex<double> *beta, std::complex<double> *c,
                     int *ldc) {
  return zgemm_(ta, tb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

// General factorization routines
inline void LAPACKgetrf(int *m, int *n, std::complex<double> *a, int *lda,
                        int *ipiv, int *info) {
  return zgetrf_(m, n, a, lda, ipiv, info);
}
inline void LAPACKgetrs(const char *c, int *n, int *nrhs,
                        std::complex<double> *a, int *lda, int *ipiv,
                        std::complex<double> *b, int *ldb, int *info) {
  return zgetrs_(c, n, nrhs, a, lda, ipiv, b, ldb, info);
}

// Factorization of packed-storage matrices
inline void LAPACKpptrf(const char *c, int *n, std::complex<double> *ap,
                        int *info) {
  return zpptrf_(c, n, ap, info);
}
inline void LAPACKpptrs(const char *c, int *n, int *nrhs,
                        std::complex<double> *ap, std::complex<double> *rhs,
                        int *ldrhs, int *info) {
  return zpptrs_(c, n, nrhs, ap, rhs, ldrhs, info);
}

}  // namespace A2D

#endif