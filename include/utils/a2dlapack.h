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

/* Function templates - never call these */

// Level 1 BLAS routines
template <typename T>
T BLASddot(int *n, T *x, int *incx, T *y, int *incy) = delete;
template <typename T>
double BLASdnrm2(int *n, T *x, int *incx) = delete;
template <typename T>
void BLASdaxpy(int *n, T *a, T *x, int *incx, T *y, int *incy) = delete;
template <typename T>
void BLASdscal(int *n, T *a, T *x, int *incx) = delete;

// Compute C := alpha*A*A**T + beta*C or C := alpha*A**T*A + beta*C
template <typename T>
void BLASsyrk(const char *uplo, const char *trans, int *n, int *k, T *alpha,
              T *a, int *lda, T *beta, T *c, int *ldc) = delete;

// Solve A*x = b or A^T*x = b where A is in packed format
template <typename T>
void BLAStpsv(const char *uplo, const char *transa, const char *diag, int *n,
              T *a, T *x, int *incx) = delete;

// Level 2 BLAS routines
// y = alpha * A * x + beta * y, for a general matrix
template <typename T>
void BLASgemv(const char *c, int *m, int *n, T *alpha, T *a, int *lda, T *x,
              int *incx, T *beta, T *y, int *incy) = delete;

// Level 3 BLAS routines
// C := alpha*op( A )*op( B ) + beta*C,
template <typename T>
void BLASgemm(const char *ta, const char *tb, int *m, int *n, int *k, T *alpha,
              T *a, int *lda, T *b, int *ldb, T *beta, T *c, int *ldc) = delete;

// General factorization routines
template <typename T>
void LAPACKdgetrf(int *m, int *n, T *a, int *lda, int *ipiv,
                  int *info) = delete;
template <typename T>
void LAPACKdgetrs(const char *c, int *n, int *nrhs, T *a, int *lda, int *ipiv,
                  T *b, int *ldb, int *info) = delete;

// Factorization of packed-storage matrices
template <typename T>
void LAPACKdpptrf(const char *c, int *n, T *ap, int *info) = delete;
template <typename T>
void LAPACKdpptrs(const char *c, int *n, int *nrhs, T *ap, T *rhs, int *ldrhs,
                  int *info) = delete;

// Level 1 BLAS routines

template <>
inline double BLASddot(int *n, double *x, int *incx, double *y,
                       int *incy) = delete;
template <>
inline double BLASdnrm2(int *n, double *x, int *incx) = delete;
template <>
inline void BLASdaxpy(int *n, double *a, double *x, int *incx, double *y,
                      int *incy) = delete;
template <>
inline void BLASdscal(int *n, double *a, double *x, int *incx) = delete;

// Compute C := alpha*A*A**double + beta*C or C := alpha*A**double*A + beta*C
template <>
inline void BLASsyrk(const char *uplo, const char *trans, int *n, int *k,
                     double *alpha, double *a, int *lda, double *beta,
                     double *c, int *ldc) = delete;

// Solve A*x = b or A^double*x = b where A is in packed format
template <>
inline void BLAStpsv(const char *uplo, const char *transa, const char *diag,
                     int *n, double *a, double *x, int *incx) = delete;

// Level 2 BLAS routines
// y = alpha * A * x + beta * y, for a general matrix
template <>
inline void BLASgemv(const char *c, int *m, int *n, double *alpha, double *a,
                     int *lda, double *x, int *incx, double *beta, double *y,
                     int *incy) = delete;

// Level 3 BLAS routines
// C := alpha*op( A )*op( B ) + beta*C,
template <>
inline void BLASgemm(const char *ta, const char *tb, int *m, int *n, int *k,
                     double *alpha, double *a, int *lda, double *b, int *ldb,
                     double *beta, double *c, int *ldc) = delete;

// General factorization routines
template <>
inline void LAPACKdgetrf(int *m, int *n, double *a, int *lda, int *ipiv,
                         int *info) = delete;
template <>
inline void LAPACKdgetrs(const char *c, int *n, int *nrhs, double *a, int *lda,
                         int *ipiv, double *b, int *ldb, int *info) = delete;

// Factorization of packed-storage matrices
template <>
inline void LAPACKdpptrf(const char *c, int *n, double *ap, int *info) = delete;
template <>
inline void LAPACKdpptrs(const char *c, int *n, int *nrhs, double *ap,
                         double *rhs, int *ldrhs, int *info) = delete;
/* Specializations - double precision real */

/* Specializations - double precision complex */

// Level 1 BLAS routines
template <>
inline std::complex<double> BLASddot(int *n, std::complex<double> *x, int *incx,
                                     std::complex<double> *y,
                                     int *incy) = delete;
template <>
inline double BLASdnrm2(int *n, std::complex<double> *x, int *incx) = delete;
template <>
inline void BLASdaxpy(int *n, std::complex<double> *a, std::complex<double> *x,
                      int *incx, std::complex<double> *y, int *incy) = delete;
template <>
inline void BLASdscal(int *n, std::complex<double> *a, std::complex<double> *x,
                      int *incx) = delete;

// Compute C := alpha*A*A**std::complex<double> + beta*C or C :=
// alpha*A**std::complex<double>*A + beta*C
template <>
inline void BLASsyrk(const char *uplo, const char *trans, int *n, int *k,
                     std::complex<double> *alpha, std::complex<double> *a,
                     int *lda, std::complex<double> *beta,
                     std::complex<double> *c, int *ldc) = delete;

// Solve A*x = b or A^std::complex<double>*x = b where A is in packed format
template <>
inline void BLAStpsv(const char *uplo, const char *transa, const char *diag,
                     int *n, std::complex<double> *a, std::complex<double> *x,
                     int *incx) = delete;

// Level 2 BLAS routines
// y = alpha * A * x + beta * y, for a general matrix
template <>
inline void BLASgemv(const char *c, int *m, int *n, std::complex<double> *alpha,
                     std::complex<double> *a, int *lda, std::complex<double> *x,
                     int *incx, std::complex<double> *beta,
                     std::complex<double> *y, int *incy) = delete;

// Level 3 BLAS routines
// C := alpha*op( A )*op( B ) + beta*C,
template <>
inline void BLASgemm(const char *ta, const char *tb, int *m, int *n, int *k,
                     std::complex<double> *alpha, std::complex<double> *a,
                     int *lda, std::complex<double> *b, int *ldb,
                     std::complex<double> *beta, std::complex<double> *c,
                     int *ldc) = delete;

// General factorization routines
template <>
inline void LAPACKdgetrf(int *m, int *n, std::complex<double> *a, int *lda,
                         int *ipiv, int *info) = delete;
template <>
inline void LAPACKdgetrs(const char *c, int *n, int *nrhs,
                         std::complex<double> *a, int *lda, int *ipiv,
                         std::complex<double> *b, int *ldb, int *info) = delete;

// Factorization of packed-storage matrices
template <>
inline void LAPACKdpptrf(const char *c, int *n, std::complex<double> *ap,
                         int *info) = delete;
template <>
inline void LAPACKdpptrs(const char *c, int *n, int *nrhs,
                         std::complex<double> *ap, std::complex<double> *rhs,
                         int *ldrhs, int *info) = delete;

}  // namespace A2D

#endif