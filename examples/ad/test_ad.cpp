#include "ad/a2dgemm.h"
#include "ad/core/a2dgemmcore.h"

int main() {
  using T = double;

  // // A: 3x3, B: 3x3, C: 3,3
  // T A[1000], B[1000], C[1000];

  // A2D::MatMatMultCore<T, 3, 4, 4, 5, 5, 3, A2D::MatOp::NORMAL,
  //                     A2D::MatOp::NORMAL, A2D::MatOp::TRANSPOSE>(
  //     A, B,
  //     C);  // calls the general  function

  A2D::Mat<T, 3, 3> A, B, C, Ab, Bb, Cb;
  A2D::ADMat<A2D::Mat<T, 3, 3>> Amat(A, Ab), Bmat(B, Bb), Cmat(C, Cb);

  A2D::MatMatMult(Amat, Bmat, Cmat);
  A2D::MatMatMult(A, Bmat, Cmat);
  A2D::MatMatMult(A, B, C);
  A2D::MatMatMult(Amat, B, Cmat);
}