#include "ad/a2dgemm.h"
#include "ad/a2dstack.h"
#include "ad/core/a2dgemmcore.h"

int main() {
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
  exprH.forward();
  exprH.reverse();
  exprH.hforward();
  exprH.hreverse();
}