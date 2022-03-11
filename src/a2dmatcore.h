#ifndef A2D_MAT_CORE_H
#define A2D_MAT_CORE_H

/*
  Compute the determinant of a 3x3 matrix

  input:
  A:        a 3x3 matrix in row-major order

  returns:  the determinant of A
*/
inline TacsScalar Mat3x3DetCore( const TacsScalar A[] ){
  return (A[8] * (A[0] * A[4] - A[3] * A[1]) -
          A[7] * (A[0] * A[5] - A[3] * A[2]) +
          A[6] * (A[1] * A[5] - A[2] * A[4]));
}

inline TacsScalar Mat3x3DetDerivForwardCore( const TacsScalar A[],
                                             TacsScalar Ad[] ){
  return (
    Ad[0] * (A[8] * A[4] - A[7] * A[5]) +
    Ad[1] * (A[6] * A[5] - A[8] * A[3]) +
    Ad[2] * (A[7] * A[3] - A[6] * A[4]) +
    Ad[3] * (A[7] * A[2] - A[8] * A[1]) +
    Ad[4] * (A[8] * A[0] - A[6] * A[2]) +
    Ad[5] * (A[6] * A[1] - A[7] * A[0]) +
    Ad[6] * (A[1] * A[5] - A[2] * A[4]) +
    Ad[7] * (A[3] * A[2] - A[0] * A[5]) +
    Ad[8] * (A[0] * A[4] - A[3] * A[1]));
}

inline void Mat3x3DetDerivReverseCore( const TacsScalar detd,
                                       const TacsScalar A[],
                                       TacsScalar Ad[] ){
  Ad[0] += (A[8] * A[4] - A[7] * A[5]) * detd;
  Ad[1] += (A[6] * A[5] - A[8] * A[3]) * detd;
  Ad[2] += (A[7] * A[3] - A[6] * A[4]) * detd;
  Ad[3] += (A[7] * A[2] - A[8] * A[1]) * detd;
  Ad[4] += (A[8] * A[0] - A[6] * A[2]) * detd;
  Ad[5] += (A[6] * A[1] - A[7] * A[0]) * detd;
  Ad[6] += (A[1] * A[5] - A[2] * A[4]) * detd;
  Ad[7] += (A[3] * A[2] - A[0] * A[5]) * detd;
  Ad[8] += (A[0] * A[4] - A[3] * A[1]) * detd;
}

inline TacsScalar Symm3x3DetCore( const TacsScalar S[] ){
  return (S[5] * (S[0] * S[3] - S[1] * S[1]) -
          S[4] * (S[0] * S[4] - S[1] * S[2]) +
          S[2] * (S[1] * S[4] - S[2] * S[3]));
}

inline TacsScalar Symm3x3DetDerivForwardCore( const TacsScalar S[],
                                              const TacsScalar Sd[] ){
  return (
    Sd[0] * (S[5] * S[3] - S[4] * S[4]) +
    2.0 * Sd[1] * (S[2] * S[4] - S[5] * S[1]) +
    2.0 * Sd[2] * (S[1] * S[4] - S[3] * S[2]) +
    Sd[3] * (S[5] * S[0] - S[2] * S[2]) +
    2.0 * Sd[4] * (S[1] * S[2] - S[0] * S[4]) +
    Sd[5] * (S[0] * S[3] - S[1] * S[1]));
}

inline void Symm3x3DetDerivReverseCore( const TacsScalar detd,
                                        const TacsScalar S[],
                                        TacsScalar Sd[] ){
  Sd[0] += (S[5] * S[3] - S[4] * S[4]) * detd;
  Sd[1] += 2.0 * (S[2] * S[4] - S[5] * S[1]) * detd;
  Sd[2] += 2.0 * (S[1] * S[4] - S[3] * S[2]) * detd;
  Sd[3] += (S[5] * S[0] - S[2] * S[2]) * detd;
  Sd[4] += 2.0 * (S[1] * S[2] - S[0] * S[4]) * detd;
  Sd[5] += (S[0] * S[3] - S[1] * S[1]) * detd;
}

inline void Symm3x3SymmMultCore( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
  C[0] = (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
  C[1] = (A[0] * B[1] + A[1] * B[3] + A[2] * B[4]);
  C[2] = (A[0] * B[2] + A[1] * B[4] + A[2] * B[5]);
  C[3] = (A[1] * B[0] + A[3] * B[1] + A[4] * B[2]);
  C[4] = (A[1] * B[1] + A[3] * B[3] + A[4] * B[4]);
  C[5] = (A[1] * B[2] + A[3] * B[4] + A[4] * B[5]);
  C[6] = (A[2] * B[0] + A[4] * B[1] + A[5] * B[2]);
  C[7] = (A[2] * B[1] + A[4] * B[3] + A[5] * B[4]);
  C[8] = (A[2] * B[2] + A[4] * B[4] + A[5] * B[5]);
}

inline void Symm3x3SymmMultScaleCore( TacsScalar scale, const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
  C[0] = scale * (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
  C[1] = scale * (A[0] * B[1] + A[1] * B[3] + A[2] * B[4]);
  C[2] = scale * (A[0] * B[2] + A[1] * B[4] + A[2] * B[5]);
  C[3] = scale * (A[1] * B[0] + A[3] * B[1] + A[4] * B[2]);
  C[4] = scale * (A[1] * B[1] + A[3] * B[3] + A[4] * B[4]);
  C[5] = scale * (A[1] * B[2] + A[3] * B[4] + A[4] * B[5]);
  C[6] = scale * (A[2] * B[0] + A[4] * B[1] + A[5] * B[2]);
  C[7] = scale * (A[2] * B[1] + A[4] * B[3] + A[5] * B[4]);
  C[8] = scale * (A[2] * B[2] + A[4] * B[4] + A[5] * B[5]);
}

inline void Symm3x3SymmMultAddCore( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
  C[0] += (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
  C[1] += (A[0] * B[1] + A[1] * B[3] + A[2] * B[4]);
  C[2] += (A[0] * B[2] + A[1] * B[4] + A[2] * B[5]);
  C[3] += (A[1] * B[0] + A[3] * B[1] + A[4] * B[2]);
  C[4] += (A[1] * B[1] + A[3] * B[3] + A[4] * B[4]);
  C[5] += (A[1] * B[2] + A[3] * B[4] + A[4] * B[5]);
  C[6] += (A[2] * B[0] + A[4] * B[1] + A[5] * B[2]);
  C[7] += (A[2] * B[1] + A[4] * B[3] + A[5] * B[4]);
  C[8] += (A[2] * B[2] + A[4] * B[4] + A[5] * B[5]);
}

inline void Symm3x3SymmMultSubCore( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
  C[0] -= (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
  C[1] -= (A[0] * B[1] + A[1] * B[3] + A[2] * B[4]);
  C[2] -= (A[0] * B[2] + A[1] * B[4] + A[2] * B[5]);
  C[3] -= (A[1] * B[0] + A[3] * B[1] + A[4] * B[2]);
  C[4] -= (A[1] * B[1] + A[3] * B[3] + A[4] * B[4]);
  C[5] -= (A[1] * B[2] + A[3] * B[4] + A[4] * B[5]);
  C[6] -= (A[2] * B[0] + A[4] * B[1] + A[5] * B[2]);
  C[7] -= (A[2] * B[1] + A[4] * B[3] + A[5] * B[4]);
  C[8] -= (A[2] * B[2] + A[4] * B[4] + A[5] * B[5]);
}

inline void Symm3x3SymmMultAddScaleCore( TacsScalar scale, const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
  C[0] += scale * (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
  C[1] += scale * (A[0] * B[1] + A[1] * B[3] + A[2] * B[4]);
  C[2] += scale * (A[0] * B[2] + A[1] * B[4] + A[2] * B[5]);
  C[3] += scale * (A[1] * B[0] + A[3] * B[1] + A[4] * B[2]);
  C[4] += scale * (A[1] * B[1] + A[3] * B[3] + A[4] * B[4]);
  C[5] += scale * (A[1] * B[2] + A[3] * B[4] + A[4] * B[5]);
  C[6] += scale * (A[2] * B[0] + A[4] * B[1] + A[5] * B[2]);
  C[7] += scale * (A[2] * B[1] + A[4] * B[3] + A[5] * B[4]);
  C[8] += scale * (A[2] * B[2] + A[4] * B[4] + A[5] * B[5]);
}

inline void Symm3x3MatMultCore( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
  C[0] = (A[0] * B[0] + A[1] * B[3] + A[2] * B[6]);
  C[1] = (A[0] * B[1] + A[1] * B[4] + A[2] * B[7]);
  C[2] = (A[0] * B[2] + A[1] * B[5] + A[2] * B[8]);
  C[3] = (A[1] * B[0] + A[3] * B[3] + A[4] * B[6]);
  C[4] = (A[1] * B[1] + A[3] * B[4] + A[4] * B[7]);
  C[5] = (A[1] * B[2] + A[3] * B[5] + A[4] * B[8]);
  C[6] = (A[2] * B[0] + A[4] * B[3] + A[5] * B[6]);
  C[7] = (A[2] * B[1] + A[4] * B[4] + A[5] * B[7]);
  C[8] = (A[2] * B[2] + A[4] * B[5] + A[5] * B[8]);
}

inline void Symm3x3MatMultScaleCore( TacsScalar scale, const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
  C[0] = scale * (A[0] * B[0] + A[1] * B[3] + A[2] * B[6]);
  C[1] = scale * (A[0] * B[1] + A[1] * B[4] + A[2] * B[7]);
  C[2] = scale * (A[0] * B[2] + A[1] * B[5] + A[2] * B[8]);
  C[3] = scale * (A[1] * B[0] + A[3] * B[3] + A[4] * B[6]);
  C[4] = scale * (A[1] * B[1] + A[3] * B[4] + A[4] * B[7]);
  C[5] = scale * (A[1] * B[2] + A[3] * B[5] + A[4] * B[8]);
  C[6] = scale * (A[2] * B[0] + A[4] * B[3] + A[5] * B[6]);
  C[7] = scale * (A[2] * B[1] + A[4] * B[4] + A[5] * B[7]);
  C[8] = scale * (A[2] * B[2] + A[4] * B[5] + A[5] * B[8]);
}

inline void Symm3x3MatMultAddCore( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
  C[0] += (A[0] * B[0] + A[1] * B[3] + A[2] * B[6]);
  C[1] += (A[0] * B[1] + A[1] * B[4] + A[2] * B[7]);
  C[2] += (A[0] * B[2] + A[1] * B[5] + A[2] * B[8]);
  C[3] += (A[1] * B[0] + A[3] * B[3] + A[4] * B[6]);
  C[4] += (A[1] * B[1] + A[3] * B[4] + A[4] * B[7]);
  C[5] += (A[1] * B[2] + A[3] * B[5] + A[4] * B[8]);
  C[6] += (A[2] * B[0] + A[4] * B[3] + A[5] * B[6]);
  C[7] += (A[2] * B[1] + A[4] * B[4] + A[5] * B[7]);
  C[8] += (A[2] * B[2] + A[4] * B[5] + A[5] * B[8]);
}

inline void Symm3x3MatMultSubCore( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
  C[0] -= (A[0] * B[0] + A[1] * B[3] + A[2] * B[6]);
  C[1] -= (A[0] * B[1] + A[1] * B[4] + A[2] * B[7]);
  C[2] -= (A[0] * B[2] + A[1] * B[5] + A[2] * B[8]);
  C[3] -= (A[1] * B[0] + A[3] * B[3] + A[4] * B[6]);
  C[4] -= (A[1] * B[1] + A[3] * B[4] + A[4] * B[7]);
  C[5] -= (A[1] * B[2] + A[3] * B[5] + A[4] * B[8]);
  C[6] -= (A[2] * B[0] + A[4] * B[3] + A[5] * B[6]);
  C[7] -= (A[2] * B[1] + A[4] * B[4] + A[5] * B[7]);
  C[8] -= (A[2] * B[2] + A[4] * B[5] + A[5] * B[8]);
}

inline void Symm3x3MatMultAddScaleCore( TacsScalar scale, const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
  C[0] += scale * (A[0] * B[0] + A[1] * B[3] + A[2] * B[6]);
  C[1] += scale * (A[0] * B[1] + A[1] * B[4] + A[2] * B[7]);
  C[2] += scale * (A[0] * B[2] + A[1] * B[5] + A[2] * B[8]);
  C[3] += scale * (A[1] * B[0] + A[3] * B[3] + A[4] * B[6]);
  C[4] += scale * (A[1] * B[1] + A[3] * B[4] + A[4] * B[7]);
  C[5] += scale * (A[1] * B[2] + A[3] * B[5] + A[4] * B[8]);
  C[6] += scale * (A[2] * B[0] + A[4] * B[3] + A[5] * B[6]);
  C[7] += scale * (A[2] * B[1] + A[4] * B[4] + A[5] * B[7]);
  C[8] += scale * (A[2] * B[2] + A[4] * B[5] + A[5] * B[8]);
}

inline void Symm3x3MatTransMultCore( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
  C[0] = (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
  C[1] = (A[0] * B[3] + A[1] * B[4] + A[2] * B[5]);
  C[2] = (A[0] * B[6] + A[1] * B[7] + A[2] * B[8]);
  C[3] = (A[1] * B[0] + A[3] * B[1] + A[4] * B[2]);
  C[4] = (A[1] * B[3] + A[3] * B[4] + A[4] * B[5]);
  C[5] = (A[1] * B[6] + A[3] * B[7] + A[4] * B[8]);
  C[6] = (A[2] * B[0] + A[4] * B[1] + A[5] * B[2]);
  C[7] = (A[2] * B[3] + A[4] * B[4] + A[5] * B[5]);
  C[8] = (A[2] * B[6] + A[4] * B[7] + A[5] * B[8]);
}

inline void Symm3x3MatTransMultScaleCore( TacsScalar scale, const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
  C[0] = scale * (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
  C[1] = scale * (A[0] * B[3] + A[1] * B[4] + A[2] * B[5]);
  C[2] = scale * (A[0] * B[6] + A[1] * B[7] + A[2] * B[8]);
  C[3] = scale * (A[1] * B[0] + A[3] * B[1] + A[4] * B[2]);
  C[4] = scale * (A[1] * B[3] + A[3] * B[4] + A[4] * B[5]);
  C[5] = scale * (A[1] * B[6] + A[3] * B[7] + A[4] * B[8]);
  C[6] = scale * (A[2] * B[0] + A[4] * B[1] + A[5] * B[2]);
  C[7] = scale * (A[2] * B[3] + A[4] * B[4] + A[5] * B[5]);
  C[8] = scale * (A[2] * B[6] + A[4] * B[7] + A[5] * B[8]);
}

inline void Symm3x3MatTransMultAddCore( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
  C[0] += (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
  C[1] += (A[0] * B[3] + A[1] * B[4] + A[2] * B[5]);
  C[2] += (A[0] * B[6] + A[1] * B[7] + A[2] * B[8]);
  C[3] += (A[1] * B[0] + A[3] * B[1] + A[4] * B[2]);
  C[4] += (A[1] * B[3] + A[3] * B[4] + A[4] * B[5]);
  C[5] += (A[1] * B[6] + A[3] * B[7] + A[4] * B[8]);
  C[6] += (A[2] * B[0] + A[4] * B[1] + A[5] * B[2]);
  C[7] += (A[2] * B[3] + A[4] * B[4] + A[5] * B[5]);
  C[8] += (A[2] * B[6] + A[4] * B[7] + A[5] * B[8]);
}

inline void Symm3x3MatTransMultSubCore( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
  C[0] -= (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
  C[1] -= (A[0] * B[3] + A[1] * B[4] + A[2] * B[5]);
  C[2] -= (A[0] * B[6] + A[1] * B[7] + A[2] * B[8]);
  C[3] -= (A[1] * B[0] + A[3] * B[1] + A[4] * B[2]);
  C[4] -= (A[1] * B[3] + A[3] * B[4] + A[4] * B[5]);
  C[5] -= (A[1] * B[6] + A[3] * B[7] + A[4] * B[8]);
  C[6] -= (A[2] * B[0] + A[4] * B[1] + A[5] * B[2]);
  C[7] -= (A[2] * B[3] + A[4] * B[4] + A[5] * B[5]);
  C[8] -= (A[2] * B[6] + A[4] * B[7] + A[5] * B[8]);
}

inline void Symm3x3MatTransMultAddScaleCore( TacsScalar scale, const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
  C[0] += scale * (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
  C[1] += scale * (A[0] * B[3] + A[1] * B[4] + A[2] * B[5]);
  C[2] += scale * (A[0] * B[6] + A[1] * B[7] + A[2] * B[8]);
  C[3] += scale * (A[1] * B[0] + A[3] * B[1] + A[4] * B[2]);
  C[4] += scale * (A[1] * B[3] + A[3] * B[4] + A[4] * B[5]);
  C[5] += scale * (A[1] * B[6] + A[3] * B[7] + A[4] * B[8]);
  C[6] += scale * (A[2] * B[0] + A[4] * B[1] + A[5] * B[2]);
  C[7] += scale * (A[2] * B[3] + A[4] * B[4] + A[5] * B[5]);
  C[8] += scale * (A[2] * B[6] + A[4] * B[7] + A[5] * B[8]);
}

inline void Mat3x3SymmMultCore( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
  C[0] = (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
  C[1] = (A[0] * B[1] + A[1] * B[3] + A[2] * B[4]);
  C[2] = (A[0] * B[2] + A[1] * B[4] + A[2] * B[5]);
  C[3] = (A[3] * B[0] + A[4] * B[1] + A[5] * B[2]);
  C[4] = (A[3] * B[1] + A[4] * B[3] + A[5] * B[4]);
  C[5] = (A[3] * B[2] + A[4] * B[4] + A[5] * B[5]);
  C[6] = (A[6] * B[0] + A[7] * B[1] + A[8] * B[2]);
  C[7] = (A[6] * B[1] + A[7] * B[3] + A[8] * B[4]);
  C[8] = (A[6] * B[2] + A[7] * B[4] + A[8] * B[5]);
}

inline void Mat3x3SymmMultScaleCore( TacsScalar scale, const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
  C[0] = scale * (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
  C[1] = scale * (A[0] * B[1] + A[1] * B[3] + A[2] * B[4]);
  C[2] = scale * (A[0] * B[2] + A[1] * B[4] + A[2] * B[5]);
  C[3] = scale * (A[3] * B[0] + A[4] * B[1] + A[5] * B[2]);
  C[4] = scale * (A[3] * B[1] + A[4] * B[3] + A[5] * B[4]);
  C[5] = scale * (A[3] * B[2] + A[4] * B[4] + A[5] * B[5]);
  C[6] = scale * (A[6] * B[0] + A[7] * B[1] + A[8] * B[2]);
  C[7] = scale * (A[6] * B[1] + A[7] * B[3] + A[8] * B[4]);
  C[8] = scale * (A[6] * B[2] + A[7] * B[4] + A[8] * B[5]);
}

inline void Mat3x3SymmMultAddCore( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
  C[0] += (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
  C[1] += (A[0] * B[1] + A[1] * B[3] + A[2] * B[4]);
  C[2] += (A[0] * B[2] + A[1] * B[4] + A[2] * B[5]);
  C[3] += (A[3] * B[0] + A[4] * B[1] + A[5] * B[2]);
  C[4] += (A[3] * B[1] + A[4] * B[3] + A[5] * B[4]);
  C[5] += (A[3] * B[2] + A[4] * B[4] + A[5] * B[5]);
  C[6] += (A[6] * B[0] + A[7] * B[1] + A[8] * B[2]);
  C[7] += (A[6] * B[1] + A[7] * B[3] + A[8] * B[4]);
  C[8] += (A[6] * B[2] + A[7] * B[4] + A[8] * B[5]);
}

inline void Mat3x3SymmMultSubCore( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
  C[0] -= (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
  C[1] -= (A[0] * B[1] + A[1] * B[3] + A[2] * B[4]);
  C[2] -= (A[0] * B[2] + A[1] * B[4] + A[2] * B[5]);
  C[3] -= (A[3] * B[0] + A[4] * B[1] + A[5] * B[2]);
  C[4] -= (A[3] * B[1] + A[4] * B[3] + A[5] * B[4]);
  C[5] -= (A[3] * B[2] + A[4] * B[4] + A[5] * B[5]);
  C[6] -= (A[6] * B[0] + A[7] * B[1] + A[8] * B[2]);
  C[7] -= (A[6] * B[1] + A[7] * B[3] + A[8] * B[4]);
  C[8] -= (A[6] * B[2] + A[7] * B[4] + A[8] * B[5]);
}

inline void Mat3x3SymmMultAddScaleCore( TacsScalar scale, const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
  C[0] += scale * (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
  C[1] += scale * (A[0] * B[1] + A[1] * B[3] + A[2] * B[4]);
  C[2] += scale * (A[0] * B[2] + A[1] * B[4] + A[2] * B[5]);
  C[3] += scale * (A[3] * B[0] + A[4] * B[1] + A[5] * B[2]);
  C[4] += scale * (A[3] * B[1] + A[4] * B[3] + A[5] * B[4]);
  C[5] += scale * (A[3] * B[2] + A[4] * B[4] + A[5] * B[5]);
  C[6] += scale * (A[6] * B[0] + A[7] * B[1] + A[8] * B[2]);
  C[7] += scale * (A[6] * B[1] + A[7] * B[3] + A[8] * B[4]);
  C[8] += scale * (A[6] * B[2] + A[7] * B[4] + A[8] * B[5]);
}

inline void MatTrans3x3SymmMultCore( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
  C[0] = (A[0] * B[0] + A[3] * B[1] + A[6] * B[2]);
  C[1] = (A[0] * B[1] + A[3] * B[3] + A[6] * B[4]);
  C[2] = (A[0] * B[2] + A[3] * B[4] + A[6] * B[5]);
  C[3] = (A[1] * B[0] + A[4] * B[1] + A[7] * B[2]);
  C[4] = (A[1] * B[1] + A[4] * B[3] + A[7] * B[4]);
  C[5] = (A[1] * B[2] + A[4] * B[4] + A[7] * B[5]);
  C[6] = (A[2] * B[0] + A[5] * B[1] + A[8] * B[2]);
  C[7] = (A[2] * B[1] + A[5] * B[3] + A[8] * B[4]);
  C[8] = (A[2] * B[2] + A[5] * B[4] + A[8] * B[5]);
}

inline void MatTrans3x3SymmMultScaleCore( TacsScalar scale, const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
  C[0] = scale * (A[0] * B[0] + A[3] * B[1] + A[6] * B[2]);
  C[1] = scale * (A[0] * B[1] + A[3] * B[3] + A[6] * B[4]);
  C[2] = scale * (A[0] * B[2] + A[3] * B[4] + A[6] * B[5]);
  C[3] = scale * (A[1] * B[0] + A[4] * B[1] + A[7] * B[2]);
  C[4] = scale * (A[1] * B[1] + A[4] * B[3] + A[7] * B[4]);
  C[5] = scale * (A[1] * B[2] + A[4] * B[4] + A[7] * B[5]);
  C[6] = scale * (A[2] * B[0] + A[5] * B[1] + A[8] * B[2]);
  C[7] = scale * (A[2] * B[1] + A[5] * B[3] + A[8] * B[4]);
  C[8] = scale * (A[2] * B[2] + A[5] * B[4] + A[8] * B[5]);
}

inline void MatTrans3x3SymmMultAddCore( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
  C[0] += (A[0] * B[0] + A[3] * B[1] + A[6] * B[2]);
  C[1] += (A[0] * B[1] + A[3] * B[3] + A[6] * B[4]);
  C[2] += (A[0] * B[2] + A[3] * B[4] + A[6] * B[5]);
  C[3] += (A[1] * B[0] + A[4] * B[1] + A[7] * B[2]);
  C[4] += (A[1] * B[1] + A[4] * B[3] + A[7] * B[4]);
  C[5] += (A[1] * B[2] + A[4] * B[4] + A[7] * B[5]);
  C[6] += (A[2] * B[0] + A[5] * B[1] + A[8] * B[2]);
  C[7] += (A[2] * B[1] + A[5] * B[3] + A[8] * B[4]);
  C[8] += (A[2] * B[2] + A[5] * B[4] + A[8] * B[5]);
}

inline void MatTrans3x3SymmMultSubCore( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
  C[0] -= (A[0] * B[0] + A[3] * B[1] + A[6] * B[2]);
  C[1] -= (A[0] * B[1] + A[3] * B[3] + A[6] * B[4]);
  C[2] -= (A[0] * B[2] + A[3] * B[4] + A[6] * B[5]);
  C[3] -= (A[1] * B[0] + A[4] * B[1] + A[7] * B[2]);
  C[4] -= (A[1] * B[1] + A[4] * B[3] + A[7] * B[4]);
  C[5] -= (A[1] * B[2] + A[4] * B[4] + A[7] * B[5]);
  C[6] -= (A[2] * B[0] + A[5] * B[1] + A[8] * B[2]);
  C[7] -= (A[2] * B[1] + A[5] * B[3] + A[8] * B[4]);
  C[8] -= (A[2] * B[2] + A[5] * B[4] + A[8] * B[5]);
}

inline void MatTrans3x3SymmMultAddScaleCore( TacsScalar scale, const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
  C[0] += scale * (A[0] * B[0] + A[3] * B[1] + A[6] * B[2]);
  C[1] += scale * (A[0] * B[1] + A[3] * B[3] + A[6] * B[4]);
  C[2] += scale * (A[0] * B[2] + A[3] * B[4] + A[6] * B[5]);
  C[3] += scale * (A[1] * B[0] + A[4] * B[1] + A[7] * B[2]);
  C[4] += scale * (A[1] * B[1] + A[4] * B[3] + A[7] * B[4]);
  C[5] += scale * (A[1] * B[2] + A[4] * B[4] + A[7] * B[5]);
  C[6] += scale * (A[2] * B[0] + A[5] * B[1] + A[8] * B[2]);
  C[7] += scale * (A[2] * B[1] + A[5] * B[3] + A[8] * B[4]);
  C[8] += scale * (A[2] * B[2] + A[5] * B[4] + A[8] * B[5]);
}

inline void Mat3x3MatMultCore( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
  C[0] = (A[0] * B[0] + A[1] * B[3] + A[2] * B[6]);
  C[1] = (A[0] * B[1] + A[1] * B[4] + A[2] * B[7]);
  C[2] = (A[0] * B[2] + A[1] * B[5] + A[2] * B[8]);
  C[3] = (A[3] * B[0] + A[4] * B[3] + A[5] * B[6]);
  C[4] = (A[3] * B[1] + A[4] * B[4] + A[5] * B[7]);
  C[5] = (A[3] * B[2] + A[4] * B[5] + A[5] * B[8]);
  C[6] = (A[6] * B[0] + A[7] * B[3] + A[8] * B[6]);
  C[7] = (A[6] * B[1] + A[7] * B[4] + A[8] * B[7]);
  C[8] = (A[6] * B[2] + A[7] * B[5] + A[8] * B[8]);
}

inline void Mat3x3MatMultScaleCore( TacsScalar scale, const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
  C[0] = scale * (A[0] * B[0] + A[1] * B[3] + A[2] * B[6]);
  C[1] = scale * (A[0] * B[1] + A[1] * B[4] + A[2] * B[7]);
  C[2] = scale * (A[0] * B[2] + A[1] * B[5] + A[2] * B[8]);
  C[3] = scale * (A[3] * B[0] + A[4] * B[3] + A[5] * B[6]);
  C[4] = scale * (A[3] * B[1] + A[4] * B[4] + A[5] * B[7]);
  C[5] = scale * (A[3] * B[2] + A[4] * B[5] + A[5] * B[8]);
  C[6] = scale * (A[6] * B[0] + A[7] * B[3] + A[8] * B[6]);
  C[7] = scale * (A[6] * B[1] + A[7] * B[4] + A[8] * B[7]);
  C[8] = scale * (A[6] * B[2] + A[7] * B[5] + A[8] * B[8]);
}

inline void Mat3x3MatMultAddCore( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
  C[0] += (A[0] * B[0] + A[1] * B[3] + A[2] * B[6]);
  C[1] += (A[0] * B[1] + A[1] * B[4] + A[2] * B[7]);
  C[2] += (A[0] * B[2] + A[1] * B[5] + A[2] * B[8]);
  C[3] += (A[3] * B[0] + A[4] * B[3] + A[5] * B[6]);
  C[4] += (A[3] * B[1] + A[4] * B[4] + A[5] * B[7]);
  C[5] += (A[3] * B[2] + A[4] * B[5] + A[5] * B[8]);
  C[6] += (A[6] * B[0] + A[7] * B[3] + A[8] * B[6]);
  C[7] += (A[6] * B[1] + A[7] * B[4] + A[8] * B[7]);
  C[8] += (A[6] * B[2] + A[7] * B[5] + A[8] * B[8]);
}

inline void Mat3x3MatMultSubCore( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
  C[0] -= (A[0] * B[0] + A[1] * B[3] + A[2] * B[6]);
  C[1] -= (A[0] * B[1] + A[1] * B[4] + A[2] * B[7]);
  C[2] -= (A[0] * B[2] + A[1] * B[5] + A[2] * B[8]);
  C[3] -= (A[3] * B[0] + A[4] * B[3] + A[5] * B[6]);
  C[4] -= (A[3] * B[1] + A[4] * B[4] + A[5] * B[7]);
  C[5] -= (A[3] * B[2] + A[4] * B[5] + A[5] * B[8]);
  C[6] -= (A[6] * B[0] + A[7] * B[3] + A[8] * B[6]);
  C[7] -= (A[6] * B[1] + A[7] * B[4] + A[8] * B[7]);
  C[8] -= (A[6] * B[2] + A[7] * B[5] + A[8] * B[8]);
}

inline void Mat3x3MatMultAddScaleCore( TacsScalar scale, const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
  C[0] += scale * (A[0] * B[0] + A[1] * B[3] + A[2] * B[6]);
  C[1] += scale * (A[0] * B[1] + A[1] * B[4] + A[2] * B[7]);
  C[2] += scale * (A[0] * B[2] + A[1] * B[5] + A[2] * B[8]);
  C[3] += scale * (A[3] * B[0] + A[4] * B[3] + A[5] * B[6]);
  C[4] += scale * (A[3] * B[1] + A[4] * B[4] + A[5] * B[7]);
  C[5] += scale * (A[3] * B[2] + A[4] * B[5] + A[5] * B[8]);
  C[6] += scale * (A[6] * B[0] + A[7] * B[3] + A[8] * B[6]);
  C[7] += scale * (A[6] * B[1] + A[7] * B[4] + A[8] * B[7]);
  C[8] += scale * (A[6] * B[2] + A[7] * B[5] + A[8] * B[8]);
}

inline void Mat3x3MatTransMultCore( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
  C[0] = (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
  C[1] = (A[0] * B[3] + A[1] * B[4] + A[2] * B[5]);
  C[2] = (A[0] * B[6] + A[1] * B[7] + A[2] * B[8]);
  C[3] = (A[3] * B[0] + A[4] * B[1] + A[5] * B[2]);
  C[4] = (A[3] * B[3] + A[4] * B[4] + A[5] * B[5]);
  C[5] = (A[3] * B[6] + A[4] * B[7] + A[5] * B[8]);
  C[6] = (A[6] * B[0] + A[7] * B[1] + A[8] * B[2]);
  C[7] = (A[6] * B[3] + A[7] * B[4] + A[8] * B[5]);
  C[8] = (A[6] * B[6] + A[7] * B[7] + A[8] * B[8]);
}

inline void Mat3x3MatTransMultScaleCore( TacsScalar scale, const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
  C[0] = scale * (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
  C[1] = scale * (A[0] * B[3] + A[1] * B[4] + A[2] * B[5]);
  C[2] = scale * (A[0] * B[6] + A[1] * B[7] + A[2] * B[8]);
  C[3] = scale * (A[3] * B[0] + A[4] * B[1] + A[5] * B[2]);
  C[4] = scale * (A[3] * B[3] + A[4] * B[4] + A[5] * B[5]);
  C[5] = scale * (A[3] * B[6] + A[4] * B[7] + A[5] * B[8]);
  C[6] = scale * (A[6] * B[0] + A[7] * B[1] + A[8] * B[2]);
  C[7] = scale * (A[6] * B[3] + A[7] * B[4] + A[8] * B[5]);
  C[8] = scale * (A[6] * B[6] + A[7] * B[7] + A[8] * B[8]);
}

inline void Mat3x3MatTransMultAddCore( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
  C[0] += (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
  C[1] += (A[0] * B[3] + A[1] * B[4] + A[2] * B[5]);
  C[2] += (A[0] * B[6] + A[1] * B[7] + A[2] * B[8]);
  C[3] += (A[3] * B[0] + A[4] * B[1] + A[5] * B[2]);
  C[4] += (A[3] * B[3] + A[4] * B[4] + A[5] * B[5]);
  C[5] += (A[3] * B[6] + A[4] * B[7] + A[5] * B[8]);
  C[6] += (A[6] * B[0] + A[7] * B[1] + A[8] * B[2]);
  C[7] += (A[6] * B[3] + A[7] * B[4] + A[8] * B[5]);
  C[8] += (A[6] * B[6] + A[7] * B[7] + A[8] * B[8]);
}

inline void Mat3x3MatTransMultSubCore( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
  C[0] -= (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
  C[1] -= (A[0] * B[3] + A[1] * B[4] + A[2] * B[5]);
  C[2] -= (A[0] * B[6] + A[1] * B[7] + A[2] * B[8]);
  C[3] -= (A[3] * B[0] + A[4] * B[1] + A[5] * B[2]);
  C[4] -= (A[3] * B[3] + A[4] * B[4] + A[5] * B[5]);
  C[5] -= (A[3] * B[6] + A[4] * B[7] + A[5] * B[8]);
  C[6] -= (A[6] * B[0] + A[7] * B[1] + A[8] * B[2]);
  C[7] -= (A[6] * B[3] + A[7] * B[4] + A[8] * B[5]);
  C[8] -= (A[6] * B[6] + A[7] * B[7] + A[8] * B[8]);
}

inline void Mat3x3MatTransMultAddScaleCore( TacsScalar scale, const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
  C[0] += scale * (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
  C[1] += scale * (A[0] * B[3] + A[1] * B[4] + A[2] * B[5]);
  C[2] += scale * (A[0] * B[6] + A[1] * B[7] + A[2] * B[8]);
  C[3] += scale * (A[3] * B[0] + A[4] * B[1] + A[5] * B[2]);
  C[4] += scale * (A[3] * B[3] + A[4] * B[4] + A[5] * B[5]);
  C[5] += scale * (A[3] * B[6] + A[4] * B[7] + A[5] * B[8]);
  C[6] += scale * (A[6] * B[0] + A[7] * B[1] + A[8] * B[2]);
  C[7] += scale * (A[6] * B[3] + A[7] * B[4] + A[8] * B[5]);
  C[8] += scale * (A[6] * B[6] + A[7] * B[7] + A[8] * B[8]);
}

inline void MatTrans3x3MatMultCore( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
  C[0] = (A[0] * B[0] + A[3] * B[3] + A[6] * B[6]);
  C[1] = (A[0] * B[1] + A[3] * B[4] + A[6] * B[7]);
  C[2] = (A[0] * B[2] + A[3] * B[5] + A[6] * B[8]);
  C[3] = (A[1] * B[0] + A[4] * B[3] + A[7] * B[6]);
  C[4] = (A[1] * B[1] + A[4] * B[4] + A[7] * B[7]);
  C[5] = (A[1] * B[2] + A[4] * B[5] + A[7] * B[8]);
  C[6] = (A[2] * B[0] + A[5] * B[3] + A[8] * B[6]);
  C[7] = (A[2] * B[1] + A[5] * B[4] + A[8] * B[7]);
  C[8] = (A[2] * B[2] + A[5] * B[5] + A[8] * B[8]);
}

inline void MatTrans3x3MatMultScaleCore( TacsScalar scale, const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
  C[0] = scale * (A[0] * B[0] + A[3] * B[3] + A[6] * B[6]);
  C[1] = scale * (A[0] * B[1] + A[3] * B[4] + A[6] * B[7]);
  C[2] = scale * (A[0] * B[2] + A[3] * B[5] + A[6] * B[8]);
  C[3] = scale * (A[1] * B[0] + A[4] * B[3] + A[7] * B[6]);
  C[4] = scale * (A[1] * B[1] + A[4] * B[4] + A[7] * B[7]);
  C[5] = scale * (A[1] * B[2] + A[4] * B[5] + A[7] * B[8]);
  C[6] = scale * (A[2] * B[0] + A[5] * B[3] + A[8] * B[6]);
  C[7] = scale * (A[2] * B[1] + A[5] * B[4] + A[8] * B[7]);
  C[8] = scale * (A[2] * B[2] + A[5] * B[5] + A[8] * B[8]);
}

inline void MatTrans3x3MatMultAddCore( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
  C[0] += (A[0] * B[0] + A[3] * B[3] + A[6] * B[6]);
  C[1] += (A[0] * B[1] + A[3] * B[4] + A[6] * B[7]);
  C[2] += (A[0] * B[2] + A[3] * B[5] + A[6] * B[8]);
  C[3] += (A[1] * B[0] + A[4] * B[3] + A[7] * B[6]);
  C[4] += (A[1] * B[1] + A[4] * B[4] + A[7] * B[7]);
  C[5] += (A[1] * B[2] + A[4] * B[5] + A[7] * B[8]);
  C[6] += (A[2] * B[0] + A[5] * B[3] + A[8] * B[6]);
  C[7] += (A[2] * B[1] + A[5] * B[4] + A[8] * B[7]);
  C[8] += (A[2] * B[2] + A[5] * B[5] + A[8] * B[8]);
}

inline void MatTrans3x3MatMultSubCore( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
  C[0] -= (A[0] * B[0] + A[3] * B[3] + A[6] * B[6]);
  C[1] -= (A[0] * B[1] + A[3] * B[4] + A[6] * B[7]);
  C[2] -= (A[0] * B[2] + A[3] * B[5] + A[6] * B[8]);
  C[3] -= (A[1] * B[0] + A[4] * B[3] + A[7] * B[6]);
  C[4] -= (A[1] * B[1] + A[4] * B[4] + A[7] * B[7]);
  C[5] -= (A[1] * B[2] + A[4] * B[5] + A[7] * B[8]);
  C[6] -= (A[2] * B[0] + A[5] * B[3] + A[8] * B[6]);
  C[7] -= (A[2] * B[1] + A[5] * B[4] + A[8] * B[7]);
  C[8] -= (A[2] * B[2] + A[5] * B[5] + A[8] * B[8]);
}

inline void MatTrans3x3MatMultAddScaleCore( TacsScalar scale, const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
  C[0] += scale * (A[0] * B[0] + A[3] * B[3] + A[6] * B[6]);
  C[1] += scale * (A[0] * B[1] + A[3] * B[4] + A[6] * B[7]);
  C[2] += scale * (A[0] * B[2] + A[3] * B[5] + A[6] * B[8]);
  C[3] += scale * (A[1] * B[0] + A[4] * B[3] + A[7] * B[6]);
  C[4] += scale * (A[1] * B[1] + A[4] * B[4] + A[7] * B[7]);
  C[5] += scale * (A[1] * B[2] + A[4] * B[5] + A[7] * B[8]);
  C[6] += scale * (A[2] * B[0] + A[5] * B[3] + A[8] * B[6]);
  C[7] += scale * (A[2] * B[1] + A[5] * B[4] + A[8] * B[7]);
  C[8] += scale * (A[2] * B[2] + A[5] * B[5] + A[8] * B[8]);
}

inline void MatTrans3x3MatTransMultCore( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
  C[0] = (A[0] * B[0] + A[3] * B[1] + A[6] * B[2]);
  C[1] = (A[0] * B[3] + A[3] * B[4] + A[6] * B[5]);
  C[2] = (A[0] * B[6] + A[3] * B[7] + A[6] * B[8]);
  C[3] = (A[1] * B[0] + A[4] * B[1] + A[7] * B[2]);
  C[4] = (A[1] * B[3] + A[4] * B[4] + A[7] * B[5]);
  C[5] = (A[1] * B[6] + A[4] * B[7] + A[7] * B[8]);
  C[6] = (A[2] * B[0] + A[5] * B[1] + A[8] * B[2]);
  C[7] = (A[2] * B[3] + A[5] * B[4] + A[8] * B[5]);
  C[8] = (A[2] * B[6] + A[5] * B[7] + A[8] * B[8]);
}

inline void MatTrans3x3MatTransMultScaleCore( TacsScalar scale, const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
  C[0] = scale * (A[0] * B[0] + A[3] * B[1] + A[6] * B[2]);
  C[1] = scale * (A[0] * B[3] + A[3] * B[4] + A[6] * B[5]);
  C[2] = scale * (A[0] * B[6] + A[3] * B[7] + A[6] * B[8]);
  C[3] = scale * (A[1] * B[0] + A[4] * B[1] + A[7] * B[2]);
  C[4] = scale * (A[1] * B[3] + A[4] * B[4] + A[7] * B[5]);
  C[5] = scale * (A[1] * B[6] + A[4] * B[7] + A[7] * B[8]);
  C[6] = scale * (A[2] * B[0] + A[5] * B[1] + A[8] * B[2]);
  C[7] = scale * (A[2] * B[3] + A[5] * B[4] + A[8] * B[5]);
  C[8] = scale * (A[2] * B[6] + A[5] * B[7] + A[8] * B[8]);
}

inline void MatTrans3x3MatTransMultAddCore( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
  C[0] += (A[0] * B[0] + A[3] * B[1] + A[6] * B[2]);
  C[1] += (A[0] * B[3] + A[3] * B[4] + A[6] * B[5]);
  C[2] += (A[0] * B[6] + A[3] * B[7] + A[6] * B[8]);
  C[3] += (A[1] * B[0] + A[4] * B[1] + A[7] * B[2]);
  C[4] += (A[1] * B[3] + A[4] * B[4] + A[7] * B[5]);
  C[5] += (A[1] * B[6] + A[4] * B[7] + A[7] * B[8]);
  C[6] += (A[2] * B[0] + A[5] * B[1] + A[8] * B[2]);
  C[7] += (A[2] * B[3] + A[5] * B[4] + A[8] * B[5]);
  C[8] += (A[2] * B[6] + A[5] * B[7] + A[8] * B[8]);
}

inline void MatTrans3x3MatTransMultSubCore( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
  C[0] -= (A[0] * B[0] + A[3] * B[1] + A[6] * B[2]);
  C[1] -= (A[0] * B[3] + A[3] * B[4] + A[6] * B[5]);
  C[2] -= (A[0] * B[6] + A[3] * B[7] + A[6] * B[8]);
  C[3] -= (A[1] * B[0] + A[4] * B[1] + A[7] * B[2]);
  C[4] -= (A[1] * B[3] + A[4] * B[4] + A[7] * B[5]);
  C[5] -= (A[1] * B[6] + A[4] * B[7] + A[7] * B[8]);
  C[6] -= (A[2] * B[0] + A[5] * B[1] + A[8] * B[2]);
  C[7] -= (A[2] * B[3] + A[5] * B[4] + A[8] * B[5]);
  C[8] -= (A[2] * B[6] + A[5] * B[7] + A[8] * B[8]);
}

inline void MatTrans3x3MatTransMultAddScaleCore( TacsScalar scale, const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
  C[0] += scale * (A[0] * B[0] + A[3] * B[1] + A[6] * B[2]);
  C[1] += scale * (A[0] * B[3] + A[3] * B[4] + A[6] * B[5]);
  C[2] += scale * (A[0] * B[6] + A[3] * B[7] + A[6] * B[8]);
  C[3] += scale * (A[1] * B[0] + A[4] * B[1] + A[7] * B[2]);
  C[4] += scale * (A[1] * B[3] + A[4] * B[4] + A[7] * B[5]);
  C[5] += scale * (A[1] * B[6] + A[4] * B[7] + A[7] * B[8]);
  C[6] += scale * (A[2] * B[0] + A[5] * B[1] + A[8] * B[2]);
  C[7] += scale * (A[2] * B[3] + A[5] * B[4] + A[8] * B[5]);
  C[8] += scale * (A[2] * B[6] + A[5] * B[7] + A[8] * B[8]);
}

inline TacsScalar Mat3x3InverseCore( const TacsScalar A[], TacsScalar Ainv[] ){
  TacsScalar det = (A[8] * (A[0] * A[4] - A[3] * A[1]) -
                    A[7] * (A[0] * A[5] - A[3] * A[2]) +
                    A[6] * (A[1] * A[5] - A[2] * A[4]));
  TacsScalar detinv = 1.0/det;

  Ainv[0] = (A[4] * A[8] - A[5] * A[7]) * detinv;
  Ainv[1] =-(A[1] * A[8] - A[2] * A[7]) * detinv;
  Ainv[2] = (A[1] * A[5] - A[2] * A[4]) * detinv;

  Ainv[3] =-(A[3] * A[8] - A[5] * A[6]) * detinv;
  Ainv[4] = (A[0] * A[8] - A[2] * A[6]) * detinv;
  Ainv[5] =-(A[0] * A[5] - A[2] * A[3]) * detinv;

  Ainv[6] = (A[3] * A[7] - A[4] * A[6]) * detinv;
  Ainv[7] =-(A[0] * A[7] - A[1] * A[6]) * detinv;
  Ainv[8] = (A[0] * A[4] - A[1] * A[3]) * detinv;

  return det;
}

inline void Mat3x3InverseDerivForwardCore( const TacsScalar Ainv[],
                                           const TacsScalar Ad[],
                                           TacsScalar Bd[] ){
  TacsScalar t[9];
  Mat3x3MatMultCore(Ainv, Ad, t);
  Mat3x3MatMultScaleCore(-1.0, t, Ainv, Bd);
}

inline void Mat3x3InverseDerivReverseCore( const TacsScalar Ainv[],
                                           const TacsScalar Bd[],
                                           TacsScalar Ad[] ){
  TacsScalar t[9];
  MatTrans3x3MatMultCore(Ainv, Bd, t);
  Mat3x3MatTransMultAddScaleCore(-1.0, t, Ainv, Ad);
}

#endif // A2D_MAT_CORE_H
