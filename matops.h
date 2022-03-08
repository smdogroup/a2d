
/*
  The following file contains the operations

  C Op(=|+=|-=) scale * Or(A) * Or(B)

  C Op(=|+=|-=) scale * Or(A) * Or(A)

  scale =
  Op = "" = assignment, +=
  Or = "" = normal, "Trans" = transpose, "Symm" = symmetric
*/

inline void Symm3x3Symm( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
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

inline void Symm3x3SymmScale( TacsScalar scale, const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
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

inline void Symm3x3SymmAdd( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
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

inline void Symm3x3SymmSub( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
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

inline void Symm3x3SymmAddScale( TacsScalar scale, const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
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

inline void Symm3x3Mat( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
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

inline void Symm3x3MatScale( TacsScalar scale, const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
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

inline void Symm3x3MatAdd( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
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

inline void Symm3x3MatSub( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
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

inline void Symm3x3MatAddScale( TacsScalar scale, const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
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

inline void Symm3x3MatTrans( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
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

inline void Symm3x3MatTransScale( TacsScalar scale, const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
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

inline void Symm3x3MatTransAdd( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
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

inline void Symm3x3MatTransSub( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
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

inline void Symm3x3MatTransAddScale( TacsScalar scale, const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
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

inline void Mat3x3Symm( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
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

inline void Mat3x3SymmScale( TacsScalar scale, const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
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

inline void Mat3x3SymmAdd( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
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

inline void Mat3x3SymmSub( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
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

inline void Mat3x3SymmAddScale( TacsScalar scale, const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
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

inline void MatTrans3x3Symm( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
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

inline void MatTrans3x3SymmScale( TacsScalar scale, const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
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

inline void MatTrans3x3SymmAdd( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
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

inline void MatTrans3x3SymmSub( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
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

inline void MatTrans3x3SymmAddScale( TacsScalar scale, const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
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

inline void Mat3x3Mat( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
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

inline void Mat3x3MatScale( TacsScalar scale, const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
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

inline void Mat3x3MatAdd( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
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

inline void Mat3x3MatSub( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
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

inline void Mat3x3MatAddScale( TacsScalar scale, const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
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

inline void Mat3x3MatTrans( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
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

inline void Mat3x3MatTransScale( TacsScalar scale, const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
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

inline void Mat3x3MatTransAdd( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
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

inline void Mat3x3MatTransSub( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
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

inline void Mat3x3MatTransAddScale( TacsScalar scale, const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
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

inline void MatTrans3x3Mat( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
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

inline void MatTrans3x3MatScale( TacsScalar scale, const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
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

inline void MatTrans3x3MatAdd( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
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

inline void MatTrans3x3MatSub( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
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

inline void MatTrans3x3MatAddScale( TacsScalar scale, const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
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

inline void MatTrans3x3MatTrans( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
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

inline void MatTrans3x3MatTransScale( TacsScalar scale, const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
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

inline void MatTrans3x3MatTransAdd( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
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

inline void MatTrans3x3MatTransSub( const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
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

inline void MatTrans3x3MatTransAddScale( TacsScalar scale, const TacsScalar A[], const TacsScalar B[], TacsScalar C[] ){
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
