#include <complex>

/*
  Use the cplx type for TacsComplex
*/
typedef std::complex<double> TacsComplex;
typedef double TacsReal;

/*
  Define the basic scalar type TacsScalar
*/
#ifdef TACS_USE_COMPLEX
typedef TacsComplex TacsScalar;
#else
typedef TacsReal TacsScalar;
#endif


/*
  AD operator class
*/
class ADOp { // AD Operator base class
public:
  ADOp(){}
  virtual void computeDeriv(){}
};

/*
  AD object class
*/
class ADObj {
public:
  ADObj(){ op = NULL; }
  ADOp *op; // Operator that produced this object (if any)
};

/*
  Scalar type
*/
class ADScalar : public ADObj {
public:
  ADScalar(){
    adjoint = 0.0;
  }
  ADScalar( const ADScalar& a ){
    value = a.value;
    adjoint = a.adjoint;
  }

  TacsScalar value;
  TacsScalar adjoint;
};

/*
  General 3x3 matrix class
*/
class ADMat3x3 : public ADObj {
public:
  ADMat3x3(){
    for ( int i = 0; i < 9; i++ ){
      Ad[i] = 0.0;
    }
  }
  ADMat3x3( const TacsScalar a[] ){
    for ( int i = 0; i < 9; i++ ){
      A[i] = a[i];
      Ad[i] = 0.0;
    }
  }
  ADMat3x3( const ADMat3x3& a ){
    for ( int i = 0; i < 9; i++ ){
      A[i] = a.A[i];
      Ad[i] = a.Ad[i];
    }
  }

  TacsScalar A[9], Ad[9];
};