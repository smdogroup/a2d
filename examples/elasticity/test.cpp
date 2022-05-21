#include "a2dtmp.h"
#include <iostream>
#include <iomanip>
#include <complex>

template <class A>
class A2DExpression {
public:
  A& cast(){ return static_cast<A&>(*this); }
  const A& cast() const { return static_cast<const A&>(*this); }

  void forward(){ cast().forward(); }
  void reverse(){ cast().reverse(); }

  void hforward(){ cast().hforward(); }
  void hproduct(){ cast().hproduct(); }
  void hreverse(){ cast().hreverse(); }
};

template<class ScalarType>
class A2DScalarType {
public:
  A2DScalarType( ScalarType& value, ScalarType& bvalue,
                 ScalarType& pvalue, ScalarType& hvalue ) :
                 value(value), bvalue(bvalue), pvalue(pvalue), hvalue(hvalue) {}

  ScalarType& value;
  ScalarType& bvalue;
  ScalarType& pvalue;
  ScalarType& hvalue;
};

template<class MatType>
class A2DMat {
public:
  A2DMat( MatType *A_, MatType *Ab_=NULL, MatType *Ad2_=NULL, MatType *Ab2_=NULL ){
    A = A_;
    Ab = Ab_;
    Ad2 = Ad2_;
    Ab2 = Ab2_;
  }

  MatType& value(){ return *A; }
  const MatType& value() const { return *A; }

  MatType& bvalue(){ return *Ab; }
  const MatType& bvalue() const { return *Ab; }

  MatType& pvalue(){ return *Ad2; }
  const MatType& pvalue() const { return *Ad2; }

  MatType& hvalue(){ return *Ab2; }
  const MatType& hvalue() const { return *Ab2; }

  MatType* A; // Matrix
  MatType* Ab; // Reverse mode derivative value
  MatType *Ad2; // Projected second derivative value
  MatType *Ab2; // Reverse mode second derivative
};


template<class SMatType, class EMatType, class ScalarType>
class Symm3x3SymmMultTrace : public A2DExpression<Symm3x3SymmMultTrace<SMatType, EMatType, ScalarType> > {
public:
  Symm3x3SymmMultTrace( A2DMat<SMatType>& SObj, A2DMat<EMatType>& EObj, A2DScalarType<ScalarType>& output ) :
    SObj(SObj), EObj(EObj), output(output) {
    const SMatType& S = SObj.value();
    const EMatType& E = EObj.value();

    output.value =
      S(0, 0) * E(0, 0) + S(1, 1) * E(1, 1) + S(2, 2) * E(2, 2) +
      2.0 * (S(0, 1) * E(0, 1) + S(0, 2) * E(0, 2) + S(1, 2) * E(1, 2));
  }

  void reverse(){
    const EMatType& E = EObj.value();
    EMatType& Eb = EObj.bvalue();
    const SMatType& S = SObj.value();
    SMatType& Sb = SObj.bvalue();

    Eb(0, 0) += output.bvalue * S(0, 0);
    Eb(1, 1) += output.bvalue * S(1, 1);
    Eb(2, 2) += output.bvalue * S(2, 2);
    Eb(0, 1) += 2.0 * output.bvalue * S(0, 1);
    Eb(0, 2) += 2.0 * output.bvalue * S(0, 2);
    Eb(1, 2) += 2.0 * output.bvalue * S(1, 2);

    Sb(0, 0) += output.bvalue * E(0, 0);
    Sb(1, 1) += output.bvalue * E(1, 1);
    Sb(2, 2) += output.bvalue * E(2, 2);
    Sb(0, 1) += 2.0 * output.bvalue * E(0, 1);
    Sb(0, 2) += 2.0 * output.bvalue * E(0, 2);
    Sb(1, 2) += 2.0 * output.bvalue * E(1, 2);
  }

  // Compute E.pvalue() = J * Ux.pvalue()
  void hforward(){
    const EMatType& E = EObj.value();
    const EMatType& Ed = EObj.pvalue();
    const SMatType& S = SObj.value();
    const SMatType& Sd = SObj.pvalue();

    output.pvalue =
      S(0, 0) * Ed(0, 0) + S(1, 1) * Ed(1, 1) + S(2, 2) * Ed(2, 2) +
      2.0 * (S(0, 1) * Ed(0, 1) + S(0, 2) * Ed(0, 2) + S(1, 2) * Ed(1, 2)) +
      Sd(0, 0) * E(0, 0) + Sd(1, 1) * E(1, 1) + Sd(2, 2) * E(2, 2) +
      2.0 * (Sd(0, 1) * E(0, 1) + Sd(0, 2) * E(0, 2) + Sd(1, 2) * E(1, 2));
  }

  // Compute df/d(trace) * d^2(trace) * (dS/dp, dE/dp)
  void hproduct(){
    const EMatType& Ed = EObj.pvalue();
    EMatType& Eh = EObj.hvalue();
    const SMatType& Sd = SObj.pvalue();
    SMatType& Sh = SObj.hvalue();

    Eh(0, 0) += output.bvalue * Sd(0, 0);
    Eh(1, 1) += output.bvalue * Sd(1, 1);
    Eh(2, 2) += output.bvalue * Sd(2, 2);
    Eh(0, 1) += 2.0 * output.bvalue * Sd(0, 1);
    Eh(0, 2) += 2.0 * output.bvalue * Sd(0, 2);
    Eh(1, 2) += 2.0 * output.bvalue * Sd(1, 2);

    Sh(0, 0) += output.bvalue * Ed(0, 0);
    Sh(1, 1) += output.bvalue * Ed(1, 1);
    Sh(2, 2) += output.bvalue * Ed(2, 2);
    Sh(0, 1) += 2.0 * output.bvalue * Ed(0, 1);
    Sh(0, 2) += 2.0 * output.bvalue * Ed(0, 2);
    Sh(1, 2) += 2.0 * output.bvalue * Ed(1, 2);
  }

  void hreverse(){
    const EMatType& Ed = EObj.pvalue();
    EMatType& Eb = EObj.hvalue();
    const SMatType& Sd = SObj.pvalue();
    SMatType& Sb = SObj.hvalue();

    Eb(0, 0) += output.hvalue * Sd(0, 0);
    Eb(1, 1) += output.hvalue * Sd(1, 1);
    Eb(2, 2) += output.hvalue * Sd(2, 2);
    Eb(0, 1) += 2.0 * output.hvalue * Sd(0, 1);
    Eb(0, 2) += 2.0 * output.hvalue * Sd(0, 2);
    Eb(1, 2) += 2.0 * output.hvalue * Sd(1, 2);

    Sb(0, 0) += output.hvalue * Ed(0, 0);
    Sb(1, 1) += output.hvalue * Ed(1, 1);
    Sb(2, 2) += output.hvalue * Ed(2, 2);
    Sb(0, 1) += 2.0 * output.hvalue * Ed(0, 1);
    Sb(0, 2) += 2.0 * output.hvalue * Ed(0, 2);
    Sb(1, 2) += 2.0 * output.hvalue * Ed(1, 2);
  }

  A2DMat<SMatType>& SObj;
  A2DMat<EMatType>& EObj;
  A2DScalarType<ScalarType>& output;
};

template<class ScalarType, class EMatType, class SMatType>
class Symm3x3IsotropicConstitutive : public A2DExpression<Symm3x3IsotropicConstitutive<ScalarType, EMatType, SMatType> >{
public:
  Symm3x3IsotropicConstitutive( const ScalarType& mu, const ScalarType& lambda,
                                A2DMat<EMatType>& EObj, A2DMat<SMatType>& SObj ) : mu(mu), lambda(lambda), EObj(EObj), SObj(SObj) {
    const EMatType& E = EObj.value();
    SMatType& S = SObj.value();
    ScalarType tr = lambda * (E(0, 0) + E(1, 1) + E(2, 2));
    ScalarType mu2 = 2.0 * mu;
    S(0, 0) = mu2 * E(0, 0) + tr;
    S(0, 1) = mu2 * E(0, 1);
    S(0, 2) = mu2 * E(0, 2);
    S(1, 1) = mu2 * E(1, 1) + tr;
    S(1, 2) = mu2 * E(1, 2);
    S(2, 2) = mu2 * E(2, 2) + tr;
  }

  void reverse(){
    const SMatType& Sb = SObj.bvalue();
    EMatType& Eb = EObj.bvalue();

    ScalarType tr = lambda * (Sb(0, 0) + Sb(1, 1) + Sb(2, 2));
    ScalarType mu2 = 2.0 * mu;
    Eb(0, 0) += mu2 * Sb(0, 0) + tr;
    Eb(0, 1) += mu2 * Sb(0, 1);
    Eb(0, 2) += mu2 * Sb(0, 2);
    Eb(1, 1) += mu2 * Sb(1, 1) + tr;
    Eb(1, 2) += mu2 * Sb(1, 2);
    Eb(2, 2) += mu2 * Sb(2, 2) + tr;
  }

  void hforward(){
    const EMatType& Ed = EObj.pvalue();
    SMatType& Sd = SObj.pvalue();

    ScalarType tr = lambda * (Ed(0, 0) + Ed(1, 1) + Ed(2, 2));
    ScalarType mu2 = 2.0 * mu;
    Sd(0, 0) = mu2 * Ed(0, 0) + tr;
    Sd(0, 1) = mu2 * Ed(0, 1);
    Sd(0, 2) = mu2 * Ed(0, 2);
    Sd(1, 1) = mu2 * Ed(1, 1) + tr;
    Sd(1, 2) = mu2 * Ed(1, 2);
    Sd(2, 2) = mu2 * Ed(2, 2) + tr;
  }

  void hproduct(){}

  void hreverse(){
    const SMatType& Sh = SObj.hvalue();
    EMatType& Eh = EObj.hvalue();

    ScalarType tr = lambda * (Sh(0, 0) + Sh(1, 1) + Sh(2, 2));
    ScalarType mu2 = 2.0 * mu;
    Eh(0, 0) += mu2 * Sh(0, 0) + tr;
    Eh(0, 1) += mu2 * Sh(0, 1);
    Eh(0, 2) += mu2 * Sh(0, 2);
    Eh(1, 1) += mu2 * Sh(1, 1) + tr;
    Eh(1, 2) += mu2 * Sh(1, 2);
    Eh(2, 2) += mu2 * Sh(2, 2) + tr;
  }

  const ScalarType& mu;
  const ScalarType& lambda;
  A2DMat<EMatType>& EObj;
  A2DMat<SMatType>& SObj;
};

template<class UxMatType, class EMatType>
class Mat3x3GreenStrain : public A2DExpression<Mat3x3GreenStrain<UxMatType, EMatType> > {
public:
  Mat3x3GreenStrain( A2DMat<UxMatType>& UxObj, A2DMat<EMatType>& EObj ) : UxObj(UxObj), EObj(EObj) {
    const UxMatType& Ux = UxObj.value();
    EMatType& E = EObj.value();
    E(0, 0) = Ux(0, 0) + 0.5*(Ux(0, 0) * Ux(0, 0) + Ux(1, 0) * Ux(1, 0) + Ux(2, 0) * Ux(2, 0));
    E(1, 1) = Ux(1, 1) + 0.5*(Ux(0, 1) * Ux(0, 1) + Ux(1, 1) * Ux(1, 1) + Ux(2, 1) * Ux(2, 1));
    E(2, 2) = Ux(2, 2) + 0.5*(Ux(0, 2) * Ux(0, 2) + Ux(1, 2) * Ux(1, 2) + Ux(2, 2) * Ux(2, 2));

    E(0, 1) = 0.5*(Ux(0, 1) + Ux(1, 0) + Ux(0, 0) * Ux(0, 1) + Ux(1, 0) * Ux(1, 1) + Ux(2, 0) * Ux(2, 1));
    E(0, 2) = 0.5*(Ux(0, 2) + Ux(2, 0) + Ux(0, 0) * Ux(0, 2) + Ux(1, 0) * Ux(1, 2) + Ux(2, 0) * Ux(2, 2));
    E(1, 2) = 0.5*(Ux(1, 2) + Ux(2, 1) + Ux(0, 1) * Ux(0, 2) + Ux(1, 1) * Ux(1, 2) + Ux(2, 1) * Ux(2, 2));
  }

  void reverse(){
    const UxMatType& Ux = UxObj.value();
    const EMatType& Eb = EObj.bvalue();
    UxMatType& Uxb = UxObj.bvalue();

    // Uxb = (I + Ux) * Eb
    Uxb(0, 0) +=       (Ux(0, 0) + 1.0) * Eb(0, 0) + 0.5 * Ux(0, 1) * Eb(0, 1) + 0.5 * Ux(0, 2) * Eb(0, 2);
    Uxb(0, 1) += 0.5 * (Ux(0, 0) + 1.0) * Eb(0, 1) +       Ux(0, 1) * Eb(1, 1) + 0.5 * Ux(0, 2) * Eb(1, 2);
    Uxb(0, 2) += 0.5 * (Ux(0, 0) + 1.0) * Eb(0, 2) + 0.5 * Ux(0, 1) * Eb(1, 2) +       Ux(0, 2) * Eb(2, 2);

    Uxb(1, 0) +=       Ux(1, 0) * Eb(0, 0) + 0.5 * (Ux(1, 1) + 1.0) * Eb(0, 1) + 0.5 * Ux(1, 2) * Eb(0, 2);
    Uxb(1, 1) += 0.5 * Ux(1, 0) * Eb(0, 1) +       (Ux(1, 1) + 1.0) * Eb(1, 1) + 0.5 * Ux(1, 2) * Eb(1, 2);
    Uxb(1, 2) += 0.5 * Ux(1, 0) * Eb(0, 2) + 0.5 * (Ux(1, 1) + 1.0) * Eb(1, 2) +       Ux(1, 2) * Eb(2, 2);

    Uxb(2, 0) +=       Ux(2, 0) * Eb(0, 0) + 0.5 * Ux(2, 1) * Eb(0, 1) + 0.5 * (Ux(2, 2) + 1.0) * Eb(0, 2);
    Uxb(2, 1) += 0.5 * Ux(2, 0) * Eb(0, 1) +       Ux(2, 1) * Eb(1, 1) + 0.5 * (Ux(2, 2) + 1.0) * Eb(1, 2);
    Uxb(2, 2) += 0.5 * Ux(2, 0) * Eb(0, 2) + 0.5 * Ux(2, 1) * Eb(1, 2) +       (Ux(2, 2) + 1.0) * Eb(2, 2);
  }

  void hforward(){
    const UxMatType& Ux = UxObj.value();
    const UxMatType& Uxd = UxObj.pvalue();
    EMatType& Ed = EObj.pvalue();

    Ed(0, 0) = Uxd(0, 0) + Ux(0, 0) * Uxd(0, 0) + Ux(1, 0) * Uxd(1, 0) + Ux(2, 0) * Uxd(2, 0);
    Ed(1, 1) = Uxd(1, 1) + Ux(0, 1) * Uxd(0, 1) + Ux(1, 1) * Uxd(1, 1) + Ux(2, 1) * Uxd(2, 1);
    Ed(2, 2) = Uxd(2, 2) + Ux(0, 2) * Uxd(0, 2) + Ux(1, 2) * Uxd(1, 2) + Ux(2, 2) * Uxd(2, 2);

    Ed(0, 1) = 0.5*(Uxd(0, 1) + Uxd(1, 0) +
                    Ux(0, 0) * Uxd(0, 1) + Ux(1, 0) * Uxd(1, 1) + Ux(2, 0) * Uxd(2, 1) +
                    Uxd(0, 0) * Ux(0, 1) + Uxd(1, 0) * Ux(1, 1) + Uxd(2, 0) * Ux(2, 1));
    Ed(0, 2) = 0.5*(Uxd(0, 2) + Uxd(2, 0) +
                    Ux(0, 0) * Uxd(0, 2) + Ux(1, 0) * Uxd(1, 2) + Ux(2, 0) * Uxd(2, 2) +
                    Uxd(0, 0) * Ux(0, 2) + Uxd(1, 0) * Ux(1, 2) + Uxd(2, 0) * Ux(2, 2));
    Ed(1, 2) = 0.5*(Uxd(1, 2) + Uxd(2, 1) +
                    Ux(0, 1) * Uxd(0, 2) + Ux(1, 1) * Uxd(1, 2) + Ux(2, 1) * Uxd(2, 2) +
                    Uxd(0, 1) * Ux(0, 2) + Uxd(1, 1) * Ux(1, 2) + Uxd(2, 1) * Ux(2, 2));
  }

  void hproduct(){
    const UxMatType& Eb = EObj.bvalue();
    const UxMatType& Uxd = UxObj.pvalue();
    UxMatType& Uxh = UxObj.hvalue();

    Uxh(0, 0) +=       Uxd(0, 0) * Eb(0, 0) + 0.5 * Uxd(0, 1) * Eb(0, 1) + 0.5 * Uxd(0, 2) * Eb(0, 2);
    Uxh(0, 1) += 0.5 * Uxd(0, 0) * Eb(0, 1) +       Uxd(0, 1) * Eb(1, 1) + 0.5 * Uxd(0, 2) * Eb(1, 2);
    Uxh(0, 2) += 0.5 * Uxd(0, 0) * Eb(0, 2) + 0.5 * Uxd(0, 1) * Eb(1, 2) +       Uxd(0, 2) * Eb(2, 2);

    Uxh(1, 0) +=       Uxd(1, 0) * Eb(0, 0) + 0.5 * Uxd(1, 1) * Eb(0, 1) + 0.5 * Uxd(1, 2) * Eb(0, 2);
    Uxh(1, 1) += 0.5 * Uxd(1, 0) * Eb(0, 1) +       Uxd(1, 1) * Eb(1, 1) + 0.5 * Uxd(1, 2) * Eb(1, 2);
    Uxh(1, 2) += 0.5 * Uxd(1, 0) * Eb(0, 2) + 0.5 * Uxd(1, 1) * Eb(1, 2) +       Uxd(1, 2) * Eb(2, 2);

    Uxh(2, 0) +=       Uxd(2, 0) * Eb(0, 0) + 0.5 * Uxd(2, 1) * Eb(0, 1) + 0.5 * Uxd(2, 2) * Eb(0, 2);
    Uxh(2, 1) += 0.5 * Uxd(2, 0) * Eb(0, 1) +       Uxd(2, 1) * Eb(1, 1) + 0.5 * Uxd(2, 2) * Eb(1, 2);
    Uxh(2, 2) += 0.5 * Uxd(2, 0) * Eb(0, 2) + 0.5 * Uxd(2, 1) * Eb(1, 2) +       Uxd(2, 2) * Eb(2, 2);
  }

  void hreverse(){
    const UxMatType& Ux = UxObj.value();
    const EMatType& Eb = EObj.hvalue();
    UxMatType& Uxb = UxObj.hvalue();

    // Uxb = (I + Ux) * Eb
    Uxb(0, 0) +=       (Ux(0, 0) + 1.0) * Eb(0, 0) + 0.5 * Ux(0, 1) * Eb(0, 1) + 0.5 * Ux(0, 2) * Eb(0, 2);
    Uxb(0, 1) += 0.5 * (Ux(0, 0) + 1.0) * Eb(0, 1) +       Ux(0, 1) * Eb(1, 1) + 0.5 * Ux(0, 2) * Eb(1, 2);
    Uxb(0, 2) += 0.5 * (Ux(0, 0) + 1.0) * Eb(0, 2) + 0.5 * Ux(0, 1) * Eb(1, 2) +       Ux(0, 2) * Eb(2, 2);

    Uxb(1, 0) +=       Ux(1, 0) * Eb(0, 0) + 0.5 * (Ux(1, 1) + 1.0) * Eb(0, 1) + 0.5 * Ux(1, 2) * Eb(0, 2);
    Uxb(1, 1) += 0.5 * Ux(1, 0) * Eb(0, 1) +       (Ux(1, 1) + 1.0) * Eb(1, 1) + 0.5 * Ux(1, 2) * Eb(1, 2);
    Uxb(1, 2) += 0.5 * Ux(1, 0) * Eb(0, 2) + 0.5 * (Ux(1, 1) + 1.0) * Eb(1, 2) +       Ux(1, 2) * Eb(2, 2);

    Uxb(2, 0) +=       Ux(2, 0) * Eb(0, 0) + 0.5 * Ux(2, 1) * Eb(0, 1) + 0.5 * (Ux(2, 2) + 1.0) * Eb(0, 2);
    Uxb(2, 1) += 0.5 * Ux(2, 0) * Eb(0, 1) +       Ux(2, 1) * Eb(1, 1) + 0.5 * (Ux(2, 2) + 1.0) * Eb(1, 2);
    Uxb(2, 2) += 0.5 * Ux(2, 0) * Eb(0, 2) + 0.5 * Ux(2, 1) * Eb(1, 2) +       (Ux(2, 2) + 1.0) * Eb(2, 2);
  }

  A2DMat<UxMatType>& UxObj;
  A2DMat<EMatType>& EObj;
};

int main( int argc, char *argv[] ){
  // typedef int32_t IndexType;
  typedef std::complex<double> ScalarType;
  typedef A2D::Mat<ScalarType, 3, 3> Mat3x3;
  typedef A2D::SymmMat<ScalarType, 3> SymmMat3x3;
  typedef A2D::Mat2ndDeriv<ScalarType, 3, 3> Mat2ndDeriv;

  double dh = 1e-30;

  Mat3x3 Ux, Uxb, Uxd, Uxh;
  SymmMat3x3 E, Eb, Ed, Eh;
  SymmMat3x3 S, Sb, Sd, Sh;
  ScalarType outputval, outputb = 0.0, outputd = 0.0, outputh = 0.0;

  A2DScalarType<ScalarType> output(outputval, outputb, outputd, outputh);
  A2DMat<Mat3x3> UxObj(&Ux, &Uxb, &Uxd, &Uxh);
  A2DMat<SymmMat3x3> EObj(&E, &Eb, &Ed, &Eh);
  A2DMat<SymmMat3x3> SObj(&S, &Sb, &Sd, &Sh);

  ScalarType mu(0.2533), lambda(0.71236);

  // Set random values for Ux
  for ( int i = 0; i < 3; i++ ){
    for ( int j = 0; j < 3; j++ ){
      Ux(i, j) = -1.0 + 2.0 * rand()/RAND_MAX;
      Uxd(i, j) = -1.0 + 2.0 * rand()/RAND_MAX;
    }
  }

  for ( int i = 0; i < 3; i++ ){
    for ( int j = 0; j < 3; j++ ){
      Ux(i, j) = Ux(i, j) + ScalarType(0.0, dh) * Uxd(i, j);
    }
  }

  Mat3x3GreenStrain<Mat3x3, SymmMat3x3> strain(UxObj, EObj);
  Symm3x3IsotropicConstitutive<ScalarType, SymmMat3x3, SymmMat3x3> constitutive(mu, lambda, EObj, SObj);
  Symm3x3SymmMultTrace<SymmMat3x3, SymmMat3x3, ScalarType> trace(SObj, EObj, output);

  output.bvalue = 1.0;

  trace.reverse();
  constitutive.reverse();
  strain.reverse();

  ScalarType result = 0.0;
  for ( int i = 0; i < 3; i++ ){
    for ( int j = 0; j < 3; j++ ){
      result += Uxb(i, j) * Uxd(i, j);
    }
  }

  // double forward = output.valueb.real();
  double res = result.real();
  double fd = output.value.imag()/dh;
  double error = (res - fd)/fd;
  std::cout << "result: " << std::setw(20) << res << " fd: " << std::setw(20) << fd
    << " error: " << std::setw(20) << error << std::endl;

  strain.hforward();
  constitutive.hforward();
  trace.hproduct();
  constitutive.hreverse();
  strain.hreverse();

  strain.hproduct();

  for ( int i = 0; i < 3; i++ ){
    for ( int j = 0; j < 3; j++ ){
      double res = Uxh(i, j).real();
      double fd = Uxb(i, j).imag()/dh;
      double error = (res - fd)/fd;

      std::cout << i << ", " << j << " result: " << std::setw(20) << res <<
        " fd: " << std::setw(20) << fd << " error: " << std::setw(20) <<  error << std::endl;
    }
  }

  // // Compute the second derivative
  // Mat2ndDeriv Ux2d;

  // // Compute the second derivatives
  // for ( int i = 0; i < 3; i++ ){
  //   for ( int j = 0; j < 3; j++ ){
  //     Mat3x3 Uxp, Uxb2;
  //     SymmMat3x3 Ep, Sp, Sb2, Eb2;

  //     Uxp.zero();
  //     Uxb2.zero();
  //     Sb2.zero();
  //     Eb2.zero();
  //     Uxp(i, j) = 1.0;

  //     strain.hforward();
  //     strain.hproduct();

  //     constitutive.hforward();
  //     constitutive.hproduct();

  //     trace.h
  //     trace.hessian_vec_product(outputb, Sp, Ep, Sb2, Eb2);
  //     constitutive.reverse(Sb2, Eb2);
  //     strain.reverse(Eb2, Uxb2);

  //     for ( int k = 0; k < 3; k++ ){
  //       for ( int l = 0; l < 3; l++ ){
  //         Ux2d(k, l, i, j) = Uxb2(k, l);
  //       }
  //     }
  //   }
  // }

  // Mat3x3 prod;
  // prod.zero();
  // for ( int i = 0; i < 3; i++ ){
  //   for ( int j = 0; j < 3; j++ ){
  //     for ( int k = 0; k < 3; k++ ){
  //       for ( int l = 0; l < 3; l++ ){
  //         prod(i, j) += Ux2d(i, j, k, l) * dUx(k, l);
  //       }
  //     }
  //   }
  // }

  // for ( int i = 0; i < 3; i++ ){
  //   for ( int j = 0; j < 3; j++ ){
  //     double res = prod(i, j).real();
  //     double fd = Uxb(i, j).imag()/dh;
  //     double error = (res - fd)/fd;

  //     std::cout << i << ", " << j << " result: " << std::setw(20) << res <<
  //       " fd: " << std::setw(20) << fd << " error: " << std::setw(20) <<  error << std::endl;
  //   }
  // }

  return (0);
}