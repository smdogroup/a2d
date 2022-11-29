template <typename T, A2D::index_t D>
class NonlinearElasticity {
 public:
  NonlinearElasticity(T mu, T lambda) : mu(mu), lambda(lambda) {}

  // Number of dimensions
  static const A2D::index_t dim = D;

  // Number of data dimensions
  static const A2D::index_t data_dim = 2;

  // Space for the finite-element data
  typedef A2D::FESpace<T, data_dim, A2D::H1Space<T, data_dim, dim>> DataSpace;

  // Space for the element geometry
  typedef A2D::FESpace<T, dim, A2D::H1Space<T, dim, dim>> FiniteElementGeometry;

  // Finite element space
  typedef A2D::FESpace<T, dim, A2D::H1Space<T, dim, dim>> FiniteElementSpace;

  // The type of matrix used to store data at each quadrature point
  typedef A2D::SymmMat<T, FiniteElementSpace::ncomp> QMatType;

  /**
   * @brief Evaluate the weak form coefficients for nonlinear elasticity
   *
   * @param wdetJ The quadrature weight times determinant of the Jacobian
   * @param data The data at the quadrature point
   * @param s The trial solution
   * @param coef Output weak form coefficients of the test space
   */
  A2D_INLINE_FUNCTION void weak(T wdetJ, const DataSpace& data,
                                const FiniteElementGeometry& geo,
                                const FiniteElementSpace& s,
                                FiniteElementSpace& coef) {
    // Get the constitutive data at the points
    T mu = data[0];
    T lambda = data[1];

    // Extract the trial solution gradient and the coefficient terms. Here
    // Uxb is the output computed as the derivative of the strain energy
    // w.r.t. Ux
    A2D::Mat<T, dim, dim> Ux0 = (s.template get<0>()).get_grad();
    A2D::Mat<T, dim, dim>& Uxb = (coef.template get<0>()).get_grad();
    A2D::ADMat<A2D::Mat<T, dim, dim>> Ux(Ux0, Uxb);

    // The Green-Langrange strain terms
    A2D::SymmMat<T, dim> E0, Eb;
    A2D::ADMat<A2D::SymmMat<T, dim>> E(E0, Eb);

    // The strain energy output
    A2D::ADScalar<T> output;

    auto strain = A2D::MatGreenStrain(Ux, E);
    auto energy = A2D::SymmIsotropicEnergy(mu, lambda, E, output);

    // Seed the output value with the wdetJ
    output.bvalue = wdetJ;

    // Reverse the derivatives through the code
    energy.reverse();
    strain.reverse();
  }

  /**
   * @brief Construct a JacVecProduct functor
   *
   * This functor computes a Jacobian-vector product of the
   *
   * @param wdetJ The quadrature weight times determinant of the Jacobian
   * @param data The data at the quadrature point
   * @param s The solution at the quadrature point
   */
  class JacVecProduct {
   public:
    A2D_INLINE_FUNCTION JacVecProduct(const NonlinearElasticity<T, D>& pde,
                                      T wdetJ, const DataSpace& data,
                                      const FiniteElementGeometry& geo,
                                      const FiniteElementSpace& s)
        :  // Initialize constitutive data
          mu(data[0]),
          lambda(data[1]),

          // Initialize the displacement gradient
          Ux(s.template get<0>().get_grad()),

          // Compute the strain from the displacement gradient
          strain(Ux, E),

          // Compute the strain energy from the strain
          energy(mu, lambda, E, output) {
      // Set the seed on the derivative
      output.bvalue = wdetJ;

      // Reverse mode for the first derivative
      energy.reverse();
      strain.reverse();
    }

    A2D_INLINE_FUNCTION void operator()(const FiniteElementSpace& p,
                                        FiniteElementSpace& Jp) {
      Ux.set_pvalue((p.template get<0>()).get_grad());

      strain.hforward();
      energy.hreverse();
      strain.hreverse();

      Ux.get_hvalue((Jp.template get<0>()).get_grad());
    }

   private:
    T mu, lambda;
    A2D::A2DMat<A2D::Mat<T, dim, dim>> Ux;
    A2D::A2DMat<A2D::SymmMat<T, dim>> E;
    A2D::A2DScalar<T> output;

    // Declare types of the operators
    decltype(A2D::MatGreenStrain(Ux, E)) strain;
    decltype(A2D::SymmIsotropicEnergy(mu, lambda, E, output)) energy;
  };

  /**
   * @brief Construct a JacVecProduct functor
   *
   * @param wdetJ The quadrature weight times determinant of the Jacobian
   * @param data The data at the quadrature point
   * @param s The solution at the quadrature point
   */
  class AdjVecProduct {
   public:
    A2D_INLINE_FUNCTION AdjVecProduct(const NonlinearElasticity<T, D>& pde,
                                      T wdetJ, const DataSpace& data,
                                      const FiniteElementGeometry& geo,
                                      const FiniteElementSpace& s)
        :  // Initialize constitutive data
          mu(data[0]),
          lambda(data[1]),

          // Initialize the displacement gradient
          Ux(s.template get<0>().get_grad()),

          // Compute the strain from the displacement gradient
          strain(Ux, E),

          // Compute the strain energy from the strain
          energy(mu, lambda, E, output) {
      // Set the seed on the derivative
      output.bvalue = wdetJ;

      // Reverse mode for the first derivative
      energy.reverse();
      strain.reverse();
    }

    A2D_INLINE_FUNCTION void operator()(const FiniteElementSpace& p,
                                        DataSpace& dfdx) {
      Ux.set_pvalue((p.template get<0>()).get_grad());

      strain.hforward();
      energy.hreverse();
      strain.hreverse();

      dfdx.set_value(0, mu.hvalue);
      dfdx.set_value(1, lambda.hvalue);
    }

   private:
    A2D::A2DScalar<T> mu, lambda;
    A2D::A2DMat<A2D::Mat<T, dim, dim>> Ux;
    A2D::A2DMat<A2D::SymmMat<T, dim>> E;
    A2D::A2DScalar<T> output;

    // Declare types of the operators
    decltype(A2D::MatGreenStrain(Ux, E)) strain;
    decltype(A2D::SymmIsotropicEnergy(mu, lambda, E, output)) energy;
  };
};