# Function space, element and basis
Solution (as well as test) functions involved in the PDEs has different
- spatial dimention ```D``` (2-dimensional, 3-dimensional, etc.),
- variable dimensin ```C``` (```C=1``` for scalar variable, ```C>1``` for
vector variable).

More importantly, each of them falls in one of the following function spaces:
- L2 space: (scalar or vector) values only
- H1 space: (scalar or vector) values and 1st order derivatives
- H(div) space: vector values and divergence
- H(curl) space: vector values and curl

Abstract function spaces above are implemented in ```multiphysics/fespace.h```.
Currently, supported function spaces are:
- ```L2Space<T, C, D>```: L2 space for ```C```-dimensional vector variable with
spatial dimension ```D```
- ```H1Space<T, C, D>```: H1 space for ```C```-dimensional vector variable with
spatial dimension ```D```
- ```Hdiv2DSpace<T>```: H1 space for 2-component vector variable with spatial
dimensional 2
- ```Hcurl2DSpace<T>```: H1 space for 2-component vector variable with spatial
dimensional 2

Numerical approximations of abstract function spaces are implemented in
```multiphysics/febasis.h```.
Each basis class corresponds to exact one type of the function space.
Currently, supported bases are:
- 2-dimensional:
    - ```LagrangeTri0<T, C>```: associated with ```L2Space<T, C, 2>```, degree-0
    Lagrange triangle
    - ```LagrangeTri1<T, C>```: associated with ```H1Space<T, C, 2>```, degree-1
    Lagrange triangle
    - ```RT2DTri1<T>```: associated with ```Hdiv2DSpace<T>```, degree-1
    Raviart–Thomas triangle


# Class template skeletons

## Function spaces: ```fespace.h```
Two types of class templates are defined:
- ```SpaceType```
- ```FESpace```

 ### ```SpaceType```
```SpaceType``` implements specific function space instances, such as
L2 space, H1 space, H(div) space, etc.
Such class templates comply with the following skeleton:
 ```c++
 class SpaceType {
    // The type of function variable associated with this function space
    using VarType = typename SpaceType::VarType;

    // number of components (e.g. number of solution variables and derivatives
    // combined)
    static const A2D::index_t ncomp;

    // spatial dimension (usually 2 or 3)
    static const A2D::index_t dim;

    // zero the variables (solution, grad, div, curl, egc.)
    void zero();

    /**
     * Get a component of the variables
     * @param comp component index
    */
    T& get_value(const A2D::index_t comp);
    const T& get_value(const A2D::index_t comp) const;

    /**
     * Get the solution variable (not including grad, div, etc.)
    */
    VarType& get_value();
    const VarType& get_value() const;


    /**
     * Transform solution variables from this (reference) element to physical
     * element s
     * @param s the function space associated with the physical element
    */
    void transform(const T& detJ, const A2D::Mat<T, D, D>& J,
                    const A2D::Mat<T, D, D>& Jinv, L2Space<T, C, D>& s) const;


    /**
     * Transform solution variables from this (physical) element to reference
     * element s
     * @param s the function space associated with the reference element
    */
    void rtransform(const T& detJ, const A2D::Mat<T, D, D>& J,
                    const A2D::Mat<T, D, D>& Jinv, L2Space<T, C, D>& s) const;


    /* Depend on the instance, could have methods below */
    grad_type& get_grad();  // get gradient
    grad_type& get_grad() const;
    div_type& get_div();  // get divergence
    div_type& get_div() const;
    curl_type& get_curl();  // get curl
    curl_type& get_curl() const;
 };
```

### ```FESpace```
```FESpace``` forms a collection of selected space class template ``instances''
for specific problems via variadic template argument.

## Element bases: ```febasis.h```
Two types of class templates are defined:
- ```BasisType```
- ```FEBasis```

### ```BasisType```
```BasisType``` implements specific element basis, such as Lagrange element,
Raviart-Thomas element, etc.
Such class templates comply with the following skeleton:
```c++
class BasisType {
    // The space type associated with this basis type
    using SpaceType = ...;

    // The type of function variable associated with this function space
    using VarType = typename SpaceType::VarType;

    // number of degrees of freedom of such element
    static const A2D::index_t ndof;

    // number of components (e.g. number of solution variables and derivatives
    // combined)
    static const A2D::index_t ncomp = SpaceType::ncomp;

    /**
     * Interpolate component(s) at a quadrature point.
     * @tparam offset number of dof combined for previous bases of this element
     * @param n quadrature index
     * @param sol element(potentially containing multiple bases)-level dof
     * @param out function space object that stores component(s) to be computed
    */
    template <class Quadrature, A2D::index_t offset, class SolnType>
    static void interp(A2D::index_t n, const SolnType sol,
                       SpaceType& out);

    /**
     * Add basis-level components to element-level dof
     * @tparam offset number of dof combined for previous bases of this element
     * @param n quadrature index
     * @param in function space object that stores component(s) to be used
     * @param res element(potentially containing multiple bases)-level dof
    */
    template <class Quadrature, A2D::index_t offset, class SolnType>
    static void add(A2D::index_t n, const SpaceType& in,
                    SolnType res);

    //
    static const A2D::index_t stride;

    // Number of basis functions, including gradient, div, etc. if supported by
    // the basis type
    static const A2D::index_t basis_size;

    // number of dof and component per stride
    static const A2D::index_t ndof_per_stride = ndof / stride;
    static const A2D::index_t ncomp_per_stride = ncomp / stride;

    /**
     * Evaluate basis functions (and derivatives, div, curl, etc.) at a quadrature point
     * @param n quadrature intex
     * @param N []-indexable array of basis function values
    */
    template <class Quadrature, class BasisType>
    static void basis(A2D::index_t n, BasisType N);
};
```

### ```FEBasis```
```FEBasis``` forms a collection of selected basis class template ``instances''
for specific problems via variadic template argument.

## Quadratures: ```fequadrature.h```
Quadrature schemes (numbering included) for different element types used for
numerical integration are defined here.
Each quadrature class template will be used as template argument in the static
sense, i.e. no object instantiation, and complies with the following skeleton:
```c++
class QuadratureType {
    static index_type get_num_points();  // return number of quadrature points
    static void get_point(const index_type n,
                          double pt[]);  // get n-th quadrature point in reference
                                         // coordinates
    static double get_weight(const index_type n);  // get weighting factor of n-th
                                                 // quadrature point
};
```


## PDEIntegrand implementations, for example: ```integrand_poisson.h```
A2D solves PDEs. Such problem instances need to be implemented to comply with
the following skeleton:
```c++
class PDEInstance {
    static const index_type dim;  // spatial dimension
    typename FiniteElementSpace;  // specialized FESpace for solution
    typename FiniteElementGeometry;  // specialized FESpace for geometry
    static T eval_weak_form(T wdetJ, FiniteElementSpace& s,
                            FiniteElementSpace& t);  // evaluate the scalar value
                                                     // of the weak form at a
                                                     // quadrature point given
                                                     // solution space s and
                                                     // test space t
    static void eval_weak_coef(T wdetJ, const FiniteElementSpace& s,
                               FiniteElementSpace& coef);  // evaluate the residuals
    static void eval_weak_jacobian_vec_product(T wdetJ,
                                               const FiniteElementSpace& s,
                                               const FiniteElementSpace& p,
                                               FiniteElementSpace& coef);
};
```

## Enabling access and manipulation of global dof on an element level: ```feelementvector.h```
This class implements interface for element-level computations.
Such interfaces are called ```ElementVector```, and must comply with the
following skeleton:
```c++
class ElementVector {

    ElementVector(...);  // takes in global solution vector and other needed
                             // information

    // A light-weight nested class that gains access to an element's local
    // degrees of freedom
    class FEDof {
        FEDof();  // must have empty argument list
        template<index_type index>
        square_bracket_indexable_type get();  // returns an object associated
                                              // with local dof, index is the
                                              // index of certain basis in the
                                              // basis collection FEBasis
        template<index_type index>
        const square_bracket_indexable_type get() const;  // const variant
    };

    // Following methods may be empty implementations depending on the specific
    // ElementVector implementation
    void get_element_values(index_type elem,
                            FEDof& dof);  // populate dof object, elem is element
                                          // index
    void add_element_values(index_type elem,
                            const FEDof& dof);  // update global solution vector
                                                // from this element
    void init_values();  // initialize local storage from global solution vector
    void get_zero_values();  // allocate local storage and set values to 0
    void add_values();  // add values to source from local storage
};
```

## The finite element problem assembler: ```feelement.h```
The finite element problem is assembled here.