# Class template skeletons

## Function spaces: ```fespace.h```
Two types of class templates are defined:
- ```XXXSpace```
- ```FESpace```

 ### ```XXXSpace```
```XXXSpace``` implements specific function function space instances, such as
L2 space, H1 space, H(div) space, etc.
Such class templates comply with the following skeleton:
 ```c++
 class XXXSpace {
    XXXSpace();
    static const index_type ncomp;  // number of components (e.g. number of
                                    // solution variables and derivatives
                                    // combined)
    static const index_type dim;  // spatial dimension (e.g. 2D or 3D)
    void set_seed(const index_type comp);  // set input seed (1.0) comp-th
                                           // component (could be solution or
                                           // derivative)
    sol_type& get_value(const index_type comp);  // get the value of comp-the
                                                 // component
    sol_type& get_value();  // get solution variable (could be scalar or vector
                            // type)
    sol_type& get_value() const;
    void transform(const T& detJ, const mat_type& J, const mat_type& Jinv,
                   XXXSpace& s);  // transform from this (reference) element to
                                  // physical element
    void rtransform(const T& detJ, const mat_type& J, const mat_type& Jinv,
                   XXXSpace& sref);  // transform from this (physical) element
                                     // to reference element

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
- ```XXXBasis```
- ```FEBasis```

### ```XXXBasis```
```XXXBasis``` implements specific element basis, such as Lagrange element,
Raviart-Thomas element, etc.
Such class templates comply with the following skeleton:
```c++
class XXXBasis {
    static const index_type ndof;  // number of degrees of freedom of such
                                   // element
    static const index_type ncomp;  // number of components (e.g. number of
                                    // solution variables and derivatives
                                    // combined)
    static void interp(index_type n, const sol_type sol,
                       function_space_type out);  // at n-th quadrature point,
                                                  // given solution to the local
                                                  // dof, evaluate components of
                                                  // function space object out
    static void add(index_type n, const function_space_type in,
                    const sol_type res);  // at n-th quadrature point, given
                                          // function space object in, evaluate
                                          // local degrees of freedom and add to
                                          // res (of a local solution type)
    static const index_type stride;
    static void basis(index_type n, basis_type N);  // evaluate basis functions
                                                    // N at n-th quadrature point
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
class XXXQuadrature {
    static index_type get_num_points();  // return number of quadrature points
    static void get_point(const index_type n,
                          double pt[]);  // get n-th quadrature point in reference
                                         // coordinates
    static double get_weight(const index_type n);  // get weighting factor of n-th
                                                 // quadrature point
};
```


## PDE implementations, for example: ```poisson.h```
A2D solves PDEs. Such problem instances need to be implemented to comply with
the following skeleton:
```c++
class XXXPDE {
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
Such interfaces are called ```ElementVector_XXX```, and must comply with the
following skeleton:
```c++
class ElementVector_XXX {

    ElementVector_XXX(...);  // takes in global solution vector and other needed
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
    // ElementVector_XXX implementation
    void get_element_values(index_type elem,
                            FEDof& dof);  // populate dof object, elem is element
                                          // index
    void add_element_values(index_type elem,
                            const FEDof& dof);  // update global solution vector
                                                // from this element
    void init_values();  // initialize local storage from global solution vector
    void init_zero_values();  // allocate local storage and set values to 0
    void add_values();  // add values to source from local storage
};
```

## The finite element problem assembler: ```feelement.h```
The finite element problem is assembled here.