# AD core routines

The core AD routines implement automatic differentiation on matrix and vector operations.

When using the AD routines, it is important to keep in mind that the first and second-order derivatives are stored in place. As a result, overwriting intermediate varaibles in an evaluation procedure will produce incorrect derivative results. It is recommended to implement tests for all of your AD code. A2D has a test class that implements tests to ensure that the derivatives are correct.

## Matrix operations

A2D implements the following set of matrix operations

### Matrix multiplication

Given $A \in \mathbb{R}^{n \times m}$ and $B \in \mathbb{R}^{m \times k}$, compute $C = A B$

```c++
MatMatMult(A, B, C);
```

For more general matrix multiplication with transpose, for instance, $C = A^{T} B$

```c++
MatMatMult<MatOp::TRANSPOSE, MatOpt::NORMAL>(A, B, C);
```

Note that the matrices must be the correct size.

### Matrix addition

Given $A, B \in \mathbb{R}^{n \times m}$, compute $C = A + B$

```c++
MatSum(A, B, C);
```

More generally matrix sums can be performed with scalar multiples $\alpha$ and $\beta$ such that $C = \alpha A + \beta B$

```c++
MatSum(alpha, A, beta, B, C);
```

### Symmetrix matrix multiplication

Given $A \in \mathbb{R}^{n \times k}$, compute the symmetric often rank-k matrix $S$ as $S = A A^{T}$

```c++
SymMatRK(A, S);
```

Similarly, $S = A^{T} A$ is

```c++
SymMatRK<MatOpt::TRANSPOSE>(A, S);
```

This operation can also be scaled as $S = \alpha A^{T} A$

```c++
SymMatRK<MatOpt::TRANSPOSE>(alpha, A, S);
```

For this operation, $\alpha$ must be a passive numeric constant.

### Symmetrix matrix addition

Given $A \in \mathbb{R}^{n \times n}$, compute $S = A + A^{T}$

```c++
SymMatSum(A, S);
```

Or, more generally, commpute $S = \alpha(A + A^{T})$ using

```c++
SymMatSum(alpha, A, S);
```

### Matrix inverse

Given $A \in \mathbb{R}^{n \times n}$, compute $B = A^{-1}$ for $n \le 3$

```c++
MatInv(A, B);
```

### Matrix determinant

Given $A \in \mathbb{R}^{n \times n}$, compute $\alpha = \text{det}(A)$

```c++
MatDet(A, alpha);
```

### Matrix trace

Given $A \in \mathbb{R}^{n \times n}$, compute $\alpha = \text{tr}(A)$

```c++
MatTrace(A, alpha);
```

### Green strain

Given $A \in \mathbb{R}^{n \times n}$, compute $E = \frac{1}{2} (A + A^{T} + A^{T} A)$

```c++
MatGreenStrain(A, E);
```

### Isotropic constitutive relationship

Given $\mu$, $\lambda$ and $E$ compute $S = 2 \mu E + \lambda \text{tr}(E) I$

```c++
SymIsotropic(mu, lambda, E, S);
```

### Symmetric multiplication with trace

Given $E, S \in \mathbb{S}^{n}$, compute $\alpha = \text{tr}(E S)$

```c++
SymMatTrace(E, S, alpha);
```

### Matrix scale

- [ ] Complete operation

Given $\alpha$ and $A \in \mathbb{R}^{n \times n}$, compute $B = \alpha A$

```c++
MatScale(alpha, A, B);
```

## Vector operations

A2D implements the following set of vector operations

### Norm

Compute the 2-norm of the vector $x \in \mathbb{R}^{n}$ such that $\alpha = ||x||_{2}$

```c++
VecNorm(x, alpha);
```

### Normalize

Normalize the input vector $x \in \mathbb{R}^{n}$ such that $y = \frac{1}{|| x ||_{2}} x$

```c++
VecNormalize(x, y);
```

### Scale

Scale the input vector $x \in \mathbb{R}^{n}$ such that $y = \alpha x$

```c++
VecScale(alpha, x, y);
```

### Dot product

Given the input vectors $x, y \in \mathbb{R}^{n}$ compute the scalar value $\alpha = x^{T} y$

```c++
VecDot(x, y, alpha);
```

### Vector sum

- [ ] Complete operation

$z = \alpha x + \beta y$

```c++
VecSum(alpha, x, beta, y, z);
```

### Cross-product ($n = 3$ only)

Given the input vectors $x, y \in \mathbb{R}^{3}$ compute the output vector $z = x \times y$

```c++
VecCross(x, y, z);
```

### Outer-product

- [ ] Complete operation

$A = x y^{T}$

```c++
VecOuterProduct(x, y, A);
```

### Symmetric outer-product

- [ ] Complete operation

$S = \alpha x x^{T}$

```c++
VecSymOuterProduct(alpha, x, S);
```

$\alpha$ is a passive numeric constant.

## Example use of A2D routines

The AD routines can be used in the following manner. Consider the computation of the strain energy given the displacement gradient in the computational coordinates $U_{\xi} \in \mathbb{R}^{n \times n}$ and the derivative of the physical coordinates with respect to the computational coordinates $J \in \mathbb{R}^{n \times n}$.

This computation would involve the following steps:

* Compute $J^{-1}$
* Compute $U_{x} = U_{\xi} J^{-1}$
* Compute $F = I + U_{x}$
* Compute $E = \frac{1}{2} F^{T} F$
* Compute $S = 2 \mu E + \lambda \text{tr}(E) I$
* Compute $U = \text{tr}(E S)$

Note that this computes two times the strain energy.

In A2D, this could be implemented with the following sequence of steps:

```c++
Mat<T, N, N> Uxi, J; // Input
T output;            // Output

Mat<T, N, N> Jinv, Ux, F, Id;
SymMat<T, N> E, S;

MatInv(J, Jinv);
MatMatMult(Uxi, Jinv, Ux);
MatSum(Ux, Id, F);
SymMatRK<MatOp::TRANSPOSE>(T(0.5), F, E);
SymIsotropic(T(0.35), T(0.51), E, S);
SymMatTrace(E, S, output);
```

Here $T$ and $N$ are template parameters. The scalar, matrix and symmetric matrix types are all A2D types designed to be used directly for computations without AD information.

To implement AD for first derivatives, A2D uses AD types within the code. The A2D operations support a mixture of AD types and regular types that represent passive variables.

The first order AD types use a container where the data for the value and the first derivative must be passed into the constructor for matrix and vector types. For instance the constructor for `Uxi` takes the value and derivative `Uxi0` and `Uxib`.

The stack of operations is created using a call to `MakeStack`. Once created, the stack can be used to execute the sequence of calls needed to perform reverse mode AD.

The seed for the reverse mode AD is set by the statement `output.bvalue = 1.0;`, and the call to `stack.reverse();` performs the reverse mode AD. After this call, the values `Uxib` and `Jb` contain the desired derivatives.

```c++
// Passive matrix (constant)
Mat<T, N, N> Id;

// Input and derivative output
ADMat<Mat<T, N, N>> Uxi(Uxi0, Uxib), J(J0, Jb);

// Output
ADScalar<T> output;

// Intermediate values
ADMat<Mat<T, N, N>> Jinv(Jinv0, Jinvb), Ux(Ux0, Uxb), F(F0, Fb);
ADMat<SymMat<T, N>> E(E0, Eb), S(S0, Sb);

auto stack = MakeStack(
    MatInv(J, Jinv),
    MatMatMult(Uxi, Jinv, Ux),
    MatSum(Ux, Id, F),
    SymMatRK<MatOp::TRANSPOSE>(T(0.5), F, E),
    SymIsotropic(T(0.35), T(0.51), E, S),
    SymMatTrace(E, S, output));

// Set the seed value
output.bvalue = 1.0;

// Reverse mode AD through the stack
stack.reverse();

// Jb and Uxib contain the derivatives
```

For second derivatives, A2D uses Hessian-vector products that require a combination of forward and reverse mode.

```c++
// Passive matrix (constant)
Mat<T, N, N> Id;

// Input
A2DMat<Mat<T, N, N>> Uxi, J;

// Outputs
A2DScalar<T> output;

// Intermediate data
A2DMat<Mat<T, N, N>> Jinv, Ux, F;
A2DMat<SymMat<T, N>> E, S;

auto stack = MakeStack(
    MatInv(J, Jinv),
    MatMatMult(Uxi, Jinv, Ux),
    MatSum(Ux, Id, F),
    SymMatRK<MatOp::TRANSPOSE>(T(0.5), F, E),
    SymIsotropic(T(0.35), T(0.51), E, S),
    SymMatTrace(E, S, output));

// Set the seed value and the second derivative value
output.bvalue = 1.0;
stack.reverse();

// Set values for the direction p
// Set Uxi.pvalue();
// Set J.pvalue();
output.hvalue = 0.0;

// Perform the forward and reverse passes
stack.hforward();
stack.hreverse();

// Second derivatives are now available in Uxi.hvalue(); and J.hvalue();
```

## Testing AD expressions