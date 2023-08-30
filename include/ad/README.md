# AD core routines

The core AD routines implement automatic differentiation on matrix and vector operations.

When using the AD routines, it is important to keep in mind that the first and second-order derivatives are stored in place. As a result, overwriting intermediate varaibles in a solution procedure will produce incorrect derivative results. It is recommended to implement tests for all of your AD code. A2D has a test class that implements tests to ensure that the derivatives are correct.

## Matrix operations

The matrix operations include the following set of operations

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

Given $\mu$, $\lambda$ and $E$ compute $S = 2 \mu E + \lambda I \text{tr}(E)$

```c++
SymIsotropic(mu, lambda, E, S);
```