# Github markdown equation glitches

Github started to support equation typesetting around 2022 [(see this
post)](https://docs.github.com/en/get-started/writing-on-github/working-with-advanced-formatting/writing-mathematical-expressions)
with a handful of bugs (some of them are features perhaps?).
Below is a short list of pitfalls I stepped into.

| Gotcha! | Reason | Solution |
| ------- | ------ | -------- |
|```${a}_{1} \sum_{x}$```| ```_..._``` is sometimes parsed as _italic_ | add space before and after ```_```|
|```${a}~{1} \sum~{x}$```|unclear|add space before and after ```~```|
|$\text{Text in equation is very messy}$|unclear | don't use ```$\text{...}$``` for long text strings|
|```$...```<br />```...$```| inline equation must stay in single line | put equation in single line: ```$... ...$```|


# So, how does A2D solve PDEs?

### Nomenclature

| Symbol | Definition |
| ------ | ---------- |
| node   |a control point in an element about which a shape function is defined|
| DOF    |degree of freedom associated with a node|
|$n _ \text{dof}$ |global number of basis nodes of the discretized problem |
|$l _ {\text{dof}, e} $|local number of basis nodes for $e$-th element|
|$x, y, z $|spatial coordinates in the physical coordinate system|
|$\xi, \eta, \zeta $|spatial coordinates in the local (computational) coordinate system|
|$\mathbf{x} $| $\mathbf{x} = (x, y, z)^T$|
|$\boldsymbol{\xi} $| $\boldsymbol{\xi} = (\xi, \eta, \zeta)^T$|
|$\Omega $| spatial domain where the PDE is defined|
|$\Omega _ e $| spatial domain for $e$-th finite element|
|$p _ \mathbf{x} \qquad $|PDE operator, derivatives are with respect to $\mathbf{x}$|
|$f _ \mathbf{x} $| integration functional derived from $p _ \mathbf{x}$, derivatives are w.r.t. $\mathbf{x}$|
|$u(\mathbf{x}) $|PDE solution function (infinite-dimensional, defined everywhere within the domain)|
|$w(\mathbf{x}) $| test function (infinite-dimensional, defined everywhere within the domain)|
|$N _ e(\mathbf{x}) $|A collection of basis function values at location $\mathbf{x}$, row vector|
|$\nabla N _ e(\mathbf{x}) $|A collection of basis function derivatives at location $\mathbf{x}$, 3-by-$l _ {\text{dof}, e}$ matrix|
|$\tilde{u} _ e(\mathbf{x}) $|solution functions approximated by bases of $e$-th finite element|
|$\tilde{w} _ e(\mathbf{x}) $|test functions approximated by bases of $e$-th finite element|
|$u _ h $|discretized nodal PDE solution (finite-dimensional, defined on basis nodes), $\in \mathbb{R}^{n _ \text{dof}}$|
|$w _ h $|discretized nodal test function (finite-dimensional, defined on basis nodes), $\in \mathbb{R}^{n _ \text{dof}}$|
|$I $| exact weak form integral|
|$\tilde{I} $|numerically approximated weak form integral|
|$J$| Jacobian coordinate transformation, $J = \dfrac{\partial\mathbf{x}}{\partial\boldsymbol{\xi}}, J^{-1} = \dfrac{\partial\boldsymbol{\xi}}{\partial \mathbf{x}}$|
|$m _ q$|quadrature weight for $q$-th quadrature|
|$P _ e $|short-and-wide selection matrix that gets DOFs for $e$-th element|
|$u _ {h, e}, w _ {h, e} $|local nodal solution and test function DOFs for $e$-th element|


| Subscripts | Definition |
| ---------- | ---------- |
|$( ~ ) _ h $|by finite element discretization|
|$( ~ ) _ e $|$e$-th element|
|$( ~ ) _ q $|$q$-th quadrature point|

### Convention for matrix calculus

- Vectors are column vectors (if not stated otherwise)
- Gradient of a scalar (i.e. $\nabla x$) is column vector
- Gradient of a row vector (i.e. $\nabla N$) is a matrix where each column is
a gradient of one entry of that row vector
- [Numerator-layout
notation](https://en.wikipedia.org/wiki/Matrix_calculus#Numerator-layout_notation)
is employed, i.e. $\dfrac{\partial x}{\partial y} \in \mathbb{R}^{m\times n}$
for $x \in \mathbb{R}^m$ and $y \in \mathbb{R}^n$


### PDE and weak form

Consider the following PDE:

$$
p _ \mathbf{x}\left(u(\mathbf{x})\right) = 0
$$

and its weak form:

$$
I = \int _ {\Omega} f _ \mathbf{x}(u, w) d \mathbf{x} = 0
$$

where $p$ is the PDE operator, $u$ is the exact solution function, $w$ is the
test function.
$f$ is some general functional derived from $p$ (e.g. via integration by parts).
A2D solves the weak form numerically by solving the following (potentially
nonlinear) discretized system:

$$
\dfrac{\partial \tilde{I}}{\partial w _ h} = 0
\tag{1}
$$

where $\tilde{I}$ is an approximation to $I$ using finite element basis
functions and numerical integration.
$w _ h \in \mathbb{R}^{n _ \text{dof}}$ is a vector of unknowns where
$n _ \text{dof}$ is total number of degrees of freedom nodes of the discretized
finite element system.

### Finite element basis

The equivalent algebraic system to equation (1) is constructed using
finite element method.
Domain $\Omega$ is divided into $n_e$ mesh elements $\Omega _ e$ s.t.
$\cup _ {e=1}^{n_e} \Omega _ e = \Omega$.
Functionals $u$ and $w$ are approximated by finite-dimensional
vectors $u$ and $w$ using predefined basis functions associated with the mesh.
Within each element $\Omega _ e$, the solution and test functions are approximated
by a weighted sum of basis function evaluations

$$
\begin{align*}
u(\mathbf{x}) &\approx \tilde{u} _ e(\mathbf{x}) = N _ e(\mathbf{x}) u _ e \\
w(\mathbf{x}) &\approx \tilde{w} _ e(\mathbf{x}) = N _ e(\mathbf{x}) w _ e \\
\end{align*}
$$

where row vector $N _ e$ is a
collection of basis functions of this element $\Omega _ e$, $u _ e$ and $w _ e$ are
nodal solution and test function DOF for this element.
$u _ e$ and $w _ e$ can be obtained from global vecotors $u$ and $w$ using
the selection matrix:

$$
u _ e = P _ e u _ {h, e}, ~  ~ w _ e = P _ e w _ {h, e}.
$$

### Numerical integration

$I$ is evaluated using numerical integration.
In each finite element, a change of variable is performed to transfer the
element from the global physical coordinates to a local computational
coordinates:

$$
I
= \sum _ {e=1}^{n _ e} \int _ {\Omega _ e} f _ \mathbf{x}(u, w) d \mathbf{x}
= \sum _ {e=1}^{n _ e} \int _ {\Omega _ \xi} \det\left(J\left(\mathbf{x}\right)\right)
f _ \mathbf{x}(u, w) d \boldsymbol{\xi}
$$

where $n _ e$ is number of subdomains, i.e. elements, $J$ is the Jacobian
transformation matrix, for 3-dimensional problem, $J$ is

$$
J = \dfrac{\partial \mathbf{x}}{\partial \boldsymbol{\xi}}
= \begin{bmatrix}
\dfrac{\partial x}{\partial \xi} & \dfrac{\partial x}{\partial \eta} &
\dfrac{\partial x}{\partial \zeta} \\
\dfrac{\partial y}{\partial \xi} & \dfrac{\partial y}{\partial \eta} &
\dfrac{\partial y}{\partial \zeta} \\
\dfrac{\partial z}{\partial \xi} & \dfrac{\partial z}{\partial \eta} &
\dfrac{\partial z}{\partial \zeta} \\
\end{bmatrix}
$$

Next, the exact integral is approximated using basis functions and numerically
evaluated using quadrature rule:

$$
\begin{align*}
I
&\approx \sum _ {e=1}^{n _ e} \int _ {\Omega _ \xi} \det\left(J(\mathbf{x})\right)
f _ \mathbf{x}\left(\tilde{u} _ e(\boldsymbol{\xi}), \tilde{w} _ e(\boldsymbol{\xi})\right) d
\boldsymbol{\xi} \\
&\approx  \sum _ {e=1}^{n _ e} \sum _ {q=1}^{n _ q} m _ q
\det\left(J(\boldsymbol{\xi} _ q)\right)
f _ \mathbf{x}\left(\tilde{u}(\boldsymbol{\xi} _ q),
\tilde{w}(\boldsymbol{\xi} _ q)\right) \\
&= \tilde{I}(u _ h, w _ h)
\end{align*}
$$

where $m _ q$ is the quadrature weights.

### Reference frame transformation

Element-wise computations are done with the help of a local (computational)
coordinate system.
For example, for a hexahedral element, its computational coordinate system is
often chosen as $\boldsymbol{\xi} \in [-1, 1]^3$.
Such auxiliary frame helps to simplify the computation as evaluations of the
function and, more importantly, derivatives can be separated.
Once those local computations are done, they need to be transferred to the
global (physical) systems by the following rules:

$$
\begin{align}
\nabla _ {\mathbf{x}} v &= J^{-T} \nabla _ {\boldsymbol{\xi}} v  \tag{1a} \\
\dfrac{\partial g}{\partial \nabla _ {\boldsymbol{\xi}}v} &=
\dfrac{\partial g}{\partial \nabla _ {\mathbf{x}}v} J^{-T} \tag{1b} \\
\dfrac{\partial^2 g}{\partial \nabla _ {\boldsymbol{\xi}}v\partial
\nabla _ {\boldsymbol{\xi}}t} &= J^{-1}\dfrac{\partial^2 g}{\partial
\nabla _ {\mathbf{x}}v\partial \nabla _ {\mathbf{x}}t} J^{-T} \tag{1c} \\
\end{align}
$$

where $v$ and $t$ are some arbitrary functions such as solution or test
function, $g$ is some scalar functional.
(1a) is straightforward by chain rule, we derive (1b) and (1c) below.
Trace identity states that for $b=f(a)$, we have $\bar{b}^T \dot{b} = \bar{a}^T \dot{a}$,
where row vector $\bar{( ~ )} = \dfrac{\partial s}{\partial ( ~ )}$ for some scalar
function $s$ by definition.
Use the trace identity, we have

$$
\begin{align*}
\bar{(\nabla _ {\mathbf{x}} v)} \dot{\nabla _ {\mathbf{x}} v} =
\bar{(\nabla _ {\mathbf{x}} v)} J^{-T} \dot{\nabla _ {\boldsymbol{\xi}} v} &=
\bar{(\nabla _ {\boldsymbol{\xi}} v)} \dot{\nabla _ {\boldsymbol{\xi}} v}  \\
\Rightarrow \bar{(\nabla _ {\boldsymbol{\xi}} v)} &= \bar{(\nabla _ {\mathbf{x}} v)}
J^{-T} \\
\end{align*}
$$

where we define $\bar{( ~ )} = \dfrac{\partial g}{\partial ( ~ )}$.
For second order derivatives of $g$ w.r.t. gradients, we first consider it's
$j$-column, we define

$$
\hat{( ~ )} = \dfrac{\partial }{\partial ( ~ )} \left[\dfrac{\partial g}{\partial
\nabla _ {\boldsymbol{\xi}} t}\right] _ j.
$$

Furthermore, (1b) gives us

$$
\left[\dfrac{\partial g}{\partial \nabla _ {\boldsymbol{\xi}} t}\right] _ j =
[J^{-1}] _ {ji} \left[\dfrac{\partial g}{\partial \nabla _ {\mathbf{x}}
t}\right] _ i.
$$

Then use the exact same analysis above, we have

$$
\begin{align*}
\hat{(\nabla _ {\boldsymbol{\xi}} v)} &= \hat{(\nabla _ {\mathbf{x}} v)} J^{-T} \\
\Rightarrow \dfrac{\partial}{\partial \nabla _ {\boldsymbol{\xi}}v}
\left[\dfrac{\partial g}{\partial \nabla _ {\boldsymbol{\xi}} t}\right] _ j &=
\dfrac{\partial}{\partial \nabla _ {\mathbf{x}}v} \left[\dfrac{\partial g}{\partial
\nabla _ {\boldsymbol{\xi}} t}\right] _ j J^{-T}\\
\Rightarrow \left[ \dfrac{\partial}{\partial \nabla _ {\boldsymbol{\xi}}v}
\left[\dfrac{\partial g}{\partial \nabla _ {\boldsymbol{\xi}} t} \right] _ j \right]^T &= [J^{-1}] _ {ji} J^{-1} \left[
\dfrac{\partial}{\partial \nabla _ {\mathbf{x}}v}
\left[\dfrac{\partial g}{\partial \nabla _ {\mathbf{x}} t}\right] _ i \right]^T.
\end{align*}
$$

As a result, the transformation for the Hessian matrix is:

$$
\dfrac{\partial^2 g}{\partial \nabla _ {\boldsymbol{\xi}}v\partial
\nabla _ {\boldsymbol{\xi}}t} = J^{-1}\dfrac{\partial^2 g}{\partial
\nabla _ {\mathbf{x}}v\partial \nabla _ {\mathbf{x}}t} J^{-T}.
$$

### Function spaces

A2D is able to handle a mixture of different function spaces, including $L _ 2$,
$H^1$, $H(\text{div})$ and $H(\text{curl})$.
A proper set of function spaces needs to be chosen such that integrand $f$ can
be evaluated (e.g. if $f$ contains $\nabla u$, then $L _ 2$ alone is not enough to
serve as the function space for $u$ because it does not support first order
derivatives).

### Evaluate the system residual

Given the expression of $I$ above, residuals can be evaluated as follows:

$$
R(u _ h) = \dfrac{\partial \tilde{I}(u _ h, w _ h)}{\partial w} = \sum _ {e=1}^{n _ e} \sum _ {q=1}^{n _ q}
m _ q \det\left(J(\boldsymbol{\xi} _ q)\right) \dfrac{\partial f _ {\mathbf{x}}\left(\tilde{u} _ e(\boldsymbol{\xi} _ q),
\tilde{w} _ e(\boldsymbol{\xi} _ q)\right)}{\partial w _ h}
\tag{2}
$$

note that $f(u, w)$ is linear with respect to $w$, so $R$ is only a function of
$u$.
Next we elaborate on $\dfrac{\partial f _ {\mathbf{x}}}{\partial w}$ by assuming
that the only $\tilde{w} _ e$-related term in $f _ {\mathbf{x}}$ contains
$\nabla \tilde{w} _ e$.
Not that despite the example we choose, general rules apply to other function
spaces.
Since we need first order derivatives, we choose $H^1$ space for the bases.
As a result, we can express $\dfrac{\partial f _ {\mathbf{x}}}{\partial w _ h}$ at each
quadarture point as:

$$
\dfrac{\partial f _ {\mathbf{x}}}{\partial w _ h}
= \dfrac{\partial f _ {\mathbf{x}}}{\partial \nabla _ {\boldsymbol{\xi}}\tilde{w} _ e}
\dfrac{\partial \nabla _ {\boldsymbol{\xi}}\tilde{w} _ e}{\partial w _ h}
= \dfrac{\partial f _ {\mathbf{x}}}{\partial \nabla _ {\boldsymbol{\xi}}\tilde{w} _ e}
\nabla _ {\boldsymbol{\xi}}N _ e(\boldsymbol{\xi}) P _ e.
\tag{3}
$$

Plug (3) into (2), we have

$$
\begin{align*}
R(u _ h)
&= \sum _ {e=1}^{n _ e} \sum _ {q=1}^{n _ q} m _ q \det\left(J(\boldsymbol{\xi} _ q)\right)
\dfrac{\partial f _ {\mathbf{x}}\left(\tilde{u} _ e(\boldsymbol{\xi} _ q),
\tilde{w} _ e(\boldsymbol{\xi} _ q)\right)}{\partial \nabla _ {\boldsymbol{\xi}}\tilde{w} _ e}
\nabla _ {\boldsymbol{\xi}} N _ e(\boldsymbol{\xi}) P _ e \\
&= \sum _ {e=1}^{n _ e} \sum _ {q=1}^{n _ q} m _ q \det\left(J(\boldsymbol{\xi} _ q)\right)
J^{-1}(\boldsymbol{\xi} _ q) \dfrac{\partial
f _ {\mathbf{x}}\left(\tilde{u} _ e(\boldsymbol{\xi} _ q),
\tilde{w} _ e(\boldsymbol{\xi} _ q)\right)}{\partial \nabla _ {\mathbf{x}}\tilde{w} _ e}
\nabla _ {\boldsymbol{\xi}} N _ e(\boldsymbol{\xi} _ q) P _ e \\
\end{align*}
$$

where the transformations between physical coordinates and computational
coordinates are used.