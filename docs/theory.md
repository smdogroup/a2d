# So, how does A2D solve PDEs?

### Nomenclature

$$
\begin{align*}
&\text{node} &&\text{a control point in an element about which a shape function
is defined} \\
&\text{DOF} &&\text{degree of freedom associated with a node} \\
&n_\text{dof} &&\text{global number of basis nodes of the discretized problem}
\\
&l_{\text{dof}, e} &&\text{local number of basis nodes for $e$-th element}
\\
&x, y, z &&\text{spatial coordinates in the physical coordinate system} \\
&\xi, \eta, \zeta &&\text{spatial coordinates in the local (computational)
coordinate system} \\
&\mathbf{x} && \mathbf{x} = (x, y, z)^T\\
&\boldsymbol{\xi} && \boldsymbol{\xi} = (\xi, \eta, \zeta)^T\\
&\Omega && \text{spatial domain where the PDE is defined} \\
&\Omega_e && \text{spatial domain for $e$-th finite element} \\
&p_\mathbf{x} \qquad &&\text{PDE operator, derivatives are with respect to
$\mathbf{x}$} \\
&f_\mathbf{x} && \text{integration functional derived from $p_\mathbf{x}$,
derivatives are w.r.t. $\mathbf{x}$} \\
&u(\mathbf{x}) &&\text{PDE solution function (infinite-dimensional, defined
everywhere within the domain)} \\
&w(\mathbf{x}) && \text{test function (infinite-dimensional, defined everywhere
within the domain)} \\
&N_e(\mathbf{x}) &&\text{A collection of basis function values at location
$\mathbf{x}$, row vector} \\
&\nabla N_e(\mathbf{x}) &&\text{A collection of basis function derivatives at
location $\mathbf{x}$, 3-by-$l_{\text{dof}, e}$ matrix} \\
&\tilde{u}_e(\mathbf{x}) &&\text{solution functions approximated by bases of
$e$-th finite element} \\
&\tilde{w}_e(\mathbf{x}) &&\text{test functions approximated by bases of $e$-th
finite element} \\
&u_h &&\text{discretized nodal PDE solution (finite-dimensional, defined on
basis nodes)}, \in \mathbb{R}^{n_\text{dof}} \\
&w_h &&\text{discretized nodal test function (finite-dimensional, defined on
basis nodes)}, \in \mathbb{R}^{n_\text{dof}} \\
&I && \text{exact weak form integral} \\
&\tilde{I} &&\text{numerically approximated weak form integral} \\
&J &&\text{Jacobian coordinate transformation, } J =
\dfrac{\partial\mathbf{x}}{\partial\boldsymbol{\xi}}, J^{-1} =
\dfrac{\partial\boldsymbol{\xi}}{\partial \mathbf{x}} \\
&m_q &&\text{quadrature weight for $q$-th quadrature} \\
&P_e &&\text{short-and-wide selection matrix that gets DOFs for $e$-th element}
\\
&u_{h, e}, w_{h, e} &&\text{local nodal solution and test function DOFs for
$e$-th element} \\
\end{align*}
$$

### Subscripts
$$
\begin{align*}
&(~)_h &&\text{by finite element discretization} \\
&(~)_e &&\text{$e$-th element} \\
&(~)_q &&\text{$q$-th quadrature point} \\
\end{align*}
$$

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
p_\mathbf{x}\left(u(\mathbf{x})\right) = 0
$$
and its weak form:
$$
I = \int_{\Omega} f_\mathbf{x}(u, w) d \mathbf{x} = 0
$$
where $p$ is the PDE operator, $u$ is the exact solution function, $w$ is the
test function.
$f$ is some general functional derived from $p$ (e.g. via integration by parts).
A2D solves the weak form numerically by solving the following (potentially
nonlinear) discretized system:
$$
\dfrac{\partial \tilde{I}}{\partial w_h} = 0
\tag{1}
$$
where $\tilde{I}$ is an approximation to $I$ using finite element basis
functions and numerical integration.
$w_h \in \mathbb{R}^{n_\text{dof}}$ is a vector of unknowns where
$n_\text{dof}$ is total number of degrees of freedom nodes of the discretized
finite element system.

### Finite element basis

The equivalent algebraic system to equation (1) is constructed using
finite element method.
Domain $\Omega$ is divided into $n$ mesh elements $\Omega_e$ s.t.
$\bigcup_{e=1}^{n} \Omega_e = \Omega$.
Functionals $u$ and $w$ are approximated by finite-dimensional
vectors $u$ and $w$ using predefined basis functions associated with the mesh.
Within each element $\Omega_e$, the solution and test functions are approximated
by a weighted sum of basis function evaluations, $u(\mathbf{x}) \approx
\tilde{u}_e(\mathbf{x}) = N_e(\mathbf{x}) u_e$ and $w(\mathbf{x}) \approx
\tilde{w}_e(\mathbf{x}) = N_e(\mathbf{x}) w_e$, where row vector $N_e$ is a
collection of basis functions of this element $\Omega_e$, $u_e$ and $w_e$ are
nodal solution and test function DOF for this element.
$u_e$ and $w_e$ can be obtained from global vecotors $u$ and $w$ using
the selection matrix:
$$
u_e = P_e u_{h, e},~~w_e = P_e w_{h, e}.
$$

### Numerical integration

$I$ is evaluated using numerical integration.
In each finite element, a change of variable is performed to transfer the
element from the global physical coordinates to a local computational
coordinates:
$$
I
= \sum_{e=1}^{n_e} \int_{\Omega_e} f_\mathbf{x}(u, w) d \mathbf{x}
= \sum_{e=1}^{n_e} \int_{\Omega_\xi} \det\left(J\left(\mathbf{x}\right)\right)
f_\mathbf{x}(u, w) d \boldsymbol{\xi}
$$
where $n_e$ is number of subdomains, i.e. elements, $J$ is the Jacobian
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
&\approx \sum_{e=1}^{n_e} \int_{\Omega_\xi} \det\left(J(\mathbf{x})\right)
f_\mathbf{x}\left(\tilde{u}_e(\boldsymbol{\xi}), \tilde{w}_e(\boldsymbol{\xi})\right) d
\boldsymbol{\xi} \\
&\approx  \sum_{e=1}^{n_e} \sum_{q=1}^{n_q} m_q
\det\left(J(\boldsymbol{\xi}_q)\right)
f_\mathbf{x}\left(\tilde{u}(\boldsymbol{\xi}_q),
\tilde{w}(\boldsymbol{\xi}_q)\right) \\
&= \tilde{I}(u_h, w_h)
\end{align*}
$$
where $m_q$ is the quadrature weights.

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
\begin{align*}
\nabla_{\mathbf{x}} v &= J^{-T} \nabla_{\boldsymbol{\xi}} v \\
\dfrac{\partial g}{\partial \nabla_{\boldsymbol{\xi}}v} &=
\dfrac{\partial g}{\partial \nabla_{\mathbf{x}}v} J^{-T}\\
\dfrac{\partial^2 g}{\partial \nabla_{\boldsymbol{\xi}}v\partial
\nabla_{\boldsymbol{\xi}}t} &= J^{-1}\dfrac{\partial^2 g}{\partial
\nabla_{\mathbf{x}}v\partial \nabla_{\mathbf{x}}t} \\
\end{align*}
$$
where $v$ and $t$ are some arbitrary functions such as solution or test
function, $g$ is some scalar functional.


Trace identity states that for $b=f(a)$, we have $\bar{b}^T \dot{b} = \bar{a}^T
\dot{a}$, where row vector $\bar{(~)} = \dfrac{\partial s}{\partial (~)}$ for some scalar
function $s$ by definition.
Use the identity above, we have
$$
\begin{align*}
\bar{(\nabla_{\mathbf{x}} v)} \dot{\nabla_{\mathbf{x}} v} =
\bar{(\nabla_{\mathbf{x}} v)} J^{-T} \dot{\nabla_{\boldsymbol{\xi}} v} &=
\bar{(\nabla_{\boldsymbol{\xi}} v)} \dot{\nabla_{\boldsymbol{\xi}} v}  \\
\Rightarrow \bar{(\nabla_{\boldsymbol{\xi}} v)} &= \bar{(\nabla_{\mathbf{x}} v)}
J^{-T} \\
\end{align*}
$$
where we define $\bar{(~)} = \dfrac{\partial g}{\partial (~)}$.
For second order derivatives of $g$ w.r.t. gradients, we first consider it's $j$-column, we define
$\hat{(~)} = \dfrac{\partial }{\partial (~)} \left[\dfrac{\partial g}{\partial
\nabla_{\boldsymbol{\xi}} t}\right]_j$.
Furthermore, we have
$$
\left[\dfrac{\partial g}{\partial \nabla_{\boldsymbol{\xi}} t}\right]_j =
[J^{-1}]_{ji} \left[\dfrac{\partial g}{\partial \nabla_{\mathbf{x}}
t}\right]_i.
$$
Then use the exact same analysis above, we have
$$
\begin{align*}
\hat{(\nabla_{\boldsymbol{\xi}} v)} &= \hat{(\nabla_{\mathbf{x}} v)} J^{-T} \\
\Rightarrow \dfrac{\partial}{\partial \nabla_{\boldsymbol{\xi}}v}
\left[\dfrac{\partial g}{\partial \nabla_{\boldsymbol{\xi}} t}\right]_j &=
\dfrac{\partial}{\partial \nabla_{\mathbf{x}}v} \left[\dfrac{\partial g}{\partial
\nabla_{\boldsymbol{\xi}} t}\right]_j J^{-T}\\
\Rightarrow \left[ \dfrac{\partial}{\partial \nabla_{\boldsymbol{\xi}}v}
\left[\dfrac{\partial g}{\partial \nabla_{\boldsymbol{\xi}} t} \right]_j \right]^T &= [J^{-1}]_{ji} J^{-1} \left[
\dfrac{\partial}{\partial \nabla_{\mathbf{x}}v}
\left[\dfrac{\partial g}{\partial \nabla_{\mathbf{x}} t}\right]_i \right]^T.
\end{align*}
$$
As a result, the transformation for the Hessian matrix is:
$$
\dfrac{\partial^2 g}{\partial \nabla_{\boldsymbol{\xi}}v\partial
\nabla_{\boldsymbol{\xi}}t} = J^{-1}\dfrac{\partial^2 g}{\partial
\nabla_{\mathbf{x}}v\partial \nabla_{\mathbf{x}}t} J^{-T}.
$$



### Function spaces

A2D is able to handle a mixture of different function spaces, including $L_2$,
$H^1$, $H(\text{div})$ and $H(\text{curl})$.
A proper set of function spaces needs to be chosen such that integrand $f$ can
be evaluated (e.g. if $f$ contains $\nabla u$, then $L_2$ alone is not enough to
serve as the function space for $u$ because it does not support first order
derivatives).

### Evaluate the system residual

Given the expression of $I$ above, residuals can be evaluated as follows:
$$
R(u_h) = \dfrac{\partial \tilde{I}(u_h, w_h)}{\partial w} = \sum_{e=1}^{n_e} \sum_{q=1}^{n_q}
m_q \det\left(J(\boldsymbol{\xi}_q)\right) \dfrac{\partial f_{\mathbf{x}}\left(\tilde{u}_e(\boldsymbol{\xi}_q),
\tilde{w}_e(\boldsymbol{\xi}_q)\right)}{\partial w_h}
\tag{2}
$$
note that $f(u, w)$ is linear with respect to $w$, so $R$ is only a function of
$u$.
Next we elaborate on $\dfrac{\partial f_{\mathbf{x}}}{\partial w}$ by assuming
that the only $\tilde{w}_e$-related term in $f_{\mathbf{x}}$ contains $\nabla
\tilde{w}_e$.
Not that despite the example we choose, general rules apply to other function
spaces.
Since we need first order derivatives, we choose $H^1$ space for the bases.
As a result, we can express $\dfrac{\partial f_{\mathbf{x}}}{\partial w_h}$ at each
quadarture point as:
$$
\dfrac{\partial f_{\mathbf{x}}}{\partial w_h}
= \dfrac{\partial f_{\mathbf{x}}}{\partial \nabla_{\boldsymbol{\xi}}\tilde{w}_e}
\dfrac{\partial \nabla_{\boldsymbol{\xi}}\tilde{w}_e}{\partial w_h}
= \dfrac{\partial f_{\mathbf{x}}}{\partial \nabla_{\boldsymbol{\xi}}\tilde{w}_e}
\nabla_{\boldsymbol{\xi}}N_e(\boldsymbol{\xi}) P_e.
\tag{3}
$$

Plug (3) into (2), we have
$$
\begin{align*}
R(u_h)
&= \sum_{e=1}^{n_e} \sum_{q=1}^{n_q} m_q \det\left(J(\boldsymbol{\xi}_q)\right)
\dfrac{\partial f_{\mathbf{x}}\left(\tilde{u}_e(\boldsymbol{\xi}_q),
\tilde{w}_e(\boldsymbol{\xi}_q)\right)}{\partial \nabla_{\boldsymbol{\xi}}\tilde{w}_e}
\nabla_{\boldsymbol{\xi}} N_e(\boldsymbol{\xi}) P_e \\
&= \sum_{e=1}^{n_e} \sum_{q=1}^{n_q} m_q \det\left(J(\boldsymbol{\xi}_q)\right)
J^{-1}(\boldsymbol{\xi}_q) \dfrac{\partial
f_{\mathbf{x}}\left(\tilde{u}_e(\boldsymbol{\xi}_q),
\tilde{w}_e(\boldsymbol{\xi}_q)\right)}{\partial \nabla_{\mathbf{x}}\tilde{w}_e}
\nabla_{\boldsymbol{\xi}} N_e(\boldsymbol{\xi}_q) P_e \\
\end{align*}
$$
Note that the transformations between physical coordinates and computational
coordinates are used.










<!-- A2D implements a set of general operations in ```feelement.h``` that construct
the (linear or nonlinear) algebraic systems derived from the governing PDEs. -->
