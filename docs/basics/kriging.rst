Kriging
=======

Intro
-----

Kriging is **"a collection of gneralized linear regression techniques for "**
**"minimizing an estimation variance defined from a priori model for a convariance".**

Consider the estimate of an unsampled value :math:`z(\mathbf{u})` from
neighboring data values :math:`z({\mathbf{u}}_{\alpha}),\alpha=1,\dots,n`.
The RF model :math:`Z(\mathbf{u})` is stationary with mean :math:`m` and
covariance :math:`C(\mathbf{h})`. In its simplest form, also known as
**Simple Kriging (SK)**, the algorithm considers the following linear estimator:

.. math::
    {Z}_{SK}^{*}(\mathbf{u})=\sum_{\alpha=1}^{n}{\lambda}_{\alpha}(\mathbf{u})
    Z({\mathbf{u}}_{\alpha})+\left(1-\sum_{\alpha=1}^{n}{\lambda}_{\alpha}
    (\mathbf{u})\right)m

The weights :math:`{\lambda}_{\alpha}(\mathbf{u})` are determined to minnimize
the error variance, also called the "estimation variance." That minimization
results in a set of normal equations:

.. math::
    \sum_{\beta=1}^{n}{\lambda}_{\beta}(\mathbf{u})C({\mathbf{u}}_{\beta}-{\mathbf{u}}_{\alpha})=C(\mathbf{u}-{\mathbf{u}}_{\alpha})\\
    \forall{\alpha}=1,\dots,n

The corresponding minimized estimation variance, or kriging variance, is:

.. math::
    {\sigma}_{SK}^{2}(\mathbf{u})=C(0)-\sum_{\alpha=1}^{n}{\lambda}_{\alpha}(\mathbf{u})C(\mathbf{u}-{\mathbf{u}}_{\alpha})\geq 0

**Ordinary Kriging (OK)** is the most commonly used variant of the previous
simple kriging algorithm, whereby the sum of the weights
:math:`\sum_{\alpha=1}^{n}{\lambda}_{\alpha}(\mathbf{u})` is considered to equal 1.
This allows building an estimator :math:`Z_{OK}^{*}(\mathbf{u})` that does not
require prior knowledge of the stationary mean m, yet remains unbiased in the
sense that :math:`E\{{Z}_{OK}^{*}(\mathbf{u})\}=E\{Z(\mathbf{u})\}`.

**None-linear kriging** is but linear kriging performed on some non-linear
transform of the z-data, e.g., the **log-transform**
:math:`\mathrm{ln}z` privided that :math:`z>0`, or the **indicator transform**
as defined in relation:

.. math::
    I(\mathbf{u};z)=\begin{cases}
    1,\quad Z(\mathbf{u})\leq z\\
    0,\quad \text{otherwise}
    \end{cases}

Traditionally, kriging (**SK** or **OK**) has been performed to provide
a "best" linear unbiased estimate (**BLUE**) for unsampled values
:math:`z(\mathbf{u})`, with the kriging variance being used to define Gaussian-type
confidence intervals, e.g.,

.. math::
    \text{Prob}\{Z(\mathbf{u})\in [{z}_{SK}^{*}(\mathbf{u})\pm 2{\sigma}_{SK}(\mathbf{u})]\cong 0.95\}

**Unfortunately, kriging variances of the type, being independent of the data**
**values, only provides a comparison of alternative geometric data configurations.**
**Kriging variances are usually not measures of local estimation accuracy.**

The kriging algorithm has two characteristic properties that allow its
use in determining **posterior ccdfs**. These two characteristic properties
are the basis for, respectively, the **multi-Gaussian (MG) approach** and
the **indicator kriging (IK) approach** to determination of ccdfs:

1. The Multi-Gaussian Approach: If the RF model :math:`Z(\mathbf{u})` is
multivariate Gaussian, then the simple kriging estimate and variance identify
the mean and variance of the posterior ccdf. In addition, since that ccdf is
Gaussian, it is fully determined by these two parameters. This remarkable
result is at the basis of multi-Gaussian (MG) kriging and simulation.
The MG approach is said to be **parametric** in the sense that it determines
the ccdfs through their parameters (mean and variance). The MG algorithm is
remarkably fast and trouble-free; its limitation is the reliance on the very
specific and sometimes inappropriate properties fo the Gaussian RF model.

2. The Indicator Kriging Approach: If the value to be estimated is the expected
value (mean) of a distribution, then least-squares (LS) regression
(i.e., kriging) is a priori the prefered algorithm. The reason is that
the LS estimator fo the variable :math:`Z(\mathbf{u})` is also the LS estimator of
its conditional expectation :math:`E\{Z(\mathbf{u})\mid(n)\}`, that is, of the
expected value of the ccdf. Instead of the variable :math:`Z(\mathbf{u})`,
consider its binary indicator transform :math:`I(\mathbf{u};z)`.
Kriging of the indicator RV :math:`I(\mathbf{u};z)` provides an estimate that
is also the best LS estimate fo the conditional expectation of :math:`I(\mathbf{u};z)`.
Now, the conditional expecation of :math:`I(\mathbf{u};z)`
is equal to the ccdf of :math:`Z(\mathbf{u})`, indeed:

.. math::
    \begin{array}{l}
    E\{I(\mathbf{u};z)\mid (n)\}&=
    \begin{array}{l}
    1\cdot \text{Prob}\{I(\mathbf{u};z)=1\mid (n)\}\\
    +0\cdot \text{Prob}\{I(\mathbf{u};z)=0\mid (n)\}
    \end{array}\\
    &=1\cdot \text{Prob}\{Z(\mathbf{u})\leq z\mid (n)\}\\\
    &\equiv F(\mathbf{u};z\mid (n))
    \end{array}

Thus the kriging algorithm applied to indicator data provides LS estimates
of the ccdf. Note that indicator kriging (IK) is not aimed at estimating the
unsampled value :math:`z(\mathbf{u})` or its indicator transform
:math:`I(\mathbf{u};z)` but at providing a ccdf model of uncertainty about
:math:`z(\mathbf{u})`. The IK algorithm is said to be **non-parametric**
in the sense that is does not approach the ccdf through its parameters
(mean and variance); rather, the ccdf values for various threshold values
:math:`z` are estimated directly.

Types of Kriging
----------------

+-----------------------------------+---------+--------------------------+--------------+
| Kriging Form                      | Mean    | Drift Model              | Prerequisite |
+===================================+=========+==========================+==============+
| Simple Kriging (SK)               | Known   | None                     | Covariance   |
+-----------------------------------+---------+--------------------------+--------------+
| Ordinary Kriging (OK)             | Unknown | Constant                 | Variogram    |
+-----------------------------------+---------+--------------------------+--------------+
| Universal Kriging (UK)            | Unknown | Functions of coordinates | Variogram    |
+-----------------------------------+---------+--------------------------+--------------+
| Kriging with external drift (KED) | Unknown | External variable        | Variogram    |
+-----------------------------------+---------+--------------------------+--------------+


Simple Kriging
--------------

In its simplist form, also known as simple kriging (SK), the algorithm considers
the following linear estimator:

.. math::
    Z_{SK}^{*}(\mathbf{u}) = \sum_{\alpha=1}^{n} \lambda_{\alpha}(\mathbf{u}) Z(\mathbf{u_{\alpha}}) + \left(1-\sum_{\alpha=1}^{n}\lambda_{\alpha}(\mathbf{u})\right) m

The weights :math:`\lambda_{\alpha}` are determined to minimize the error
variance, also called the "estimation vairiance." That minimization result
in a set of normal equations known as *Simple Kriging System*:

.. math::
    \sum_{\beta=1}^{n} \lambda_{\beta}(\mathbf{u}) C(\mathbf{u_{\beta}}-\mathbf{u_{\alpha}})=C(\mathbf{u}-\mathbf{{u}_{\alpha}}),\\\forall \alpha=1, ... , n

In matrix notation, we have

.. math::
    \boldsymbol{\Sigma}\boldsymbol{\lambda}=\boldsymbol{\sigma_{0}}

where :math:`\boldsymbol{\Sigma}=[{\sigma}_{\alpha\beta}]` is the :math:`N\times N`
matrix of data-to-data covariances, :math:`\boldsymbol{\sigma_{0}}=[{\sigma}_{\alpha0}]`
is the N-vector of covariances between the data and the target, and :math:`\boldsymbol{\lambda}=[\lambda_\alpha]`
is the N-vector of solutions.

The corresponding minimized estimation variance, or kirging variance, is:

.. math::
    \sigma_{SK}^{2}(\mathbf{u}) = C(0) - \sum_{\lambda=1}^{n}\lambda_{\alpha}(\mathbf{u}) C(\mathbf{u}-\mathbf{u_{\alpha}}) \geq 0


Ordinary Kriging (OK)
---------------------

Ordinary Kriging (OK) filters the mean from the SK estimator by requiring that
the kriging weights sum to one. This results in the following ordinary kriging
estimator:

.. math::
    {Z}_{OK}^{*}(\mathbf{u})=\sum_{\alpha=1}^{n}{{\lambda}_{\alpha}^{(OK)}(\mathbf{u})Z({\mathbf{u}}_{\alpha})}

and the sationary OK system:

.. math::
    \begin{cases}
    \sum_{\beta=1}^{n}{{\lambda}_{\beta}^{(OK)}(\mathbf{u}) C({\mathbf{u}}_{\beta}-{\mathbf{u}}_{\alpha})}+\mu(\mathbf{u})=C(\mathbf{u}-{\mathbf{u}}_{\alpha}),\quad \alpha=1,\dots,n \\
    \sum_{\beta=1}^{n}{{\lambda}_{\beta}^{(OK)}(\mathbf{u})}=1\\
    \end{cases}

In matrix notation, the above linear equations correspond to:

.. math::
    \begin{bmatrix}{C}_{11} & {C}_{12} & \cdots & {C}_{1N} & 1 \\{C}_{21} & {C}_{22} & \cdots & {C}_{2N} & 1 \\ \vdots & \vdots & \cdots & \vdots & 1 \\ {C}_{N1} & {C}_{N2} & \cdots & {C}_{NN} & 1 \\ 1 & 1 & 1 & 1 & 0\end{bmatrix} \times \begin{bmatrix}{\lambda}_{1}\\{\lambda}_{2}\\ \vdots \\{\lambda}_{N}\\ \mu \end{bmatrix}
    = \begin{bmatrix}{C}_{10}\\{C}_{20}\\ \vdots \\{C}_{N0}\\ 1 \end{bmatrix}


The kriging variance is obtained by multiplying the first N equations of the
kriging system by :math:`\lambda_\alpha`, summing over :math:`\alpha`, and
then using the last equations. The result is the OK variance:

.. math::
    {\sigma}_{OK}^{2}=E{({Z}^{*}-{Z}_{0})}^{2}={\sigma}_{00}-\sum\limits_{\alpha}{{\lambda}_{\alpha}{\sigma}_{\alpha0}}-\mu

The linear system has a unique solution if and only if the covarance matrix
:math:`\boldsymbol{\Sigma}[{\sigma}_{\alpha\beta}]` is strictly positive
definite, which is the case if we use strictly positive definite covariance
function model and if all data are distinct.

Universal Kriging (UK) or Kriging with a Trend Model (KT)
---------------------------------------------------------

The general model, which Matheron(1969) named the *universal kriging* model
for reasons explained below, assumes that the mean function can be represented
as a reponse surface function

.. math::
    m(x)=\sum\limits_{\mathscr{l}=0}^{L}{{a}_{\mathscr{l}}{f}^{\mathscr{l}}(x)}

where the :math:`{f}^{\mathscr{l}}(x)` are kown basis functions and :math:`{a}_{\mathscr{l}}`
are fixed but unknown coefficients. Usually the first basis function
(case :math:`\mathscr{l}=0`) is the constant function identically equal to 1,
which guarantees that the constant-mean case is included in the model.
The other functions are typically monomials of low degree in the cooridinates
of x (in practice, the degree does not exceed two). In the case of monomials,
the superscript :math:`\mathscr{l}`, which is an index, has the meaning of a
power (in 1D, :math:`{f}^{\mathscr{l}}(x)={x}^{\mathscr{l}}`). Note that the
above function may be regarded as a local approximation to :math:`m(x)`; that
is, the coefficients :math:`{a}_{\mathscr{l}}` may vary in space but sufficiently
slowly to be considered constant within estimation neighborhoods.

The universal kriging model is the decomposition of the variable :math:`Z(x)`
into the sum:

.. math::
    Z(x)=m(x)+Y(x)

of a smooth deterministic function :math:`m(x)`, describing the systematic aspect of
the phenomenon, and called the drift, and a zero-mean random function :math:`Y(x)`,
called the residual and capturing its erratic fluctuations. Note that the drift
refers to a technically precise notion (the mean of the RF :math:`Z`),
whereas *trend* is a generic term designating a general tendency, a systematic
effect (besides, "trend" may imply an underlying driving force).

In order to minimize :math:`E{({Z}^{*}-{Z}_{0})}^{2}:math:`, we have to make
:math:`{[E({Z}^{*}-{Z}_{0})]}^{2}` zero whatever the unknown coefficients
:math:`{a}_{\mathscr{l}}`, which implies annihilating their factors in the above.
This leads to the set of L+1 conditions:

.. math::
    \sum\limits_{\alpha}{\lambda}_{\alpha}{f}_{\alpha}^{\mathscr{l}}={f}_{0}^{\mathscr{l}}, \quad \mathscr{l}=0,1,\dots,L

that Matheron(1969) called universality conditions, hence the name universal
kriging (UK). *They express that the estimator :math:`{Z}^{*}` is unbiased for *
*all values of :math:`{\alpha}_{\mathscr{l}}`*.

The Universal Kriging System can be expressed as:

.. math::
    \begin{cases}
    \sum\limits_{\beta}{{\lambda}_{\beta}{\sigma}_{\alpha\beta}}+\sum\limits_{\mathscr{l}}{{\mu}_{\mathscr{l}}{f}_{\alpha}^{\mathscr{l}}}={\sigma}_{\alpha0}, &\quad \alpha=1,\dots,N\\
    \sum\limits_{\alpha}{{\lambda}_{\alpha}{f}_{\alpha}^{\mathscr{l}}}={f}_{0}^{\mathscr{l}}, &\quad \mathscr{l}=0,\dots,L
    \end{cases}

In matrix notation the system is of the form :math:`\mathbf{Aw=b}` with the following structure:

.. math::
    \begin{bmatrix}
    \boldsymbol{\Sigma} & \mathbf{F} \\
    {\mathbf{F}}^{'} & 0
    \end{bmatrix}
    \begin{bmatrix}
    \boldsymbol{\lambda} \\
    \boldsymbol{\mu}
    \end{bmatrix}
    =
    \begin{bmatrix}
    {\boldsymbol{\sigma}}_{0}\\
    {\mathbf{f}}_{0}
    \end{bmatrix}

where :math:`\boldsymbol{\Sigma}`, :math:`\boldsymbol{\lambda}` and
:math:`{\boldsymbol{\sigma}}_{0}` are defined as for simple kriging and where

.. math::
    \mathbf{F}=
    \begin{bmatrix}
    1&{f}_{1}^{1}&.&{f}_{1}^{L}\\
    1&{f}_{1}^{1}&.&{f}_{1}^{L}\\
    .&.&.&.\\
    .&.&.&.\\
    .&.&.&.\\
    1&{f}_{1}^{1}&.&{f}_{1}^{L}
    \end{bmatrix}, \quad
    \boldsymbol{\mu}=
    \begin{bmatrix}
    {\mu}_{0}\\
    {\mu}_{1}\\
    .\\
    .\\
    .\\
    {\mu}_{L}
    \end{bmatrix}, \quad
    {\mathbf{f}}_{0}=
    \begin{bmatrix}
    1 \\
    {f}_{0}^{1} \\
    .\\
    .\\
    .\\
    {f}_{0}^{L}
    \end{bmatrix}

Those :math:`1`s in :math:`\mathbf{F}` correspond to OK.


Kriging with an External Drift
------------------------------

Kriging with an external drift variable is an extention of UK. The trend model
is limited to two terms :math:`m(\mathbf{u})={a}_{0}+{a}_{1}{f}_{1}(\mathbf{u})` with
the term :math:`{f}_{1}(\mathbf{u})` set equal to a secondary (external) variable.
The smooth variability of the second variable is deemed related to that of the
primary variable :math:`Z(\mathbf{u})` being estimated.

Let :math:`y(\mathbf{u})` be the secondary variable; the trend model is then:

.. math::
    E\{Z(\mathbf{u})\}=m(\mathbf{u})={a}_{0}+{a}_{1}y(\mathbf{u})

:math:`y(\mathbf{u})` is assumed to reflect the spacial trends of the :math:`z` variability
up to a linear rescalling of units (corresponding to the two parameters :math:`{a}_{0}`
and :math:`{a}_{1}`)

The estimate of the :math:`z` variable and the corresponding system of equations are
identical to the UK estimate and system with K=1, and
:math:`{f}_{1}(\mathbf{u})={y}(\mathbf{u})`

.. math::
    Z_{UK}^{*}(\mathbf{u})=\sum_{\alpha=1}^{n}{{\lambda}_{\alpha}^{UK}(\mathbf{u})Z({\mathbf{u}}_{\alpha})}

.. math::
    \begin{cases}
    \sum_{\beta=1}^{n}{{\lambda}_{\beta}^{UK}(\mathbf{u})C({\mathbf{u}}_{\alpha}-{\mathbf{u}}_{\alpha})} + {\mu}_{0}(\mathbf{u}) + {\mu}_{a}(\mathbf{u})y({\mathbf{u}}_{\alpha}) = C(\mathbf{u}-{\mathbf{u}}_{\alpha}) &\alpha=1,\dots,n\\
    \sum_{\beta=1}^{n}{{\lambda}_{\beta}^{UK}}=1\\
    \sum_{\beta=1}^{n}{{\lambda}_{\beta}^{UK}y({\mathbf{u}}_{\beta})}=y(\mathbf{u})
    \end{cases}

The fundamental (hypothesis) relation must make physical sense.

Two conditions must be met before applying the external drift algorithm:
(1) The external variable must vary smoothly in space, otherwise the resulting
UK system may be unstable; and (2) the external variable must be known at all
locations :math:`{\mathbf{u}}_{\alpha}` of the primary data values and at all locations
:math:`\mathbf{u}` to be estimated.

Block Kriging
-------------

The linearity of the kriging algorithm allows direct estimation of *linear*
averages of the attributes :math:`z(\mathbf{u})`. For example, consider the
estimation of the block average defined as:

.. math::
    z_{V}(\mathbf{u})=\frac{1}{|V|}\int_{V(\mathbf{u})}{z({\mathbf{u}}')d{\mathbf{u}}'}\approx \frac{1}{N}\sum_{j=1}^{N}{z({\mathbf{u}}_{j}^{'})}

where :math:`V(\mathbf{u})` is a block of measure :math:`|V|` centered at u, and the
:math:`{\mathbf{u}}_{j}^{'}` are N points discretizing the volume :math:`V(\mathbf{u})`.

Doing point kriging or block krging only affect the right handside of the
kriging system. Each element in the right hand side matrix is the average
covariance between the sample point and all points in the target block instead
of just the covariance between the sample point and the target point.