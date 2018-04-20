CoKriging
=========

Intro
-----

The term kriging is traditionally reserved for liear regression using data on
the same attribute as that being estimated. For example, an unsampled porosity
value :math:`z(u)` is estimated from neighboring porosity sample values defined on
the same volume support.

The term cokriging is reserved for linear regression that also uses data defined
on different attributes. For example, the porosity values :math:`z(u)` may be estimated
from combination of porosity samples and related acoustic data values.

In the case of a single secondary variable (:math:`Y`), the ordinary cokriging
estimator of :math:`Z(\mathbf{u})` is written:

.. math::
    Z_{COK}^{*}(\mathbf{u})
    =\sum_{{\alpha}_{1}=1}^{{n}_{1}}{{\lambda}_{{\alpha}_{1}}(\mathbf{u})Z({\mathbf{u}}_{{\alpha}_{1}})}
    +\sum_{{\alpha}_{2}=1}^{{n}_{2}}{{\lambda}_{{\alpha}_{2}}^{'}(\mathbf{u})Y({\mathbf{u}}_{{\alpha}_{2}}^{'})}

where the :math:`{\lambda}_{{\alpha}_{1}}` are the weights applied to the :math:`{n}_{1}`
:math:`z` samples and the :math:`{\lambda}_{{\alpha}_{2}}^{'}` are the weights applied to
the :math:`n_2` `y` samples.

Kriging requires a model for the :math:`Z` covariance. Cokriging requires a joint
model for the matrix of covariance functions including the :math:`Z` covariance
:math:`C_{Z}(\mathbf{h})`, the :math:`Y` covariance :math:`C_{Y}(\mathbf{h})`, the cross :math:`Z-Y`
covariance :math:`C_{ZY}(\mathbf{h})=Cov\{Z(\mathbf{u}),Y(\mathbf{u+h})\}`, and the
cross :math:`Y-Z` covariance :math:`C_{YZ}(\mathbf{h})`

The covariance matrix requires :math:`K^2` covariance functions when :math:`K` different
variables are considered in a cokriging exercise. The inference becomes
extremely demanding in terms of data and the subsequent joint modeling is
particularly tedious. This is the main reason why cokriging has not been
extensively used in practice. Algorithms such as kriging with an external
drift and collocated cokriging have been developed to shortcut the tedious
inference and modeling process required by cokriging.

Ordinary Cokriging
------------------

The sum of the weights applied to the primary variable is set to one, and
the sum of the weigths applied to any other variable is set to zero. In the
case of two variables, these two conditions are:

.. math::
    \begin{cases}
    \sum\limits_{{\alpha}_{1}}^{}{{\lambda}_{{\alpha}_{1}}(\mathbf{u})}=1\\
    \sum\limits_{{\alpha}_{2}}^{}{{\lambda}_{{\alpha}_{2}}(\mathbf{u})}=0
    \end{cases}

The problem with this traditional formalism is that the second condition tends
to limit severely the influence of the secondary variables.

Standardized Ordinary Cokriging
-------------------------------

Often, a better approach consists of creating new secondary variables with the
same mean as the primary variable. Then all the weights are constrained to
sum to one.

In the case of two variables, the expression could be written as:

.. math::
    Z_{COK}^{*}(\mathbf{u})
    =\sum_{{\alpha}_{1}=1}^{{n}_{1}}{{\lambda}_{{\alpha}_{1}}(\mathbf{u})Z({\mathbf{u}}_{{\alpha}_{1}})}
    +\sum_{{\alpha}_{2}=1}^{{n}_{2}}{{\lambda}_{{\alpha}_{2}}^{'}(\mathbf{u})[Y({\mathbf{u}}_{{\alpha}_{2}}^{'})+{m}_{Z}-{m}_{Y}]}


with a single condition:

.. math::
    \sum_{{\alpha}_{1}=1}^{{n}_{1}}{{\lambda}_{{\alpha}_{1}}}(\mathbf{u})+\sum_{{\alpha}_{2}=1}^{{n}_{2}}{{\lambda}_{{\alpha}_{2}}}(\mathbf{u})=1

where :math:`m_Z=E\{Z(u)\}` and :math:`m_Y=E\{Y(u)\}` are stationary means
of :math:`Z` and :math:`Y`.

Simple Cokriging
----------------

There is no constraint on the weights. Just like simple kriging, this version
of cokriging requires working on data residuals or equivalently, on variables
whose means have all been standardized to zero. This is the case when applying
simple cokriging in an MG approach (the normal score transforms of each variable
have a stationary mean of zero).

Collocated Cokriging
--------------------

A reduced form of cokriging consists of retaining only the collocated
variable :math:`y(\mathbf{u})`, provided that it is availible at all locations
:math:`\mathbf{u}` being estimated. The cokriging estimator is written as:

.. math::
    Z_{COK}^{*}(\mathbf{u})
    =\sum_{{\alpha}_{1}=1}^{{n}_{1}}{{\lambda}_{{\alpha}_{1}}(\mathbf{u})Z({\mathbf{u}}_{{\alpha}_{1}})}
    +{\lambda}^{'}(\mathbf{u})Y(\mathbf{u})

The corresponding cokriging system requires knowledge of only the :math:`Z`
covariance :math:`C_{Z}(\mathbf{h})` and the :math:`Z-Y` cross-covariance
:math:`C_{ZY}(\mathbf{h})`. The latter can be approximated through the following model:

.. math::
    C_{ZY}(\mathbf{h})=B\cdot C_{Z}(\mathbf{h}),\quad\forall \mathbf{h}

where :math:`B=\sqrt{C_Y(0)/C_Z(0)}\cdot{\rho}_{ZY}(0)`, :math:`C_Z(0)`, :math:`C_Y(0)` are
the variances of Z and Y, and :math:`{\rho}_{ZY}(0)` is the linear coefficient
of correlation of collocated z-y data.
