Simulation
==========

Simulation differs from kriging or any interpolation algorithm, in two major aspects:

1. In most interpolation algorithms, including kriging, the goal is to provide
a "best", hence unique, local estimate of the variable or any of its trend
components without specific regard to the resulting spatial statistics of the
estimates taken together. **In simulation, reproduction of global features **
**(texture) and statistics (histogram, covariance) take precedence over local **
**accuracy.** Kriging provides a set of local representations, say
:math:`z^{*}(\mathbf{u}),\mathbf{u}\in A`, where local accuracy prevails.
Simulation provides alternative global representations, :math:`z^{(l)}(u),u\in A`,
where reproduction of patterns of spatial continuity prevails.

2. Except if a Gaussian model for errors is assumed, kriging provides only an
incomplete measure of local accuracy, and no appreciation of joint accuracy
when several locations are considered together.
**Simulations are designed specifically to provide such measures of accuracy, both local and involving several locations.**
These measures are given by the differences between :math:`L` alternative
simulated vlaues at any location (local accuracy) or the :math:`L` alternative
simulated fields (global or joint accuracy).

Different simulation algorithms impart different global statistics and spatial
features on each realization, For example, simulated categorical values can be
made to honor specific geometrical patterns as in *object-based simulation* or
the covariance of simulated continuous values can be made to honor a prior
covariance model as for *Gaussian-related simulations*. A hybrid approach
could be considered to generate numerical models that reflect widely different
types of features. For example, one may start with an object-based process or
categorical *indicator simulation* to generate the geometric architecture of
the various lithofacies, following with a Gaussian algorithm to simulate the
distribution of continuous petrophysical properties within each sperate lithofacies,
then a simulated annealing process could be used to modify locally the petrophysical
properties to match, say, well test data.