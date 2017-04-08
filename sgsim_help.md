# Description:

Sequential Gaussian simulation program

# Parameters:

- `datafl`: the input data in a simplified Geo-EAS formatted file. If this file does not exist then an unconditional simulation will be generated.

- `icolx`, `icoly`, `icolz`, `icolvr`, `icolwt` and `icolsec`: the column numbers for
the x, y and z coordinates, the variable to be simulated, the declustering
weight, and the secondary variable (e.g., for external drift if used).
One or two of the coordinate column numbers can be set to zero which indicates
that the simulation is 2-D or 1-D. For equal weighting, set `icolwt` to zero.

- `tmin` and `tmax`: all values strictly less than `tmin`
and strictly greater than `tmax` are ignored.

- `itrans`: if set to 0 then no transformation will be performed;
the variable is assumed already standard normal (the simulation results
will also be left unchanged). If `itrans`=1, transformations are performed.

- `transfl`: output file for the transformation table if transformation
is required (igauss=0).

- `ismooth`: if set to 0, then the data histogram, possibly with declustering
weights is used for transformation, if set to 1, then the data are transformed
according to the values in another file (perhaps from histogram smoothing).

- `smthfl`: file with the values to use for transformation to normal scores
 (if ismooth is set to 0).

- `icolvr` and `icolwt`: columns in smthfl for the variable and the
 declustering weight (set to 1 and 2 if smthfl is the output from histsmth).

- `zmin` and `zmax`: the minimum and maximum allowable data values.
These are used in the back transformation procedure.

- `ltail` and `ltpar` specify the back transformation implementation in
the lower tail of the distribution:
    - `ltail`=1 implements linear interpolation to the lower limit zmin,
    - `ltail`=2 implements power model interpolation, with w=`ltpar`,
    to the lower limit `zmin`.
    - The middle class interpolation is linear.

- `utail` and `utpar` specify the back transformation implementation in the
upper tail of the distribution:
    - `utail`=1 implements linear interpolation to the upper limit `zmax`,
    - `utail`=2 implements power model interpolation, with w=`utpar`,
    to the upper limit `zmax`,
    - `utail`=4 implements hyperbolic model extrapolation with w=`utpar`.
    - The hyperbolic tail extrapolation is limited by zmax.

- `idbg`: an integer debugging level between 0 and 3.
The larger the debugging level the more information written out.

- `dbgfl`: the file for the debugging output.

- `outfl`: the output grid is written to this file. The output file will
contain the results, cycling fastest on x then y then z then simulation by simulation.

- `nsim`: the number of simulations to generate.

- `nx`, `xmn`, `xsiz`: definition of the grid system (x axis).

- `ny`, `ymn`, `ysiz`: definition of the grid system (y axis).

- `nz`, `zmn`, `zsiz`: definition of the grid system (z axis).

- `seed`: random number seed (a large odd integer).

- `ndmin` and `ndmax`: the minimum and maximum number of original data that
should be used to simulate a grid node. If there are fewer than
`ndmin` data points the node is not simulated.

- `ncnode`: the maximum number of previously simulated nodes to use
for the simulation of another node.

- `sstrat`: if set to 0, the data and previously simulated grid nodes are
searched separately: the data are searched with a *super block* search and
the previously simulated nodes are searched with a *spiral search*.
If set to 1, the data are relocated to grid nodes and a spiral search is used
and the parameters `ndmin` and `ndmax` are not considered.

- `multgrid`: a multiple grid simulation will be performed if this is set to 1
(otherwise a standard spiral search for previously simulated nodes is considered).

- `nmult`: the number of multiple grid refinements to consider
(used only if multgrid is set to 1).

- `noct`: the number of original data to use per octant. If this parameter is
set less than or equal to 0, then it is not used; otherwise, it overrides
the ndmax parameter and the data is partitioned into octants and
the closest `noct` data in each octant is retained for the simulation of a grid node.

- `radius_hmax`, `radius_hmin` and `radius_vert`:
the search radii in the maximum horizontal direction,
minimum horizontal direction, and vertical direction (see angles below).

- `sang1`, `sang2` and `sang3`:
the angle parameters that describe the orientation of the search ellipsoid.

- `ktype`: the kriging type (
    0 = simple kriging,
    1 = ordinary kriging,
    2 = simple kriging with a locally varying mean,
    3 = kriging with an external drift, or
    4 = collocated cokriging with one secondary variable)
    used throughout the loop over all nodes.
    SK is required by theory; only in cases where the number of original data
    found in the neighborhood is large enough can OK be used without
    the risk of spreading data values beyond their range of influence.

- `rho`: correlation coefficient to use for collocated cokriging
(used only if ktype = 4).

- `secfl`: the file for the locally varying mean, the external drift variable,
or the secondary variable for collocated cokriging (the secondary variable
must be gridded at the same resolution as the model being constructed by sgsim).

- `nst` and `c0`: the number of semivariogram structures and the isotropic nugget constant.
For each of the nst nested structures one must define
    - `it`, the type of structure;
    - `cc`, the c parameter;
    - `ang1`, `ang2`, `ang3`, the angles defining the geometric anisotropy;
    - `aa_hmax`, the maximum horizontal range;
    - `aa_hmin`, the minimum horizontal range; and
    - `aa_vert`, the vertical range.

# Application notes:

This program requires standard normal data and writes standard normal simulated values. Normal score transform and back transform are to be performed outside of this program
Recall that the power model is not a legitimate model for a multiGaussian phenomenon and it is not allowed in `sgsim`
The semivariogram model is that of the normal scores. The kriging variance is directly interpreted as the variance of the conditional distribution; consequently, the nugget constant `c0` and `c` (sill) parameters should add to 1.0.
