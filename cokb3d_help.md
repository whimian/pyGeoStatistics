# Description:

A cokriging program for a points or blocks on a regular grid.

# Parameters:

- `datafl`: the input data in a simplified Geo-EAS formatted file.

- `nvar`: the number of variables (primary plus all secondary).
For example, `nvar`=2 if there is only one secondary variable.

- `icolx`, `icoly`, `icolz` and `icolvr()`: the columns for the x, y and
z coordinates, the primary variable to be kriged, and all secondary variables.

- `tmin` and `tmax`: all values (for all variables) strictly less than `tmin`
and greater than or equal to `tmax` are ignored.

- `coloc`: set to 1 if performing co-located cokriging with a gridded secondary
 variable, otherwise set to 0.

- `secfl`: if co-located cokriging, the file with
the gridded secondary variable.

- `icolsec`: if co-located cokriging, the column number for the secondary
variable in `secfl`

- `idbg`: an integer debugging level between 0 and 3. The higher the
debugging level the more output. Normally level 0 or 1 should be chosen.

- `dbgfl`: the debugging output is written to this file.

- `outfl`: the output grid is written to this file. The output file will
contain both the kriging estimate and the kriging variance for all
points/blocks. The output grid cycles fastest on x then y then z.

- `nx`, `xmn`, `xsiz`: definition of the grid system (x axis).

- `ny`, `ymn`, `ysiz`: definition of the grid system (y axis).

- `nz`, `zmn`, `zsiz`: definition of the grid system (z axis).

- `nxdis`, `nydis` and `nzdis`: the number of discretization points for
a block. If `nxdis`, `nydis` and `nzdis` are set to 1 then point cokriging
is performed.

- `ndmin`, `ndmaxp` and `ndmaxs`: the minimum and maximum number of primary
data, and the maximum number of secondary data (regardless of which secondary
variable) to use for kriging a block.

- `pradius_hmax', 'pradius_hmin` and `pradius_vert`:
search radii for primary data

- `sradius_hmax`, `sradius_hmin` and `sradius_vert`:
search radii for secondary data (same for all types)

- `sangle`, `sangle1` and `sangle2`: the angles defining the common
orientation of the search ellipsoids for primary and secondary data

- `ktype: the kriging type must be specified:
    - 0 = simple cokriging;
    - 1 = standardized ordinary cokriging with re-centered variables
    and a single unbiasedness constraint;
    - 2 = traditional ordinary cokriging.

- `mean()`: the mean of the primary and all secondary variables are required
 input if either simple cokriging or standardized ordinary cokriging are used.
 The program calculates the data residuals from these means.

The direct and cross variograms may be specified in any order; they are
specified according to the variable number Variable `1` is the primary
(regardless of its column ordering in the input data files)
and the secondary variables are numbered from `2` depending on their
 ordered specification in `icolvr()`. It is unnecessary to specify
 the `j` to `i` cross variogram if the `i` to `j` cross variogram
 has been specified; the cross variogram is expected to be symmetric
 (as from theory). For each `i` to `j` variogram the following are required:

- `nst` and `c0`: the number of variogram structures and the isotropic
nugget constant. The nugget constant does not count as a structure.

- For each of the `nst` nested structures one must define `it` the type of
structure (the power model is not allowed);
    - `cc`, the c parameter;
    - `ang1`, `ang2`, `ang3` the angles defining the geometric anisotropy;
    - `aa_hmax`, the maximum horizontal range;
    - `aa_hmin`, the minimum horizontal range; and
    - `aa_vert`, the vertical range.

# Application notes:

The construction of the cokriging matrix requires the **linear model of
coregionalization**. The input variogram parameters are checked for
positive definiteness. *The power model is not allowed.*

A specific search is done for secondary data (same for all secondary)
allowing the option of collocated cokriging.

A cokriging program for scattered points and cross validation is not provided;
programs `cokb3d` and `kt3d` could be combined for this purpose.
