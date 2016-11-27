`datafl`: the input data in a simplified Geo_EAS formatted file. The data are
ordered rowwise (X cycles fastest, then Y, then Z).

`nvar` and `ivar(1)` ... `ivar(nvar)`: the number of variables and their
columns in the data file.

`tmin` and `tmax`: all values, regardless of which variable, strictly less
than tmin and greater than or equal to tmax are ignored.

`outfl`: the output variograms are written to a single output file named outfl.
 The output file contains the variograms ordered by direction and then variogram
 type specified in the parameter file (the directions cycle fastest then the
 variogram number.) For each variogram there is a one-line description and then
 nlag lines each with the following:

 1. lag number (increasing from 1 to nlag).
 2. average separation distance for the lag.
 3. the *semivariogram* value (whatever type was specified).
 4. number of pairs for the lag.
 5. mean of the data contributing to the tail.
 6. mean of the data contributing to the head.
 7. the tail and head variances (for the correlogram).

`igrid`: the grid or realization number. Recall that realization or grids are
written on after another; therefore, if igrid=2 the input file must contain
at least 2 nx ny nz values and the second set of nx ny nz values will be taken
as the second grid.

`nx`, `xmn`, `xsiz`: definition of the grid system (x axis)

`ny`, `ymn`, `ysiz`: definition of the grid system (y axis)

`mz`, `zmn`, `zsiz`: definition of the grid system (z axis)