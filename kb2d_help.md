# Description:

This is a straightforward 2-D simple and ordinary kriging subroutine that can
be used as is or as a basis for custom kriging programs.


# Paramters:

- `datafl`: the input data in a simplified Geo-EAS formatted file.

- `icolx`, `icoly` and `icolvr`: the columns for the x and y coordinates, and
the variable to be kriged.

- `tmin` and `tmax`: all values strictly less than `tmin` and greater than or
equal to `tmax` are ignored.

- `idbg`: an iteger debugging level between 0 and 3. The higher the debugging
level, the more output. Normally level 0 or 1 should be chosen. If there are
suspected problems, or if you would like to see the actual kriging matrices,
level 2 or 3 can be chosen. It is advisable to restrict the actual number of
points being estimated when the debugging level is high (the debugging file
can become extremely large.)

- `dbgfl`: the debugging output is written to this file.

- `outfl`: the output grid is written to this file. The output file will contain
both the kriging estimates and the kriging variance for all points/blocks. The
output grid cycles fastest on *x* and then *y*.

- `nx`, `xmn`, `xsiz`: definition of the grid system (*x* axis).

- `ny`, `ymn`, `ysiz`: definition of the grid system (*y* axis).

- `nxdis` and `nydis`: the number of discretization points for a block. If both
`nxdis` and `nydis` are set to 1, then point kriging is performed.

- `ndmin` and `ndmax`: the minimum and maximum number of data points to use for
kriging a block.

- `radius`: the maximum isotopic search radius.

- `isk` and `skmean`: if `isk`=0, then simple kriging will be performed with a
mean of `skmean`.

- `nst` and `c0`: the number of variogram structures and the isotopic nugget
constant. The nugget constant does not count as a structure.

- For each of the `nst` nested structures one must define `it`, the type of
structures; `cc`, the *c* parameter; `azm`, the maximum range; `a_max`, the
maximum range; and `a_min`, the minimum range. A detailed description of these
parameters is given in section II.3.
