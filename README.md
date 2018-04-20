# pyGeoStatistics

![status](https://img.shields.io/badge/status-alpha-green.svg)
[![Documentation Status](https://readthedocs.org/projects/pygeostatistics/badge/?version=latest)](http://pygeostatistics.readthedocs.io/en/latest/?badge=latest)

A collection of python routines (accelerated with [Numba](https://github.com/numba/numba))
and jupyter notebooks for geostatistics,
which is immensely inspired by gslib (in Fortran).

# Usage

Every routine reads its parameters from a parameter file written in `json`.
All parameters including input/output file path need to be specified in these parameter
files.

I've created scripts that assist in creating parameter files, they could be
found in `\parameters` folder.

I tried to adhere to the naming convention of `gslib` when it comes to parameter
names.

Markdown files describing parameters needed for each routine are in
`\gslib_help`.

## Example:

```Python
from pygeostatistics import Sgsim

sgsimulator = Sgsim("testData/test_sgsim.par")
sgsimulator.simulate()
```

# Routines

- `eda.py`: exploratory data anaylysis.

- [`nst.py`](#normal-score-transform-nstpy): apply normal score transform to data.

- `gam.py`: calculate variogram for regular data.

- `gamv.py`: calculate variogram for irregular data.

- `sa.ipynb`: interactive structural analysis.

- [`krige2d.py`](#2d-kriging-krige2dpy): kriging 2d data.

  - Simple Kriging
  - Ordinary Kriging

- [`krige3d.py`](#3d-kriging-krige3dpy): kriging 3d data.

  - Simple Kriging
  - Ordinary Kriging
  - Universal Kriging (Kriging with a Trend)
  - Kriging the Trend
  - Kriging with External drift
  - SK with non-stationary drift

- [`sgsim.py`](#sequential-gaussian-simulation-sgsimpy): Sequential Gaussian Simulation.

# Other Utilities

- `super_block.py`: Class for performing super block search used in kriging.
  - used in `krige3d.py`
  - used in `sgsim.py`

- `normal_score_transform.py`: Class for NST used in Gaussian Simulation.
  - used in `sgsim.py`

# Documentation

For full documentation, including installation, tutorials and PDF documents, please see http://pygeostatistics.readthedocs.io/.
