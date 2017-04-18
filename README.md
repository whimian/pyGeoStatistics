# pyGeoStatistics

A collection of python routines and jupyter notebooks for geostatistics,
which is immensely inspired by gslib (in Fortran).

**Caution:** its development just started, hence super buggy.

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

- ['sgsim.py'](#sgsimpy): Sequential Gaussian Simulation.

# Other Utilities

- `super_block.py`: Class for performing super block search used in kriging.
    - used in `krige3d.py`

# Normal Score Transform (nst.py)

# 2D Kriging (krige2d.py)

## Instructions:

This program supports the basic Simple Kriging and Ordinary Kriging.

1. Simple Kriging (SK)
    - set `ikrige = 0`
2. Ordinary Kriging (OK)
    - set `ikrige = 1`

# 3D Kriging (krige3d.py)

## Instructions:

This program can perform all four major kinds of kriging:

1. Simple Kriging (SK)
    - set `ikrige = 0`
2. Ordinary Kriging (OK)
    - set `ikrige = 1`
3. Universal Kriging (UK) a.k.a. Kriging with a Trend
    - set `ikrige = 1` and set `idrift`
4. Kriging with External Drift (KED)
    - set `ikrige = 3`
    - KED requires external variable values at both sample points (in `datafl`)
    and points being estimated (in `secfl`)

and two kinds of additional functionalities:

1. non-stationary simple kriging with known means
    - set `ikrige = 2`
    - requires mean values at each point being estimated (in `secfl`)
2. kriging the trend
    - set `itrend = True`

*A detailed mathematical explanation of kriging could be found in [kriging.ipynb](kriging.ipynb)*
## Validation

Two validation methods are supported in krige3d, which are cross-validation
with data in import data file `datafl` and jackknife validation with data
in another file `jackfl`.

Validation data are read on the fly when doing cross-validation or jackknife
validation.

## Paramters
- `ikrige` and `skmean`:
    - if `ikrige` is set to 0 then stationary simple kriging with (`skmean`)
    will be performed,
    - if `ikrige` is set to 1 then ordinary kriging will be performed,
    - if `ikrige` is set to 2 then nonz-stationary simple kriging with means
    taken from `secfile` will be performed,
    - if `ikrige` is set to 3 then kriging with an external drift will be
    performed.
    - Note that power law variogram models (`it`=4) are not allowed with
    simple kriging.

# Sequential Gaussian Simulation (sgsim.py)
