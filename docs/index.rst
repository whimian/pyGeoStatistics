.. Pore Pressure Prediction documentation master file, created by
   sphinx-quickstart on Thu Oct 26 10:58:51 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

===============
pyGeoStatistics
===============

Overview
========

A collection of python routines (accelerated with Numba) and jupyter notebooks
for geostatistics, which is immensely inspired by gslib (in Fortran).


Usage
=====
Every routine reads its parameters from a parameter file written in json. All parameters including input/output file path need to be specified in these parameter files.

I've created scripts that assist in creating parameter files, they could be found in \parameters folder.

I tried to adhere to the naming convention of gslib when it comes to parameter names.

Markdown files describing parameters needed for each routine are in \gslib_help.

Example:
========
::

    from pygeostatistics import Sgsim

    sgsimulator = Sgsim("testData/test_sgsim.par")
    sgsimulator.simulate()

Contribute
==========
- Issue Tracker: https://github.com/whimian/pyGeoStatistics/issues

- Source Code: https://github.com/whimian/pyGeoStatistics

License
=======
`MIT <https://github.com/whimian/pyGeoStatistics/blob/master/LICENSE>`_

.. toctree::
    :maxdepth: 1
    :caption: Getting Started
    :hidden:

    install

.. toctree::
  :maxdepth: 1
  :caption: Geostatistics Basics
  :hidden:

  basics/kriging
  basics/cokriging
  basics/simulation

.. toctree::
  :maxdepth: 1
  :caption: Tutorials
  :hidden:

  tutorials/EDA
  tutorials/NST
