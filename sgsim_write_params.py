# -*- coding: utf-8 -*-
"""
Created on Sun Apr 2 2017
"""
import json

__author__ = "yuhao"

params = {
    'datafl': 'testData/test.gslib',
    'icolx': 0,
    'icoly': 1,
    'icolz': -1,
    'icolvr': 2,
    'icolsec': -1,  # for external drift if used
    'icolwt': -1,  # declustering weights

    # data limits
    'tmin': -1.0e21,
    'tmax': 1.0e21,
    'itrans': True,  # boolean
    # output file for transformation table if transformation is needed
    'transfl': 'sgsim.trn',
    'ismooth': False,  # boolean
    # file with values used for transformation to normal scores
    'smthfl': 'histsmth.out',
    'icolsvr': 0,
    'icolswt': 1,
    # allowable data values used for backtransform
    'zmin': 0,
    'zmax': 30,
    # lower and upper tail model specification for backtransform
    'ltail': 1,
    'ltpar': 0,
    'utail': 1,
    'utpar': 15,
    # debug and output data file
    'idbg': 3,
    'dbgfl': 'sgsim.dbg',
    'outfl': 'sgsim.out',
    'nsim': 1,  # number of simulation
    # Grid definition
    'nx': 98,
    'xmn': 100,
    'xsiz': 200,
    'ny': 79,
    'ymn': 100,
    'ysiz': 200,
    'nz': 1,
    'zmn': 0,
    'zsiz': 200,
    'seed': 1,  # random seed
    # maximum and minimum data points used in kriging
    'ndmin': 1,
    'ndmax': 30,
    'nodmax': 12,  # previously simulated nodes to use
    'sstrat': 0,  # search strategy
    'multgrid': False,  # boolean
    'nmult': 2,  # scalar
    'noct': 0,  # maximum number to retain from an octant
    # search radii
    'radius_hmax': 4000,
    'radius_hmin': 4000,
    'radius_vert': 0,
    # search anisotropy angles
    'sang1' : 0,
    'sang2' : 0,
    'sang3' : 0,
    # size of covariance lookup table size
    'mxctx': 30,
    'mxcty': 30,
    'mxctz': 30,
    # kriging type
    'ikrige': 0,
    # self.skmean = params['skmean']
    'rho': 0.7, # correlation coefficient for COCOK
    'varred': 0.1, # variance reduction factor for COCOK
    'secfl': 'ydata.dat',
    'icollvm': 4,
    # Vairography definition
    'nst': 1,
    'c0': 0.05,
    'it': [1],
    'cc': [0.65],
    'ang1': [0],
    'ang2': [0],
    'ang3': [0],
    'aa_hmax': [3715.9],
    'aa_hmin': [3715.9],
    'aa_vert': [3715.9]
}

with open('testData/test_sgsim.par', 'w') as fout:
    fout.write(json.dumps(params, sort_keys=True, indent=4))