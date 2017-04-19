# -*- coding: utf-8 -*-
"""
Created on Sun Nov 06 18:19:28 2016
"""
__author__ = "yuhao"

import os
import json

PARAMS = {
    'datafl': 'testData/test.gslib',
    'icolx': 1,
    'icoly': 2,
    'icolz': 0,
    'nvar': 1,
    'ivar': [3, 4],
    'tmin': -1.0e21,
    'tmax': 1.0e21,
    'outfl': 'out.dat',
    'nlag': 20,
    'xlag': 500.0,
    'xltol': 300.0,
    'ndir': 1,
    'azm': [0.0],  # [0.0, 0.0, 90.],
    'atol': [90.0],  # [90.0, 22.5, 22.5],
    'bandwh': [200.0],   # [200.0, 200.0, 200.0],
    'dip': [0.0],  # [0.0, 0.0, 0.0],
    'dtol': [90.0],  # [90.0, 22.5, 22.5],
    'bandwd': [200.0],  # [200.0, 200.0, 200.0],
    'standardize': False,
    'nvarg': 3,
    'ivtail': [1, 1, 2],
    'ivhead': [1, 1, 2],
    'ivtype': [1, 3, 1]
}

PARENT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.path.pardir)
PARAM_DIR = os.path.join(PARENT_DIR, 'testData')

with open(os.path.join(PARAM_DIR, 'xihuSmall_sparse_gamv.par'), 'w') as fout:
    fout.write(json.dumps(PARAMS, sort_keys=True, indent=4))
