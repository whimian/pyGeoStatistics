# -*- coding: utf-8 -*-
"""
Created on Thu Nov 03 20:53:02 2016

@author: yuhao
"""

import json

params = {
    'datafl': 'testData/xihu_sparse.gslib',
    'nvar': 1,
    'ivar': [1, 2],
    'tmin': -1.0e21,
    'tmax': 1.0e21,
    'outfl': 'gam.out',
    'igrid': 1,
    'nx': 15,
    'xmn': 0.5,
    'xsiz': 5,
    'ny': 23,
    'ymn': 0.5,
    'ysiz': 5,
    'nz': 161,
    'zmn': 0.5,
    'zsiz': 0.5,
    'ndir': 2,
    'nlag': 10,
    'ixd': [1, 0],
    'iyd': [0, 1],
    'izd': [0, 0],
    'standardize': True,
    'nvarg': 5,
    'ivtail': [1, 1, 2, 2, 1],
    'ivhead': [1, 1, 2, 2, 1],
    'ivtype': [1, 3, 1, 3, 9]
    }
with open('xihuSmall_sparse_gam.par', 'w') as fout:
    fout.write(json.dumps(params, sort_keys=True, indent=4))

#with open('abc.txt', 'r') as fin:
#    kkk = json.load(fin)
