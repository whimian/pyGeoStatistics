# -*- coding: utf-8 -*-
"""
Created on Tue Nov 2016
"""
import json

__author__ = "yuhao"

params = {
    'datafl': 'testData/test.gslib',
    'icolx': 0,
    'icoly': 1,
    'icolz': 2,
    'icolvr': 3,
    'icolsec': 4,

    'tmin': -1.0e21,
    'tmax': 1.0e21,

    'option': 0,
    'jackfl': 'jackfl.dat',
    'jicolx': 0,
    'jicoly': 1,
    'jicolz': 2,
    'jicolvr': 3,
    'jicolsec': 4,

    'idbg': 3,
    'dbgfl': 'kb2d.dbg',
    'outfl': 'out.dat',

    'nx': 98,
    'xmn': 100,
    'xsiz': 200,
    'ny': 79,
    'ymn': 100,
    'ysiz': 200,
    'nz': 79,
    'zmn': 200,
    'zsiz': 200,

    'nxdis': 1,
    'nydis': 1,
    'nzdis': 1,

    'ndmin': 1,
    'ndmax': 10,

    'noct': 0,
    'radius_hmax': 4000,
    'radius_hmin': 4000,
    'radius_vert': 4000,
    'sang1' : 30,
    'sang2' : 30,
    'sang3' : 30,

    'ikrige': 0,
    'skmean': 14.69588,

    'idrift': [0, 0, 0, 0, 0, 0, 0, 0, 0],
    'itrend': 0,
    'secfl': 'secfl.dat',
    'iseccol': 3,

    'nst': 1,
    'c0': 0.05,
    'it': [0],
    'cc': [0.65],
    'ang1': [0],
    'ang2': [0],
    'ang3': [0],
    'aa_hmax': [3715.9],
    'aa_hmin': [3715.9],
    'aa_vert': [3715.9],
}

with open('testData/test_krige3d.par', 'w') as fout:
    fout.write(json.dumps(params, sort_keys=True, indent=4))
