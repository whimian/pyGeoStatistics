# -*- coding: utf-8 -*-
"""
A straightforward 2D kriging program

Created on Fri Nov 11 2016
"""
__author__ = "yuhao"

import yaml
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
from itertools import product
import time

from pygeostatistics.yaml_patch import loader_patched


class Krige2d():
    def __init__(self, param_file):
        self.param_file = param_file
        self._read_params()
        self._check_params()
        self.property_name = None
        self.vr = None
        self.maxcov = None
        self.rotmat = None
        self.estimation = None
        self.estimation_variance = None

        self.xdb = None
        self.ydb = None

        self._block_covariance = None
        self._unbias = None

        self._2d = False

    def _read_params(self):
        with open(self.param_file) as fin:
            params = yaml.load(fin, Loader=loader_patched())
            self.datafl = params['datafl']  #: 'testData/test.gslib',
            self.icolx = params['icolx']  #: 1,
            self.icoly = params['icoly']  #: 2,
            self.icolvr = params['icolvr']  #: 0,
            self.tmin = params['tmin']  #: -1.0e21,
            self.tmax = params['tmax']  #: 1.0e21,
            self.idbg = params['idbg']  #: 3,
            self.dbgfl = params['dbgfl']  #: 'kb2d.dbg',
            self.outfl = params['outfl']  #: 'out.dat',
            self.nx = params['nx']  #: 50,
            self.xmn = params['xmn']  #: 0.5,
            self.xsiz = params['xsiz']  #: 1.0,
            self.ny = params['ny']  #: 50,
            self.ymn = params['ymn']  #: 0.5,
            self.ysiz = params['ysiz']  #: 1.0,
            self.nxdis = params['nxdis']  #: 1,
            self.nydis = params['nydis']  #: 1,
            self.ndmin = params['ndmin']  #: ,
            self.ndmax = params['ndmax']  #: ,
            self.radius = params['radius']  #: ,
            self.ktype = params['isk']  #: ,
            self.skmean = params['skmean']  #: ,
            self.nst = params['nst']  #: 1,
            self.c0 = params['c0']  #: 0,
            self.it = params['it']  #: [],
            self.cc = params['cc']  #: [],
            self.azm = params['azm']  #: [],
            self.a_max = params['a_max']  #:[],
            self.a_min = params['a_min']  #: []

    def read_data(self):
        data_list = None
        with open(self.datafl, 'r') as fin:
            data_list = fin.readlines()
        name = data_list[0].strip()
        ncols = int(data_list[1].strip())
        column_name = [item.strip() for item in data_list[2: ncols+2]]
        self.property_name = [item for item in column_name
                              if item not in ['x', 'y', 'z']]
        if 'z' not in column_name:
            self._2d = True
            column_name.append('z')
            data_list = [tuple(item.strip().split() + ['0'])
                         for item in data_list[ncols+2:]]
        else:
            data_list = [tuple(item.strip().split())
                         for item in data_list[ncols+2:]]
        data_dtype = np.dtype({
            'names': column_name,
            'formats': ['f8'] * len(column_name)})
        self.vr = np.array(data_list, dtype=data_dtype)

    def _check_params(self):
        for vtype, a_range in zip(self.it, self.a_max):
            if vtype not in np.arange(1, 6):
                raise ValueError("INVALID variogram number {}".format(vtype))
            if vtype == 4:
                if a_range < 0:
                    raise ValueError("INVALID power variogram")
                elif a_range > 2.0:
                    raise ValueError("INVALID power variogram")
            if vtype == 5:
                raise ValueError("Cannot handle this type of variogram.")

    def _rotation_matirx(self):
        azumth = np.deg2rad(90.0 - np.array(self.ang))
        self.rotmat = np.zeros((4, self.nst))
        self.rotmat[0] = np.cos(azumth)
        self.rotmat[1] = np.sin(azumth)
        self.rotmat[2] = -np.sin(azumth)
        self.rotmat[3] = np.cos(azumth)

    def _max_covariance(self):
        PMX = 9999.0  # max value used for power model
        self.maxcov = self.c0
        for kind, contri in zip(self.it, self.cc):
            if kind == 4:
                self.maxcov += PMX
            else:
                self.maxcov += contri

    def _cova2(self, x1, y1, x2, y2):
        "calculte covariance using provided variogram model"
        PMX = 9999.0  # max value used for power model
        dx = x2 - x1
        dy = y2 - y1
        # check for small distance
        if (dx*dx + dy*dy) < np.finfo("float").eps:
            return self.maxcov
        # for non-zero distance
        cova = 0.0
        for iss in range(self.nst):
            dx1 = dx*self.rotmat[0, iss] + dy*self.rotmat[1, iss]
            dy1 = (dx*self.rotmat[2, iss] + dy*self.rotmat[3, iss]) / \
                  self.anis[iss]
            h = np.sqrt(np.maximum(dx1*dx1 + dy1*dy1, 0))
            if self.it[iss] == 1:  # spherical model
                hr = h/self.a_max[iss]
                if hr < 1:
                    cova += self.cc[iss] * (1 - hr * (1.5 - 0.5 * hr * hr))
            elif self.it[iss] == 2:  # exponential model
                cova += self.cc[iss]*np.exp(-3.0*h/self.a_max[iss])
            elif self.it[iss] == 3:  # gaussian model
                cova += self.cc*np.exp(-3.0 * h * h / \
                                        (self.a_max[iss] * self.a_max[iss]))
            elif self.it[iss] == 4:  # power model
                cova += PMX - self.cc[iss]*(h**(self.a_max[iss]))
        return cova

    def _block_discretization(self):
        """
        Set up the discretization points per block. Figure out how many are
        needed, the spacing, and fill the xdb and ydb arrays with the
        offsets relative to the block center
        """
        xdis = self.xsiz / np.maximum(self.nxdis, 1.0)
        ydis = self.ysiz / np.maximum(self.nydis, 1.0)
        xloc = -0.5*(self.xsiz + xdis)
        yloc = -0.5*(self.ysiz + ydis)
        xdb_temp = np.arange(1, self.nxdis+1, 1) * xdis + xloc
        ydb_temp = np.arange(1, self.nydis+1, 1) * ydis + yloc
        xdb, ydb = np.meshgrid(xdb_temp, ydb_temp)
        self.xdb, self.ydb = xdb.flat, ydb.flat
        # xdb and ydb are nxdis * nydis array

    @property
    def unbias(self):
        "the unbiasedness constraint"
        if self._unbias is None:
            self._unbias = self._cova2(self.xdb[0], self.ydb[0],
                                      self.xdb[0], self.ydb[0])
        return self._unbias

    @property
    def block_covariance(self):
        "the block covariance"
        if self._block_covariance is None:
            self._block_covariance = 0
            if self.ndb <= 1:  # point kriging
                self._block_covariance = self.unbias
            else:  # block kriging
                cov = list()
                for x1, y1 in zip(self.xdb, self.ydb):
                    for x2, y2 in zip(self.xdb, self.ydb):
                        cov.append(self._cova2(x1, y1, x2, y2))
                cov = np.array(cov).reshape((self.ndb, self.ndb))
                cov[np.diag_indices_from(cov)] -= self.c0
                self._block_covariance = np.mean(cov)
        return self._block_covariance

    def _preprocess(self):
        self._read_params()
        # number of points in discretization block
        self.ndb = self.nxdis * self.nydis

        self.anis = np.array(self.a_min)/np.array(self.a_max)
        self.ang = np.array(self.azm)

        self._rotation_matirx()
        self._max_covariance()
        self._block_discretization()
        if self.nxdis == 1 and self.nydis == 1:
            self.block_kriging = False

    def kd2d(self):
        self._preprocess()
        print("Start kriging...")
        # For each target point on the grid
        xloc_temp = np.arange(self.nx) * self.xsiz + self.xmn
        yloc_temp = np.arange(self.ny) * self.ysiz + self.ymn
        yloc_mesh, xloc_mesh = np.meshgrid(yloc_temp, xloc_temp)
        self.estimation = list()
        self.estimation_variance = list()
        num_of_points = self.nx*self.ny
        t1 = time.time()
        ts = 0
        percent_od = 0
        for idx, (xloc, yloc) in enumerate(zip(xloc_mesh.flat, yloc_mesh.flat)):
            ts_1 = time.time()
            # Find the nearest samples within each octant:
            nums, dist = self._search(xloc, yloc)

            ts += time.time() - ts_1
            # is there enough samples?
            if len(dist) < self.ndmin:
                print("Block {},{} not estimated.".format(
                    (xloc-self.xmn)/self.xsiz,
                    (yloc-self.ymn)/self.ysiz))
                self.estimation.append(np.nan)
                self.estimation_variance.append(np.nan)
                continue
            na = dist.shape[0]

            # Put coordinates and values of neighborhood samples into xa,ya,vra
            xa = self.vr['x'][nums]
            ya = self.vr['y'][nums]
            vra = self.vr[self.property_name[0]][nums]
            # handle the situation of only one sample:
            if na == 1:
                est, estv = self._one_sample(xloc, yloc, xa, ya, vra)
                self.estimation.append(est)
                self.estimation_variance.append(estv)
            else:  # many samples
                est, estv = self._many_sample(xloc, yloc, xa, ya, vra)
                self.estimation.append(est)
                self.estimation_variance.append(estv)

            percent = np.round(idx/num_of_points*100, decimals=0)
            dtime = time.time() - t1
            if percent != percent_od:
                print("{}% ".format(percent) +\
                  "."*20 + "{}s elapsed.".format(np.round(dtime, decimals=3)))
            percent_od = percent
        print("Kriging finished.")
        print("Time used for searching: {}s".format(ts))
        self.estimation = np.array(self.estimation).reshape((self.nx, self.ny))
        self.estimation_variance = np.array(
            self.estimation_variance).reshape((self.nx, self.ny))

    def _search(self, xloc, yloc):
        "Search all points return point index and distance to (xloc,yloc)"
        dist = list()
        nums = list()
        # Scan all the samples:
        for idd in range(self.vr.shape[0]):
            dx = self.vr['x'][idd] - xloc
            dy = self.vr['y'][idd] - yloc
            h2 = dx*dx + dy*dy
            if h2 > self.radius*self.radius:
                continue
            # do not consider this sample if there are enough close ones:
            if len(nums) == self.ndmax:
                if h2 >= dist[-1]:
                    continue
                elif h2 < dist[-1]:
                    del nums[-1]
                    del dist[-1]
            # consider this sample (it will be added in the correct location):
            if len(nums) < self.ndmax:
                nums.append(idd)
                dist.append(h2)
        if len(dist) == 0:
            return np.array([]), np.array([])
        else:
            # Sort samples found thus far in increasing order of distance:
            dist = np.array(dist)
            nums = np.array(nums)
            sort_index = np.argsort(dist)
            dist = dist[sort_index]
            nums = nums[sort_index]
            return nums, dist

    def _one_sample(self, xloc, yloc, xa, ya, vra):
        # Left Hand Side Covariance:
        left = self._cova2(xa[0], ya[0], xa[0], ya[0])

        # Right Hand Side Covariance:
        xx = xa[0] - xloc
        yy = ya[0] - yloc
        if not self.block_kriging:  # point kriging
            right = self._cova2(xx, yy, self.xdb[0], self.ydb[0])
        else:  # block kriging
            right = 0.0
            # cb_list = list()
            for i in range(self.ndb):
                right = self._cova2(xx, yy, self.xdb[i], self.ydb[i])
                dx = xx - self.xdb[i]
                dy = yy - self.ydb[i]
                if dx*dx + dy*dy < np.finfo('float').eps:
                    right -= self.c0
            right /= self.ndb

        # Estimation
        if self.ktype == 0:  # Simple kriging
            # Solve for lambda
            s = right / self.block_covariance

            est = s * vra[0] + (1.0 - s) * self.skmean
            estv = self.block_covariance - s * right
            return est, estv

        else:  # Ordinary kriging
            est = vra[0]
            estv = self.block_covariance - 2.0 * right + left
            return est, estv

    def _many_sample(self, xloc, yloc, xa, ya, vra):
        "Solve the Kriging System with more than one sample"
        na = len(vra)
        # number of equations, for simple kriging there're na,
        # for ordinary there're na + 1
        neq = na + self.ktype

        # Establish left hand side covariance matrix:
        left = np.full((neq, neq), np.nan)
        for i, j in product(range(na), range(na)):
            if np.isnan(left[j, i]):
                left[i, j] = self._cova2(xa[i], ya[i], xa[j], ya[j])
            else:
                left[i, j] = left[j, i]

        # Establish the Right Hand Side Covariance:
        right = list()

        for j in range(na):
            xx = xa[j] - xloc
            yy = ya[j] - yloc
            if not self.block_kriging:
                cb = self._cova2(xx, yy, self.xdb[0], self.ydb[0])
            else:
                cb = 0.0
                for i in range(self.ndb):
                    cb += self._cova2(xx, yy, self.xdb[i], self.ydb[i])
                    dx = xx - self.xdb[i]
                    dy = yy - self.ydb[i]
                    if dx*dx + dy*dy < np.finfo('float').eps:
                        cb -= self.c0
                cb /= self.ndb
            right.append(cb)

        if self.ktype == 1:  # for ordinary kriging
            # Set the unbiasedness constraint
            left[neq-1, :-1] = self.unbias
            left[:-1, neq-1] = self.unbias
            left[-1, -1] = 0
            right.append(self.unbias)

        # Solve the kriging system
        s = None
        try:
            s = linalg.solve(left, right)
        except linalg.LinAlgError as inst:
            print("Warning kb2d: singular matrix for block " + \
                  "{},{}".format((xloc-self.xmn)/self.xsiz,
                                 (yloc-self.ymn)/self.ysiz))
            return np.nan, np.nan

        estv = self.block_covariance
        if self.ktype == 1:  # ordinary kriging
            estv -= s[-1]*self.unbias  # s[-1] is mu
        est = np.sum(s[:na]*vra[:na])
        estv -= np.sum(s[:na]*right[:na])
        if self.ktype == 0:  # simple kriging
            est += (1 - np.sum(s[:na])) * self.skmean
        return est, estv

    def view(self, pname=None):
        pname = self.property_name[0] if pname is None else pname
        fig, ax = plt.subplots()
        im = ax.imshow(self.estimation.T, interpolation='nearest',
                       origin='lower',
                       extent=[self.xmn,
                               self.xmn + (self.nx - 1)*self.xsiz,
                               self.ymn,
                               self.ymn + (self.ny - 1)*self.ysiz])
        ax.set_xlabel("X (m)")
        ax.set_ylabel("Y (m)")
        ax.set_title("Estimation")
        ax.set_aspect('equal')
        fig.colorbar(im)
        fig.show()

if __name__ == '__main__':
    test_krige = Krige2d("testData/test_krige2d.par")
    test_krige.read_data()
    test_krige.kd2d()
    test_krige.view()
