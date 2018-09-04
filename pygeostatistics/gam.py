# -*- coding: utf-8 -*-
"""
Compute Variogram on regularly spaced data

Created on Tue Nov 03 2016
"""
from __future__ import division, print_function
import json
import numpy as np
import matplotlib.pyplot as plt

__author__ = "yuhao"


class Gam():
    def __init__(self, param_file):
        self.param_file = param_file
        self._read_params()
        self._check_params()
        self.vr = None
        self.gam = None
        self.npair = None
        self.directions = None
        self.mean = None
        self.variance = None

    def _read_params(self):
        with open(self.param_file) as fin:
            params = json.load(fin)
            self.datafl = params['datafl']
            self.nvar = params['nvar']
            #   'ivar'
            self.tmin = params['tmin']
            self.tmax = params['tmax']
            #   'outfl'
            #   'igrid'
            self.nx = params['nx']
            self.xmn = params['xmn']
            self.xsiz = params['xsiz']
            self.ny = params['ny']
            self.ymn = params['ymn']
            self.ysiz = params['ysiz']
            self.nz = params['nz']
            self.zmn = params['zmn']
            self.zsiz = params['zsiz']
            self.ndir = params['ndir']
            self.nlag = params['nlag']
            self.ixd = params['ixd']
            self.iyd = params['iyd']
            self.izd = params['izd']
            self.standardize = params['standardize']
            self.nvarg = params['nvarg']
            self.ivtail = params['ivtail']
            self.ivhead = params['ivhead']
            self.ivtype = params['ivtype']

    def read_data(self):
        data_list = list()
        name = None
        ncols = None
        column_name = None
        data_dtype = None
        with open(self.datafl, 'r') as fin:
            name = fin.readline()
            ncols = int(fin.readline())
            column_name = list()
            for i in range(ncols):
                column_name.append(fin.readline().rstrip('\n'))
            data_dtype = np.dtype({
                'names': column_name,
                'formats': ['f8'] * ncols})
            for line in fin:
                data = line.split()
                for i in range(ncols):
                    data[i] = float(data[i])
                data = tuple(data)
                data_list.append(data)
        input_data = np.array(data_list, dtype=data_dtype)
        self.vr = input_data[column_name[-1]]
        self.vr = self.vr.reshape((self.nx, self.ny, self.nz))

    def _check_params(self):
        try:
            # self.datafl
            if not isinstance(type(self.nvar), int) != 'int' or self.nvar < 1:
                raise ValueError("wrong value with number of variables")
            #   'ivar'
            # self.tmin = params['tmin']
            # self.tmax = params['tmax']
            # #   'outfl'
            # #   'igrid'
            # self.nx = params['nx']
            # self.xmn = params['xmn']
            # self.xsiz = ['xsiz']
            # self.ny = params['ny']
            # self.ymn = params['ymn']
            # self.ysiz = ['ysiz']
            # self.nz = params['nz']
            # self.zmn = params['zmn']
            # self.zsiz = ['zsiz']
            # self.ndir = ['ndir']
            # self.nlag = ['nlag']
            # self.ixd = params['ixd']
            # self.iyd = params['iyd']
            # self.izd = params['izd']
            # #   'standardize':True,
            # self.nvarg = ['nvarg']
            # self.ivtail = ['ivtail']
            # self.ivhead = ['ivhead']
            # self.ivtype = ['ivtype']
        except Exception as inst:
            print(inst)

    def _preprocess(self):
        # put three input directional vector into a direction list
        self.directions = list()
        for ix, iy, iz in zip(self.ixd, self.iyd, self.izd):
            self.directions.append((ix, iy, iz))
        self.mean = np.mean(self.vr)
        self.variance = np.var(self.vr)

    def gamma(self):
        """
        This subroutine computes any of eight different measures of spatial
        continuity for regular spaced 3-D data.  Missing values are allowed
        and the grid need not be cubic.
        """
        # initialize the summation arrays for each direction, variogram and lag
        self._preprocess()

        head_data = list()
        tail_data = list()
        for i in range(self.ndir):
            head_data.append(list())
            tail_data.append(list())
            for j in range(self.nlag):
                head_data[i].append(list())
                tail_data[i].append(list())
        # loop over all points on the grid
        coordination = np.meshgrid(
            np.arange(self.nx), np.arange(self.ny), np.arange(self.nz))
        # loop over all points on the grid
        for ix, iy, iz in zip(
                coordination[0].flatten(),
                coordination[1].flatten(),
                coordination[2].flatten()):
            # loop over each direction
            for idir, (ixd, iyd, izd) in enumerate(self.directions):
                ix1, iy1, iz1 = ix + ixd, iy + iyd, iz + izd
                ilag = 0
                # add indexes of every eligible pair to
                # index_pairs[nlag][npairs] list
                while ix1 < self.nx and iy1 < self.ny and iz1 < self.nz:
                    if ilag < self.nlag:
                        head_value = self.vr[(ix, iy, iz)]
                        tail_value = self.vr[(ix1, iy1, iz1)]
                        # if head_value > s
                        if not np.isnan(head_value) and \
                                not np.isnan(tail_value):
                            head_data[idir][ilag].append(head_value)
                            tail_data[idir][ilag].append(tail_value)
                    else:
                        break
                    ilag = ilag + 1
                    ix1 = ix1 + ixd
                    iy1 = iy1 + iyd
                    iz1 = iz1 + izd
        # after figure out all the index pairs, use the index to get data

        self.gam = np.zeros((self.ndir, self.nlag))
        self.npair = np.zeros((self.ndir, self.nlag), dtype='>i4')
        # hm = np.zeros((self.ndir, self.nlag))
        # tm = np.zeros((self.ndir, self.nlag))
        # hv = np.zeros((self.ndir, self.nlag))
        # tv = np.zeros((self.ndir, self.nlag))

        for i, (head_lags, tail_lags) in enumerate(zip(head_data, tail_data)):
            for j, (head_lag, tail_lag) in enumerate(zip(head_lags, tail_lags)):
                self.npair[i][j] = len(head_lag)
                for hd, td in zip(head_lag, tail_lag):
                    self.gam[i][j] = self.gam[i][j] + (hd - td)**2
        self.gam /= self.npair
        self.gam *= 0.5  # * self.gam
        if self.standardize is True:
            self.gam /= self.variance

    def graph(self):
        fig, axes = plt.subplots(nrows=self.ndir, ncols=1)
        for ax, data, (idx, idy, idz) in zip(axes, self.gam, self.directions):
            real_lag_value = np.sqrt((idx * self.xsiz)**2 +
                                     (idy * self.ysiz)**2 +
                                     (idz * self.zsiz)**2)
            temp = np.ones(self.nlag) * real_lag_value
            abscissa = temp * np.arange(1, self.nlag + 1)
            ax.plot(abscissa, data, linestyle='--', marker='s',
                    fillstyle='none', color='black', markeredgecolor='blue')
#            ax.step(abscissa, data, label='pre (default)')
#            ax.scatter(abscissa, data, marker='D')
            ax.set_xlim(left=0, right=real_lag_value*(self.nlag+1))
            ax.set_ylim(bottom=0)
            ax.set_title("directional vector [{}, {}, {}]".format(
                                                              idx, idy, idz))
            ax.set_ylabel("Variogram $\gamma(h)$")
            ax.set_xlabel("Distance")
        fig.tight_layout()
        plt.draw()

if __name__ == "__main__":
    data_analysis = Gam("testData/xihuSmall_sparse_gam.par")
    data_analysis.read_data()
    data_analysis.gamma()
    data_analysis.graph()
