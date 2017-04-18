# -*- coding: utf-8 -*-
"""
Compute Variogram on irregularly spaced data

Created on Sun Nov 06 2016
"""
from __future__ import division, print_function
import json
import numpy as np
import matplotlib.pyplot as plt
# import pandas as pd
#from scipy.spatial.distance import pdist, squareform

__author__ = "yuhao"

class Gamv():
    def __init__(self, param_file):
        self.param_file = param_file
        self._read_params()
        self._check_params()
        self.vr = None
        self.gam = None
        self.npair = None
        self.distance = None
        self.mean = None
        self.variance = None
        self.hm = None
        self.tm = None
        self.property_name = None
        self.lag_interval = None

    def _read_params(self):
        with open(self.param_file) as fin:
            params = json.load(fin)
            self.datafl = params['datafl']
            self.icolx = params['icolx']  #: 1,
            self.icoly = params['icoly']  #: 2,
            self.icolz = params['icolz']  #: 0,
            self.nvar = params['nvar']  # : 2,
            self.ivar = params['ivar']  # : [3, 4],
            self.tmin = params['tmin']  # : -1.0e21,
            self.tmax = params['tmax']  # : 1.0e21,
            self.outfl = params['outfl']  # : 'out.dat',
            self.nlag = params['nlag']  # : 10,
            self.xlag = params['xlag']  # : 5.0,
            self.xltol = params['xltol']  # : 3.0,
            self.ndir = params['ndir']  # : 3,
            self.azm = params['azm']  # : [0.0, 0.0, 90.],
            self.atol = params['atol']  # : [90.0, 22.5, 22.5],
            self.bandwh = params['bandwh']  # : [50.0, 25.0, 25.0],
            self.dip = params['dip']  # : [0.0, 0.0, 0.0],
            self.dtol = params['dtol']  # : [90.0, 22.5, 22.5],
            self.bandwd = params['bandwd']  # : [50.0, 25.0, 25.0],
            self.standardize = params['standardize']  # : False,
            self.nvarg = params['nvarg']  # : 3,
            self.ivtail = params['ivtail']  # : [1, 1, 2],
            self.ivhead = params['ivhead']  # : [1, 1, 2],
            self.ivtype = params['ivtype']  # : [1, 3, 1]

    def _check_params(self):
        # # check lag tolerance
        # if self.xltol <= 0:
        #     self.xltol = 0.5 * self.xlag
        # # check azimuth tolerance
        # for i, item in enumerate(self.atol):
        #     if item <= 0:
        #         self.atol[i] = 45.0
        # # check dip tolerance
        # for i, item in enumerate(self.dtol):
        #     if item <= 0:
        #         self.dtol[i] = 45.0
        if self.ndir != len(self.azm):
            raise ValueError('number of directions does not match provided \
                              azimuth.')

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

    def _preprocess(self):
        self.mean = np.mean(self.vr[self.property_name[0]])
        self.variance = np.var(self.vr[self.property_name[0]])
        # check lag tolerance
        if self.xltol <= 0:
            self.xltol = 0.5 * self.xlag
        # check azimuth tolerance
        # for i, item in enumerate(self.atol):
        #     if item <= 0:
        #         self.atol[i] = 45.0
        self.atol = [45.0 if item <= 0 else item for item in self.atol]
        # check dip tolerance
        # for i, item in enumerate(self.dtol):
        #     if item <= 0:
        #         self.dtol[i] = 45.0
        self.dtol = [45.0 if item <= 0 else item for item in self.dtol]

    def gamv(self):
        self._check_params()
        self._preprocess()
        self.gam = np.zeros((self.ndir, self.nlag+2))
        self.npair = np.zeros((self.ndir, self.nlag+2), dtype='>i4')
        self.distance = np.zeros((self.ndir, self.nlag+2))
        self.tm = np.zeros((self.ndir, self.nlag+2))
        self.hm = np.zeros((self.ndir, self.nlag+2))


#        tail_data = list()
#        head_data = list()
#        for i in xrange(self.ndir):
#            tail_data.append(list())
#            head_data.append(list())
#            for j in xrange(self.nlag+2):
#                tail_data[i].append(list())
#                head_data[i].append(list())
        # calculate the value interval of each lag
        # method 1
#        temp = np.arange(0, self.nlag+2, 1.)
#        temp -= 1
#        temp[0] = 0
#        temp *= self.xlag
#        length, = temp.shape
#        temp.reshape((length, 1))
#        upper = temp + self.xltol
#        upper[0][0] = 0
#        lower = temp - self.xltol
#        lower[0][0] = 0
#        self.lag_interval = np.concatenate((lower, upper), axis=1)
        # method 2
        self.lag_interval = np.ones((self.nlag+2, 2))
        self.lag_interval[:2, :] = 0
        self.lag_interval *= self.xlag
        scale = np.arange(0, self.nlag+2, 1) - 1
        self.lag_interval[:, 0] *= scale
        self.lag_interval[:, 1] *= scale
        self.lag_interval[:, 0] -= self.xltol
        self.lag_interval[0, 0] = 0.0
        self.lag_interval[:, 1] += self.xltol
        self.lag_interval[0, 1] = 0.0
        # The mathematical azimuth is measured counterclockwise from EW and
        # not clockwise from NS as the conventional azimuth is:
        # name it azmuth instead of azimuth
        azmuth = np.deg2rad(90.0 - np.array(self.azm))
        uvxazm = np.cos(azmuth)
        uvyazm = np.sin(azmuth)
        csatol = np.cos(np.deg2rad(self.atol))
#        for i, az in enumerate(self.atol):
#            if az == 90:
#                csatol[i] = 0
        # The declination is measured positive from vertical (up) rather than
        # negative down from horizontal:
        declin = np.deg2rad(90.0 - np.array(self.dip))
        uvzdec = np.cos(declin)
        uvhdec = np.sin(declin)
        csdtol = np.cos(np.deg2rad(self.dtol))
#        for i, di in enumerate(self.dtol):
#            if di == 90:
#                csdtol[i] = 0
        # square of maxmium distance
        dismxs = ((self.nlag + 0.5 - np.finfo('float').eps) * self.xlag)**2
        num_of_data = int(self.vr.shape[0])
        ijlist = list()
        for i in xrange(num_of_data - 1):
            for j in xrange(i+1, num_of_data):
                ijlist.append((i, j))
        for ijtuple in ijlist:
            i, j = ijtuple
            # i_coor = list(self.vr[i])[:-1]
            # j_coor = list(self.vr[j])[:-1]
            # calculate the lag corresponding to the current pair
            # dx = j_coor[0] - i_coor[0]
            # dy = j_coor[1] - i_coor[1]
            # dz = j_coor[2] - i_coor[2]
            dx = self.vr['x'][j] - self.vr['x'][i]
            dy = self.vr['y'][j] - self.vr['y'][i]
            dz = self.vr['z'][j] - self.vr['z'][i]
            # square of lag h
            hs = dx**2 + dy**2 + dz**2
            if hs > dismxs:
                # print("skip pair {},{} for maxdistance".format(i, j))
                continue  # skip to next pair
            h = np.sqrt(hs)
            # determine which lag
            lag_num = list()  # could be in two lags
            for k, lg_in in enumerate(self.lag_interval):
                if h >= lg_in[0] and h <= lg_in[1]:
                    lag_num.append(k)
            if len(lag_num) == 0:
                # print("skip pair {},{} for no lag".format(i, j))
                continue  # skip if cannot find which lag

            for idir in xrange(self.ndir):
                omni = self.atol[idir] >= 90.0
                # check for an acceptable azimuth angle:
                dxy = np.sqrt(max(dx**2 + dy**2, 0.0))
                if dxy < np.finfo('float').eps:
                    dcazm = 1.0
                else:
                    dcazm = (dx*uvxazm[idir] + dy*uvyazm[idir])/dxy
                if np.abs(dcazm) < csatol[idir]:
                    # print("skip pair {},{} for az".format(i, j))
                    continue
                # check for the horizontal bandwidth criteria (maximium
                # deviation perpendicular to the specified direction azimuth):
                band = uvxazm[idir]*dy - uvyazm[idir]*dx
                if np.abs(band) > self.bandwh[idir] and not omni:
                    # print("skip pair {},{} for az bwd".format(i, j))
                    continue
                # check for an acceptable dip angle:
                if dcazm < 0:
                    dxy = -dxy
                if lag_num[0] == 0:
                    dcdec = 0
                else:
                    dcdec = (dxy*uvhdec[idir] + dz*uvzdec[idir])/h
                    if np.abs(dcdec) < csdtol[idir]:
                        # print("skip pair {},{} for dip".format(i, j))
                        continue
                # check the vertical bandwidth criteria (maximium deviation
                # perpendicular to the specified dip direction):
                band = uvhdec[idir]*dz - uvzdec[idir]*dxy
                if np.abs(band) > self.bandwd[idir]  and not omni:
                    # print("skip pair {},{} for dip bwd".format(i, j))
                    continue
                # check whether or not an omi-directional variogram is being
                # computed:
                # omni = False
                # if self.atol[idir] >= 90.0:
                #     omni = True
                omni = self.atol[idir] >= 90.0
                # then this pair is acceptable, proceed to compute variogram
                # sort out which is tail and head:
                if dcazm >= 0 and dcdec <= 0:
                    # vrh = list(self.vr[i])[-1]
                    # vrt = list(self.vr[j])[-1]
                    vrh = self.vr[self.property_name[0]][i]
                    vrt = self.vr[self.property_name[0]][j]
                    if omni:
                        # vrtpr = list(self.vr[i])[-1]
                        # vrhpr = list(self.vr[j])[-1]
                        vrtpr = self.vr[self.property_name[0]][i]
                        vrhpr = self.vr[self.property_name[0]][j]
                else:
                    # vrh = list(self.vr[j])[-1]
                    # vrt = list(self.vr[i])[-1]
                    vrh = self.vr[self.property_name[0]][j]
                    vrt = self.vr[self.property_name[0]][i]
                    if omni:
                        # vrtpr = list(self.vr[j])[-1]
                        # vrhpr = list(self.vr[i])[-1]
                        vrtpr = self.vr[self.property_name[0]][j]
                        vrhpr = self.vr[self.property_name[0]][i]
                # reject this pair on the basis of missing values:
                if vrt < self.tmin or vrh < self.tmin or\
                        vrt > self.tmax or vrh > self.tmax:
                    continue
                # COMPUTE THE APPRORIATE "VARIOGRAM" MEASURE
                #         ***Semivariogram***
                for ilag in lag_num:
                    self.npair[idir][ilag] += 1
                    self.distance[idir][ilag] += h
                    self.tm[idir][ilag] += vrt
                    self.hm[idir][ilag] += vrh
                    self.gam[idir][ilag] += (vrh - vrt)**2
                    if omni:
                        if vrtpr >= self.tmin or vrhpr >= self.tmin or\
                                vrtpr < self.tmax or vrhpr < self.tmax:
                            self.npair[idir][ilag] += 1
                            self.distance[idir][ilag] += h
                            self.tm[idir][ilag] += vrtpr
                            self.hm[idir][ilag] += vrhpr
                            self.gam[idir][ilag] += (vrhpr - vrtpr)**2
        self.gam /= self.npair
        if self.standardize is True:
            self.gam /= self.variance
        self.gam /= 2

        self.distance /= self.npair
        self.tm /= self.npair
        self.hm /= self.npair

    def graph(self):
        abscissa = np.arange(0, self.nlag+1, 1.0)
        abscissa *= self.xlag
        fig, axes = plt.subplots(nrows=self.ndir, ncols=1)
        if isinstance(axes, list):
            for i, ax in enumerate(axes):
                ordinate = np.insert(self.gam[i][2:], 0, None)
                ax.scatter(abscissa, ordinate)
                ax.set_title(r"Azimuth: {}$^\circ$({}$^\circ$), ".format(
                    self.azm[i], self.atol[i]) +
                             r"Dip: {}$^\circ$({}$^\circ$)".format(
                                 self.dip[i], self.dtol[i]))
                ax.set_ylim(bottom=0)
                ax.set_xlim(left=0)
                ax.grid()
        else:
            ordinate = np.insert(self.gam[0][2:], 0, None)
            axes.scatter(abscissa, ordinate)
            axes.set_title(r"Azimuth: {}$^\circ$({}$^\circ$), ".format(
                self.azm[0], self.atol[0]) +
                           r"Dip: {}$^\circ$({}$^\circ$)".format(
                               self.dip[0], self.dtol[0]))
            axes.set_ylim(bottom=0)
            axes.set_xlim(left=0)
            axes.grid()
        fig.tight_layout()
        # plt.draw()
        return fig, axes


if __name__ == "__main__":
    data_analysis = Gamv("testData/xihuSmall_sparse_gamv.par")
    data_analysis.read_data()
    data_analysis.gamv()
    data_analysis.graph()
