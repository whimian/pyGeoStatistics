# -*- coding: utf-8 -*-
"""
A 3D kriging program utilizing block search scheme

which supports simple kriging (SK), ordinary kriging (OK), or kriging with
a polynomial trend model (KT) with up to nine monomial terms.

Created on Tue Nov 22 2016
"""
from __future__ import absolute_import, division, print_function

__author__ = "yuhao"

import json
import time
from collections import namedtuple
from itertools import izip, product

import matplotlib.pyplot as plt
import numpy as np
from numba import jit
from scipy import linalg

from .super_block import SuperBlockSearcher


class Krige3d(object):
    """
    Performing 3d Kriging with super block search
    """
    def __init__(self, param_file):
        self.param_file = param_file
        self._read_params()
        self._check_params()
        self.property_name = None
        self.vr = None
        self.rotmat = None
        self.estimation = None
        self.estimation_variance = None

        self.xdb = None
        self.ydb = None
        self.zdb = None

        self._2d = False
        self.searcher = None
        self.const = None

        self._block_covariance = None
        self._unbias = None
        self.maxcov = None
        self._mdt = None

        self.resc = None

    def _read_params(self):
        with open(self.param_file) as fin:
            params = json.load(fin)
            # data file definition
            self.datafl = params['datafl']  #: 'testData/test.gslib',
            # self.idhl = None  # ????
            self.ixl = params['icolx']  #: 1,
            self.iyl = params['icoly']  #: 2,
            self.izl = params['icolz']
            self.ivrl = params['icolvr']  #: 0,
            self.iextv = params['icolsec']  # scalar, used for cross-validation
            # data limits
            self.tmin = params['tmin']  #: -1.0e21,
            self.tmax = params['tmax']  #: 1.0e21,
            # Validation Options: 0:no, 1:crossvalidation, 2:jackknife
            self.koption = params['option']
            # definition of jackknife data file
            self.jackfl = params['jackfl']
            self.ixlj = params['jicolx']  #: 1,
            self.iylj = params['jicoly']  #: 2,
            self.izlj = params['jicolz']
            self.ivrlj = params['jicolvr']  #: 0,
            self.iextvj = params['jicolsec']
            # debug and output data file
            self.idbg = params['idbg']  #: 3,
            self.dbgfl = params['dbgfl']  #: 'kt3d.dbg',
            self.outfl = params['outfl']  #: 'out.dat',
            # Grid definition
            self.nx = params['nx']  #: 50,
            self.xmn = params['xmn']  #: 0.5,
            self.xsiz = params['xsiz']  #: 1.0,
            self.ny = params['ny']  #: 50,
            self.ymn = params['ymn']  #: 0.5,
            self.ysiz = params['ysiz']  #: 1.0,
            self.nz = params['nz']  #: 50,
            self.zmn = params['zmn']  #: 0.5,
            self.zsiz = params['zsiz']  #: 1.0,
            # discretization definition
            self.nxdis = params['nxdis']  #: 1,
            self.nydis = params['nydis']  #: 1,
            self.nzdis = params['nzdis']  #: 1,
            # maximum and minimum data points used in kriging
            self.ndmin = params['ndmin']  #: ,
            self.ndmax = params['ndmax']  #: ,
            # maximum number to retain from an octant
            # (an octant search is not used if noct=0)
            self.noct = params['noct']
            # search radii
            self.radius_hmax = params['radius_hmax']  # scalar
            self.radius_hmin = params['radius_hmin']  # scalar
            self.radius_vert = params['radius_vert']  # scalar
            # search ellipsoid
            self.sang1 = params['sang1']  # scalar
            self.sang2 = params['sang2']  # scalar
            self.sang3 = params['sang3']  # scalar
            # kriging type
            self.ktype = params['ikrige']
            self.skmean = params['skmean']
            # external drift definition
            self.idrift = params['idrift']  # list of boolean
            self.itrend = params['itrend']  # boolean
            self.extfl = params['secfl']
            self.iextve = params['iseccol']  # scalar
            # Vairography definition
            self.nst = params['nst']
            self.c0 = params['c0']
            self.it = params['it']
            self.cc = params['cc']
            self.ang1 = params['ang1']
            self.ang2 = params['ang2']
            self.ang3 = params['ang3']
            self.aa_hmax = params['aa_hmax']
            self.aa_hmin = params['aa_hmin']
            self.aa_vert = params['aa_vert']

    def _check_params(self):
        # Check search radius
        if self.radius_hmax <= 0:
            raise ValueError("radius_hmax should be larger than zero.")
        # Check variograms
        if self.nst <= 0:
            raise ValueError("nst must be at least 1.")
        for vtype, a_range in zip(self.it, self.aa_hmax):
            if vtype not in np.arange(1, 6):
                raise ValueError("INVALID variogram number {}".format(vtype))
            if vtype == 4:
                if a_range < 0:
                    raise ValueError("INVALID power variogram")
                elif a_range > 2.0:
                    raise ValueError("INVALID power variogram")
        # Check data file definition
        if self.ktype == 3 and self.iextv <= 0:
            raise ValueError("Must have exteranl variable")
        if self.ixl < 0 and self.nx > 1:
            raise ValueError("WARNING: ixl=0 and nx>1 !")
        if self.iyl < 0 and self.ny > 1:
            raise ValueError("WARNING: iyl=0 and ny>1 !")
        if self.izl < 0 and self.nz > 1:
            raise ValueError("WARNING: izl=0 and nz>1 !")
        # check Trend term
        for item in self.idrift:
            if not isinstance(item, bool):
                raise ValueError("Invalid drift term {}".format(item))

    def read_data(self):
        "Read a simplified Geo-EAS formatted file."
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

    def _preprocess(self):
        """create variables needed before performing kriging"""
        # calculate dimensional constants
        krige3d_const = namedtuple('Krige3d_const',
                                   ['PMX', 'MAXNST', 'MAXDT', 'MAXSB',
                                    'MAXDIS', 'MAXSAM', 'UNEST'])
        maxsbx = 1
        if self.nx > 1:
            maxsbx = int(self.nx/2)
            if maxsbx > 50:
                maxsbx = 50
        maxsby = 1
        if self.ny > 1:
            maxsby = int(self.ny/2)
            if maxsby > 50:
                maxsby = 50
        maxsbz = 1
        if self.nz > 1:
            maxsbz = int(self.nz/2)
            if maxsbz > 50:
                maxsbz = 50
        self.const = krige3d_const(
            PMX=999,
            MAXNST=4,
            MAXDT=9,
            MAXSB=(maxsbx, maxsby, maxsbz),
            MAXDIS=self.nxdis * self.nydis * self.nzdis,
            MAXSAM=self.ndmax + 1,
            UNEST=np.nan
            )
        # Calculate needed programing variables from input parameters
        self.radsqd = self.radius_hmax * self.radius_hmax
        self.sanis1 = self.radius_hmin / self.radius_hmax
        self.sanis2 = self.radius_vert / self.radius_hmax

        self.anis1 = np.array(self.aa_hmin) / \
                     np.maximum(self.aa_hmax, np.finfo(float).eps)
        self.anis2 = np.array(self.aa_vert) / \
                     np.maximum(self.aa_hmax, np.finfo(float).eps)

        # set up for validation, if cross-validation, set jackfl as datafl
        if self.koption == 1:
            self.jackfl = self.datafl
            # self.idhlj = self.idhl
            self.ixlj = self.ixl
            self.iylj = self.iyl
            self.izlj = self.izl
            self.ivrlj = self.ivrl
            self.iextvj = self.iextv

    def _set_rotation(self):
        """
        Set up rotation matrix for both anisotropy and searching.
        with self.rotmat being an array of 3*3 rotation matrix, the last matrix
        in the array are the searching matrix
        """
        ang1 = np.append(self.ang1, self.sang1)
        ang2 = np.append(self.ang2, self.sang2)
        ang3 = np.append(self.ang3, self.sang3)
        anis1 = np.append(self.anis1, self.sanis1)
        anis2 = np.append(self.anis2, self.sanis2)

        self.rotmat = np.full((ang1.shape[0], 3, 3), np.nan)
        def convert_ang1(ang):
            if ang <= 0 and ang < 270:
                alpha = np.deg2rad(90 - ang)
            else:
                alpha = np.deg2rad(450 - ang)
            return alpha
        v_convert = np.vectorize(convert_ang1)

        alpha = v_convert(ang1)
        beta = np.deg2rad(-ang2)
        theta = np.deg2rad(ang3)

        sina = np.sin(alpha)
        sinb = np.sin(beta)
        sint = np.sin(theta)
        cosa = np.cos(alpha)
        cosb = np.cos(beta)
        cost = np.cos(theta)

        afac1 = 1.0 / np.maximum(anis1, np.finfo(float).eps)
        afac2 = 1.0 / np.maximum(anis2, np.finfo(float).eps)
        self.rotmat[:, 0, 0] = cosb * cosa
        self.rotmat[:, 0, 1] = cosb * sina
        self.rotmat[:, 0, 2] = -sinb
        self.rotmat[:, 1, 0] = afac1 * (-cost * sina + sint * sinb * cosa)
        self.rotmat[:, 1, 1] = afac1 * (cost * cosa + sint * sinb * sina)
        self.rotmat[:, 1, 2] = afac1 * (sint * cosb)
        self.rotmat[:, 2, 0] = afac2 * (sint * sina + cost * sinb * cosa)
        self.rotmat[:, 2, 1] = afac2 * (-sint * cosa + cost * sinb * sina)
        self.rotmat[:, 2, 2] = afac2 * (cost * cosb)

    def kt3d(self):
        self._preprocess()
        # Set up the rotation/anisotropy matrices needed for variogram
        # and searching
        self._set_rotation()
        # compute maximum covariance for the rescaling factor:
        self._max_covariance()
        # compute rescaling factor:
        self._rescaling()
        # Set up for super block searching:
        print("Setting up Super Block Search...")
        self._create_searcher()
        # Set up discretization points per block
        self._block_discretization()
        # Find unbias value
        self.unbias = self.maxcov

        # mean values of the drift function
        self.bv = np.zeros((9,))
        self.bv[0] = np.mean(self.xdb) * self.resc
        self.bv[1] = np.mean(self.ydb) * self.resc
        self.bv[2] = np.mean(self.zdb) * self.resc
        self.bv[3] = np.mean(self.xdb * self.xdb) * self.resc
        self.bv[4] = np.mean(self.ydb * self.ydb) * self.resc
        self.bv[5] = np.mean(self.zdb * self.zdb) * self.resc
        self.bv[6] = np.mean(self.xdb * self.ydb) * self.resc
        self.bv[7] = np.mean(self.xdb * self.zdb) * self.resc
        self.bv[8] = np.mean(self.ydb * self.zdb) * self.resc
        # report on progress from time to time:
        nd = self.vr.shape[0]
        if self.koption == 0:  # kriging
            nxy = self.nx * self.ny
            nxyz = self.nx * self.ny * self.nz
            nloop = nxyz
            # irepo = max(1, min((nxyz/10), 10000))
        else:  # Validation
            nloop = 10000000
            irepo = max(1, min(nd/10, 10000))
        print("Start working on the kriging...")
        # time
        t1 = time.time()
        ts = 0
        percent_od = 0

        self.estimation = np.full((nloop,), np.nan)
        self.estimation_variance = np.full((nloop,), np.nan)
        # MAIN LOOP OVER ALL THE BLOCKS IN THE GRID:
        for index in xrange(nloop):
            if self.koption == 0:
                self.iz = index // nxy
                self.iy = (index - self.iz * nxy) // self.nx
                self.ix = index - self.iz * nxy - self.iy * self.nx
                xloc = self.xmn + self.ix * self.xsiz
                yloc = self.ymn + self.iy * self.ysiz
                zloc = self.zmn + self.iz * self.zsiz
            else:  # crossvalidation or jackknife
                # read(ljack,*,err=96,end=2) (var(i),i=1,nvarij)
                var = list()
                # ddh = 0.0
                xloc = self.xmn
                yloc = self.ymn
                zloc = self.zmn
                true = self.const.UNEST
                # secj = self.const.UNEST
                # if self.idhlj > 0:
                #     ddh = var[idhlj]
                if self.ixlj > 0:
                    xloc = var[self.ixlj]
                if self.iylj > 0:
                    yloc = var[self.iylj]
                if self.izlj > 0:
                    zloc = var[self.izlj]
                if self.ivrlj > 0:
                    true = var[self.ivrlj]
                if self.iextvj > 0:
                    self.extest = var[self.iextvj]
                if true < self.tmin or true >= self.tmax:
                    true = self.const.UNEST
            # read in the external drift variable if needed.
            if self.ktype == 2 or self.ktype == 3:  # non-SK or KED
                if self.koption == 0:
                    # read(lext,*) (var(i),i=1,iextve)
                    var = list()
                    self.extest = var[self.iextve]  # colocated external value
                if self.extest < self.tmin or self.extest >= self.tmax:
                    self.estimation[index] = self.const.UNEST
                    self.estimation_variance[index] = self.const.UNEST
                    continue
                # rescalling factor for external drift variable
                self.resce = self.maxcov / max(self.extest, 0.0001)
            # Search for proximity data
            # ts_1 = time.time()
            self.searcher.search(xloc, yloc, zloc)
            # ts += time.time() - ts_1
            # load nearest data in xa, ya, za, vra, vea
            xa = list()
            ya = list()
            za = list()
            vra = list()
            vea = list()  # colocated external drift value
            na = 0
            for i in xrange(self.searcher.nclose):
                ind = self.searcher.close_samples[i]
                accept = True
                if self.koption != 0 and \
                    abs(self.vr['x'][ind] - xloc) + \
                    abs(self.vr['y'][ind] - yloc) + \
                    abs(self.vr['z'][ind] - zloc) < np.finfo(float).eps:
                    accept = False
                # if self.koption != 0 and \
                #     abs(dh[ind] - ddh) < np.finfo(float).eps:
                #     accept = False
                if accept:
                    if na < self.ndmax:
                        xa.append(self.vr['x'][ind] - xloc)
                        ya.append(self.vr['y'][ind] - yloc)
                        za.append(self.vr['z'][ind] - zloc)
                        vra.append(self.vr[self.property_name[0]][ind])
                        if self.ktype == 3:  # KED
                            vea.append(self.vr[self.property_name[1]])
                        na += 1
            # check number of samples found
            if na < self.ndmin:
                self.estimation[index] = self.const.UNEST
                self.estimation_variance[index] = self.const.UNEST
                print("not enough data.")
                continue
            # Test if there are enough samples to estimate all drift terms:
            if na >= 1 and na <= self.mdt:
                # if firon:
                #     firon = False
                self.estimation[index] = self.const.UNEST
                self.estimation_variance[index] = self.const.UNEST
                print("not enough data to estimate all drift terms")
                continue
            xa = np.array(xa)
            ya = np.array(ya)
            za = np.array(za)
            vra = np.array(vra)
            vea = np.array(vea)
            # Enough data, proceed with estimation
            if na <= 1:
                est, estv = self._one_sample(xa, ya, za, vra)
                self.estimation[index] = est
                self.estimation_variance[index] = estv
            else:
                est, estv = self._many_samples(xa, ya, za, vra, vea)
                self.estimation[index] = est
                self.estimation_variance[index] = estv
            # print working percentage
            percent = np.round(index/nloop*100, decimals=0)
            dtime = time.time() - t1
            if percent != percent_od and percent % 10 == 0:
                print("{}% ".format(percent) +\
                  "."*20 + "{}s elapsed.".format(np.round(dtime, decimals=3)))
            percent_od = percent
        print("Kriging Finished.")
        # print("Time used for searching: {}s".format(ts))

    def _rescaling(self):
        if self.radsqd < 1:
            self.resc = 2 * self.radius_hmax / max(self.maxcov, 0.0001)
        else:
            self.resc = (4 * self.radsqd) / max(self.maxcov, 0.0001)
        if self.resc <= 0:
            raise ValueError("rescaling value is wrong, {}".format(self.resc))
        self.resc = 1 / self.resc

    def _one_sample(self, xa, ya, za, vra):
        """
        If only one sample, perform SK or OK

        Parameters
        ----------
        xa, ya, za, vra: 1-D ndarray
        """
        # Left hand side
        left = cova3(
            (xa[0], ya[0], za[0]), (xa[0], ya[0], za[0]),
            self.rotmat, self.maxcov, self.nst, self.it, self.cc, self.aa_hmax)
        #Right hand side
        if self.ndb <= 1:
            right = cova3(
                (xa[0], ya[0], za[0]), (self.xdb[0], self.ydb[0], self.zdb[0]),
                self.rotmat, self.maxcov,
                self.nst, self.it, self.cc, self.aa_hmax)
        else:
            right = 0
            for i in xrange(self.ndb):
                # cov = self._cova3((xa[0], ya[0], za[0]),
                #                   (self.xdb[i], self.ydb[i], self.zdb[i]))
                cov = cova3(
                    (xa[0], ya[0], za[0]),
                    (self.xdb[i], self.ydb[i], self.zdb[i]),
                    self.rotmat, self.maxcov,
                    self.nst, self.it, self.cc, self.aa_hmax)
                right += cov
                dx = xa[0] - self.xdb[i]
                dy = ya[0] - self.ydb[i]
                dz = za[0] - self.zdb[i]
                if dx*dx + dy*dy + dz*dz < np.finfo(float).eps:
                    right -= self.c0
            right /= self.ndb
        if self.ktype == 2: # non-sationary SK
            self.skmean = self.extest
        if self.ktype == 0 or self.ktype == 2: # SK
            # solve for lambda
            s = right / self.block_covariance
            est = s * vra[0] + (1-s) * self.skmean
            estv = self.block_covariance - s * right
        else:  # OK
            est = vra[0]
            estv = self.block_covariance - 2 * right + left
        return est, estv

    def _many_samples(self, xa, ya, za, vra, vea):
        """
        More than one sample

        Parameters
        ----------
        xa, ya, za, vra: 1-D ndarray
        """
        na = len(xa)
        neq = self.mdt + na  # number of equations
        # Left Hand Side
        # first = False
        left = left_side(
            xa, ya, za, neq, self.unbias, self.rotmat, self.maxcov, self.nst,
            self.it, self.cc, self.aa_hmax)

        # Right Hand Side
        right = right_side(
            xa, ya, za, self.xdb, self.ydb, self.zdb, neq, self.unbias,
            self.rotmat, self.maxcov, self.nst, self.it, self.cc, self.aa_hmax,
            self.c0)

        # Add the additional unbiasedness constraints:
        im = na + 1
        # First drift term (linear in 'x')
        if self.idrift[0] is True:
            im += 1
            left[im, :im] = xa * self.resc
            left[:im, im] = xa * self.resc
            left[im, im:] = 0
            left[im:, im] = 0
            # right.append(self.bv[0])
            right[im] = self.bv[0]
        # Second drift term (linear in 'y'):
        if self.idrift[1] is True:
            im += 1
            left[im, :im] = ya * self.resc
            left[:im, im] = ya * self.resc
            left[im, im:] = 0
            left[im:, im] = 0
            # right.append(self.bv[1])
            right[im] = self.bv[1]
        # Third drift term (linear in 'z')
        if self.idrift[2] is True:
            im += 1
            left[im, :im] = za * self.resc
            left[:im, im] = za * self.resc
            left[im, im:] = 0
            left[im:, im] = 0
            # right.append(self.bv[2])
            right[im] = self.bv[2]
        # Fourth drift term (quadratic in 'x')
        if self.idrift[3] is True:
            im += 1
            left[im, :im] = xa * xa * self.resc
            left[:im, im] = xa * xa * self.resc
            left[im, im:] = 0
            left[im:, im] = 0
            # right.append(self.bv[3])
            right[im] = self.bv[3]
        # Fifth drift term (quadratic in 'y')
        if self.idrift[4] is True:
            im += 1
            left[im, :im] = ya * ya * self.resc
            left[:im, im] = ya * ya * self.resc
            left[im, im:] = 0
            left[im:, im] = 0
            # right.append(self.bv[4])
            right[im] = self.bv[4]
        # Sixth drift term (quadratic in 'z')
        if self.idrift[5] is True:
            im += 1
            left[im, :im] = za * za * self.resc
            left[:im, im] = za * za * self.resc
            left[im, im:] = 0
            left[im:, im] = 0
            # right.append(self.bv[5])
            right[im] = self.bv[5]
        # Seventh drift term (quadratic in 'xy')
        if self.idrift[6] is True:
            im += 1
            left[im, :im] = xa * ya * self.resc
            left[:im, im] = xa * ya * self.resc
            left[im, im:] = 0
            left[im:, im] = 0
            # right.append(self.bv[6])
            right[im] = self.bv[6]
        # Eighth drift term (quadratic in 'xz')
        if self.idrift[7] is True:
            im += 1
            left[im, :im] = xa * za * self.resc
            left[:im, im] = xa * za * self.resc
            left[im, im:] = 0
            left[im:, im] = 0
            # right.append(self.bv[7])
            right[im] = self.bv[7]
        # Ninth drift term (quadratic in 'yz')
        if self.idrift[8] is True:
            im += 1
            left[im, :im] = ya * za * self.resc
            left[:im, im] = ya * za * self.resc
            left[im, im:] = 0
            left[im:, im] = 0
            # right.append(self.bv[8])
            right[im] = self.bv[8]
        # External drift term
        if self.ktype == 3:  # KED
            im += 1
            left[im, :im] = vea * self.resce
            left[:im, im] = vea * self.resce
            left[im, im:] = 0
            left[im:, im] = 0
            # right.append(self.extest * self.resce)
            right[im] = self.extest * self.resce

        # if estimating the trend then reset the right terms all to 0.0
        if self.itrend == True:
            right = np.full((neq,), np.nan)

        # Solve the kriging system
        s = None
        try:
            s = linalg.solve(left, right)
        except linalg.LinAlgError as inst:
            print("Warning kt3d: Singular matrix " + \
                    "{}, {}, {}".format(self.ix, self.iy, self.iz))
            return np.nan, np.nan
        # Estimate and estimation variance
        if self.ktype == 2:  # non-stationary SK
            self.skmean = self.extest
        # Variance
        estv = self.block_covariance - np.sum(s * right)
        # Estimate
        est = 0
        if self.ktype == 0:  # SK
            est = np.sum(s[:na] * (vra[:na] - self.skmean))
        elif self.ktype == 2:  # non-stationary SK
            est = np.sum(s[:na] * (vra - vea))
        else:  # OK, KED
            est = np.sum(s[:na] * vra)

        if self.ktype == 0 or self.ktype == 2:  # SK or non-stationary SK
            est += self.skmean

        return est, estv

    @property
    def block_covariance(self):
        "return average covariance within block"
        if self._block_covariance is None:
            if self.ndb <= 1:  # point kriging
                self._block_covariance = self.unbias
            else:
                cov = list()
                for x1, y1, z1 in izip(self.xdb, self.ydb, self.zdb):
                    for x2, y2, z2 in izip(self.xdb, self.ydb, self.zdb):
                        # cov.append(self._cova3((x1, y1, z1), (x2, y2, z2)))
                        cov.append(cova3(
                            (x1, y1, z1), (x2, y2, z2),
                            self.rotmat, self.maxcov, self.nst,
                            self.it, self.cc, self.aa_hmax))
                cov = np.array(cov).reshape((self.ndb, self.ndb))
                cov[np.diag_indices_from(cov)] -= self.c0
                self._block_covariance = np.mean(cov)
        return self._block_covariance

    @property
    def mdt(self):
        """
        The number of drift terms.

        If an external drift is being considered then there is one more
        drift term other than those less than nine drift terms

        And if SK is being considered,
        then we will set all drift terms off and mdt to 0.

        The property is used to determine how many extra rows/columns there
        are in the kriging matrix.
        """
        if self._mdt is None:
            self._mdt = 1
            for i in xrange(9):
                if self.ktype == 0 or self.ktype == 2:
                    self.idrift[i] = 0
                self._mdt += self.idrift[i]
            if self.ktype == 3:  #KED
                self._mdt += 1
            elif self.ktype == 0:  # SK
                self._mdt = 0
            elif self.ktype == 2:  # non-stationary SK
                self._mdt = 0
        return self._mdt

    def _create_searcher(self):
        "Help create and initialize the searcher object"
        self.searcher = SuperBlockSearcher()
        # initialize required atrributes
        # grid definition
        self.searcher.nx = self.nx
        self.searcher.xmn = self.xmn
        self.searcher.xsiz = self.xsiz
        self.searcher.ny = self.ny
        self.searcher.ymn = self.ymn
        self.searcher.ysiz = self.ysiz
        self.searcher.nz = self.nz
        self.searcher.zmn = self.zmn
        self.searcher.zsiz = self.zsiz
        # data
        self.searcher.vr = self.vr
        self.searcher.MAXSB = self.const.MAXSB
        # rotation matrix
        self.searcher.rotmat = self.rotmat[-1]
        self.searcher.radsqd = self.radsqd
        # octant search
        self.searcher.noct = self.noct
        # Setup
        self.searcher.setup()
        self.searcher.pickup()
        # sort data according to superblock number
        self.vr = self.vr[self.searcher.sort_index]

    def _block_discretization(self):
        self.nxdis = 1 if self.nxdis < 1 else self.nxdis
        self.nydis = 1 if self.nydis < 1 else self.nydis
        self.nzdis = 1 if self.nzdis < 1 else self.nzdis
        self.ndb = self.nxdis * self.nydis * self.nzdis
        if self.ndb > self.const.MAXDIS:
            raise ValueError("Too many discretization points")
        xdis = self.xsiz / max(self.nxdis, 1)
        ydis = self.ysiz / max(self.nydis, 1)
        zdis = self.zsiz / max(self.nzdis, 1)
        self.xdb = np.arange(0, self.nxdis, 1) * xdis + \
                   (-0.5 * self.xsiz + 0.5 * xdis)
        self.ydb = np.arange(0, self.nydis, 1) * ydis + \
                   (-0.5 * self.ysiz + 0.5 * ydis)
        self.zdb = np.arange(0, self.nzdis, 1) * zdis + \
                   (-0.5 * self.zsiz + 0.5 * zdis)


    def _max_covariance(self):
        '''
        Calculate the maximum covariance value (used for zero distances and
        for power model covariance):
        '''
        self.maxcov = self.c0
        for ist in xrange(self.nst):
            if self.it[ist] == 4:
                self.maxcov += self.const.PMX
            else:
                self.maxcov += self.cc[ist]

    def view2d(self):
        "View 2D data using matplotlib"
        if self._2d is False:
            print("3D data, use view3d() instead.")
        else:
            fig, ax = plt.subplots()
            im = ax.imshow(self.estimation.reshape(self.ny, self.nx),
                           interpolation='nearest',
                           origin='lower',
                           extent=[self.xmn,
                                   self.xmn + (self.nx - 1)*self.xsiz,
                                   self.ymn,
                                   self.ymn + (self.ny - 1)*self.ysiz],
                           cmap='jet')
            ax.set_xlabel("X (m)")
            ax.set_ylabel("Y (m)")
            ax.set_title("Estimation")
            ax.set_aspect('equal')
            fig.colorbar(im)
            fig.show()

    def view3d(self):
        "View 3D data using mayavi"
        pass

@jit(nopython=True)
def left_side(xa, ya, za, neq, unbias, rotmat, maxcov, nst, it, cc, aa_hmax):
    na = len(xa)
    left = np.full((neq, neq), np.nan)
    # fill the kriging matrix:
    # for i, j in product(xrange(na), xrange(na)):
    for i in xrange(na):
        for j in xrange(na):
            if np.isnan(left[j, i]):
                # left[i, j] = self._cova3((xa[i], ya[i], za[i]),
                #                          (xa[j], ya[j], za[j]))
                left[i, j] = cova3(
                    (xa[i], ya[i], za[i]), (xa[j], ya[j], za[j]),
                    rotmat, maxcov, nst, it, cc, aa_hmax)
            else:
                left[i, j] = left[j, i]
    if neq > na:  # fill for OK
        left[na, :na] = unbias
        left[:na, na] = unbias
        left[na, na] = 0
    return left

@jit(nopython=True)
def right_side(xa, ya, za, xdb, ydb, zdb, neq, unbias, rotmat, maxcov,
               nst, it, cc, aa_hmax, c0):
    na = len(xa)
    ndb = len(xdb)
    right = np.full((neq,), np.nan)
    for i in xrange(na):
        if ndb <= 1:
            cb = cova3(
                (xa[i], ya[i], za[i]),
                (xdb[0], ydb[0], zdb[0]),
                rotmat, maxcov, nst,
                it, cc, aa_hmax)
        else:
            cb = 0
            for j in xrange(ndb):
                cov = cova3(
                    (xa[i], ya[i], za[i]),
                    (xdb[j], ydb[j], zdb[j]),
                    rotmat, maxcov,
                    nst, it, cc, aa_hmax)
                cb += cov
                dx = xa[i] - xdb[j]
                dy = ya[i] - ydb[j]
                dz = za[i] - zdb[j]
                if dx*dx + dy*dy + dz*dz < 7./3-4./3-1:
                    cb -= c0
        cb /= ndb
        right[i] = cb
    if neq > na:
        right[na:] = unbias
    return right

@jit(nopython=True)
def sqdist(point1, point2, rotmat):
    """
    This routine calculates the anisotropic distance between two points
    given the coordinates of each point and a definition of the
    anisotropy.

    This method only consider a single anisotropy senario.

    Parameters
    ----------
    point1 : tuple
        Coordinates of first point (x1,y1,z1)
    point2 : tuple
        Coordinates of second point (x2,y2,z2)
    rotmat : 3*3 ndarray
        matrix of rotation for this structure

    Returns
    -------
    sqdist : scalar
        The squared distance accounting for the anisotropy
        and the rotation of coordinates (if any).
    """
    dx = point1[0] - point2[0]
    dy = point1[1] - point2[1]
    dz = point1[2] - point2[2]
    sqdist = 0.0
    for i in xrange(3):
        cont = rotmat[i, 0] * dx + \
                rotmat[i, 1] * dy + \
                rotmat[i, 2] * dz
        sqdist += cont * cont
    return sqdist

@jit(nopython=True)
def cova3(point1, point2, rotmat, maxcov, nst, it, cc, aa_hmax):
    """
    Parameters
    ----------
    point1, point2: tuple of 3
        coordinates of two points

    Returns
    -------
    cova: scalar
        covariance between (x1,y1,z1) and (x2,y2,z2)
    """
    # check for 'zero' distance, return maxcov if so:
    hsqd = sqdist(point1, point2, rotmat[0])
    if hsqd < 7./3-4./3-1:
        cova = maxcov
        return cova
    # loop over all structures
    cova = 0
    for ist in xrange(nst):
        if ist != 0:
            hsqd = sqdist(point1, point2, rotmat[ist])
        h = np.sqrt(hsqd)
        if it[ist] == 1:  # Spherical
            hr = h / aa_hmax[ist]
            if hr < 1:
                cova += cc[ist] * (1 - hr * (1.5 - 0.5 * hr * hr))
        elif it[ist] == 2:  # Exponential
            cova += cc[ist] * np.exp(-3.0 * h / aa_hmax[ist])
        elif it[ist] == 3:  # Gaussian
            cova += cc[ist] * \
                np.exp(-3.0 * (h / aa_hmax[ist]) * (h/aa_hmax[ist]))
        elif it[ist] == 4:  # Power
            cova += maxcov - cc[ist] * (h**(aa_hmax[ist]))
        elif it[ist] == 5:  # Hole Effect
            cova += cc[ist] * np.cos(h / aa_hmax[ist] * np.pi)
    return cova
