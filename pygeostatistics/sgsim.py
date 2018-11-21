# -*- coding: utf-8 -*-
"""
Sequential Gaussian Simulation

Created on Sun Dec 4 2016
"""
from __future__ import absolute_import, division, print_function

__author__ = "yuhao"

# import json
import os
import time
from collections import namedtuple
from itertools import product

import yaml
import matplotlib.pyplot as plt
import numpy as np
from numba import jit
from scipy import interpolate, linalg

from pygeostatistics.normal_score_transform import NormalScoreTransform, gauinv
from pygeostatistics.super_block import SuperBlockSearcher
from pygeostatistics.yaml_patch import loader_patched


class Sgsim(object):
    """
    Performing 3d Kriging with super block search
    """
    def __init__(self, param_file):
        self.param_file = param_file
        self._read_params()
        self._check_params()
        self.property_name = None
        self.var = None
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

        self.tmpfl = None
        self.icolvr = None
        self.icolwt = None

    def _read_params(self):
        with open(self.param_file, 'r') as fin:
            params = yaml.load(fin, Loader=loader_patched())
            # data file definition
            self.datafl = params['datafl']  #: 'testData/test.gslib',
            # self.idhl = None  # ????
            self.ixl = params['icolx']  #: 1,
            self.iyl = params['icoly']  #: 2,
            self.izl = params['icolz']
            self.ivrl = params['icolvr']  #: 0,
            self.iwt = params['icolwt']  # declustering weights
            self.isecvr = params['icolsec']  # for external drift if used
            # data limits
            self.tmin = params['tmin']  #: -1.0e21,
            self.tmax = params['tmax']  #: 1.0e21,
            self.itrans = params['itrans']  # boolean
            # output file for transformation table if transformation is needed
            self.transfl = params['transfl']
            self.ismooth = params['ismooth']  # boolean
            # file with values used for transformation to normal scores
            self.smthfl = params['smthfl']
            self.isvr = params['icolsvr']
            self.iswt = params['icolswt']
            # allowable data values used for backtransform
            self.zmin = params['zmin']
            self.zmax = params['zmax']
            # lower and upper tail model specification for backtransform
            self.ltail = params['ltail']
            self.ltpar = params['ltpar']
            self.utail = params['utail']
            self.utpar = params['utpar']
            # debug and output data file
            self.idbg = params['idbg']  #: 3,
            self.dbgfl = params['dbgfl']  #: 'kt3d.dbg',
            self.outfl = params['outfl']  #: 'out.dat',
            self.nsim = params['nsim']  # number of simulation
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
            self.ixv = params['seed']  # random seed
            # maximum and minimum data points used in kriging
            self.ndmin = params['ndmin']  #: ,
            self.ndmax = params['ndmax']  #: ,
            self.nodmax = params['nodmax']  # previously simulated nodes to use
            self.sstrat = params['sstrat']  # search strategy
            self.mults = params['multgrid']  # boolean
            self.nmult = params['nmult']  # scalar
            self.noct = params['noct']  # maximum number to retain from an octant
            # search radii
            self.radius_hmax = params['radius_hmax']  # scalar
            self.radius_hmin = params['radius_hmin']  # scalar
            self.radius_vert = params['radius_vert']  # scalar
            # search anisotropy angles
            self.sang1 = params['sang1']  # scalar
            self.sang2 = params['sang2']  # scalar
            self.sang3 = params['sang3']  # scalar
            # size of maximum covariance lookup table size
            self.mxctx = params['mxctx']
            self.mxcty = params['mxcty']
            self.mxctz = params['mxctz']
            # kriging type
            self.ktype = params['ikrige']
            # self.skmean = params['skmean']
            self.colocorr = params['rho']  # correlation coefficient for COCOK
            self.varred = params['varred']  # variance reduction factor for COCOK
            self.lvmfl = params['secfl']
            self.icollvm = params['icollvm']
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
                raise ValueError("Power model is not allowed.")
        # check for sill
        self.sill = self.c0 + np.sum(self.cc)
        if self.sill != 1:
            print("Warning: Sill is not 1. ({})".format(self.sill))
        # check tail interpolation parameters
        if self.ltail not in [1, 2]:
            raise ValueError("invalid lower tail option")
        if self.utail not in [1, 2, 4]:
            raise ValueError("invalid upper tail option")
        if self.utail == 4 and self.utpar < 0:
            raise ValueError("invalid power for hyperbolic tail")
        if self.ltail == 2 and self.ltpar < 0:
            raise ValueError("invalid power for power model")
        if self.utail == 2 and self.utpar < 0:
            raise ValueError("invalid power for power model")

    def _preprocess(self):
        # calculate dimensional constants
        # sgsim_const = namedtuple('sgsim_const',
        #                            ['PMX', 'MAXNST', 'MAXDT', 'MAXSB',
        #                             'MAXDIS', 'MAXSAM', 'UNEST'])
        sgsim_const = namedtuple('sgsim_const',
                                 ['PMX', 'MAXNST', 'MAXDT', 'MAXSB',
                                  'MAXSAM', 'UNEST'])
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
        self.const = sgsim_const(
            PMX=999,
            MAXNST=4,
            MAXDT=9,
            MAXSB=(maxsbx, maxsby, maxsbz),
#            MAXDIS=self.nxdis * self.nydis * self.nzdis,
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

    def _read_data(self):
        if self.itrans is True:
            # Establish the reference histogram for simulation
            print("Setting up transformation table...")
            # decide which file to use for establishing the transformation table
            if self.ismooth is True:
                self.tmpfl = self.smthfl
                self.icolvr = self.isvr
                self.icolwt = self.iswt
            else:  # data histogram, possibly with declustering weights is used for transformation.
                self.tmpfl = self.datafl
                self.icolvr = self.ivrl
                self.icolwt = self.iwt
            # read in file for transformation table
            var = self._read_file(self.tmpfl)
            # check
            column_names = var.dtype.names
            vrtr = var[column_names[self.icolvr]]
            if self.icolwt < 0:
                vrgtr = np.ones_like(vrtr)
            else:
                vrgtr = var[column_names[self.icolwt]]
            # remove illegal data
            mask1 = vrtr >= self.tmin
            mask2 = vrtr < self.tmax
            mask3 = vrgtr > 0
            mask = mask1 * mask2 * mask3
            vrtr = vrtr[mask]
            vrgtr = vrgtr[mask]
            # sort data by value
            sort_index = np.argsort(vrtr)
            vrtr = vrtr[sort_index]
            vrgtr = vrgtr[sort_index]
            # perform transform
            self.nst_primary = NormalScoreTransform(
                vrtr, vrgtr, self.zmin, self.zmax, self.ltail, self.ltpar,
                self.utail, self.utpar)
            self.nst_primary.create_transform_func()

            # write transformation table to file

        # read in hard data for simulation
        print("Reading Data...")
        self.var = self._read_file(self.datafl)
        column_names = self.var.dtype.names
        self.x = self.var['x']
        self.y = self.var['y']
        self.z = self.var['z']
        if self.ivrl >= 0:
            self.vr = self.var[column_names[self.ivrl]]
        if self.iwt >= 0:
            self.wt = self.var[column_names[self.iwt]]
        if self.isecvr >= 0:
            self.sec = self.var[column_names[self.isecvr]]
        # normal score transform
        if self.itrans is True:
            vrg = self.nst_primary.transform(self.vr)
            self.vr = vrg  # give nst value to data.

        # read in secondary attribute model if needed:
        if self.ktype >= 2:
            print("Reading secondary data...")
            var = self._read_file(self.lvmfl)
            column_names = var.dtype.names
            self.lvm = var[column_names[self.icollvm]]  # locally varying mean
            self.sim = np.arange(self.nx * self.ny * self.nz)
        # transform the secondary variable for a local mean?
        # if trans is True and self.ktype == 2 and self.itrans is True:
        if self.ktype == 2 and self.itrans is True:
            # self.lvm = interp_func(self.lvm)
            self.lvm = self.nst_primary.transform(self.lvm)
        # if need to work with data residuals? (locally varying mean)
        if self.ktype == 2:
            x_index, y_index, z_index = self._getindx(self.x, self.y, self.z)
            index = x_index + y_index * self.nx + z_index * self.nx * self.ny
            self.sec = self.lvm[index]
            # calculation of residual moved to krige subroutine.
        # if need to get an external drift attribute for the data?
        if self.ktype == 3:
            x_index, y_index, z_index = self._getindx(self.x, self.y, self.z)
            index = x_index + y_index * self.nx + z_index * self.nx * self.ny
            self.sec = self.lvm[index]
        if self.ktype == 4:
            print("transforming secondary data with variance reduction"+\
                  " of {}".format(self.varred))
            # sort according to lvm
            sort_index = np.argsort(self.lvm)
            self.lvm = self.lvm[sort_index]
            self.sim = self.sim[sort_index]
            nxyz = self.nx * self.ny * self.nz
            oldcp = np.arange(nxyz) * (1/nxyz)
            cp = oldcp + (1/nxyz)
            w = (cp + oldcp)*0.5
            self.lvm = gauinv(w) * self.varred
            # sort back according to sim
            sort_index = np.argsort(self.sim)
            self.lvm = self.lvm[sort_index]
            self.sim = self.sim[sort_index]

    def _getindx(self, xloc, yloc, zloc):
        """
        determine which superblock are the given point or list of points in

        Parameters
        ----------
        xloc, yloc, zloc: scalar or 1-D ndarray
        """
        x_block = np.arange(self.xmn - 0.5 * self.xsiz,
                            self.xmn + (self.nx + 1) * self.xsiz + 1,
                            self.xsiz)
        x_index = np.searchsorted(x_block, xloc) - 1

        y_block = np.arange(self.ymn - 0.5 * self.ysiz,
                            self.ymn + (self.ny + 1) * self.ysiz + 1,
                            self.ysiz)
        y_index = np.searchsorted(y_block, yloc) - 1

        z_block = np.arange(self.zmn - 0.5 * self.zsiz,
                            self.zmn + (self.nz + 1) * self.zsiz + 1,
                            self.zsiz)
        z_index = np.searchsorted(z_block, zloc) - 1

        return (x_index, y_index, z_index)

    def _read_file(self, filename):
        "Read a simplified Geo-EAS formatted file."
        data_list = None
        with open(filename, 'r') as fin:
            data_list = fin.readlines()
        name = data_list[0].strip()
        ncols = int(data_list[1].strip())
        column_name = [item.strip() for item in data_list[2: ncols+2]]
        if 'z' not in column_name:
            self.izl = -1
            column_name.append('z')
            data_list = [tuple(item.strip().split() + ['0'])
                         for item in data_list[ncols+2:]]
        else:
            data_list = [tuple(item.strip().split())
                         for item in data_list[ncols+2:]]
        data_dtype = np.dtype({
            'names': column_name,
            'formats': ['f8'] * len(column_name)})
        return np.array(data_list, dtype=data_dtype)

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
        self.searcher.vr = self.var
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
        self.var = self.var[self.searcher.sort_index]

    def _ctable(self):
        """
        Create covariance lookup table

        The idea is to establish a 3-D network that contains the covariance
        value for a range of grid node offsets that should be at as large
        as twice the search radius in each direction.  The reason it has to
        be twice as large as the search radius is because we want to use it
        to compute the data covariance matrix as well as the data-point
        covariance matrix.
        Secondly, we want to establish a search for nearby nodes that
        in order of closeness as defined by the variogram.

        OUTPUT: self.covtab
        """
        tiny = 1.0e-10
        # definition of covariance table dimension (radius)
        self.nctx = min(((self.mxctx-1)//2), (self.nx-1))
        self.ncty = min(((self.mxcty-1)//2), (self.ny-1))
        self.nctz = min(((self.mxctz-1)//2), (self.nz-1))

        self.covtab = np.full(
            (2 * self.nctx + 1, 2 * self.ncty + 1, 2 * self.nctz + 1), np.nan)
        tmp = list()
        order = list()
        for i, j, k in product(range(-self.nctx, self.nctx+1),
                               range(-self.ncty, self.ncty+1),
                               range(-self.nctz, self.nctz+1)):
            xx = i * self.xsiz
            ic = self.nctx + i
            yy = j * self.ysiz
            jc = self.ncty + j
            zz = k * self.zsiz
            kc = self.nctz + k
            # cova = self._cova3((0, 0, 0), (xx, yy, zz))
            cova = cova3(
                (0, 0, 0), (xx, yy, zz), self.rotmat,
                self.maxcov, self.nst, self.it, self.cc, self.aa_hmax)
            self.covtab[ic, jc, kc] = cova
            # hsqd = self._sqdist((0, 0, 0), (xx, yy, zz), self.rotmat[-1])
            hsqd = sqdist((0, 0, 0), (xx, yy, zz), self.rotmat[-1])
            if hsqd <= self.radsqd:
                tmp.append(-(cova - tiny * hsqd))
                order.append((ic, jc, kc))
        # sort according to variogram distance tmp
        sort_index = np.argsort(tmp)
        order = np.array(order,
                         dtype=np.dtype({
                             'names': ['i', 'j', 'k'],
                             'formats': ['i4'] * 3}))
        order = order[sort_index]
        self.ixnode = order['i']
        self.iynode = order['j']
        self.iznode = order['k']

    def _srchnd(self, ix, iy, iz, sim):
        """
        Search for nearby Simulated Grid nodes in a Spiral fashion.

        The idea is to spiral away from the node being simulated and note all
        the nearby nodes that have been simulated.

        INPUT VARIABLES:
          ix,iy,iz:        index of the point currently being simulated
          sim:             the realization so far
          nodmax:          the maximum number of nodes that we want
          nlookup:          the number of nodes in the look up table
          ixnode, iynode, iznode:    the relative indices of those nodes.
          xmn, ymn, zmn:       the origin of the global grid netwrok
          xsiz, ysiz, zsiz:      the spacing of the grid nodes.

        OUTPUT VARIABLES:
          ncnode          the number of close nodes
          icnode()        the number in the look up table
          cnodex, cnodey, cnodez: list
            the location of the nodes
          cnodev()        the values at the nodes
        """
        nx = self.nx
        ny = self.ny
        nz = self.nz
        nxy = nx * ny
        nxyz = nx * ny * nz
        UNEST = -99
        nlookup = self.ixnode.shape[0]

        self.ncnode, self.icnode, self.cnodex, \
        self.cnodey, self.cnodez, self.cnodev = \
            search(ix, iy, iz,
                   self.xmn, self.ymn, self.zmn,
                   self.xsiz, self.ysiz, self.zsiz,
                   nx, ny, nz, nxy,
                   self.ixnode, self.iynode, self.iznode,
                   self.nctx, self.ncty, self.nctz,
                   sim,
                   self.noct,
                   nlookup, self.nodmax)

    def _krige(self, ix, iy, iz, xx, yy, zz, lktype):
        """
        INPUT VARIABLES:
        ix, iy, iz:
            index of the point currently being simulated
        xx, yy, zz:
            location of the point currently being simulated
        OUTPUT VARIABLES:
        cmean:
            kriged estimate
        cstdev:
            kriged standard deviation
        """
        # Size of the kriging system
        nxy = self.nx * self.ny
        nx = self.nx
        # first = False
        na = self.searcher.nclose + self.ncnode  # sample points + simulated points
        if lktype == 0:  # simple kriging
            neq = na
        elif lktype == 1:  # ordinary kriging
            neq = na + 1
        elif lktype == 2:  # simple kriging with a locally varying mean
            neq = na
        elif lktype == 3:  # kriging with external drift
            neq = na + 2
        elif neq == 4:
            neq = na + 1
        else:
            print("unsupported kriging type")
        if lktype >= 3:
            ind = ix + (iy-1)*nx + (iz-1)*nxy
            if self.lvm[ind] <= -6.0 or self.lvm[ind] >= 6.0:
                lktype = 0
        # Set up kriging matrices:
        # initialize coordination, data, secondary data lists
        vra = np.full((na,), np.nan)
        left, right = krige_matrix(
            lktype, ix, iy, iz,
            xx, yy, zz,
            vra,
            na, neq, nx, nxy,
            self.searcher.nclose, self.searcher.close_samples,
            self.x, self.y, self.z, self.vr,
            self.cnodex, self.cnodey, self.cnodez, self.cnodev, self.icnode,
            self.ixnode, self.iynode, self.iznode,
            self.nctx, self.ncty, self.nctz,
            self.mxctx, self.mxcty, self.mxctz,
            self.rotmat, self.maxcov, self.nst, self.it, self.cc, self.aa_hmax,
            self.covtab)
        # vrea = np.full((na,), np.nan)  # secondary data

        # Solve kriging system
        s = None
        try:
            s = linalg.solve(left, right)
        except linalg.LinAlgError as inst:
            print("Warning kt3d: Singular matrix " + \
                    "{}, {}, {}".format(ix, iy, iz))
            return np.nan, np.nan
        except ValueError as e:
            print(e.message)
            return np.nan, np.nan

        # Compute the estimate and kriging variance.
        cmean = 0
        cestdev = self.maxcov
#        sumwts = 0

        cmean = np.sum(s * vra)
        cstdev = cestdev - np.sum(s * right)
#        sumwts = np.sum(s)
        if lktype == 1: # OK
            cstdev -= 2*s[-1]

        return cmean, cstdev

    def _max_covariance(self):
        '''
        Calculate the maximum covariance value (used for zero distances and
        for power model covariance):
        '''
        self.maxcov = self.c0
        for ist in range(self.nst):
            if self.it[ist] == 4:
                self.maxcov += self.const.PMX
            else:
                self.maxcov += self.cc[ist]

    def simulate(self):
        """
        Simulate
        """
        print("Preparing...")
        self._preprocess()
        self._read_data()
        self._set_rotation()
        self._max_covariance()
        if self.sstrat == 0:  # data and points are searched separately
            print("Setting up super block search...")
            self._create_searcher()
        print("Setting up Covariance Lookup Table and Spiral Search...")
        self._ctable()
        nxyz = self.nx * self.ny * self.nz
        nxy = self.nx * self.ny
        nx = self.nx
        # doing `nsim` times of simulation
        for isim in range(self.nsim):
            if isim > 0 and self.ktype == 4:
                pass
            # random path for this realization
            sim = np.random.rand(nxyz,)  # random numbers in sim[]
            order = np.arange(nxyz)

            # multiple grid search, 4 is arbitrary
            if self.mults is True:
                for imult in range(self.nmult):
                    nnz = max(1, self.nz / (imult * 4))
                    nny = max(1, self.ny / (imult * 4))
                    nnx = max(1, self.nx / (imult * 4))
                    for iz, iy, ix in product(range(nnz), range(nny),
                                              range(nnx)):
                        jz = iz * imult * 4
                        jy = iy * imult * 4
                        jx = ix * imult * 4
                        index = jx + jy*self.nx + jz*nxy
                        sim[index] -= imult
                        # larger grid nodes are pushed forward by giving them
                        # negative values
            # sort according to random number order
            sort_index = np.argsort(sim)
            order = order[sort_index]  # simulation path

            print("Working on realization {}".format(isim + 1))
            # GRID NODE data assignment
            TINY = 0.0001
            UNEST = -99
            EPSILON = np.finfo(float).eps
            sim = np.full_like(sim, UNEST)
            ix, iy, iz = self._getindx(self.x, self.y, self.z)
            ind = ix + iy * nx + iz * nxy
            xx = self.xmn + ix * self.xsiz
            yy = self.ymn + iy * self.ysiz
            zz = self.zmn + iz * self.zsiz
            test = np.abs(xx - self.x) + np.abs(yy - self.y) + \
                   np.abs(zz - self.z)
            nd = self.x.shape[0]

            for idx in range(nd):
                # assign this data to the node unless there is a closer data.
                if self.sstrat == 1:
                    if sim[ind[idx]] >= 0:
                        id2 = int(sim[ind[idx]]+0.5)
                        test2 = abs(xx[idx] - self.x[id2]) + \
                                abs(yy[idx] - self.y[id2]) + \
                                abs(zz[idx] - self.z[id2])
                        if test[idx] <= test2:
                            sim[ind[idx]] = idx
                            print("data values are both assigned to the same node.")
                        else:
                            sim[ind[idx]] = idx

                if self.sstrat == 0 and test[idx] <= TINY:
                    # too close, has to use hard data on this node
                    # assign a flag so that this node does not get simulated.
                    sim[ind[idx]] = 10 * UNEST

            # enter data values into the simulated grid:
            for i in range(nxyz):
                idx = sim[i]
                if idx > 0:
                    sim[i] = self.vr[idx]

            t1 = time.time()
            ts = 0
            percent_od = 0
            # MAIN LOOP
            for igrid in range(nxyz):
                index = order[igrid]
                if sim[index] > UNEST + EPSILON or sim[index] < UNEST * 2:
                    continue
                iz = index // nxy
                iy = (index - iz * nxy) // nx
                ix = index - iz * nxy - iy * nx
                xcoor = self.xmn + ix * self.xsiz
                ycoor = self.ymn + iy * self.ysiz
                zcoor = self.zmn + iz * self.zsiz

                # simulate this (ix, iy, iz)
                if self.sstrat == 0:
                    self.searcher.search(xcoor, ycoor, zcoor)
                    if self.searcher.nclose < self.ndmin:
                        continue
                    if self.searcher.nclose > self.ndmax:
                        self.searcher.nclose = self.ndmax
                self._srchnd(ix, iy, iz, sim) # spiral search on the covariance lookup table
                # calculate the conditional mean and standard deviation
                # this will be done with kriging if there are data,
                # otherwise, the glocal mean and standard deviation will be used.
                if self.ktype == 2:
                    gmean = self.lvm[index]
                else:
                    gmean = 0
                if self.searcher.nclose + self.ncnode < 1:
                    cmean = gmean
                    cstdev = 1
                else:
                    # kriging
                    lktype = self.ktype
                    if self.ktype == 1 and self.searcher.nclose + self.ncnode < 4:
                        lktype = 0
                    cmean, cstdev = self._krige(ix, iy, iz, xcoor, ycoor, zcoor, lktype)
                # draw a random number and assign a value to this node
                p = np.random.rand()
                xp = gauinv(p)
                sim[index] = xp * cstdev + cmean
                percent = np.round(igrid/nxyz*100, decimals=0)
                dtime = time.time() - t1
                if percent != percent_od and percent % 10 == 0:
                    print("{}% ".format(percent) +\
                      "."*20 + "{}s elapsed.".format(np.round(dtime, decimals=3)))
                percent_od = percent
            # END MAIN LOOP

            # do we need to reassign the data to the grid nodes?
            if self.sstrat == 0:
                for idx in range(nd):
                    if test[idx] <= TINY:
                        sim[ind[idx]] = self.vr[idx]
            # back transform each value and write results:
            if self.itrans is True:
                mask = sim > UNEST + EPSILON
                temp = sim[mask]
                sim[mask] = self.nst_primary.back_func(temp)
                mask = sim < self.zmin
                sim[mask] = self.zmin
                mask = sim > self.zmax
                sim[mask] = self.zmax

            # save this simulation in file
            parent_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.path.pardir)
            sim_dir = os.path.join(parent_dir, 'simulations')
            np.save(os.path.join(sim_dir, 'sim_{}'.format(isim+1)), sim)

    def view2d(self):
        "View 2D data using matplotlib"
        sim_dir = os.path.join(os.path.join(os.path.dirname(
            os.path.abspath(__file__)), os.path.pardir), 'simulations')
        for file_name in os.listdir("simulations/"):
            estimation = np.load("simulations/{}".format(file_name))
            fn, ext = os.path.splitext(file_name)
            fig, ax = plt.subplots()
            im = ax.imshow(estimation.reshape(self.ny, self.nx),
                           interpolation='nearest',
                           origin='lower',
                           extent=[self.xmn,
                                   self.xmn + (self.nx - 1)*self.xsiz,
                                   self.ymn,
                                   self.ymn + (self.ny - 1)*self.ysiz],
                           cmap='jet')
            ax.set(xlabel="X (m)", ylabel="Y (m)",
                   title="Realization {}".format(fn), aspect='equal')
            fig.colorbar(im)
            fig.show()

            fig2, ax2 = plt.subplots()
            hist, bin_edges = np.histogram(estimation, bins=20)
            ax2.set_title("Histogram {}".format(fn))
            ax2.bar(bin_edges[:-1], hist, width=bin_edges[1]-bin_edges[0],
                    color='red', alpha=0.5)
            fig2.show()

    def average(self):
        result = np.zeros((self.nx * self.ny, ))
        for file_name in os.listdir("simulations/"):
            estimation = np.load("simulations/{}".format(file_name))
            result += estimation
        result /= 100
        fn, ext = os.path.splitext(file_name)
        fig, ax = plt.subplots()
        im = ax.imshow(result.reshape(self.ny, self.nx),
                       interpolation='nearest',
                       origin='lower',
                       extent=[self.xmn,
                               self.xmn + (self.nx - 1)*self.xsiz,
                               self.ymn,
                               self.ymn + (self.ny - 1)*self.ysiz],
                       cmap='jet')
        ax.set(xlabel="X (m)", ylabel="Y (m)",
               title="Average", aspect='equal')
        fig.colorbar(im)
        fig.show()

        fig2, ax2 = plt.subplots()
        hist, bin_edges = np.histogram(result, bins=20)
        ax2.set_title("Histogram")
        ax2.bar(bin_edges[:-1], hist, width=bin_edges[1]-bin_edges[0],
                color='red', alpha=0.5)
        fig2.show()

    def view3d(self):
        "View 3D data using mayavi"
        pass

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
    for i in range(3):
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
    for ist in range(nst):
        if ist != 0:
            # hsqd = self._sqdist(point1, point2, self.rotmat[ist])
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

@jit(nopython=True)
def search(ix, iy, iz,
           xmn, ymn, zmn,
           xsiz, ysiz, zsiz,
           nx, ny, nz, nxy,
           ixnode, iynode, iznode,
           nctx, ncty, nctz,
           sim,
           noct,
           nlookup, nodmax):
    ncnode = 0
    icnode = [] # the number in the lookup table
    cnodex = []
    cnodey = []
    cnodez = []
    cnodev = []
    ninoct = np.zeros((8,))

    for il in range(0, nlookup):
        if ncnode == nodmax:
            return
        # calculate location of point in simulation grid
        i = ix + ixnode[il] - nctx - 1
        j = iy + iynode[il] - ncty - 1
        k = iz + iznode[il] - nctz - 1
        # covaraince lookup table is large so we need to make sure
        # the point is within the simulation grid.
        if i < 0 or j < 0 or k < 0:
            continue
        if i > nx or j > ny or k > nz:
            continue
        ind = i + j * nx + k * nxy
        if not np.isnan(sim[ind]):
            # check the number of data already taken from this octant:
            if noct > 0:
                idx = ix - i
                idy = iy - j
                idz = iz - k
                if idz > 0:
                    iq = 3
                    if idx <= 0 and idy > 0:
                        iq = 0
                    if idx > 0 and idy >= 0:
                        iq = 1
                    if idx < 0 and idy <= 0:
                        iq = 2
                else:
                    iq = 7
                    if idx <= 0 and idy > 0:
                        iq = 4
                    if idx > 0 and idy >= 0:
                        iq = 5
                    if idx < 0 and idy <= 0:
                        iq = 6
                ninoct[iq] += 1
                if ninoct[iq] > noct:
                    continue
            ncnode += 1
            icnode.append(il)
            cnodex.append(xmn + i * xsiz)
            cnodey.append(ymn + j * ysiz)
            cnodez.append(zmn + k * zsiz)
            cnodev.append(sim[ind])
    icnode_a = np.array(icnode)
    cnodex_a = np.array(cnodex)
    cnodey_a = np.array(cnodey)
    cnodez_a = np.array(cnodez)
    cnodev_a = np.array(cnodev)
    return ncnode, icnode_a, cnodex_a, cnodey_a, cnodez_a, cnodev_a

@jit(nopython=True)
def krige_matrix(lktype, ix, iy, iz,
                 xx, yy, zz,
                 vra, na, neq, nx, nxy,
                 nclose, close_samples,
                 x, y, z, vr,
                 cnodex, cnodey, cnodez, cnodev, icnode,
                 ixnode, iynode, iznode,
                 nctx, ncty, nctz,
                 mxctx, mxcty, mxctz,
                 rotmat, maxcov, nst, it, cc, aa_hmax,
                 covtab):
    left = np.full((neq, neq), np.nan)
    right = np.full((neq,), np.nan)
    for j in range(na):
        # sort out the actual location of point "j"
        # if j <= self.nclose:
        if j < nclose:
            index = close_samples[j]
            x1 = x[index]
            y1 = y[index]
            z1 = z[index]
            vra[j] = vr[index]
            # if lktype >= 2:
            #     vrea[j] = self.sec[index]
            # if lktype == 2:
            #     vra[j] -= vrea[j]
        # if it is a previously simulated node (keep index for table look-up)
        else:
            index = j - nclose
            x1 = cnodex[index]
            y1 = cnodey[index]
            z1 = cnodez[index]
            vra[j] = cnodev[index]
            ind = icnode[index]
            ix1 = ix + ixnode[ind] - nctx - 1
            iy1 = iy + iynode[ind] - ncty - 1
            iz1 = iz + iznode[ind] - nctz - 1
            index = ix1 + (iy1-1)*nx + (iz1-1)*nxy
            # if lktype >= 2:
            #     vrea[j] = lvm[index]
            # if lktype == 2:
            #     vra[j] -= vrea[j]
        for i in range(j+1):
            if i < nclose:
                index = close_samples[i]
                x2 = x[index]
                y2 = y[index]
                z2 = z[index]
            else:  # it's a previously simulated node (keep index for table look-up)
                index = i - nclose
                x2 = cnodex[index]
                y2 = cnodex[index]
                z2 = cnodez[index]
                ind = icnode[index]
                ix2 = ix + ixnode[ind] - nctx - 1
                iy2 = iy + iynode[ind] - ncty - 1
                iz2 = iz + iznode[ind] - nctz -1

            if j <= nclose or i <= nclose:
                # calculate covariance value
                # left[i, j] = self._cova3((x1, y1, z1), (x2, y2, z2))
                left[i, j] = cova3(
                    (x1, y1, z1), (x2, y2, z2), rotmat,
                    maxcov, nst, it, cc, aa_hmax)
            else:
                # try to use lookup table if distance is in range.
                # ii = self.nctx + 1 + (ix1 - ix2)
                # jj = self.ncty + 1 + (iy1 - iy2)
                # kk = self.nctz + 1 + (iz1 - iz2)
                ii = nctx + (ix1 - ix2)
                jj = ncty + (iy1 - iy2)
                kk = nctz + (iz1 - iz2)
                if ii < 0 or ii > mxctx or \
                    jj < 0 or jj > mxcty or \
                    kk < 0 or kk > mxctz:
                    # left[i, j] = self._cova3((x1, y1, z1), (x2, y2, z2))
                    left[i, j] = cova3(
                        (x1, y1, z1), (x2, y2, z2), rotmat,
                        maxcov, nst, it, cc,
                        aa_hmax)
                else:
                    left[i, j] = covtab[ii, jj, kk]
        # calculate right side value:
        if j < nclose:
            # right.append(_cova3((xx, yy, zz), (x1, y1, z1)))
            right[j] = cova3(
                (xx, yy, zz), (x1, y1, z1), rotmat,
                maxcov, nst, it, cc, aa_hmax)
        else:
            # try to use lookup table if distance is in range.
            # ii = nctx + 1 + (ix1 - ix2)
            # jj = ncty + 1 + (iy1 - iy2)
            # kk = nctz + 1 + (iz1 - iz2)
            ii = nctx + (ix1 - ix2)
            jj = ncty + (iy1 - iy2)
            kk = nctz + (iz1 - iz2)
            if ii < 0 or ii > mxctx or \
                jj < 0 or jj > mxcty or \
                kk < 0 or kk > mxctz:
                # right.append(_cova3((xx, yy, zz), (x1, y1, z1)))
                right[i] = cova3(
                    (xx, yy, zz), (x1, y1, z1), rotmat,
                    maxcov, nst, it, cc,
                    aa_hmax)

            else:
                right[i] = covtab[ii, jj, kk]
    # fill void elements of left matrix:
    # for i, j in product(range(na), range(na)):
    for i in range(na):
        for j in range(na):
            if np.isnan(left[i, j]):
                left[i, j] = left[j, i]

    # Addition of OK constraint:
    if lktype == 1 or lktype == 3:
        right[neq:] = 1
        left[neq, :neq] = 1  # self.unbias
        left[:neq, neq] = 1  # self.unbias
        left[neq, neq] = 0
    # Addition of the External Drift Constraint:
    # not implemented
    # Addition of Collocated Cosimulation Constraint:
    # not implemented
    return left, right
