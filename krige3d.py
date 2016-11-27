# -*- coding: utf-8 -*-
"""
A 3D kriging program utilizing block search scheme

which supports simple kriging (SK), ordinary kriging (OK), or kriging with
a polynomial trend model (KT) with up to nine monomial terms.

Created on Tue Nov 22 2016
"""
from __future__ import division, print_function
import json
from itertools import izip
import time
from collections import namedtuple
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
from .super_block import SuperBlockSearcher

__author__ = "yuhao"

class Krige3d(object):
    def __init__(self, param_file):
        Krige3d_const = namedtuple('Krige3d_const', ['PMX', 'MAXNST', 'MAXDT'])
        self.const = Krige3d_const(999, 4, 9)
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
        self.searcher = None

        self.MAXSB = None

    def _read_params(self):
        with open(self.param_file) as fin:
            params = json.load(fin)
            self.datafl = params['datafl']  #: 'testData/test.gslib',

            self.idhl = None  # ????

            self.ixl = params['icolx']  #: 1,
            self.iyl = params['icoly']  #: 2,
            self.izl = params['icolz']
            self.ivrl = params['icolvr']  #: 0,
            self.iextv = params['icolsec']

            self.tmin = params['tmin']  #: -1.0e21,
            self.tmax = params['tmax']  #: 1.0e21,

            self.koption = params['option']
            self.jackfl = params['jackfl']
            self.ixlj = params['jicolx']  #: 1,
            self.iylj = params['jicoly']  #: 2,
            self.izlj = params['jicoloz']
            self.ivrlj = params['jicolvr']  #: 0,
            self.iextvj = params['jicolsec']

            self.idbg = params['idbg']  #: 3,
            self.dbgfl = params['dbgfl']  #: 'kb2d.dbg',
            self.outfl = params['outfl']  #: 'out.dat',

            self.nx = params['nx']  #: 50,
            self.xmn = params['xmn']  #: 0.5,
            self.xsiz = params['xsiz']  #: 1.0,

            self.ny = params['ny']  #: 50,
            self.ymn = params['ymn']  #: 0.5,
            self.ysiz = params['ysiz']  #: 1.0,

            self.nz = params['nz']  #: 50,
            self.zmn = params['zmn']  #: 0.5,
            self.zsiz = params['zsiz']  #: 1.0,

            self.nxdis = params['nxdis']  #: 1,
            self.nydis = params['nydis']  #: 1,
            self.nzdis = params['nzdis']  #: 1,

            self.ndmin = params['ndmin']  #: ,
            self.ndmax = params['ndmax']  #: ,
            # the maximum number to retain from an octant (an octant search
            # is not used if noct=0
            self.noct = params['noct']
            # search radii
            self.radius_hmax = params['radius_hmax']  # scalar
            self.radius_hmin = params['radius_hmin']  # scalar
            self.radius_vert = params['radius_vert']  # scalar
            # search ellipsoid
            self.sang1 = params['sang1']  # scalar
            self.sang2 = params['sang2']  # scalar
            self.sang3 = params['sang3']  # scalar

            self.ktype = params['ikrige']
            self.skmean = params['skmean']

            self.idrift = params['idrift']
            self.itrend = params['itrend']
            self.extfl = params['secfl']
            self.iextve = params['iseccol']

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
        if self.radius_hmax <= 0:
            raise ValueError("radius_hmax should be larger than zero.")
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
        if self.ktype == 3 and self.iextv <= 0:
            raise ValueError("Must have exteranl variable")
        if self.ixl <= 0 and self.nx > 1:
            raise ValueError("WARNING: ixl=0 and nx>1 !")
        if self.iyl <= 0 and self.ny > 1:
            raise ValueError("WARNING: iyl=0 and ny>1 !")
        if self.izl <= 0 and self.nz > 1:
            raise ValueError("WARNING: izl=0 and nz>1 !")
        # check MAXSB
        if not isinstance(self.MAXSB, list):
            raise ValueError("SuperBlock Definition should be set with a list")
        elif len(self.MAXSB) != 3:
            raise ValueError("SuperBlock Search Definition need 3 elements")

        # check idrift
        for item in self.idrift:
            if item < 0 or item > 1:
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
        """
        Calculate needed programing variables from input parameters
        """
        self.radsqd = self.radius_hmax * self.radius_hmax
        self.sanis1 = self.radius_hmin / self.radius_hmax
        self.sanis2 = self.radius_vert / self.radius_hmax

        self.anis1 = np.array(self.aa_hmin) / \
                     np.maximum(self.aa_hmax, np.finfo(float).eps)
        self.anis2 = np.array(self.aa_vert) / \
                     np.maximum(self.aa_hmax, np.finfo(float).eps)

    def setrot(self):
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
        self.setrot()
        # compute maximum covariance for the rescaling factor:
        covmax = self.c0[0]
        for ist in xrange(self.nst):
            if self.it[ist] == 4:
                covmax += self.const.PMX
            else:
                covmax += self.cc[ist]
        # compute rescaling factor:
        if self.radsqd < 1:
            resc = 2 * self.radius_hmax / max(covmax, 0.0001)
        else:
            resc = (4 * self.radsqd) / max(covmax, 0.0001)
        if resc <= 0:
            raise ValueError("rescaling value is wrong, {}".format(resc))
        resc = 1 / resc
        # Set up for super block searching:
        self._create_searcher()
        # compute the number of drift terms, if an external drift is being
        # considered then it is one more drift term, if SK is being considered
        # then we will set all drift terms off and mdt to 0.
        mdt = 1
        for i in xrange(9):
            if self.ktype == 0 or self.ktype == 2:
                self.idrift[i] = 0
            mdt += self.idrift[i]
        if self.ktype == 3:  #KED
            mdt += 1
        elif self.ktype == 0:  # SK
            mdt = 0
        elif self.ktype == 2:  # UK
            mdt = 0
        # Set up discretization points per block
        self._block_discretization()

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
        self.searcher.MAXSB = self.MAXSB
        # rotation matrix
        self.searcher.rotmat = self.rotmat[-1]
        self.searcher.radsqd = self.radsqd
        # octant search
        self.searcher.noct = self.noct

        self.searcher.setup()
        self.searcher.pickup()

    def _block_discretization(self):
        self.nxdis = 1 if self.nxdis < 1 else self.nxdis
        self.nydis = 1 if self.nydis < 1 else self.nydis
        self.nzdis = 1 if self.nzdis < 1 else self.nzdis
        ndb = self.nxdis * self.nydis * self.nzdis
        if ndb > self.const.MAXDIS:
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

    def _cova3(self, point1, point2):
        """
        INPUT VARIABLES:
          x1,y1,z1         coordinates of first point
          x2,y2,z2         coordinates of second point
          nst(ivarg)       number of nested structures (maximum of 4)
          ivarg            variogram number (set to 1 unless doing cokriging
                              or indicator kriging)
          MAXNST           size of variogram parameter arrays
          c0(ivarg)        isotropic nugget constant
          it(i)            type of each nested structure:
                             1. spherical model of range a;
                             2. exponential model of parameter a;
                                  i.e. practical range is 3a
                             3. gaussian model of parameter a;
                                  i.e. practical range is a*sqrt(3)
                             4. power model of power a (a must be gt. 0  and
                                  lt. 2).  if linear model, a=1,c=slope.
                             5. hole effect model
          cc(i)            multiplicative factor of each nested structure.
                             (sill-c0) for spherical, exponential,and gaussian
                             slope for linear model.
          aa(i)            parameter "a" of each nested structure.

        OUTPUT VARIABLES:
          cmax             maximum covariance
          cova             covariance between (x1,y1,z1) and (x2,y2,z2)
        """
        # Calculate the maximum covariance value (used for zero distances and
        # for power model covariance):
        cmax = self.c0
        for ist in xrange(self.nst):
            if self.it[ist] == 4:
                cmax += self.const.PMX
            else:
                cmax += self.cc[ist]
        # check for 'zero' distance, return cmax if so:
        hsqd = self.sqdist(point1, point2)
        if hsqd < np.finfo(float).eps:
            cova = cmax
            return cova
        # loop over all structures
        cova = 0
        for ist in xrange(self.nst):
            h = np.sqrt(hsqd)
            if self.it[ist] == 1:  # Spherical
                hr = h / self.aa_hmax[ist]
                if hr < 1:
                    cova += self.cc[ist] * (1 - hr * (1.5 - 0.5 * hr * hr))
            elif self.it[ist] == 2:  # Exponential
                cova += self.cc[ist] * np.exp(-3.0 * h / self.aa_hmax[ist])
            elif self.it[ist] == 3:  # Gaussian
                cova += self.cc[ist] * \
                        np.exp(-3.0 * (h / self.aa_hmax[ist]) *
                               (h/self.aa_hmax[ist]))
            elif self.it[ist] == 4:  # Power
                cova += cmax - self.cc[ist] * (h**(self.aa_hmax[ist]))
            elif self.it[ist] == 5:  # Hole Effect
                cova += self.cc[ist] * np.cos(h / self.aa_hmax[ist] * np.pi)

    def sqdist(self, point1, point2):
        """
        This routine calculates the anisotropic distance between two points
        given the coordinates of each point and a definition of the
        anisotropy.

        Parameters
        ----------
        point1 : tuple
            Coordinates of first point (x1,y1,z1)
        point2 : tuple
            Coordinates of second point (x2,y2,z2)

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
            cont = self.rotmat[i, 0] * dx + \
                   self.rotmat[i, 1] * dy + \
                   self.rotmat[i, 2] * dz
            sqdist += cont * cont
        return sqdist

    def view2d(self):
        "View 2D data using matplotlib"
        if self._2d is False:
            print("3D data, use view3d() instead.")
        else:
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

    def view3d(self):
        "View 3D data using mayavi"
        pass

if __name__ == "__main__":
    test_krige3d = Krige3d("testData/test_krige3d.par")
    test_krige3d.read_data()
    test_krige3d.kt3d()
    test_krige3d.view2d()
    test_krige3d.view3d()

