# -*- coding: utf-8 -*-
"""
A cokriging program for a points or blocks on a regular grid.

Created on Fri Dec 2 2016
"""
from __future__ import division, print_function, absolute_import
import json
from itertools import product
import time
from collections import namedtuple
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
from super_block import SuperBlockSearcher

__author__ = "yuhao"

class Cokrige(object):
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

        self.nst = list()
        self.c0 = list()
        self.it = list()
        self.cc = list()
        self.ang1 = list()
        self.ang2 = list()
        self.ang3 = list()
        self.aa_hmax = list()
        self.aa_hmin = list()
        self.aa_vert = list()

    def _read_params(self):
        with open(self.param_file) as fin:
            params = json.load(fin)
            # data file definition
            self.datafl = params['datafl']  #: 'testData/test.gslib',
            self.nvr = params['nvar']  # number (primary + secondary)
            self.ixl = params['icolx']  #: 1,
            self.iyl = params['icoly']  #: 2,
            self.izl = params['icolz']
            self.ivrl = params['icolvr']  # list
            # data limits
            self.tmin = params['tmin']  #: -1.0e21,
            self.tmax = params['tmax']  #: 1.0e21,
            # collocated cokriging or not
            self.icolloc = params['icolloc'] # boolean

            # definition of collocated data file
            self.secfl = params['secfl']
            self.iclcol = params['iclcol']

            self.idbg = params['idbg']  #: 3,
            self.dbgfl = params['dbgfl']  #: 'kb2d.dbg',
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
            self.ndmin = params['ndmin']  # for both
            self.ndmaxp = params['ndmaxp']  # primary
            self.ndmaxs = params['ndmaxs']  # secondary
            # search radii for primary variable
            self.pradius_hmax = params['radius_hmax']  # scalar
            self.pradius_hmin = params['radius_hmin']  # scalar
            self.pradius_vert = params['radius_vert']  # scalar
            # search radii for secondary variables
            self.sradius_hmax = params['radius_hmax']  # scalar
            self.sradius_hmin = params['radius_hmin']  # scalar
            self.sradius_vert = params['radius_vert']  # scalar
            # search ellipsoid
            self.sang1 = params['sang1']  # scalar
            self.sang2 = params['sang2']  # scalar
            self.sang3 = params['sang3']  # scalar
            # kriging type
            self.ktype = params['ikrige']
            # mean values for primary and secondary variables
            self.vmean = params['mean']  # list
            # Vairography definition
            self.vario = params['vario']  # list of dictionaries

    def _fill_check_covariance(self):
        self.variography = [dict()] * self.nvr * self.nvr
        for var in self.vario:
            self.variography[(var['i']-1) * self.nvr + (var['j']-1)] = var
        # try fill in symmetric covariance element
        for i, j in product(range(self.nvr), range(self.nvr)):
            idx1 = i + j * self.nvr
            idx2 = j + i * self.nvr
            if idx1 == {} and idx2 == {}:
                raise ValueError("need variogram between {},{}".format(i, j))
            elif idx1 == {}:
                self.variography[idx1] = self.variography[idx2]
            elif idx2 == {}:
                self.variography[idx2] = self.variography[idx1]
        for var in self.variography:
            self.nst.append(var['nst'])
            self.c0.append(var['c0'])
            self.it.append(var['it'])
            for idx in range(var['nst']):
                self.cc.append(var['cc'][idx])
                self.ang1.append(var['ang1'][idx])
                self.ang2.append(var['ang2'][idx])
                self.ang3.append(var['ang3'][idx])
                self.aa_hmax.append(var['aa_hmax'][idx])
                self.aa_hmin.append(var['aa_hmin'][idx])
                self.aa_vert.append(var['aa_vert'][idx])

        # check linear model of coregionalization
        # check definite positiveness

    def _check_params(self):
        # Check search radius
        if self.pradius_hmax <= 0:
            raise ValueError("pradius_hmax should be larger than zero.")
        if self.sradius_hmax <= 0:
            raise ValueError("sradius_hmax should be larger than zero.")
        # Check data file definition
        if self.ixl < 0 and self.nx > 1:
            raise ValueError("WARNING: ixl=0 and nx>1 !")
        if self.iyl < 0 and self.ny > 1:
            raise ValueError("WARNING: iyl=0 and ny>1 !")
        if self.izl < 0 and self.nz > 1:
            raise ValueError("WARNING: izl=0 and nz>1 !")
        if self.ndmin <= 0:
            raise ValueError("ndmin too small")
        if self.ndmaxs/2 <= self.nvr and self.ktype == 2:
            print('WARNING: with traditional ordinary cokriging the '+\
                  'sum of the weights applied to EACH secondary data'+\
                  'is zero.  With ndmaxs set low and nvr large the'+\
                  'secondary data will not contribute to the estimate')

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
        cokrige_const = namedtuple('Cokrige_const',
                                   ['PMX', 'MAXNST', 'MAXSB', 'MAXDIS',
                                    'MAXSAM', 'UNEST', 'MAXVAR', 'MAXARG',
                                    'MAXCOK'])
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
        self.const = cokrige_const(
            PMX=999,
            MAXNST=4,
            MAXSB=(maxsbx, maxsby, maxsbz),
            MAXDIS=self.nxdis * self.nydis * self.nzdis,
            MAXSAM=self.ndmaxp + self.ndmaxs,
            UNEST=np.nan,
            MAXVAR=self.nvr,
            MAXARG=self.nvr*self.nvr,
            MAXCOK=(self.ndmaxp + self.ndmaxs)*self.nvr + self.nvr
            )
        # Calculate needed programing variables from input parameters
        self.pradsqd = self.pradius_hmax * self.pradius_hmax
        self.psanis1 = self.pradius_hmin / self.pradius_hmax
        self.psanis2 = self.pradius_vert / self.pradius_hmax

        self.sradsqd = self.sradius_hmax * self.sradius_hmax
        self.ssanis1 = self.sradius_hmin / self.sradius_hmax
        self.ssanis2 = self.sradius_vert / self.sradius_hmax

        self.anis1 = np.array(self.aa_hmin) / \
                     np.maximum(self.aa_hmax, np.finfo(float).eps)
        self.anis2 = np.array(self.aa_vert) / \
                     np.maximum(self.aa_hmax, np.finfo(float).eps)

        self._fill_check_covariance()

    def _set_rotation(self):
        """
        Set up rotation matrix for both anisotropy and searching.
        with self.rotmat being an array of 3*3 rotation matrix, the last matrix
        in the array are the searching matrix
        """
        ang1 = np.append(self.ang1, self.sang1)
        ang2 = np.append(self.ang2, self.sang2)
        ang3 = np.append(self.ang3, self.sang3)
        anis1 = np.append(self.anis1, self.psanis1)
        anis2 = np.append(self.anis2, self.psanis2)
        anis1 = np.append(anis1, self.ssanis1)
        anis2 = np.append(anis2, self.ssanis2)
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

    def krige(self):
        self._fill_check_covariance()
        self._preprocess()
        # Set up the rotation/anisotropy matrices needed for variogram
        # and searching
        self._set_rotation()
        # compute maximum covariance for the rescaling factor:
        self._max_covariance()
        # Set up for super block searching:
        print("Setting up Super Block Search...")
        self._create_searcher()
        # Set up discretization points per block
        self._block_discretization()
        # Find unbias value
        self.unbias = self.maxcov

        nxy = self.nx * self.ny
        nloop = self.nx * self.ny * self.nz
        print("Start working on the kriging...")
        # time
        t1 = time.time()
        ts = 0
        percent_od = 0
        self.estimation = np.full((nloop,), np.nan)
        self.estimation_variance = np.full((nloop,), np.nan)
        # MAIN LOOP
        for index in range(nloop):
            self.iz = index // nxy
            self.iy = (index - self.iz * nxy) // self.nx
            self.ix = index - self.iz * nxy - self.iy * self.nx
            xloc = self.xmn + self.ix * self.xsiz
            yloc = self.ymn + self.iy * self.ysiz
            zloc = self.zmn + self.iz * self.zsiz
            # Search for proximity data
            ts_1 = time.time()
            self.searcher.search(xloc, yloc, zloc)
            ts += time.time() - ts_1
            # load nearest data in xa, ya, za, vra, vea
            xa = list()
            ya = list()
            za = list()
            vra = list()
            iva = list()  # which variable
            npri = 0  # number of primary data
            nsec = 0  # number of secondary data
            na = 0  # number of both kinds
            for i in range(self.searcher.nclose):
                if npri == self.ndmaxp and nsec == self.ndmaxs:
                    continue
                idx = self.searcher.close_samples[i]
                # Load primary data
                prim = self.vr[self.property_name[0]][idx]
                if prim <= self.tmin and prim > self.tmax and \
                        npri < self.ndmaxp:
                    npri += 1
                    na += 1
                    xa.append(self.vr['x'][idx] - xloc)
                    ya.append(self.vr['y'][idx] - yloc)
                    za.append(self.vr['z'][idx] - zloc)
                    vra.append(prim)
                    iva.append(0)
                # Load secondary data
                sec1 = self.vr[self.property_name[1]][idx]
                if sec1 <= self.tmin and sec1 > self.tmax and \
                        nsec < self.ndmaxs:
                    nsec += 1
                    na += 1
                    xa.append(self.vr['x'][idx] - xloc)
                    ya.append(self.vr['y'][idx] - yloc)
                    za.append(self.vr['z'][idx] - zloc)
                    if self.ktype != 2:
                        vra.append(sec1 - self.vmean[1] - self.vmean[0])
                    else:
                        vra.append(sec1)
                    iva.append(1)
                sec2 = self.vr[self.property_name[2]][idx]
                if sec2 <= self.tmin and sec2 > self.tmax and \
                        nsec < self.ndmaxs:
                    nsec += 1
                    na += 1
                    xa.append(self.vr['x'][idx] - xloc)
                    ya.append(self.vr['y'][idx] - yloc)
                    za.append(self.vr['z'][idx] - zloc)
                    if self.ktype != 2:
                        vra.append(sec1 - self.vmean[2] - self.vmean[0])
                    else:
                        vra.append(sec1)
                    iva.append(2)
                sec3 = self.vr[self.property_name[3]][idx]
                if sec3 <= self.tmin and sec3 > self.tmax and \
                        nsec < self.ndmaxs:
                    nsec += 1
                    na += 1
                    xa.append(self.vr['x'][idx] - xloc)
                    ya.append(self.vr['y'][idx] - yloc)
                    za.append(self.vr['z'][idx] - zloc)
                    if self.ktype != 2:
                        vra.append(sec1 - self.vmean[3] - self.vmean[0])
                    else:
                        vra.append(sec1)
                    iva.append(3)

            est, estv = self._many_samples(xa, ya, za, vra, na)
            self.estimation[index] = est
            self.estimation_variance[index] = estv
            # print working percentage
            percent = np.round(index/nloop*100, decimals=0)
            dtime = time.time() - t1
            if percent != percent_od:
                print("{}% ".format(percent) +\
                  "."*20 + "{}s elapsed.".format(np.round(dtime, decimals=3)))
            percent_od = percent
        print("Kriging Finished.")
        print("Time used for searching: {}s".format(ts))

    def _many_samples(self, xa, ya, za, vra, na):
        if self.ktype == 0:
            neq = na
        elif self.ktype == 1:
            neq = na + 1
        elif self.ktype == 2:
            neq = na + self.nvr
        if (neq - na) > na or na < self.ndmin:
            print("not enough data.")
            return np.nan, np.nan
        # left side
        left = np.full((neq, neq), np.nan)
        # fill the kriging matrix:
        for i, j in product(range(na), range(na)):
            if np.isnan(left[j, i]):
                left[i, j] = self._cova3((xa[i], ya[i], za[i]),
                                         (xa[j], ya[j], za[j]))
            else:
                left[i, j] = left[j, i]


    @property
    def block_covariance(self):
        "return average covariance within block"
        if self._block_covariance is None:
            if self.ndb <= 1:  # point kriging
                self._block_covariance = self.unbias
            else:
                cov = list()
                for x1, y1, z1 in zip(self.xdb, self.ydb, self.zdb):
                    for x2, y2, z2 in zip(self.xdb, self.ydb, self.zdb):
                        cov.append(self._cova3((x1, y1, z1), (x2, y2, z2)))
                cov = np.array(cov).reshape((self.ndb, self.ndb))
                cov[np.diag_indices_from(cov)] -= self.c0
                self._block_covariance = np.mean(cov)
        return self._block_covariance

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
        for ist in range(self.nst):
            if self.it[ist] == 4:
                self.maxcov += self.const.PMX
            else:
                self.maxcov += self.cc[ist]

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

    def _cova3(self, point1, point2, ivarg):
        """
        Parameters
        ----------
        point1, point2: tuple of 3
            coordinates of two points
        ivarg: 0, 1, 2, 3
            0 for primary, 1,2,3 for secondary
        Returns
        -------
        cova: scalar
            covariance between (x1,y1,z1) and (x2,y2,z2)
        """
        # Calculate the maximum covariance
        istart = sum(self.nst[:ivarg])
        cmax = self.c0[ivarg]
        for iss in range(self.nst[ivarg]):
            ist = istart + iss
            if self.it[ist] == 4:
                cmax += self.const.PMX
            else:
                cmax += self.cc[ist]
        # check for 'zero' distance, return maxcov if so:
        hsqd = self._sqdist(point1, point2, self.rotmat[istart])
        if hsqd < np.finfo(float).eps:
            cova = cmax
            return cova

        # loop over all structures
        cova = 0
        for ist in range(istart, self.nst[ivarg]):
            if ist != 1:
                hsqd = self._sqdist(point1, point2, self.rotmat[ist])
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
                cova += self.maxcov - self.cc[ist] * (h**(self.aa_hmax[ist]))
            elif self.it[ist] == 5:  # Hole Effect
                cova += self.cc[ist] * np.cos(h / self.aa_hmax[ist] * np.pi)
        return cova

    def _sqdist(self, point1, point2, rotmat):
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

if __name__ == '__main__':
    test_cokrige = Cokrige("testData/test_cokrige.par")
    test_cokrige.read_data()
    test_cokrige.krige()
