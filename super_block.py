# -*- coding: utf-8 -*-
"""
class for performing Super Block Search

Created on Tue Nov 22 2016
"""
from __future__ import division, print_function
from itertools import izip, product
import numpy as np
# from numba import jit


class SuperBlockSearcher(object):
    """
    Class for performing Super Block Search

    This subroutine sets up a 3-D "super block" model and orders the data
    by super block number.  The limits of the super block is set to the
    minimum and maximum limits of the grid; data outside are assigned to
    the nearest edge block.

    The idea is to establish a 3-D block network that contains all the
    relevant data. The data are then sorted by their index location in
    the search network, i.e., the index location is given after knowing
    the block index in each coordinate direction (ix,iy,iz):
            ii = (iz-1)*nxsup*nysup + (iy-1)*nxsup + ix
    An array, the same size as the number of super blocks, is constructed
    that contains the cumulative number of data in the model. With this
    array it is easy to quickly check what data are located near any given
    location.

    Parameters
    ----------
    nx,xmn,xsiz      Definition of the X grid being considered
    ny,ymn,ysiz      Definition of the Y grid being considered
    nz,zmn,zsiz      Definition of the Z grid being considered
    vr(nd)           x, y, z, other variables
    MAXSB[X,Y,Z]     Maximum size of super block network

    """
    def __init__(self):
        # grid definition
        self.nx = None
        self.xmn = None
        self.xsiz = None
        self.ny = None
        self.ymn = None
        self.ysiz = None
        self.nz = None
        self.zmn = None
        self.zsiz = None
        # data
        self.vr = None  # with x,y,z and variable and other secondary variable
        self.MAXSB = []

        # rotation matrix
        self.rotmat = None  # rotation matrix for searching!!!
        self.radsqd = None  # squared search radius

        # octant search
        self.noct = None  #  the number of data noct to retain from each octant

        #To be calculated
        self.nisb = None  # array with cumulative number of data in each super block.
        # super block definitions
        self.nxsup = None
        self.xmnsup = None
        self.xsizsup = None
        self.nysup = None
        self.ymnsup = None
        self.ysizsup = None
        self.nzsup = None
        self.zmnsup = None
        self.zsizsup = None
        # superblocks to search
        self.nsbtosr = None    # Number of super blocks to search
        self.ixsbtosr = None   # X offsets for super blocks to search
        self.iysbtosr = None   # Y offsets for super blocks to search
        self.izsbtosr = None   # Z offsets for super blocks to search
        # points found within nearby super blocks
        self.nclose = None
        self.close_samples = None
        self.infoct = None
        # output sort_index
        self.sort_index = None

    def _super_flat_index(self, ixsup, iysup, izsup):
        return ixsup + iysup * self.nxsup + izsup * self.nxsup * self.nysup

    def setup(self):
        """
        Variables estimated
        -------------------
        nisb()                Array with cumulative number of data in each
                                super block.
        nxsup,xmnsup,xsizsup  Definition of the X super block grid
        nysup,ymnsup,ysizsup  Definition of the Y super block grid
        nzsup,zmnsup,zsizsup  Definition of the Z super block grid
        """
        # Establish super block definition
        self.nxsup = min(self.nx, self.MAXSB[0])
        self.nysup = min(self.ny, self.MAXSB[1])
        self.nzsup = min(self.nz, self.MAXSB[2])

        self.xsizsup = self.nx * self.xsiz / self.nxsup
        self.ysizsup = self.ny * self.ysiz / self.nysup
        self.zsizsup = self.nz * self.zsiz / self.nzsup

        self.xmnsup = (self.xmn - 0.5 * self.xsiz) + 0.5 * self.xsizsup
        self.ymnsup = (self.ymn - 0.5 * self.ysiz) + 0.5 * self.ysizsup
        self.zmnsup = (self.zmn - 0.5 * self.zsiz) + 0.5 * self.zsizsup

        # partition data into each super block
        x_block = np.arange(self.xmnsup - 0.5 * self.xsizsup,
                            self.xmnsup + (self.nxsup + 1) * self.xsizsup + 1,
                            self.xsizsup)
        x_index = np.searchsorted(x_block, self.vr['x']) - 1

        y_block = np.arange(self.ymnsup - 0.5 * self.ysizsup,
                            self.ymnsup + (self.nysup + 1) * self.ysizsup + 1,
                            self.ysizsup)
        y_index = np.searchsorted(y_block, self.vr['y']) - 1

        z_block = np.arange(self.zmnsup - 0.5 * self.zsizsup,
                            self.zmnsup + (self.nzsup + 1) * self.zsizsup + 1,
                            self.zsizsup)
        z_index = np.searchsorted(z_block, self.vr['z']) - 1

        # self.super_block = np.full((self.nxsup, self.nysup, self.nzsup), [])
        temp = np.zeros_like(self.vr['x'])
        self.nisb = np.zeros((self.nxsup*self.nysup*self.nzsup,))
        for idx, (ix, iy, iz) in enumerate(izip(x_index, y_index, z_index)):
            ii = self._super_flat_index(ix, iy, iz)
            temp[idx] = ii
            self.nisb[ii] += 1

        # sort data by asceding super block number:
        self.sort_index = np.argsort(temp)
        self.vr = self.vr[self.sort_index]
        # set up nisb
        self.nisb = np.cumsum(self.nisb, dtype=np.int)

    def pickup(self):
        """
        This subroutine establishes which super blocks must be searched given
        that a point being estimated/simulated falls within a super block
        centered at 0,0,0.

        Variables estimated
        -------------------
        nsbtosr          Number of super blocks to search
        ixsbtosr         X offsets for super blocks to search
        iysbtosr         Y offsets for super blocks to search
        izsbtosr         Z offsets for super blocks to search
        """
        self.nsbtosr = 0
        self.ixsbtosr = list()
        self.iysbtosr = list()
        self.izsbtosr = list()
        for i, j, k in product(xrange(-(self.nxsup-1), self.nxsup),
                               xrange(-(self.nysup-1), self.nysup),
                               xrange(-(self.nzsup-1), self.nzsup)):
            xo = i * self.xsizsup
            yo = j * self.ysizsup
            zo = k * self.zsizsup
            shortest = np.finfo(float).max
            for i1, j1, k1 in product([-1, 1], [-1, 1], [-1, 1]):
                for i2, j2, k2 in product([-1, 1], [-1, 1], [-1, 1]):
                    xdis = (i1 - i2) * 0.5 * self.xsizsup + xo
                    ydis = (j1 - j2) * 0.5 * self.ysizsup + yo
                    zdis = (k1 - k2) * 0.5 * self.zsizsup + zo
                    hsqd = sqdist((0, 0, 0), (xdis, ydis, zdis), self.rotmat)
                    shortest = hsqd if hsqd < shortest else shortest
            if shortest <= self.radsqd:
                self.nsbtosr += 1
                self.ixsbtosr.append(i)
                self.iysbtosr.append(j)
                self.izsbtosr.append(k)

    def search(self, xloc, yloc, zloc):
        """
        Variables estimated
        -------------------
        nclose           Number of close data
        close()          Index of close data
        infoct           Number of informed octants (only computes if
                             performing an octant search)
        """
        ix, iy, iz = self.getindx(xloc, yloc, zloc)
        self.nclose = 0
        self.close_samples = list()
        distance = list()
        # loop over all super blocks
        for isup in xrange(self.nsbtosr):
            ixsup = ix + self.ixsbtosr[isup]
            iysup = iy + self.iysbtosr[isup]
            izsup = iz + self.izsbtosr[isup]
            if ixsup < 0 or ixsup >= self.nxsup or \
                    iysup < 0 or iysup >= self.nysup or \
                    izsup < 0 or izsup >= self.nzsup:
                continue
            # find number of points within this super block
            ii = self._super_flat_index(ixsup, iysup, izsup)
            i = None
            if ii == 0:
                nums = self.nisb[ii]
                i = 0
            else:
                nums = self.nisb[ii] - self.nisb[ii - 1]
                i = self.nisb[ii - 1]
            # loop over all the data within this super block
            for k in xrange(0, nums):
                # hsqd = self.sqdist((xloc, yloc, zloc),
                #                    (self.vr['x'][i], self.vr['y'][i],
                #                     self.vr['z'][i]))
                hsqd = sqdist((xloc, yloc, zloc),
                                   (self.vr['x'][i], self.vr['y'][i],
                                    self.vr['z'][i]), self.rotmat)
                if hsqd > self.radsqd:
                    continue
                self.nclose += 1
                self.close_samples.append(i)
                distance.append(i)
                i += 1
        # sort nearby samples by distance
        distance = np.array(distance)
        self.close_samples = np.array(self.close_samples)
        sort_index = np.argsort(distance)
        self.close_samples = self.close_samples[sort_index]
        if self.noct <= 0:
            return
        else:  # partition the data into octant
            inoct = np.zeros((8,))
            # pick up the closes samples in each octant
            nt = self.noct * 8
            na = 0
            for j in xrange(self.nclose):
                i = int(self.close_samples[j])
                h = distance[j]
                dx = self.vr['x'][i] - xloc
                dy = self.vr['y'][i] - yloc
                dz = self.vr['z'][i] - zloc
                if dz >= 0:
                    iq = 3
                    if dx <= 0 and dy > 0:
                        iq = 0
                    if dx > 0 and dy >= 0:
                        iq = 1
                    if dx < 0 and dy <= 0:
                        iq = 2
                else:
                    iq = 7
                    if dx <= 0 and dy > 0:
                        iq = 4
                    if dx > 0 and dy >= 0:
                        iq = 5
                    if dx < 0 and dy >= 0:
                        iq = 6
                inoct[iq] += 1

                if inoct[iq] <= self.noct:
                    self.close_samples[na] = i
                    distance[na] = h
                    if na == nt:
                        self.nclose = na
                        break
                    na += 1
        # how many octants from which samples are drawn
        self.infoct = np.count_nonzero(inoct)

    def getindx(self, xloc, yloc, zloc):
        """
        determine which superblock are the given point or list of points in

        Parameters
        ----------
        xloc, yloc, zloc: scalar or 1-D ndarray
        """
        x_block = np.arange(self.xmnsup - 0.5 * self.xsizsup,
                            self.xmnsup + (self.nxsup + 1) * self.xsizsup + 1,
                            self.xsizsup)
        x_index = np.searchsorted(x_block, xloc) - 1

        y_block = np.arange(self.ymnsup - 0.5 * self.ysizsup,
                            self.ymnsup + (self.nysup + 1) * self.ysizsup + 1,
                            self.ysizsup)
        y_index = np.searchsorted(y_block, yloc) - 1

        z_block = np.arange(self.zmnsup - 0.5 * self.zsizsup,
                            self.zmnsup + (self.nzsup + 1) * self.zsizsup + 1,
                            self.zsizsup)
        z_index = np.searchsorted(z_block, zloc) - 1

        return (x_index, y_index, z_index)

# @jit
def sqdist(point1, point2, rotmat):
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
        cont = rotmat[i, 0] * dx + \
               rotmat[i, 1] * dy + \
               rotmat[i, 2] * dz
        sqdist += cont * cont
    return sqdist
