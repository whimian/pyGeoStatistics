# -*- coding: utf-8 -*-
"""
Normal Score Transform

Created on Tue Dec 6 2016
"""
from __future__ import division, print_function, absolute_import
import json
from itertools import izip, product
import time
from collections import namedtuple
import numpy as np
from numba import jit
from scipy import linalg, interpolate


class NormalScoreTransform(object):
    def __init__(self, data, weights, zmin, zmax, ltail, ltpar, utail, utpar):
        """
        Perform Normal score transform of data.

        Atrributes
        ----------
        data: 1-d ndarray
            data to be transformed
        weights: 1-d ndarray
            declustering weights for transform
        zmin, zmax: float, float
            allowable values for back-transform
        ltail: {1, 2}
            option to handle values smaller than minimun data, lower tail
        ltpar: float
            parameter for lower tail option
        utail: int
            potion to handle values larger than maximum data, upper tail
        utpar: float
            parameter for upper tail option
        """
        self.data = np.array(data) # input data ndarray
        self.weights = np.array(weights) # input declustering weight ndarray

        self.transform_table = None
        self.zmin = zmin # allowable value for backtransform
        self.zmax = zmax # allowable value for backtransform
        self.ltail = ltail  # option to handle values less than vrg[0]:
        self.ltpar = ltpar  # parameter required for option ltail
        self.utail = utail  # option to handle values greater than vrg[-1]:
        self.utpar = utpar  # parameter required for option utail

    def _create_table(self):
        "create transformation lookup table"
        # sort input data by value
        sort_index = np.argsort(self.data)
        sorted_data = self.data[sort_index]
        sorted_weight = self.weights[sort_index]
        # compute cumulative probabilities
        weight_sum = np.sum(sorted_weight)
        cum_weight = np.cumsum(sorted_weight / weight_sum)
        cum_weight_old = np.append(np.array([0]), cum_weight[:-1])
        average = 0.5 * (cum_weight + cum_weight_old)
        # calculate normal score value:
        score = [gauinv(element) for element in average]
        # create lookup table
        table = [(da, sc) for da, sc in zip(sorted_data, score)]
        self.transform_table = np.array(table, dtype=np.dtype({
            'names': ['value', 'score'],
            'formats': ['f8'] * 2
        }))

    def create_transform_func(self):
        self._create_table()
        nrows = self.transform_table['value'].shape[0]
        self.forward_func = interpolate.interp1d(
            self.transform_table['value'].reshape((nrows,)),
            self.transform_table['score'].reshape((nrows,)),
            kind='linear',
            fill_value="extrapolate")

        self.back_func = interpolate.interp1d(
            self.transform_table['score'].reshape((nrows,)),
            self.transform_table['value'].reshape((nrows,)),
            kind='linear',
            fill_value="extrapolate")

    def transform(self, values):
        "transform data to nore score"
        return self.forward_func(values)

    def back_transform(self, scores):
        "transform nore score back to orginal data"
        values = np.full_like(scores, np.nan)

        lo_value = self.transform_table['value'][0]
        up_value = self.transform_table['value'][-1]
        lo_score = self.transform_table['score'][0]
        up_score = self.transform_table['score'][-1]
        # scores in normal range
        normal_mask = np.logical_and(scores <= up_score, scores >= lo_score)
        normal_scores = scores[normal_mask]
        values[normal_mask] = self.back_func(normal_scores)
        # scores in lower tail: 1=linear, 2=power
        lower_mask = scores < lo_score
        lower_scores = scores[lower_mask]
        temp = list()
        for sc in lower_scores:
            backtr = lo_value
            cdflo = gcum(lo_score)
            cdfbt = gcum(sc)
            if self.ltail == 1:  # linear
                backtr = powint(0, cdflo, self.zmin, lo_value, cdfbt, 1)
                temp.append(backtr)
            elif self.ltail == 2:  # power
                cpow = 1.0 / self.ltpar
                backtr = powint(0, cdflo, self.zmin, lo_value, cdfbt, cpow)
                temp.append(backtr)
        values[lower_mask] = temp
        # scores in upper tail: 1=linear, 2=power, 4=hyperbolic
        upper_mask = scores > up_score
        upper_scores = scores[upper_mask]
        temp = list()
        for sc in up_score:
            backtr = up_value
            cdfhi = gcum(up_score)
            cdfbt = gcum(sc)  # cdf value of the score to be back-transformed
            if self.utail == 1:  # linear
                backtr = powint(cdfhi, 1.0, up_value, self.zmax, cdfbt, 1)
                temp.append(backtr)
            elif self.utail == 2:  # power
                cpow = 1.0 / self.utpar
                backtr = powint(cdfhi, 1.0, up_value, self.zmax, cdfbt, cpow)
                temp.append(backtr)
            elif self.utail == 4:  # hyperbolic
                l = (up_value**self.utpar) * (1 - gcum(up_score))
                backtr = (l / (1 - gcum(sc)))**(1 / self.utpar)
                temp.append(backtr)
        values[upper_mask] = temp
        return values

@jit(nopython=True)
def gauinv(p):
    """
    Computes the inverse of the standard normal cumulative distribution
    function with a numerical approximation.

    Parameters
    ----------
    p : scalar, ndarray
        Cumulative probability funciton value

    Returns
    -------
    xp : scalar, ndarray
        Quantile function value

    Notes
    -----
    .. [1] Statistical Computing, by W.J. Kennedy, Jr. and James E. Gentle,
            1980, p. 95.
    """
    lim = 1.0e-10
    p0 = -0.322232431088
    p1 = -1.0
    p2 = -0.342242088547
    p3 = -0.0204231210245
    p4 = -0.0000453642210148
    q0 = 0.0993484626060
    q1 = 0.588581570495
    q2 = 0.531103462366
    q3 = 0.103537752850
    q4 = 0.0038560700634
    # check for an error situation
    if p < lim:
        return -1e10
    if p > 1 - lim:
        return 1e10
    pp = p
    if p > 0.5:
        pp = 1 - pp
    if p == 0.5:
        return 0
    y = np.sqrt(np.log(1.0/(pp*pp)))
    xp = y + ((((y*p4 + p3)*y + p2)*y + p1)*y + p0) / \
        ((((y*q4 + q3)*y + q2)*y + q1)*y + q0)
    if p == pp:
        xp = -xp

    return xp

@jit(nopython=True)
def gcum(x):
    """
    Evaluate the standard normal cdf given a normal deviate x.  gcum is
    the area under a unit normal curve to the left of x.  The results are
    accurate only to about 5 decimal places.
    """
    z = -x if x < 0 else x
    t = 1. / (1. + 0.2316419 * z)
    gcum = t*(0.31938153 + t*(-0.356563782 + t*(1.781477937 + \
           t*(-1.821255978 + t*1.330274429))))
    e2 = np.exp(-z*z/2.)*0.3989422803 if z <= 6 else 0
    gcum = 1.0 - e2 * gcum
    if x >= 0:
        return gcum
    else:
        return 1.0 - gcum

@jit(nopython=True)
def powint(xlow, xhigh, ylow, yhigh, value, power):
    "power interpolation"
    if xhigh-xlow < np.finfo(float).eps:
        return (yhigh + ylow) / 2.0
    else:
        return ylow + (yhigh - ylow) * \
               (((value - xlow) / (xhigh - xlow))**power)
