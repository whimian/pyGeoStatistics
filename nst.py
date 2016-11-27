# -*- coding: utf-8 -*-
"""
Created on Mon Jul 04 23:42:04 2016

@author: yuhao
"""
from __future__ import division
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def nst(data):
    """
    Perform Normal Score Transformation

    Parameters
    ----------
    data : ndarray

    Returns
    -------

    """
    data_sorted = np.sort(data)
    cdf_original = np.arange(1, data_sorted.size + 1) / data_sorted.size

    cdf0 = cdf_original[:-1]
    cdf0 = np.insert(cdf0, 0, 0)
    cdf_norm = (cdf_original + cdf0) / 2
#    plt.plot(data_sorted, cdf_original, data_sorted, cdf_norm, 'r*')
    y = list()
    for c in cdf_norm:
        y.append(gauinv(c))
    return y


def gauinv(p):
    """
    Computes the inverse of the standard normal cumulative distribution
    function with a numerical approximation.

    Parameters
    ----------
    p : float
        Cumulative probability funciton value

    Returns
    -------
    xp : float
        Quantile function value

    Notes
    -----
    .. [1] Statistical Computing, by W.J. Kennedy, Jr. and James E. Gentle,
            1980, p. 95.
    """
    # lim = 1.0e-10
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
    pp = p
    if p > 0.5:
        pp = 1 - pp
    xp = 0
    if p == 0.5:
        return xp

    y = np.sqrt(np.log(1.0/(pp*pp)))
    xp = y + ((((y*p4 + p3)*y + p2)*y + p1)*y + p0) / \
        ((((y*q4 + q3)*y + q2)*y + q1)*y + q0)
    if p == pp:
        xp = -xp

    return xp

if __name__ == "__main__":
    # import data
    z = open('testData/ZoneA.dat', 'r').readlines()
    z = [i.strip().split() for i in z[10:]]
    z = np.array(z, dtype=np.float)
    z = pd.DataFrame(z, columns=['x', 'y', 'thk', 'por', 'perm',
                                 'lperm', 'lpermp', 'lpermr'])
    points = np.array(z[['x', 'y', 'por']])
    por = points[:, 2]

#    fig, axes = plt.subplots(nrows=1, ncols=2)
    # pdf
#    hist, bin_edges = np.histogram(por, bins=15, normed=True)
#    axes[0].plot(bin_edges[1:], hist)
    # cdf
    # the traditional way of making probabilities suming to n/(n+1) less than 1
#    por_sorted = np.sort(por)
#    yvals = np.arange(por_sorted.size) / por_sorted.size
#    axes[1].plot(por_sorted, yvals)
#
#    points_f = pd.DataFrame(points, columns=['x', 'y', 'por'])
#
#    sorted_points = points_f.sort_values('por')
    y = nst(por)
    ycdf = np.arange(1, len(y) + 1) / len(y)
    plt.figure()
    plt.plot(y, ycdf)
