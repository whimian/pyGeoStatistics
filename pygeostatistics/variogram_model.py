# -*- coding: utf-8 -*-
"""
Created on Fri Nov 2016

Five different kinds of variogram model
-----
use parameter name crange instead of range which is a built-in python function
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


def spherical(lag, sill, crange):
    if lag <= crange:
        return sill*(1.5*(lag/crange) - 0.5*(lag/crange)**3)
    else:
        return sill


def exponential(lag, sill, crange):
    return sill*(1 - np.exp(-(3*lag/crange)))


def gaussian(lag, sill, crange):
    return sill*(1 - np.exp(-(3*lag**2/crange**2)))


def power(lag, sill, omega):
    return sill*lag**omega


def hole_effect(lag, sill, crange):
    return sill*(1-np.cos((lag/crange)*np.pi))

if __name__ == '__main__':
    func = np.vectorize(exponential)
    abscissa = np.arange(0, 100, 0.1)
    ordinate = func(abscissa, 1, 40)
    plt.plot(abscissa, ordinate)
