# -*- coding: utf-8 -*-
"""
Read gslib file format

Created on Wen Sep 5th 2018
"""
from __future__ import absolute_import, division, print_function

__author__ = "yuhao"

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist


class SpatialData(object):
    def __init__(self, file_path):
        self.datafl = file_path
        self.vr = None
        self.property_name = None
        self._2d = False
        self._read_data()

    def _read_data(self):
        """
        read gslib file
        """
        column_name = []
        with open(self.datafl, 'r') as fin:
            _ = fin.readline().strip()
            ncols = int(fin.readline().strip())
            for _ in range(ncols):
                column_name.append(fin.readline().strip())
        self.property_name = [item for item in column_name
                              if item not in ['x', 'y', 'z']]
        df = pd.read_csv(self.datafl, sep='\t', header=None, names=column_name,
                         skiprows=ncols+2)
        if 'z' not in column_name:
            self._2d = True
            column_name.append('z')
            df['z'] = 0
        self.df = df

        data_dtype = np.dtype({
            'names': column_name,
            'formats': ['f8'] * len(column_name)})

        self.vr = np.core.records.fromarrays(
            df.values.transpose(), dtype=data_dtype)

    def preview(self):
        return self.vr.head(20)

    def pdf(self, ax, bins=15):
        hist, bin_edges = np.histogram(self.vr[self.property_name[0]],
                                       bins=bins)
        ax.set_title("pdf")
        ax.bar(bin_edges[:-1], hist, width=bin_edges[1]-bin_edges[0],
               color='red', alpha=0.5)

    def cdf(self, ax):
        data = self.vr[self.property_name[0]]
        data = np.sort(data)
        cdf = np.arange(1, len(data) + 1) / len(data)
        ax.set_title("cdf")
        ax.plot(data, cdf)

    @property
    def maximum(self):
        return self.df[self.property_name[0]].max()

    @property
    def minimum(self):
        return self.df[self.property_name[0]].min()

    @property
    def mean(self):
        return self.df[self.property_name[0]].mean()

    @property
    def variance(self):
        return self.df[self.property_name[0]].var()

    @property
    def meadian(self):
        return np.median(self.vr[self.property_name[0]])

    @property
    def upper_quartile(self):
        return self.df[self.property_name[0]].quantile(0.75)

    @property
    def lower_quartile(self):
        return self.df[self.property_name[0]].quantile(0.25)

    @property
    def num(self):
        return self.vr.shape[0]

    def distance(self):
        num = self.vr.shape[0]
        return pdist(np.concatenate((self.vr['x'].reshape((num, 1)),
                                     self.vr['y'].reshape((num, 1))), axis=1))

    @property
    def summary(self):
        return (
            "Summary\n"
            "-------\n"
            "Number of Data: {}\n"
            "Mean: {}\n"
            "Variance: {}\n"
            "Minimum: {}\n"
            "Lower Quartile: {}\n"
            "Median: {}\n"
            "Upper Quartile: {}\n"
            "Maximum: {}\n").format(
                self.num,
                self.mean,
                self.variance,
                self.minimum,
                self.lower_quartile,
                self.meadian,
                self.upper_quartile,
                self.maximum)
