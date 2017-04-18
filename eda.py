# -*- coding: utf-8 -*-
"""
Exploratory Data Analysis

Created on Mon Nov 07 2016
"""
from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist


class EDA():
    def __init__(self, datafile):
        self.datafl = datafile
        self.vr = None
        self.property_name = None
        self._2d = False

    def preview(self):
        print("Data File")
        print("---------")
        with open(self.datafl, 'r') as fin:
            for line in fin.readlines(20):
                print(line.strip())

    def read(self):
        data_list = None
        with open(self.datafl, 'r') as fin:
            data_list = fin.readlines()
        self.name = data_list[0].strip()
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

    def pdf(self, bins=15):
        hist, bin_edges = np.histogram(self.vr[self.property_name[0]],
                                       bins=bins)
        fig, ax = plt.subplots()
        ax.set_title("pdf")
        plt.bar(bin_edges[:-1], hist, width=bin_edges[1]-bin_edges[0],
                color='red', alpha=0.5)
        fig.show()

    def cdf(self):
        data = self.vr[self.property_name[0]]
        data = np.sort(data)
        cdf = np.arange(1, len(data) + 1) / len(data)
        fig, ax = plt.subplots()
        ax.set_title("cdf")
        ax.plot(data, cdf)
        fig.show()

    @property
    def maximum(self):
        return np.max(self.vr[self.property_name[0]])

    @property
    def minimum(self):
        return np.min(self.vr[self.property_name[0]])

    @property
    def mean(self):
        return np.mean(self.vr[self.property_name[0]])

    @property
    def variance(self):
        return np.var(self.vr[self.property_name[0]])

    @property
    def meadian(self):
        return np.median(self.vr[self.property_name[0]])

    @property
    def upper_quartile(self):
        index = None
        even = False
        length = self.vr.shape[0]
        if length / 2 == 0:
            even = True
        if even:
            index = int((3*length + 2)/4) + 1
        else:
            index = int((3*length + 3)/4) + 1
        return self.vr[self.property_name[0]][index]

    @property
    def lower_quartile(self):
        index = None
        even = False
        length = self.vr.shape[0]
        if length / 2 == 0:
            even = True
        if even:
            index = int((length + 2)/4)
        else:
            index = int((length + 1)/4)
        return self.vr[self.property_name[0]][index]
    @property
    def num(self):
        return self.vr.shape[0]

    def statistics(self):
        print("\nStatistics")
        print("-"*10)
        print("Number of Data: {}".format(self.num))
        print("Mean: {}".format(self.mean))
        print("Variance: {}".format(self.variance))
        print("-"*10)
        print("Minimum: {}".format(self.minimum))
        print("Lower Quartile: {}".format(self.lower_quartile))
        print("Median: {}".format(self.meadian))
        print("Upper Quartile: {}".format(self.upper_quartile))
        print("Maximum: {}".format(self.maximum))
        print("-"*10)

    def distance(self):
        num = self.vr.shape[0]
        dist = pdist(np.concatenate((self.vr['x'].reshape((num, 1)),
                                     self.vr['y'].reshape((num, 1))), axis=1))
        print("\nDistance\n"+"-"*8)
        print("Max distance: {}\nMin distance: {}".format(np.max(dist),
                                                          np.min(dist)))
        print("-"*8)
        print("X: {} - {}".format(np.min(self.vr['x']), np.max(self.vr['x'])))
        print("Y: {} - {}".format(np.min(self.vr['y']), np.max(self.vr['y'])))
        if self._2d is False:
            print("Z: {} - {}".format(np.min(self.vr['z']),
                                      np.max(self.vr['z'])))

    def view2d(self, pname=None):
        pname = self.property_name[0] if pname is None else pname
        if self._2d is False:
            print("3D data, use view3d() instead.")
        else:
            fig, ax = plt.subplots()
            abscissa = self.vr['x']
            ordinate = self.vr['y']
            sc = ax.scatter(abscissa, ordinate, c=self.vr[pname],
                            cmap='jet')
            ax.set(xlim=(np.min(abscissa), np.max(abscissa)),
                   ylim=(np.min(ordinate), np.max(ordinate)),
                   xlabel="X (m)", ylabel="Y (m)",
                   title="Data Scatter", aspect='equal',
                   facecolor='grey')
            fig.colorbar(sc)
            fig.show()

if __name__ == "__main__":
    eda = EDA("testData/test.gslib")
#    eda.preview()
    eda.read()
    eda.pdf()
    eda.cdf()
    eda.statistics()
    eda.distance()
    eda.view2d()
