# -*- coding: utf-8 -*-
"""
Test

Created on Sep. 4th 2018
"""
__author__ = "yuhao"

import pytest
from pygeostatistics.eda import EDA


def test__EDA():
    eda = EDA("testData/test.gslib")
    eda.read()
    assert eda.maximum == 16.9583
    assert eda.minimum == 12.1491
    assert float("{:.4f}".format(eda.mean)) == 14.6959
    assert float("{:.4f}".format(eda.variance)) == 0.7776
    assert eda.meadian == 14.6515
