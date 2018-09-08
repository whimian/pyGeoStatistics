#!/usr/bin/env python
"""
Created on Sep 8th 2018
"""
from distutils.core import setup
from setuptools import find_packages
import versioneer

CLASSIFIERS = [
    'Development Status :: 4 - Beta',
    'Intended Audience :: Developers',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Mathematics',
    'Topic :: Scientific/Engineering :: Physics',
    'Operating System :: Microsoft :: Windows',
    'Operating System :: POSIX',
    'Operating System :: Unix',
    'Operating System :: MacOS',
    'Natural Language :: English',
]

with open("README.md") as f:
    LONG_DESCRIPTION = ''.join(f.readlines())

setup(
    name="pyGeoStatistics",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    install_requires=[
        'scipy',
        'pandas',
        'numba',
        'matplotlib'
    ],
    packages=find_packages(exclude=['tests', 'testData']),
    author="Yu Hao",
    author_email="yuhao@live.cn",
    description="pyGeoStatistics: Geostatistics with Python",
    long_description=LONG_DESCRIPTION,
    license="MIT",
    keywords="geostatistics",
    url="https://github.com/whimian/pyGeoStatistics",
    download_url="https://github.com/whimian/pyGeoStatistics",
    classifiers=CLASSIFIERS,
    platforms=["Windows", "Linux", "Solaris", "Mac OS-X", "Unix"],
    zip_safe=False
)
