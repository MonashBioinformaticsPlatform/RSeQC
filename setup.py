#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

import os
from setuptools import setup

readme = open('README.rst', 'r')
README_TEXT = readme.read()

setup(
    name = "RSeQC",
    version = "2.6.4",
    author = "Liguo Wang",
    author_email ="wangliguo78@gmail.com",
    description = ("RNA-seq QC Package"),
    long_description=README_TEXT,
    license = "GPLv2",
    keywords='RNA-seq, RNAseq, QC, metrics',
    url = "https://github.com/MonashBioinformaticsPlatform/RSeQC",
    scripts = ['scripts/rseqc'],
    install_requires = [
        'numpy',
        'cython>=0.17',
        'pysam',
        'bx-python'
        ], 
    packages=['rseqc']
    )

#platforms = ['Linux','MacOS'],
#classifiers=[
#    "Development Status :: 3 - Alpha",
#    "Topic :: Utilities",
#    "License :: OSI Approved :: BSD License",
