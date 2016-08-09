#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

import os
from setuptools import setup

setup(
    name = "RSeQC",
    version = "2.6.4",
    author = "Liguo Wang",
    author_email ="wangliguo78@gmail.com",
    description = ("RNA-seq QC Package"),
    long_description=read('README'),
    license = "GPLv2",
    keywords='RNA-seq, RNAseq, QC, metrics',
    url = "https://github.com/MonashBioinformaticsPlatform/RSeQC"
    scripts = ['scripts/multiqc'],
    install_requires = [
        'numpy',
        'cython>=0.17',
        'pysam',
        'bx-python'
        ], 
    packages=['RSeQC'],

    platforms = ['Linux','MacOS'],
    requires = ['cython (>=0.17)'],

    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Utilities",
        "License :: OSI Approved :: BSD License",
    ], 
        )
