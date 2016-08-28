#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

import os
from setuptools import setup
from setuptools import find_packages

setup(
    name = "RSeQC",
    version = "2.6.4",
    author = "Liguo Wang",
    author_email ="wangliguo78@gmail.com",
    description = "RNA-seq QC Package",
    license = "GPLv2",
    keywords = 'RNA-seq, RNAseq, QC, metrics',
    url = "https://github.com/MonashBioinformaticsPlatform/RSeQC",
    scripts = ['scripts/rseqc'],
    install_requires = [
        'pysam',
        'bx-python'
        ], 
    packages=find_packages(),
    zip_safe=False
    )

#setup_requires=['numpy'],
#platforms = ['Linux','MacOS'],
#classifiers=[
#    "Development Status :: 3 - Alpha",
#    "Topic :: Utilities",
#    "License :: OSI Approved :: BSD License",
