"""
Setup script for RSeQC -- Comprehensive QC package for RNA-seq data
"""

import sys, os, platform, glob
from distutils.core import setup
from setuptools import *


if sys.version_info[0] != 2 or sys.version_info[1] < 7:
	print >> sys.stderr, "ERROR: RSeQC requires Python 2.7"
	sys.exit()
	

def main():
	setup(  name = "RSeQC",
            version = "2.6.4",
            packages = find_packages( 'lib' ),
            package_dir = { '': 'lib' },
            package_data = { '': ['*.ps'] },
            scripts = glob.glob( "scripts/*.py"),
            ext_modules = [],
            py_modules = [ 'psyco_full' ],
            test_suite = 'nose.collector',
            setup_requires = ['nose>=0.10.4','cython>=0.12'],
            author = "Liguo Wang",
			author_email ="wangliguo78@gmail.com",
			platforms = ['Linux','MacOS'],
			requires = ['cython (>=0.17)'],
			install_requires = ['cython>=0.17','pysam','bx-python','numpy'], 
            description = "RNA-seq QC Package",
            url = "http://rseqc.sourceforge.net/",
            zip_safe = False,
            dependency_links = [],
			classifiers=[
              'Development Status :: 5 - Production/Stable',
              'Environment :: Console',
              'Intended Audience :: Science/Research',
              'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
              'Operating System :: MacOS :: MacOS X',
              'Operating System :: POSIX',
              'Programming Language :: Python',
              'Topic :: Scientific/Engineering :: Bio-Informatics',
              ],
			keywords='RNA-seq, QC',
        )


if __name__ == "__main__":
	main()
