from setuptools import setup, find_packages

setup(
    name = "RSeQC",
    version = "2.6.4",
    url = "https://github.com/MonashBioinformaticsPlatform/RSeQC",
    license = "GPLv2",
    author = "kizza_a",
    author_email ="kirill.tsyganov@monash.edu",
    description = "RNAseq QCs suite",
    packages=find_packages(exclude=['test']),
    zip_safe=False
    keywords = 'RNA-seq, RNAseq, QC, metrics',
    scripts = ['scripts/rseqc'],
    install_requires = ['pysam',
                        'bx-python'
                        ], 
    )

#setup_requires=['numpy'],
#platforms = ['Linux','MacOS'],
#classifiers=[
#    "Development Status :: 3 - Alpha",
#    "Topic :: Utilities",
#    "License :: OSI Approved :: BSD License",
