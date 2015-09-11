import os
import shutil
import sys
import glob
from setuptools import setup, find_packages

setup(
    name='primer3tools',
    version='0.0.2',
    description='primer3tools: Design PCR primers and check uniqueness across multiple genomes',
    packages = find_packages(),
    author='Martin Hunt',
    author_email='path-help@sanger.ac.uk',
    url='https://github.com/sanger-pathogens/primer3tools',
    scripts=glob.glob('scripts/*'),
    test_suite='nose.collector',
    tests_require=['nose >= 1.3'],
    install_requires=[
        'pyfastaq >= 3.7.0',
        'pysam >= 0.8.3',
    ],
    license='GPLv3',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3 :: Only',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    ],
)
