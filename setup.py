#!/usr/bin/env python
# -*- coding: utf-8 -*-

import glob
import os
from setuptools import setup, find_packages

with open('README.rst', 'r') as fh:
    long_description = fh.read()

SCRIPTS_DIR = os.path.join('xena_gdc_etl', 'scripts')
scripts = glob.glob(os.path.join(SCRIPTS_DIR, '*.py'))
scripts += glob.glob(os.path.join(SCRIPTS_DIR, '*.sh'))

setup(
    author='Yunhai Luo',
    author_email='yunhail@stanford.edu',
    classifiers=[
        'Programming Language :: Python :: 2.7',
        'Operating System :: OS Independent',
    ],
    description='Scripts for importing GDC data into UCSC Xena.',
    install_requires=[
        'jinja2>=2.7.2',
        'lxml>=4.0.0',
        'numpy>=1.13.0',
        'pandas>=0.20.2',
        'requests>=2.20.0',
        'xlrd>=1.1.0',
    ],
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    long_description=long_description,
    long_description_content_type='text/x-rst',
    name='xena_gdc_etl',
    packages=find_packages(),
    package_data={'': ['resources/*']},
    scripts=scripts,
    url='https://github.com/yunhailuo/xena-GDC-ETL',
    version='0.2.0',
)
