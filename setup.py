#!/usr/bin/env python
# -*- coding: utf-8 -*-

import glob
import os
from setuptools import setup, find_packages

with open('requirements.txt') as f_requires:
    install_requires = [req for req in f_requires if not req.startswith('#')]
with open('README.rst') as f_readme:
    long_description = f_readme.read()
SCRIPTS_DIR = os.path.join('xena_gdc_etl', 'scripts')
scripts = glob.glob(os.path.join(SCRIPTS_DIR, '*.py'))
scripts += glob.glob(os.path.join(SCRIPTS_DIR, '*.sh'))

setup(
    author='Yunhai Luo',
    author_email='yunhail@stanford.edu',
    classifiers=[
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'License :: OSI Approved :: Apache Software License',
        'Operating System :: OS Independent',
    ],
    description='Scripts for importing GDC data into UCSC Xena.',
    install_requires=install_requires,
    setup_requires=['pytest-runner'],
    tests_require=['pytest', 'pytest-cov', 'flake8'],
    long_description=long_description,
    long_description_content_type='text/x-rst',
    name='xena_gdc_etl',
    packages=find_packages(),
    package_data={'': ['resources/*']},
    scripts=scripts,
    entry_points={
        'console_scripts': [
            'xge = xena_gdc_etl.main:main',
        ],
    },
    license='Apache License 2.0',
    keywords='Xena Genomic Data Commons GDC',
    url='https://github.com/yunhailuo/xena-GDC-ETL',
    version='0.2.0',
)
