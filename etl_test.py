#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""A test module for ETL pipelines importing GDC data into Xena
"""

# Ensure Python 2 and 3 compatibility
from __future__ import division
from __future__ import print_function

import os
import timeit

import xena

def main():
    start_time = timeit.default_timer()
    script_dir = os.path.dirname(os.path.abspath(__file__))
    work_dir = os.path.join(script_dir, 'gitignore', 'test')
    projects = ['TCGA-CHOL', 'TCGA-UCS']
    dataset_types = ['rna_counts', 'rna_fpkm', 'rna_fpkmuq', 'mirna', 
                     'masked_cnv', 'muse', 'mutect2', 'somaticsniper', 
                     'varscan2']
    probemap_dict = {
            'rna': os.path.join(work_dir, 'probeMaps', 
                                'gencode.v22.annotation.gene.probeMap')
        }
    for project in projects:
        for dtype in dataset_types:
            xena.import_gdc(project, dtype, work_path=work_dir, 
                            probemap_path=probemap_dict)
    end_time = timeit.default_timer()

    print('Finish in {} sec.'.format(end_time - start_time))
    test_size_mb = 149.27
    total_size_mb = 23870
    est_time = total_size_mb / test_size_mb * (end_time - start_time) / 60 / 60
    print('Expect to ETL all data in {} hours.'.format(est_time))

if __name__ == '__main__':
    main()
