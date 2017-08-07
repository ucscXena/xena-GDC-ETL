#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""A test module for ETL pipelines importing GDC data into Xena
"""

# Ensure Python 2 and 3 compatibility
from __future__ import division
from __future__ import print_function

import os
import timeit

from xena_dataset import XenaDataset

def main():
    start_time = timeit.default_timer()
    script_dir = os.path.dirname(os.path.abspath(__file__))
    root_dir = os.path.join(script_dir, 'gitignore', 'test')
    projects = ['TCGA-CHOL', 'TCGA-UCS']
    xena_dtypes = ['htseq.counts', 'htseq.fpkm', 'htseq.fpkm-uq', 'mirna', 
                   'masked.cnv', 'muse.snv', 'mutect2.snv', 
                   'somaticsniper.snv', 'varscan2.snv']
#    probemap_dict = {
#            'rna': os.path.join(work_dir, 'probeMaps', 
#                                'gencode.v22.annotation.gene.probeMap')
#        }
    for project in projects:
        for dtype in xena_dtypes:
            try:
                cohort = XenaDataset(project, dtype, root_dir)
                cohort.download_gdc().transform().metadata()
            except:
                print('No {} data for cohort {}.'.format(dtype, project))
    end_time = timeit.default_timer()

    print('Finish in {} sec.'.format(end_time - start_time))
    test_size_mb = 149.27
    total_size_mb = 23870
    est_time = total_size_mb / test_size_mb * (end_time - start_time) / 60 / 60
    print('Expect to ETL all data in {} hours.'.format(est_time))

if __name__ == '__main__':
    main()
