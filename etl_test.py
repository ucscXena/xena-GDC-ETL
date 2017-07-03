#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Ensure Python 2 and 3 compatibility
from __future__ import division
from __future__ import print_function

import os
import timeit

import xena

def main():
    start_time = timeit.default_timer()
    script_dir = os.path.dirname(os.path.abspath(__file__))
    work_dir = os.path.join(script_dir, 'gitignore')
    matrix_dir = None
#    work_dir = os.path.join(matrix_dir, 'TCGA-BRCA', 'HTSeq-Counts')
    mode = 'both'
#    xena.import_gdc(project='TCGA-BRCA', dataset_type='rna_counts', 
#                    work_path=work_dir, mode=mode, matrix_dir=matrix_dir)
#    test_size_mb = 310.76 # TCGA-BRCA rna_couts
#    xena.import_gdc(project='TCGA-BRCA', dataset_type='muse', 
#                    work_path=work_dir, mode=mode, matrix_dir=matrix_dir)
#    xena.import_gdc(project='TCGA-BRCA', dataset_type='mutect2', 
#                    work_path=work_dir, mode=mode, matrix_dir=matrix_dir)
#    xena.import_gdc(project='TCGA-BRCA', dataset_type='somaticsniper', 
#                    work_path=work_dir, mode=mode, matrix_dir=matrix_dir)
#    xena.import_gdc(project='TCGA-BRCA', dataset_type='varscan2', 
#                    work_path=work_dir, mode=mode, matrix_dir=matrix_dir)
#    test_size_mb = 81.33 # TCGA-BRCA muse, mutect2, somaticsniper and varscan2
    xena.import_gdc(project='TCGA-CHOL', dataset_type='rna_counts', 
                    work_path=work_dir, mode=mode, matrix_dir=matrix_dir)
#    test_size_mb = 11.22 # TCGA-CHOL rna_counts
    xena.import_gdc(project='TCGA-CHOL', dataset_type='rna_fpkm', 
                    work_path=work_dir, mode=mode)
    xena.import_gdc(project='TCGA-CHOL', dataset_type='rna_fpkmuq', 
                    work_path=work_dir, mode=mode)
    xena.import_gdc(project='TCGA-CHOL', dataset_type='mirna', 
                    work_path=work_dir, mode=mode)
    xena.import_gdc(project='TCGA-CHOL', dataset_type='mirna_isoform', 
                    work_path=work_dir, mode=mode)
    xena.import_gdc(project='TCGA-CHOL', dataset_type='masked_cnv', 
                    work_path=work_dir, mode=mode)
    test_size_mb = 75.18 # TCGA-CHOL 6 types
#    xena.import_gdc(project='TCGA-DLBC', dataset_type='muse', 
#                    work_path=work_dir, mode=mode, matrix_dir=matrix_dir)
#    test_size_mb = 1.29 # TCGA-DLBC muse
    end_time = timeit.default_timer()
    print('Finish in {} sec.'.format(end_time - start_time))
    total_size_mb = 23870
    est_time = total_size_mb / test_size_mb * (end_time - start_time) / 60 / 60
    print('Expect to process ({}) all data in {} hours.'.format(mode, est_time))

if __name__ == '__main__':
    main()
