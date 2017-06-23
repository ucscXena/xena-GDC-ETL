#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Ensure Python 2 and 3 compatibility
from __future__ import division
from __future__ import print_function

import os
import timeit

import xena

start_time = timeit.default_timer()
script_dir = os.path.dirname(os.path.abspath(__file__))
work_dir = os.path.join(script_dir, r'gitignore')
mode = 'both'
xena.import_gdc(project='TCGA-BRCA', dataset_type='rna_counts', 
                work_dir=work_dir, mode=mode)
#xena.import_gdc(project='TCGA-CHOL', dataset_type='rna_counts', 
#                work_dir=work_dir, mode=mode)
#xena.import_gdc(project='TCGA-CHOL', dataset_type='rna_fpkm', 
#                work_dir=work_dir, mode=mode)
#xena.import_gdc(project='TCGA-CHOL', dataset_type='rna_fpkmuq', 
#                work_dir=work_dir, mode=mode)
#xena.import_gdc(project='TCGA-CHOL', dataset_type='mirna', 
#                work_dir=work_dir, mode=mode)
#xena.import_gdc(project='TCGA-CHOL', dataset_type='mirna_isoform', 
#                work_dir=work_dir, mode=mode)
#xena.import_gdc(project='TCGA-CHOL', dataset_type='masked_cnv', 
#                work_dir=work_dir, mode=mode)
end_time = timeit.default_timer()
print('Finish in {} sec.'.format(end_time - start_time))
#test_size_mb = 75.18 # TCGA-CHOL 6 types
test_size_mb = 310.76 # TCGA-BRCA rna_couts
total_size_mb = 23870
est_time = total_size_mb / test_size_mb * end_time / 60 / 60
print('Expect to process ({}) all data in {} hours.'.format(mode, est_time))
