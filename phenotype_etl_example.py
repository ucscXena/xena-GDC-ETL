#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""A example script for importing phenotype data on GDC into Xena

For a simple complete import, the script can be adapted from this example by 
modifying 3 variables: root_dir, projects and xena_dtypes. For example::

    script_dir = os.path.dirname(os.path.abspath(__file__))
    root_dir = os.path.join(script_dir, 'gitignore', 'test')
    projects = ['TARGET-RT']
    xena_dtypes = ['htseq.counts', 'htseq.fpkm', 'htseq.fpkm-uq', 'mirna']

"""

# Ensure Python 2 and 3 compatibility
from __future__ import division
from __future__ import print_function

import logging
import os
import timeit

import pandas as pd

from xena_dataset import XenaDataset

def main():
    start_time = timeit.default_timer()
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    root_dir = os.path.join(script_dir, 'gitignore', 'test')
    projects = ['TCGA-CHOL']
    
    counts = 0
    total_projects = len(projects)
    log_format = '%(asctime)-15s [%(levelname)s]: %(message)s'
    logging.basicConfig(level=logging.WARNING, format=log_format, 
                        datefmt='%Y-%m-%d %H:%M:%S',
                        filename=os.path.join(root_dir, 'etl.err'),
                        filemode='w')
    logger = logging.getLogger('Xena-GDC-ETL')
    for project in projects:
        counts += 1
        msg = '[{}/{}] Importing phenotype data for projects: {}'
        print(msg.format(counts, total_projects, project))
        try:
            biospecimen = XenaDataset(project, 'biospecimen', root_dir)
            clinical = XenaDataset(project, 'clinical', root_dir)
            biospecimen_matrix = biospecimen.download_gdc().transform().matrix
            clinical_matrix = clinical.download_gdc().transform().matrix
            phenotype = pd.merge(pd.read_table(biospecimen_matrix),
                                 pd.read_table(clinical_matrix),
                                 how='outer', on='bcr_patient_barcode')
            phenotype.set_index('bcr_sample_barcode', inplace=True)
            phenotype.to_csv(os.path.join(root_dir, project, 'Xena_Matrices', 
                                          '{}.phenotype.tsv'.format(project)), 
                             sep='\t')
        except Exception:
            msg = 'No phenotype data for cohort {}.'.format(project)
            logger.warning(msg, exc_info=True)
            print(msg)

    end_time = timeit.default_timer()
    print('Finish in {} hours.'.format((end_time - start_time) / 60 / 60))

if __name__ == '__main__':
    main()
