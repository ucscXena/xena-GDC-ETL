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

from xena_dataset import XenaTCGAPhenoset, XenaTARGETPhenoset

def main():
    start_time = timeit.default_timer()
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    root_dir = os.path.join(script_dir, 'gitignore', 'test')
    projects = ['TCGA-CHOL', 'TARGET-AML']
    
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
            if project.startswith('TCGA'):
                phenoset = XenaTCGAPhenoset(project, root_dir)
            if project.startswith('TARGET'):
                phenoset = XenaTARGETPhenoset(project, root_dir)
                print(phenoset.gdc_download_dict)
            phenoset.download_gdc().transform().metadata()
        except Exception:
            msg = 'No phenotype data for cohort {}.'.format(project)
            logger.warning(msg, exc_info=True)
            print(msg)
    logging.shutdown()
    
    end_time = timeit.default_timer()
    print('Finish in {} hours.'.format((end_time - start_time) / 60 / 60))

if __name__ == '__main__':
    main()
