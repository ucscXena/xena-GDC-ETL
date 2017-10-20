#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""A example script for importing genomic data on GDC into Xena

For a simple complete import, the script can be adapted from this example by 
modifying 3 variables: root_dir, projects and xena_dtypes. For example::

    script_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    root_dir = os.path.join(script_dir, 'gitignore', 'test')
    projects = ['TCGA-CHOL', 'TARGET-RT']
    xena_dtypes = ['htseq_counts', 'htseq_fpkm', 'htseq_fpkm-uq', 'mirna']

"""

# Ensure Python 2 and 3 compatibility
from __future__ import division
from __future__ import print_function

import logging
import os
import timeit

def main():
    start_time = timeit.default_timer()
    
    root_dir = os.path.abspath('/home/yunhai/gdc/xena/files')
    # Get all project_ids on GDC, and convert unicode to str for python 2
    projects = [str(x) for x in gdc.get_project_info().index]
    # Selected types of datasets for Xena
    xena_dtypes = ['htseq_counts', 'htseq_fpkm', 'htseq_fpkm-uq', 'mirna', 
                   'masked_cnv', 'muse_snv', 'mutect2_snv', 
                   'somaticsniper_snv', 'varscan2_snv']
    
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
        print('Importing [{}/{}] projects: {}'.format(counts, total_projects, 
                                                      project))
        for dtype in xena_dtypes:
            try:
                dataset = GDCOmicset(project, dtype, root_dir)
                dataset.download().transform().metadata()
            except Exception:
                msg = 'No {} data for cohort {}.'.format(dtype, project)
                logger.warning(msg, exc_info=True)
                print(msg)
    logging.shutdown()
    
    end_time = timeit.default_timer()
    print('Finish in {} hours.'.format((end_time - start_time) / 60 / 60))

if __name__ == '__main__':
    if __package__ is None:
        import sys
        sys.path.insert(
            0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        )
        import gdc
        from xena_dataset import GDCOmicset
    else:
        from .. import gdc
        from ..xena_dataset import GDCOmicset
    main()
