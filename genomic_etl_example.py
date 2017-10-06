#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""A example script for importing genomic data on GDC into Xena

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

import gdc
from xena_dataset import XenaDataset

def main():
    start_time = timeit.default_timer()
    
    root_dir = os.path.abspath('/home/yunhai/gdc/imported_GDC')
    # Get all project_ids on GDC, and convert unicode to str for python 2
    projects = [str(x) for x in gdc.get_all_project_info().index]
    # Selected types of datasets for Xena
    xena_dtypes = ['htseq.counts', 'htseq.fpkm', 'htseq.fpkm-uq', 'mirna', 
                   'masked.cnv', 'muse.snv', 'mutect2.snv', 
                   'somaticsniper.snv', 'varscan2.snv']
    
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
                dataset = XenaDataset(project, dtype, root_dir)
                dataset.download_gdc().transform().metadata()
            except Exception:
                msg = 'No {} data for cohort {}.'.format(dtype, project)
                logger.warning(msg, exc_info=True)
                print(msg)
    logging.shutdown()
    
    end_time = timeit.default_timer()
    print('Finish in {} hours.'.format((end_time - start_time) / 60 / 60))

if __name__ == '__main__':
    main()
