#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""A example script for importing GDC data into Xena

For a simple complete import, the script can be adapted from this example by 
modifying 3 variables: root_dir, projects and xena_dtypes. For example::

    import os
    script_dir = os.path.dirname(os.path.abspath(__file__))
    root_dir = os.path.join(script_dir, 'gitignore', 'test')
    projects = ['TARGET-RT']
    xena_dtypes = ['htseq.counts', 'htseq.fpkm', 'htseq.fpkm-uq', 'mirna']

"""

# Ensure Python 2 and 3 compatibility
from __future__ import division
from __future__ import print_function

import sys
import timeit
import traceback

import gdc
from xena_dataset import XenaDataset

def main():
    start_time = timeit.default_timer()
    
    root_dir = '/home/yunhai/gdc/imported_GDC'
    # Get all project_ids on GDC, and convert unicode to str for python 2
    projects = [str(x) for x in gdc.get_all_project_info().index]
    # Selected types of datasets for Xena
    xena_dtypes = ['htseq.counts', 'htseq.fpkm', 'htseq.fpkm-uq', 'mirna', 
                   'masked.cnv', 'muse.snv', 'mutect2.snv', 
                   'somaticsniper.snv', 'varscan2.snv']
    
    counts = 0
    total_projects = len(projects)
    for project in projects:
        counts += 1
        print('Importing [{}/{}] projects: {}'.format(counts, total_projects, 
                                                      project))
        for dtype in xena_dtypes:
            try:
                cohort = XenaDataset(project, dtype, root_dir)
                cohort.download_gdc().transform().metadata()
            except Exception:
                print('\nNo {} data for cohort {}.'.format(dtype, project), 
                      file=sys.stderr)
                traceback.print_exc(file=sys.stderr)

    end_time = timeit.default_timer()
    print('Finish in {} hours.'.format((end_time - start_time) / 60 / 60))

if __name__ == '__main__':
    main()
