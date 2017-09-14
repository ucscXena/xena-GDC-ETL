#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""A script for importing GDC data into Xena
"""

# Ensure Python 2 and 3 compatibility
from __future__ import division
from __future__ import print_function

import os
import timeit

import gdc
from xena_dataset import XenaDataset

def main():
    start_time = timeit.default_timer()
    root_dir = '/home/yunhai/gdc/imported_GDC'
    projects = gdc.get_all_project_info().index.tolist()
    # Convert unicode to str for python 2
    projects = [str(x) for x in projects]
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
            except:
                print('No {} data for cohort {}.'.format(dtype, project))
    end_time = timeit.default_timer()

    print('Finish in {} hours.'.format((end_time - start_time) / 60 / 60))

if __name__ == '__main__':
    main()
