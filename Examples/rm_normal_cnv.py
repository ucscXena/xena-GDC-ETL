#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""A script for removing all normal samples from TCGA'S CNV data.

Normal samples in TCGA have sample codes 10, 11, 12, 13 or 14 in their 
barcodes, i.e. sample submitter ID.
"""

# Ensure Python 2 and 3 compatibility
from __future__ import division
from __future__ import print_function

import logging
import os
import shutil
import timeit

def main():
    start_time = timeit.default_timer()
    
    root_dir = os.path.abspath('/home/yunhai/gdc/xena/files')
    # Get all project_ids on GDC, and convert unicode to str for python 2
    projects = [str(x) for x in gdc.get_all_project_info().index]
    # Remove TARGET projects which don't have CNV data
    for project_id in ['TARGET-OS', 'TARGET-CCSK', 'TARGET-WT', 'TARGET-AML', 
                       'TARGET-RT', 'TARGET-NBL']:
        projects.remove(project_id)
    
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
        print('Processing [{}/{}] projects: {}'.format(counts, total_projects, 
                                                       project))
        try:
            dataset = GDCOmicset(project, 'masked.cnv', root_dir)
            raw_data_dir = os.path.join(root_dir, project, 'GDC_Raw_Data', 
                                        'masked_cnv')
            normal_data_dir = os.path.join('/home/yunhai/gdc/TCGA_normal', 
                                           project, 'GDC_Raw_Data', 
                                           'masked_cnv')
            mkdir_p(normal_data_dir)
            dataset.raw_data_list = []
            for f in os.listdir(raw_data_dir):
                sample_id = f.split('.', 1)[0]
                if int(sample_id[-3:-1]) in range(10, 15):
                    shutil.move(os.path.join(raw_data_dir, f), normal_data_dir)
                    continue
                dataset.raw_data_list.append(os.path.join(raw_data_dir, f))
            dataset.transform().metadata()
        except Exception:
            msg = 'Fail to process CNV data for cohort {}.'.format(project)
            logger.warning(msg, exc_info=True)
            print(msg)

    end_time = timeit.default_timer()
    print('Finish in {} hours.'.format((end_time - start_time) / 60 / 60))

if __name__ == '__main__':
    if __package__ is None:
        import sys
        sys.path.insert(
            0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        )
        import gdc
        from xena_dataset import GDCOmicset, mkdir_p
    else:
        from .. import gdc
        from ..xena_dataset import GDCOmicset, mkdir_p
    main()
