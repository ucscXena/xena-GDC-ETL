#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""A example script for importing phenotype data through GDC API into Xena

For a simple complete import, the script can be adapted from this example by 
modifying 3 variables: root_dir, projects and xena_dtypes. For example::

    script_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    root_dir = os.path.join(script_dir, 'gitignore', 'test')
    projects = ['TARGET-RT']

"""

# Ensure Python 2 and 3 compatibility
from __future__ import division
from __future__ import print_function

import logging
import os
import timeit

import numpy as np

def main():
    start_time = timeit.default_timer()
    
    root_dir = os.path.abspath('/home/yunhai/gdc/xena/files')
    # Get all project_ids on GDC, and convert unicode to str for python 2
    projects = [str(x) for x in gdc.get_project_info().index]
    script_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
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
        msg = '[{}/{}] Querying GDC phenotype data for projects: {}'
        print(msg.format(counts, total_projects, project))
        matrix = os.path.join(root_dir, project, 'Xena_Matrices', 
                              '{}.GDC_phenotype.tsv'.format(project))
        try:
            df = gdc.get_clinical_samples(project)
            sample_mask = df['samples.submitter_id'].map(
                    lambda s: s[-3:-1] not in ['10']
                )
            df = df[sample_mask].set_index('samples.submitter_id')
            df['TCGA_study'] = df['project.project_id'].apply(
                    lambda x: x[5:]
                ).map(gdc.TCGA_STUDY_ABBR)
            df['year_at_diagnosis'] = df['diagnoses.age_at_diagnosis'].apply(
                    lambda x: x if np.isnan(x) else int(round(x/365))
                )
            gdc.mkdir_p(os.path.dirname(matrix))
            df.to_csv(matrix, sep='\t', encoding='utf-8')
        except Exception:
            msg = 'Fail to get GDC phenotype data for cohort {}.'
            msg = msg.format(project)
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
    else:
        from .. import gdc
    main()
