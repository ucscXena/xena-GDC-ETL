#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""A example script for importing phenotype data through GDC API into Xena

For a simple complete import, the script can be adapted from this example by 
modifying 3 variables: root_dir, projects and xena_dtypes. For example::

    root_dir = os.path.abspath('/home/user/gdc/xena/files')
    projects = [str(x) for x in gdc.get_project_info().index]

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
            # Query GDC API for phenotype info
            df = gdc.get_clinical_samples(project)
            # Remove code 10, Blood Derived Normal, sample:
            # https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes
            sample_mask = df['samples.submitter_id'].map(
                    lambda s: s[-3:-1] not in ['10']
                )
            df = df[sample_mask].set_index('samples.submitter_id')
            # Get more readable study names:
            # https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations
            df['TCGA_study'] = df['project.project_id'].apply(
                    lambda x: x[5:]
                ).map(gdc.TCGA_STUDY_ABBR)
            # Convert age from days to years
            df['year_at_diagnosis'] = df['diagnoses.age_at_diagnosis'].apply(
                    lambda x: x if np.isnan(x) else int(round(x/365))
                )
            # Save matrix
            gdc.mkdir_p(os.path.dirname(matrix))
            df.to_csv(matrix, sep='\t', encoding='utf-8')
        except Exception:
            msg = 'Fail to get GDC phenotype data for cohort {}.'
            msg = msg.format(project)
            logger.warning(msg, exc_info=True)
            print(msg)
    logging.shutdown()
    
    end_time = timeit.default_timer()
    m, s = divmod(int(end_time - start_time), 60)
    h, m = divmod(m, 60)
    print('Finish in {:d}:{:02d}:{:02d}.'.format(h, m, s))

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
