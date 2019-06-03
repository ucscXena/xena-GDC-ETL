#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""A specialized script for importing phenotype data for project "TARGET-CCSK"
from GDC to Xena

Clinical data for "TARGET-CCSK" has the column "TARGET USI" not matching the
"cases.submitter_id" on GDC. All IDs in that column miss prefixes. As of
10/20/2017, 13 out of 32 are present in GDC. All 13 should have a prefix of
"TARGET-51-".
"""

# Ensure Python 2 and 3 compatibility
from __future__ import division
from __future__ import print_function

import logging
import os
import timeit

import pandas as pd


def ccsk_clin_dfs2matrix(df_list):
    print('Merging into one matrix ...', end='')
    clin_df = pd.concat(df_list)
    clin_df = clin_df.rename(index=lambda x: 'TARGET-51-' + x)
    print('\rMapping clinical info to individual samples...', end='')
    cases = gdc.search(
        'cases',
        in_filter={'project.project_id': 'TARGET-CCSK'},
        fields=['submitter_id', 'samples.submitter_id'],
        typ='json',
    )
    cases_samples = [c for c in cases if 'samples' in c]
    from pandas.io.json import json_normalize

    cases_samples_map = json_normalize(
        cases_samples, 'samples', ['submitter_id'], meta_prefix='cases.'
    ).rename(
        columns={
            'submitter_id': 'sample_id',
            'cases.submitter_id': 'TARGET USI',
        }
    )
    return pd.merge(
        clin_df.reset_index(), cases_samples_map, how='inner', on='TARGET USI'
    ).set_index('sample_id')


def main():
    start_time = timeit.default_timer()

    root_dir = os.path.abspath('/home/yunhai/gdc/xena/files')

    log_format = '%(asctime)-15s [%(levelname)s]: %(message)s'
    logging.basicConfig(
        level=logging.WARNING,
        format=log_format,
        datefmt='%Y-%m-%d %H:%M:%S',
        filename=os.path.join(root_dir, 'etl.err'),
        filemode='w',
    )
    logger = logging.getLogger('Xena-GDC-ETL')
    print('Importing phenotype data for TARGET-CCSK')
    try:
        phenoset = TARGETPhenoset('TARGET-CCSK', root_dir)
        phenoset.raws2matrix = ccsk_clin_dfs2matrix
        phenoset.download().transform().metadata()
    except Exception:
        logger.warning('No phenotype data for "TARGET-CCSK".', exc_info=True)
        print('No phenotype data for "TARGET-CCSK".')
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
        from xena_dataset import TARGETPhenoset
    else:
        from .. import gdc
        from ..xena_dataset import TARGETPhenoset
    main()
