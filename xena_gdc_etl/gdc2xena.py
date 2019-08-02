#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This script provides a command line tool for importing data from GDC to
Xena.

This module also provides an independent function ``gdc2xena`` for building an
importing pipeline directly by providing 3 arguments: root_dir, projects and
xena_dtypes. For example::

    script_dir = os.path.dirname(os.path.abspath(__file__))
    root_dir = os.path.join(script_dir, 'gitignore', 'test')
    projects = ['TCGA-CHOL', 'TARGET-RT']
    xena_dtypes = ['htseq_counts', 'phenotype', 'survival']
    gdc2xena(root_dir, projects, xena_dtypes)

Supported types of data are genomic data, phenotype data and survival data.
Please check help message with "-h" option for details.
"""

# Ensure Python 2 and 3 compatibility
from __future__ import division
from __future__ import print_function

import os
import logging
import timeit
import json
import time
import shutil

from .xena_dataset import (
    GDCOmicset,
    GDCPhenoset,
    GDCAPIPhenoset,
    GDCSurvivalset,
)


def gdc2xena(root_dir, projects, xena_dtypes, delete_raw_data=False):
    """Start a pipeline for importing data from GDC to Xena.

    Data will be imported on a dataset basis, which is defined by project
    and one specific data type. Data will be downloaded and transformed into
    the ``root_dir`` directory. Each project in projects will have its own
    directory under ``root_dir``. A project directory will have two
    subdirectories, "Raw_Data" and "Xena_Matrices", which hold downloaded
    data and transformed Xena matrices (with corresponding metadata)
    respectively. Errors raised during the importing of a dataset will be
    recorded in "etl.err" under ``root_dir`` and will not affect the importing
    of other datasets. Importing process, including the time consumption of
    the whole process, will be output to "stdout".

    Args:
        root_dir (str): An existing directory for saving imported data and log
            of errors.
        projects (list of str): A list of valid GDC project_id(s).
        xena_dtypes (list of str): A list of supported data type codes, which
            are (without quotation marks): "htseq_counts", "htseq_fpkm",
            "htseq_fpkm-uq", "mirna", "masked_cnv", "muse_snv", "mutect2_snv",
            "somaticsniper_snv", "varscan2_snv", "phenotype", "survival",
            "methylation27", "methylation450".
        delete_raw_data (optional, bool): Delete raw data upon generation
            of xena_matrix.
    """
    start_time = timeit.default_timer()
    counts = 0
    unfinished = {}
    total_projects = len(projects)
    log_format = '%(asctime)-15s [%(levelname)s]: %(message)s'
    logging.basicConfig(
        level=logging.WARNING,
        format=log_format,
        datefmt='%Y-%m-%d %H:%M:%S',
        filename=os.path.join(
            root_dir, 'etl_' + time.strftime("%Y%m%d-%H%M%S") + '.err',
        ),
        filemode='w',
    )
    logger = logging.getLogger('Xena-GDC-ETL')
    for project in projects:
        counts += 1
        msg = 'Importing [{}/{}] projects: {}'
        print(msg.format(counts, total_projects, project))
        for dtype in xena_dtypes:
            if dtype == 'survival':
                dataset = GDCSurvivalset(project, root_dir)
            elif dtype == 'raw_phenotype':
                if project.startswith('TCGA'):
                    dataset = GDCPhenoset(project, 'raw_phenotype', root_dir)
                if project.startswith('TARGET'):
                    dataset = GDCPhenoset(project, 'clinical', root_dir)
            elif dtype == 'GDC_phenotype':
                dataset = GDCPhenoset(project, 'GDC_phenotype', root_dir)
            elif dtype == 'Xena_phenotype':
                dataset = GDCAPIPhenoset(project, root_dir)
            else:
                dataset = GDCOmicset(project, dtype, root_dir)
            try:
                dataset.download().transform().metadata()
                if delete_raw_data:
                    print("Deleting raw data ...")
                    shutil.rmtree(dataset.raw_data_dir)
            except Exception:
                if project not in unfinished:
                    unfinished[project] = [dtype]
                else:
                    unfinished[project].append(dtype)
                with open(
                    os.path.join(root_dir, "unfinished.json"),
                    "w",
                ) as outfile:
                    json.dump(unfinished, outfile)
                msg = 'No {} data for cohort {}.'.format(dtype, project)
                logger.warning(msg, exc_info=True)
                print(msg)
    logging.shutdown()
    end_time = timeit.default_timer()
    m, s = divmod(int(end_time - start_time), 60)
    h, m = divmod(m, 60)
    print('Finish in {:d}:{:02d}:{:02d}.'.format(h, m, s))
