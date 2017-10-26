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

import argparse
import logging
import os
import timeit

import gdc
from xena_dataset import (GDCOmicset,
                          TCGAPhenoset, TARGETPhenoset,
                          GDCSurvivalset)

def gdc2xena(root_dir, projects, xena_dtypes):
    start_time = timeit.default_timer()
    
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
        msg = 'Importing [{}/{}] projects: {}'
        print(msg.format(counts, total_projects, project))
        for dtype in xena_dtypes:
            if dtype == 'survival':
                dataset = GDCSurvivalset(project, root_dir)
            elif dtype == 'phenotype':
                if project.startswith('TCGA'):
                    dataset = TCGAPhenoset(project, root_dir)
                if project.startswith('TARGET'):
                    dataset = TARGETPhenoset(project, root_dir)
            else:
                dataset = GDCOmicset(project, dtype, root_dir)
            try:
                dataset.download().transform().metadata()
            except Exception:
                msg = 'No {} data for cohort {}.'.format(dtype, project)
                logger.warning(msg, exc_info=True)
                print(msg)
    logging.shutdown()
    
    end_time = timeit.default_timer()
    m, s = divmod(int(end_time - start_time), 60)
    h, m = divmod(m, 60)
    print('Finish in {:d}:{:02d}:{:02d}.'.format(h, m, s))


def main():
    valid_dtype = ['htseq_counts', 'htseq_fpkm', 'htseq_fpkm-uq', 'mirna',
                   'masked_cnv', 'muse_snv', 'mutect2_snv', 
                   'somaticsniper_snv', 'varscan2_snv', 'phenotype',
                   'survival']
    parser = argparse.ArgumentParser(
            description='Pipeline for importing data from GDC to Xena.',
            epilog='Supported data types are: {}'.format(str(valid_dtype))
        )
    parser.add_argument('-r', '--root', type=str,
                        help='Root directory for imported data.', default='.')
    projects_group = parser.add_mutually_exclusive_group()
    projects_group.add_argument('-p', '--projects', type=str, nargs='+',
                                help='GDC project ID(s) to be imported; or '
                                     '"all" if all projects on GDC are going '
                                     'to be imported. Defaults to "all".',
                                default=['all'])
    projects_group.add_argument('-P', '--not-projects', type=str, nargs='+',
                                help='Import all projects on GDC except '
                                     'projects specified by this option. '
                                     'This option and the "-p" option are '
                                     'mutually exclusive.',
                                default=[])
    datatype_group = parser.add_mutually_exclusive_group()
    datatype_group.add_argument('-t', '--datatype', type=str, nargs='+',
                                help='Data type code(s) to be imported; or '
                                     '"all" if all supported types are going '
                                     'to be imported. Defaults to "all".',
                                default=['all'])
    datatype_group.add_argument('-T', '--not-datatype', type=str, nargs='+',
                                help='Import all supported types except '
                                     'projects specified by this option. '
                                     'This option and the "-t" option are '
                                     'mutually exclusive.',
                                default=[])
    args = parser.parse_args()
    root_dir = os.path.abspath(args.root)
    projects = args.projects
    if 'all' in [p.lower() for p in projects]:
        projects = [str(x) for x in gdc.get_project_info().index]
    for p in args.not_projects:
        projects.remove(p)
    xena_dtypes = args.datatype
    if 'all' in [t.lower() for t in xena_dtypes]:
        xena_dtypes = valid_dtype
    for t in args.not_datatype:
        xena_dtypes.remove(t)
    print('#### GDC to Xena Importing Settings ####')
    total_projects = len(projects)
    print('Import the following {} projects:'.format(total_projects))
    print(repr(projects), end='\n\n')
    print('for the following {} types of data:'.format(len(xena_dtypes)))
    print(str(xena_dtypes), end='\n\n')
    print('into this directory: {}'.format(root_dir))
    print('########################################', end='\n\n')
    gdc2xena(root_dir, projects, xena_dtypes)


if __name__ == '__main__':
    main()
