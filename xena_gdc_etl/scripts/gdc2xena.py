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

from xena_gdc_etl import gdc
from xena_gdc_etl.xena_dataset import GDCOmicset, GDCPhenoset, GDCSurvivalset


def gdc2xena(root_dir, projects, xena_dtypes):
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
    """
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
            elif dtype == 'raw_phenotype':
                if project.startswith('TCGA'):
                    dataset = GDCPhenoset(project, 'raw_phenotype', root_dir)
                if project.startswith('TARGET'):
                    dataset = GDCPhenoset(project, 'clinical', root_dir)
            elif dtype == 'GDC_phenotype':
                dataset = GDCPhenoset(project, 'GDC_phenotype', root_dir)
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
                   'somaticsniper_snv', 'varscan2_snv', 'raw_phenotype',
                   'GDC_phenotype', 'survival', 'methylation27',
                   'methylation450']
    parser = argparse.ArgumentParser(
            description='Pipeline for importing data from GDC to Xena.'
        )
    subparsers = parser.add_subparsers(title='Subcommands', dest='subcomm',
                                       metavar='')

    # Subcommand for full ETL (download, transform, and metadata)
    etlparser = subparsers.add_parser(
            'etl',
            help='Download and transform GDC data into Xena matrix, '
                 'and generate corresponding metadata.',
            epilog='Supported data types are: {}'.format(str(valid_dtype))
        )
    etlparser.add_argument('-r', '--root', type=str, default='.',
                           help='Root directory for imported data.')
    projects_group = etlparser.add_mutually_exclusive_group()
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
    datatype_group = etlparser.add_mutually_exclusive_group()
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

    # Subcommand for making metadata
    metaparser = subparsers.add_parser(
            'metadata',
            help='Generate metadata for a Xena matrix',
            epilog='Supported data types are: {}'.format(str(valid_dtype))
        )
    metaparser.add_argument('-p', '--project', type=str, required=True,
                            help='The project of the matrix.')
    metaparser.add_argument('-t', '--datatype', type=str, required=True,
                            help='One data type code for the matrix.')
    metaparser.add_argument('-m', '--matrix', type=str, required=True,
                            help='Path to a Xena matrix')
    metaparser.add_argument('-r', '--release', type=float, required=True,
                            help='GDC data release number.')

    args = parser.parse_args()
    if args.subcomm == 'etl':
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
    elif args.subcomm == 'metadata':
        root_dir = os.path.dirname(args.matrix)
        if args.datatype == 'survival':
            dataset = GDCSurvivalset(args.project, root_dir)
        elif args.datatype == 'raw_phenotype':
            if args.project.startswith('TCGA'):
                dataset = GDCPhenoset(args.project, 'raw_phenotype',
                                      root_dir)
            if args.project.startswith('TARGET'):
                dataset = GDCPhenoset(args.project, 'clinical', root_dir)
        elif args.datatype == 'GDC_phenotype':
            dataset = GDCPhenoset(args.project, 'GDC_phenotype', root_dir)
        else:
            dataset = GDCOmicset(args.project, args.datatype, root_dir)
        dataset.matrix = args.matrix
        dataset.gdc_release = (
                'https://docs.gdc.cancer.gov/Data/Release_Notes/Data_Release_Notes/#data-release-'
                + str(args.release).replace('.', '')
            )
        dataset.metadata()


if __name__ == '__main__':
    main()
