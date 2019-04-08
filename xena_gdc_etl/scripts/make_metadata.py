#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This script provides a command line tool for generating Xena metadata for a
Xena matrix.

Supported types of data are genomic data, phenotype data and survival data.
Please check help message with "-h" option for details.
"""

# Ensure Python 2 and 3 compatibility
from __future__ import print_function

import argparse
import os
import sys
import time

import jinja2

from ..constants import METADATA_TEMPLATE, METADATA_VARIABLES


def metadata(matrix, xena_dtypes):
    """Generating Xena metadata for a Xena matrix.

    Args:
        matrix (str): The path, including the filename, of the Xena matrix.
        xena_dtypes (str): One data type code indication the data type in
            matrices to be merged. Supported data type codes include (without
            quotation marks): "htseq_counts", "htseq_fpkm", "htseq_fpkm-uq",
            "mirna", "masked_cnv", "muse_snv", "mutect2_snv",
            "somaticsniper_snv", "varscan2_snv", "GDC_phenotype", "survival",
            "methylation27".
    """

    # Generate metadata.
    print('Creating metadata file ...', end='')
    sys.stdout.flush()
    jinja2_env = jinja2.Environment(
            loader=jinja2.PackageLoader('xena_gdc_etl', 'resources')
        )
    metadata_template = jinja2_env.get_template(METADATA_TEMPLATE[xena_dtypes])
    matrix_date = time.strftime("%m-%d-%Y",
                                time.gmtime(os.path.getmtime(matrix)))
    variables = {
            'project_id': 'GDC-PANCAN',
            'date': matrix_date,
            'gdc_release': 'https://docs.gdc.cancer.gov/Data/Release_Notes/Data_Release_Notes/#data-release-90', # noqa
            'xena_cohort': 'GDC Pan-Cancer (PANCAN)'
        }
    try:
        variables.update(METADATA_VARIABLES[xena_dtypes])
    except KeyError:
        pass
    outmetadata = matrix + '.json'
    with open(outmetadata, 'w') as f:
        f.write(metadata_template.render(**variables))
    print('\rMetadata JSON is saved at {}.'.format(outmetadata))


def main():
    valid_dtype = ['htseq_counts', 'htseq_fpkm', 'htseq_fpkm-uq', 'mirna',
                   'masked_cnv', 'muse_snv', 'mutect2_snv',
                   'somaticsniper_snv', 'varscan2_snv', 'GDC_phenotype',
                   'survival', 'methylation27', 'methylation450']
    parser = argparse.ArgumentParser(
            description='Generating Xena metadata for a Xena matrix.'
        )
    parser.add_argument('-m', '--matrix', type=str, required=True,
                        help='The path, including the filename, of the Xena '
                             'matrix.')
    parser.add_argument('-t', '--datatype', type=str, required=True,
                        help='One data type code indication the data type in '
                             'matrices to be merged. Supported data type '
                             'codes include: {}'.format(str(valid_dtype)))
    args = parser.parse_args()
    metadata(args.matrix, args.datatype)


if __name__ == '__main__':
    main()
