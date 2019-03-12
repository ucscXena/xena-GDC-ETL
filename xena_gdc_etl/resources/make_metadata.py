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

    # Map xena_dtype to corresponding metadata template.
    meta_templates = {'htseq_counts': 'template.rna.meta.json',
                      'htseq_fpkm': 'template.rna.meta.json',
                      'htseq_fpkm-uq': 'template.rna.meta.json',
                      'mirna': 'template.mirna.meta.json',
                      'mirna_isoform': 'template.mirna_isoform.meta.json',
                      'cnv': 'template.cnv.meta.json',
                      'masked_cnv': 'template.cnv.meta.json',
                      'muse_snv': 'template.snv.meta.json',
                      'mutect2_snv': 'template.snv.meta.json',
                      'somaticsniper_snv': 'template.snv.meta.json',
                      'varscan2_snv': 'template.snv.meta.json',
                      'GDC_phenotype': 'template.phenotype.meta.json',
                      'survival': 'template.survival.meta.json',
                      'methylation27': 'template.methylation.meta.json',
                      'methylation450': 'template.methylation.meta.json'}
    meta_templates_dir = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            'Resources')
    meta_templates = {
            k: os.path.join(meta_templates_dir, v)
            for k, v in meta_templates.items()
        }
    # Jinja2 template variables for corresponding "xena_dtype".
    meta_vars = {
            'htseq_counts': {'gdc_type': 'HTSeq - Counts',},
            'htseq_fpkm': {'gdc_type': 'HTSeq - FPKM',
                           'unit': 'fpkm'},
            'htseq_fpkm-uq': {'gdc_type': 'HTSeq - FPKM-UQ',
                              'unit': 'fpkm-uq'},
            'mirna': {'gdc_type': 'miRNA Expression Quantification'},
            'mirna_isoform': {'gdc_type': 'Isoform Expression Quantification'},
            'cnv': {'gdc_type': 'Copy Number Segment'},
            'masked_cnv': {'gdc_type': 'Masked Copy Number Segment'},
            'muse_snv': {'gdc_type': 'MuSE Variant Aggregation and Masking'},
            'mutect2_snv': {
                    'gdc_type': 'MuTect2 Variant Aggregation and Masking'
                },
            'somaticsniper_snv': {
                    'gdc_type': 'SomaticSniper Variant Aggregation and Masking'
                },
            'varscan2_snv': {
                    'gdc_type': 'VarScan2 Variant Aggregation and Masking'
                },
            'methylation27': {'platform_num': '27'},
            'methylation450': {'platform_num': '450'}
        }
    
    # Generate metadata.
    print('Creating metadata file ...', end='')
    sys.stdout.flush()
    template_json = meta_templates[xena_dtypes]
    file_dir = os.path.dirname(template_json)
    file_name = os.path.basename(template_json)
    jinja2_env = jinja2.Environment(
            loader=jinja2.FileSystemLoader(file_dir)
        )
    metadata_template = jinja2_env.get_template(file_name)
    matrix_date = time.strftime("%m-%d-%Y",
                                time.gmtime(os.path.getmtime(matrix)))
    variables = {
            'project_id': 'GDC-PANCAN',
            'date': matrix_date,
            'gdc_release': 'https://docs.gdc.cancer.gov/Data/Release_Notes/Data_Release_Notes/#data-release-90',
            'xena_cohort': 'GDC Pan-Cancer (PANCAN)'
        }
    try:
        variables.update(meta_vars[xena_dtypes])
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
