#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This script provides a command line tool for merging Xena matrices of the
same data type.

This module also provides an stand-alone function ``merge`` which can be used
in scripts. For example::

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
from datetime import date
import glob
import os
import time
import timeit
import sys

import jinja2
import pandas as pd

def merge(filelist, xena_dtypes, out_matrix):
    """Merge a list (``filelist``) of Xena matrices of the same data type
    (``xena_dtypes``) into a single Xena matrix (``out_matrix``).
    
    Args:
        filelist (list of path): A list of Xena matrices files to be merged,
            which will be read by pandas.read_table.
        xena_dtypes (str): One data type code indication the data type in
            matrices to be merged. Supported data type codes include (without
            quotation marks): "htseq_counts", "htseq_fpkm", "htseq_fpkm-uq",
            "mirna", "masked_cnv", "muse_snv", "mutect2_snv",
            "somaticsniper_snv", "varscan2_snv", "GDC_phenotype", "survival",
            "methylation27".
        out_matrix (str): The path, including the filename, for merged Xena
            matrix.
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
                      'survival': 'template.survival.meta.json'}
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
                }
        }
    gdc_release = 'https://docs.gdc.cancer.gov/Data/Release_Notes/Data_Release_Notes/#data-release-90'
    
    start_time = timeit.default_timer()
    
    if xena_dtypes in ['htseq_counts', 'htseq_fpkm', 'htseq_fpkm-uq',
                       'mirna', 'methylation27']:
        merge_axis = 1
    elif xena_dtypes in ['masked_cnv', 'muse_snv', 'mutect2_snv',
                         'somaticsniper_snv', 'varscan2_snv',
                         'GDC_phenotype', 'survival']:
        merge_axis = 0
    else:
        msg = 'Invalid datatype: {}\nSupported data types are: {}'
        valid_dtype = ['htseq_counts', 'htseq_fpkm', 'htseq_fpkm-uq',
                       'mirna', 'masked_cnv', 'muse_snv', 'mutect2_snv',
                       'somaticsniper_snv', 'varscan2_snv', 'GDC_phenotype',
                       'survival', 'methylation27']
        raise ValueError(msg.format(xena_dtypes, valid_dtype))
        return
    # Merge by growing merged matrix one matrix at a time.
    merged = pd.DataFrame()
    count = 0
    total = len(filelist)
    for path in filelist:
        count += 1
        print('\r[{}/{}] Merging {} ...'.format(count, total, path), end='')
        sys.stdout.flush()
        merged = pd.concat(
                [merged, pd.read_table(path, header=0, index_col=0)],
                axis=merge_axis, copy=False
            )
    print('\rSaving merged matrix to {} ...'.format(out_matrix), end='')
    sys.stdout.flush()
    merged.to_csv(out_matrix, sep='\t', encoding='utf-8')
    del merged # Prevent from doubling memory usage
    print('\rMerged "{}" matrix is ready at {}'.format(xena_dtypes,
                                                       out_matrix))
    return
    



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
    matrix_date = time.strftime(
                    "%m-%d-%Y", time.gmtime(os.path.getmtime(out_matrix))
                )
    variables = {'date': matrix_date,
                 'gdc_release': gdc_release,
                 'xena_cohort': 'GDC Pan-Cancer (PANCAN)'}
    try:
        variables.update(meta_vars[xena_dtypes])
    except KeyError:
        pass
    outmetadata = out_matrix + '.json'
    with open(outmetadata, 'w') as f:
        f.write(metadata_template.render(**variables))
    print('\rMetadata JSON is saved at {}.'.format(outmetadata))



    end_time = timeit.default_timer()
    m, s = divmod(int(end_time - start_time), 60)
    h, m = divmod(m, 60)
    print('Finish in {:d}:{:02d}:{:02d}.'.format(h, m, s))


def main():
    valid_dtype = ['htseq_counts', 'htseq_fpkm', 'htseq_fpkm-uq', 'mirna',
                   'masked_cnv', 'muse_snv', 'mutect2_snv',
                   'somaticsniper_snv', 'varscan2_snv', 'GDC_phenotype',
                   'survival', 'methylation27']
    parser = argparse.ArgumentParser(
            description='Pipeline for merging Xena matrices of the same data '
                        'type.'
        )
    parser.add_argument('-f', '--files', type=str, nargs='+', required=True,
                        help='A list of paths for Xena matrices files to be '
                             'merged. All paths in this list support UNIX '
                             'style pathname pattern expansion with "glob". '
                             'Files will be read by pandas.read_table.')
    parser.add_argument('-t', '--datatype', type=str, required=True,
                        help='One data type code indication the data type in '
                             'matrices to be merged. Supported data type '
                             'codes include: {}'.format(str(valid_dtype)))
    parser.add_argument('-o', '--outdir', type=str, default='.',
                        help='A directory to put the merged matrix. Defaults '
                             'to the current working directory of python.')
    parser.add_argument('-n', '--name', type=str, default=None,
                        help='Filename for the merged matrix. Defaults to '
                             'None. If None, the filename will be derived '
                             'from the cohort name and the data type. Check '
                             '"-t" and "-c" options for details.')
    parser.add_argument('-c', '--cohort', type=str, default=None,
                        help='A cohort name for the merged matrix. Defaults '
                             'to None. If None, it will be set to a format of '
                             '"MergedCohort<date>" by default. For example, '
                             '"MergedCohort{}".'
                             ''.format(date.today().strftime('%m%d%Y')))
    args = parser.parse_args()
    print('Checking matrices to be merged ...')
    matrix_list = []
    for path in args.files:
        for f in glob.glob(path):
            f = os.path.abspath(f)
            if os.path.isfile(f):
                print('\r{}'.format(f), end='')
                matrix_list.append(f)
    print('\r{} matrices to be merged'.format(len(matrix_list)))
    if args.name is None:
        if args.cohort is None:
            cohort = 'MergedCohort{}'.format(date.today().strftime('%m%d%Y'))
        else:
            cohort = args.cohort
        matrix_name = '{}.{}.tsv'.format(cohort, args.datatype)
    else:
        matrix_name = args.name
    matrix = os.path.join(os.path.abspath(args.outdir), matrix_name)
    merge(matrix_list, args.datatype, matrix)


if __name__ == '__main__':
    main()
