#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This script merges, for each type of data, matrices from all TCGA projects
and generate a pan-cancer matrix of corresponding datatype for GDC TCGA
Pan-Cancer (PANCAN) cohort. It was written specifically for the directory
structure and file nomenclature used by ``gdc2xena.py`` ETL process.
In general, this script should be used right after importing/updating any
individual TCGA project.

To use/customize this script, change the assignment of ``root_dir``,
``datatypes`` and ``gdc_release`` accordingly. ``root_dir`` is the parent
directory for all TCGA projects to be merged. Merging process will iterate
through the list of ``datatype`` and generates merged pan-cancer matrix one by
one. Supported types of data are ['htseq_counts', 'htseq_fpkm',
'htseq_fpkm-uq', 'mirna', 'masked_cnv', 'muse_snv', 'mutect2_snv',
'somaticsniper_snv', 'varscan2_snv', 'raw_phenotype', 'GDC_phenotype',
'survival']. ``gdc_release`` is the URL to GDC's data release note of the data
included in this matrix.
"""

# Ensure Python 2 and 3 compatibility
from __future__ import division
from __future__ import print_function

import glob
import os
import sys
import time

import jinja2
import pandas as pd

from ..constants import METADATA_TEMPLATE, METADATA_VARIABLES  # noqa


def main():
    root_dir = r'/mnt/gdc/updates'
    out_dir = r'/mnt/gdc/updates/GDC-PANCAN/Xena_Matrices'
    datatypes = ['htseq_counts', 'htseq_fpkm', 'htseq_fpkm-uq', 'mirna',
                 'masked_cnv', 'muse_snv', 'mutect2_snv', 'somaticsniper_snv',
                 'varscan2_snv', 'survival']
    gdc_release = 'https://docs.gdc.cancer.gov/Data/Release_Notes/Data_Release_Notes/#data-release-100'
    meta_templates_dir = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            'Resources')
    meta_templates = {
            k: os.path.join(meta_templates_dir, v)
            for k, v in METADATA_TEMPLATE.items()
        }
    
    for dtype in datatypes:
        if dtype in ['htseq_counts', 'htseq_fpkm', 'htseq_fpkm-uq', 'mirna']:
            merge_axis = 1
        elif dtype in ['masked_cnv', 'muse_snv', 'mutect2_snv',
                       'somaticsniper_snv', 'varscan2_snv',
                       'raw_phenotype', 'GDC_phenotype', 'survival']:
            merge_axis = 0
        else:
            msg = 'Invalid datatype: {}\nSupported data types are: {}'
            valid_dtype = ['htseq_counts', 'htseq_fpkm', 'htseq_fpkm-uq',
                           'mirna', 'masked_cnv', 'muse_snv', 'mutect2_snv',
                           'somaticsniper_snv', 'varscan2_snv',
                           'raw_phenotype', 'GDC_phenotype', 'survival']
            raise ValueError(msg.format(dtype, valid_dtype))
            return
        print('\n########################################')
        # Gather the list of matrices to be merged
        pathpattern = os.path.join(root_dir, 'TCGA-*', 'Xena_Matrices',
                                   '*.{}.tsv'.format(dtype))
        matrices = []
        for path in glob.glob(pathpattern):
            print('\rReading {} ...'.format(path), end='')
            sys.stdout.flush()
            matrices.append(pd.read_table(path, header=0, index_col=0))
        # Merge matrices
        print('\rMerging {} "{}" matrices ...'.format(len(matrices), dtype))
        merged = pd.concat(matrices, axis=merge_axis)
        outmatrix = os.path.join(out_dir, 'GDC-PANCAN.{}.tsv'.format(dtype))
        print('Saving merged matrix to {} ...'.format(outmatrix), end='')
        sys.stdout.flush()
        merged.to_csv(outmatrix, sep='\t', encoding='utf-8')
        del merged # Prevent from doubling memory usage
        print('\rMerged "{}" matrix is ready at {}'.format(dtype, outmatrix))
        # Generate metadata.
        print('Creating metadata file ...', end='')
        sys.stdout.flush()
        template_json = meta_templates[dtype]
        file_dir = os.path.dirname(template_json)
        file_name = os.path.basename(template_json)
        jinja2_env = jinja2.Environment(
                loader=jinja2.FileSystemLoader(file_dir)
            )
        metadata_template = jinja2_env.get_template(file_name)
        variables = {
                'project_id': 'GDC-PANCAN',
                'date': time.strftime(
                        "%m-%d-%Y", time.gmtime(os.path.getmtime(outmatrix))
                    ),
                'gdc_release': gdc_release,
                'xena_cohort': 'GDC Pan-Cancer (PANCAN)'
            }
        try:
            variables.update(METADATA_VARIABLES[dtype])
        except KeyError:
            pass
        outmetadata = outmatrix + '.json'
        with open(outmetadata, 'w') as f:
            f.write(metadata_template.render(**variables))
        print('\rMetadata JSON is saved at {}.'.format(outmetadata))


if __name__ == '__main__':
    main()
