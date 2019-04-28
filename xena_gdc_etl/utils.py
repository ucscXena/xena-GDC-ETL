from __future__ import print_function
import os
import sys
import time

import pandas as pd
from pandas.util.testing import assert_frame_equal
import jinja2

from .constants import METADATA_TEMPLATE, METADATA_VARIABLES


def mkdir_p(dir_name):
    """Make the directory as needed: no error if existing.

    Args:
        dir_name (str): Directory name or path.

    Returns:
        str: The absolute path for the directory.
    """

    dir_path = os.path.abspath(dir_name)
    try:
        os.makedirs(dir_path)
    except OSError:
        if not os.path.isdir(dir_path):
            raise
    return dir_path


def equal_matrices(df1, df2):
    """This function is used for testing the equality of 2 Xena matrices

    Args:
        xena matrices df1 and df2 whose equality is to be checked
    """
    df1 = pd.read_table(df1, header=0,
                        index_col=0).sort_index(axis=0).sort_index(axis=1)
    df2 = pd.read_table(df2, header=0,
                        index_col=0).sort_index(axis=0).sort_index(axis=1)
    try:
        assert_frame_equal(df1, df2)
        print('Equal.')
    except:  # appeantly AssertionError doesn't catch all
        print('Not equal.')


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
