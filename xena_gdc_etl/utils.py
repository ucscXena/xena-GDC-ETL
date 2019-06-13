from __future__ import print_function
from datetime import date
import glob
import os
import sys
import time
import timeit

import pandas as pd
import jinja2
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from .constants import (
    METADATA_TEMPLATE,
    METADATA_VARIABLES,
    valid_dtype,
    GDC_RELEASE_URL,
)


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


def handle_merge_xena(name, files, cohort, datatype, outdir):
    print('Checking matrices to be merged ...')
    matrix_list = []
    for path in files:
        for f in glob.glob(path):
            f = os.path.abspath(f)
            if os.path.isfile(f):
                print('\r{}'.format(f), end='')
                matrix_list.append(f)
    print('\r{} matrices to be merged'.format(len(matrix_list)))
    if name is None:
        if cohort is None:
            cohort = 'MergedCohort{}'.format(date.today().strftime('%m%d%Y'))
        else:
            cohort = cohort
        matrix_name = '{}.{}.tsv'.format(cohort, datatype)
    else:
        matrix_name = name
    matrix = os.path.join(os.path.abspath(outdir), matrix_name)
    merge(matrix_list, datatype, matrix)


def merge(filelist, xena_dtypes, out_matrix):
    """Merge a list (``filelist``) of Xena matrices of the same data type
    (``xena_dtypes``) into a single Xena matrix (``out_matrix``).

    Args:
        filelist (list of path): A list of Xena matrices files to be merged,
            which will be read by pandas.read_csv.
        xena_dtypes (str): One data type code indication the data type in
            matrices to be merged. Supported data type codes include (without
            quotation marks): "htseq_counts", "htseq_fpkm", "htseq_fpkm-uq",
            "mirna", "masked_cnv", "muse_snv", "mutect2_snv",
            "somaticsniper_snv", "varscan2_snv", "GDC_phenotype", "survival",
            "methylation27".
        out_matrix (str): The path, including the filename, for merged Xena
            matrix.
    """

    meta_templates_dir = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        'Resources',
    )
    meta_templates = {
        k: os.path.join(meta_templates_dir, v)
        for k, v in METADATA_TEMPLATE.items()
    }
    gdc_release = GDC_RELEASE_URL + '#data-release-90'
    start_time = timeit.default_timer()
    if xena_dtypes in [
        'htseq_counts',
        'htseq_fpkm',
        'htseq_fpkm-uq',
        'mirna',
        'methylation27',
    ]:
        merge_axis = 1
    elif xena_dtypes in [
        'masked_cnv',
        'muse_snv',
        'mutect2_snv',
        'somaticsniper_snv',
        'varscan2_snv',
        'GDC_phenotype',
        'survival',
    ]:
        merge_axis = 0
    else:
        msg = 'Invalid datatype: {}\nSupported data types are: {}'
        raise ValueError(msg.format(xena_dtypes, valid_dtype))
    # Merge by growing merged matrix one matrix at a time.
    merged = pd.DataFrame()
    count = 0
    total = len(filelist)
    for path in filelist:
        count += 1
        print('\r[{}/{}] Merging {} ...'.format(count, total, path), end='')
        sys.stdout.flush()
        merged = pd.concat(
            [merged, pd.read_csv(path, sep='\t', header=0, index_col=0)],
            axis=merge_axis,
            copy=False,
        )
    print('\rSaving merged matrix to {} ...'.format(out_matrix), end='')
    sys.stdout.flush()
    merged.to_csv(out_matrix, sep='\t', encoding='utf-8')
    del merged  # Prevent from doubling memory usage
    print(
        '\rMerged "{}" matrix is ready at {}'.format(xena_dtypes, out_matrix)
    )
    return
    # Generate metadata.
    print('Creating metadata file ...', end='')
    sys.stdout.flush()
    template_json = meta_templates[xena_dtypes]
    file_dir = os.path.dirname(template_json)
    file_name = os.path.basename(template_json)
    jinja2_env = jinja2.Environment(loader=jinja2.FileSystemLoader(file_dir))
    metadata_template = jinja2_env.get_template(file_name)
    matrix_date = time.strftime(
        "%m-%d-%Y", time.gmtime(os.path.getmtime(out_matrix))
    )
    variables = {
        'date': matrix_date,
        'gdc_release': gdc_release,
        'xena_cohort': 'GDC Pan-Cancer (PANCAN)',
    }
    try:
        variables.update(METADATA_VARIABLES[xena_dtypes])
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


def reduce_json_array(j):
    """Recursively go over a JSON and unpack arrays which have only one item,
    i.e. remove unnecessary arrays (brackets).

    Args:
        j (list of dict): A JSON to be reduced.

    Returns:
        list or dict: a reduced JSON with unnecessary array removed.
    """

    if isinstance(j, list):
        if len(j) == 1:
            reduced = reduce_json_array(j[0])
        else:
            reduced = [reduce_json_array(e) for e in j]
    elif isinstance(j, dict):
        reduced = {k: reduce_json_array(v) for k, v in j.items()}
    else:
        reduced = j
    return reduced


def requests_retry_session(
    retries=5,
    backoff_factor=0.3,
    status_forcelist=(500, 502, 504),
    session=None,
):
    session = session or requests.Session()
    retry = Retry(
        total=retries,
        read=retries,
        connect=retries,
        backoff_factor=backoff_factor,
        status_forcelist=status_forcelist,
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount('http://', adapter)
    session.mount('https://', adapter)
    return session
