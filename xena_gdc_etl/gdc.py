#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This module provides basic and minimum necessary functions for carrying out
data query and download for Xena GDC ETL pipelines.
"""

# Ensure Python 2 and 3 compatibility
from __future__ import division
from __future__ import print_function

import json
import os
import sys
import warnings

import pandas as pd
import requests

from .utils import (
    mkdir_p,
    reduce_json_array,
    extract_leaves,
    get_intersection,
)

GDC_API_BASE = 'https://api.gdc.cancer.gov'
_SUPPORTED_FILE_TYPES = {
    'txt',
    'vcf',
    'bam',
    'tsv',
    'xml',
    'maf',
    'xlsx',
    'tar',
    'gz',
    'md5',
    'xls',
}
_SUPPORTED_DATASETS = [
    {'data_type': 'Copy Number Segment'},
    {'data_type': 'Masked Copy Number Segment'},
    {'data_type': 'Isoform Expression Quantification'},
    {'data_type': 'miRNA Expression Quantification'},
    {'data_type': 'Methylation Beta Value'},
    {'analysis.workflow_type': 'HTSeq - Counts'},
    {'analysis.workflow_type': 'HTSeq - FPKM'},
    {'analysis.workflow_type': 'HTSeq - FPKM-UQ'},
    {'analysis.workflow_type': 'MuSE Variant Aggregation and Masking'},
    {'analysis.workflow_type': 'MuTect2 Variant Aggregation and Masking'},
    {
        'analysis.workflow_type':
        'SomaticSniper Variant Aggregation and Masking'
    },
    {'analysis.workflow_type': 'VarScan2 Variant Aggregation and Masking'},
    {'data_type': 'Biospecimen Supplement'},
    {'data_type': 'Clinical Supplement'},
]
# https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations
TCGA_STUDY_ABBR = {
    'LAML': 'Acute Myeloid Leukemia',
    'ACC': 'Adrenocortical carcinoma',
    'BLCA': 'Bladder Urothelial Carcinoma',
    'LGG': 'Brain Lower Grade Glioma',
    'BRCA': 'Breast invasive carcinoma',
    'CESC': (
        'Cervical squamous cell carcinoma and endocervical adenocarcinoma'
    ),
    'CHOL': 'Cholangiocarcinoma',
    'LCML': 'Chronic Myelogenous Leukemia',
    'COAD': 'Colon adenocarcinoma',
    'CNTL': 'Controls',
    'ESCA': 'Esophageal carcinoma',
    'FPPP': 'FFPE Pilot Phase II',
    'GBM': 'Glioblastoma multiforme',
    'HNSC': 'Head and Neck squamous cell carcinoma',
    'KICH': 'Kidney Chromophobe',
    'KIRC': 'Kidney renal clear cell carcinoma',
    'KIRP': 'Kidney renal papillary cell carcinoma',
    'LIHC': 'Liver hepatocellular carcinoma',
    'LUAD': 'Lung adenocarcinoma',
    'LUSC': 'Lung squamous cell carcinoma',
    'DLBC': 'Lymphoid Neoplasm Diffuse Large B-cell Lymphoma',
    'MESO': 'Mesothelioma',
    'MISC': 'Miscellaneous',
    'OV': 'Ovarian serous cystadenocarcinoma',
    'PAAD': 'Pancreatic adenocarcinoma',
    'PCPG': 'Pheochromocytoma and Paraganglioma',
    'PRAD': 'Prostate adenocarcinoma',
    'READ': 'Rectum adenocarcinoma',
    'SARC': 'Sarcoma',
    'SKCM': 'Skin Cutaneous Melanoma',
    'STAD': 'Stomach adenocarcinoma',
    'TGCT': 'Testicular Germ Cell Tumors',
    'THYM': 'Thymoma',
    'THCA': 'Thyroid carcinoma',
    'UCS': 'Uterine Carcinosarcoma',
    'UCEC': 'Uterine Corpus Endometrial Carcinoma',
    'UVM': 'Uveal Melanoma',
}


def simple_and_filter(in_dict={}, exclude_dict={}):
    """Make a simple GDC API compatible query filter from a dict, in which
    individual conditions are joint by the "and" logic.

    In the return filter, individual conditions in the ``in_dict`` and
    ``exclude_dict`` will be joint by the "and" operator, meaning a hit has to
    match all conditions. Here, a condition can use either a "in" operator
    (specified in the ``in_dict``) or a "exclude" operator (specified in the
    ``exclude_dict``).
    See details at
    https://docs.gdc.cancer.gov/API/Users_Guide/Search_and_Retrieval/#filters-specifying-the-query

    Args:
        in_dict (dict): A dict describing query conditions with the "in"
            operator. Each (key, value) pair represents for one condition. The
            "key" is the 'field' operand. Operator between "key" and "value"
            is "in".
        exclude_dict (dict): A dict describing query conditions with the
            "exclude" operator. Each (key, value) pair represents for one
            condition. The "key" is the 'field' operand. Operator between
            "key" and "value" is "exclude_dict".
    Returns:
        dict: A dict of filter conforming to GDC API's format. It should then
        be converted to JSON format and used in the following http request.
    """

    if not in_dict and not exclude_dict:
        return in_dict
    operation_list = []
    for key in in_dict:
        value = in_dict[key]
        if not isinstance(value, list):
            value = [value]
        operation_list.append(
            {"op": "in", "content": {"field": key, "value": value}}
        )
    for key in exclude_dict:
        value = exclude_dict[key]
        if not isinstance(value, list):
            value = [value]
        operation_list.append(
            {"op": "exclude", "content": {"field": key, "value": value}}
        )
    return {"op": "and", "content": operation_list}


def search(
    endpoint,
    in_filter={},
    exclude_filter={},
    fields=[],
    expand=[],
    typ='dataframe',
    method='GET',
):
    """Search one GDC endpoints and return searching results in a pandas
    DataFrame if possible.

    When searching results cannot be safely converted to a pandas DataFrame,
    results will be returned in the JSON format as it is returned from GDC
    API.

    Args:
        endpoint (str): One string of GDC API supported endpoint. See:
            https://docs.gdc.cancer.gov/API/Users_Guide/Getting_Started/#api-endpoints
        in_filter (dict, optional): A dict of query conditions which will be
            used to perform the query. Each (key, value) pair represents for
            one ondition. It will be passed to ``simple_and_filter`` for
            making a query filter compatible with GDC API. Please check
            ``simple_and_filter`` function for details.
        exclude_filter (dict, optional): An optional dict of query conditions
            which will be used to perform the query. Each (key, value) pair
            represents for one condition. It will be passed to
            ``simple_and_filter`` for making a query filter compatible with
            GDC API. Please check ``simple_and_filter`` function for details.
        fields (list or str, optional): One or more fields to be queried. Each
            field will be used as a column name in the returned DataFrame. It
            can be a comma separated string or a list of field strings or a
            combination of both.
        expand (list or str, optional): One or more field groups to be
            queried. It can be a comma separated string or a list of field
            strings or a combination of both.
        typ (str): type of search result to return (JSON or dataframe).
            Defaults to 'dataframe'.
        method (str): HTTP method for the search. Defaults to 'GET'.
            ..

    Returns:
        pandas.core.frame.DataFrame or str: A search result in form of a
            pandas DataFrame or a JSON formatted string, depending on the
            value of ``typ`` and the DataFrame convertibility of JSON.
    """

    try:
        assert typ.lower() in ['json', 'dataframe']
    except (AttributeError, AssertionError):
        raise ValueError(
            'typ should be a string of either JSON or dataframe, '
            'not {}'.format(typ)
        )
    filters = simple_and_filter(in_dict=in_filter, exclude_dict=exclude_filter)
    if isinstance(fields, str):
        fields = [fields]
    if isinstance(expand, str):
        expand = [expand]
    payload = {'size': 1}
    if filters:
        payload['filters'] = json.dumps(filters)
    if fields:
        payload['fields'] = ','.join(fields)
    if expand:
        payload['expand'] = ','.join(expand)
    url = '{}/{}'.format(GDC_API_BASE, endpoint)
    if method.upper() == 'POST':
        response = requests.post(url, data=payload)
    elif method.upper() == 'GET':
        response = requests.get(url, params=payload)
    else:
        raise ValueError(
            'Invalid method: {}\n method must be either "GET" '
            'or "POST".'.format(method)
        )
    try:
        payload['size'] = response.json()['data']['pagination']['total']
    except KeyError:
        payload.pop('size')
        response = requests.get(url, params=payload)
        if typ.lower() == 'json':
            return response.json()
        else:
            warnings.warn(
                'Fail to get a table of results. JSON returned. '
                'Please check the result carefully.',
                stacklevel=2,
            )
            return response.json()
    if method.upper() == 'POST':
        response = requests.post(url, data=payload)
    else:
        response = requests.get(url, params=payload)
    if response.status_code == 200:
        results = response.json()['data']['hits']
        if typ.lower() == 'json':
            return results
        try:
            return pd.io.json.json_normalize(reduce_json_array(results))
        except Exception:
            warnings.warn(
                'Fail to convert searching results into table. '
                'JSON will be returned.',
                stacklevel=2,
            )
            return results
    else:
        warnings.warn(
            'Searching failed with HTTP status code: '
            '{}'.format(response.status_code),
            stacklevel=2,
        )
        return None


def get_ext(file_name):
    """Get all extensions supported by this module in the file name.

    Supported extensions are defined in the constant "_SUPPORTED_FILE_TYPES".
    Multiple extensions will be separated by ".".

    Args:
        file_name (str): The filename will be split by "." and checked from
            left to right. Extensions will be kept starting from the first
            (and left most) supported extension.

    Returns:
        str: A string of extensions joint by ".".
    """

    # https://github.com/broadinstitute/gdctools/blob/master/gdctools/lib/meta.py
    name_list = file_name.split('.')
    for i in range(len(name_list)):
        if name_list[i] in _SUPPORTED_FILE_TYPES:
            break
    return '.'.join(name_list[i:])


def download(uuids, download_dir='.', chunk_size=4096):
    """Download GDC's open access data according to UUID(s).

    Args:
        uuids (str, list or dict): A single UUID (str), a list of UUIDs (list)
            or a dict whose keys are UUIDs for target file(s). If "uuids" is
            str or list, downloaded file(s) will be renamed to
            "UUID.extension" where "extension" is extracted by "get_ext()"
            from the original filename. Renamed file(s) will be saved at
            "download_dir". If "uuids" is a dict, the argument "download_dir"
            will be ignored; values of dict will be paths for saving
            corresponding downloaded files.
        download_dir (str, optional): The directory for saving downloaded
            file(s) when "uuids" is str or list. It will be ignored if "uuids"
            is dict. Defaults to ".".
        chunk_size (int, optional): The chunk size is the number of bytes it
            should read into memory, when the response is got with
            "stream=True". Check the documentation of "requests" module for
            details. Defaults to 4096.

    Returns:
        list: a list of paths for downloaded files.
    """

    if isinstance(uuids, str):
        uuids = {uuids: None}
    elif isinstance(uuids, list):
        uuids = {uuid: None for uuid in uuids}
    elif not isinstance(uuids, dict):
        raise TypeError(
            'uuids is a {}; it should be a string, a list or a '
            'dict'.format(type(uuids))
        )
    total = len(uuids)
    count = 0
    download_list = []
    data_endpt = '{}/data/'.format(GDC_API_BASE)
    for uuid in uuids:
        count += 1
        response = requests.get(data_endpt + uuid, stream=True)
        if response.status_code == 200:
            file_size = int(response.headers['Content-Length'])
            if uuids[uuid] is None:
                content_disp = response.headers['Content-Disposition']
                ori_name = content_disp[content_disp.find('filename=') + 9 :]  # noqa: E203, E501
                new_filename = uuid + '.' + get_ext(ori_name)
                path = os.path.join(
                    os.path.abspath(download_dir), new_filename
                )
            else:
                path = os.path.abspath(uuids[uuid])
            status = '\r[{:d}/{:d}] Download to "{}": {:4.0%}'
            mkdir_p(os.path.dirname(path))
            with open(path, 'wb') as f:
                downloaded = 0
                print(status.format(count, total, path, 0), end='')
                sys.stdout.flush()
                for chunk in response.iter_content(chunk_size):
                    f.write(chunk)
                    downloaded = downloaded + chunk_size
                    print(
                        status.format(
                            count, total, path, min(1, downloaded / file_size)
                        ),
                        end='',
                    )
                    sys.stdout.flush()
            download_list.append(path)
        else:
            print('\rFail to download file {}.'.format(uuid))
    print('')
    return download_list


def get_project_info(projects=None):
    """Get info for project(s) of interest through GDC API.

    Args:
        projects (list or str): one (str) or a list of GDC "project_id"(s),
            whose info will be returned. If None, projects will not be
            filtered, i.e. info for all GDC projects will be returned.
            Defaults to None.

    Returns:
        pandas.core.frame.DataFrame: A DataFrame of project info including
        "project ID", "project name", "primary site" and "program name".
    """

    in_filter = {}
    if projects is not None:
        if isinstance(projects, list):
            in_filter = {'projects.project_id': projects}
        else:
            in_filter = {'projects.project_id': [projects]}
    project_df = search(
        'projects',
        in_filter=in_filter,
        fields=['name', 'primary_site', 'project_id', 'program.name'],
    )
    return project_df.set_index('id')


def get_samples_clinical(projects=None):
    """Get info for all samples of ``projects`` and clinical info for all
    cases of ``projects`` through GDC API.

    Args:
        projects (list or str): one (str) or a list of GDC "project_id"(s),
            whose info will be returned. If None, projects will not be
            filtered, i.e. info for all GDC projects will be returned.
            Defaults to None.

    Returns:
        pandas.core.frame.DataFrame: A DataFrame organized by samples, having
        info for all samples of ``projects``, as well as corresponding
        clinical info.
    """

    in_filter = {}
    if projects is not None:
        if isinstance(projects, list):
            in_filter = {'project.project_id': projects}
        else:
            in_filter = {'project.project_id': [projects]}
    fields = [
        'case_id',
        'created_datetime',
        'disease_type',
        'id',
        'primary_site',
        'state',
        'submitter_id',
        'updated_datetime',
    ]
    expand = [
        'demographic',
        'diagnoses',
        'diagnoses.treatments',
        'exposures',
        'family_histories',
        'project',
        'samples',
        'tissue_source_site',
    ]
    res = search(
        'cases', in_filter=in_filter, fields=fields, expand=expand, typ='json'
    )
    reduced_no_samples_json = reduce_json_array(
        [{k: v for k, v in d.items() if k != 'samples'} for d in res]
    )
    cases_df = pd.io.json.json_normalize(reduced_no_samples_json)
    # In the list of reduced json, "samples" fields for each case are not
    # consistently ``list`` (if there is only 1 sample for the case, it will
    # be reduced into "naked" ``dict``). Therefore, it cannot be normalized
    # correctly with ``record_path `` "samples". Use the raw json instead.
    # Besides, there are cases (34 by 12/11/2017) which doesn't have any
    # samples and thus doesn't have the key "samples". Ignore them.
    #    for r in res:
    #        r.setdefault('samples', [{}])
    #        samples_json.append(r)
    samples_df = pd.io.json.json_normalize(
        [r for r in res if 'samples' in r],
        'samples',
        'id',
        record_prefix='samples.',
    )
    return pd.merge(cases_df, samples_df, how='inner', on='id')


def gdc_check_new(new_file_uuids):
    """
    This function help check a list of GDC's updated files and summarize
    impacted project(s), data_type(s) and analysis.workflow_type(s).
    """

    df_list = []
    for uuids in (
        new_file_uuids[i : i + 20000]  # noqa: E203
        for i in range(0, len(new_file_uuids), 20000)
    ):
        df = search(
            'files',
            in_filter={'access': 'open', 'file_id': uuids},
            fields=[
                'cases.project.project_id',
                'data_type',
                'analysis.workflow_type',
            ],
            method='POST',
        )
        try:
            df['cases'] = df['cases'].map(
                lambda c: ', '.join({p['project']['project_id'] for p in c})
            )
        except:  # noqa: E722
            pass
        df_list.append(df)
    df = pd.concat(df_list, axis=0)
    try:
        df = df.drop('id', axis=1)
    except KeyError:
        pass
    try:
        df = df.drop_duplicates()
    except:  # noqa: E722
        pass
    df.to_csv(sys.stdout, sep='\t', index=False)


def get_mapping(
    results,
    object_path,
    mapping,
    input_field_name,
    output_field_name,
    input_field_items,
):
    """Gets mapping between ``input_field_name`` and ``output_field_name``.
    This function is intended to be a utility function for ``map_fields``.

    Args:
        results (list or list of dicts): The output from the API.
        object_path (str): Path leading to the same object
            (``input_field_name`` and ``output_field_name``).
        mapping (dict): Mapping between ``input_field_name`` and
            ``output_field_name``. Keys are ``input_field_name`` and values are
            mapped ``output_field_name``.
        input_field_name (str): The field which is to be mapped to
            ``output_field_name``.
        output_field_name (str): ``input_field_name`` is to be mapped with
            ``output_field_name``.
        input_field_items (list of str or str): A list of input ids which is
            to mapped.
    For detailed explanation of the args please see ``map_fields``
    function.

    """

    for hit in results:
        for index in range(len(object_path)):
            if isinstance(hit, list):
                get_mapping(
                    hit,
                    object_path[index:],
                    mapping,
                    input_field_name,
                    output_field_name,
                    input_field_items,
                )
                return
            elif isinstance(hit, dict):
                if object_path[index] in hit:
                    hit = hit[object_path[index]]
        if not isinstance(hit, list):
            hit = [hit]
        for data in hit:
            field_1_leaves = []
            field_2_leaves = []
            extract_leaves(data, input_field_name, field_1_leaves)
            extract_leaves(data, output_field_name, field_2_leaves)
            for leaf in field_1_leaves:
                if leaf in input_field_items:
                    mapping[leaf] = field_2_leaves
    return mapping


def map_fields(
    endpoint,
    input_field_items,
    input_field_name,
    output_field_name,
):
    """Maps two fields - ``input_field_name`` and ``output_field_name`` from
    the GDC API in a dictionary. The fields should be within same object
    defined by the argument ``object_path``. The ``input_field_name`` and
    ``output_field_name`` should not be nested..


    Examples:
        Vaild examples would be:
        1.  object_path: "samples.portions.analytes.aliquots"
            input_field_name: "aliquot_id"
            output_field_name: "submitter_id"

        2.  object_path: "samples.portions"
            input_field_name: "analytes.aliquots.submitter_id"
            output_field_name: "submitter_id"

    Args:
        endpoint (str): One string of GDC API supported endpoint.
        object_path (str): Path leading to the same object (input_field_name
            and output_field_name). It should point to an object and should
            not be a property itself. For example,
            "samples.portions.analytes.aliquots.aliquot_id" is not a valid
            ``object_path`` since it leads to the property "aliquot_id". But,
            "samples.portions.analytes.aliquots" is a valid object_path since
            it consists of one or more property ("aliquot_id").
        input_field_items (list of str or str): A list of items that are of
            input field type which are to mapped. The items should be hashable
            since they will be used as key in the returned map.
        input_field_name (str): A direct property of the object pointed by
            ``object_path`` the field which is to be mapped to
            ``output_field_name``. List of all available fields can be found
            here
            https://docs.gdc.cancer.gov/API/Users_Guide/\
Appendix_A_Available_Fields/
        output_field_name (str): ``input_field_name`` is to be mapped with
            ``output_field_name``. If the field is not in the object pointed
            by the ``object_path``, the mapped values will be None. The mapped
            values can even be of ``list`` or ``dict`` datatype for proper
            ``output_field_name``.

    Returns:
        dict: The mapping between ``input_field_items`` and
            ``output_field_name``. Keys are ``input_field_items`` and the
            corresponding values are mapped ``output_field_items``.
    """

    if isinstance(input_field_items, str):
        input_field_items = [input_field_items]
    in_filter = {input_field_name: input_field_items}
    result = search(
        endpoint=endpoint,
        in_filter=in_filter,
        fields=[input_field_name, output_field_name],
        typ="json",
    )
    mapping = {}
    field_1 = input_field_name.split(".")
    field_2 = output_field_name.split(".")
    object_path = get_intersection(field_1, field_2)
    object_path_length = len(object_path)
    input_field_name = ".".join(
        input_field_name.split(".")[object_path_length:]
    )
    output_field_name = ".".join(
        output_field_name.split(".")[object_path_length:]
    )
    get_mapping(
        result,
        object_path,
        mapping,
        input_field_name,
        output_field_name,
        input_field_items,
    )
    for item in input_field_items:
        if item not in mapping:
            mapping[item] = None
    return mapping
