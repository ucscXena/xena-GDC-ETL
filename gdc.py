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

_GDC_API_BASE = 'https://api.gdc.cancer.gov'
_SUPPORTED_FILE_TYPES = {'txt', 'vcf', 'bam', 'tsv', 'xml', 'maf', 'xlsx', 
                         'tar', 'gz', 'md5', 'xls'}
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
        {'analysis.workflow_type': 
            'SomaticSniper Variant Aggregation and Masking'},
        {'analysis.workflow_type': 'VarScan2 Variant Aggregation and Masking'},
        {'data_type': 'Biospecimen Supplement'},
        {'data_type': 'Clinical Supplement'}
    ]

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
        operation_list.append({"op":"in", 
                               "content":{"field":key, "value":value}})
    for key in exclude_dict:
        value = exclude_dict[key]
        if not isinstance(value, list):
            value = [value]
        operation_list.append({"op":"exclude", 
                               "content":{"field":key, "value":value}})
    return {"op":"and", "content":operation_list}

def traverse_field_json(data, field=None):
    """Assuming the nested JSON having a single (set of) data, walk into/down 
    this nested JSON and return the value.
    
    A lot of times, a single GDC data is kept in a highly nested JSON objects. 
    For example, in a search of the files endpoint, a sample's submitter ID 
    "TCGA-C8-A133-01A" is kept as::
    
        "cases": [
          {
            "samples": [
              {
                "submitter_id": "TCGA-C8-A133-01A"
              }
            ]
          }
        ]
    
    If multiple sets of data are found in the nested JSON, a "ValueError" will 
    be raised.
    
    Args:
        data (str): A string of nested JSON. This JSON should contain only one 
            set of data, meaning arrays nested in this JSON should be size 1. 
            Otherwise, a "ValueError" will be raised.
        field (str, optional): A string representing the structure of the 
            nested JSON. Keys for each level are separated by ".". If None, 
            the nested JSON must contain only one data (rather than one set of 
            data),  i.e. only one key at each level of nesting. Otherwise, a 
            "ValueError" will be raised. Defaults to None.
        
    Returns: 
        str or else: One data specified by the "field" or the only data 
        contained in the nested JSON. It should be str most of the time but 
        can be any types except list and dict.

    """
    
    if field is None:
        while isinstance(data, list) or isinstance(data, dict):
            if len(data) > 1:
                raise ValueError('Multiple data found in a list!')
            elif isinstance(data, list):
                data = data[0]
            else:
                data = data.values()[0]
    else:
        keys = field.split('.')
        for k in keys:
            if isinstance(data, list):
                if len(data) > 1:
                    raise ValueError('Multiple data found in a list!')
                else:
                    data = data[0]
            data = data[k]
    return data

def search(endpoint, fields, in_filter, exclude_filter={}, expand=[]):
    """Search one GDC endpoints and return searching results in a pandas 
    DataFrame if possible.
    
    When searching results cannot be safely converted to a pandas DataFrame, 
    results will be returned in the JSON format as it is returned from GDC 
    API.
    
    Args:
        endpoint (str): One string of GDC API supported endpoint. See:
            https://docs.gdc.cancer.gov/API/Users_Guide/Getting_Started/#api-endpoints
        fields (list or str): One or more fields to be queried. Each field 
            will be used as a column name in the returned DataFrame. It can be 
            a comma separated string or a list of field strings or a 
            combination of both.
        in_filter (dict): A dict of query conditions which will be used to 
            perform the query. Each (key, value) pair represents for one 
            condition. It will be passed to ``simple_and_filter`` for making a 
            query filter compatible with GDC API. Please check 
            ``simple_and_filter`` function for details.
        exclude_filter (dict): An optional dict of query conditions which will 
            be used to perform the query. Each (key, value) pair represents 
            for one condition. It will be passed to ``simple_and_filter`` for 
            making a query filter compatible with GDC API. Please check 
            ``simple_and_filter`` function for details.
        expand (list or str): One or more field groups to be queried. It can 
            be a comma separated string or a list of field strings or a 
            combination of both.
    
    Returns:
        pandas.core.frame.DataFrame or str: A search result in form of a 
            pandas DataFrame or a JSON formatted string.
    """
    
    url = '{}/{}'.format(_GDC_API_BASE, endpoint)
    filters =  simple_and_filter(in_dict=in_filter, 
                                 exclude_dict=exclude_filter)
    if isinstance(fields, str):
        fields = [fields]
    if isinstance(expand, str):
        expand = [expand]
    params = {'filters':json.dumps(filters), 
              'size':1, 
              'fields':','.join(fields),
              'expand':','.join(expand)}
    response = requests.post(url, data=params)
    params['size'] = response.json()['data']['pagination']['total']
    response = requests.get(url, params=params)
    if response.status_code == 200:
        results = response.json()['data']['hits']
        try:
            df = pd.read_json(json.dumps(results), orient='records', 
                              typ='frame')
            col_to_del = []
            for field in [f for f in fields if '.' in f]:
                col, keys = field.split('.', 1)
                df[field] = df[col].apply(
                        lambda x: traverse_field_json(x, keys)
                    )
                col_to_del.append(col)
            for col in set(col_to_del):
                df = df.drop(col, axis=1)
            results = df
        except Exception:
            warnings.warn('Fail to convert searching results into table. '
                          'JSON will be returned.', stacklevel=2)
    else:
        results = None
    return results

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
        raise TypeError('uuids is a {}; it should be a string, a list or a '
                        'dict'.format(type(uuids)))
    total = len(uuids)
    count = 0
    download_list = []
    data_endpt = '{}/data/'.format(_GDC_API_BASE)
    for uuid in uuids:
        count += 1
        response = requests.get(data_endpt + uuid, stream=True)        
        if response.status_code == 200:
            file_size = int(response.headers['Content-Length'])
            if uuids[uuid] is None:
                content_disp = response.headers['Content-Disposition']
                ori_name = content_disp[content_disp.find('filename=') + 9:]
                new_filename = uuid + '.' + get_ext(ori_name)
                path = os.path.join(os.path.abspath(download_dir), 
                                    new_filename)
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
                    print(status.format(count, total, path, 
                                        min(1, downloaded/file_size)), end='')
                    sys.stdout.flush()
            download_list.append(path)
        else:
            print('\rFail to download file {}.'.format(uuid))
    print('')
    return download_list

def get_all_project_info():
    """Get project info for all projects on GDC.
    
    Returns:
        dict: A dict of project info for all projects in GDC. The key 
        will be project_id and the corresponding value is a dict contains 
        "project name", "primary_site" and "program name" info.
    """
    
    fields = ['name', 'primary_site', 'project_id', 'program.name']
    project_df = search('projects', fields, {})
    return project_df.set_index('id')

def main():
    print('A simple python module providing selected GDC API functionalities.')
    
    # Simple test
    print(get_all_project_info().head())
    
if __name__ == '__main__':
    main()
