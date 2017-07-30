#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This module provides basic and minimum necessary functions for carrying out 
data query and download for Xena GDC ETL pipelines
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
_SUPPORTED_FILE_TYPES = {'xml', 'txt', 'tar', 'gz', 'md5', 'xlsx', 'xls'}
_SUPPORTED_DATASETS = [{'data_type': 'Copy Number Segment'},
                       {'data_type': 'Masked Copy Number Segment'},
                       {'data_type': 'Isoform Expression Quantification'},
                       {'data_type': 'miRNA Expression Quantification'},
                       {'data_type': 'Methylation Beta Value'},
                       {'analysis.workflow_type': 'HTSeq - Counts'},
                       {'analysis.workflow_type': 'HTSeq - FPKM'},
                       {'analysis.workflow_type': 'HTSeq - FPKM-UQ'},
                       {'analysis.workflow_type': 'MuSE Variant Aggregation and Masking'},
                       {'analysis.workflow_type': 'MuTect2 Variant Aggregation and Masking'},
                       {'analysis.workflow_type': 'SomaticSniper Variant Aggregation and Masking'},
                       {'analysis.workflow_type': 'VarScan2 Variant Aggregation and Masking'},
                       {'data_type': 'Biospecimen Supplement'},
                       {'data_type': 'Clinical Supplement'}]

def and_in_filter_constructor(filter_dict):
    """A simple constructor converting a query dictionary into GDC API 
    specific filters.
    
    Convert a dict of query condition into a diction conforming to GDC 
    API's format. Conditions in the input dict will be combined together by 
    an AND relation. Every (key, value) pair will be joined together with "in" 
    operator. If value is not a list, it will be converted to a list first. 
    See details at 
    https://docs.gdc.cancer.gov/API/Users_Guide/Search_and_Retrieval/#filters-specifying-the-query
    
    Args:
        filter_dict: dict
            A dict describing query conditions. Each (key, value) pair 
            represents for one condition. Operator between conditions is 
            "AND". Within a condition, the "key" is the 'field' operand. 
            Operator between "key" and "value" is "in".
    Returns: 
        query_dict: A dict of filter conforming to GDC API's format. It should 
        then be converted to JSON format and used in the following http 
        request.
    """
    
    if not filter_dict:
        return filter_dict
    operands_list = []
    for key in filter_dict:
        value = filter_dict[key]
        if not isinstance(value, list):
            value = [value]
        operands_list.append({"op":"in", 
                              "content":{"field":key, "value":value}})
    return {"op":"and", "content":operands_list}

def traverse_field_json(data, field=None):
    """Assuming the nested JSON having a single (set of) data, walk into/down 
    this nested JSON and return the value.
    
    A lot of times, a single GDC data is kept in a highly nested JSON objects. 
    For example, in a search of the files endpoint, a sample's submitter ID 
    "TCGA-C8-A133-01A" is kept as:
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
        data: str
            A string of nested JSON. This JSON should contain only one set of 
            data, meaning arrays nested in this JSON should be size 1. 
            Otherwise, a "ValueError" will be raised.
        field: str, default None
            A string representing the structure of the nested JSON. Keys for 
            each level are separated by ".". If None, the nested JSON must 
            contain only one data (rather than one set of data),  i.e. only 
            one key at each level of nesting. Otherwise, a "ValueError" will 
            be raised.
        
    Returns: 
        data: One data specified by the "field" or the only data contained in 
        the nested JSON. It should be str most of the time but can any types 
        except list and dict.

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

def search(endpoint, filter_dict, fields):
    """Search one GDC endpoints and return searching results in a pandas 
    DataFrame if possible.
    
    When searching results cannot be safely converted to a pandas DataFrame, 
    results will be returned in the JSON format as it is returned from GDC 
    API.
    
    Args:
        endpoint: str
            One string of GDC API supported endpoint. See:
            https://docs.gdc.cancer.gov/API/Users_Guide/Getting_Started/#api-endpoints
        filter_dict: str
            A dict of query conditions which will be used to perform the 
            query. Each (key, value) pair represents for one condition. 
            It will be passed to "and_in_filter_constructor" for making a 
            query filter compatible with GDC API. Please check 
            "and_in_filter_constructor" function for details.
        fields: list or str
            One or more fields to be queried. Each field will be used as a 
            column name in the returned DataFrame. It can be a comma separated 
            string or a list of field strings or a combination of both.
    
    Returns:
        results: A pandas DataFrame or a JSON formatted string.
    """
    
    url = '{}/{}'.format(_GDC_API_BASE, endpoint)
    filters =  and_in_filter_constructor(filter_dict)
    if isinstance(fields, str):
        fields = [fields]
    params = {'filters':json.dumps(filters), 
              'size':1, 
              'fields':','.join(fields)}
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
        file_name: str
            The filename will be split by "." and checked from left to right. 
            Extensions will be kept starting from the first (and left most) 
            supported extension.
    
    Returns:
        exts: A string of extensions joint by ".".
    """
    
    # https://github.com/broadinstitute/gdctools/blob/master/gdctools/lib/meta.py
    name_list = file_name.split('.')
    for i in range(len(name_list)):
        if name_list[i] in _SUPPORTED_FILE_TYPES:
            break
    return '.'.join(name_list[i:])

def get_file_dict(filter_dict, label_field=None):
    """Get a dictionary of files matching query conditions.
    
    Args:
        filter_dict: dict
            A dict of query conditions which will be used to search for files 
            of interest. Each (key, value) pair represents for one condition. 
            It will be passed to "and_in_filter_constructor" for making a 
            query filter compatible with GDC API. Please check 
            "and_in_filter_constructor" function for details.
        label_field: str, default None
            A single GDC available file field whose value will be 
            used for renaming downloaded files. Default is None which makes 
            the final file name become "UUID.file_extension". If provided, the 
            final file name will be "label.UUID.file_extension". GDC available 
            file fields can be found at 
            https://docs.gdc.cancer.gov/API/Users_Guide/Appendix_A_Available_Fields/
            File_extensions supported by this module are defined by the 
            constant set _SUPPORTED_FILE_TYPES.
        
    Returns: 
        file_dict: A dict of files matching query conditions. The key will be 
        a file's UUID, while the value is the file name formatted accordingly, 
        which can be used during downloading.
    """
    
    fields = ['file_id', 'file_name']
    if label_field is not None:
        fields.append(label_field)
    try:
        file_df = search('files', filter_dict, fields)
        file_df.set_index('file_id', drop=False, inplace=True)
        if label_field is None:
            file_s = (file_df['file_id'].astype(str) + '.' 
                      + file_df['file_name'].apply(get_ext).astype(str))
        else:
            file_s = (file_df[label_field].astype(str) + '.' 
                      + file_df['file_id'].astype(str) + '.' 
                      + file_df['file_name'].apply(get_ext).astype(str))
        file_dict = file_s.to_dict()
    except:
        file_dict = {}
    return file_dict

def download_data(uuid, path=None):
    """ Download a single file from GDC.
    
    Args:
        uuid: str
            UUID for the target file.
        path: str, default None.
            Path for saving the downloaded file. If None, the downloaded file 
            will be saved under current python work directory with the 
            original filename from GDC.
    
    Returns:
        status: boolean
    """
    
    data_endpt = '{}/data/'.format(_GDC_API_BASE)
    chunk_size = 4096
    response = requests.get(data_endpt + uuid, stream=True)
    if response.status_code == 200:
        file_size = int(response.headers['Content-Length'])
        if path is None:
            content_disp = response.headers['Content-Disposition']
            path = content_disp[content_disp.find('filename=') + 9:]
        with open(path, 'wb') as f:
            downloaded = 0
            print('  0%', end='')
            for chunk in response.iter_content(chunk_size):
                f.write(chunk)
                downloaded = downloaded + chunk_size
                print('{}{:4.0%}'.format('\b' * 4, downloaded/file_size), 
                      end='')
                sys.stdout.flush()
            print('\b\b\b\b100%')
        return True
    else:
        return False

def download(uuids, download_dir='.', chunk_size=4096):
    """Download GDC's open access data according to UUID(s).
    
    Args:
        uuids: str, list, or dict
            A single UUID (str), a list of UUIDs (list) or a dict whose keys 
            are UUIDs for target file(s). If "uuids" is str or list, 
            downloaded file(s) will be renamed to "UUID.extension" where 
            "extension" is extracted by "get_ext()" from the original 
            filename. Renamed file(s) will be saved at "download_dir". If 
            "uuids" is a dict, the argument "download_dir" will be ignored; 
            values of dict will be paths for saving corresponding downloaded 
            files.
        download_dir: str, default '.'
            The directory for saving downloaded file(s) when "uuids" is str or 
            list. It will be ignored if "uuids" is dict.
        chunk_size: int, default 4096
            The chunk size is the number of bytes it should read into memory, 
            when the response is got with "stream=True". Check the 
            documentation of "requests" module for details.
    
    Returns:
        download_list: a list of paths for downloaded files.
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
            with open(path, 'wb') as f:
                downloaded = 0
                print(status.format(count, total, path, 0), end='')
                sys.stdout.flush()
                for chunk in response.iter_content(chunk_size):
                    f.write(chunk)
                    downloaded = downloaded + chunk_size
                    print(status.format(count, total, path, 
                                        downloaded/file_size), end='')
                    sys.stdout.flush()
            download_list.append(path)
        else:
            print('\rFail to download file {}.'.format(uuid))
    print('')
    return download_list

def get_all_project_info():
    """Get project info for all projects on GDC.
    
    Returns:
        project_dict: A dict of project info for all projects in GDC. The key 
        will be project_id and the corresponding value is a dict contains 
        "project name", "primary_site" and "program name" info.
    """
    
    fields = ['name', 'primary_site', 'project_id', 'program.name']
    project_df = search('projects', {}, fields)
    return project_df.set_index('id')

def main():
    print('A simple python module providing selected GDC API functionality.')
    
if __name__ == '__main__':
    main()
