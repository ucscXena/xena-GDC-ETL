#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This module provides basic and minimum necessary functions for carrying out 
data query and download for Xena GDC ETL pipelines
"""

# Ensure Python 2 and 3 compatibility
from __future__ import division
from __future__ import print_function

import json
import sys

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
    
    operands_list = []
    for key in filter_dict:
        value = filter_dict[key]
        if not isinstance(value, list):
            value = [value]
        operands_list.append({"op":"in", 
                              "content":{"field":key, "value":value}})
    return {"op":"and", "content":operands_list}

def get_file_dict(query_filter, label_field=None):
    """Get a dictionary of files matching query conditions.
    
    Args:
        query_filter: dict
            A dict of query conditions which will be used to search for files 
            of interest. It should follow GDC API's query format. See:
            https://docs.gdc.cancer.gov/API/Users_Guide/Search_and_Retrieval/#filters-specifying-the-query
        label_field: A single GDC available file field whose value will be 
            used for renaming downloaded files. Default is None which makes 
            the final file name become "UUID.file_extension". If provided, 
            the final file name will be "label.UUID.file_extension". GDC 
            available file fields can be found at 
            https://docs.gdc.cancer.gov/API/Users_Guide/Appendix_A_Available_Fields/
            File_extensions supported by this module are defined by the 
            constant set _SUPPORTED_FILE_TYPES.
        
    Returns: 
        file_dict: A dict of files matching query conditions. The key will be 
        a file's UUID, while the value is the file name formatted accordingly, 
        which can be used during downloading.
    """
    
    files_endpt = '{}/files'.format(_GDC_API_BASE)
    params = {'filters':json.dumps(query_filter), 'size':1}
    if label_field is None:
        params['fields'] = 'file_id,file_name'
    else:
        params['fields'] = 'file_id,file_name,{}'.format(label_field)
    response = requests.post(files_endpt, data=params)
    params['size'] = response.json()['data']['pagination']['total']
    response = requests.get(files_endpt, params=params)

    file_dict = {}
    if response.status_code == 200:
        if label_field is None:
            root_field_group = 'file_id'
        else:
            root_field_group = label_field.split('.', 1)[0]
        for hit in response.json()['data']['hits']:
            label = hit[root_field_group]
            while isinstance(label, list) or isinstance(label, dict):
                if isinstance(label, list):
                    label = label[0]
                elif isinstance(label, dict):
                    label = label.values()
            file_id = hit['file_id']
            # https://github.com/broadinstitute/gdctools/blob/master/gdctools/lib/meta.py
            name_list = hit['file_name'].split('.')
            for i in range(len(name_list)):
                if name_list[i] in _SUPPORTED_FILE_TYPES:
                    break
            file_name = '.'.join([label, file_id] + name_list[i:])
            file_dict[file_id] = file_name
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

def get_all_project_info():
    """Get project info for all projects on GDC.
    
    Returns:
        project_dict: A dict of project info for all projects in GDC. The key 
        will be project_id and the corresponding value is a dict contains 
        "project name", "primary_site" and "program name" info.
    """
    
    projects_endpt = '{}/projects'.format(_GDC_API_BASE)
    projects_r = requests.get(projects_endpt)
    total_projects = projects_r.json()['data']['pagination']['total']
    fields = 'name,primary_site,project_id,program.name'
    url = '{}?size={}&fields={}'.format(projects_endpt, total_projects, fields)
    projects_r = requests.get(url)
    project_dict = {}
    for project in projects_r.json()['data']['hits']:
        project['program'] = project['program']['name']
        project['primary_site'] = project['primary_site'][0]
        project_dict[project.pop('project_id')] = project
    return project_dict

def main():
    print('A simple python module providing selected GDC API functionality.')

if __name__ == '__main__':
    main()
