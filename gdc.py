#!/usr/bin/env python

# Ensure Python 2 and 3 compatibility
from __future__ import division
from __future__ import print_function

import gzip
import json
import os.path
import sys
import timeit

import numpy as np
import pandas as pd
import requests

GDC_API_BASE = 'https://api.gdc.cancer.gov'

def get_all_project_ids():
    """Get project ids for all projects in GDC
    
    Return: A list of project id for all projects in GDC
    """
    
    projects_endpt = '{}/projects'.format(GDC_API_BASE)
    projects_r = requests.get(projects_endpt)
    total_projects = projects_r.json()['data']['pagination']['total']
    url = '{}?size={}&fields=project_id'.format(projects_endpt, 
                                                total_projects)
    projects_r = requests.get(url)
    projects_list = projects_r.json()['data']['hits']
    return [p['project_id'] for p in projects_list]

def get_files_uuids(query_filter):
    """Get a list of file UUIDs matching query conditions
    
    Args:
        query_filter: A dict for query conditions conforming to GDC API's 
            query format. See:
            https://docs.gdc.cancer.gov/API/Users_Guide/Search_and_Retrieval/#filters-specifying-the-query
        
    Return: A list of file UUIDs matching query conditions
    """
    
    files_endpt = '{}/files'.format(GDC_API_BASE)
    params = {'filters':json.dumps(query_filter), 'fields':'file_id'}
    files_r = requests.post(files_endpt, data=params)
    size = files_r.json()['data']['pagination']['total']
    params['size'] = size
    files_r = requests.get(files_endpt, params=params)
    files_list = files_r.json()['data']['hits']
    return [f['file_id'] for f in files_list]

def and_eq_filter_constructor(filter_dict):
    """ A simple constructor for GDC API's query filter
    
    Convert a diction of query condition into a diction conforming to GDC 
    API's format. Conditions in the input dict will be combined together with 
    an AND relation. Every (key, value) pair will be joined together with "=" 
    operator.
    
    Args:
        filter_dict: A dict describing query conditions. Each (key, value) 
            pair represents for one condition. Operator between conditions is 
            "AND". Within a condition, the "key" is the 'field' operand. 
            Operator between "key" and "value" is "=".
    Return: A diction of filter conforming to GDC API's format. It should then 
        be converted to JSON format and used in the following http request.
    """
    
    operands_list = []
    for key in filter_dict:
        operands_list.append({"op":"=", "content":{"field":key,
                                                   "value":filter_dict[key]}})
    return {"op":"and", "content":operands_list}

def download_data(ids_dict, dir_path="", rename_ext_num=1):
    """Download GDC open access files
    
    Use GDC's data endpoint to download open access files stored in the GDC by 
    specifying file UUID(s). UUID(s) should be provided as a list of string. 
    Files will be "GET" from GDC one by one with progress information.
    
    Args:
        ids_list: A dict of UUIDs for files to be downloaded. Values are 
            UUIDs. Keys will be used in return dict.
        dir_path: Optinal. A string for downloading directory. Default 
            download directory is the current working directory.
        rename_ext_num: Optional. An int specify how many file name extensions 
            should be kept when renaming is required for name collision. 
            Default is to keep one last extension.
    Returns:
        A dict of absolute paths for downloaded files. Keys are corresponding
        keys in the input ids_list.
    """
    
    # Decide to download files one by one because multiple file downloading 
    # through GDC's API doesn't give file size which make download progress 
    # report impossible.
    
    total_count = len(ids_dict)
    file_count = 0
    data_endpt = '{}/data/'.format(GDC_API_BASE)
    chunk_size = 1024
    status = '\r[{{}}/{:d}] Downloading "{{}}": {{:3.0%}}'.format(total_count)
    file_dict = {}
    print('Download data to "{}"'.format(os.path.abspath(dir_path)))
    
    for key in ids_dict:
        file_count = file_count + 1
        downloaded = 0
        response = requests.get(data_endpt + ids_dict[key], stream=True)
        if response.status_code == 200:
            # Assume file name provided in the response header.
            content_disp = response.headers['Content-Disposition']
            file_name = content_disp[content_disp.find('filename=') + 9:]
            file_path = os.path.join(dir_path, file_name)
            # Handle file name collision
            if os.path.isfile(file_path):
                file_name_split = file_name.rsplit('.', rename_ext_num)
                file_name_split[0] = ids_dict[key]
                file_name = ".".join(file_name_split)
                file_path = os.path.join(dir_path, file_name)
            # Assume file size provided in the response header.
            file_size = int(response.headers['Content-Length'])
            with open(file_path, 'wb') as file:
                for chunk in response.iter_content(chunk_size):
                    file.write(chunk)
                    downloaded = downloaded + chunk_size
                    print(status.format(file_count, 
                                        file_name, 
                                        downloaded/file_size), 
                          end='')
                    sys.stdout.flush()
                file_dict[key] = os.path.abspath(file_path)
        print('\r[{}/{}] "{}" downloaded.      '.format(str(file_count),
                                                        str(total_count), 
                                                        file_name))
    return file_dict

def xena_matrix_update_RNA(file_dict, init_matrix = pd.DataFrame()):
    """Data transformation and assembly for RNAseq data
    
    Args:
        file_dict: A dict of file path.
        init_matrix: A matrix to be updated with new data.
    Returns:
        A pandas data frame for updated data matrix.
    """
    # TODO: Move to a separated module.
    # TODO: Name/Doc the function with its actual functionality so that it 
    #       doesn't has to be restricted to a specific data type?
    
    total_files = len(file_dict)
    new_matrix = init_matrix
    i = 0
    for key in file_dict:
        i = i + 1
        print('\rProcessing {}/{} file ...'.format(i, total_files), end='')
        with gzip.open(file_dict[key], 'r') as f:
            raw_col = pd.read_table(f, 
                                    sep='\t', 
                                    header=None, 
                                    names=[key], 
                                    index_col=0)
            transformed_col = np.log2(raw_col + 1)
            new_matrix = new_matrix.merge(transformed_col, 
                                          how='outer', 
                                          left_index=True, 
                                          right_index=True)
    print('\rAll {} files have been processed. '.format(total_files), end='')
    print('Data Transformation done.')
    return new_matrix

def label_files(uuids, label_field):
    """Query GDC with files' UUIDs for a proper sample label.
    
    Args:
        uuids: A list of UUIDs for data files to be labeled.
        label_field: A single GDC available file field whose value will be 
            used as the sample label.
    Return: A dict where each value is one input UUID and the key is the
        corresponding label.
    """
    
    label_key = label_field.split('.', 1)[0]
    files_endpt = '{}/files'.format(GDC_API_BASE)
    file_ids_filter = {"op":"in",
                       "content":{
                               "field":"file_id",
                               "value":uuids
                               }
                       }
    params = {'filters':json.dumps(file_ids_filter), 
              'fields':'file_id,{}'.format(label_field),
              'size':len(uuids)}
    response = requests.post(files_endpt, data=params)
    if response.status_code == 200:
        uuids_dict = {}
        for hit in response.json()['data']['hits']:
            label = hit[label_key]
            while isinstance(label, list) or isinstance(label, dict):
                if isinstance(label, list):
                    label = label[0]
                elif isinstance(label, dict):
                    label = label.values()
            uuids_dict[label] = hit['file_id']
    return uuids_dict

def main():
    start_time = timeit.default_timer()
    print('Script starts')
    query_test = {'cases.project.project_id': 'TCGA-BRCA',
                  'data_category': 'Transcriptome Profiling',
                  'files.analysis.workflow_type': 'HTSeq - FPKM-UQ'}

    file_ids_list = get_files_uuids(and_eq_filter_constructor(query_test))

    label_field = 'cases.samples.portions.analytes.aliquots.submitter_id'
    file_ids_dict = label_files(file_ids_list, label_field)

    print('Start to download at {} sec.'.format(timeit.default_timer() - start_time))
    file_dict = download_data(file_ids_dict, rename_ext_num=3)
    print('Download finishes at {} sec.'.format(timeit.default_timer() - start_time))

    test = xena_matrix_update_RNA(file_dict)
    print('Data transformation finishes at {} sec.'.format(timeit.default_timer() - start_time))
    test.to_csv('test.tsv', sep='\t')
    print('Script finishes at {} sec.'.format(timeit.default_timer() - start_time))

if __name__ == '__main__':
    main()
