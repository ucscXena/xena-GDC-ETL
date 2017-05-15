#!/usr/bin/env python

import requests
import json
import gzip
import pandas as pd
import os.path
import sys

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
    files_r = requests.get(files_endpt, params=params)
    size = files_r.json()['data']['pagination']['total']
    params['size'] = size
    files_r = requests.get(files_endpt, params=params)
    files_list = files_r.json()['data']['hits']
    return [f['file_id'] for f in files_list]

def and_eq_filter_constructor(filter_dict):
    operands_list = []
    for key in filter_dict:
        operands_list.append({"op":"=", "content":{"field":key,
                                                   "value":filter_dict[key]}})
    return {"op":"and", "content":operands_list}

def download_data(ids_list, dir_path="", rename_ext_num=1):
    """Download GDC open access files
    
       Use GDC's data endpoint to download open access files stored in the GDC 
       by specifying file UUID(s). UUID(s) should be provided as a list of 
       string. Files will be "GET" from GDC one by one with progress 
       information.
       
       Args:
           ids_list: A list of string for file's UUID.
           dir_path: Optinal. A string for downloading directory. Default 
               download directory is the current working directory.
           rename_ext_num: Optional. An int specify how many file name 
               extensions should be kept when renaming is required for name 
               collision. Default is to keep one last extension.
    """
    
    # Decide to download files one by one because multiple file downloading 
    # through GDC's API doesn't give file size which make download progress 
    # report impossible.
    
    total_count = len(ids_list)
    file_count = 0
    data_endpt = '{}/data/'.format(GDC_API_BASE)
    chunk_size = 1024
    status = '\r[{{}}/{:d}] Downloading "{{}}": {{:3.0%}}'.format(total_count)
    print('Download data to "{}"'.format(os.path.abspath(dir_path)))
    
    for uuid in ids_list:
        file_count = file_count + 1
        downloaded = 0
        response = requests.get(data_endpt + uuid, stream=True)
        if response.status_code == 200:
            # Assume file name provided in the response header.
            content_disp = response.headers['Content-Disposition']
            file_name = content_disp[content_disp.find('filename=') + 9:]
            file_path = os.path.join(dir_path, file_name)
            # Handle file name collision
            if os.path.isfile(file_path):
                file_name_split = file_name.rsplit('.', rename_ext_num)
                file_name_split[0] = uuid
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
        print('\r[{}/{}] "{}" downloaded.      '.format(str(file_count),
                                                        str(total_count), 
                                                        file_name))

def main():
    query_test = {'cases.project.project_id': 'TCGA-BRCA',
                  'data_category': 'Transcriptome Profiling',
                  'files.analysis.workflow_type': 'HTSeq - FPKM-UQ'}
    file_ids_list = get_files_uuids(and_eq_filter_constructor(query_test))
    print(len(file_ids_list))
    download_data(file_ids_list[0:10], rename_ext_num=3)
    with gzip.open('0d604955-ce10-4fa9-ba8b-9417c5719653.FPKM-UQ.txt.gz', 'r') as f:
        test = pd.read_table(f, names=["ensemble_id", "FPKM_UQ"])
    print(test.head(5))

if __name__ == '__main__':
    main()
