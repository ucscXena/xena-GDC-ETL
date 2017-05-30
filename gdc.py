#!/usr/bin/env python

# Ensure Python 2 and 3 compatibility
from __future__ import division
from __future__ import print_function

import json
import os.path
import sys

import pandas as pd
import requests

_GDC_API_BASE = 'https://api.gdc.cancer.gov'
_SUPPORTED_FILE_TYPES = {'xml', 'txt', 'tar', 'gz', 'md5', 'xlsx', 'xls'}

def get_all_project_ids():
    """Get project ids for all projects on GDC.
    
    Return: A list of project id for all projects in GDC
    """
    
    projects_endpt = '{}/projects'.format(_GDC_API_BASE)
    projects_r = requests.get(projects_endpt)
    total_projects = projects_r.json()['data']['pagination']['total']
    url = '{}?size={}&fields=project_id'.format(projects_endpt, 
                                                total_projects)
    projects_r = requests.get(url)
    projects_list = projects_r.json()['data']['hits']
    return [p['project_id'] for p in projects_list]

def get_files_uuids(query_filter):
    """Get UUIDs for files matching query conditions.
    
    Args:
        query_filter: A dict for query conditions conforming to GDC API's 
            query format. See:
            https://docs.gdc.cancer.gov/API/Users_Guide/Search_and_Retrieval/#filters-specifying-the-query
        
    Return: A list of file UUIDs matching query conditions
    """
    
    files_endpt = '{}/files'.format(_GDC_API_BASE)
    params = {'filters':json.dumps(query_filter), 'fields':'file_id'}
    files_r = requests.post(files_endpt, data=params)
    size = files_r.json()['data']['pagination']['total']
    params['size'] = size
    files_r = requests.get(files_endpt, params=params)
    files_list = files_r.json()['data']['hits']
    return [f['file_id'] for f in files_list]

def and_eq_filter_constructor(filter_dict):
    """ A simple constructor converting a query dictionary into GDC API 
    specific filters.
    
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

def download_data(ids_dict, dir_path=""):
    """Download open access data from GDC according to input UUIDs.
    
    Use GDC's data endpoint to download open access files stored in the GDC by 
    specifying file UUID(s). UUID(s) should be provided as a list of string. 
    Files will be "GET" from GDC one by one with progress information.
    
    Args:
        ids_list: A dict of UUIDs for files to be downloaded. Values are 
            UUIDs. Keys will be used in return dict.
        dir_path: Optinal. A string for downloading directory. Default 
            download directory is the current working directory.
    Returns:
        A dict of absolute paths for downloaded files. Keys are corresponding
        keys in the input ids_list.
    """
    
    # Decide to download files one by one because multiple file downloading 
    # through GDC's API doesn't give file size which make download progress 
    # report impossible.
    
    total_count = len(ids_dict)
    file_count = 0
    data_endpt = '{}/data/'.format(_GDC_API_BASE)
    chunk_size = 1024
    status = '\r[{{}}/{:d}] Download "{{}}": {{:3.0%}}'.format(total_count)
    file_dict = {}
    print('Download data to "{}"'.format(os.path.abspath(dir_path)))
    
    for key in ids_dict:
        file_count = file_count + 1
        downloaded = 0
        response = requests.get(data_endpt + ids_dict[key], stream=True)
        if response.status_code == 200:
            # Assume file name provided in the response header.
            content_disp = response.headers['Content-Disposition']
            gdc_file_name = content_disp[content_disp.find('filename=') + 9:]
            # https://github.com/broadinstitute/gdctools/blob/master/gdctools/lib/meta.py
            name_list = gdc_file_name.split('.')
            for i in range(len(name_list)):
                if name_list[i] in _SUPPORTED_FILE_TYPES:
                    break
            download_file_name = '.'.join([key, ids_dict[key]] + name_list[i:])
            file_path = os.path.join(dir_path, download_file_name)
            # Assume file size provided in the response header.
            file_size = int(response.headers['Content-Length'])
            with open(file_path, 'wb') as file:
                for chunk in response.iter_content(chunk_size):
                    file.write(chunk)
                    downloaded = downloaded + chunk_size
                    print(status.format(file_count, 
                                        download_file_name, 
                                        downloaded/file_size), 
                          end='')
                    sys.stdout.flush()
                file_dict[key] = os.path.abspath(file_path)
        print('')
    return file_dict

def label_files(uuids, label_field):
    """Label a list of file UUIDs with their corresponding field of selection. 
    
    Query GDC with files' UUIDs for proper sample labels, such as their 
    corresponding "aliquots.submitter_id".
    
    Args:
        uuids: A list of UUIDs for data files to be labeled.
        label_field: A single GDC available file field whose value will be 
            used as the sample label.
    Return: A dict where each value is one input UUID with its key being the
        corresponding label.
    """
    
    label_key = label_field.split('.', 1)[0]
    files_endpt = '{}/files'.format(_GDC_API_BASE)
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
    uuids_dict = {}
    if response.status_code == 200:
        for hit in response.json()['data']['hits']:
            label = hit[label_key]
            while isinstance(label, list) or isinstance(label, dict):
                if isinstance(label, list):
                    label = label[0]
                elif isinstance(label, dict):
                    label = label.values()
            uuids_dict[label] = hit['file_id']
    return uuids_dict

def get_all_case_info():
    """Get some basic information for all cases on GDC.
    """
    
    cases_endpt = '{}/cases'.format(_GDC_API_BASE)
    fields_list = ['case_id',
                   'project.project_id',
                   'project.primary_site',
                   'project.disease_type',
                   'submitter_id']
    expand_list = ['demographic',
                   'diagnoses']
    size = 10
    params = {'fields':','.join(fields_list), 
              'expand':','.join(expand_list),
              'size':size}
    cur_page = 0
    total_pages = 1
    cases_dataframe = pd.DataFrame()
    # Cases Endpoint on GDC fail to handle big size query; have to iterate 
    # over pages.
    while (total_pages > cur_page):
        from_record = size * cur_page + 1
        params['from'] = from_record
        cases_r = requests.post(cases_endpt, data=params)
        cur_page = cases_r.json()['data']['pagination']['page']
        total_pages = cases_r.json()['data']['pagination']['pages']
        print('\rProcessing page {}/{}...'.format(cur_page, total_pages), 
              end='')
        for hit in cases_r.json()['data']['hits']:
            case_record = {}
            if 'diagnoses' in hit:
                case_record.update({('diagnoses', key): value 
                                    for key, value 
                                    in hit['diagnoses'][0].items()})
            if 'demographic' in hit:
                case_record.update({('demographic', key): value 
                                    for key, value 
                                    in hit['demographic'].items()})
            if 'project' in hit:
                case_record.update({('project', key): value 
                                    for key, value 
                                    in hit['project'].items()})
            case_record[('cases', 'submitter_id')] = hit['submitter_id']
            case_df = pd.DataFrame(case_record, 
                                   index=[hit['case_id']]).sort_index(axis=1)
            cases_dataframe = cases_dataframe.append(case_df)
    cases_dataframe.to_csv('cases.tsv', sep='\t')

def main():
    print('A simple python module for GDC API functionality.')
    #get_all_case_info()

if __name__ == '__main__':
    main()
