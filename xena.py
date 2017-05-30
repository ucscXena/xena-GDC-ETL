#!/usr/bin/env python

# Ensure Python 2 and 3 compatibility
from __future__ import division
from __future__ import print_function

import gzip
import os
import timeit

import numpy as np
import pandas as pd

import gdc

def xena_matrix_update_RNA(file_dict, init_matrix=pd.DataFrame(), 
                           update=True, deletes=None, index_cols=None, 
                           data_col=None):
    """Transformation (log2(x+1)) and data matrix assembly for RNA-seq data 
    (HTSeq - Counts, HTSeq - FPKM and HTSeq - FPKM-UQ).
    
    Args:
        file_dict: A dict of file path, with the key being sample label and 
            the value being one absolute path for the data.
        init_matrix: A Xena matrix loaded as a Pandas DataFrame which be 
            updated with new data. Default is an empty DataFrame, i.e. the 
            function will return a new matrix in Pandas DataFrame format.
        update: A boolean switch decides whether a new column should overwrite 
            a column in the matrix if it has an overlapping column name. 
            Default is to overwrite the old column. If given False, the 
            funciton will rename both old and new columns.
        deletes: A list of strings which correspond to the name of columns to 
            be deleted.
    Returns:
        A Pandas DataFrame for Xena matrix.
    """
    # TODO: Name/Doc the function with its actual functionality so that it 
    #       doesn't has to be restricted to a specific data type?
    
    total_files = len(file_dict)
    new_matrix = init_matrix
    if isinstance(deletes, list):
        for column in deletes:
            if column in new_matrix.columns:
                new_matrix.drop(column, axis=1, inplace=True)
            else:
                message = ('Sample {} was not found in matrix. '
                           'Nothing to delete.')
                print(message.format(column))
    i = 0
    for key in file_dict:
        i = i + 1
        print('\rProcessing {}/{} file...'.format(i, total_files), end='')
        with gzip.open(file_dict[key], 'r') as f:
            raw_col = pd.read_table(f, sep='\t', header=None, names=[key], 
                                    index_col=0)
            transformed_col = np.log2(raw_col + 1)
            if update and (key in new_matrix.columns):
                new_matrix.drop(key, axis=1, inplace=True)
            new_matrix = new_matrix.merge(transformed_col, how='outer', 
                                          left_index=True, right_index=True)
    print('\rAll {} files have been processed. '.format(total_files), end='')
    print('Data Transformation done.')
    return new_matrix

def main():
    """The main pipeline for importing GDC data into Xena.
    """

    xena_root = os.path.abspath('.')
    print('Current Xena root directory is {}'.format(xena_root))
    label_field = 'cases.samples.portions.analytes.aliquots.submitter_id'
    data_type_tuples = [('Gene Expression Quantification', 'HTSeq - Counts'),
                        ('Gene Expression Quantification', 'HTSeq - FPKM'),
                        ('Gene Expression Quantification', 'HTSeq - FPKM-UQ')]

    project_ids_list = gdc.get_all_project_ids()
    # Remove the line below to get RNA-Seq data for all 39 projects on GDC.
    project_ids_list = ['TCGA-BRCA']

    start_time = timeit.default_timer()
    for project_id in project_ids_list:
        status = 'Start to work on project "{}".'
        print(status.format(project_id))
        project_dir = os.path.join(xena_root, project_id)
        try: 
            os.makedirs(project_dir)
        except OSError:
            if not os.path.isdir(project_dir):
                raise
        for (data_type, workflow_type) in data_type_tuples:
            data_dir_name = '{}_{}'.format(data_type.replace(' ', ''),
                                           workflow_type.replace(' ', ''))
            data_dir = os.path.join(project_dir, data_dir_name)
            try: 
                os.makedirs(data_dir)
            except OSError:
                if not os.path.isdir(data_dir):
                    raise
            query_dict = {'cases.project.project_id': project_id,
                          'data_type': data_type,
                          'analysis.workflow_type': workflow_type}
            gdc_filter = gdc.and_eq_filter_constructor(query_dict)
            file_ids_list = gdc.get_files_uuids(gdc_filter)
            if not file_ids_list:
                message = 'No {} - {} data for "{}" project.'
                print(message.format(data_type, workflow_type, project_id))
                continue
            file_ids_dict = gdc.label_files(file_ids_list, label_field)
            if not file_ids_dict:
                message = 'Fail to label {} - {} data for "{}" project.'
                print(message.format(data_type, workflow_type, project_id))
                continue
            
            status = 'Start to download {} - {} for {} at {} sec.'
            print(status.format(data_type, workflow_type, project_id, 
                                timeit.default_timer() - start_time))
            files_dict = gdc.download_data(file_ids_dict, dir_path=data_dir)
            status = 'Download finishes at {} sec.'
            if not files_dict:
                print('Download failed.')
                message = 'No downloaded {} - {} data for "{}" project.'
                print(message.format(data_type, workflow_type, project_id))
                continue
            print(status.format(timeit.default_timer() - start_time))
            
            status = 'Start to transform {} - {} data for {} at {} sec.'
            print(status.format(data_type, workflow_type, project_id,
                                timeit.default_timer() - start_time))
            new_matrix = xena_matrix_update_RNA(files_dict, update=True)
            status = 'Data transformation finishes at {} sec.'
            print(status.format(timeit.default_timer() - start_time))
            
            status = 'Start to save new matrix at {} sec.'
            print(status.format(timeit.default_timer() - start_time))
            new_matrix_name = '{}_{}_{}.tsv'.format(project_id, data_type, 
                                                    workflow_type)
            new_matrix_path = os.path.join(project_dir, new_matrix_name)
            new_matrix.to_csv(new_matrix_path, sep='\t')
            status = 'New matrix saved at {} sec.'
            print(status.format(timeit.default_timer() - start_time))

if __name__ == '__main__':
    main()
