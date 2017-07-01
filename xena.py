#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Ensure Python 2 and 3 compatibility
from __future__ import division
from __future__ import print_function

import datetime
import os
import sys

import jinja2
import numpy as np
import pandas as pd

import gdc

def mkdir_p(dir_name):
    """Make the directory as needed: no error if existing.
    
    Args:
        dir_name: str
            directory name or path
    
    Returns:
        path: absolute path for the directory.
    """
    
    dir_path = os.path.abspath(dir_name)
    try:
        os.makedirs(dir_path)
    except OSError:
        if not os.path.isdir(dir_path):
            raise
    return dir_path

def download_dataset(dataset_type, projects=None, work_dir='.'):
    """Download open access data based on "project_id(s)" and "dataset_type" 
    and put them together under one directory as a single dataset.
    
    Files will be renamed into 
    "<cases.samples.submitter_id[:-1]>.<UUID>.<file extension>"
    
    Args:
        dataset_type: dict
            A dict specifying dataset(s) of interest on GDC. It will be used 
            directly for querying GDC's file endpoint. The key of dict should 
            be GDC available fields for the file endpoint.
        projects: str or list, default None
            One project_id (str) or a list of project_id(s) of interest. 
            Default is to download data for all projects on GDC.
        work_dir: str, default '.'
            A directory structure will be built under "work_dir" according to 
            "project_id(s)" and "dataset_type":
                work_dir
                └── project_id(s) joint by "_"
                    └── dataset_type value(s) joint by "_" with whitespace 
                        removed
                        ├── data1
                        ├── data2
                        ├── ...
                        └── dataN
    
    Returns:
        data_dict: a dict with on entry. The key is the directory path for the 
        aggregated dataset and the valuse is a list of abspaths for downloaded 
        files belonging to the dataset.
    """
    
    sampleid_field='cases.samples.submitter_id'
    smapleid_processing = lambda x: x[:-1]
    
    if projects is None:
        projects = gdc.get_all_project_info().keys()
    if isinstance(projects, str):
        projects = [projects]
    total_projects = len(projects)
    project_count = 0
    
    work_dir = os.path.abspath(work_dir)
    project_dir = os.path.join(work_dir, '_'.join(projects))
    mkdir_p(project_dir)
    dataset_dir = '_'.join(sorted(dataset_type.values())).replace(' ', '')
    dataset_dir = os.path.join(project_dir, dataset_dir)
    mkdir_p(dataset_dir)
    
    print('Datasets will be downloaded to "{}".'.format(dataset_dir))
    dataset_description = ' - '.join(sorted(dataset_type.values()))
    dataset_status = 'Get "{}" dataset for project "[{}/{:d}] {}".'
    dataset_status = dataset_status.format(dataset_description, '{}', 
                                           total_projects, '{}')

    for project_id in projects:
        project_count += 1
        print(dataset_status.format(project_count, project_id))
        query_dict = {'cases.project.project_id': project_id,
                      'access': 'open'}
        query_dict.update(dataset_type)
        query_filter = gdc.and_eq_filter_constructor(query_dict)
        file_dict = gdc.get_file_dict(query_filter, 
                                      label_field=sampleid_field,
                                      label_func=smapleid_processing)
        if not file_dict:
            message = 'No {} data for "{}" project.'
            print(message.format(dataset_description, project_id))
            continue
        data_dict = {dataset_dir: []}
        total_files = len(file_dict)
        file_count = 0
        for uuid in file_dict:
            file_count += 1
            download_status = '[{:d}/{:d}] Download "{}": '
            print(download_status.format(file_count, total_files, 
                                         file_dict[uuid]), end='')
            file_path = os.path.join(dataset_dir, file_dict[uuid])
            if gdc.download_data(uuid, file_path):
                data_dict[dataset_dir].append(file_path)
    return data_dict

def read_by_ext(filename, mode='r'):
    """Automatically decide how to open a file which might be compressed.
    
    Uses the code from the "hook_compressed" function in the fileinput module.
    
    Args:
        filename: Must contain proper extension indicating the compression 
            condition.
        mode: Default is 'r'.
    
    Returns:
        A filehandle to be used with `with`.
    """
    
    ext = os.path.splitext(filename)[1]
    if ext == '.gz':
        import gzip
        return gzip.open(filename, mode)
    elif ext == '.bz2':
        import bz2
        return bz2.BZ2File(filename, mode)
    else:
        return open(filename, mode)

def xena_matrix_merge(file_list, read_table_kwargs={}, merge_axis=1, 
                      average=True, log_transform=True):
    total_files = len(file_list)
    file_count = 0
    df_list = []
    for file_path in file_list:
        file_count = file_count + 1
        print('\rProcessing {}/{} file...'.format(file_count, total_files), 
              end='')
        sys.stdout.flush()
        with read_by_ext(file_path) as f:
            filename = os.path.basename(file_path)
            sample_id = filename.split('.', 1)[0]
            pd_kwargs = {'index_col': 0, 'usecols': [0, 1]}
            pd_kwargs.update(read_table_kwargs)
            raw_df = pd.read_table(f, **pd_kwargs)
            raw_df.index.name = 'sampleid'
            if merge_axis == 1:
                raw_df.columns = [sample_id]
            elif merge_axis == 0:
                raw_df.set_index([[sample_id] * raw_df.shape[0]], inplace=True)
            else:
                error_msg = 'Unrecognized merge_axis: {}'.format(merge_axis)
                raise ValueError(error_msg)
            df_list.append(raw_df)
    print('\rAll {} files have been processed. '.format(total_files))
    print('Merging into one matrix ...')
    xena_matrix = pd.concat(df_list, axis=merge_axis)
    if average:
        print('Averaging duplicated samples ...')
        xena_matrix = xena_matrix.groupby(xena_matrix.columns, axis=1).mean()
    if log_transform:
        print('Transforming data matrix ...')
        xena_matrix = np.log2(xena_matrix + 1)
    print('Xena matrix ready.')
    return xena_matrix

def xena_maf_transform(input_path):
    if isinstance(input_path, list):
        if len(input_path) == 1:
            input_file = input_path[0]
        else:
            msg = ('MAF file should be transformed one at a time. {} files '
                + 'found in the input list')
            raise ValueError(msg)
    elif isinstance(input_path, str):
        input_file = input_path
    else:
        msg = 'Unrecognized type for input_path: {}'.format(type(input_path))
        raise ValueError(msg)

    print('Reading {} ...'.format(input_file))
    with read_by_ext(input_file) as f:
        xena_matrix = pd.read_table(
                f, 
                header=0, 
                usecols=[12, 36, 4, 5, 6, 39, 41, 50, 0, 10, 15, 109], 
                comment='#'
            )
    rename_dict = {'Hugo_Symbol': 'gene', 
                   'Chromosome': 'chrom', 
                   'Start_Position': 'chromstart', 
                   'End_Position': 'chromend', 
                   'Reference_Allele': 'ref', 
                   'Tumor_Seq_Allele2': 'alt', 
                   'Tumor_Sample_Barcode': 'sampleid', 
                   'HGVSp_Short': 'Amino_Acid_Change', 
                   'Consequence': 'effect',
                   'FILTER': 'filter'}
    xena_matrix['dna_vaf'] = (xena_matrix['t_alt_count'] 
                           / xena_matrix['t_depth'])
    xena_matrix = (
            xena_matrix.drop(['t_alt_count', 't_depth'], axis=1)
                       .rename(columns=rename_dict)
                       .set_index('sampleid', drop=False)
        )
    print('Xena matrix ready.')
    return xena_matrix

def import_gdc(project, dataset_type, work_path='.', matrix_dir=None,
               mode='both'):
    """An interface to the main pipeline for importing GDC RNA-seq data into 
    Xena.
    
    Args:
        project: str
            One project_id for a GDC project.
        dataset_type: str in ['rna_counts', 'rna_fpkm', 'rna_fpkmuq', 'mirna', 
            'mirna_isoform', 'masked_cnv']
        work_path: str, default '.'
            For 'both' or 'download' mode, it will be used for building the 
            download directory structure and saving download files. For 
            'transform' mode, all files under this directory will be treated 
            as data and used for building Xena compatible matrix. For 
            'transform' mode, work_path can also point to a file to be 
            transformed rather than a directory.
        matrix_dir: str, default None
            One directory to save the final Xena matrix.
        mode: str in ['both', 'download', 'transform'], default 'both'
            Action(s) to take for importing data. For just 'download' data, a 
            directory structure will be built under "work_path" according to 
            "projects" and "data_type_list". Data will be downloaded to 
            corresponding directories. For just 'transform' data, all files 
            under "data_dir" will be treated as one set of data. These data 
            will be transformed into one single Xena compatible matrix and 
            saved under work_dir. When 'both' actions are performed, 
            downloaded data will be automatically organized into data sets 
            according to "projects" and "data_type_list". Every data set will 
            be transformed into one single Xena compatible matrix and saved 
            together with its corresponding data. "data_dir" will be ignored 
            under 'both' mode.
    """
    
    if dataset_type not in ['rna_counts', 'rna_fpkm', 'rna_fpkmuq', 'mirna',
                            'mirna_isoform', 'masked_cnv','muse', 'mutect2', 
                            'somaticsniper', 'varscan2']:
        raise ValueError('Unrecognized dataset_type: {}.'.format(dataset_type))

    gdc_rna_counts = {'analysis.workflow_type': 'HTSeq - Counts'}
    gdc_rna_fpkm = {'analysis.workflow_type': 'HTSeq - FPKM'}
    gdc_rna_fpkmuq = {'analysis.workflow_type': 'HTSeq - FPKM-UQ'}
    gdc_mirna = {'data_type': 'miRNA Expression Quantification'}
    gdc_mirna_isoform = {'data_type': 'Isoform Expression Quantification'}
    gdc_masked_cnv = {'data_type': 'Masked Copy Number Segment'}
    gdc_snv_muse = {
            'analysis.workflow_type': 
                'MuSE Variant Aggregation and Masking'
        }
    gdc_snv_mutect2 = {
            'analysis.workflow_type': 
                'MuTect2 Variant Aggregation and Masking'
        }
    gdc_snv_somaticsniper = {
            'analysis.workflow_type': 
                'SomaticSniper Variant Aggregation and Masking'
        }
    gdc_snv_varscan2 = {
            'analysis.workflow_type': 
                'VarScan2 Variant Aggregation and Masking'
        }

    trans_args_rna = {'read_table_kwargs': {'header': None, 'comment': '_'}}
    trans_args_mirna = {'read_table_kwargs': {'header': 0, 'usecols': [0, 2]}}
    trans_args_mirna_isoform = {'read_table_kwargs': {'header': 0, 
                                                      'usecols': [1, 3]}}
    trans_args_masked_cnv = {'read_table_kwargs': {'header': 0,
                                                   'index_col': None,
                                                   'usecols': [1, 2, 3, 5]},
                             'merge_axis': 0,
                             'average': False,
                             'log_transform': False}

    dataset_map = {'rna_counts': {'gdc_type': gdc_rna_counts,
                                  'trans_func': xena_matrix_merge,
                                  'trans_args': trans_args_rna,
                                  'matrix_name': '{}.htseq.counts.tsv'},
                   'rna_fpkm': {'gdc_type': gdc_rna_fpkm,
                                'trans_func': xena_matrix_merge,
                                'trans_args': trans_args_rna,
                                'matrix_name': '{}.htseq.fpkm.tsv'},
                   'rna_fpkmuq': {'gdc_type': gdc_rna_fpkmuq,
                                  'trans_func': xena_matrix_merge,
                                  'trans_args': trans_args_rna,
                                  'matrix_name': '{}.htseq.fpkm-uq.tsv'},
                   'mirna': {'gdc_type': gdc_mirna,
                             'trans_func': xena_matrix_merge,
                             'trans_args': trans_args_mirna,
                             'matrix_name': '{}.mirna.tsv'},
                   'mirna_isoform': {'gdc_type': gdc_mirna_isoform,
                                     'trans_func': xena_matrix_merge,
                                     'trans_args': trans_args_mirna_isoform,
                                     'matrix_name': '{}.mirna.isoform.tsv'}, 
                   'masked_cnv': {'gdc_type': gdc_masked_cnv,
                                  'trans_func': xena_matrix_merge,
                                  'trans_args': trans_args_masked_cnv,
                                  'matrix_name': '{}.masked.cnv.tsv'},
                   'muse': {'gdc_type': gdc_snv_muse,
                            'trans_func': xena_maf_transform,
                            'trans_args': {},
                            'matrix_name': '{}.muse.snv.tsv'},
                   'mutect2': {'gdc_type': gdc_snv_mutect2,
                               'trans_func': xena_maf_transform,
                               'trans_args': {},
                               'matrix_name': '{}.mutect2.snv.tsv'},
                   'somaticsniper': {
                           'gdc_type': gdc_snv_somaticsniper,
                           'trans_func': xena_maf_transform,
                           'trans_args': {},
                           'matrix_name': '{}.somaticsniper.snv.tsv'
                        },
                   'varscan2': {'gdc_type': gdc_snv_varscan2,
                                'trans_func': xena_maf_transform,
                                'trans_args': {},
                                'matrix_name': '{}.varscan2.snv.tsv'}}


    if mode == 'both' or mode == 'download':
        gdc_type = dataset_map[dataset_type]['gdc_type']
        dataset = download_dataset(projects=project, dataset_type=gdc_type,
                                   work_dir=work_path)
        file_list = dataset.values()[0]

    if mode == 'transform':
        if os.path.isdir(work_path):
            file_list = []
            for f in os.listdir(work_path):
                file_list.append(os.path.join(work_path, f))
        else:
            file_list = [work_path]
            work_path = os.path.dirname(work_path)
    
    if mode == 'both' or mode == 'transform':
        trans_func = dataset_map[dataset_type]['trans_func']
        trans_args = dataset_map[dataset_type]['trans_args']
        xena_matrix = trans_func(file_list, **trans_args)
        if matrix_dir is None:
            matrix_dir = work_path
        matrix_name = dataset_map[dataset_type]['matrix_name']
        matrix_path = os.path.join(matrix_dir, matrix_name.format(project))
        print('Saving matrix to {} ...'.format(matrix_path))
        xena_matrix.to_csv(matrix_path, sep='\t')

def render_rna_counts_metadata(matrix_dir, matrix_name, keywords):
    """Make "metadata.json" for Xena importing
    
    Args:
        matrix_dir: str
            Directory for corresponding data matrix. Generated metadata file 
            will be saved in the same directory.
        matrix_name: str
            Filename for corresponding data matrix. Its extension will be 
            changed into '.json' and then used as the name for generated 
            metadata file.
        keywords: dict
            Must contain the following keywords in the template for proper 
            rendering:
                program
                project
                data_url
                tissue
    
    Returns:
        metadata: JSON formatted string. ready to be written into a file.
    """

    metadata_name = os.path.splitext(matrix_name)[0]+'.json'
    metadata_path = os.path.join(matrix_dir, metadata_name)
    keywords.update({'data_type': 'HTSeq - Counts', 
                     'date': datetime.date.today().isoformat()})
    
    template_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 
                                 'Resources')
    jinja2_env = jinja2.Environment(
            loader=jinja2.FileSystemLoader(template_dir))
    template = jinja2_env.get_template('template.rna.meta.json')
    with open(metadata_path, 'w') as f:
        f.write(template.render(**keywords))

def main():
    print('A python module of Xena specific importing pipeline for GDC data.')

#    print(gdc.get_all_project_info())

#    rna_counts_metadata = {'program': 'TCGA',
#                           'project': 'TCGA-BRCA',
#                           'data_url': 'https://api.gdc.cancer.gov/data/',
#                           'tissue': 'Breast'}
#    render_rna_counts_metadata('gitignore', 'test.uuid.tsv', 
#                               rna_counts_metadata)

if __name__ == '__main__':
    main()
