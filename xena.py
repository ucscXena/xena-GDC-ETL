#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This module provides functions/wrappers necessary for quickly assembling an 
ETL pipeline importing GDC data into Xena.
"""

# Ensure Python 2 and 3 compatibility
from __future__ import division
from __future__ import print_function

import time
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

def download_dataset(dataset_type, projects=None, download_dir='.', 
                     label_field='cases.samples.submitter_id'):
    """Download GDC's open access data.
    
    Data was selected based on "project_id(s)" and "dataset_type". As long as 
    "project_id(s)" and "dataset_type" are valid, retrieved data will be put 
    together under one directory as a single dataset. Tough this function 
    won't check, putting different types of data into one dataset is NOT 
    recommended. By default, data files will be renamed into 
    "<cases.samples.submitter_id>.<UUID>.<file extension>" and saved under 
    current python work directory.
    
    Args:
        dataset_type: dict
            A dict specifying dataset(s) of interest on GDC. It will be used 
            directly for querying GDC's file endpoint. Key(s) of dict should 
            be GDC available fields for the file endpoint. See details at
            https://docs.gdc.cancer.gov/API/Users_Guide/Appendix_A_Available_Fields/#file-fields
        projects: str or list, default None
            One project_id (str) or a list of project_id(s) of interest. If a 
            list of project_id is provided, those projects will be aggregated 
            into one cohort which will be named as "_".join(projects). 
            Default is to download the specific type ("dataset_type") of data 
            for all projects on GDC as one single cohort, which will be named 
            using today's date as "All_GDC_projects_%m-%d-%Y".
        download_dir: str, default '.'
            A directory for saving downloaded data. Default "download_dir" is 
            the current working directory of the python script. It is 
            recommended that the "download_dir" is an empty directory 
            dedicated for the new cohort of dataset. 
        label_field:  str, default 'cases.samples.submitter_id'
            A single GDC available file field whose value will be 
            used for renaming downloaded files. By default, Xena uses 
            'cases.samples.submitter_id' as sample ID. Therefore filenames for 
            downloaded files will be "submitter_id.UUID.file_extension".
    
    Returns:
        data_dict: a dict with one entry. The key is the directory path for 
        this new (aggregated) dataset and the valuse is a list of abspaths for 
        downloaded files.
    """
    
    if isinstance(projects, str):
        projects = [projects]
    print('Datasets will be downloaded to\n"{}".'.format(download_dir))
    query_dict = {'access': 'open'}
    if projects is not None:
        query_dict['cases.project.project_id'] = projects
    query_dict.update(dataset_type)
    query_filter = gdc.and_in_filter_constructor(query_dict)
    file_dict = gdc.get_file_dict(query_filter, label_field=label_field)
    if not file_dict:
        message = 'No {} data found for project {}.'
        dataset_description = ' - '.join(sorted(dataset_type.values()))
        # raise ValueError(message.format(dataset_description, projects))
        print(message.format(dataset_description, projects))
        return {}
    data_dict = {download_dir: []}
    total_files = len(file_dict)
    file_count = 0
    for uuid in file_dict:
        file_count += 1
        download_status = '[{:d}/{:d}] Download "{}": '
        print(download_status.format(file_count, total_files, 
                                     file_dict[uuid]), end='')
        file_path = os.path.join(download_dir, file_dict[uuid])
        if gdc.download_data(uuid, file_path):
            data_dict[download_dir].append(file_path)
    return data_dict

def read_by_ext(filename, mode='r'):
    """Automatically decide how to open a file which might be compressed.
    
    Uses the code from the "hook_compressed" function in python's fileinput 
    module.
    
    Args:
        filename: str
            Must contain proper extension indicating the compression condition.
        mode: str, default 'r'
    
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

def process_average_log(df):
    """Process Xena data matrix by first averaging columns having the same 
    name and then transform it by log(x + 1)
    
    Args:
        df: pandas.DataFrame
    
    Returns:
        df_transformed: Transformed pandas.DataFrame
    """
    
    print('Averaging duplicated samples ...')
    df_avg = df.groupby(df.columns, axis=1).mean()
    print('Log transforming data ...')
    return np.log2(df_avg + 1)

def process_maf(df):
    """ Transform pre-sliced GDC's MAF data into Xena data matrix.
    
    A new column of DNA variant allele frequencies named "dna_vaf" will 
    calculated by division "t_alt_count" / "t_depth". Columns "t_alt_count" 
    and "t_depth" will then be dropped. At last column will be renamed 
    accordingly and row index will be set as sample ID.
    
    Args:
        df: pandas.DataFrame
    
    Returns:
        df_transformed: Transformed pandas.DataFrame
    """
    print('Calculating "dna_vaf" ...')
    df['dna_vaf'] = df['t_alt_count'] / df['t_depth']
    print('Trim "Tumor_Sample_Barcode" into Xena sample ID ...')
    trim_func = lambda x: '-'.join(x.split('-', 4)[0:4])
    df['Tumor_Sample_Barcode'] = df['Tumor_Sample_Barcode'].apply(trim_func)
    print('Re-organizing matrix ...')
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
    return (
            df.drop(['t_alt_count', 't_depth'], axis=1)
              .rename(columns=rename_dict)
              .set_index('sampleid')
        )

def xena_matrix_merge(file_list, read_table_kwargs={}, merge_axis=1, 
                      matrix_process=None, index_name=None):
    """Transform a GDC dataset into Xena matrix
    
    Args:
        file_list: list
            A list of pathes pointing to data file(s) of a GDC dataset which 
            will be assembled and/or transformed into Xena matrix
        read_table_kwargs: dict, default {}
            Keyword arguments which will be passed to pandas.read_table for 
            reading each GDC data file as a pandas.DataFrame. By default, the 
            first and second columns of the GDC data file will be read in 
            ('usecols': [0, 1]) and the first column will be used as the row 
            index ('index_col': 0).
        merge_axis: {0/’index’, 1/’columns’}, default 1
            The axis argument which will be passed to pandas.concat for 
            assembling multiple GDC data files into one pandas.DataFrame.
        matrix_process: function, default None
            A function which takes in only one required argument, 
            pandas.DataFrame, and returns one pandas.DataFrame. It is used for 
            further processing the Xena matrix after it is sliced and/or 
            assembled from GDC data. By default, no further processing will 
            happen.
        index_name: str, default None
            Name used for renaming row/column index (i.e. cell(1,1)). By 
            default, the index name will be set by pandas.read_table and left 
            as it is.
    
    Returns:
        xena_matrix: assembled/transformed pandas.DataFrame
    """
    
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
            if merge_axis == 1:
                pd_kwargs = {'index_col': 0, 'usecols': [0, 1]}
            elif merge_axis == 0:
                pd_kwargs = {}
            else:
                error_msg = 'Unrecognized merge_axis: {}'.format(merge_axis)
                raise ValueError(error_msg)
            pd_kwargs.update(read_table_kwargs)
            raw_df = pd.read_table(f, **pd_kwargs)
            if merge_axis == 1:
                raw_df.columns = [sample_id]
            elif merge_axis == 0:
                raw_df.set_index([[sample_id] * raw_df.shape[0]], inplace=True)
            df_list.append(raw_df)
    print('\rAll {} files have been processed. '.format(total_files))
    print('Merging into one matrix ...')
    xena_matrix = pd.concat(df_list, axis=merge_axis)
    if matrix_process is not None:
        xena_matrix = matrix_process(xena_matrix)
    if index_name is not None:
        xena_matrix.index.name = index_name
    print('Xena matrix ready.')
    return xena_matrix

def render_metadata(matrix_path, project=None, dataset_type=None, 
                    probemap=None, keywords={}):
    """Make "metadata.json" for Xena importing
    
    Args:
        matrix_path: str
            Path to the file of Xena matrix. Generated metadata file 
            will be saved in the same directory, with a ".json" postfix 
            appended to the filename of Xena matrix. If project and/or 
            dataset_type is None, project and/or dataset_type info will be 
            inferrd from the filename of Xena matrix, assuming the filename is 
            in the format of "project_id.dtype_string.tsv".
        project: str
            One GDC project_id.
        dataset_type: str in ['rna_counts', 'rna_fpkm', 'rna_fpkmuq', 'mirna', 
            'mirna_isoform', 'cnv', 'masked_cnv', 'muse', 'mutect2', 
            'somaticsniper', 'varscan2']
        probemap: str or dict
            A string of path pointing to the probeMap file corresponding to 
            the dataset_type. A probeMap can be sepecified in metadata for 
            'rna_counts', 'rna_fpkm', 'rna_fpkmuq', and 'mirna'. If a dict is 
            provided, it can have 2 entries, 'rna' or 'mirna': probemap['rna'] 
            is the path of probeMap for dataset_type 'rna_counts', 'rna_fpkm' 
            and 'rna_fpkmuq'; probemap['mirna'] is the path of probeMap for 
            dataset_type 'mirna'. The pointed probeMap file won't be verified.
        keywords: dict
            Optional keywords to be rendered from jinja2 templates. Keywords 
            set by this dict have priority and won't be override. Supported 
            keywords are "xena_cohort", "project_id", "gdc_type", 
            "date", "notes", "probemap_url", "maf_uuid".
    
    Returns:
        metadata: JSON formatted string. ready to be written into a file.
    """
    
    # Map GDC project_id to Xena specific cohort name.
    cohort_map = {'TCGA-BRCA': 'TCGA Breast Cancer (BRCA)',
                  'TCGA-LUAD': 'TCGA Lung Adenocarcinoma (LUAD)',
                  'TCGA-UCEC': 'TCGA Endometrioid Cancer (UCEC)',
                  'TCGA-LGG': 'TCGA Lower Grade Glioma (LGG)',
                  'TCGA-HNSC': 'TCGA Head and Neck Cancer (HNSC)',
                  'TCGA-PRAD': 'TCGA Prostate Cancer (PRAD)',
                  'TCGA-LUSC': 'TCGA Lung Squamous Cell Carcinoma (LUSC)',
                  'TCGA-THCA': 'TCGA Thyroid Cancer (THCA)',
                  'TCGA-SKCM': 'TCGA Melanoma (SKCM)',
                  'TCGA-OV': 'TCGA Ovarian Cancer (OV)',
                  'TCGA-STAD': 'TCGA Stomach Cancer (STAD)',
                  'TCGA-COAD': 'TCGA Colon Cancer (COAD)',
                  'TCGA-BLCA': 'TCGA Bladder Cancer (BLCA)',
                  'TCGA-GBM': 'TCGA Glioblastoma (GBM)',
                  'TCGA-LIHC': 'TCGA Liver Cancer (LIHC)',
                  'TCGA-KIRC': 'TCGA Kidney Clear Cell Carcinoma (KIRC)',
                  'TCGA-CESC': 'TCGA Cervical Cancer (CESC)',
                  'TCGA-KIRP': 'TCGA Kidney Papillary Cell Carcinoma (KIRP)',
                  'TCGA-SARC': 'TCGA Sarcoma (SARC)',
                  'TCGA-ESCA': 'TCGA Esophageal Cancer (ESCA)',
                  'TCGA-PAAD': 'TCGA Pancreatic Cancer (PAAD)',
                  'TCGA-PCPG': 'TCGA Pheochromocytoma & Paraganglioma (PCPG)',
                  'TCGA-READ': 'TCGA Rectal Cancer (READ)',
                  'TCGA-TGCT': 'TCGA Testicular Cancer (TGCT)',
                  'TCGA-LAML': 'TCGA Acute Myeloid Leukemia (LAML)',
                  'TCGA-THYM': 'TCGA Thymoma (THYM)',
                  'TCGA-ACC': 'TCGA Adrenocortical Cancer (ACC)',
                  'TCGA-MESO': 'TCGA Mesothelioma (MESO)',
                  'TCGA-UVM': 'TCGA Ocular melanomas (UVM)',
                  'TCGA-KICH': 'TCGA Kidney Chromophobe (KICH)',
                  'TCGA-UCS': 'TCGA Uterine Carcinosarcoma (UCS)',
                  'TCGA-CHOL': 'TCGA Bile Duct Cancer (CHOL)',
                  'TCGA-DLBC': 'TCGA Large B-cell Lymphoma (DLBC)'}
    # Map for inferring dataset_type from matrix filename.
    nominal_dtype_map = {'htseq.counts.tsv': 'rna_counts', 
                         'htseq.fpkm.tsv': 'rna_fpkm', 
                         'htseq.fpkm-uq.tsv': 'rna_fpkmuq', 
                         'mirna.tsv': 'mirna', 
                         'mirna.isoform.tsv': 'mirna_isoform', 
                         'cnv.tsv': 'cnv', 
                         'masked.cnv.tsv': 'masked_cnv', 
                         'muse.snv.tsv': 'muse', 
                         'mutect2.snv.tsv': 'mutect2', 
                         'somaticsniper.snv.tsv': 'somaticsniper', 
                         'varscan2.snv.tsv': 'varscan2'}
    # Map dataset_map to GDC data type or GDC workflow type.
    gdc_type_map = {'rna_counts': 'HTSeq - Counts',  
                    'rna_fpkm': 'HTSeq - FPKM', 
                    'rna_fpkmuq': 'HTSeq - FPKM-UQ', 
                    'mirna': 'miRNA Expression Quantification', 
                    'mirna_isoform': 'Isoform Expression Quantification', 
                    'cnv': 'Copy Number Segment', 
                    'masked_cnv': 'Masked Copy Number Segment', 
                    'muse': 'MuSE Variant Aggregation and Masking', 
                    'mutect2': 'MuTect2 Variant Aggregation and Masking', 
                    'somaticsniper': 
                        'SomaticSniper Variant Aggregation and Masking', 
                    'varscan2': 'VarScan2 Variant Aggregation and Masking'}
    # Map dataset_type to corresponding metadata template.
    template_map = {'rna_counts': 'template.rna.meta.json', 
                    'rna_fpkm': 'template.rna.meta.json', 
                    'rna_fpkmuq': 'template.rna.meta.json', 
                    'mirna': 'template.mirna.meta.json', 
                    'mirna_isoform': 'template.mirna.isoform.meta.json', 
                    'cnv': 'template.cnv.meta.json', 
                    'masked_cnv': 'template.cnv.meta.json', 
                    'muse': 'template.snv.meta.json', 
                    'mutect2': 'template.snv.meta.json', 
                    'somaticsniper': 'template.snv.meta.json', 
                    'varscan2': 'template.snv.meta.json'}

    # Start to generate metadata.
    if not os.path.isfile(matrix_path):
        raise ValueError('{} is not a file.'.format(matrix_path))

    print('Creating metadata file ...')
    matrix_path = os.path.abspath(matrix_path)
    matrix_date = time.strftime("%m-%d-%Y", 
                                time.gmtime(os.path.getmtime(matrix_path)))
    matrix_dir, matrix_filename = os.path.split(matrix_path)
    metadata_path = os.path.join(matrix_dir, matrix_filename + '.json')

    if (project is None) or (dataset_type is None):
        nominal_project, nominal_dtype = matrix_filename.split('.', 1)
        if project is None:
            project = nominal_project
        if dataset_type is None:
            dataset_type = nominal_dtype_map[nominal_dtype]

    if project in cohort_map:
        xena_cohort = cohort_map[project]
    else:
        xena_cohort = project
    gdc_type = gdc_type_map[dataset_type]
    variables = {'xena_cohort': xena_cohort, 
                 'project_id': project, 
                 'gdc_type': gdc_type, 
                 'date': matrix_date}
    if dataset_type in ['rna_counts', 'rna_fpkm', 'rna_fpkmuq', 'mirna']:
        if isinstance(probemap, str):
            probemap_url = os.path.relpath(probemap, 
                                           matrix_dir).replace('\\', '/')
            variables['probemap_url'] = probemap_url
        elif isinstance(probemap, dict):
            if dataset_type.startswith('rna') and 'rna' in probemap:
                probemap_url = os.path.relpath(probemap['rna'], 
                                               matrix_dir).replace('\\', '/')
                variables['probemap_url'] = probemap_url
            elif dataset_type.startswith('mirna') and 'mirna' in probemap:
                probemap_url = os.path.relpath(probemap['mirna'], 
                                               matrix_dir).replace('\\', '/')
                variables['probemap_url'] = probemap_url
    if dataset_type in ['muse', 'mutect2', 'somaticsniper', 'varscan2']:
        query_dict = {'access': 'open',
                      'cases.project.project_id': project,
                      'analysis.workflow_type': gdc_type}
        query_filter = gdc.and_in_filter_constructor(query_dict)
        file_dict = gdc.get_file_dict(query_filter)
        if len(file_dict) == 1:
            variables['maf_uuid'] = file_dict.keys()[0]
    variables.update(keywords)

    template_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                'Resources')
    jinja2_env = jinja2.Environment(
            loader=jinja2.FileSystemLoader(template_dir))
    template_json = template_map[dataset_type]
    template = jinja2_env.get_template(template_json)
    with open(metadata_path, 'w') as f:
        f.write(template.render(**variables))
    print('Metadata ready.')
    return metadata_path

def import_gdc(project, dataset_type, work_path='.', matrix_path=None, 
               probemap_path=None, mode='all'):
    """An wrapper for building a main pipeline to import GDC RNA-seq data into 
    Xena.
    
    The default directory structure is like this:
        work_dir
        └── project
            ├── "GDC_Raw_Data"
            │   └── dataset_type value(s) joint by "_" with whitespace 
            │       removed
            │       ├── data1
            │       ├── data2
            │       ├── ...
            │       └── dataN
            └── "Xena_Matrices"
                ├── project_id.dataset_type_string1.tsv
                ├── project_id.dataset_type_string1.tsv.json
                ├── project_id.dataset_type_string2.tsv
                ├── project_id.dataset_type_string2.tsv.json
                ├── ...
                ├── project_id.dataset_type_stringN.tsv
                └── project_id.dataset_type_stringN.tsv.json
    
    Args:
        project: str
            One GDC project_id.
        dataset_type: str in ['rna_counts', 'rna_fpkm', 'rna_fpkmuq', 'mirna', 
            'mirna_isoform', 'cnv', 'masked_cnv', 'muse', 'mutect2', 
            'somaticsniper', 'varscan2']
        work_path: str, default '.'
            For 'all' or 'download' mode, it will be used for building the 
            default directory structure shown above to save download files. 
            If "work_path" is not a directory, its directory 
            (os.path.dirname(work_path)) will be used. 
            For 'transform' mode, if "work_path" is a directory, all files 
            under this directory will be treated as data and used for building 
            Xena compatible matrix. If "work_path" is not a directory, it will 
            be assumed to be a file and this single file will be tranformed. 
            For 'metadata' mode, "work_path" will be ignored.
        matrix_path: str, default None
            For 'download' mode, "matrix_path" will be ignored.
            If "matrix_path" is None:
                For 'all' mode, the Xena matrix will be saved under 
                "Xena_Matrices" of the default directory structure shown above.
                For 'transform' mode, the default directory structure is 
                assumed and inferred by looking for "project" in "work_path". 
                If the "project" is found, the Xena matrix will be saved under 
                "Xena_Matrices" of the default directory structure shown above.
                If "project" is not found, i.e. the default directory 
                structure doesn't exist, the Xena matrix will be saved under 
                "work_path" directory or the directory of "work_path" 
                (os.path.dirname(work_path)).
                For 'metadata' mode, a NoneType "matrix_path" will trigger an 
                error.
            If "matrix_path" is a directory:
                For 'all' and 'transform' mode, the Xena matrix will be saved 
                under this "matrix_path" directory.
                For 'metadata' mode, "matrix_path" being a directory will 
                trigger an error in this "import_gdc" wrapper.
            If "matrix_path" is not a directory:
                For 'all' and 'transform' mode, the Xena matrix will be saved 
                as "matrix_path". 
                For 'metadata' mode, the "matrix_path" file will be treated as 
                a Xena matrix and the corresponding metadata will be generated.
        probemap_path: str, default None
            Path to one probeMap file corresponding to the dataset_type. For 
            'download' and 'transform' mode, "probemap_path" will be ignored.
        mode: str in ['all', 'download', 'transform', 'metadata'], default 
              'all'
            Action(s) to be taken for the data importing pipeline. 
            For 'all' mode, a default directory structure (shown above) will 
            be created under the "work_path". Both "project" and 
            "dataset_type" are required for defining a single a dataset which 
            will be downloaded to one corresponding directory in the default 
            directory structure. Then the dataset will be transformed into one 
            Xena matrix and saved accordingly, depending on the "matrix_path". 
            Finally, one metadata will be generated for the Xena matrix.
            For 'download' mode, a default directory structure (shown above) 
            will be created under the "work_path". Both "project" and 
            "dataset_type" are required for defining a single a dataset which 
            will be downloaded to one corresponding directory in the default 
            directory structure.
            For 'transform' mode, all data file(s) determined according to the 
            "work_path" will be treated as one dataset, which will be 
            transformed into one Xena matrix and saved accordingly, depending 
            on the "matrix_path".
            For 'metadata' mode, metadata will generated for the single Xena 
            matrix defined by the "matrix_path".
    """

    # Download settings for each type of dataset.
    download_args_rna_counts = {
            'dataset_type': {'analysis.workflow_type': 'HTSeq - Counts'}
        }
    download_args_rna_fpkm = {
            'dataset_type': {'analysis.workflow_type': 'HTSeq - FPKM'}
        }
    download_args_rna_fpkmuq = {
            'dataset_type': {'analysis.workflow_type': 'HTSeq - FPKM-UQ'}
        }
    download_args_mirna = {
            'dataset_type': {'data_type': 'miRNA Expression Quantification'}
        }
    download_args_mirna_isoform = {
            'dataset_type': {'data_type': 'Isoform Expression Quantification'}
        }
    download_args_cnv = {
            'dataset_type': {'data_type': 'Copy Number Segment'}
        }
    download_args_masked_cnv = {
            'dataset_type': {'data_type': 'Masked Copy Number Segment'}
        }
    download_args_snv_muse = {
            'dataset_type': {
                    'analysis.workflow_type': 
                        'MuSE Variant Aggregation and Masking'
                },
            'label_field': 'cases.project.project_id'
        }
    download_args_snv_mutect2 = {
            'dataset_type': {
                    'analysis.workflow_type': 
                        'MuTect2 Variant Aggregation and Masking'
                },
            'label_field': 'cases.project.project_id'
        }
    download_args_snv_somaticsniper = {
            'dataset_type': {
                    'analysis.workflow_type': 
                        'SomaticSniper Variant Aggregation and Masking'
                },
            'label_field': 'cases.project.project_id'
        }
    download_args_snv_varscan2 = {
            'dataset_type': {
                    'analysis.workflow_type': 
                        'VarScan2 Variant Aggregation and Masking'
                },
            'label_field': 'cases.project.project_id'
        }
    # Xena transformation settings for each type of dataset.
    trans_args_rna = {'read_table_kwargs': {'header': None, 'comment': '_'},
                      'matrix_process': process_average_log,
                      'index_name': 'Ensembl_ID'}
    trans_args_mirna = {'read_table_kwargs': {'header': 0, 'usecols': [0, 2]},
                        'matrix_process': process_average_log}
    trans_args_mirna_isoform = {'read_table_kwargs': {'header': 0, 
                                                      'usecols': [1, 3]},
                                'matrix_process': process_average_log}
    trans_args_cnv = {'read_table_kwargs': {'header': 0, 
                                            'usecols': [1, 2, 3, 5]},
                             'merge_axis': 0}
    trans_args_snv = {'read_table_kwargs': {'header': 0,
                                            'usecols': [12, 36, 4, 5, 6, 39, 
                                                        41, 51, 0, 10, 15, 
                                                        110], 
                                            'comment': '#'},
                      'merge_axis': 0,
                      'matrix_process': process_maf,
                      'index_name': 'Sample_ID'}
    # Map for dataset specific download settings, transformation settings, and 
    # Xena matrix name boilerplates.
    dataset_map = {
            'rna_counts': {
                    'download_args': download_args_rna_counts,
                    'trans_args': trans_args_rna,
                    'matrix_name': '{}.htseq.counts.tsv'
                },
            'rna_fpkm': {
                    'download_args': download_args_rna_fpkm,
                    'trans_args': trans_args_rna,
                    'matrix_name': '{}.htseq.fpkm.tsv'
                },
            'rna_fpkmuq': {
                    'download_args': download_args_rna_fpkmuq,
                    'trans_args': trans_args_rna,
                    'matrix_name': '{}.htseq.fpkm-uq.tsv'
                },
            'mirna': {
                    'download_args': download_args_mirna,
                    'trans_args': trans_args_mirna,
                    'matrix_name': '{}.mirna.tsv'
                },
            'mirna_isoform': {
                    'download_args': download_args_mirna_isoform,
                    'trans_args': trans_args_mirna_isoform,
                    'matrix_name': '{}.mirna.isoform.tsv'
                },
            'cnv': {
                    'download_args': download_args_cnv,
                    'trans_args': trans_args_cnv,
                    'matrix_name': '{}.cnv.tsv'
                },
            'masked_cnv': {
                    'download_args': download_args_masked_cnv,
                    'trans_args': trans_args_cnv,
                    'matrix_name': '{}.masked.cnv.tsv'
                },
            'muse': {
                    'download_args': download_args_snv_muse,
                    'trans_args': trans_args_snv,
                    'matrix_name': '{}.muse.snv.tsv'
                },
            'mutect2': {
                    'download_args': download_args_snv_mutect2,
                    'trans_args': trans_args_snv,
                    'matrix_name': '{}.mutect2.snv.tsv'
                },
            'somaticsniper': {
                    'download_args': download_args_snv_somaticsniper,
                    'trans_args': trans_args_snv,
                    'matrix_name': '{}.somaticsniper.snv.tsv'
                },
            'varscan2': {
                    'download_args': download_args_snv_varscan2,
                    'trans_args': trans_args_snv,
                    'matrix_name': '{}.varscan2.snv.tsv'
                }
        }

    # Start importing process here.
    if dataset_type not in dataset_map:
        raise ValueError('Unrecognized dataset_type: {}.'.format(dataset_type))

    work_path = os.path.abspath(work_path)
    if os.path.isdir(work_path):
        work_dir = work_path
    else:
        work_dir = os.path.dirname(work_path)
    
    if mode == 'all' or mode == 'download':
        download_args = dataset_map[dataset_type]['download_args']
        gdc_type = download_args['dataset_type']
        print('Get {} dataset for {}.'.format(gdc_type.values(), project))
        dataset_dirname = gdc_type.values()[0].replace(' ', '')
        dataset_dir = os.path.join(work_dir, project, 'GDC_Raw_Data', 
                                   dataset_dirname)
        mkdir_p(dataset_dir)
        dataset = download_dataset(projects=project, download_dir=dataset_dir, 
                                   **download_args)
        if not dataset:
            return
        file_list = dataset.values()[0]

    if mode == 'transform':
        if os.path.isdir(work_path):
            file_list = []
            for f in os.listdir(work_dir):
                f_path = os.path.abspath(os.path.join(work_dir, f))
                if os.path.isfile(f_path):
                    file_list.append(f_path)
        else:
            file_list = [work_path]
    
    if mode == 'all' or mode == 'transform':
        trans_args = dataset_map[dataset_type]['trans_args']
        xena_matrix = xena_matrix_merge(file_list, **trans_args)
        if matrix_path is None:
            matrix_dir = os.path.dirname(file_list[0])
            while project in os.path.dirname(matrix_dir):
                matrix_dir = os.path.dirname(matrix_dir)
            matrix_dir = os.path.join(matrix_dir, 'Xena_Matrices')
        elif os.path.isdir(matrix_path):
            matrix_dir = matrix_path
        else:
            matrix_dir = os.path.dirname(matrix_path)
        mkdir_p(matrix_dir)
        if matrix_path is None or os.path.isdir(matrix_path):
            matrix_name = dataset_map[dataset_type]['matrix_name']
            matrix_path = os.path.join(matrix_dir, matrix_name.format(project))
        print('Saving matrix to {} ...'.format(matrix_path))
        xena_matrix.to_csv(matrix_path, sep='\t')
    
    if mode == 'all' or mode == 'metadata':
        if os.path.isfile(matrix_path):
            render_metadata(matrix_path, project=project,
                            dataset_type=dataset_type, probemap=probemap_path)
        else:
            raise ValueError('{} is not a Xena matrix.'.format(matrix_path))
    
    print('Data transformation is finished.')

def main():
    print('A python module of Xena specific importing pipeline for GDC data.')

if __name__ == '__main__':
    main()
