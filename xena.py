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

def download_dataset(dataset_type, projects=None, work_dir='.', 
                     label_field='cases.samples.submitter_id'):
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
            One project_id (str) or a list of project_id(s) of interest. If a 
            list of project_id is provided, those projects will be aggregated 
            into one cohort. Default is to download specific type 
            ("dataset_type") of data for all projects on GDC as one single 
            cohort.
        work_dir: str, default '.'
            A directory structure will be built under "work_dir" according to 
            "project_id(s)" and "dataset_type":
                work_dir
                └── project_id(s) joint by "_" or "All_GDC_projects_%m-%d-%Y"
                    └── GDC_Raw_Data
                        └── dataset_type value(s) joint by "_" with whitespace 
                            removed
                            ├── data1
                            ├── data2
                            ├── ...
                            └── dataN
            Default is the current working directory of the python script.
        label_field:  str, default 'cases.samples.submitter_id'
            A single GDC available file field whose value will be 
            used for renaming downloaded files. By default Xena uses 
            'cases.samples.submitter_id' as sample ID. Therefore filenames for 
            downloaded files are in the format of 
            "submitter_id.UUID.file_extension".
    
    Returns:
        data_dict: a dict with on entry. The key is the directory path for the 
        aggregated dataset and the valuse is a list of abspaths for downloaded 
        files belonging to the dataset.
    """
    
    dataset_description = ' - '.join(sorted(dataset_type.values()))
    if isinstance(projects, str):
        projects = [projects]
    work_dir = os.path.abspath(work_dir)
    if projects is None:
        cohort_dirname = datetime.date.today().strftime(
                'All_GDC_projects_%m-%d-%Y'
            )
        print('Get "{}" dataset for {}.'.format(dataset_description, 
                                                cohort_dirname))
    else:
        cohort_dirname = '_'.join(projects)
        print('Get "{}" dataset for {}.'.format(dataset_description, projects))
    dataset_dirname = '_'.join(sorted(dataset_type.values())).replace(' ', '')
    dataset_dir = os.path.join(work_dir, cohort_dirname, 'GDC_Raw_Data', 
                               dataset_dirname)
    mkdir_p(dataset_dir)
    print('Datasets will be downloaded to\n"{}".'.format(dataset_dir))
    query_dict = {'access': 'open'}
    if projects is not None:
        query_dict['cases.project.project_id'] = projects
    query_dict.update(dataset_type)
    query_filter = gdc.and_in_filter_constructor(query_dict)
    file_dict = gdc.get_file_dict(query_filter, label_field=label_field)
    if not file_dict:
        message = 'No {} data for project {}.'
        raise ValueError(message.format(dataset_description, projects))
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

def process_average_log(df):
    # The following if is for RNA-seq data
    # TODO: Move to settings for rna data type
    if df.index.name == 0:
        df.index.name = 'Ensemble_ID'
    print('Averaging duplicated samples ...')
    df_avg = (
            df.rename(columns=lambda x: x[:-1])
              .groupby(df.columns, axis=1)
              .mean()
        )
    print('Log transforming data ...')
    return np.log2(df_avg + 1)

def process_maf(df):
    df.index.name = 'Sample_ID'
    print('Calculating "dna_vaf" ...')
    df['dna_vaf'] = df['t_alt_count'] / df['t_depth']
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
    print('Re-organizing matrix ...')
    return (
            df.drop(['t_alt_count', 't_depth'], axis=1)
              .rename(columns=rename_dict)
              .set_index('sampleid', drop=False)
        )

def xena_matrix_merge(file_list, read_table_kwargs={}, merge_axis=1, 
                      matrix_process=None, average=True, log_transform=True):
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
#    if average:
#        print('Averaging duplicated samples ...')
#        xena_matrix = xena_matrix.groupby(xena_matrix.columns, axis=1).mean()
#    if log_transform:
#        print('Transforming data matrix ...')
#        xena_matrix = np.log2(xena_matrix + 1)
    print('Xena matrix ready.')
    return xena_matrix

def render_metadata(matrix_path, project=None, dataset_type=None, 
                    keywords = {}):
    """Make "metadata.json" for Xena importing
    
    Args:
        matrix_path: str
            Path to the file of Xena matrix. Generated metadata file 
            will be saved in the same directory, with a ".json" postfix 
            appended to the filename of Xena matrix. If project and/or 
            dataset_type is None, project and/or dataset_type info will be 
            extracted from the filename of Xena matrix.
        project: str
            GDC project_id.
        dataset_type: str in ['rna_counts', 'rna_fpkm', 'rna_fpkmuq', 'mirna', 
            'mirna_isoform', 'cnv', 'masked_cnv', 'muse', 'mutect2', 
            'somaticsniper', 'varscan2']
        keywords: dict
            Optional keywords to be rendered from jinja2 templates. Keywords 
            set by this dict have priority and won't be override. Supported 
            keywords are "xena_cohort", "project_id", "gdc_type", 
            "date", "notes", "maf_url".
    
    Returns:
        metadata: JSON formatted string. ready to be written into a file.
    """

    print('Creating metadata file ...')
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
    #       {{ maf_url}}
    
    if os.path.isdir(matrix_path):
        message = ('Require a real file; '
                +  '"matrix_path" now points to a directory: {}')
        raise ValueError(message.format(matrix_path))
    elif os.path.islink(matrix_path):
        message = ('Require a real file; '
                +  '"matrix_path" now points to a symbolic link: {}')
        raise ValueError(message.format(matrix_path))
    matrix_date = os.path.getmtime(matrix_path)
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
    template_json = template_map[dataset_type]
    variables = {'xena_cohort': xena_cohort, 
                 'project_id': project, 
                 'gdc_type': gdc_type, 
                 'date': matrix_date}
    variables.update(keywords)
    
    template_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 
                                 'Resources')
    jinja2_env = jinja2.Environment(
            loader=jinja2.FileSystemLoader(template_dir))
    template = jinja2_env.get_template(template_json)
    with open(metadata_path, 'w') as f:
        f.write(template.render(**variables))
    return metadata_path

def import_gdc(project, dataset_type, work_path='.', matrix_dir=None,
               mode='all'):
    """An interface to the main pipeline for importing GDC RNA-seq data into 
    Xena.
    
    Args:
        project: str
            One project_id for a GDC project.
        dataset_type: str in ['rna_counts', 'rna_fpkm', 'rna_fpkmuq', 'mirna', 
            'mirna_isoform', 'cnv', 'masked_cnv', 'muse', 'mutect2', 
            'somaticsniper', 'varscan2']
        work_path: str, default '.'
            For 'all' or 'download' mode, it will be used for building the 
            download directory structure and saving download files. If 
            "work_path" is not a directory, its directory 
            (os.path.dirname(work_path)) will be used. For 'transform' mode, 
            if "work_path" is a directory, all files under this directory will 
            be treated as data and used for building Xena compatible matrix. 
            If "work_path" is not a directory, it will be assumed to be a file 
            and this single file will be tranformed.
        matrix_dir: str, default None
            One directory to save the final Xena matrix.
        mode: str in ['all', 'download', 'transform'], default 'all'
            Action(s) to take for importing data. For just 'download' data, a 
            directory structure will be built under "work_path" according to 
            "projects" and "data_type_list". Data will be downloaded to 
            corresponding directories. For just 'transform' data, all files 
            under "data_dir" will be treated as one set of data. These data 
            will be transformed into one single Xena compatible matrix and 
            saved under work_dir. When 'all' actions are performed, 
            downloaded data will be automatically organized into data sets 
            according to "projects" and "data_type_list". Every data set will 
            be transformed into one single Xena compatible matrix and saved 
            together with its corresponding data. "data_dir" will be ignored 
            under 'all' mode.
    """
    
    if dataset_type not in ['rna_counts', 'rna_fpkm', 'rna_fpkmuq', 'mirna',
                            'mirna_isoform', 'cnv', 'masked_cnv', 'muse', 
                            'mutect2', 'somaticsniper', 'varscan2']:
        raise ValueError('Unrecognized dataset_type: {}.'.format(dataset_type))

    gdc_rna_counts = {'analysis.workflow_type': 'HTSeq - Counts'}
    gdc_rna_fpkm = {'analysis.workflow_type': 'HTSeq - FPKM'}
    gdc_rna_fpkmuq = {'analysis.workflow_type': 'HTSeq - FPKM-UQ'}
    gdc_mirna = {'data_type': 'miRNA Expression Quantification'}
    gdc_mirna_isoform = {'data_type': 'Isoform Expression Quantification'}
    gdc_cnv = {'data_type': 'Copy Number Segment'}
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

    trans_args_rna = {'read_table_kwargs': {'header': None, 'comment': '_'},
                      'matrix_process': process_average_log}
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
                      'matrix_process': process_maf}

    dataset_map = {'rna_counts': {'gdc_type': gdc_rna_counts,
                                  'trans_args': trans_args_rna,
                                  'matrix_name': '{}.htseq.counts.tsv'},
                   'rna_fpkm': {'gdc_type': gdc_rna_fpkm,
                                'trans_args': trans_args_rna,
                                'matrix_name': '{}.htseq.fpkm.tsv'},
                   'rna_fpkmuq': {'gdc_type': gdc_rna_fpkmuq,
                                  'trans_args': trans_args_rna,
                                  'matrix_name': '{}.htseq.fpkm-uq.tsv'},
                   'mirna': {'gdc_type': gdc_mirna,
                             'trans_args': trans_args_mirna,
                             'matrix_name': '{}.mirna.tsv'},
                   'mirna_isoform': {'gdc_type': gdc_mirna_isoform,
                                     'trans_args': trans_args_mirna_isoform,
                                     'matrix_name': '{}.mirna.isoform.tsv'}, 
                   'cnv': {'gdc_type': gdc_cnv, 
                           'trans_args': trans_args_cnv, 
                           'matrix_name': '{}.cnv.tsv'},
                   'masked_cnv': {'gdc_type': gdc_masked_cnv,
                                  'trans_args': trans_args_cnv,
                                  'matrix_name': '{}.masked.cnv.tsv'},
                   'muse': {'gdc_type': gdc_snv_muse,
                            'trans_args': trans_args_snv,
                            'matrix_name': '{}.muse.snv.tsv'},
                   'mutect2': {'gdc_type': gdc_snv_mutect2,
                               'trans_args': trans_args_snv,
                               'matrix_name': '{}.mutect2.snv.tsv'},
                   'somaticsniper': {
                           'gdc_type': gdc_snv_somaticsniper,
                           'trans_args': trans_args_snv,
                           'matrix_name': '{}.somaticsniper.snv.tsv'
                        },
                   'varscan2': {'gdc_type': gdc_snv_varscan2,
                                'trans_args': trans_args_snv,
                                'matrix_name': '{}.varscan2.snv.tsv'}}

    if os.path.isdir(work_path):
        work_dir = work_path
    else:
        work_dir = os.path.dirname(work_path)

    if mode == 'all' or mode == 'download':
        gdc_type = dataset_map[dataset_type]['gdc_type']
        dataset = download_dataset(projects=project, dataset_type=gdc_type,
                                   work_dir=work_dir)
        file_list = dataset.values()[0]

    if mode == 'transform':
        if os.path.isdir(work_path):
            for f in os.listdir(work_dir):
                file_list.append(os.path.join(work_dir, f))
        else:
            file_list = [work_path]
    
    if mode == 'all' or mode == 'transform':
        trans_args = dataset_map[dataset_type]['trans_args']
        xena_matrix = xena_matrix_merge(file_list, **trans_args)
        if matrix_dir is None:
            matrix_dir = work_dir
            while project in os.path.dirname(matrix_dir):
                matrix_dir = os.path.dirname(matrix_dir)
        matrix_dir = os.path.join(matrix_dir, 'Xena_Matrices')
        mkdir_p(matrix_dir)
        matrix_name = dataset_map[dataset_type]['matrix_name']
        matrix_path = os.path.join(matrix_dir, matrix_name.format(project))
        print('Saving matrix to {} ...'.format(matrix_path))
        xena_matrix.to_csv(matrix_path, sep='\t')
    
    if mode == 'all':
        render_metadata(matrix_path)
    
    print('Data transformation is finished.')

def main():
    print('A python module of Xena specific importing pipeline for GDC data.')

#    l = os.listdir(r'gitignore\TCGA-CHOL\Xena_Matrices')
#    for f in l:
#        render_metadata(os.path.join(r'gitignore\TCGA-CHOL\Xena_Matrices', f))

if __name__ == '__main__':
    main()
