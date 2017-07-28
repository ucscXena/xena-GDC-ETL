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
import warnings

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

def read_by_ext(filename, mode='r'):
    """Automatically decide how to open a file which might be compressed.
    
    Leveraged codes from the "hook_compressed" function in python's fileinput 
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
    
    print('\rAveraging duplicated samples ...', end='')
    df_avg = df.groupby(df.columns, axis=1).mean()
    print('\rLog transforming data ...', end='')
    return np.log2(df_avg + 1)

def process_maf(df):
    """Transform pre-sliced GDC's MAF data into Xena data matrix.
    
    A new column of DNA variant allele frequencies named "dna_vaf" will 
    calculated by division "t_alt_count" / "t_depth". Columns "t_alt_count" 
    and "t_depth" will then be dropped. At last column will be renamed 
    accordingly and row index will be set as sample ID.
    
    Args:
        df: pandas.DataFrame
    
    Returns:
        df_transformed: Transformed pandas.DataFrame
    """
    print('\rCalculating "dna_vaf" ...', end='')
    df['dna_vaf'] = df['t_alt_count'] / df['t_depth']
    print('\rTrim "Tumor_Sample_Barcode" into Xena sample ID ...', end='')
    trim_func = lambda x: '-'.join(x.split('-', 4)[0:4])
    df['Tumor_Sample_Barcode'] = df['Tumor_Sample_Barcode'].apply(trim_func)
    print('\rRe-organizing matrix ...', end='')
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

class XenaDataset(object):
    
    # Map Xena dtype code to GDC data query dict
    __XENA_GDC_DTYPE = {
            'htseq.counts': {'data_type': 'Gene Expression Quantification',
                             'analysis.workflow_type': 'HTSeq - Counts'},
            'htseq.fpkm': {'data_type': 'Gene Expression Quantification',
                           'analysis.workflow_type': 'HTSeq - FPKM'},
            'htseq.fpkm-uq': {'data_type': 'Gene Expression Quantification',
                              'analysis.workflow_type': 'HTSeq - FPKM-UQ'},
            'mirna': {'data_type': 'miRNA Expression Quantification',
                      'analysis.workflow_type': 'BCGSC miRNA Profiling'},
            'mirna.isoform': {
                    'data_type': 'Isoform Expression Quantification',
                    'analysis.workflow_type': 'BCGSC miRNA Profiling'
                },
            'cnv': {'data_type': 'Copy Number Segment',
                    'analysis.workflow_type': 'DNAcopy'},
            'masked.cnv': {'data_type': 'Masked Copy Number Segment',
                           'analysis.workflow_type': 'DNAcopy'},
            'muse.snv': {
                    'data_type': 'Masked Somatic Mutation',
                    'analysis.workflow_type': 
                        'MuSE Variant Aggregation and Masking'
                },
            'mutect2.snv': {
                    'data_type': 'Masked Somatic Mutation',
                    'analysis.workflow_type': 
                        'MuTect2 Variant Aggregation and Masking'
                },
            'somaticsniper.snv': {
                    'data_type': 'Masked Somatic Mutation',
                    'analysis.workflow_type': 
                        'SomaticSniper Variant Aggregation and Masking'
                },
            'varscan2.snv': {
                    'data_type': 'Masked Somatic Mutation',
                    'analysis.workflow_type': 
                        'VarScan2 Variant Aggregation and Masking'
                }
        }

    # Settings for making Xena matrix from raw data
    __RNA_TRANSFORM_ARGS = {
            'read_table_kwargs': {'header': None, 'comment': '_'},
            'merge_axis': 1,
            'matrix_process': process_average_log,
            'index_name': 'Ensembl_ID'
        }
    __MIRNA_TRANSFORM_ARGS = {
            'read_table_kwargs': {'header': 0, 'usecols': [0, 2]},
            'merge_axis': 1,
            'matrix_process': process_average_log
        }
    __MIRNA_ISOFORM_TRANSFORM_ARGS = {
            'read_table_kwargs': {'header': 0, 'usecols': [1, 3]},
            'merge_axis': 1,
            'matrix_process': process_average_log
        }
    __CNV_TRANSFORM_ARGS = {
            'read_table_kwargs': {'header': 0, 'usecols': [1, 2, 3, 5]},
            'merge_axis': 0
        }
    __SNV_TRANSFORM_ARGS = {
            'read_table_kwargs': {'header': 0, 
                                  'usecols': [12, 36, 4, 5, 6, 39, 41, 51, 0, 
                                              10, 15, 110],
                                  'comment': '#'},
            'merge_axis': 0,
            'matrix_process': process_maf,
            'index_name': 'Sample_ID'
        }
    __TRANSFORM_ARGS = {
            'htseq.counts': __RNA_TRANSFORM_ARGS,
            'htseq.fpkm': __RNA_TRANSFORM_ARGS,
            'htseq.fpkm-uq': __RNA_TRANSFORM_ARGS,
            'mirna': __MIRNA_TRANSFORM_ARGS,
            'mirna.isoform': __MIRNA_ISOFORM_TRANSFORM_ARGS,
            'cnv': __CNV_TRANSFORM_ARGS,
            'masked.cnv': __CNV_TRANSFORM_ARGS,
            'muse.snv': __SNV_TRANSFORM_ARGS,
            'mutect2.snv': __SNV_TRANSFORM_ARGS,
            'somaticsniper.snv': __SNV_TRANSFORM_ARGS,
            'varscan2.snv': __SNV_TRANSFORM_ARGS
        }

    # Map xena_dtype to corresponding metadata template.
    __METADATA_TEMPLATE_DIR = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), 'Resources'
        )
    __METADATA_TEMPLATE = {'htseq.counts': 'template.rna.meta.json',
                           'htseq.fpkm': 'template.rna.meta.json',
                           'htseq.fpkm-uq': 'template.rna.meta.json',
                           'mirna': 'template.mirna.meta.json',
                           'mirna.isoform': 'template.mirna.isoform.meta.json',
                           'cnv': 'template.cnv.meta.json',
                           'masked.cnv': 'template.cnv.meta.json',
                           'muse.snv': 'template.snv.meta.json',
                           'mutect2.snv': 'template.snv.meta.json',
                           'somaticsniper.snv': 'template.snv.meta.json',
                           'varscan2.snv': 'template.snv.meta.json'}
    # Map GDC project_id to Xena specific cohort name.
    __XENA_COHORT = {'TCGA-BRCA': 'TCGA Breast Cancer (BRCA)',
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
    # Map dataset_map to (unique) GDC data type or GDC workflow type.
    __METADATA_GDC_TYPE = {
            'htseq.counts': 'HTSeq - Counts', 
            'htseq.fpkm': 'HTSeq - FPKM', 
            'htseq.fpkm-uq': 'HTSeq - FPKM-UQ', 
            'mirna': 'miRNA Expression Quantification', 
            'mirna.isoform': 'Isoform Expression Quantification', 
            'cnv': 'Copy Number Segment', 
            'masked.cnv': 'Masked Copy Number Segment', 
            'muse.snv': 'MuSE Variant Aggregation and Masking', 
            'mutect2.snv': 'MuTect2 Variant Aggregation and Masking', 
            'somaticsniper.snv': 
                'SomaticSniper Variant Aggregation and Masking', 
            'varscan2.snv': 'VarScan2 Variant Aggregation and Masking'
        }
    
    def __init__(self, projects, xena_dtype, root_dir='.', 
                 raw_data_dir=None, matrix_dir=None):
        self.projects = projects
        self.xena_dtype = xena_dtype
        self.root_dir = root_dir
        if raw_data_dir is not None:
            self.raw_data_dir = raw_data_dir
        if matrix_dir is not None:
            self.matrix_dir = matrix_dir
    
    # Define projects property
    def __get_projects(self):
        return self.__projects
    
    def __set_projects(self, projects):
        if isinstance(projects, str):
            self.__projects = [projects]
        elif isinstance(projects, list):
            self.__projects = projects
        else:
            raise ValueError('"projects" must be either str or list.')
    
    projects = property(__get_projects, __set_projects,
                        doc="""A list of GDC's project_id(s). Corresponding
                            projects will be included in this dataset.""")
    
    # Define xena_dtype property
    def __get_xena_dtype(self):
        return self.__xena_dtype
    
    def __set_xena_dtype(self, xena_dtype):
        if xena_dtype in self.__XENA_GDC_DTYPE:
            self.__xena_dtype = xena_dtype
            self.gdc_dtype = self.__XENA_GDC_DTYPE[xena_dtype]
            self.transform_args = self.__TRANSFORM_ARGS[xena_dtype]
            self.metadata_template = self.__METADATA_TEMPLATE[xena_dtype]
        else:
            raise ValueError("Unsupported data type: {}".format(xena_dtype))
    
    xena_dtype = property(__get_xena_dtype, __set_xena_dtype, 
                          doc="""A string of one dataset type supported by 
                              this class. To get a list of supported types, 
                              use "get_supported_dtype()".""")
    
    def get_supported_dtype(self):
        return self.__XENA_GDC_DTYPE.keys()
    
    # Setup default directory structure from "root_dir"
    def __get_root_dir(self):
        """The "root_dir" property defines the root directory for this 
        dataset.
        
        By default, all files related to this class, such as raw data, 
        Xena matrix, metadata, should and highly recommended to be organized 
        and saved under this directory. The default directory structure is:
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

        When instantiating this class, default "root_dir" points to current 
        python work directory. Setting "root_dir" will not only set the 
        "root_dir" property but also set "raw_data_dir" and/or "matrix_dir" 
        properties if the directory for GDC data ("raw_data_dir") and/or Xena 
        matrix ("matrix_dir") is not set yet. However, no actual directories 
        will be made by just setting these properties. Directories will only 
        be made when needed.
        
        If you want to reset the default directory structure while 
        "raw_data_dir" and/or "matrix_dir" properties are set already, you can 
        override it using the "set_default_dir_tree" method with 
        "reset_default=True".
        """
        
        return self.__root_dir
    
    def set_default_dir_tree(self, root_dir, reset_default=False):
        if os.path.isdir(root_dir):
            self.__root_dir = os.path.abspath(root_dir)
            self.__dataset_dir = os.path.join(self.__root_dir,
                                              '_'.join(self.projects))
            if (not hasattr(self, 'raw_data_dir')) or reset_default:
                self.raw_data_dir = os.path.join(
                        self.__dataset_dir,
                        'GDC_Raw_Data',
                        self.xena_dtype.replace('.', '_')
                    )
            if (not hasattr(self, 'matrix_dir')) or reset_default:
                self.matrix_dir = os.path.join(self.__dataset_dir,
                                               'Xena_Matrices')
        else:
            raise OSError('{} is not a valid directory.'.format(root_dir))
    
    root_dir = property(__get_root_dir, set_default_dir_tree)

    def download(self, label_field='cases.samples.submitter_id'):
        """Download GDC's open access data for this dataset.
        
        Data was selected based on "projects" and "gdc_dtype" properties. 
        "gdc_dtype" can be assigned directly or set indirectly through 
        "xena_dtype". All open access data matches the two criteria 
        will be put together under one directory as a single dataset. Though 
        this method won't check, putting different types of data into one 
        dataset is NOT recommended. 
        
        By default, data files will be renamed to 
        "<cases.samples.submitter_id>.<UUID>.<file extension>" and saved under 
        the directory defined by the "raw_data_dir" property. Check the 
        "root_dir" property for details about the default directory structure. 
        
        A list of paths for downloaded files will be assigned to the 
        "raw_data_list" property which can be used for Xena matrix "transform" 
        processing. Check the "transform" method for details.
        
        Args:
            label_field: str, default 'cases.samples.submitter_id'
                A single GDC available file field whose value will be used for 
                renaming downloaded files. By default, Xena uses 
                'cases.samples.submitter_id' as sample ID. Therefore filenames 
                for downloaded files will be 
                "submitter_id.UUID.file_extension".
        
        Returns:
            self: allow method chaining.
        """
        
        assert self.projects is not None
        assert hasattr(self, 'gdc_dtype') and isinstance(self.gdc_dtype, dict)
        query_dict = {'access': 'open', 
                      'cases.project.project_id': self.projects}
        query_dict.update(self.gdc_dtype)
        print('Searching for raw data ...', end='')
        file_dict = gdc.get_file_dict(query_dict, label_field=label_field)
        if not file_dict:
            message = '\rNo {} data found for project {}.'
            dataset_description = ' - '.join(sorted(self.gdc_dtype.values()))
            print(message.format(dataset_description, str(self.projects)))
            return self
        print('\r{} files found for {} data of {}.'.format(len(file_dict),
                                                           self.xena_dtype,
                                                           self.projects))
        assert hasattr(self, 'raw_data_dir')
        self.raw_data_dir = os.path.abspath(self.raw_data_dir)
        mkdir_p(self.raw_data_dir)
        print('Datasets will be downloaded to')
        print(self.raw_data_dir)
        file_dict = {
                k: os.path.join(self.raw_data_dir, v) 
                for k, v in file_dict.items()
            }
        self.raw_data_list = gdc.download(file_dict)
        print('Raw {} data for {} is ready.'.format(self.projects, 
                                                    self.xena_dtype))
        return self

    def transform(self, raw_data_list=None, matrix_name=None):
        """Transform raw data in a dataset into Xena matrix
        
        Args:
            raw_data_list: str or list, default None
                One (str) or a list of raw data file(s) used for building Xena 
                matrix. This "raw_data_list" argument has a higher priority 
                than the "raw_data_list" property, i.e. if "raw_data_list" is 
                not None, it will be used to define or overwrite the 
                "raw_data_list" property. If this "raw_data_list" argument is 
                None, the "raw_data_list" property will be checked first. If 
                "raw_data_list" property is defined and is not None, it will 
                be used. The "raw_data_list" property can be assigned directly 
                or be set by the "download" method. Check the "download" 
                method for details. If the "raw_data_list" property is not 
                usable, the "raw_data_dir" property will be checked. All files 
                under this directory will be treated as data and used for 
                building Xena matrix.
            matrix_name: str, default None
                File name of the Xena matrix. The matrix will be save under 
                the directory defined by the "matrix_dir" property. If None, 
                default filename will be used, which has a pattern as 
                "projects.xena_type.tsv". Full path of this matrix will be 
                assigned to the "matrix" property which can be used for making 
                metadata. Check the "metadata" method for details.
        
        Returns:
            self: allow method chaining.
        """
        
        if raw_data_list is not None:
            if isinstance(raw_data_list, str):
                self.raw_data_list = [raw_data_list]
            elif isinstance(raw_data_list, list):
                self.raw_data_list = raw_data_list
            else:
                message = '"raw_data_list" should be str or list, not {}.'
                raise TypeError(message.format(type(raw_data_list)))
        if ((not hasattr(self, 'raw_data_list')) 
            or (self.raw_data_list is None)):
            try:
                raw_data_dir = os.path.abspath(self.raw_data_dir)
                raw_data_list = []
                for f in os.listdir(raw_data_dir):
                    f_path = os.path.join(raw_data_dir, f)
                    if os.path.isfile(f_path):
                        raw_data_list.append(f_path)
                if not raw_data_list:
                    raise ValueError
                self.raw_data_list = raw_data_list
            except:
                raise ValueError('Cannot find raw data for Xena matrix '
                                 'transformation.')

        # Start transformation
        merge_axis = self.transform_args['merge_axis']
        total = len(self.raw_data_list)
        count = 0
        df_list = []
        for path in self.raw_data_list:
            count = count + 1
            print('\rProcessing {}/{} file...'.format(count, total), end='')
            sys.stdout.flush()
            with read_by_ext(path) as f:
                sample_id = os.path.basename(path).split('.', 1)[0]
                if merge_axis == 1:
                    pd_kwargs = {'index_col': 0, 'usecols': [0, 1]}
                elif merge_axis == 0:
                    pd_kwargs = {}
                else:
                    message = 'Invalid merge_axis: {}'.format(merge_axis)
                    raise ValueError(message)
                pd_kwargs.update(self.transform_args['read_table_kwargs'])
                df = pd.read_table(f, **pd_kwargs)
                if merge_axis == 1:
                    df.columns = [sample_id]
                elif merge_axis == 0:
                    df.set_index([[sample_id] * df.shape[0]], inplace=True)
                df_list.append(df)
        print('\rAll {} files have been processed. '.format(total))
        print('Merging into one matrix ...', end='')
        xena_matrix = pd.concat(df_list, axis=merge_axis)
        if 'matrix_process' in self.transform_args:
            xena_matrix = self.transform_args['matrix_process'](xena_matrix)
        if 'index_name' in self.transform_args:
            xena_matrix.index.name = self.transform_args['index_name']
        # Transformation done
        if matrix_name is None:
            matrix_name = '{}.{}.tsv'.format('_'.join(self.projects),
                                             self.xena_dtype)
        self.matrix = os.path.join(os.path.abspath(self.matrix_dir), 
                                   matrix_name)
        print('\rSaving matrix to {} ...'.format(self.matrix), end='')
        mkdir_p(os.path.abspath(self.matrix_dir))
        xena_matrix.to_csv(self.matrix, sep='\t')
        print('\rXena matrix is saved at {}.'.format(self.matrix))
        return self
    
    def metadata(self):
        """Make "metadata.json" for Xena data loading
        
        One metadata will be created for one Xena matrix defined by the 
        "matrix" property. The metadata JSON file will be saved under the same 
        directory as the matrix file and named with a ".json" postfix appended 
        to the filename of Xena matrix. The "transform" method will generate a 
        Xena matrix from raw data and assign it to "matrix"; or "matrix" can 
        be assigned directly. JSON templates for metatdata are defined by 
        "__METADATA_TEMPLATE_DIR" and "__METADATA_TEMPLATE" constants; 
        pre-defined Xena cohort names are defined by the "__XENA_COHORT" 
        constant; data type labels in the metadata are defined by the 
        "__METADATA_GDC_TYPE" constant. Specific raw data (MAF) urls for 
        "Masked Somatic Mutation" data will be queried according to "projects" 
        and "gdc_dtype" properties.
        
        Returns:
            self: allow method chaining.
        """
        
        if ((not hasattr(self, 'matrix')) or (self.matrix is None) 
            or (not os.path.isfile(self.matrix))):
            raise ValueError('Cannot find Xena matrix for this dataset; '
                             'please create a matrix or assign a matrix with '
                             'the "matrix" property before making metadata.')
        else:
            self.matrix = os.path.abspath(self.matrix)

        # Start to generate metadata.
        # General jinja2 Expressions
        print('Creating metadata file ...', end='')
        matrix_date = time.strftime("%m-%d-%Y", 
                                    time.gmtime(os.path.getmtime(self.matrix)))
        projects = ','.join(self.projects)
        variables = {'project_id': projects, 
                     'gdc_type': self.__METADATA_GDC_TYPE[self.xena_dtype], 
                     'date': matrix_date}
        if projects in self.__XENA_COHORT:
            variables['xena_cohort'] = self.__XENA_COHORT[projects]
        else:
            variables['xena_cohort'] = projects
        # Data type specific jinja2 Expression
        if self.xena_dtype == 'htseq.fpkm':
            variables['unit'] = 'fpkm'
        if self.xena_dtype == 'htseq.fpkm-uq':
            variables['unit'] = 'fpkm-uq'
        if self.xena_dtype in ['muse.snv', 'mutect2.snv', 'somaticsniper.snv', 
                               'varscan2.snv']:
            query_dict = {'access': 'open',
                          'cases.project.project_id': self.projects}
            query_dict.update(self.gdc_dtype)
            try:
                print('\rSearching the specific URL for raw MAF data ...', 
                      end='')
                res_df = gdc.search('files', query_dict, 'file_id')
                if res_df['file_id'].shape == (1,):
                    variables['maf_uuid'] = str(res_df['file_id'][0])
            except:
                message = ('Fail to get a specific URL for the MAF file for: ' 
                           'matrix "{}"; "{}" data of cohort "{}".')
                warnings.warn(message.format(self.matrix, 
                                             variables['project_id'], 
                                             variables['xena_cohort']), 
                              stacklevel=2)
        
        # Render jinja2 template
        jinja2_env = jinja2.Environment(
                loader=jinja2.FileSystemLoader(self.__METADATA_TEMPLATE_DIR))
        template_json = self.__METADATA_TEMPLATE[self.xena_dtype]
        template = jinja2_env.get_template(template_json)
        self.metadata = os.path.join(self.matrix + '.json')
        with open(self.metadata, 'w') as f:
            f.write(template.render(**variables))
        print('\rMetadata JSON is saved at {}.'.format(self.metadata))
        return self
    
def main():
    print('A python module of Xena specific importing pipeline for GDC data.')
    
#    test = XenaDataset('TCGA-CHOL', 'htseq.fpkm', r'E:\test')
#    test.download().transform().metadata()
#    test.transform().metadata()
    
if __name__ == '__main__':
    main()
