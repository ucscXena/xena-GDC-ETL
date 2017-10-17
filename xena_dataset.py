#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This module mainly provides a ``XenaDataset`` class representing one Xena 
matrix in a Xena cohort. Three class, ``GDCOmicset``, ``TCGAPhenoset`` and 
``TARGETPhenoset`` are derived from ``XenaDataset``, representing genomic 
data, phenotype info of TCGA and phenotype info of TARGET respectively.

In general, a ``XenaDataset`` class contains 3 methods, ``download``, 
``transform`` and ``metadata``, which can be used for quickly assembling an 
ETL pipeline importing data into Xena.
"""

# Ensure Python 2 and 3 compatibility
from __future__ import division
from __future__ import print_function

import time
import os
import sys
import warnings

import jinja2
from lxml import etree
import numpy as np
import pandas as pd
import requests

import gdc

def mkdir_p(dir_name):
    """Make the directory as needed: no error if existing.
    
    Args:
        dir_name (str): Directory name or path.
    
    Returns:
        str: The absolute path for the directory.
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
        filename (str): Must contain proper extension indicating the 
            compression condition.
        mode (str, optional): To specify the mode in which the file is opened. 
            It will be passed to corresponding open function (``open``, 
            ``gzip.open`` or ``bz2.BZ2File``); please check them for details. 
            Defaults to 'r'.
    
    Returns:
        file object: A filehandle to be used with `with`.
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


def read_biospecimen(fileobj):
    """Extract info from GDC's biospecimen supplement and re-organize them 
    into a pandas DataFrame.
    
    Args:
        fileobj (file or path): XML file of GDC's biospecimen supplement.
    
    Returns:
        pandas.core.frame.DataFrame: Transformed pandas DataFrame.
    """
    
    if isinstance(fileobj, file):
        filename = fileobj.name
    else:
        filename = os.path.basename(fileobj)
    ext = os.path.splitext(filename)[1]
    if ext == '.xlsx':
        # Design specifically for TARGET biospecimen
        df = pd.read_excel(filename, header=None)
        df.iloc[0].fillna(method='ffill', inplace=True)
        df.columns = df.iloc[0:2].apply(lambda x: x.str.cat(sep='.'))
        return df.drop(df.index[0:2]).set_index(df.columns[0])
    elif ext != '.xml':
        raise IOError('Unknown file type for biospecimen data: {}'.format(ext))
    
    disease_dict = {
            'LAML': 'Acute Myeloid Leukemia',
            'ACC': 'Adrenocortical carcinoma',
            'BLCA': 'Bladder Urothelial Carcinoma',
            'LGG': 'Brain Lower Grade Glioma',
            'BRCA': 'Breast invasive carcinoma',
            'CESC': ('Cervical squamous cell carcinoma and endocervical '
                     'adenocarcinoma'),
            'CHOL': 'Cholangiocarcinoma',
            'LCML': 'Chronic Myelogenous Leukemia',
            'COAD': 'Colon adenocarcinoma',
            'CNTL': 'Controls',
            'ESCA': 'Esophageal carcinoma',
            'FPPP': 'FFPE Pilot Phase II',
            'GBM': 'Glioblastoma multiforme',
            'HNSC': 'Head and Neck squamous cell carcinoma',
            'KICH': 'Kidney Chromophobe',
            'KIRC': 'Kidney renal clear cell carcinoma',
            'KIRP': 'Kidney renal papillary cell carcinoma',
            'LIHC': 'Liver hepatocellular carcinoma',
            'LUAD': 'Lung adenocarcinoma',
            'LUSC': 'Lung squamous cell carcinoma',
            'DLBC': 'Lymphoid Neoplasm Diffuse Large B-cell Lymphoma',
            'MESO': 'Mesothelioma',
            'MISC': 'Miscellaneous',
            'OV': 'Ovarian serous cystadenocarcinoma',
            'PAAD': 'Pancreatic adenocarcinoma',
            'PCPG': 'Pheochromocytoma and Paraganglioma',
            'PRAD': 'Prostate adenocarcinoma',
            'READ': 'Rectum adenocarcinoma',
            'SARC': 'Sarcoma',
            'SKCM': 'Skin Cutaneous Melanoma',
            'STAD': 'Stomach adenocarcinoma',
            'TGCT': 'Testicular Germ Cell Tumors',
            'THYM': 'Thymoma',
            'THCA': 'Thyroid carcinoma',
            'UCS': 'Uterine Carcinosarcoma',
            'UCEC': 'Uterine Corpus Endometrial Carcinoma',
            'UVM': 'Uveal Melanoma',
        }
    root = etree.parse(fileobj).getroot()
    ns = root.nsmap
    samples_common = {}
    for child in root.find('admin:admin', ns):
        try:
            samples_common[child.tag.split('}', 1)[-1]] = child.text.strip()
        except AttributeError:
            samples_common[child.tag.split('}', 1)[-1]] = ''
    for child in root.find('bio:patient', ns):
        try:
            samples_common[child.tag.split('}', 1)[-1]] = child.text.strip()
        except AttributeError:
            samples_common[child.tag.split('}', 1)[-1]] = ''
    # Add 'primary_diagnosis' according to
    # https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations
    samples_common['primary_diagnosis'] = disease_dict[
            samples_common['disease_code']
        ]
    
    samples = {}
    for sample in root.find('bio:patient/bio:samples', ns):
        record = {}
        for child in sample:
            if child.text and child.text.strip():
                record[child.tag.split('}', 1)[-1]] = child.text.strip()
        record.update(samples_common)
        samples[record['bcr_sample_barcode']] = record
    df = pd.DataFrame(samples).T
    sample_mask = df.bcr_sample_barcode.map(
            lambda s: s[-3:-1] not in ['10']
        )
    df = df[sample_mask]
    df['bcr_patient_barcode'] = root.find(
            'bio:patient/shared:bcr_patient_barcode', ns
        ).text
    return df


def read_clinical(fileobj):
    """Extract info from GDC's clinical supplement and re-organize them into a 
    pandas DataFrame.
    
    Args:
        fileobj (file or path): XML file of GDC's clinical supplement.
    
    Returns:
        pandas.core.frame.DataFrame: Transformed pandas DataFrame.
    """
    
    if isinstance(fileobj, file):
        filename = fileobj.name
    else:
        filename = os.path.basename(fileobj)
    ext = os.path.splitext(filename)[1]
    if ext == '.xlsx':
        return pd.read_excel(filename, sheetname='Clinical Data', index_col=0)
    elif ext != '.xml':
        raise IOError('Unknown file type for clinical data: {}'.format(ext))
    
    root = etree.parse(fileobj).getroot()
    ns = root.nsmap
    patient = {}
    # "Dirty" extraction
    for child in root.xpath('.//*[not(*)]'):
        try:
            patient[child.tag.split('}', 1)[-1]] = child.text.strip()
        except AttributeError:
            patient[child.tag.split('}', 1)[-1]] = ''
    # Redo 'race'
    if 'race_list' in patient:
        del patient['race_list']
    try:
        patient['race'] = ','.join(
                [child.text.strip()
                 for child in root.find('.//clin_shared:race_list', ns) 
                 if child.text and child.text.strip()]
            )
    except:
        patient['race'] = ''
    # Redo the most recent "follow_up" and update the patient dict if there is 
    # an overlapped key.
    follow_ups = root.xpath('.//*[local-name()="follow_up"]')
    if follow_ups:
        most_recent = follow_ups[0]
        for follow_up in follow_ups:
            if follow_up.attrib['version'] > most_recent.attrib['version']:
                most_recent = follow_up
        for child in most_recent:
            try:
                patient[child.tag.split('}', 1)[-1]] = child.text.strip()
            except AttributeError:
                patient[child.tag.split('}', 1)[-1]] = ''
    return pd.DataFrame({patient['bcr_patient_barcode']: patient}).T


def rna_columns_matrix(df_list):
    """Merge and process a list of dataframes to make a Xena data matrix.
    
    Every dataframe contains data for individual sample. They are merged 
    horizontally (``axis=1``; growing with more and more columns). For merged 
    matrix, first average columns having the same name and then transform its 
    data by log(x + 1).
    
    Args:
        df_list (list of pandas.core.frame.DataFrame): Input raw data. One 
            dataframe is expected to have data for one sample, with column 
            header being sample ID (not strictly required).
    
    Returns:
        pandas.core.frame.DataFrame: Ready to load Xena matrix.
    """
    
    print('\rMerging into one matrix ...', end='')
    df = pd.concat(df_list, axis=1)
    # ``pandas.concat`` Will lose index name when there are rows outside 
    # intersection during outter join.
    if df.index.name is None:
        df.index.name = df_list[0].index.name
    print('\rAveraging duplicated samples ...', end='')
    df_avg = df.groupby(df.columns, axis=1).mean()
    print('\rLog transforming data ...       ', end='')
    return np.log2(df_avg + 1)


def snv_maf_matrix(df_list):
    """Transform pre-sliced GDC's MAF data into Xena data matrix.
    
    A new column of DNA variant allele frequencies named "dna_vaf" will 
    calculated by division "t_alt_count" / "t_depth". Columns "t_alt_count" 
    and "t_depth" will then be dropped. In the end, columns will be renamed 
    accordingly and row index will be set as sample ID.
    
    Args:
        df (pandas.core.frame.DataFrame): Input raw matrix for mutation data 
           (from MAF file).
    
    Returns:
        pandas.core.frame.DataFrame: Transformed pandas DataFrame.
    """
    
    assert len(df_list) == 1
    df = df_list[0]
    print('\rCalculating "dna_vaf" ...', end='')
    df['dna_vaf'] = df['t_alt_count'] / df['t_depth']
    print('\rTrim "Tumor_Sample_Barcode" into Xena sample ID ...', end='')
    trim_func = lambda x: '-'.join(x.split('-', 4)[0:4])
    df['Tumor_Sample_Barcode'] = df['Tumor_Sample_Barcode'].apply(trim_func)
    print('\rRe-organizing matrix ...', end='')
    df = (df.drop(['t_alt_count', 't_depth'], axis=1)
            .set_index('Tumor_Sample_Barcode')
            .rename(columns={'Hugo_Symbol': 'gene', 
                             'Chromosome': 'chrom', 
                             'Start_Position': 'start', 
                             'End_Position': 'end', 
                             'Reference_Allele': 'ref', 
                             'Tumor_Seq_Allele2': 'alt', 
                             'Tumor_Sample_Barcode': 'sampleid', 
                             'HGVSp_Short': 'Amino_Acid_Change', 
                             'Consequence': 'effect',
                             'FILTER': 'filter'}))
    df.index.name = 'Sample_ID'
    return df


class XenaDataset(object):
    """XenaDataset represents for one Xena matrix in a Xena cohort.

    This class provides a set of method for downloading and transforming data 
    into Xena matrix, as well as generating associated metadata. It also 
    provides a set of attributes to control these processes.

    Attributes:
        projects (str or list): One (string) or a list of project IDs 
            describing study projects included in this dataset.
        xena_dtype (str): A short string (liek an ID) describing the type of 
            data in this dataset. It is highly recommended, though not 
            required, that one dataset has a single type of data.
        root_dir (str): Defines the root directory for this dataset.
            By default, all files related to this class, such as raw data, 
            Xena matrix, metadata, should and highly recommended to be 
            organized and saved under this directory. However, this is not 
            required. You can set directories for raw data (through the 
            ``raw_data_dir`` property) and Xena matrix (through the 
            ``matrix_dir``property) specifically. No actual directories will 
            be made by just setting these properties. Directories will only be 
            made when needed.The default directory structure is::
            
                root_dir
                └── projects
                    ├── "Raw_Data"
                    │   └── xena_dtype
                    │       ├── data1
                    │       ├── data2
                    │       ├── ...
                    │       └── dataN
                    └── "Xena_Matrices"
                        ├── projects.xena_dtype(1).tsv
                        ├── projects.xena_dtype(1).tsv.json
                        ├── projects.xena_dtype(2).tsv
                        ├── projects.xena_dtype(2).tsv.json
                        ├── ...
                        ├── projects.xena_dtype(N).tsv
                        └── projects.xena_dtype(N).tsv.json
            
            Defaults to "." which points to current python work directory.
        
        raw_data_dir (str): A directory for saving raw data. Check the 
            ``root_dir`` property for the default directory structure. 
            Defaults to None.
        download_map (dict): A dict with the key being a URL for one raw data 
            to be downloaded and the value being a path for saving downloaded 
            raw data.
        raw_data_list (list): A list of file path(s) for all raw data related 
            to this dataset. It will be automatically set by the ``download`` 
            method; or it can be assigned directly as a public attribute. This 
            ``raw_data_list`` attribute will be used by the ``transform`` 
            method for making a Xena matrix from raw data. If the 
            ``raw_data_list`` is not available when trying to get this 
            attribute, the ``raw_data_dir`` property will be checked. All 
            files under ``raw_data_dir`` will be treated as data and used for 
            creating a ``raw_data_list``.
        read_raw (func): A function used for reading a raw data. A file object 
            will be passed to it.
        raws2matrix (func): A function used for merging multiple raw data read 
            by ``read_raw`` into one Xena matrix. A list of data read by 
            ``read_raw`` will be passed to this function.
        matrix_dir (str): A path for saving the Xena matrix for this 
            dataset. If the ``matrix_dir`` is not available when trying to get 
            this attribute, the ``matrix`` property will be checked first. If 
            the ``matrix`` property is available, its directory will be 
            assigned to the ``matrix_dir``. If the ``matrix`` property is not 
            available, a default path will be assigned according to the 
            default directory structure. Check the ``root_dir`` property for 
            the default directory structure. Defaults to None.
        matrix (str): A path for the Xena matrix of this dataset. 
            This attribute will be used but not validated by the ``transform`` 
            method for saving newly generated Xena matrix. This attribute will 
            also be used yet will be validated (i.e. it has to point to a 
            valid existing file) by the ``metadata`` method for making 
            metadata associated with the Xena matrix and with this dataset. By 
            default, when used by the ``transform`` method, this ``matrix`` 
            attribute will adapte a pattern as "projects.xena_type.tsv".
        metadata_template (jinja2.environment.Template or str): A Jinja2 
            template for rendering metadata of this dataset. When using a 
            string to set ``metadata_template``, it will be used as a path to 
            the template file and the corresponding template will be got and 
            assigned to this attribute.
        metadata_vars: A dict of variables which will be used (by \*\* 
            unpacking) for rendering the ``metadata_template``.
    """
    
    @property
    def projects(self):
        return self.__projects
    
    @projects.setter
    def projects(self, projects):
        if isinstance(projects, str):
            self.__projects = [projects]
        elif isinstance(projects, list):
            self.__projects = projects
        else:
            raise ValueError('"projects" must be either str or list.')
    
    @property
    def root_dir(self):
        """A path of an existing directory for keeping files (raw data, matrix 
        and metadata) of this dataset.
        """
        
        return self.__root_dir
    
    @root_dir.setter
    def root_dir(self, path):
        if os.path.isdir(path):
            self.__root_dir = os.path.abspath(path)
        else:
            raise IOError('{} is not an existing directory.'.format(path))
    
    @property
    def raw_data_dir(self):
        try:
            return self.__raw_data_dir
        except AttributeError:
            self.__raw_data_dir = os.path.join(
                    self.root_dir, '_'.join(self.projects), 'Raw_Data', 
                    str(self.xena_dtype)
                )
            return self.__raw_data_dir
    
    @raw_data_dir.setter
    def raw_data_dir(self, path):
        self.__raw_data_dir = os.path.abspath(path)
    
    @property
    def download_map(self):
        assert self.__download_map and isinstance(self.__download_map, dict)
        return self.__download_map
    
    @download_map.setter
    def download_map(self, d):
        if isinstance(d, dict):
            self.__download_map = d
        else:
            raise TypeError('download_map should be a dict, '
                            'not {}'.format(type(d)))
    
    # Raw data list: try to get the list from ``raw_data_dir`` if unavailable.
    @property
    def raw_data_list(self):
        try:
            return self.__raw_data_list
        except AttributeError:
            try:
                raw_data_dir = os.path.abspath(self.raw_data_dir)
                raw_data = []
                for f in os.listdir(raw_data_dir):
                    f_path = os.path.join(raw_data_dir, f)
                    if os.path.isfile(f_path):
                        raw_data.append(f_path)
                if raw_data:
                    self.__raw_data_list = raw_data
                else:
                    raise ValueError
            except Exception:
                raise ValueError('Cannot find raw data.')
            return self.__raw_data_list
    
    @raw_data_list.setter
    def raw_data_list(self, raw_data):
        self.__raw_data_list = raw_data
    
    @property
    def matrix_dir(self):
        try:
            return self.__matrix_dir
        except AttributeError:
            try:
                self.__matrix_dir = os.path.dirname(self.__matrix)
                return self.__matrix_dir
            except AttributeError:
                self.__matrix_dir = os.path.join(
                        self.root_dir, '_'.join(self.projects), 
                        'Xena_Matrices'
                    )
                return self.__matrix_dir
    
    @matrix_dir.setter
    def matrix_dir(self, path):
        self.__matrix_dir = os.path.abspath(path)
    
    @property
    def matrix(self):
        try:
            assert self.__matrix
            return self.__matrix
        except (AttributeError, AssertionError):
            self.__matrix = os.path.join(
                    self.matrix_dir,
                    '{}.{}.tsv'.format('_'.join(self.projects), 
                                       self.xena_dtype)
                )
            return self.__matrix
    
    @matrix.setter
    def matrix(self, path):
        self.__matrix = os.path.abspath(path)
        self.matrix_dir = os.path.dirname(self.__matrix)
    
    @property
    def metadata_template(self):
        assert isinstance(self.__metadata_template, 
                          jinja2.environment.Template)
        return self.__metadata_template
    
    @metadata_template.setter
    def metadata_template(self, template):
        if isinstance(template, jinja2.environment.Template):
            self.__metadata_template = template
        elif isinstance(template, str):
            assert os.path.isfile(template)
            file_dir = os.path.dirname(template)
            file_name = os.path.basename(template)
            jinja2_env = jinja2.Environment(
                    loader=jinja2.FileSystemLoader(file_dir)
                )
            self.__metadata_template = jinja2_env.get_template(file_name)
        else:
            raise TypeError('metadata_template should be a jinja2 template '
                            'or an existing path to a template JSON file, '
                            'not a :{}.'.format(type(template)))
    
    def __init__(self, projects, xena_dtype, root_dir='.', 
                 raw_data_dir=None, matrix_dir=None):
        self.projects = projects
        self.xena_dtype = xena_dtype
        self.root_dir = root_dir
        if raw_data_dir is not None:
            self.raw_data_dir = raw_data_dir
        if matrix_dir is not None:
            self.matrix_dir = matrix_dir
    
    def download(self, chunk_size=4096):
        """Download file(s) according to the ``download_map`` property.
        
        A list of paths for downloaded files will be assigned to the 
        ``raw_data_list`` property which can be used for Xena matrix 
        ``transform`` processing. Check the ``transform`` method for details.
        
        Args:
            chunk_size (int, optional): The chunk size is the number of bytes 
                it should read into memory, when the response is got with 
                "stream=True". Check the documentation of "requests" module 
                for details. Defaults to 4096.
        
        Returns:
            self: allow method chaining.
        """
        
        print('Starts to download...', end='')
        total = len(self.download_map)
        count = 0
        download_list = []
        for url, path in self.download_map.items():
            count += 1
            response = requests.get(url, stream=True)        
            if response.status_code == 200:
                file_size = int(response.headers['Content-Length'])
                path = os.path.abspath(path)
                status = '\r[{:d}/{:d}] Download to "{}": {:4.0%}'
                mkdir_p(os.path.dirname(path))
                with open(path, 'wb') as f:
                    downloaded = 0
                    print(status.format(count, total, path, 0), end='')
                    sys.stdout.flush()
                    for chunk in response.iter_content(chunk_size):
                        f.write(chunk)
                        downloaded = downloaded + chunk_size
                        print(status.format(count, total, path, 
                                            min(1, downloaded/file_size)), 
                              end='')
                        sys.stdout.flush()
                download_list.append(path)
            else:
                print('\rFail to download file {}.'.format(url))
        print('')
        self.raw_data_list = download_list
        print('Raw {} data for {} is ready.'.format(self.projects, 
                                                    self.xena_dtype))
        return self
    
    def transform(self):
        """Transform raw data in a dataset into Xena matrix.
        
        A path for the generated matrix will be assigned to the ``matrix`` 
        property which can be used for making metadata by the ``metadata`` 
        method. Check the ``metadata`` method for details.
        
        The transformation process 1) takes in a list of path for raw data; 2) 
        open each data based on its file extension; 3) read the file object by 
        ``read_func`` and append the readout to a list, which will be 4) 
        assembled into a Xena matrix by ``raws2matrix``. The generated Xena 
        matrix will be saved at the path defined by the ``matrix`` property.
        
        Returns:
            self: allow method chaining.
        """
        
        message = 'Make Xena matrix for {} data of {}.'
        print(message.format(self.xena_dtype, self.projects))
        total = len(self.raw_data_list)
        count = 0
        df_list = []
        for path in self.raw_data_list:
            count = count + 1
            print('\rProcessing {}/{} file...'.format(count, total), end='')
            sys.stdout.flush()
            with read_by_ext(path) as f:
                df = self.read_raw(f)
            df_list.append(df)
        print('\rAll {} files have been processed. '.format(total))
        xena_matrix = self.raws2matrix(df_list)
        # Transformation done
        print('\rSaving matrix to {} ...'.format(self.matrix), end='')
        mkdir_p(self.matrix_dir)
        xena_matrix.to_csv(self.matrix, sep='\t', encoding='utf-8')
        print('\rXena matrix is saved at {}.'.format(self.matrix))
        return self
    
    def metadata(self):
        """Make "metadata.json" for Xena data loading.
        
        One metadata will be created for one Xena matrix defined by the 
        ``matrix`` property. The ``matrix`` property has to point to an 
        existing file when this ``metadata`` method is being called. The 
        metadata JSON file will be saved under the same directory as the 
        matrix file and named with a ".json" postfix appended to the filename 
        of Xena matrix. JSON templates for making metatdata are defined by the 
        ``metadata_template`` property, and variables for rendering the 
        template are defined by the ``metadata_vars`` property.
        
        Returns:
            self: allow method chaining.
        """
        
        message = 'Create metadata for {} data matrix of {}.'
        print(message.format(self.xena_dtype, self.projects))
        try:
            assert os.path.isfile(self.matrix)
        except AttributeError:
            raise IOError('Xena matrix for this dataset is unknown; please '
                          'create a matrix or assign an existing matrix file '
                          'to the "matrix" property before making metadata.')
        except AssertionError:
            raise IOError('{} is not an existing file.'.format(self.matrix))

        # Start to generate metadata.
        # General jinja2 Variables
        print('Creating metadata file ...', end='')
        self._metadata = self.matrix + '.json'
        with open(self._metadata, 'w') as f:
            f.write(self.metadata_template.render(**self.metadata_vars))
        print('\rMetadata JSON is saved at {}.'.format(self._metadata))
        return self


class GDCOmicset(XenaDataset):
    """GDCOmicset is derived from the ``XenaDataset`` class and represents for 
    a Xena matrix whose data is genomic data from GDC.

    This class provides a set of default configurations for downloading and 
    transforming GDC data, as well as generating associated metadata for the 
    transformed Xena matrix.

    Attributes:
        projects (str or list): One (string) or a list of GDC's 
            "cases.project.project_id". All corresponding projects will be 
            included in this dataset.
        xena_dtype (str): A dataset type supported by this class. To get a 
            list of supported types, use ``GDCOmicset.get_supported_dtype()``.
        gdc_filter (dict): A filter dict specifying GDC files relevant to this 
            dataset. Its key is one GDC API available field and the 
            corresponding value should be a string or a list of strings.
        root_dir (str, optional): Defines the root directory for this dataset.
            By default, all files related to this class, such as raw data, 
            Xena matrix, metadata, should and highly recommended to be 
            organized and saved under this directory. The default directory 
            structure is::
            
                root_dir
                └── projects
                    ├── "GDC_Raw_Data"
                    │   └── xena_dtype
                    │       ├── data1
                    │       ├── data2
                    │       ├── ...
                    │       └── dataN
                    └── "Xena_Matrices"
                        ├── projects.xena_dtype(1).tsv
                        ├── projects.xena_dtype(1).tsv.json
                        ├── projects.xena_dtype(2).tsv
                        ├── projects.xena_dtype(2).tsv.json
                        ├── ...
                        ├── projects.xena_dtype(N).tsv
                        └── projects.xena_dtype(N).tsv.json
            
            Defaults to "." which points to current python work directory. 
            Setting "root_dir" will not only set the "root_dir" property but 
            also set "raw_data_dir" and/or "matrix_dir" properties if the 
            directory for GDC data ("raw_data_dir") and/or Xena matrix 
            ("matrix_dir") is not set yet. However, no actual directories will 
            be made by just setting these properties. Directories will only be 
            made when needed.
        
            If you want to reset the default directory structure while 
            "raw_data_dir" and/or "matrix_dir" properties are set already, you 
            can override it using the "set_default_dir_tree" method with 
            "reset_default=True".
        raw_data_dir (str, optional): A path for saving raw data downloaded 
            from GDC. Defaults to None. By default, it will be set by the 
            ``set_default_dir_tree`` method, based on the ``root_dir``.
        raw_data_list (list): A list of file path(s) for all GDC raw data 
            related to this dataset. It will be automatically set by the 
            ``download`` method; or it can be assigned directly as a public 
            attribute. This ``raw_data_list`` attribute will be used by 
            ``transform`` method for making a Xena matrix from GDC raw data. 
            If the ``raw_data_list`` is not usable when trying to get this 
            attribute, the ``raw_data_dir`` property will be checked. All 
            files under ``raw_data_dir`` will be treated as data and used for 
            creating a ``raw_data_list``. 
        matrix_dir (str, optional): A path for saving the Xena matrix for this 
            dataset. Defaults to None. By default, it will be set by the 
            ``set_default_dir_tree`` method, based on the ``root_dir``.
        matrix (str, optional): A path for the Xena matrix of this dataset. 
            This attribute will be used but not validated (i.e. it can be any 
            desired directory and filename) by the ``transform`` method for 
            saving newly generated Xena matrix. This attribute will also be 
            used and validated (i.e. it has to point to a valid file) by the 
            ``metadata`` method for making metadata assciated with the Xena 
            matrix and with this dataset.
            
            By default, when used by the ``transform`` method, this ``matrix`` 
            attribute will adapte a pattern as "projects.xena_type.tsv". There
            are two ways to pass custom matrix name: 1) you can set this 
            ``matrix`` attribute before calling ``transform``. This allows you 
            to customize both directory and filename for the Xena matrix. 
            2) you can pass a ``matrix_name`` argument to the ``transform`` 
            method. You can only customize the name (not its directory) of the 
            Xena matrix.
    """
    
    # Map Xena dtype code to GDC data query dict
    _XENA_GDC_DTYPE = {
            'htseq_counts': {'data_type': 'Gene Expression Quantification',
                             'analysis.workflow_type': 'HTSeq - Counts'},
            'htseq_fpkm': {'data_type': 'Gene Expression Quantification',
                           'analysis.workflow_type': 'HTSeq - FPKM'},
            'htseq_fpkm-uq': {'data_type': 'Gene Expression Quantification',
                              'analysis.workflow_type': 'HTSeq - FPKM-UQ'},
            'mirna': {'data_type': 'miRNA Expression Quantification',
                      'analysis.workflow_type': 'BCGSC miRNA Profiling'},
            'mirna_isoform': {
                    'data_type': 'Isoform Expression Quantification',
                    'analysis.workflow_type': 'BCGSC miRNA Profiling'
                },
            'cnv': {'data_type': 'Copy Number Segment',
                    'analysis.workflow_type': 'DNAcopy'},
            'masked_cnv': {'data_type': 'Masked Copy Number Segment',
                           'analysis.workflow_type': 'DNAcopy',
                           'cases.samples.sample_type_id': 
                               (['01', '02', '03', '04', '05', '06', '07', 
                                 '08', '09', '15', '16', '20', '40', '50', 
                                 '60', '61', '99'])},
            'muse_snv': {
                    'data_type': 'Masked Somatic Mutation',
                    'analysis.workflow_type': 
                        'MuSE Variant Aggregation and Masking'
                },
            'mutect2_snv': {
                    'data_type': 'Masked Somatic Mutation',
                    'analysis.workflow_type': 
                        'MuTect2 Variant Aggregation and Masking'
                },
            'somaticsniper_snv': {
                    'data_type': 'Masked Somatic Mutation',
                    'analysis.workflow_type': 
                        'SomaticSniper Variant Aggregation and Masking'
                },
            'varscan2_snv': {
                    'data_type': 'Masked Somatic Mutation',
                    'analysis.workflow_type': 
                        'VarScan2 Variant Aggregation and Masking'
                },
            'biospecimen': {'data_type': 'Biospecimen Supplement'},
            'clinical': {'data_type': 'Clinical Supplement'}
        }

    # Prefix in filenames for downloaded files
    _GDC_PREFIX = {
            'htseq_counts': 'cases.samples.submitter_id', 
            'htseq_fpkm': 'cases.samples.submitter_id', 
            'htseq_fpkm-uq': 'cases.samples.submitter_id', 
            'mirna': 'cases.samples.submitter_id', 
            'mirna_isoform': 'cases.samples.submitter_id', 
            'cnv': 'cases.samples.submitter_id', 
            'masked_cnv': 'cases.samples.submitter_id', 
            'muse_snv': 'submitter_id', 
            'mutect2_snv': 'submitter_id', 
            'somaticsniper_snv': 'submitter_id', 
            'varscan2_snv': 'submitter_id',
            'biospecimen': 'cases.submitter_id', 
            'clinical': 'cases.submitter_id'
        }
    
    # Settings for making Xena matrix from GDC data
    _READ_RAW_FUNCS = dict.fromkeys(
            ['htseq_counts', 'htseq_fpkm', 'htseq_fpkm-uq'], 
            lambda x: pd.read_table(
                    x, header=None, 
                    names=['Ensembl_ID', 
                           os.path.basename(x.name).split('.', 1)[0]], 
                    index_col=0, usecols=[0, 1], comment='_'
                )
        )
    _READ_RAW_FUNCS['mirna'] = lambda x: pd.read_table(
            x, header=0,
            names=['miRNA_ID', os.path.basename(x.name).split('.', 1)[0]],
            index_col=0, usecols=[0, 2]
        )
    _READ_RAW_FUNCS['mirna_isoform'] = lambda x: pd.read_table(
            x, header=0,
            names=['isoform_coords',
                   os.path.basename(x.name).split('.', 1)[0]],
            index_col=0, usecols=[1, 3]
        )
    _READ_RAW_FUNCS.update(dict.fromkeys(
            ['cnv', 'masked_cnv'], 
            lambda x: pd.read_table(
                    x, header=0, usecols=[1, 2, 3, 5]
                ).assign(sample=os.path.basename(x.name).split('.', 1)[0])
        ))
    _READ_RAW_FUNCS.update(dict.fromkeys(
            ['muse_snv', 'mutect2_snv', 'somaticsniper_snv', 'varscan2_snv'],
            lambda x: pd.read_table(
                    x, header=0, 
                    usecols=[12, 36, 4, 5, 6, 39, 41, 51, 0, 10, 15, 110], 
                    comment='#'
                )
        ))
    _READ_RAW_FUNCS['clinical'] = read_clinical
    _READ_RAW_FUNCS['biospecimen'] = read_biospecimen
    _RAWS2MATRIX_FUNCS = {
            'htseq_counts': rna_columns_matrix,
            'htseq_fpkm': rna_columns_matrix,
            'htseq_fpkm-uq': rna_columns_matrix,
            'mirna': rna_columns_matrix,
            'mirna_isoform': rna_columns_matrix,
            'cnv': 
                lambda x: (pd.concat(x, axis=0)
                             .rename(columns={'Chromosome': 'Chrom',
                                              'Segment_Mean': 'value'})
                             .set_index('sample')),
            'masked_cnv': 
                lambda x: (pd.concat(x, axis=0)
                             .rename(columns={'Chromosome': 'Chrom',
                                              'Segment_Mean': 'value'})
                             .set_index('sample')),
            'muse_snv': snv_maf_matrix,
            'mutect2_snv': snv_maf_matrix,
            'somaticsniper_snv': snv_maf_matrix,
            'varscan2_snv': snv_maf_matrix,
            'biospecimen': 
                lambda x: (pd.concat(x, axis=0)
                             .replace(r'^\s*$', np.nan, regex=True)
                             .dropna(axis=1, how='all')
                             .set_index('bcr_sample_barcode')),
            'clinical': 
                lambda x: (pd.concat(x, axis=0)
                             .replace(r'^\s*$', np.nan, regex=True)
                             .dropna(axis=1, how='all')
                             .set_index('bcr_patient_barcode'))
        }
    
    # Map xena_dtype to corresponding metadata template.
    _METADATA_TEMPLATE = {'htseq_counts': 'template.rna.meta.json',
                          'htseq_fpkm': 'template.rna.meta.json',
                          'htseq_fpkm-uq': 'template.rna.meta.json',
                          'mirna': 'template.mirna.meta.json',
                          'mirna_isoform': 'template.mirna_isoform.meta.json',
                          'cnv': 'template.cnv.meta.json',
                          'masked_cnv': 'template.cnv.meta.json',
                          'muse_snv': 'template.snv.meta.json',
                          'mutect2_snv': 'template.snv.meta.json',
                          'somaticsniper_snv': 'template.snv.meta.json',
                          'varscan2_snv': 'template.snv.meta.json'}
    _METADATA_TEMPLATE = {
            k: os.path.join(os.path.dirname(os.path.abspath(__file__)), 
                            'Resources', v) 
            for k, v in _METADATA_TEMPLATE.items()
        }
    # Map GDC project_id to Xena specific cohort name.
    _XENA_COHORT = {
            'TCGA-BRCA': 'TCGA Breast Cancer (BRCA)',
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
            'TCGA-DLBC': 'TCGA Large B-cell Lymphoma (DLBC)'
        }
    # Jinja2 template variables for corresponding "xena_dtype".
    _METADATA_VARIABLES = {
            'htseq_counts': {'gdc_type': 'HTSeq - Counts',}, 
            'htseq_fpkm': {'gdc_type': 'HTSeq - FPKM',
                           'unit': 'fpkm'}, 
            'htseq_fpkm-uq': {'gdc_type': 'HTSeq - FPKM-UQ',
                              'unit': 'fpkm-uq'}, 
            'mirna': {'gdc_type': 'miRNA Expression Quantification'}, 
            'mirna_isoform': {'gdc_type': 'Isoform Expression Quantification'},
            'cnv': {'gdc_type': 'Copy Number Segment'}, 
            'masked_cnv': {'gdc_type': 'Masked Copy Number Segment'}, 
            'muse_snv': {'gdc_type': 'MuSE Variant Aggregation and Masking'}, 
            'mutect2_snv': {
                    'gdc_type': 'MuTect2 Variant Aggregation and Masking'
                }, 
            'somaticsniper_snv': {
                    'gdc_type': 'SomaticSniper Variant Aggregation and Masking'
                }, 
            'varscan2_snv': {
                    'gdc_type': 'VarScan2 Variant Aggregation and Masking'
                }
        }
    
    @property
    def xena_dtype(self):
        return self.__xena_dtype
    
    @xena_dtype.setter
    def xena_dtype(self, xena_dtype):
        if xena_dtype in self._XENA_GDC_DTYPE:
            self.__xena_dtype = xena_dtype
        else:
            raise ValueError("Unsupported data type: {}".format(xena_dtype))
    
    @classmethod
    def get_supported_dtype(cls):
        """Return a list of dataset type codes supported by this class."""
        
        return cls._XENA_GDC_DTYPE.keys()
    
    # Set default query filter dict for GDC API if it hasn't been set yet.
    @property
    def gdc_filter(self):
        """A filter dict which will be used on GDC's API for querying data 
        files belonging to this dataset.
        
        If the filter dict hasn't been defined, a default filter querying for
        only open access data will be returned. The default filter also 
        contains extra conditions according to "projects" and "xena_dtype" 
        properties.
        """
        
        try:
            assert self.__gdc_filter
            return self.__gdc_filter
        except (AttributeError, AssertionError):
            self.__gdc_filter = {'access': 'open', 
                                 'cases.project.project_id': self.projects}
            self.__gdc_filter.update(self._XENA_GDC_DTYPE[self.xena_dtype])
            return self.__gdc_filter
    
    @gdc_filter.setter
    def gdc_filter(self, filter_dict):
        self.__gdc_filter = filter_dict
    
    # Set default GDC field for prefixing filename of downloaded files.
    @property
    def gdc_prefix(self):
        """A single GDC field whost value will be used as the prefix of the 
        filename for downloaded files.
        
        If this GDC field hasn't been defined, a default field will be used 
        according to the "xena_dtype" property.
        """
        
        try:
            assert self.__gdc_prefix
            return self.__gdc_prefix
        except (AttributeError, AssertionError):
            self.__gdc_prefix = self._GDC_PREFIX[self.xena_dtype]
            return self.__gdc_prefix
    
    @gdc_prefix.setter
    def gdc_prefix(self, gdc_field):
        self.__gdc_prefix = gdc_field
    
    @XenaDataset.download_map.getter
    def download_map(self):
        """Get a dictionary of GDC files to be downloaded for this dataset.
        
        If the dict hasn't been defined, data files will be queried through 
        GDC's API based on the "gdc_filter" property. Check the "gdc_filter" 
        property for details about querying conditions.
        
        Keys of this dict are UUIDs of data files, while values are filenames 
        to be used for saving corresponding files. If the dict hasn't been 
        defined, data files for each sample, by default, will be renamed to 
        "<prefix>.<UUID>.<file extension>" (without brackets), where
        prefix is defined by the "gdc_prefix" property. 
        
        It is worth noting that the data transformation process may need a ID 
        for every data files. The ID, if needed, will be extracted from the 
        filename (the first substring when splitting the filename by "."). For 
        example, Xena uses GDC's "cases.samples.submitter_id" for sample ID. 
        Therefore, data files for each sample will be renamed to 
        "<cases.samples.submitter_id>.<UUID>.<file extension>".
        Please keep that in mind when defining your own download dict. Please 
        check the "gdc_prefix" property and the "transform" method for details.
        """
        
        try:
            assert self.__download_map
            return self.__download_map
        except (AttributeError, AssertionError):
            fields = ['file_id', 'file_name', self.gdc_prefix]
            try:
                print('Searching for raw data ...', end='')
                file_df = gdc.search('files', in_filter=self.gdc_filter, 
                                     fields=fields)
            except Exception:
                file_dict = {}
            else:
                file_df.set_index('file_id', drop=False, inplace=True)
                file_dict = (
                        file_df[self.gdc_prefix].astype(str)
                        + '.' + file_df['file_id'].astype(str)
                        + '.' + file_df['file_name'].apply(gdc.get_ext)
                    ).to_dict()
            if not file_dict:
                msg = '\rNo {} data found for project {}.'
                gdc_dtype = self._XENA_GDC_DTYPE[self.xena_dtype]
                print(msg.format(' - '.join(sorted(gdc_dtype.values())), 
                                 str(self.projects)))
                return None
            file_dict = {
                    '{}/data/{}'.format(gdc.GDC_API_BASE, uuid): 
                        os.path.join(self.raw_data_dir, name) 
                    for uuid, name in file_dict.items()
                }
            self.__download_map = file_dict
            msg = '\r{} files found for {} data of {}.'
            print(msg.format(len(file_dict), self.xena_dtype, self.projects))
            return self.__download_map
    
    @property
    def read_raw(self):
        try:
            return self.__read_raw
        except:
            self.__read_raw = self._READ_RAW_FUNCS[self.xena_dtype]
            return self.__read_raw
    
    @read_raw.setter
    def read_raw(self, func):
        self.__read_raw = func
    
    @property
    def raws2matrix(self):
        try:
            return self.__raws2matrix
        except:
            self.__raws2matrix = self._RAWS2MATRIX_FUNCS[self.xena_dtype]
            return self.__raws2matrix
    
    @raws2matrix.setter
    def raws2matrix(self, func):
        self.__raws2matrix = func
    
    @XenaDataset.metadata_template.getter
    def metadata_template(self):
        try:
            assert isinstance(self.__metadata_template, 
                              jinja2.environment.Template)
            return self.__metadata_template
        except (AttributeError, AssertionError):
            template_json = self._METADATA_TEMPLATE[self.xena_dtype]
            file_dir = os.path.dirname(template_json)
            file_name = os.path.basename(template_json)
            jinja2_env = jinja2.Environment(
                    loader=jinja2.FileSystemLoader(file_dir)
                )
            self.__metadata_template = jinja2_env.get_template(file_name)
            return self.__metadata_template
    
    @property
    def metadata_vars(self):
        try:
            assert (self.__metadata_vars 
                    and isinstance(self.__metadata_vars, dict))
            return self.__metadata_vars
        except (AttributeError, AssertionError):
            matrix_date = time.strftime(
                    "%m-%d-%Y", time.gmtime(os.path.getmtime(self.matrix))
                )
            projects = ','.join(self.projects)
            variables = {'project_id': projects, 
                         'date': matrix_date}
            if projects in self._XENA_COHORT:
                variables['xena_cohort'] = self._XENA_COHORT[projects]
            else:
                variables['xena_cohort'] = projects
            try:
                variables.update(self._METADATA_VARIABLES[self.xena_dtype])
            except KeyError:
                pass
            # Data type specific jinja2 Variables
            if self.xena_dtype in ['muse_snv', 'mutect2_snv', 
                                   'somaticsniper_snv', 'varscan2_snv']:
                try:
                    print('\rSearching the specific URL for raw MAF data ...', 
                          end='')
                    res_df = gdc.search('files', in_filter=self.gdc_filter, 
                                        fields='file_id')
                    if res_df['file_id'].shape == (1,):
                        variables['maf_uuid'] = str(res_df['file_id'][0])
                except Exception:
                    message = ('Fail to get a specific URL for the MAF file '
                               'for: matrix "{}"; "{}" data of cohort "{}".')
                    warnings.warn(message.format(self.matrix, 
                                                 variables['project_id'], 
                                                 variables['xena_cohort']), 
                                  stacklevel=2)
            self.__metadata_vars = variables
            return self.__metadata_vars
    
    @metadata_vars.setter
    def metadata_vars(self, variables):
        self.__metadata_vars = variables


class TCGAPhenoset(XenaDataset):
    """TCGAPhenoset is derived from the ``XenaDataset`` class and represents 
    for a Xena matrix whose data is phenotype data of a TCGA project.

    This class provides a set of default configurations for downloading and 
    transforming TCGA clinical and biospecimen data, as well as generating 
    associated metadata for the transformed Xena matrix.
    """
    
    @property
    def metadata_vars(self):
        try:
            assert (self.__metadata_vars 
                    and isinstance(self.__metadata_vars, dict))
            return self.__metadata_vars
        except (AttributeError, AssertionError):
            matrix_date = time.strftime(
                    "%m-%d-%Y", time.gmtime(os.path.getmtime(self.matrix))
                )
            projects = ','.join(self.projects)
            variables = {'project_id': projects, 
                         'date': matrix_date}
            if projects in GDCOmicset._XENA_COHORT:
                variables['xena_cohort'] = GDCOmicset._XENA_COHORT[projects]
            else:
                variables['xena_cohort'] = projects
            self.__metadata_vars = variables
            return self.__metadata_vars
    
    @metadata_vars.setter
    def metadata_vars(self, variables):
        self.__metadata_vars = variables
    
    def __init__(self, projects, root_dir='.', raw_data_dir=None, 
                 matrix_dir=None):
        self.projects = projects
        self.xena_dtype = 'phenotype'
        self.root_dir = root_dir
        self.clin_dataset = GDCOmicset(self.projects, 'clinical', root_dir, 
                                        raw_data_dir, matrix_dir)
        self.bio_dataset = GDCOmicset(self.projects, 'biospecimen', root_dir,
                                       raw_data_dir, matrix_dir)
        if matrix_dir is not None:
            self.matrix_dir = matrix_dir
        
        self.metadata_template = os.path.join(
                os.path.dirname(os.path.abspath(__file__)),
                'Resources', 'template.phenotype.meta.json'
            )
    
    def download(self):
        """Download GDC's open access phenotype data for projects in this 
        dataset.
        
        There are two types of phenotype data on GDC, the clinical data and 
        the biospecimen data. Data is selected based solely on the "projects" 
        property. Both clinical and biospecimen data will be downloaded. They 
        will be saved under two separated directories. Data from different 
        projects will be put together under one directory as a single dataset, 
        following the same rule as in the ``GDCOmicset`` class. In fact, both 
        downloads are performed by corresponding ``GDCOmicset`` classes. 
        Check the ``download`` method of the ``GDCOmicset`` class for details.
                        
        Returns:
            self: allow method chaining.
        """
        
        # Download clinical data
        self.clin_dataset.download()
        # Download biospecimen data
        self.bio_dataset.download()
        return self
    
    def transform(self):
        """Transform raw phenotype data into Xena matrix.
        
        Raw clinical data and biospecimen data will be transformed separately. 
        Then the clinical matrix and biospecimen matrix will be merged on 
        "cases.submitter_id" and processed properly.
        
        Returns:
            self: allow method chaining.
        """
        
        self.raw_data_list = [self.clin_dataset.transform().matrix, 
                              self.bio_dataset.transform().matrix]
        self.read_raw = pd.read_table
        self.raws2matrix = (
                lambda dfs: pd.merge(
                        dfs[0],
                        dfs[1][dfs[1].columns.difference(dfs[0].columns)
                                     .insert(0, 'bcr_patient_barcode')],
                        how='outer', on='bcr_patient_barcode'
                    ).replace(r'^\s*$', np.nan, regex=True)
                     .fillna(dfs[1])
                     .set_index('bcr_sample_barcode')
            )
        super(TCGAPhenoset, self).transform()
        os.remove(self.clin_dataset.matrix)
        os.remove(self.bio_dataset.matrix)
        return self


class TARGETPhenoset(XenaDataset):
    """TARGETPhenoset is derived from the ``XenaDataset`` class and represents 
    for a Xena matrix whose data is phenotype data of a TARGET project.

    This class provides a set of default configurations for downloading and 
    transforming TARGET clinical, as well as generating associated metadata 
    for the transformed Xena matrix.
    """
    
    @XenaDataset.download_map.getter
    def download_map(self):
        """Get a dictionary of GDC files to be downloaded for this dataset.
        
        If the dict hasn't been defined, data files will be queried through 
        GDC's API based on the "gdc_filter" property. Check the "gdc_filter" 
        property for details about querying conditions.
        
        Keys of this dict are UUIDs of data files, while values are filenames 
        to be used for saving corresponding files. If the dict hasn't been 
        defined, data files for each sample, by default, will be renamed to 
        "<prefix>.<UUID>.<file extension>" (without brackets), where
        prefix is defined by the "gdc_prefix" property. 
        
        It is worth noting that the data transformation process may need a ID 
        for every data files. The ID, if needed, will be extracted from the 
        filename (the first substring when splitting the filename by "."). For 
        example, Xena uses GDC's "cases.samples.submitter_id" for sample ID. 
        Therefore, data files for each sample will be renamed to 
        "<cases.samples.submitter_id>.<UUID>.<file extension>".
        Please keep that in mind when defining your own download dict. Please 
        check the "gdc_prefix" property and the "transform" method for details.
        """
        
        try:
            assert self.__download_map
            return self.__download_map
        except (AttributeError, AssertionError):
            filters = {'access': 'open', 
                       'cases.project.project_id': self.projects, 
                       'data_type': 'Clinical Supplement'}
            fields = ['file_id', 'file_name']
            try:
                print('Searching for raw data ...', end='')
                file_df = gdc.search('files', in_filter=filters, 
                                     fields=fields)
            except Exception:
                file_dict = {}
            else:
                file_df.set_index('file_id', drop=False, inplace=True)
                file_dict = (
                        file_df['file_id'].astype(str)
                        + '.' + file_df['file_name'].astype(str)
                    ).to_dict()
            if not file_dict:
                msg = '\rNo {} data found for project {}.'
                gdc_dtype = self._XENA_GDC_DTYPE[self.xena_dtype]
                print(msg.format(' - '.join(sorted(gdc_dtype.values())), 
                                 str(self.projects)))
                return None
            file_dict = {
                    '{}/data/{}'.format(gdc.GDC_API_BASE, uuid): 
                        os.path.join(self.raw_data_dir, name) 
                    for uuid, name in file_dict.items()
                }
            self.__download_map = file_dict
            msg = '\r{} files found for {} data of {}.'
            print(msg.format(len(file_dict), self.xena_dtype, self.projects))
            return self.__download_map
    
    @property
    def metadata_vars(self):
        try:
            assert (self.__metadata_vars 
                    and isinstance(self.__metadata_vars, dict))
            return self.__metadata_vars
        except (AttributeError, AssertionError):
            matrix_date = time.strftime(
                    "%m-%d-%Y", time.gmtime(os.path.getmtime(self.matrix))
                )
            projects = ','.join(self.projects)
            variables = {'project_id': projects, 
                         'date': matrix_date,
                         'xena_cohort': projects}
            self.__metadata_vars = variables
            return self.__metadata_vars
    
    @metadata_vars.setter
    def metadata_vars(self, variables):
        self.__metadata_vars = variables
    
    def __target_clin_dfs2matrix(self, df_list):
        print('Merging into one matrix ...', end='')
        clin_df = pd.concat(df_list)
        print('\rMapping clinical info to individual samples...', end='')
        cases = gdc.search(
                'cases', in_filter={'project.project_id': self.projects}, 
                fields=['submitter_id', 'samples.submitter_id']
            )
        cases_samples = [c for c in cases if 'samples' in c]
        from pandas.io.json import json_normalize
        cases_samples_map = json_normalize(
                cases_samples, 'samples', ['submitter_id'], 
                meta_prefix='cases.'
            ).rename(columns={'submitter_id': 'sample_id', 
                              'cases.submitter_id': 'TARGET USI'})
        return pd.merge(clin_df.reset_index(), cases_samples_map, how='inner', 
                        on='TARGET USI').set_index('sample_id')
    
    def __init__(self, projects, root_dir='.', raw_data_dir=None, 
                 matrix_dir=None):
        self.xena_dtype = 'clinical'
        self.read_raw = read_clinical
        self.raws2matrix = self.__target_clin_dfs2matrix
        super(TARGETPhenoset, self).__init__(projects, 'clinical', 
                                                 root_dir, raw_data_dir, 
                                                 matrix_dir)
        
        self.metadata_template = os.path.join(
                os.path.dirname(os.path.abspath(__file__)),
                'Resources', 'template.phenotype.meta.json'
            )


def main():
    print('A python module of Xena specific importing pipeline for GDC data.')


if __name__ == '__main__':
    main()
