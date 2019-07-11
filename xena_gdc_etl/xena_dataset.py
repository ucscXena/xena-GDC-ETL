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

import functools
import time
import os
import re
import sys
import warnings

import jinja2
from lxml import etree
import numpy as np
import pandas as pd

from xena_gdc_etl import gdc
from .utils import mkdir_p, requests_retry_session, reduce_json_array
from .constants import (
    GDC_XENA_COHORT,
    METADATA_TEMPLATE,
    METADATA_VARIABLES,
    GDC_RELEASE_URL,
)


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

    if hasattr(fileobj, 'name'):
        filename = fileobj.name
    else:
        filename = fileobj
    ext = os.path.splitext(filename)[1]
    if ext == '.xlsx':
        # Design specifically for TARGET biospecimen
        try:
            df = pd.read_excel(
                filename, sheet_name='Sample Names', header=None
            )
        except Exception:
            try:
                df = pd.read_excel(
                    filename, sheet_name='SampleNames', header=None
                )
            except Exception:
                raise
        df.iloc[0].fillna(method='ffill', inplace=True)
        df.columns = df.iloc[0:2].apply(lambda x: x.str.cat(sep='.'))
        return df.drop(df.index[0:2]).set_index(df.columns[0])
    elif ext != '.xml':
        raise IOError('Unknown file type for biospecimen data: {}'.format(ext))

    root = etree.parse(fileobj).getroot()
    ns = root.nsmap
    assert (
        'biospecimen'
        in root.xpath('@xsi:schemaLocation', namespaces=ns)[0].lower()
    )
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
    samples_common['primary_diagnosis'] = gdc.TCGA_STUDY_ABBR[
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
    sample_mask = df.bcr_sample_barcode.map(lambda s: s[-3:-1] not in ['10'])
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

    if hasattr(fileobj, 'name'):
        filename = fileobj.name
    else:
        filename = fileobj
    ext = os.path.splitext(filename)[1]
    if ext == '.xlsx':
        xl_file = pd.ExcelFile(filename)
        sheets = xl_file.sheet_names
        if sheets[0] == "Clinical Data":
            return xl_file.parse(sheets[0], index_col=0)
        else:
            print("Clincal Data not found, skipping this file ...")
            return pd.DataFrame()
    elif ext != '.xml':
        raise IOError('Unknown file type for clinical data: {}'.format(ext))

    root = etree.parse(fileobj).getroot()
    ns = root.nsmap
    assert (
        'clinical'
        in root.xpath('@xsi:schemaLocation', namespaces=ns)[0].lower()
    )
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
            [
                child.text.strip()
                for child in root.find('.//clin_shared:race_list', ns)
                if child.text and child.text.strip()
            ]
        )
    except Exception:
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


def merge_cnv(filelist):
    """Transform GDC's CNV data into Xena data matrix.

    Args:
        filelist (list of path): The list of input raw MAF file for mutation
            data.

    Returns:
        pandas.core.frame.DataFrame: Transformed pandas DataFrame.
    """

    xena_matrix = pd.DataFrame()
    total = len(filelist)
    count = 0
    for path in filelist:
        xena_matrix = xena_matrix.append(
            pd.read_csv(path, sep="\t", header=0, usecols=[1, 2, 3, 5]).assign(
                sample=os.path.basename(path).split('.', 1)[0]
            )
        )
        count += 1
        print('\rProcessed {}/{} file...'.format(count, total), end='')
        sys.stdout.flush()
    print('\rAll {} files have been processed. '.format(total))
    return xena_matrix.rename(
        columns={'Chromosome': 'Chrom', 'Segment_Mean': 'value'}
    ).set_index('sample')


def snv_maf_matrix(filelist):
    """Transform GDC's MAF data into Xena data matrix.

    A new column of DNA variant allele frequencies named "dna_vaf" will
    calculated by division "t_alt_count" / "t_depth". Columns "t_alt_count"
    and "t_depth" will then be dropped. In the end, columns will be renamed
    accordingly and row index will be set as sample ID.

    Args:
        filelist (list of path): The list of input raw MAF file for mutation
            data.

    Returns:
        pandas.core.frame.DataFrame: Transformed pandas DataFrame.
    """

    assert len(filelist) == 1
    print('\rProcessing 1/1 file...', end='')
    df = pd.read_csv(
        filelist[0],
        sep="\t",
        header=0,
        comment='#',
        usecols=[12, 36, 4, 5, 6, 39, 41, 51, 0, 10, 15, 110],
    )
    print('\rCalculating "dna_vaf" ...', end='')
    df['dna_vaf'] = df['t_alt_count'] / df['t_depth']
    print('\rTrim "Tumor_Sample_Barcode" into Xena sample ID ...', end='')
    df['Tumor_Sample_Barcode'] = df['Tumor_Sample_Barcode'].apply(
        lambda x: '-'.join(x.split('-', 4)[0:4])
    )
    print('\rRe-organizing matrix ...', end='')
    df = (
        df.drop(['t_alt_count', 't_depth'], axis=1)
        .set_index('Tumor_Sample_Barcode')
        .rename(
            columns={
                'Hugo_Symbol': 'gene',
                'Chromosome': 'chrom',
                'Start_Position': 'start',
                'End_Position': 'end',
                'Reference_Allele': 'ref',
                'Tumor_Seq_Allele2': 'alt',
                'Tumor_Sample_Barcode': 'sampleid',
                'HGVSp_Short': 'Amino_Acid_Change',
                'Consequence': 'effect',
                'FILTER': 'filter',
            }
        )
    )
    df.index.name = 'Sample_ID'
    return df


def merge_sample_cols(
    filelist,
    header='infer',
    index_col=0,
    usecols=[0, 1],
    comment=None,
    index_name='id',
    get_sid=lambda f: os.path.basename(f).split('.', 1)[0],
    log2TF=True,
):
    """Merge and process a list of raw data files to make a Xena data matrix.

    Each file will be considered as data from one sample and will be extracted
    as one column in the merged matrix. Data from the same sample (identified
    by sample ID) will be averaged before being put into the matrix. By
    default (``log2TF=True``), merged will be transformed by log2(x + 1).

    Args:
        filelist (list of path): The list of input raw data.
        header:
        index_col:
        usecols:
        comment:
        index_name:
        get_sid:
        log2TF:

    Returns:
        pandas.core.frame.DataFrame: Ready to load Xena matrix.
    """

    # Group data by sample (sid), especially for samples with more than 1 data
    # files (repeats). If repeats within a sample need to be averaged, this
    # will be very helpful.
    sample_dict = {}
    for path in filelist:
        sample_id = get_sid(path)  # os.path.basename(path).split('.', 1)[0]
        if sample_id not in sample_dict:
            sample_dict[sample_id] = []
        sample_dict[sample_id].append(path)
    # Read raw data sample by sample, average data for each sample if needed,
    # and put data for each sample into the matrix
    xena_matrix = pd.DataFrame()
    total = len(filelist)
    count = 0
    for sample_id in sample_dict:
        df_list = [
            pd.read_csv(
                f,
                sep="\t",
                header=header,
                index_col=index_col,
                usecols=usecols,
                comment=comment,
                names=[index_name, sample_id],
            )
            for f in sample_dict[sample_id]
        ]
        if len(sample_dict[sample_id]) > 1:
            df_list = [
                xena_matrix,
                pd.concat(df_list, axis=1, copy=False)
                .mean(1)
                .rename(sample_id),
            ]
        else:
            df_list.insert(0, xena_matrix)
        xena_matrix = pd.concat(df_list, axis=1, copy=False)
        count += len(sample_dict[sample_id])
        print('\rProcessed {}/{} file...'.format(count, total), end='')
        sys.stdout.flush()
    print('\rAll {} files have been processed. '.format(total))
    xena_matrix.index.name = index_name
    if log2TF:
        return np.log2(xena_matrix + 1)
    else:
        return xena_matrix


def handle_gistic(filelist):
    """Handles GISTIC Data type.
    Args:
        filelist (list of path): The list of input raw data.
    Returns:
        pandas.core.frame.DataFrame: Ready to load Xena matrix.
    """

    assert len(filelist) == 1
    print('\rProcessing file {}'.format(filelist[0]), end='')
    df = pd.read_csv(
        filelist[0],
        sep="\t",
        header=0,
        comment='#',
        index_col=0,
    )
    df = df.drop(["Gene ID", "Cytoband"], axis=1)
    columns = list(df)
    mapping = gdc.map_two_fields(
        endpoint="cases",
        input_field="samples.portions.analytes.aliquots.aliquot_id",
        output_field="samples.submitter_id",
        input_values=columns,
    )
    mapping = reduce_json_array(mapping)
    df = df.rename(columns=mapping)
    return df


class XenaDataset(object):
    r"""XenaDataset represents for one Xena matrix in a Xena cohort.

    This class provides a set of method for downloading and transforming data
    into Xena matrix, as well as generating associated metadata. It also
    provides a set of attributes to control these processes.

    Attributes:
        projects (str or list of str): One or a list of project IDs describing
            study projects included in this dataset.
        xena_dtype (str): A short string (like an ID) describing the type of
            data in this dataset. It is highly recommended, though not
            required, that one dataset has a single type of data.
        root_dir (str): Defines the root directory for this dataset. The
            XenaDataset and the importing process can be highly customized,
            with directories for every data and each step explicitly assigned.
            You can set directories for raw data (through the ``raw_data_dir``
            property) and Xena matrix (through the ``matrix_dir``property)
            specifically. The ``root_dir`` will be essentially useless under
            such situation.
            If some or all directories remain unassigned when being used, a
            default directory tree will be used, with a structure like this::

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

            By default, all files related to a dataset, such as raw data, Xena
            matrix, metadata, should and highly recommended to be organized
            and saved under the root directory. However, this is neither
            required nor checked. Setting directory related properties
            (including ``root_dir`` and some properties listed below) will not
            trigger creations of any directories. Directories, if not exist,
            will only be made right before being needed.

            Defaults to "." which points to current python work directory.
        raw_data_dir (str): A directory for raw data. Please Check the
            ``raw_data_list`` property for its potential usage for defining
            raw data for Xena matrix ``transform``, and check the ``root_dir``
            property for the default "Raw_Data" directory structure. Defaults
            to None.
        matrix_dir (str): A path for saving the transformed Xena matrix for
            this dataset. If the ``matrix_dir`` is not available at the time
            of being used, the ``matrix`` property will be checked first. If
            the ``matrix`` property is available, its directory will be
            assigned to the ``matrix_dir``. If the ``matrix`` property is not
            available, a default path will be assigned according to the
            default directory structure. Check the ``root_dir`` property for
            the default directory structure. Defaults to None.
        download_map (dict): A dict with the key being a URL for one raw data
            to be downloaded and the value being a path for saving downloaded
            raw data.
        raw_data_list (list): A list of file path(s) for all raw data related
            to this dataset. It will be automatically set by the ``download``
            method; or it can be assigned directly as a public attribute. This
            ``raw_data_list`` attribute will be used by the ``transform``
            method for making a Xena matrix from raw data. If the
            ``raw_data_list`` is not available at the time of being used, the
            ``raw_data_dir`` property will be checked. All files under
            ``raw_data_dir`` will be treated as data and used for creating a
            ``raw_data_list``.
        raws2matrix (callable): A function used for merging multiple raw data
            in the ``raw_data_list`` into one Xena matrix, as well as
            processing the merged matrix if needed. A valid ``raws2matrix``
            must accept only one argument, which is ``raw_data_list``.
        matrix (str): A path for the Xena matrix of this dataset.
            This attribute will be used but not validated by the ``transform``
            method for saving newly generated Xena matrix. This attribute will
            also be used yet will be validated (i.e. it has to point to a
            valid existing file) by the ``metadata`` method before making
            metadata associated with the Xena matrix and with this dataset. If
            ``matrix`` is not available at the time of being used by the
            ``transform`` method, it will use ``matrix_dir`` as its directory
            and will adapte a filename with the "projects.xena_type.tsv"
            pattern. Please check the ``matrix_dir`` property to see how it is
            determined.
        metadata_template (jinja2.environment.Template or str): A Jinja2
            template for rendering metadata of this dataset. When using a
            string to set ``metadata_template``, it will be used as a path to
            the template file and the corresponding template will be retrieved
            and assigned to this attribute.
        metadata_vars (dict): A dict of variables which will be used (by \*\*
            unpacking) for rendering the ``metadata_template``.
    """

    @property
    def projects(self):
        return self._projects

    @projects.setter
    def projects(self, projects):
        if isinstance(projects, str):
            self._projects = [projects]
        elif isinstance(projects, list):
            self._projects = projects
        else:
            raise ValueError('"projects" must be either str or list.')

    @property
    def root_dir(self):
        """A path of an existing directory for keeping files (raw data, matrix
        and metadata) of this dataset.
        """

        return self._root_dir

    @root_dir.setter
    def root_dir(self, path):
        if os.path.isdir(path):
            self._root_dir = os.path.abspath(path)
        else:
            raise IOError('{} is not an existing directory.'.format(path))

    @property
    def raw_data_dir(self):
        try:
            return self._raw_data_dir
        except AttributeError:
            self._raw_data_dir = os.path.join(
                self.root_dir,
                '_'.join(self.projects),
                'Raw_Data',
                str(self.xena_dtype),
            )
            return self._raw_data_dir

    @raw_data_dir.setter
    def raw_data_dir(self, path):
        self._raw_data_dir = os.path.abspath(path)

    @property
    def matrix_dir(self):
        try:
            return self._matrix_dir
        except AttributeError:
            try:
                self._matrix_dir = os.path.dirname(self._matrix)
                return self._matrix_dir
            except AttributeError:
                self._matrix_dir = os.path.join(
                    self.root_dir, '_'.join(self.projects), 'Xena_Matrices'
                )
                return self._matrix_dir

    @matrix_dir.setter
    def matrix_dir(self, path):
        self._matrix_dir = os.path.abspath(path)

    @property
    def download_map(self):
        assert self._download_map and isinstance(self._download_map, dict)
        return self._download_map

    @download_map.setter
    def download_map(self, d):
        if isinstance(d, dict):
            self._download_map = d
        else:
            raise TypeError(
                'download_map should be a dict, ' 'not {}'.format(type(d))
            )

    # Raw data list: try to get the list from ``raw_data_dir`` if unavailable.
    @property
    def raw_data_list(self):
        try:
            return self._raw_data_list
        except AttributeError:
            try:
                raw_data_dir = os.path.abspath(self.raw_data_dir)
                raw_data = []
                for f in os.listdir(raw_data_dir):
                    f_path = os.path.join(raw_data_dir, f)
                    if os.path.isfile(f_path):
                        raw_data.append(f_path)
                if raw_data:
                    self._raw_data_list = raw_data
                else:
                    raise ValueError
            except Exception:
                raise ValueError('Cannot find raw data.')
            return self._raw_data_list

    @raw_data_list.setter
    def raw_data_list(self, raw_data):
        self._raw_data_list = raw_data

    @property
    def matrix(self):
        try:
            assert self._matrix
            return self._matrix
        except (AttributeError, AssertionError):
            self._matrix = os.path.join(
                self.matrix_dir,
                '{}.{}.tsv'.format('_'.join(self.projects), self.xena_dtype),
            )
            return self._matrix

    @matrix.setter
    def matrix(self, path):
        self._matrix = os.path.abspath(path)
        self.matrix_dir = os.path.dirname(self._matrix)

    @property
    def metadata_template(self):
        assert isinstance(self._metadata_template, jinja2.environment.Template)
        return self._metadata_template

    @metadata_template.setter
    def metadata_template(self, template):
        if isinstance(template, jinja2.environment.Template):
            self._metadata_template = template
        elif isinstance(template, str):
            jinja2_env = jinja2.Environment(
                loader=jinja2.PackageLoader('xena_gdc_etl', 'resources')
            )
            self._metadata_template = jinja2_env.get_template(template)
        else:
            raise TypeError(
                'metadata_template should be a jinja2 template or an existing '
                'path to a template JSON file, not a :{}.'.format(
                    type(template)
                )
            )

    def __init__(
        self,
        projects,
        xena_dtype,
        root_dir='.',
        raw_data_dir=None,
        matrix_dir=None,
    ):
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
            response = requests_retry_session().get(url, stream=True)
            if response.ok:
                path = os.path.abspath(path)
                status = '\r[{:d}/{:d}] Downloading to "{}" ...'
                print(status.format(count, total, path), end='')
                sys.stdout.flush()
                mkdir_p(os.path.dirname(path))
                with open(path, 'wb') as f:
                    for chunk in response.iter_content(chunk_size):
                        f.write(chunk)
                download_list.append(path)
            else:
                raise IOError(
                    '\nFail to download file {}. Response {}'.format(
                        url, response.status_code
                    )
                )
        print('')
        self.raw_data_list = download_list
        print(
            'Raw {} data for {} is ready.'.format(
                self.projects, self.xena_dtype
            )
        )
        return self

    def transform(self):
        """Transform raw data in a dataset into Xena matrix.

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
        xena_matrix = self.raws2matrix(self.raw_data_list)
        # Transformation done
        print('\rSaving matrix to {} ...'.format(self.matrix), end='')
        mkdir_p(self.matrix_dir)
        xena_matrix.to_csv(self.matrix, sep='\t', encoding='utf-8')
        print('\rXena matrix is saved at {}.'.format(self.matrix))
        return self

    def metadata(self):
        """Make "metadata.json" for Xena data loading.

        A JSON of metadata will be created for the Xena matrix defined by the
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
            raise IOError(
                'Xena matrix for this dataset is unknown; please create a '
                'matrix or assign an existing matrix file to the "matrix" '
                'property before making metadata.'
            )
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
    r"""GDCOmicset is derived from the ``XenaDataset`` class and represents for
    a Xena matrix whose data is genomic data from GDC.

    This class provides a set of default configurations for downloading and
    transforming GDC data, as well as generating associated metadata for the
    transformed Xena matrix. These default configurations are stored as
    private constants, and they can be checked and/or changed through the
    following attributes: ``gdc_release``, ``gdc_filter``, ``gdc_prefix``,
    ``download_map``, ``raws2matrix``, ``metadata_template``, and
    ``metadata_vars``.

    Attributes:
        projects (str or list): One (string) or a list of GDC's
            "cases.project.project_id". All corresponding projects will be
            included in this dataset.
        xena_dtype (str): A dataset type supported by this class. To get a
            list of supported types, use ``GDCOmicset.get_supported_dtype()``.
        gdc_release (str): URL to the data release note for the dataset. It
            will be used by the ``metadata`` method when making the metadata
            for this dataset. It is highly recommended that this attribute is
            set explicitly by the user so that it is guaranteed to match the
            data (raw data) underlying this dataset. If it is not available,
            the most recent data release will be queried and used.
        gdc_filter (dict): A filter for querying GDC data underlying this
            dataset. Each item of this dict means to be an "in" operation,
            with its key being one GDC API available field and its value being
            a string or a list of strings. It can be automatically derived
            from ``projects`` and ``xena_dtype`` if it is not assigned
            explicitly by the user when being used. Please check `GDC API
            documentation
            <https://docs.gdc.cancer.gov/API/Users_Guide/Search_and_Retrieval/#filters-specifying-the-query>`_
            for details.
        gdc_prefix (str): A GDC available file field whost value will be used
            in the filename of corresponding download file. It will be used by
            ``download_map`` for making default download map. It can be
            automatically mapped from ``xena_dtype`` if it is not assigned
            explicitly by the user when being used. Please check
            ``download_map`` and `GDC API documentation
            <https://docs.gdc.cancer.gov/API/Users_Guide/Search_and_Retrieval/#filters-specifying-the-query>`_
            for details.
        download_map (dict): A dict with the key being a URL for one raw data
            to be downloaded and the value being a path for saving downloaded
            raw data. If it hasn't been assigned explicitly by the user when
            being used, it can be automatically generated by querying through
            GDC API according to ``gdc_filter`` and ``gdc_prefix`` which are
            based on ``projects`` and ``xena_dtype``. Please check
            ``gdc_filter`` for details about querying conditions. Filename of
            data files, by default, will adapt a pattern of
            "<value of gdc_prefix>.<GDC file UUID>.<file extension>"

            It is worth noting that the data transformation process may need
            an ID for every data files. The ``raws2matrix`` functions may
            extract the ID from the filename (the first substring when
            splitting the filename by "."). For example, Xena uses GDC's
            "cases.samples.submitter_id" for sample ID. Therefore,
            ``gdc_prefix`` should be set to "cases.samples.submitter_id" so
            that data files for each sample will be renamed to
            "<cases.samples.submitter_id>.<file UUID>.<file extension>",
            allowing the desired sample ID to be extracted correctly. Please
            keep that in mind when trying to define your own download dict but
            use default transformation settings (``raws2matrix``). Please
            check the ``raws2matrix`` properties, as well as the ``transform``
            method for details.
        raws2matrix (callable): A function which accepts only one argument of
            ``raw_data_list``, merges them into one Xena matrix, and processes
            the merged matrix if needed. Defaults, if needed, can be mapped
            from ``xena_dtype``.
        metadata_template (jinja2.environment.Template or str): A Jinja2
            template for rendering metadata of this dataset. When setting this
            attribute with a string, it will be taken as a path to the
            template file and the corresponding template will be retrieved and
            assigned to this attribute. Defaults, if needed, can be mapped
            from ``xena_dtype``.
        metadata_vars (dict): A dict of variables which will be used (by \*\*
            unpacking) when rendering the ``metadata_template``. Defaults, if
            needed, can be derived from corresponding matrix and ``projects``
            and ``xena_dtype`` properties.
    """

    # Map Xena dtype code to GDC data query dict
    _XENA_GDC_DTYPE = {
        'htseq_counts': {
            'data_type': 'Gene Expression Quantification',
            'analysis.workflow_type': 'HTSeq - Counts',
        },
        'htseq_fpkm': {
            'data_type': 'Gene Expression Quantification',
            'analysis.workflow_type': 'HTSeq - FPKM',
        },
        'htseq_fpkm-uq': {
            'data_type': 'Gene Expression Quantification',
            'analysis.workflow_type': 'HTSeq - FPKM-UQ',
        },
        'mirna': {
            'data_type': 'miRNA Expression Quantification',
            'analysis.workflow_type': 'BCGSC miRNA Profiling',
        },
        'mirna_isoform': {
            'data_type': 'Isoform Expression Quantification',
            'analysis.workflow_type': 'BCGSC miRNA Profiling',
        },
        'cnv': {
            'data_type': 'Copy Number Segment',
            'analysis.workflow_type': 'DNAcopy',
        },
        'masked_cnv': {
            'data_type': 'Masked Copy Number Segment',
            'analysis.workflow_type': 'DNAcopy',
            'cases.samples.sample_type_id': (
                [
                    '01',
                    '02',
                    '03',
                    '04',
                    '05',
                    '06',
                    '07',
                    '08',
                    '09',
                    '15',
                    '16',
                    '20',
                    '40',
                    '50',
                    '60',
                    '61',
                    '99',
                ]
            ),
        },
        'muse_snv': {
            'data_type': 'Masked Somatic Mutation',
            'analysis.workflow_type': 'MuSE Variant Aggregation and Masking',
        },
        'mutect2_snv': {
            'data_type': 'Masked Somatic Mutation',
            'analysis.workflow_type':
            'MuTect2 Variant Aggregation and Masking',
        },
        'somaticsniper_snv': {
            'data_type': 'Masked Somatic Mutation',
            'analysis.workflow_type':
            'SomaticSniper Variant Aggregation and Masking',
        },
        'varscan2_snv': {
            'data_type': 'Masked Somatic Mutation',
            'analysis.workflow_type':
            'VarScan2 Variant Aggregation and Masking',
        },
        'gistic': {
            'data_type': 'Gene Level Copy Number Scores',
            'analysis.workflow_type': 'GISTIC - Copy Number Score',
        },
        'star_counts': {
            'analysis.workflow_type': 'STAR - Counts',
        },
        'methylation27': {
            'data_type': 'Methylation Beta Value',
            'platform': 'Illumina Human Methylation 27',
        },
        'methylation450': {
            'data_type': 'Methylation Beta Value',
            'platform': 'Illumina Human Methylation 450',
        },
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
        'gistic': 'submitter_id',
        'star_counts': 'cases.samples.submitter_id',
        'methylation27': 'cases.samples.submitter_id',
        'methylation450': 'cases.samples.submitter_id',
    }

    # Settings for making Xena matrix from GDC data
    _RAWS2MATRIX_FUNCS = dict.fromkeys(
        ['htseq_counts', 'htseq_fpkm', 'htseq_fpkm-uq'],
        functools.partial(
            merge_sample_cols, header=None, index_name='Ensembl_ID'
        ),
    )
    _RAWS2MATRIX_FUNCS['mirna'] = functools.partial(
        merge_sample_cols, header=0, usecols=[0, 2], index_name='miRNA_ID'
    )
    _RAWS2MATRIX_FUNCS['mirna_isoform'] = functools.partial(
        merge_sample_cols,
        header=0,
        usecols=[1, 3],
        index_name='isoform_coords',
    )
    _RAWS2MATRIX_FUNCS.update(dict.fromkeys(['cnv', 'masked_cnv'], merge_cnv))
    _RAWS2MATRIX_FUNCS.update(
        dict.fromkeys(
            ['muse_snv', 'mutect2_snv', 'somaticsniper_snv', 'varscan2_snv'],
            snv_maf_matrix,
        )
    )
    _RAWS2MATRIX_FUNCS['gistic'] = handle_gistic
    _RAWS2MATRIX_FUNCS['star_counts'] = functools.partial(
        merge_sample_cols,
        header=0,
        index_name='Ensembl_ID',
    )
    _RAWS2MATRIX_FUNCS.update(
        dict.fromkeys(
            ['methylation27', 'methylation450'],
            functools.partial(
                merge_sample_cols,
                header=0,
                log2TF=False,
                index_name='Composite Element REF',
            ),
        )
    )

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

    @property
    def gdc_release(self):
        try:
            return self.__gdc_release
        except AttributeError:
            data_release = gdc.search('status', typ='json')['data_release']
            anchor = (
                re.match(r'(Data Release [^\s]+)\s', data_release)
                .group(1)
                .replace(' ', '-')
                .replace('.', '')
                .lower()
            )
            self.__gdc_release = GDC_RELEASE_URL + '#' + anchor
            return self.__gdc_release

    @gdc_release.setter
    def gdc_release(self, url):
        self.__gdc_release = url

    # Set default query filter dict for GDC API if it hasn't been set yet.
    @property
    def gdc_filter(self):
        try:
            assert self.__gdc_filter
            return self.__gdc_filter
        except (AttributeError, AssertionError):
            self.__gdc_filter = {
                'access': 'open',
                'cases.project.project_id': self.projects,
            }
            self.__gdc_filter.update(self._XENA_GDC_DTYPE[self.xena_dtype])
            return self.__gdc_filter

    @gdc_filter.setter
    def gdc_filter(self, filter_dict):
        self.__gdc_filter = filter_dict

    # Set default GDC field for prefixing filename of downloaded files.
    @property
    def gdc_prefix(self):
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
        try:
            assert self._download_map
            return self._download_map
        except (AttributeError, AssertionError):
            fields = ['file_id', 'file_name', self.gdc_prefix]
            try:
                print('Searching for raw data ...', end='')
                file_df = gdc.search(
                    'files', in_filter=self.gdc_filter, fields=fields
                )
            except Exception:
                file_dict = {}
            else:
                file_df.set_index('file_id', drop=False, inplace=True)
                file_dict = (
                    file_df[self.gdc_prefix].astype(str)
                    + '.'
                    + file_df['file_id'].astype(str)
                    + '.'
                    + file_df['file_name'].apply(gdc.get_ext)
                ).to_dict()
            if not file_dict:
                msg = '\rNo {} data found for project {}.'
                gdc_dtype = self._XENA_GDC_DTYPE[self.xena_dtype]
                print(
                    msg.format(
                        ' - '.join(sorted(gdc_dtype.values())),
                        str(self.projects),
                    )
                )
                return file_dict
            file_dict = {
                '{}/data/{}'.format(gdc.GDC_API_BASE, uuid): os.path.join(
                    self.raw_data_dir, name
                )
                for uuid, name in file_dict.items()
            }
            self._download_map = file_dict
            msg = '\r{} files found for {} data of {}.'
            print(msg.format(len(file_dict), self.xena_dtype, self.projects))
            return self._download_map

    @property
    def raws2matrix(self):
        try:
            return self.__raws2matrix
        except Exception:
            self.__raws2matrix = self._RAWS2MATRIX_FUNCS[self.xena_dtype]
            return self.__raws2matrix

    @raws2matrix.setter
    def raws2matrix(self, func):
        self.__raws2matrix = func

    @XenaDataset.metadata_template.getter
    def metadata_template(self):
        try:
            assert isinstance(
                self._metadata_template, jinja2.environment.Template
            )
            return self._metadata_template
        except (AttributeError, AssertionError):
            template_json = METADATA_TEMPLATE[self.xena_dtype]
            jinja2_env = jinja2.Environment(
                loader=jinja2.PackageLoader('xena_gdc_etl', 'resources')
            )
            self._metadata_template = jinja2_env.get_template(template_json)
            return self._metadata_template

    @property
    def metadata_vars(self):
        try:
            assert self.__metadata_vars and isinstance(
                self.__metadata_vars, dict
            )
            return self.__metadata_vars
        except (AttributeError, AssertionError):
            matrix_date = time.strftime(
                "%m-%d-%Y", time.gmtime(os.path.getmtime(self.matrix))
            )
            projects = ','.join(self.projects)
            variables = {
                'project_id': projects,
                'date': matrix_date,
                'gdc_release': self.gdc_release,
            }
            if projects in GDC_XENA_COHORT:
                variables['xena_cohort'] = GDC_XENA_COHORT[projects]
            else:
                variables['xena_cohort'] = 'GDC ' + projects
            try:
                variables.update(METADATA_VARIABLES[self.xena_dtype])
            except KeyError:
                pass
            # Data type specific jinja2 Variables
            if self.xena_dtype in [
                'muse_snv',
                'mutect2_snv',
                'somaticsniper_snv',
                'varscan2_snv',
            ]:
                try:
                    print(
                        '\rSearching the specific URL for raw MAF data ...',
                        end='',
                    )
                    res_df = gdc.search(
                        'files', in_filter=self.gdc_filter, fields='file_id'
                    )
                    if res_df['file_id'].shape == (1,):
                        variables['maf_uuid'] = str(res_df['file_id'][0])
                except Exception:
                    message = (
                        'Fail to get a specific URL for the MAF file '
                        'for: matrix "{}"; "{}" data of cohort "{}".'
                    )
                    warnings.warn(
                        message.format(
                            self.matrix,
                            variables['project_id'],
                            variables['xena_cohort'],
                        ),
                        stacklevel=2,
                    )
            self.__metadata_vars = variables
            return self.__metadata_vars

    @metadata_vars.setter
    def metadata_vars(self, variables):
        self.__metadata_vars = variables


class GDCPhenoset(XenaDataset):
    r"""GDCPhenoset is derived from the ``XenaDataset`` class and represents for
    a Xena matrix whose data is phenotype data from GDC.

    This class provides a set of default configurations for downloading and
    transforming GDC data, as well as generating associated metadata for the
    transformed Xena matrix. These default configurations are stored as
    private constants, and they can be checked and/or changed through the
    following attributes: ``gdc_release``, ``gdc_filter``, ``download_map``,
    ``raws2matrix``, ``metadata_template``, and ``metadata_vars``.

    Attributes:
        projects (str or list): One (string) or a list of GDC's
            "cases.project.project_id". All corresponding projects will be
            included in this dataset.
        xena_dtype (str): One dataset type of "biospecimen", "clinical",
            "raw_phenotype" or "GDC_phenotype". Defaults to None, for which
            the class will guess the correct type to use from ``projects``.
        gdc_release (str): URL to the data release note for the dataset. It
            will be used by the ``metadata`` method when making the metadata
            for this dataset. It is highly recommended that this attribute is
            set explicitly by the user so that it is guaranteed to match the
            data (raw data) underlying this dataset. If it is not available,
            the most recent data release will be queried and used.
        gdc_filter (dict): A filter for querying GDC data underlying this
            dataset. Each item of this dict means to be an "in" operation,
            with its key being one GDC API available field and its value being
            a string or a list of strings. It can be automatically derived
            from ``projects`` and ``xena_dtype`` if it is not assigned
            explicitly by the user when being used. Please check `GDC API
            documentation
            <https://docs.gdc.cancer.gov/API/Users_Guide/Search_and_Retrieval/#filters-specifying-the-query>`_
            for details.
        download_map (dict): A dict with the key being a URL for one raw data
            to be downloaded and the value being a path for saving downloaded
            raw data. If it hasn't been assigned explicitly by the user when
            being used, it can be automatically generated by querying through
            GDC API according to ``gdc_filter`` which are based on
            ``projects`` and ``xena_dtype``. Filename of data files, by
            default, will adapt a pattern of
            "<data_category>.<GDC file UUID>.<file extension>"

            It is worth noting the "<data_category>" prefix can be useful or
            even necessary for ``transform`` method to apply correct
            transformation to the file. "<data_category>" is closely related
            to the format of the file.
        metadata_template (jinja2.environment.Template or str): A Jinja2
            template for rendering metadata of this dataset. When setting this
            attribute with a string, it will be taken as a path to the
            template file and the corresponding template will be retrieved and
            assigned to this attribute. Defaults, if needed, can be mapped
            from ``xena_dtype``.
        metadata_vars (dict): A dict of variables which will be used (by \*\*
            unpacking) when rendering the ``metadata_template``. Defaults, if
            needed, can be derived from corresponding matrix and ``projects``
            and ``xena_dtype`` properties.
    """

    # Map Xena dtype code to GDC data query dict
    _XENA_GDC_DTYPE = {
        'biospecimen': {
            'data_category': 'Biospecimen',
            'data_format': ['BCR XML', 'XLSX'],
        },
        'clinical': {
            'data_category': 'Clinical',
            'data_format': ['BCR XML', 'XLSX'],
        },
        'raw_phenotype': {
            'data_category': ['Biospecimen', 'Clinical'],
            'data_format': ['BCR XML', 'XLSX'],
        },
        'GDC_phenotype': {
            'data_category': ['Biospecimen', 'Clinical'],
            'data_format': ['BCR XML', 'XLSX'],
        },
    }
    # To resovle overlapping between raw data and API data, remove columns
    # according to the following lists.
    _API_DROPS = [
        'id',
        'case_id',
        'state',
        'created_datetime',
        'updated_datetime',
        'demographic_id.demographic',
        'submitter_id.demographic',
        'state.demographic',
        'created_datetime.demographic',
        'updated_datetime.demographic',
        'diagnosis_id.diagnoses',
        'submitter_id.diagnoses',
        'state.diagnoses',
        'created_datetime.diagnoses',
        'updated_datetime.diagnoses',
        'treatment_id.treatments.diagnoses',
        'submitter_id.treatments.diagnoses',
        'state.treatments.diagnoses',
        'created_datetime.treatments.diagnoses',
        'updated_datetime.treatments.diagnoses',
        'exposure_id.exposures',
        'submitter_id.exposures',
        'state.exposures',
        'created_datetime.exposures',
        'updated_datetime.exposures',
        'pathology_report_uuid.samples',
        'state.project',
        'released.project',
        'sample_id.samples',
        'created_datetime.samples',
        'updated_datetime.samples',
        'tissue_source_site_id.tissue_source_site',
    ]
    _RAW_DROPS = [
        'alcohol_history_documented',
        'bcr_patient_barcode',
        'bcr_patient_uuid',
        'bcr_sample_uuid',
        'composition',
        'current_weight',
        'days_to_birth',
        'days_to_collection',
        'days_to_death',
        'days_to_last_followup',
        'days_to_sample_procurement',
        'ethnicity',
        'freezing_method',
        'gender',
        'height',
        'icd_10',
        'icd_o_3_histology',
        'icd_o_3_site',
        'initial_weight',
        'intermediate_dimension',
        'is_ffpe',
        'longest_dimension',
        'oct_embedded',
        'pathologic_stage',
        'pathology_report_uuid',
        'preservation_method',
        'primary_diagnosis',
        'race',
        'sample_type',
        'sample_type_id',
        'shortest_dimension',
        'state',
        'time_between_clamping_and_freezing',
        'time_between_excision_and_freezing',
        'tissue_type',
        'tumor_descriptor',
        'tumor_tissue_site',
        'vital_status',
    ]

    @property
    def xena_dtype(self):
        return self.__xena_dtype

    @xena_dtype.setter
    def xena_dtype(self, xena_dtype):
        if xena_dtype in self._XENA_GDC_DTYPE:
            self.__xena_dtype = xena_dtype
        else:
            raise ValueError("Unsupported data type: {}".format(xena_dtype))

    @property
    def gdc_release(self):
        try:
            return self.__gdc_release
        except AttributeError:
            data_release = gdc.search('status', typ='json')['data_release']
            anchor = (
                re.match(r'(Data Release [^\s]+)\s', data_release)
                .group(1)
                .replace(' ', '-')
                .replace('.', '')
                .lower()
            )
            self.__gdc_release = GDC_RELEASE_URL + '#' + anchor
            return self.__gdc_release

    @gdc_release.setter
    def gdc_release(self, url):
        self.__gdc_release = url

    # Set default query filter dict for GDC API if it hasn't been set yet.
    @property
    def gdc_filter(self):
        try:
            assert self.__gdc_filter
            return self.__gdc_filter
        except (AttributeError, AssertionError):
            self.__gdc_filter = {
                'access': 'open',
                'cases.project.project_id': self.projects,
            }
            self.__gdc_filter.update(self._XENA_GDC_DTYPE[self.xena_dtype])
            return self.__gdc_filter

    @gdc_filter.setter
    def gdc_filter(self, filter_dict):
        self.__gdc_filter = filter_dict

    @XenaDataset.download_map.getter
    def download_map(self):
        try:
            assert self._download_map
            return self._download_map
        except (AttributeError, AssertionError):
            fields = ['file_id', 'file_name', 'data_category']
            try:
                print('Searching for raw data ...', end='')
                file_df = gdc.search(
                    'files', in_filter=self.gdc_filter, fields=fields
                )
            except Exception:
                file_dict = {}
            else:
                file_df.set_index('file_id', drop=False, inplace=True)
                file_dict = (
                    file_df['data_category'].astype(str)
                    + '.'
                    + file_df['file_id'].astype(str)
                    + '.'
                    + file_df['file_name'].apply(gdc.get_ext)
                ).to_dict()
            if not file_dict:
                msg = '\rNo {} data found for project {}.'
                gdc_dtype = self._XENA_GDC_DTYPE[self.xena_dtype]
                print(
                    msg.format(
                        ' - '.join(sorted(gdc_dtype.values())),
                        str(self.projects),
                    )
                )
                return file_dict
            file_dict = {
                '{}/data/{}'.format(gdc.GDC_API_BASE, uuid): os.path.join(
                    self.raw_data_dir, name
                )
                for uuid, name in file_dict.items()
            }
            self._download_map = file_dict
            msg = '\r{} files found for {} data of {}.'
            print(msg.format(len(file_dict), self.xena_dtype, self.projects))
            return self._download_map

    @property
    def metadata_vars(self):
        try:
            assert self.__metadata_vars and isinstance(
                self.__metadata_vars, dict
            )
            return self.__metadata_vars
        except (AttributeError, AssertionError):
            matrix_date = time.strftime(
                "%m-%d-%Y", time.gmtime(os.path.getmtime(self.matrix))
            )
            projects = ','.join(self.projects)
            variables = {
                'project_id': projects,
                'date': matrix_date,
                'gdc_release': self.gdc_release,
            }
            if projects in GDC_XENA_COHORT:
                variables['xena_cohort'] = GDC_XENA_COHORT[projects]
            else:
                variables['xena_cohort'] = 'GDC ' + projects
            self.__metadata_vars = variables
            return self.__metadata_vars

    @metadata_vars.setter
    def metadata_vars(self, variables):
        self.__metadata_vars = variables

    def __init__(
        self,
        projects,
        xena_dtype=None,
        root_dir='.',
        raw_data_dir=None,
        matrix_dir=None,
    ):
        self.projects = projects
        if xena_dtype is not None:
            self.xena_dtype = xena_dtype
        elif all([i.startswith('TCGA-') for i in self.projects]):
            self.xena_dtype = 'GDC_phenotype'
        elif all([i.startswith('TARGET-') for i in self.projects]):
            self.xena_dtype = 'clinical'
        else:
            warnings.warn(
                'Caution: fail to guess phenotype data type for project '
                '{}; use "raw_phenotype" as default.'.format(self.projects)
            )
            self.xena_dtype = 'raw_phenotype'
        self.root_dir = root_dir
        if matrix_dir is not None:
            self.matrix_dir = matrix_dir
        jinja2_env = jinja2.Environment(
            loader=jinja2.PackageLoader('xena_gdc_etl', 'resources')
        )
        self.metadata_template = jinja2_env.get_template(
            'template.phenotype.meta.json'
        )

    def transform(self):
        """Transform raw phenotype data into Xena matrix.

        Raw clinical data and/or biospecimen data will first be transformed
        separately. Then if needed (i.e. for TCGA projects) the clinical
        matrix and biospecimen matrix will be merged on "cases.submitter_id"
        and processed properly.

        Returns:
            self: allow method chaining.
        """

        message = 'Make Xena matrix for {} data of {}.'
        print(message.format(self.xena_dtype, self.projects))
        total = len(self.raw_data_list)
        count = 0
        bio_dfs = []
        clin_dfs = []
        for path in self.raw_data_list:
            count = count + 1
            print('\rProcessing {}/{} file...'.format(count, total), end='')
            sys.stdout.flush()
            # `read_biospecimen` and `read_clinical` will check file format
            try:
                df = read_clinical(path)
                if not df.empty:
                    clin_dfs.append(df)
            except Exception:
                try:
                    df = read_biospecimen(path)
                    bio_dfs.append(df)
                except Exception:
                    raise TypeError('Fail to process file {}.'.format(path))
        print('\rAll {} files have been processed. '.format(total))
        try:
            bio_matrix = (
                pd.concat(bio_dfs, axis=0)
                .replace(r'\r\n', ' ', regex=True)
                .replace(r'^\s*$', np.nan, regex=True)
                .dropna(axis=1, how='all')
            )
        except Exception:
            bio_matrix = pd.DataFrame()
        try:
            clin_matrix = (
                pd.concat(clin_dfs, axis=0)
                .replace(r'\r\n', ' ', regex=True)
                .replace(r'^\s*$', np.nan, regex=True)
                .dropna(axis=1, how='all')
            )
        except Exception:
            clin_matrix = pd.DataFrame()
        if self.xena_dtype == 'clinical':
            try:
                xena_matrix = clin_matrix.set_index('bcr_patient_barcode')
            except Exception:
                xena_matrix = clin_matrix
            print('\rMapping clinical info to individual samples...', end='')
            cases = gdc.search(
                'cases',
                in_filter={'project.project_id': self.projects},
                fields=['submitter_id', 'samples.submitter_id'],
                typ='json',
            )
            cases_samples = [c for c in cases if 'samples' in c]
            from pandas.io.json import json_normalize

            cases_samples_map = json_normalize(
                cases_samples,
                'samples',
                ['submitter_id'],
                meta_prefix='cases.',
            )
            if all([i.startswith('TCGA-') for i in self.projects]):
                cases_samples_map = cases_samples_map.rename(
                    columns={
                        'submitter_id': 'sample_id',
                        'cases.submitter_id': 'bcr_patient_barcode',
                    }
                )
                xena_matrix = pd.merge(
                    xena_matrix.reset_index(),
                    cases_samples_map,
                    how='inner',
                    on='bcr_patient_barcode',
                ).set_index('sample_id')
            elif all([i.startswith('TARGET-') for i in self.projects]):
                cases_samples_map = cases_samples_map.rename(
                    columns={
                        'submitter_id': 'sample_id',
                        'cases.submitter_id': 'TARGET USI',
                    }
                )
                xena_matrix = pd.merge(
                    xena_matrix.reset_index(),
                    cases_samples_map,
                    how='inner',
                    on='TARGET USI',
                ).set_index('sample_id')
            else:
                warnings.warn('Fail to get per sample based clinical matrix.')
        elif self.xena_dtype == 'biospecimen':
            try:
                xena_matrix = bio_matrix.set_index('bcr_sample_barcode')
            except Exception:
                xena_matrix = bio_matrix
        elif self.xena_dtype in ['raw_phenotype', 'GDC_phenotype'] and all(
            [i.startswith('TCGA-') for i in self.projects]
        ):
            bio_columns = bio_matrix.columns.difference(
                clin_matrix.columns
            ).insert(0, 'bcr_patient_barcode')
            xena_matrix = (
                pd.merge(
                    clin_matrix,
                    bio_matrix[bio_columns],
                    how='outer',
                    on='bcr_patient_barcode',
                )
                .replace(r'^\s*$', np.nan, regex=True)
                .set_index('bcr_sample_barcode')
                .fillna(bio_matrix.set_index('bcr_sample_barcode'))
            )
        if self.xena_dtype == 'GDC_phenotype':
            # Query GDC API for GDC harmonized phenotype info
            api_clin = gdc.get_samples_clinical(self.projects)
            # Revert hierarchy order in column names
            api_clin = api_clin.rename(
                columns={
                    n: '.'.join(reversed(n.split('.')))
                    for n in api_clin.columns
                }
            )
            if all([i.startswith('TCGA-') for i in self.projects]):
                # Remove code 10, Blood Derived Normal, sample:
                # https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes
                sample_mask = api_clin['submitter_id.samples'].map(
                    lambda s: s[-3:-1] not in ['10']
                )
                api_clin = api_clin[sample_mask].set_index(
                    'submitter_id.samples'
                )
                # Remove all empty columns
                api_clin = api_clin.dropna(axis=1, how='all')
                # For overlapping columns between raw data matrix and GDC'S
                # API data matrix, use API data.
                for c in self._API_DROPS:
                    try:
                        api_clin.drop(c, axis=1, inplace=True)
                    except Exception:
                        pass
                for c in self._RAW_DROPS:
                    try:
                        xena_matrix.drop(c, axis=1, inplace=True)
                    except Exception:
                        pass
                # Merge phenotype matrix from raw data and that from GDC's API
                xena_matrix = xena_matrix.reset_index().rename(
                    columns={'bcr_sample_barcode': 'submitter_id.samples'}
                )
                xena_matrix = (
                    pd.merge(
                        xena_matrix,
                        api_clin.reset_index(),
                        how='outer',
                        on='submitter_id.samples',
                    )
                    .replace(r'^\s*$', np.nan, regex=True)
                    .set_index('submitter_id.samples')
                )
            elif all([i.startswith('TARGET-') for i in self.projects]):
                xena_matrix = api_clin.dropna(axis=1, how='all').set_index(
                    'submitter_id.samples'
                )
            else:
                raise ValueError(
                    'Getting "GDC_phenotype" for a cohort with mixed TCGA and '
                    'TARGET projects is not currently suppported.'
                )
        # Transformation done
        print('\rSaving matrix to {} ...'.format(self.matrix), end='')
        mkdir_p(self.matrix_dir)
        xena_matrix.to_csv(self.matrix, sep='\t', encoding='utf-8')
        print('\rXena matrix is saved at {}.'.format(self.matrix))
        return self


class GDCSurvivalset(XenaDataset):
    r"""GDCSurvivalset is derived from the ``XenaDataset`` class and represents
    for a Xena matrix of GDC survival data for project(s) of interest.

    In general, survival data is retrieved from GDC API's "analysis/survival"
    endpoint. This class provides two default configurations, which can be
    checked and/or changed through ``gdc_release`` and ``metadata_vars``, for
    generating metadata for the transformed Xena matrix. The ``download`` and
    ``transform`` methods are overridden by methods specific for GDC survival
    data.

    Attributes:
        gdc_release (str): URL to the data release note for the dataset. It
            will be used by the ``metadata`` method when making the metadata
            for this dataset. It is highly recommended that this attribute is
            set explicitly by the user so that it is guaranteed to match the
            GDC data underlying this dataset. If it is not available, the most
            recent data release will be queried and used.
        metadata_vars (dict): A dict of variables which will be used (by \*\*
            unpacking) when rendering the ``metadata_template``. Defaults, if
            needed, can be derived from corresponding matrix and the
            ``projects`` property.
    """

    @property
    def gdc_release(self):
        try:
            return self.__gdc_release
        except AttributeError:
            data_release = gdc.search('status', typ='json')['data_release']
            anchor = (
                re.match(r'(Data Release [^\s]+)\s', data_release)
                .group(1)
                .replace(' ', '-')
                .replace('.', '')
                .lower()
            )
            self.__gdc_release = GDC_RELEASE_URL + '#' + anchor
            return self.__gdc_release

    @gdc_release.setter
    def gdc_release(self, url):
        self.__gdc_release = url

    @property
    def metadata_vars(self):
        try:
            assert self.__metadata_vars and isinstance(
                self.__metadata_vars, dict
            )
            return self.__metadata_vars
        except (AttributeError, AssertionError):
            matrix_date = time.strftime(
                "%m-%d-%Y", time.gmtime(os.path.getmtime(self.matrix))
            )
            projects = ','.join(self.projects)
            variables = {
                'project_id': projects,
                'date': matrix_date,
                'gdc_release': self.gdc_release,
            }
            if projects in GDC_XENA_COHORT:
                variables['xena_cohort'] = GDC_XENA_COHORT[projects]
            else:
                variables['xena_cohort'] = 'GDC ' + projects
            self.__metadata_vars = variables
            return self.__metadata_vars

    @metadata_vars.setter
    def metadata_vars(self, variables):
        self.__metadata_vars = variables

    def __init__(
        self, projects, root_dir='.', raw_data_dir=None, matrix_dir=None
    ):
        super(GDCSurvivalset, self).__init__(
            projects, 'survival', root_dir, raw_data_dir, matrix_dir
        )

        jinja2_env = jinja2.Environment(
            loader=jinja2.PackageLoader('xena_gdc_etl', 'resources')
        )
        self.metadata_template = jinja2_env.get_template(
            'template.survival.meta.json'
        )

    def download(self):
        """Retrieve GDC API's survival data for project(s) in this dataset.

        The survival data is queried and retrieved through GDC API's
        "analysis/survival" endpoint for project(s) belonging to this dataset.
        JSON query results are converted to a pandas DataFrame and saved as a
        single tab-separated values ("<projects>.GDC_survival.tsv") file under
        ``raw_data_dir``.

        Returns:
            self: allow method chaining.
        """

        survival = gdc.search(
            'analysis/survival',
            in_filter={'project.project_id': self.projects},
            typ='json',
        )['results'][0]['donors']
        mkdir_p(self.raw_data_dir)
        path = os.path.join(
            self.raw_data_dir,
            '{}.GDC_survival.tsv'.format(','.join(self.projects)),
        )
        pd.DataFrame(survival).set_index('id').to_csv(path, sep='\t')
        self.raw_data_list = [path]
        print(
            'Raw {} data for {} is ready.'.format(
                self.projects, self.xena_dtype
            )
        )
        return self

    def transform(self):
        """Transform GDC survival data according to Xena survival data spec

        Only 1 GDC raw survival data (i.e. ``raw_data_list[0]``) will be read
        and used by this transformation. Xena survival data has 4 columns,
        which are "sample", "OS", "OS.time" and "_PATIENT". "OS"
        corresponds to the "censored" column in GDC survival data; "OS.time"
        corresponds to the "time" column in GDC survival data;
        "_PATIENT" corresponds to the "submitter_id" column in GDC survival
        data which is the case(patient)'s submitter ID; "sample" contains
        "samples.submitter_id" for corresponding case(patient).

        Returns:
            self: allow method chaining.
        """

        raw_df = pd.read_csv(self.raw_data_list[0], sep="\t")
        # Transform GDC data according to Xena survival data spec
        survival_df = raw_df.drop(
            ['project_id', 'survivalEstimate'], axis=1
        ).rename(
            columns={
                'censored': 'OS',
                'time': 'OS.time',
                'submitter_id': '_PATIENT',
            }
        )
        survival_df['OS'] = (~survival_df['OS']).map(int)
        # Get samples to case map
        case_samples = gdc.search(
            'cases',
            in_filter={'project.project_id': self.projects},
            fields='submitter_sample_ids',
            typ='json',
        )
        case_samples = [c for c in case_samples if 'submitter_sample_ids' in c]
        samples_df = pd.io.json.json_normalize(
            case_samples, 'submitter_sample_ids', 'id'
        ).rename(columns={0: 'sample'})
        sample_mask = samples_df['sample'].map(
            lambda s: s[-3:-1] not in ['10']
        )
        samples_df = samples_df[sample_mask]
        # Make sample indexed survival matrix
        df = (
            pd.merge(survival_df, samples_df, how='inner', on='id')
            .drop('id', axis=1)
            .set_index('sample')
        )
        mkdir_p(os.path.dirname(self.matrix))
        df.to_csv(self.matrix, sep='\t')
        print('\rXena matrix is saved at {}.'.format(self.matrix))
        return self


def main():
    print('A python module of Xena specific importing pipeline for GDC data.')
    start = time.time()
    dataset = GDCOmicset(
        projects='TCGA-BRCA',
        root_dir=r'/mnt/e/GitHub/xena-GDC-ETL/gitignore/test',
        xena_dtype='methylation450',
    )
    dataset.download().transform().metadata()
    print(time.time() - start)


if __name__ == '__main__':
    main()
