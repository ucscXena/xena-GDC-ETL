xena-GDC-ETL
============

Extract, transform and load `GDC data <https://portal.gdc.cancer.gov/>`__ onto `UCSC Xena <https://xenabrowser.net/>`__.

**Table of Contents**

- `Dependencies <#dependencies>`__
- `Installation <#installation>`__
- `Usage <#usage>`__
- `ETL process details <#etl-process-details>`__
- `API <#api>`__

Dependencies
------------

Specific versions mentioned below have been tested. Eariler versions may still work but not guaranteed. 

1. Python 2.7

   This pipeline has been tested with python 2.7. It may also work with python 3 since it was originally designed to be `single-source Python 2/3 compatible <https://docs.python.org/3/howto/pyporting.html#the-short-explanation>`__.

2. `Jinja2 <http://jinja.pocoo.org/docs/2.9/>`__ v2.7.2
3. `Numpy <http://www.numpy.org/>`__ v1.13.0
4. `Pandas <http://pandas.pydata.org/>`__ v0.20.2
5. `Requests <http://docs.python-requests.org/en/master/>`__ v1.2.3

Installation
------------

-  ``git clone https://github.com/yunhailuo/xena-GDC-ETL.git``
-  Dependencies can be installed either before or after cloning this repository.

Usage
-----

- **Build basic importing pipeline**

  1. A pipeline starts from a ``XenaDataset`` object which describes one dataset of a cohort on Xena. A ``XenaDataset`` object is defined by "projects" and "xena\_dtype". For example, 

     .. code:: python

       from xena_dataset import XenaDataset
       cohort = XenaDataset('TCGA-BRCA', 'htsep.counts')

     A "projects" is one (as string) or a list of GDC's "`cases.project.project\_id <https://docs.gdc.cancer.gov/API/Users_Guide/Appendix_A_Available_Fields/#file-fields>`__\ ". You can get some basic info about projects on GDC by 

     .. code:: python

       import gdc
       projects = gdc.get_all_project_info()
       print('List of GDC's project-ids:')
       print(projects.keys())

     A "xena\_dtype" describes the specific type of data in this dataset. It is critical for downloading corresponding GDC data, for transforming raw GDC data into a valid Xena matrix and for generating metadata. You can get a list of "xena\_dtype" codes by 

     .. code:: python

       from xena_dataset import XenaDataset
       print(XenaDataset.get_supported_dtype())

  2. An importing pipeline can do 3 jobs either separately or combinatorially. They are:

     +------------------------+------------------------------------------------------------------------------+
     | ``XenaDataset`` method | Description                                                                  |
     +========================+==============================================================================+
     | ``download_gdc``       | Query for and download GDC open access data by "projects" and "xena\_dtype". |
     +------------------------+------------------------------------------------------------------------------+
     | ``transform``          | Transform raw GDC data into a Xena compatible matrix.                        |
     +------------------------+------------------------------------------------------------------------------+
     | ``metadata``           | Generate a Xena compatible metadata for a Xena matrix of a dataset.          |
     +------------------------+------------------------------------------------------------------------------+

     All 3 methods return the original "XenaDataset" object which allows method chaining. For example, the following code will perform an import from scatch and save all data (within a default directory structure) under the current working directory: 

     .. code:: python

       cohort.download_gdc().transform().metadata()

     The 3 methods can be used separately or in specific combinations, as long as corresponding properties are set properly. The next section will describe details about the 3 ``XenaDataset`` methods.

- **Details about** ``XenaDataset``

  - **Default directory structure for data**

    There are 3 types of data for in a typical Xena dataset: raw data (GDC), Xena matrix and metadata of this dataset. By default, they will be organized like this:

    ::

       root_dir
       └── projects
           ├── "GDC_Raw_Data"
           │   └── xena_dtype with "." replaced by "_"
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

    You may change this by setting "root\_dir", "raw\_data\_dir" and "matrix\_dir" properties of the "XenaDataset" object. Default "root\_dir" is the current working directory.

  - **The** ``download_gdc`` **method**

    This method looks for GDC data relevant to this dataset, filtering by "projects" and "xena\_dtype". Files containing data for one individual sample are renamed as ".<UUID>.<file extension>". Files containing data for the whole dataset are renamed as "<UUID>.<file extension>" They will be saved under the directory defined by "raw\_data\_dir", and the "raw\_data\_list" property of this "XenaDataset" object will be set to a list of paths for downloaded files.

  - **The** ``transform`` **method**

    This method works on a list of data defined by the "raw\_data\_list" property of this "XenaDataset" object. "raw\_data\_list" can be set directly; or it will be set by the ``download_gdc`` method if it succeed. Data in this list will be merged and/or transformed into a valid Xena matrix based on the "xena\_dtype" of this dataset. The filename and location for the final Xena matrix is defined by the "matrix" property of this "XenaDataset" object. By default, Xena matrix will be saved under the "matrix\_dir" as "..tsv".

  - **The** ``metadata`` **method**

    This method works on the Xena matrix defined by the "matrix" property of this "XenaDataset" object. "matrix" can be set directly; or it will be set by the ``transform`` method if it succeed. Based on the "xena\_dtype" and "projects", a jinja2 template will be selected and variables in the template will be set accordingly. Generated metadata is a JSON file. Its name will be derived from the name of Xena matrix by adding a ".json" postfix. This metadata will be saved under the same directory as the Xena matrix.

ETL process details
-------------------

- **Settings for downloading raw data from GDC**

  +-------------------+-----------------------------------+-----------------------------------------------+------------------+
  |                   |                               GDC filter                                          |                  |
  +    xena_dtype     +-----------------------------------+-----------------------------------------------+ File count/Level +
  |                   | data_type                         | analysis.workflow_type                        |                  |
  +===================+===================================+===============================================+==================+
  | htseq.counts      | Gene Expression Quantification    | HTSeq - Counts                                | 1/Sample vial    |
  +-------------------+-----------------------------------+-----------------------------------------------+------------------+
  | htseq.fpkm        | Gene Expression Quantification    | HTSeq - FPKM                                  | 1/Sample vial    |
  +-------------------+-----------------------------------+-----------------------------------------------+------------------+
  | htseq.fpkm-uq     | Gene Expression Quantification    | HTSeq - FPKM-UQ                               | 1/Sample vial    |
  +-------------------+-----------------------------------+-----------------------------------------------+------------------+
  | mirna             | miRNA Expression Quantification   | BCGSC miRNA Profiling                         | 1/Sample vial    |
  +-------------------+-----------------------------------+-----------------------------------------------+------------------+
  | mirna.isoform     | Isoform Expression Quantification | BCGSC miRNA Profiling                         | 1/Sample vial    |
  +-------------------+-----------------------------------+-----------------------------------------------+------------------+
  | cnv               | Copy Number Segment               | DNAcopy                                       | 1/Sample vial    |
  +-------------------+-----------------------------------+-----------------------------------------------+------------------+
  | masked.cnv        | Masked Copy Number Segment        | DNAcopy                                       | 1/Sample vial    |
  +-------------------+-----------------------------------+-----------------------------------------------+------------------+
  | muse.snv          | Masked Somatic Mutation           | MuSE Variant Aggregation and Masking          | 1/Project        |
  +-------------------+-----------------------------------+-----------------------------------------------+------------------+
  | mutect2.snv       | Masked Somatic Mutation           | MuTect2 Variant Aggregation and Masking       | 1/Project        |
  +-------------------+-----------------------------------+-----------------------------------------------+------------------+
  | somaticsniper.snv | Masked Somatic Mutation           | SomaticSniper Variant Aggregation and Masking | 1/Project        |
  +-------------------+-----------------------------------+-----------------------------------------------+------------------+
  | varscan2.snv      | Masked Somatic Mutation           | VarScan2 Variant Aggregation and Masking      | 1/Project        |
  +-------------------+-----------------------------------+-----------------------------------------------+------------------+

- **Settings for transform raw data into Xena matrix**

  +-------------------+----------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------+-----------------------+-------------------------------+-----------------------------------------------------------------------------+
  |  xena_dtype       | Raw data has header? | Select columns (in order)                                                                                                                                                  | Row index       | Skip rows start with? | Merge into matrix as          | Process matrix                                                              |
  +===================+======================+============================================================================================================================================================================+=================+=======================+===============================+=============================================================================+
  | htseq.counts      | No                   | 1, 2                                                                                                                                                                       | Ensembl_ID      | _                     | 1 new column based on index   | 1. Average if there are multiple data from the same sample vial;            |
  |                   |                      | [Ensembl_ID, Counts]                                                                                                                                                       |                 |                       |                               | 2. log2(counts + 1)                                                         |
  +-------------------+----------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------+-----------------------+-------------------------------+-----------------------------------------------------------------------------+
  | htseq.fpkm        | No                   | 1, 2                                                                                                                                                                       | Ensembl_ID      | _                     | 1 new column based on index   | 1. Average if there are multiple data from the same sample vial;            |
  |                   |                      | [Ensembl_ID, Counts]                                                                                                                                                       |                 |                       |                               | 2. log2(counts + 1)                                                         |
  +-------------------+----------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------+-----------------------+-------------------------------+-----------------------------------------------------------------------------+
  | htseq.fpkm-uq     | No                   | 1, 2                                                                                                                                                                       | Ensembl_ID      | _                     | 1 new column based on index   | 1. Average if there are multiple data from the same sample vial;            |
  |                   |                      | [Ensembl_ID, Counts]                                                                                                                                                       |                 |                       |                               | 2. log2(counts + 1)                                                         |
  +-------------------+----------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------+-----------------------+-------------------------------+-----------------------------------------------------------------------------+
  | mirna             | Yes                  | 1, 3                                                                                                                                                                       | miRNA_ID        | N/A                   | 1 new column based on index   | 1. Average if there are multiple data from the same sample vial;            |
  |                   |                      | [miRNA_ID, RPM]                                                                                                                                                            |                 |                       |                               | 2. log2(counts + 1)                                                         |
  +-------------------+----------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------+-----------------------+-------------------------------+-----------------------------------------------------------------------------+
  | mirna.isoform     | Yes                  | 2, 4                                                                                                                                                                       | isoform_coords  | N/A                   | 1 new column based on index   | 1. Average if there are multiple data from the same sample vial;            |
  |                   |                      | [isoform_coords, RPM]                                                                                                                                                      |                 |                       |                               | 2. log2(counts + 1)                                                         |
  +-------------------+----------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------+-----------------------+-------------------------------+-----------------------------------------------------------------------------+
  | cnv               | Yes                  | 2, 3, 4, 6                                                                                                                                                                 | sample          | N/A                   | New rows based on column name | 1. Rename columns as::                                                      |
  |                   |                      | [Chromosome, Start, End, Segment_Mean]                                                                                                                                     |                 |                       |                               |                                                                             |
  |                   |                      |                                                                                                                                                                            |                 |                       |                               |     {                                                                       |
  |                   |                      |                                                                                                                                                                            |                 |                       |                               |         'Chromosome': 'Chrom',                                              |
  |                   |                      |                                                                                                                                                                            |                 |                       |                               |         'Segment_Mean': 'value'                                             |
  |                   |                      |                                                                                                                                                                            |                 |                       |                               |     }                                                                       |
  +-------------------+----------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------+-----------------------+-------------------------------+-----------------------------------------------------------------------------+
  | masked.cnv        | Yes                  | 1, 2, 3, 5                                                                                                                                                                 | sample          | N/A                   | New rows based on column name | 1. Rename columns as::                                                      |
  |                   |                      | [Chromosome, Start, End, Segment_Mean]                                                                                                                                     |                 |                       |                               |                                                                             |
  |                   |                      |                                                                                                                                                                            |                 |                       |                               |     {                                                                       |
  |                   |                      |                                                                                                                                                                            |                 |                       |                               |         'Chromosome': 'Chrom',                                              |
  |                   |                      |                                                                                                                                                                            |                 |                       |                               |         'Segment_Mean': 'value'                                             |
  |                   |                      |                                                                                                                                                                            |                 |                       |                               |     }                                                                       |
  +-------------------+----------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------+-----------------------+-------------------------------+-----------------------------------------------------------------------------+
  | muse.snv          | Yes                  | 13, 37, 5, 6, 7, 40, 42, 52, 1, 11, 16, 111                                                                                                                                | N/A             | #                     | N/A                           | 1. Calculate variant allele frequency (dna_vaf) by "t_alt_count"/"t_depth"; |
  | mutect2.snv       |                      | [Tumor_Seq_Allele2, HGVSp_Short, Chromosome, Start_Position, End_Position, t_depth, t_alt_count, Consequence, Hugo_Symbol, Reference_Allele, Tumor_Sample_Barcode, FILTER] |                 |                       |                               | 2. Delete "t_alt_count" and "t_depth" columns;                              |
  | somaticsniper.snv |                      |                                                                                                                                                                            |                 |                       |                               | 3. Trim "Tumor_Sample_Barcode" to sample vial level;                        |
  | varscan2.snv      |                      |                                                                                                                                                                            |                 |                       |                               | 4. Rename columns as::                                                      |
  |                   |                      |                                                                                                                                                                            |                 |                       |                               |                                                                             |
  |                   |                      |                                                                                                                                                                            |                 |                       |                               |     {                                                                       |
  |                   |                      |                                                                                                                                                                            |                 |                       |                               |         'Hugo_Symbol': 'gene',                                              |
  |                   |                      |                                                                                                                                                                            |                 |                       |                               |         'Chromosome': 'chrom',                                              |
  |                   |                      |                                                                                                                                                                            |                 |                       |                               |         'Start_Position': 'start',                                          |
  |                   |                      |                                                                                                                                                                            |                 |                       |                               |         'End_Position': 'end',                                              |
  |                   |                      |                                                                                                                                                                            |                 |                       |                               |         'Reference_Allele': 'ref',                                          |
  |                   |                      |                                                                                                                                                                            |                 |                       |                               |         'Tumor_Seq_Allele2': 'alt',                                         |
  |                   |                      |                                                                                                                                                                            |                 |                       |                               |         'Tumor_Sample_Barcode': 'sampleid',                                 |
  |                   |                      |                                                                                                                                                                            |                 |                       |                               |         'HGVSp_Short': 'Amino_Acid_Change',                                 |
  |                   |                      |                                                                                                                                                                            |                 |                       |                               |         'Consequence': 'effect',                                            |
  |                   |                      |                                                                                                                                                                            |                 |                       |                               |         'FILTER': 'filter'                                                  |
  |                   |                      |                                                                                                                                                                            |                 |                       |                               |     }                                                                       |
  +-------------------+----------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------+-----------------------+-------------------------------+-----------------------------------------------------------------------------+

API
---

Check documentation for GDC module and Xena Dataset module `here <API.rst>`__.

