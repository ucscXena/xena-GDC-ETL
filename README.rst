xena-GDC-ETL
============

.. image:: https://travis-ci.org/ucscXena/xena-GDC-ETL.svg?branch=master
    :target: https://travis-ci.org/ucscXena/xena-GDC-ETL

Extract, transform and load `GDC data <https://portal.gdc.cancer.gov/>`__ onto `UCSC Xena <https://xenabrowser.net/>`__.

**Table of Contents**

- `Dependencies`_
- `Installation`_
- `Basic usage with command line tools`_
- `Advanced usage with XenaDataset and its subclasses`_
- `GDC ETL settings`_
- `Documentation`_

Dependencies
------------

Specific versions mentioned below have been tested. Eariler versions may still work but not guaranteed. 

1. Python 2.7, 3.5+

   This pipeline has been tested with python 2.7, 3.5, 3.6 and 3.7. It may also
   work with other python 3 versions since it was originally designed to be
   `single-source Python 2/3 compatible <https://docs.python.org/3/howto/pyporting.html#the-short-explanation>`_.

2. `Requests <http://python-requests.org>`_ v1.2.3
3. `Numpy <https://www.numpy.org/>`_ v1.15.0
4. `Pandas <https://pandas.pydata.org/>`_ v0.23.2
5. `Jinja2 <http://jinja.pocoo.org/>`_ v2.10.1: used for generating metadata JSON.
6. `lxml <https://lxml.de/>`_ v4.2.0: used for parsing TCGA phenotype data
7. `xlrd <https://xlrd.readthedocs.io/en/latest/>`_ v1.1.0: used for reading TARGET phenotype data

Installation
------------

- First clone the repository from GitHub by running
  ``git clone https://github.com/ucscXena/xena-GDC-ETL.git``. Now, ``cd`` into
  ``xena-GDC-ETL`` directory and install the package using pip: ``pip install .``

    If you are developing the package, you can use ``pip``'s edit mode for
    installation: ``pip install -e .``.
- You can also directly use ``pip`` to install the package. To get the latest
  code from GitHub master branch, run:
  ``pip install git+https://github.com/ucscXena/xena-GDC-ETL``.
  To get the latest stable version, run: 
  ``pip install xena-GDC-ETL``.
- Dependencies can be installed either before or after cloning this repository.
  You can install them by running ``pip install -r requirements.txt``.
- In general,

  - ``gdc.py`` contains functionalities related to GDC API, which requires no other modules in this package;
  - ``xena_dataset.py`` contains core functionalities for importing data from GDC to Xena and needs the ``gdc.py`` module in this package;
  - ``gdc2xena.py`` defines `a command line tool`__ which requires both ``gdc.py`` module and ``xena_dataset.py`` module in this package;

    __ gdc2xena_

  - ``gdc_check_new.py`` defines `a command line tool`__ which requiresthe ``gdc.py`` module in this package.

    __ gdc_check_new_

Basic usage with command line tools
-----------------------------------

.. _gdc2xena:

- **Import selected project(s) and selected type(s) of data from GDC to Xena**

  .. code:: bash

    xge etl [-h] [-r ROOT]
            [-p PROJECTS [PROJECTS ...] | -P NOT_PROJECTS [NOT_PROJECTS ...]]
            [-t DATATYPE [DATATYPE ...] | -T NOT_DATATYPE [NOT_DATATYPE ...]]
            [-D DELETE]

  This tool will perform a full import of dataset(s) into the root directory (specified by the ``-r`` option) with a default directory tree. In general, a full import has 3 steps: downloading raw data, making Xena matrix from raw data and generating matrix associated metadata. Data from each step will be saved to corresponding directories, whose structure is like this:

  ::

    root
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

  A dataset is defined by its project and data type. Projects of interest are provided through ``-p`` or ``-P`` option, and data types of interest are provided through the ``-t`` or ``-T`` option. Multiple inputs separated by whitespace are allowed and will be treated separately with all possible combinations. Valid projects should be valid project_id on GDC. Valid data types includes (without quotation marks): 'htseq_counts', 'htseq_fpkm', 'htseq_fpkm-uq', 'mirna', 'masked_cnv', 'muse_snv', 'mutect2_snv', 'somaticsniper_snv', 'varscan2_snv', 'phenotype', and 'survival'. Upper case options (``-P`` or ``-T``) are mutually exclusive with corresponding lower case options, and they are used to define datasets of interest by excluding selections from either all projects on GDC or all supported data types. For example, the following command line imports 3 types of RNA-seq data for all but FM-AD projects from GDC to /home/user/xena_root:

  .. code:: bash

    mkdir -p /home/user/xena_root
    xge etl -P FM-AD -t htseq_counts htseq_fpkm htseq_fpkm-uq

  Notes:

  1. Root directory must be existing
  2. Please check the next section for `advanced usage with XenaDataset and its subclasses`_, if you want to customize the importing process with selected (rather than all possible) combinations of your input projects and data types or selected (rather than all 3) importing step(s).

- **Generate metadata of a xena matrix**

  .. code:: bash

    xge metadata --project TCGA-BRCA --datatype htseq_counts --matrix path/to/matrix.tsv --release 10

  This tool generates metadata for a xena matrix. For the shown example, metadata
  is generated for the matrix ``matrix.tsv`` for release ``10``, project
  ``TCGA-BRCA`` and datatype ``htsep_counts``. Note that, metadata JSON file is
  saved at the same directory as the ``matrix.tsv`` file.

.. _gdc_check_new:

- **Check against a list of updated files for affected dataset(s)**

  .. code:: bash

    xge gdc_check_new [-h] URL

  This tool takes in a file (either a URL or a local file readable by ``pandas.read_csv``) of table and read one of its columns named as "New File UUID". It then checks all file UUIDs in this table on GDC and summarize all their associated project(s), data type(s) and analysis workflow type(s). Such tables are usually provided in GDC's data release note. With the summarized info, you can design specific imports to just update datasets which are updated on GDC. For example, the following command:

  .. code:: bash

    xge gdc_check_new https://docs.gdc.cancer.gov/Data/Release_Notes/DR9.0_files_swap.txt.gz

  should give you:

  .. code:: bash

    analysis.workflow_type    cases.project.project_id    data_type
    HTSeq - FPKM              TARGET-NBL                  Gene Expression Quantification
    HTSeq - FPKM-UQ           TARGET-NBL                  Gene Expression Quantification
    HTSeq - Counts            TARGET-NBL                  Gene Expression Quantification

.. _version:

- **Shows the current version of xena_gdc_etl**

  .. code:: bash

    xge --version

.. _xena-eql:

- **Check equality of two xena matrices**

  .. code:: bash

    xge xena-eql path/to/matrix1.tsv path/to/matrix2.tsv

  This tool takes path to two xena matrices and output if they are equal or not.

.. _merge-xena:

- **Merge xena matrices**

  .. code:: bash

    xge merge-xena -f path/to/matrix1.tsv path/to/matrix2.tsv -t htseq_counts -o path/to/output -n new_name.tsv -c TCGA-BRCA

  This tool merges xena matrices and outputs the merged matrix. For the given
  example the tool will merge ``matrix1.tsv`` and ``matrix2.tsv`` matrices and
  store the merged matrix in ``path/to/output`` directory with the name
  ``new_name.tsv``. Note that, had the argument ``-n`` not been
  specified, the merged matrix would have been saved as
  ``TCGA-BRCA.htseq_counts.tsv``.


Advanced usage with XenaDataset and its subclasses
--------------------------------------------------

- **The** ``XenaDataset`` **class**

  Though this is not an abstract class, it is designed as a generalized class representing one Xena dataset and its importing process. For doing an import of GDC data, use its subclasses_, which have preloaded with some default settings, might be simpler.
  
  A Xena dataset is defined by its study project (cohort) and the type of data in this dataset. A typical importing process has the following 3 steps:
  
  1. Download raw data from the source.
  
    The ``download_map`` property defines a dict of raw data to be downloaded, with the key being the URL and the value being the path, including the filename, for saving corresponding downloaded file. The ``download`` method will read the ``download_map`` and perform the downloading, creating non-existing directories as needed. After downloading all files, a list of paths for downloaded files will be recorded in the ``raw_data_list`` property. The ``download`` method needs only a valid ``download_map``. It will return the object itself, therefore can be chained with ``transform``.
  
  2. Transform raw data into valid Xena matrix.
  
    One assumption for data transformation is that there might be multiple raw data (in the ``raw_data_list``) supporting the single Xena matrix in a dataset. Therefore, the ``transform`` method will first merge data and then process merged matrix as needed. It will open the file one by one accordingly (by extension), and read the file object and transform its data with a function defined by ``read_raw``. The list of transformed single data will be merged and processed by a function defined by ``raws2matrix``, which gives the finalized Xena matrix. The ``transform`` method requires a valid list of raw data, besides ``read_raw`` and ``raws2matrix``. A valid list of raw data can be either explicitly defined by ``raw_data_list`` or can be derived from ``raw_data_dir`` with all files under ``raw_data_dir`` being treated as raw data. It will return the object itself, therefore can be chained with ``metadata``.
  
  3. Generate metadata for the new Xena matrix.
  
    Metadata for Xena matrix is a JSON file rendered by the ``metadata`` method with ``metadata_vars`` (dict) through Jinja2 from ``metadata_template``. This JSON file will be saved under the same directory as the matrix, with a filename being the matrix name plus the '.json' postfix. The ``metadata`` method requires an existing file of Xena matrix.
  
  .. _directory related settings:
  
  ``root_dir`` is both an optional instantiation arguments and a property. By default, it points to the current working directory. It is worth mentioning that the default directory structure mentioned above is implemented in the class. However, you are free to changed the setting with the following properties:
  
  - Pass an argument for ``root_dir`` during instantiation or set the ``root_dir`` property explicitly after instantiation.
  - Downloaded raw data will be saved under ``raw_data_dir``.
  - Newly transformed Xena matrix will be saved as ``matrix`` under ``matrix_dir``. The directory path in ``matrix`` has the priority over ``matrix_dir``. By default, Xena matrix will be saved under the "matrix_dir" as "<projects>.<xena_dtype>.tsv".
  - Metadata will always have the specific pattern of name and be together with ``matrix`` (i.e. no way to change this behavior).

.. _subclasses:

- **Build GDC importing pipelines with** ``GDCOmicset``, ``GDCPhenoset`` **or** ``GDCSurvivalset`` **classes**

  ``GDCOmicset``, ``GDCPhenoset`` and ``GDCSurvivalset`` are subclasses of ``XenaDataset`` and are preloaded with settings for importing GDC genomic data, TCGA phenotype data on GDC, TARGET phenotype data on GDC and GDC's survival data respecitively. These settings can be customized by setting corresponding properties described below. For more details, please check the `next section <#gdc-etl-settings>`__ and the `documentation <https://github.com/ucscXena/xena-GDC-ETL/blob/master/docs/API.rst>`_.
  
  The script for ``gdc2xena.py`` command line is a good example for basic usage of these classes. Similar to ``XenaDataset``, a GDC dataset is defined by ``projects``, which is one or a list of valid GDC "project_id". For ``GDCOmicset``, a dataset should also be defined with one of the supported ``xena_dtype`` (find out with the class method ``GDCOmicset.get_supported_dtype()``). The ``xena_dtype`` is critical for a ``GDCOmicset`` object selecting correct default settings. For ``GDCPhenoset`` and ``GDCSurvivalset``, data type are self-explanatory and cannot be changed. Therefore, you can instantiate these classes like this:
  
  .. code:: python
  
    from xena_dataset import GDCOmicset, GDCPhenoset, GDCSurvivalset, GDCAPIPhenoset
    
    gdc_omic_cohort = GDCOmicset('TCGA-BRCA', 'htsep_counts')
    
    # Won't check if the ID is of TCGA program or not.
    tcga_pheno_cohort = GDCPhenoset('TCGA-BRCA')
    
    # Won't check if the ID is of TARGET program or not.
    target_pheno_cohort = GDCPhenoset('TARGET-NBL')
    
    gdc_survival_cohort = GDCSurvivalset('TCGA-BRCA')

    gdc_api_pheno_cohort = GDCAPIPhenoset('CPTAC-3')
  
  With such a dataset object, it is fine to call ``download``, ``transform`` and/or ``metadata`` method(s). These methods will use preloaded settings and save files under ``root_dir`` accordingly. You are free to call/chain some but not all 3 methods; just keep in mind the pre-requisites for each method and set related properties properly. Aside from `directory related settings`_ described above, you can change some default importing settings through the following properties.
  
  .. _Customize GDCOmicset:
  
  - **Customize** ``GDCOmicset``
  
  |
  
  +-------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
  |    Attributes     |                                                                         Usage                                                                          |                                                                                                        Type and Format\ :sup:`1`                                                                                                         |                                                                                                               Default settings                                                                                                                 |
  +===================+========================================================================================================================================================+==========================================================================================================================================================================================================================================+================================================================================================================================================================================================================================================+
  | gdc_filter        | Used for deriving default ``download_map`` as the GDC search filters.                                                                                  | ``dict``: the key is 1 GDC available file field and the value is either a string or a list, meaning the value of the file field matches a string or number in (a list)                                                                   | Check `GDC download settings`_ for details.                                                                                                                                                                                                    |
  +-------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
  | gdc_prefix        | Used for deriving default ``download_map`` as the GDC search fields.                                                                                   | ``str``: 1 GDC available file field whose value will be the prefix of the filename of corresponding downloaded file.                                                                                                                     | Check `GDC download settings`_ for details.                                                                                                                                                                                                    |
  +-------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
  | download_map      | Used by the ``download`` method for downloading GDC raw data supporting this dataset.                                                                  | ``dict``: the key is download URL and the value is the desired path for saving the downloaded file.                                                                                                                                      | Download URLs are in the pattern of "https://api.gdc.cancer.gov/data/<FILE UUID>", and paths are in the pattern of "<``raw_data_dir``>/<value of gdc_prefix>.<GDC file UUID>.<file extension>".                                                |
  +-------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
  | read_raw          | Used by the ``transform`` method when reading a single GDC raw data.                                                                                   | ``callable``: takes only 1 file object as its argument and returns an arbitrary result which will be put in a list and passed on to ``raws2matrix``.                                                                                     | Check `GDC genomic transform settings`_ for details                                                                                                                                                                                            |
  +-------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
  | raws2matrix       | Used by the ``transform`` method and responsible for both merging multiple GDC raw data into one Xena matrix and processing new Xena matrix as needed. | ``callable``: takes only 1 list of ``read_raw`` returns as its argument and returns an object (usually a pandas DataFrame) which has a ``to_csv`` method for saving as a file.                                                           | Check `GDC genomic transform settings`_ for details                                                                                                                                                                                            |
  +-------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
  | metadata_template | Used by the ``metadata`` method for rendering metadata by Jinja2.                                                                                      | ``jinja2.environment.Template`` or ``str``: a ``jinja2.environment.Template`` used directly by Jinja2; if it's a string, it is a path to the template file which will be silently read and converted to ``jinja2.environment.Template``. | `Resources <https://github.com/ucscXena/xena-GDC-ETL/blob/master/Resources>`_                                                                                                                                                                  |
  +-------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
  | metadata_vars     | Used by the ``metadata`` method for rendering metadata by Jinja2.                                                                                      | ``dict``: used directly by Jinja2 which should match variables in ``metadata_template``.                                                                                                                                                 | ::                                                                                                                                                                                                                                             |
  |                   |                                                                                                                                                        |                                                                                                                                                                                                                                          |                                                                                                                                                                                                                                                |
  |                   |                                                                                                                                                        |                                                                                                                                                                                                                                          |   {                                                                                                                                                                                                                                            |
  |                   |                                                                                                                                                        |                                                                                                                                                                                                                                          |       'project_id': <``projects``>,                                                                                                                                                                                                            |
  |                   |                                                                                                                                                        |                                                                                                                                                                                                                                          |       'date': <the time of last modification of ``matrix``>,                                                                                                                                                                                   |
  |                   |                                                                                                                                                        |                                                                                                                                                                                                                                          |       'gdc_release': <``gdc_release``>,                                                                                                                                                                                                        |
  |                   |                                                                                                                                                        |                                                                                                                                                                                                                                          |       'xena_cohort': <Xena specific cohort name for TCGA data or GDC project_id for TARGET data, with (for both) "GDC " prefix>                                                                                                                |
  |                   |                                                                                                                                                        |                                                                                                                                                                                                                                          |   }                                                                                                                                                                                                                                            |
  |                   |                                                                                                                                                        |                                                                                                                                                                                                                                          |                                                                                                                                                                                                                                                |
  |                   |                                                                                                                                                        |                                                                                                                                                                                                                                          | \* The first element of the "url" field in metadata will be "gdc_release" URL, and the second will be specific URL for raw data file if there is only 1 raw data file for this dataset; or it will be just "https://api.gdc.cancer.gov/data/". |
  +-------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
  | gdc_release       | Used by the ``metadata`` method for rendering metadata, showing the GDC data release of this dataset.                                                  | ``str``: an URL pointing to corresponding GDC Data Release Note.                                                                                                                                                                         | Current data release version when the ``gdc_release`` is being used/called, queried through "https://api.gdc.cancer.gov/status".                                                                                                               |
  +-------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
  
  \1. GDC API Available File Fields: https://docs.gdc.cancer.gov/API/Users_Guide/Appendix_A_Available_Fields/#file-fields
  
  - **Customize** ``GDCPhenoset`` **for TCGA projects**
  
    TCGA phenotype data for Xena includes both clinical data and biospecimen data, as `detailed below <#transform-phenotype>`_. Downloading and transformation of clinical data and biospecimen data are in fact delegated by two independent ``GDCOmicset`` object respecitively. Corresponding subdatasets can be accessed through ``clin_dataset`` and ``bio_dataset`` attributes and hence can be customized as mentioned above. Because of such complexity of TCGA phenotype data, the ``download`` and ``transform`` methods are coded specifically and overrode corresponding methods of the base class, ``XenaDataset``. Customization for downloading and matrix transformation is very limited and should be done in the following steps:
    
    1. Instantiate a ``GDCPhenoset``;
    2. Instantiate and customize one or two ``GDCOmicset`` objects for clinical data and/or biospecimen data as needed;
    3. Assign customized ``GDCOmicset`` objects to corresponding attributes, ``clin_dataset`` and ``bio_dataset``;
    4. Call desired method(s) (``download`` and/or ``transform``).
    
    - Customize ``download`` step
    
      This step can be customized only through customized ``clin_dataset`` and ``bio_dataset``, since the whole downloading process is delegated by these two GDCOmicset objects.
      
    - Customize ``transform`` step
    
      The first part of ``transform`` is delegated by ``transform`` methods of ``clin_dataset`` and ``bio_dataset``. Therefore, the only way to customized this process is to customize ``clin_dataset`` and ``bio_dataset``. How the two matrices are then merged into one Xena phenotype matrix is hard coded and cannot be customized. It is worth noting that if you want to call ``transfrom`` but skip the downloading step, you will need to define ``clin_dataset`` and ``bio_dataset`` before calling ``transform``.
      
    - Customize ``metadata`` step
    
      Different from ``download`` and ``transform``, there is no special settings for the ``metadata`` method of ``GDCPhenoset``. Therefore, similar to that of ``GDCOmicset``, this step can be customized through ``metadata_template``, ``metadata_vars`` and ``gdc_release`` properties. And to call just the ``metadata`` method, an existing ``matrix`` is enough.

  - **Customize** ``GDCPhenoset`` for **TARGET projects**

    TARGET phenotype data for Xena contains only the clinical data (no biospecimen data), as `detailed below <#transform-phenotype>`_. The importing process is quite similar to that of a ``GDCOmicset``. You can customize ``TARGETPhenoset`` with ``download_map``, ``read_raw``, ``raws2matrix``, ``metadata_template``, ``metadata_vars`` and ``gdc_release`` in the same way as that of `GDCOmicset <#customize-gdcomicset>`_.

  - **Customize** ``GDCSurvivalset``
  
    GDC data supporting Xena survival matrix does not come any GDC files. It comes from the "analysis/survival" endpoint of GDC API. Therefore, the ``download`` and ``transform`` methods are re-designed, overriding those of the base class, ``XenaDataset``. Aside from redefining ``download`` and ``transform`` methods, there is no simple way to customize ``download`` and ``transform`` steps. You can still call ``transform`` without ``download`` by just defining a valid list of raw data with ``raw_data_list`` or ``raw_data_dir``. However, only this first file in the list will be read and used.
    
    Different from ``download`` and ``transform``, there is no special settings for the ``metadata`` method of ``GDCSurvivalset``. Therefore, similar to that of ``GDCOmicset``, this step can be customized through ``metadata_template``, ``metadata_vars`` and ``gdc_release`` properties. To call just the ``metadata`` method, an existing ``matrix`` is enough.

  - **Customize** ``GDCAPIPhenoset``

    The data for this class comes from GDC API only. Therefore, the ``download`` 
    and ``transform`` methods are re-designed, overriding those of the base class, 
    ``XenaDataset``. Aside from redefining ``download`` and ``transform`` methods, 
    there is no simple way to customize ``download`` and ``transform`` steps.

    Different from ``download`` and ``transform``, there is no special settings 
    for the ``metadata`` method of ``GDCAPIPhenoset``. Therefore, similar to that 
    of ``GDCOmicset``, this step can be customized through ``metadata_template``, 
    ``metadata_vars`` and ``gdc_release`` properties. To call just the ``metadata`` 
    method, an existing ``matrix`` is enough.


GDC ETL settings
-------------------

.. _GDC download settings:

- **Settings for downloading/getting raw data (files) from GDC**

  +-------------------+-------------------+-----------------------------------+-----------------------------------------------+--------------------------+------------------------------------------------------+
  |                   |                   |                               GDC data filter                                     |                          |                                                      |
  +    xena_dtype     + GDC API endpoint  +-----------------------------------+-----------------------------------------------+ File count/Level         + GDC file field for filename prefix                   +
  |                   |                   | data_type                         | analysis.workflow_type                        |                          |                                                      |
  +===================+===================+===================================+===============================================+==========================+======================================================+
  | htseq_counts      | data              | Gene Expression Quantification    | HTSeq - Counts                                | 1/Sample vial            | cases.samples.submitter_id                           |
  +-------------------+-------------------+-----------------------------------+-----------------------------------------------+--------------------------+------------------------------------------------------+
  | htseq_fpkm        | data              | Gene Expression Quantification    | HTSeq - FPKM                                  | 1/Sample vial            | cases.samples.submitter_id                           |
  +-------------------+-------------------+-----------------------------------+-----------------------------------------------+--------------------------+------------------------------------------------------+
  | htseq_fpkm-uq     | data              | Gene Expression Quantification    | HTSeq - FPKM-UQ                               | 1/Sample vial            | cases.samples.submitter_id                           |
  +-------------------+-------------------+-----------------------------------+-----------------------------------------------+--------------------------+------------------------------------------------------+
  | mirna             | data              | miRNA Expression Quantification   | BCGSC miRNA Profiling                         | 1/Sample vial            | cases.samples.submitter_id                           |
  +-------------------+-------------------+-----------------------------------+-----------------------------------------------+--------------------------+------------------------------------------------------+
  | mirna_isoform     | data              | Isoform Expression Quantification | BCGSC miRNA Profiling                         | 1/Sample vial            | cases.samples.submitter_id                           |
  +-------------------+-------------------+-----------------------------------+-----------------------------------------------+--------------------------+------------------------------------------------------+
  | cnv               | data              | Copy Number Segment               | DNAcopy                                       | 1/Sample vial            | cases.samples.submitter_id                           |
  +-------------------+-------------------+-----------------------------------+-----------------------------------------------+--------------------------+------------------------------------------------------+
  | masked_cnv        | data              | Masked Copy Number Segment        | DNAcopy                                       | 1/Sample vial            | cases.samples.submitter_id                           |
  +-------------------+-------------------+-----------------------------------+-----------------------------------------------+--------------------------+------------------------------------------------------+
  | gistic            | data              | Gene Level Copy Number Scores     | GISTIC - Copy Number Score                    | 1/Project                | submitter_id                                         |
  +-------------------+-------------------+-----------------------------------+-----------------------------------------------+--------------------------+------------------------------------------------------+
  | star_counts       | data              | STARCounts                        | STAR - Counts                                 | 1/Sample vial            | cases.samples.submitter_id                           |
  +-------------------+-------------------+-----------------------------------+-----------------------------------------------+--------------------------+------------------------------------------------------+
  | muse_snv          | data              | Masked Somatic Mutation           | MuSE Variant Aggregation and Masking          | 1/Project                | submitter_id                                         |
  +-------------------+-------------------+-----------------------------------+-----------------------------------------------+--------------------------+------------------------------------------------------+
  | mutect2_snv       | data              | Masked Somatic Mutation           | MuTect2 Variant Aggregation and Masking       | 1/Project                | submitter_id                                         |
  +-------------------+-------------------+-----------------------------------+-----------------------------------------------+--------------------------+------------------------------------------------------+
  | somaticsniper_snv | data              | Masked Somatic Mutation           | SomaticSniper Variant Aggregation and Masking | 1/Project                | submitter_id                                         |
  +-------------------+-------------------+-----------------------------------+-----------------------------------------------+--------------------------+------------------------------------------------------+
  | varscan2_snv      | data              | Masked Somatic Mutation           | VarScan2 Variant Aggregation and Masking      | 1/Project                | submitter_id                                         |
  +-------------------+-------------------+-----------------------------------+-----------------------------------------------+--------------------------+------------------------------------------------------+
  | clinical          | data              | Clinical Supplement               | N/A                                           | 0 or 1/Case              | cases.submitter_id                                   |
  +-------------------+-------------------+-----------------------------------+-----------------------------------------------+--------------------------+------------------------------------------------------+
  | biospecimen       | data              | Biospecimen Supplement            | N/A                                           | 1/Case                   | cases.submitter_id                                   |
  +-------------------+-------------------+-----------------------------------+-----------------------------------------------+--------------------------+------------------------------------------------------+
  | survival          | analysis/survival | N/A (filtered by just the "project.project_id")                                   | 1 Record/Case (Non-file) | N/A (filename will be "<projects>.GDC_survival.tsv") |
  +-------------------+-------------------+-----------------------------------+-----------------------------------------------+--------------------------+------------------------------------------------------+

.. _GDC genomic transform settings:

- **Settings for transform "Omic" data into Xena matrix**

  +-------------------+----------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------+-----------------------+-------------------------------+-----------------------------------------------------------------------------------------------------------------------+
  |  xena_dtype       | Raw data has header? | Select columns (in order)                                                                                                                                                  | Row index       | Skip rows start with? | Merge into matrix as          | Process matrix                                                                                                        |
  +===================+======================+============================================================================================================================================================================+=================+=======================+===============================+=======================================================================================================================+
  | htseq_counts      | No                   | 1, 2                                                                                                                                                                       | Ensembl_ID      | _                     | 1 new column based on index   | 1. Average if there are multiple data from the same sample vial;                                                      |
  |                   |                      | [Ensembl_ID, Counts]                                                                                                                                                       |                 |                       |                               | 2. log2(counts + 1)                                                                                                   |
  +-------------------+----------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------+-----------------------+-------------------------------+-----------------------------------------------------------------------------------------------------------------------+
  | htseq_fpkm        | No                   | 1, 2                                                                                                                                                                       | Ensembl_ID      | _                     | 1 new column based on index   | 1. Average if there are multiple data from the same sample vial;                                                      |
  |                   |                      | [Ensembl_ID, Counts]                                                                                                                                                       |                 |                       |                               | 2. log2(counts + 1)                                                                                                   |
  +-------------------+----------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------+-----------------------+-------------------------------+-----------------------------------------------------------------------------------------------------------------------+
  | htseq_fpkm-uq     | No                   | 1, 2                                                                                                                                                                       | Ensembl_ID      | _                     | 1 new column based on index   | 1. Average if there are multiple data from the same sample vial;                                                      |
  |                   |                      | [Ensembl_ID, Counts]                                                                                                                                                       |                 |                       |                               | 2. log2(counts + 1)                                                                                                   |
  +-------------------+----------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------+-----------------------+-------------------------------+-----------------------------------------------------------------------------------------------------------------------+
  | mirna             | Yes                  | 1, 3                                                                                                                                                                       | miRNA_ID        | N/A                   | 1 new column based on index   | 1. Average if there are multiple data from the same sample vial;                                                      |
  |                   |                      | [miRNA_ID, RPM]                                                                                                                                                            |                 |                       |                               | 2. log2(counts + 1)                                                                                                   |
  +-------------------+----------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------+-----------------------+-------------------------------+-----------------------------------------------------------------------------------------------------------------------+
  | mirna_isoform     | Yes                  | 2, 4                                                                                                                                                                       | isoform_coords  | N/A                   | 1 new column based on index   | 1. Average if there are multiple data from the same sample vial;                                                      |
  |                   |                      | [isoform_coords, RPM]                                                                                                                                                      |                 |                       |                               | 2. log2(counts + 1)                                                                                                   |
  +-------------------+----------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------+-----------------------+-------------------------------+-----------------------------------------------------------------------------------------------------------------------+
  | cnv               | Yes                  | 2, 3, 4, 6                                                                                                                                                                 | sample          | N/A                   | New rows based on column name | 1. Rename columns as::                                                                                                |
  |                   |                      | [Chromosome, Start, End, Segment_Mean]                                                                                                                                     |                 |                       |                               |                                                                                                                       |
  |                   |                      |                                                                                                                                                                            |                 |                       |                               |     {                                                                                                                 |
  |                   |                      |                                                                                                                                                                            |                 |                       |                               |         'Chromosome': 'Chrom',                                                                                        |
  |                   |                      |                                                                                                                                                                            |                 |                       |                               |         'Segment_Mean': 'value'                                                                                       |
  |                   |                      |                                                                                                                                                                            |                 |                       |                               |     }                                                                                                                 |
  +-------------------+----------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------+-----------------------+-------------------------------+-----------------------------------------------------------------------------------------------------------------------+
  | masked_cnv        | Yes                  | 1, 2, 3, 5                                                                                                                                                                 | sample          | N/A                   | New rows based on column name | 1. Rename columns as::                                                                                                |
  |                   |                      | [Chromosome, Start, End, Segment_Mean]                                                                                                                                     |                 |                       |                               |                                                                                                                       |
  |                   |                      |                                                                                                                                                                            |                 |                       |                               |     {                                                                                                                 |
  |                   |                      |                                                                                                                                                                            |                 |                       |                               |         'Chromosome': 'Chrom',                                                                                        |
  |                   |                      |                                                                                                                                                                            |                 |                       |                               |         'Segment_Mean': 'value'                                                                                       |
  |                   |                      |                                                                                                                                                                            |                 |                       |                               |     }                                                                                                                 |
  +-------------------+----------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------+-----------------------+-------------------------------+-----------------------------------------------------------------------------------------------------------------------+
  | gistic            | Yes                  | 1                                                                                                                                                                          | Ensembl_ID      | _                     | N/A                           | 1. Drop "Gene ID" and "Cytoband" column;                                                                              |
  |                   |                      | [Ensembl_ID]                                                                                                                                                               |                 |                       |                               | 2. Map "samples.portions.analytes.aliquots.aliquot_id" into "samples.submitter_id" using GDC API and use it as index. |
  +-------------------+----------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------+-----------------------+-------------------------------+-----------------------------------------------------------------------------------------------------------------------+
  | muse_snv          | Yes                  | 13, 37, 5, 6, 7, 40, 42, 52, 1, 11, 16, 111                                                                                                                                | N/A             | #                     | N/A                           | 1. Calculate variant allele frequency (dna_vaf) by "t_alt_count"/"t_depth";                                           |
  | mutect2_snv       |                      | [Tumor_Seq_Allele2, HGVSp_Short, Chromosome, Start_Position, End_Position, t_depth, t_alt_count, Consequence, Hugo_Symbol, Reference_Allele, Tumor_Sample_Barcode, FILTER] |                 |                       |                               | 2. Delete "t_alt_count" and "t_depth" columns;                                                                        |
  | somaticsniper_snv |                      |                                                                                                                                                                            |                 |                       |                               | 3. Trim "Tumor_Sample_Barcode" to sample vial level;                                                                  |
  | varscan2_snv      |                      |                                                                                                                                                                            |                 |                       |                               | 4. Rename columns as::                                                                                                |
  |                   |                      |                                                                                                                                                                            |                 |                       |                               |                                                                                                                       |
  |                   |                      |                                                                                                                                                                            |                 |                       |                               |     {                                                                                                                 |
  |                   |                      |                                                                                                                                                                            |                 |                       |                               |         'Hugo_Symbol': 'gene',                                                                                        |
  |                   |                      |                                                                                                                                                                            |                 |                       |                               |         'Chromosome': 'chrom',                                                                                        |
  |                   |                      |                                                                                                                                                                            |                 |                       |                               |         'Start_Position': 'start',                                                                                    |
  |                   |                      |                                                                                                                                                                            |                 |                       |                               |         'End_Position': 'end',                                                                                        |
  |                   |                      |                                                                                                                                                                            |                 |                       |                               |         'Reference_Allele': 'ref',                                                                                    |
  |                   |                      |                                                                                                                                                                            |                 |                       |                               |         'Tumor_Seq_Allele2': 'alt',                                                                                   |
  |                   |                      |                                                                                                                                                                            |                 |                       |                               |         'Tumor_Sample_Barcode': 'sampleid',                                                                           |
  |                   |                      |                                                                                                                                                                            |                 |                       |                               |         'HGVSp_Short': 'Amino_Acid_Change',                                                                           |
  |                   |                      |                                                                                                                                                                            |                 |                       |                               |         'Consequence': 'effect',                                                                                      |
  |                   |                      |                                                                                                                                                                            |                 |                       |                               |         'FILTER': 'filter'                                                                                            |
  |                   |                      |                                                                                                                                                                            |                 |                       |                               |     }                                                                                                                 |
  +-------------------+----------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------+-----------------------+-------------------------------+-----------------------------------------------------------------------------------------------------------------------+
  | star_counts       | Yes                  | 1, 2                                                                                                                                                                       | Ensembl_ID      | _                     | 1 new column based on index   | 1. Average if there are multiple data from the same sample vial;                                                      |
  |                   |                      | [Ensembl_ID, Counts]                                                                                                                                                       |                 |                       |                               | 2. log2(counts + 1)                                                                                                   |
  +-------------------+----------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------+-----------------------+-------------------------------+-----------------------------------------------------------------------------------------------------------------------+

.. _transform phenotype:

- **Settings for transform phenotype data into Xena matrix**

  +-------------+------------------------------------------------+-----------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
  | GDC program |                  GDC raw data                  | Raw data format | Single data file transformation                                                                                                                                                                                                                                                                                                                                                                                                                                                         | Merge and matrix processing                                                                                                                                                                                            |
  +=============+================================================+=================+=========================================================================================================================================================================================================================================================================================================================================================================================================================================================================================+========================================================================================================================================================================================================================+
  | TCGA        | Clinical Supplement and Biospecimen Supplement | BCR XML         | For clincial data, info will be extracted and organized into a per patient based pandas DataFrame. It will have a column named "bcr_patient_barcode" which will be used to join with biospecimen matrix later on.                                                                                                                                                                                                                                                                       | 1. Multiple clinical data are concatenated directly by row with all empty columns removed.                                                                                                                             |
  |             |                                                |                 |                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         | 2. Multiple biospecimen data are concatenated directly by row with all empty columns removed.                                                                                                                          |
  |             |                                                |                 | The XML scheme are quite different for different projects. Therefore, to get as much info as possible while still keeping things clear, texts, if any, from all elements that have non-element children are extracted first. After such a "dirty" extraction, two clean ups will be done:                                                                                                                                                                                               | 3. Merged clinical matrix and merged biospecimen matrix are further merged on "bcr_patient_barcode". For conflict/overlapping columns, non-empty value from the clinical data has the priority.                        |
  |             |                                                |                 |                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |                                                                                                                                                                                                                        |
  |             |                                                |                 | 1. For "race" info, it will be converted into a comma separated list of races, in case there are more than one entry in <clin_shared:race_list> in the clinical XML.                                                                                                                                                                                                                                                                                                                    |                                                                                                                                                                                                                        |
  |             |                                                |                 | 2. When there is one or more follow ups, the most recent follow up will be find out. All info in the most recent follow up will be used to replace/add to previously extracted matrix.                                                                                                                                                                                                                                                                                                  |                                                                                                                                                                                                                        |
  |             |                                                |                 |                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |                                                                                                                                                                                                                        |
  |             |                                                |                 | For biospecimen data, there is one coherent XML scheme for all TCGA projects. There are two parts to be considered for biospecimen data: per sample/sample specific data and patient data (which is common for all samples). Info from both parts will be extracted and finally organized into a per sample based matrix, having a column named "bcr_patient_barcode", which will be used to join with clinical matrxi later on. In general, info extraction has the following 3 steps: |                                                                                                                                                                                                                        |
  |             |                                                |                 |                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |                                                                                                                                                                                                                        |
  |             |                                                |                 | 1. Common patient data will be extracted first, including texts from direct children of <admin:admin> and <bio:patient>. A new field of "primary_diagnosis" will be added by mapping "disease_code" to `TCGA study name <https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations>`_.                                                                                                                                                                      |                                                                                                                                                                                                                        |
  |             |                                                |                 | 2. Samples from <bio:patient/bio:samples> will be processed and have comman patient data attached one by one. Non-empty texts from direct children of sample will be extracted, i.e. details from nodes like <bio:portions> will be dropped. Samples having `type code 10 <https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes>`_ are dropped.                                                                                                               |                                                                                                                                                                                                                        |
  |             |                                                |                 | 3. A column of "bcr_patient_barcode" from <bio:patient/shared:bcr_patient_barcode> will be added to the final biospecimen matrix (same for the whole table).                                                                                                                                                                                                                                                                                                                            |                                                                                                                                                                                                                        |
  +-------------+------------------------------------------------+-----------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
  | TARGET      | Clinical Supplement only                       | XLSX            | The excel file is converted to a pandas DataFrame.                                                                                                                                                                                                                                                                                                                                                                                                                                      | 1. Multiple DataFrames will be concatenated directly by row, and arriage return and line feed are replaced by a single space.                                                                                          |
  |             |                                                |                 |                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         | 2. Clinical data is per case(patient) based, while Xena phenotype matrix is per sample based. All related samples for each case/patient will be identified and phenotype data will be mapped to corresponding samples. |
  +-------------+------------------------------------------------+-----------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

- **Settings for transform survival data into Xena matrix**

  GDC survival data is returned as JSON from GDC API. During the downloading process, it can and will be converted directly to pandas DataFrame and saved as tab delimited table. During transformation, columns in "primary" Xena survival matrix can be mapped directly (without further processing/calculation) from the raw table like this:

  +---------------------+-------------------+
  | Primary Xena column | GDC source column |
  +=====================+===================+
  | OS.time             | time              |
  +---------------------+-------------------+
  | OS                  | censored          |
  +---------------------+-------------------+
  | _PATIENT            | submitter_id      |
  +---------------------+-------------------+

  GDC survival data is per case(patient) based and so is "primary" Xena survival matrix, while Xena survival matrix is per sample based. All related samples for each case/patient will be identified and survival data will be mapped to corresponding samples.

- **CPTAC-3 Cohort**

  CPTAC-3 data consists of RNAseq data (as discussed in ``GDCOmicset``) and
  clinical data from the API. The cases and expand for clinical data are
  defined in the ``constants.py`` file.

Documentation
-------------

Check documentation for GDC module and Xena Dataset module `here <https://github.com/ucscXena/xena-GDC-ETL/blob/master/docs/API.rst>`_.
