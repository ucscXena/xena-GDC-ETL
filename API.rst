
xena-GDC-ETL
************


API
===

* **gdc module**

  This module provides basic and minimum necessary functions for carrying out
  data query and download for Xena GDC ETL pipelines.

  **gdc.download(uuids, download_dir='.', chunk_size=4096)**

     Download GDC’s open access data according to UUID(s).

     Args:
        uuids (str, list or dict): A single UUID (str), a list of UUIDs (list)
           or a dict whose keys are UUIDs for target file(s). If “uuids” is str
           or list, downloaded file(s) will be renamed to “UUID.extension” where
           “extension” is extracted by “get_ext()” from the original filename.
           Renamed file(s) will be saved at “download_dir”. If “uuids” is a
           dict, the argument “download_dir” will be ignored; values of dict
           will be paths for saving corresponding downloaded files.

        download_dir (str, optional): The directory for saving downloaded
           file(s) when “uuids” is str or list. It will be ignored if “uuids” is
           dict. Defaults to “.”.

        chunk_size (int, optional): The chunk size is the number of bytes it
           should read into memory, when the response is got with “stream=True”.
           Check the documentation of “requests” module for details. Defaults to
           4096.

     Returns:
        list: a list of paths for downloaded files.

  **gdc.get_clinical_samples(projects=None)**

     Get info for all samples of ``projects`` and clinical info for all cases of
     ``projects`` through GDC API.

     Args:
        projects (list or str): one (str) or a list of GDC “project_id”(s),
           whose info will be returned. If None, projects will not be filtered,
           i.e. info for all GDC projects will be returned. Defaults to None.

     Returns:
        pandas.core.frame.DataFrame: A DataFrame organized by samples, having
        info for all samples of ``projects``, as well as corresponding clinical
        info.

  **gdc.get_ext(file_name)**

     Get all extensions supported by this module in the file name.

     Supported extensions are defined in the constant “_SUPPORTED_FILE_TYPES”.
     Multiple extensions will be separated by “.”.

     Args:
        file_name (str): The filename will be split by “.” and checked from
           left to right. Extensions will be kept starting from the first (and
           left most) supported extension.

     Returns:
        str: A string of extensions joint by “.”.

  **gdc.get_project_info(projects=None)**

     Get info for project(s) of interest through GDC API.

     Args:
        projects (list or str): one (str) or a list of GDC “project_id”(s),
           whose info will be returned. If None, projects will not be filtered,
           i.e. info for all GDC projects will be returned. Defaults to None.

     Returns:
        pandas.core.frame.DataFrame: A DataFrame of project info including
        “project ID”, “project name”, “primary site” and “program name”.

  **gdc.mkdir_p(dir_name)**

     Make the directory as needed: no error if existing.

     Args:
        dir_name (str): Directory name or path.

     Returns:
        str: The absolute path for the directory.

  **gdc.reduce_json_array(j)**

     Recursively go over a JSON and unpack arrays which have only one item, i.e.
     remove unnecessary arrays (brackets).

     Args:
        j (list of dict): A JSON to be reduced.

     Returns:
        list or dict: a reduced JSON with unnecessary array removed.

  **gdc.search(endpoint, in_filter={}, exclude_filter={}, fields=[], expand=[],
  typ='dataframe', method='GET')**

     Search one GDC endpoints and return searching results in a pandas DataFrame
     if possible.

     When searching results cannot be safely converted to a pandas DataFrame,
     results will be returned in the JSON format as it is returned from GDC API.

     Args:
        endpoint (str): One string of GDC API supported endpoint. See:
           https://docs.gdc.cancer.gov/API/Users_Guide/Getting_Started/#api-endpoints

        in_filter (dict, optional): A dict of query conditions which will be
           used to perform the query. Each (key, value) pair represents for one
           ondition. It will be passed to ``simple_and_filter`` for making a
           query filter compatible with GDC API. Please check
           ``simple_and_filter`` function for details.

        exclude_filter (dict, optional): An optional dict of query conditions
           which will be used to perform the query. Each (key, value) pair
           represents for one condition. It will be passed to
           ``simple_and_filter`` for making a query filter compatible with GDC
           API. Please check ``simple_and_filter`` function for details.

        fields (list or str, optional): One or more fields to be queried. Each
           field will be used as a column name in the returned DataFrame. It can
           be a comma separated string or a list of field strings or a
           combination of both.

        expand (list or str, optional): One or more field groups to be
           queried. It can be a comma separated string or a list of field
           strings or a combination of both.

        typ (str): type of search result to return (JSON or dataframe).
           Defaults to ‘dataframe’.

        method (str): HTTP method for the search. Defaults to ‘GET’.
     Returns:
        pandas.core.frame.DataFrame or str: A search result in form of a
           pandas DataFrame or a JSON formatted string, depending on the value
           of ``typ`` and the DataFrame convertibility of JSON.

  **gdc.simple_and_filter(in_dict={}, exclude_dict={})**

     Make a simple GDC API compatible query filter from a dict, in which
     individual conditions are joint by the “and” logic.

     In the return filter, individual conditions in the ``in_dict`` and
     ``exclude_dict`` will be joint by the “and” operator, meaning a hit has to
     match all conditions. Here, a condition can use either a “in” operator
     (specified in the ``in_dict``) or a “exclude” operator (specified in the
     ``exclude_dict``). See details at
     https://docs.gdc.cancer.gov/API/Users_Guide/Search_and_Retrieval/#filters-specifying-the-query

     Args:
        in_dict (dict): A dict describing query conditions with the “in”
           operator. Each (key, value) pair represents for one condition. The
           “key” is the ‘field’ operand. Operator between “key” and “value” is
           “in”.

        exclude_dict (dict): A dict describing query conditions with the
           “exclude” operator. Each (key, value) pair represents for one
           condition. The “key” is the ‘field’ operand. Operator between “key”
           and “value” is “exclude_dict”.

     Returns:
        dict: A dict of filter conforming to GDC API’s format. It should then be
        converted to JSON format and used in the following http request.

* **xena_dataset module**

  This module mainly provides a ``XenaDataset`` class representing one Xena
  matrix in a Xena cohort. Three class, ``GDCOmicset``, ``TCGAPhenoset`` and
  ``TARGETPhenoset`` are derived from ``XenaDataset``, representing genomic
  data, phenotype info of TCGA and phenotype info of TARGET respectively.

  In general, a ``XenaDataset`` class contains 3 methods, ``download``,
  ``transform`` and ``metadata``, which can be used for quickly assembling an
  ETL pipeline importing data into Xena.

  **class xena_dataset.GDCOmicset(projects, xena_dtype, root_dir='.',
  raw_data_dir=None, matrix_dir=None)**

     Bases: ``xena_dataset.XenaDataset``

     GDCOmicset is derived from the ``XenaDataset`` class and represents for a
     Xena matrix whose data is genomic data from GDC.

     This class provides a set of default configurations for downloading and
     transforming GDC data, as well as generating associated metadata for the
     transformed Xena matrix. These default configurations are stored as private
     constants, and they can be checked and/or changed through the following
     attributes: ``gdc_release``, ``gdc_filter``, ``gdc_prefix``,
     ``download_map``, ``read_raw``, ``raws2matrix``, ``metadata_template``, and
     ``metadata_vars``.

     Attributes:
        projects (str or list): One (string) or a list of GDC’s
           “cases.project.project_id”. All corresponding projects will be
           included in this dataset.

        xena_dtype (str): A dataset type supported by this class. To get a
           list of supported types, use ``GDCOmicset.get_supported_dtype()``.

        gdc_release (str): URL to the data release note for the dataset. It
           will be used by the ``metadata`` method when making the metadata for
           this dataset. It is highly recommended that this attribute is set
           explicitly by the user so that it is guaranteed to match the data
           (raw data) underlying this dataset. If it is not available, the most
           recent data release will be queried and used.

        gdc_filter (dict): A filter for querying GDC data underlying this
           dataset. Each item of this dict means to be an “in” operation, with
           its key being one GDC API available field and its value being a
           string or a list of strings. It can be automatically derived from
           ``projects`` and ``xena_dtype`` if it is not assigned explicitly by
           the user when being used. Please check GDC API documentation for
           details.

        gdc_prefix (str): A GDC available file field whost value will be used
           in the filename of corresponding download file. It will be used by
           ``download_map`` for making default download map. It can be
           automatically mapped from ``xena_dtype`` if it is not assigned
           explicitly by the user when being used. Please check ``download_map``
           and GDC API documentation for details.

        download_map (dict): A dict with the key being a URL for one raw data
           to be downloaded and the value being a path for saving downloaded raw
           data. If it hasn’t been assigned explicitly by the user when being
           used, it can be automatically generated by querying through GDC API
           according to ``gdc_filter`` and ``gdc_prefix`` which are based on
           ``projects`` and ``xena_dtype``. Please check ``gdc_filter`` for
           details about querying conditions. Filename of data files, by
           default, will adapt a pattern of “<value of gdc_prefix>.<GDC file
           UUID>.<file extension>”

           It is worth noting that the data transformation process may need an
           ID for every data files. Default ``read_raw`` functions may extract
           the ID from the filename (the first substring when splitting the
           filename by “.”). For example, Xena uses GDC’s
           “cases.samples.submitter_id” for sample ID. Therefore, ``gdc_prefix``
           should be set to “cases.samples.submitter_id” so that data files for
           each sample will be renamed to “<cases.samples.submitter_id>.<file
           UUID>.<file extension>”, allowing the desired sample ID to be
           extracted correctly. Please keep that in mind when trying to define
           your own download dict but use default transformation settings
           (``read_raw`` and ``raws2matrix``). Please check ``read_raw`` and
           ``raws2matrix`` properties, as well as the ``transform`` method for
           details.

        read_raw (callable): A function which accepts only one argument of
           file object and reads it during Xena matrix ``transform``. Its
           return, which is usually a pandas DataFrame, will be put in a list
           which will then be passed to ``raws2matrix``. Defaults, if needed,
           can be mapped from ``xena_dtype``.

        raws2matrix (callable): A function which accepts only one argument of
           a list of ``read_raw`` return(s), merges them into one Xena matrix,
           and processes the merged matrix if needed. Defaults, if needed, can
           be mapped from ``xena_dtype``.

        metadata_template (jinja2.environment.Template or str): A Jinja2
           template for rendering metadata of this dataset. When setting this
           attribute with a string, it will be taken as a path to the template
           file and the corresponding template will be retrieved and assigned to
           this attribute. Defaults, if needed, can be mapped from
           ``xena_dtype``.

        metadata_vars (dict): A dict of variables which will be used (by **
           unpacking) when rendering the ``metadata_template``. Defaults, if
           needed, can be derived from corresponding matrix and ``projects`` and
           ``xena_dtype`` properties.

     ``classmethod get_supported_dtype()``

        Return a list of dataset type codes supported by this class.

  **class xena_dataset.GDCSurvivalset(projects, root_dir='.', raw_data_dir=None,
  matrix_dir=None)**

     Bases: ``xena_dataset.XenaDataset``

     GDCSurvivalset is derived from the ``XenaDataset`` class and represents for
     a Xena matrix of GDC survival data for project(s) of interest.

     In general, survival data is retrieved from GDC API’s “analysis/survival”
     endpoint. This class provides two default configurations, which can be
     checked and/or changed through ``gdc_release`` and ``metadata_vars``, for
     generating metadata for the transformed Xena matrix. The ``download`` and
     ``transform`` methods are overridden by methods specific for GDC survival
     data.

     Attributes:
        gdc_release (str): URL to the data release note for the dataset. It
           will be used by the ``metadata`` method when making the metadata for
           this dataset. It is highly recommended that this attribute is set
           explicitly by the user so that it is guaranteed to match the GDC data
           underlying this dataset. If it is not available, the most recent data
           release will be queried and used.

        metadata_vars (dict): A dict of variables which will be used (by **
           unpacking) when rendering the ``metadata_template``. Defaults, if
           needed, can be derived from corresponding matrix and the ``projects``
           property.

     **download()**

        Retrieve GDC API’s survival data for project(s) in this dataset.

        The survival data is queried and retrieved through GDC API’s
        “analysis/survival” endpoint for project(s) belonging to this dataset.
        JSON query results are converted to a pandas DataFrame and saved as a
        single tab-separated values (“<projects>.GDC_survival.tsv”) file under
        ``raw_data_dir``.

        Returns:
           self: allow method chaining.

     **transform()**

        Transform GDC survival data according to Xena survival data spec

        Only 1 GDC raw survival data (i.e. ``raw_data_list[0]``) will be read
        and used by this transformation. Xena survival data has 6 columns, which
        are “sample”, “_EVENT”, “_TIME_TO_EVENT”, “_OS_IND”, “_OS” and
        “_PATIENT”. “_EVENT” and “_OS_IND” are the same and corresponds to the
        “censored” column in GDC survival data; “_TIME_TO_EVENT” and “_OS” are
        the same and corresponds to the “time” column in GDC survival data;
        “_PATIENT” corresponds to the “submitter_id” column in GDC survival data
        which is the case(patient)’s submitter ID; “sample” contains
        “samples.submitter_id” for corresponding case(patient).

        Returns:
           self: allow method chaining.

  **class xena_dataset.TARGETPhenoset(projects, root_dir='.', raw_data_dir=None,
  matrix_dir=None)**

     Bases: ``xena_dataset.XenaDataset``

     TARGETPhenoset is derived from the ``XenaDataset`` class and represents for
     a Xena matrix whose data is phenotype data of a TARGET project.

     A Xena matrix for TARGET phenotype data is transformed from just the
     clinical data (i.e. without biospecimen data), which is different from that
     of TCGA phenotype data. This class provides a set of default configurations
     for downloading and transforming TARGET clinical, as well as generating
     associated metadata for the transformed Xena matrix. Default configurations
     can be checked and/or changed through ``gdc_release``, ``download_map`` and
     ``metadata_vars``. There is also a protected default method for
     transforming GDC data into Xena matrix. Records in TARGET clinical data are
     per case (patient) based. All samples (IDs) of patients will be queried
     through GDC API and will be merged with clinical data so that transforming
     the data to per sample based.

     Attributes:
        gdc_release (str): URL to the data release note for the dataset. It
           will be used by the ``metadata`` method when making the metadata for
           this dataset. It is highly recommended that this attribute is set
           explicitly by the user so that it is guaranteed to match the data
           (raw data) underlying this dataset. If it is not available, the most
           recent data release will be queried and used.

        download_map (dict): A dict with the key being a URL for one raw data
           to be downloaded and the value being a path for saving downloaded raw
           data. If it hasn’t been assigned explicitly by the user when being
           used, it will, by default, generated by querying through GDC API for
           open access clincial data belong to the ``project`` (won’t check if
           it’s a TARGET project) of this dataset. Files will be rename in the
           pattern of “<GDC file UUID>.<original filename>”.

        metadata_vars (dict): A dict of variables which will be used (by **
           unpacking) when rendering the ``metadata_template``. Defaults, if
           needed, can be derived from corresponding matrix and the ``projects``
           property.

  **class xena_dataset.TCGAPhenoset(projects, root_dir='.', raw_data_dir=None,
  matrix_dir=None)**

     Bases: ``xena_dataset.XenaDataset``

     TCGAPhenoset is derived from the ``XenaDataset`` class and represents for a
     Xena matrix whose data is phenotype data of a TCGA project.

     This class provides a set of default configurations and methods for
     downloading and transforming TCGA clinical and biospecimen data, as well as
     generating associated metadata for the transformed Xena matrix. Default
     configurations can be checked and/or changed through ``gdc_release`` and
     ``metadata_vars``. The ``download`` and ``transform`` methods are
     overridden by methods specific for TCGA phenotype data.

     Attributes:
        gdc_release (str): URL to the data release note for the dataset. It
           will be used by the ``metadata`` method when making the metadata for
           this dataset. It is highly recommended that this attribute is set
           explicitly by the user so that it is guaranteed to match the data
           (raw data) underlying this dataset. If it is not available, the most
           recent data release will be queried and used.

        metadata_vars (dict): A dict of variables which will be used (by **
           unpacking) when rendering the ``metadata_template``. Defaults, if
           needed, can be derived from corresponding matrix and the ``projects``
           property.

     **download()**

        Download GDC’s open access phenotype data for projects in this dataset.

        There are two types of phenotype data on GDC, the clinical data and the
        biospecimen data. Data is selected based solely on the ``projects``
        property. Both clinical and biospecimen data will be downloaded. They
        will be saved under two separated directories. In fact, both downloads
        are delegated by corresponding ``GDCOmicset`` classes. Check the
        ``GDCOmicset`` class for details.

        Returns:
           self: allow method chaining.

     **transform()**

        Transform raw phenotype data into Xena matrix.

        Raw clinical data and biospecimen data will first be transformed
        separately. Then the clinical matrix and biospecimen matrix will be
        merged on “cases.submitter_id” and processed properly.

        Returns:
           self: allow method chaining.

  **class xena_dataset.XenaDataset(projects, xena_dtype, root_dir='.',
  raw_data_dir=None, matrix_dir=None)**

     Bases: ``object``

     XenaDataset represents for one Xena matrix in a Xena cohort.

     This class provides a set of method for downloading and transforming data
     into Xena matrix, as well as generating associated metadata. It also
     provides a set of attributes to control these processes.

     Attributes:
        projects (str or list of str): One or a list of project IDs describing
           study projects included in this dataset.

        xena_dtype (str): A short string (like an ID) describing the type of
           data in this dataset. It is highly recommended, though not required,
           that one dataset has a single type of data.

        root_dir (str): Defines the root directory for this dataset. The
           XenaDataset and the importing process can be highly customized, with
           directories for every data and each step explicitly assigned. You can
           set directories for raw data (through the ``raw_data_dir`` property)
           and Xena matrix (through the ``matrix_dir``property) specifically.
           The ``root_dir`` will be essentially useless under such situation. If
           some or all directories remain unassigned when being used, a default
           directory tree will be used, with a structure like this:

           ::

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
           matrix, metadata, should and highly recommended to be organized and
           saved under the root directory. However, this is neither  required
           nor checked. Setting directory related properties (including
           ``root_dir`` and some properties listed below) will not trigger
           creations of any directories. Directories, if not exist, will only be
           made right before being needed.

           Defaults to “.” which points to current python work directory.

        raw_data_dir (str): A directory for raw data. Please Check the
           ``raw_data_list`` property for its potential usage for defining raw
           data for Xena matrix ``transform``, and check the ``root_dir``
           property for the default “Raw_Data” directory structure. Defaults to
           None.

        matrix_dir (str): A path for saving the transformed Xena matrix for
           this dataset. If the ``matrix_dir`` is not available at the time of
           being used, the ``matrix`` property will be checked first. If the
           ``matrix`` property is available, its directory will be assigned to
           the ``matrix_dir``. If the ``matrix`` property is not available, a
           default path will be assigned according to the default directory
           structure. Check the ``root_dir`` property for the default directory
           structure. Defaults to None.

        download_map (dict): A dict with the key being a URL for one raw data
           to be downloaded and the value being a path for saving downloaded raw
           data.

        raw_data_list (list): A list of file path(s) for all raw data related
           to this dataset. It will be automatically set by the ``download``
           method; or it can be assigned directly as a public attribute. This
           ``raw_data_list`` attribute will be used by the ``transform`` method
           for making a Xena matrix from raw data. If the ``raw_data_list`` is
           not available at the time of being used, the ``raw_data_dir``
           property will be checked. All files under ``raw_data_dir`` will be
           treated as data and used for creating a ``raw_data_list``.

        read_raw (callable): A function used for reading a raw data during
           Xena matrix ``transform``. A valid ``read_raw`` function must accept
           only one argument, which is a file object. Its return will be put in
           a list which will then be passed to ``raws2matrix`` function.

        raws2matrix (callable): A function used for merging multiple raw data
           read by ``read_raw`` into one Xena matrix, as well as processing the
           merged matrix if needed. A valid ``raws2matrix`` must accept only one
           argument, which is a list of ``read_raw`` return(s).

        matrix (str): A path for the Xena matrix of this dataset.
           This attribute will be used but not validated by the ``transform``
           method for saving newly generated Xena matrix. This attribute will
           also be used yet will be validated (i.e. it has to point to a valid
           existing file) by the ``metadata`` method before making metadata
           associated with the Xena matrix and with this dataset. If ``matrix``
           is not available at the time of being used by the ``transform``
           method, it will use ``matrix_dir`` as its directory and will adapte a
           filename with the “projects.xena_type.tsv” pattern. Please check the
           ``matrix_dir`` property to see how it is determined.

        metadata_template (jinja2.environment.Template or str): A Jinja2
           template for rendering metadata of this dataset. When using a string
           to set ``metadata_template``, it will be used as a path to the
           template file and the corresponding template will be retrieved and
           assigned to this attribute.

        metadata_vars (dict): A dict of variables which will be used (by **
           unpacking) for rendering the ``metadata_template``.

     **download(chunk_size=4096)**

        Download file(s) according to the ``download_map`` property.

        A list of paths for downloaded files will be assigned to the
        ``raw_data_list`` property which can be used for Xena matrix
        ``transform`` processing. Check the ``transform`` method for details.

        Args:
           chunk_size (int, optional): The chunk size is the number of bytes
              it should read into memory, when the response is got with
              “stream=True”. Check the documentation of “requests” module for
              details. Defaults to 4096.

        Returns:
           self: allow method chaining.

     **metadata()**

        Make “metadata.json” for Xena data loading.

        A JSON of metadata will be created for the Xena matrix defined by the
        ``matrix`` property. The ``matrix`` property has to point to an existing
        file when this ``metadata`` method is being called. The metadata JSON
        file will be saved under the same directory as the matrix file and named
        with a “.json” postfix appended to the filename of Xena matrix. JSON
        templates for making metatdata are defined by the ``metadata_template``
        property, and variables for rendering the template are defined by the
        ``metadata_vars`` property.

        Returns:
           self: allow method chaining.

     ``root_dir``

        A path of an existing directory for keeping files (raw data, matrix and
        metadata) of this dataset.

     **transform()**

        Transform raw data in a dataset into Xena matrix.

        The transformation process 1) takes in a list of path for raw data; 2)
        open each data based on its file extension; 3) read the file object by
        ``read_func`` and append the readout to a list, which will be 4)
        assembled into a Xena matrix by ``raws2matrix``. The generated Xena
        matrix will be saved at the path defined by the ``matrix`` property.

        Returns:
           self: allow method chaining.

  **xena_dataset.mkdir_p(dir_name)**

     Make the directory as needed: no error if existing.

     Args:
        dir_name (str): Directory name or path.

     Returns:
        str: The absolute path for the directory.

  **xena_dataset.read_biospecimen(fileobj)**

     Extract info from GDC’s biospecimen supplement and re-organize them into a
     pandas DataFrame.

     Args:
        fileobj (file or path): XML file of GDC’s biospecimen supplement.

     Returns:
        pandas.core.frame.DataFrame: Transformed pandas DataFrame.

  **xena_dataset.read_by_ext(filename, mode='r')**

     Automatically decide how to open a file which might be compressed.

     Leveraged codes from the “hook_compressed” function in python’s fileinput
     module.

     Args:
        filename (str): Must contain proper extension indicating the
           compression condition.

        mode (str, optional): To specify the mode in which the file is opened.
           It will be passed to corresponding open function (``open``,
           ``gzip.open`` or ``bz2.BZ2File``); please check them for details.
           Defaults to ‘r’.

     Returns:
        file object: A filehandle to be used with *with*.

  **xena_dataset.read_clinical(fileobj)**

     Extract info from GDC’s clinical supplement and re-organize them into a
     pandas DataFrame.

     Args:
        fileobj (file or path): XML file of GDC’s clinical supplement.

     Returns:
        pandas.core.frame.DataFrame: Transformed pandas DataFrame.

  **xena_dataset.rna_columns_matrix(df_list)**

     Merge and process a list of dataframes to make a Xena data matrix.

     Every dataframe contains data for individual sample. They are merged
     horizontally (``axis=1``; growing with more and more columns). For merged
     matrix, first average columns having the same name and then transform its
     data by log(x + 1).

     Args:
        df_list (list of pandas.core.frame.DataFrame): Input raw data. One
           dataframe is expected to have data for one sample, with column header
           being sample ID (not strictly required).

     Returns:
        pandas.core.frame.DataFrame: Ready to load Xena matrix.

  **xena_dataset.snv_maf_matrix(df_list)**

     Transform pre-sliced GDC’s MAF data into Xena data matrix.

     A new column of DNA variant allele frequencies named “dna_vaf” will
     calculated by division “t_alt_count” / “t_depth”. Columns “t_alt_count” and
     “t_depth” will then be dropped. In the end, columns will be renamed
     accordingly and row index will be set as sample ID.

     Args:
        df (pandas.core.frame.DataFrame): Input raw matrix for mutation data
           (from MAF file).

     Returns:
        pandas.core.frame.DataFrame: Transformed pandas DataFrame.
