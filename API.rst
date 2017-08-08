
xena-GDC-ETL
************


API
===

* **gdc module**

  This module provides basic and minimum necessary functions for carrying out
  data query and download for Xena GDC ETL pipelines.

  **gdc.and_in_filter_constructor(filter_dict)**

     A simple constructor converting a query dictionary into GDC API  specific
     filters.

     Convert a dict of query condition into a diction conforming to GDC  API's
     format. Conditions in the input dict will be combined together by  an AND
     relation. Every (key, value) pair will be joined together with "in"
     operator. If value is not a list, it will be converted to a list first.
     See details at
     https://docs.gdc.cancer.gov/API/Users_Guide/Search_and_Retrieval/#filters-specifying-the-query

     Args:
        filter_dict (dict): A dict describing query conditions. Each (key,
           value) pair represents for one condition. Operator between
           conditions is "AND". Within a condition, the "key" is the 'field'
           operand. Operator between "key" and "value" is "in".

     Returns:
        dict: A dict of filter conforming to GDC API's format. It should then
        be converted to JSON format and used in the following http request.

  **gdc.download(uuids, download_dir='.', chunk_size=4096)**

     Download GDC's open access data according to UUID(s).

     Args:
        uuids (str, list or dict): A single UUID (str), a list of UUIDs (list)
           or a dict whose keys are UUIDs for target file(s). If "uuids" is  str
           or list, downloaded file(s) will be renamed to  "UUID.extension"
           where "extension" is extracted by "get_ext()"  from the original
           filename. Renamed file(s) will be saved at  "download_dir". If
           "uuids" is a dict, the argument "download_dir"  will be ignored;
           values of dict will be paths for saving  corresponding downloaded
           files.

        download_dir (str, optional): The directory for saving downloaded
           file(s) when "uuids" is str or list. It will be ignored if "uuids"
           is dict. Defaults to ".".

        chunk_size (int, optional): The chunk size is the number of bytes it
           should read into memory, when the response is got with
           "stream=True". Check the documentation of "requests" module for
           details. Defaults to 4096.

     Returns:
        list: a list of paths for downloaded files.

  **gdc.get_all_project_info()**

     Get project info for all projects on GDC.

     Returns:
        dict: A dict of project info for all projects in GDC. The key  will be
        project_id and the corresponding value is a dict contains  "project
        name", "primary_site" and "program name" info.

  **gdc.get_ext(file_name)**

     Get all extensions supported by this module in the file name.

     Supported extensions are defined in the constant "_SUPPORTED_FILE_TYPES".
     Multiple extensions will be separated by ".".

     Args:
        file_name (str): The filename will be split by "." and checked from
           left to right. Extensions will be kept starting from the first  (and
           left most) supported extension.

     Returns:
        str: A string of extensions joint by ".".

  **gdc.search(endpoint, filter_dict, fields)**

     Search one GDC endpoints and return searching results in a pandas
     DataFrame if possible.

     When searching results cannot be safely converted to a pandas DataFrame,
     results will be returned in the JSON format as it is returned from GDC
     API.

     Args:
        endpoint (str): One string of GDC API supported endpoint. See:
           https://docs.gdc.cancer.gov/API/Users_Guide/Getting_Started/#api-endpoints

        filter_dict (str): A dict of query conditions which will be used to
           perform the query. Each (key, value) pair represents for one
           condition. It will be passed to "and_in_filter_constructor" for
           making a query filter compatible with GDC API. Please check
           "and_in_filter_constructor" function for details.

        fields (list or str): One or more fields to be queried. Each field
           will be used as a column name in the returned DataFrame. It can be  a
           comma separated string or a list of field strings or a  combination
           of both.

     Returns:
        pandas.core.frame.DataFrame or str: A search result in form of a
           pandas DataFrame or a JSON formatted string.

  **gdc.traverse_field_json(data, field=None)**

     Assuming the nested JSON having a single (set of) data, walk into/down
     this nested JSON and return the value.

     A lot of times, a single GDC data is kept in a highly nested JSON objects.
     For example, in a search of the files endpoint, a sample's submitter ID
     "TCGA-C8-A133-01A" is kept as:

     ::

        "cases": [
          {
            "samples": [
              {
                "submitter_id": "TCGA-C8-A133-01A"
              }
            ]
          }
        ]

     If multiple sets of data are found in the nested JSON, a "ValueError" will
     be raised.

     Args:
        data (str): A string of nested JSON. This JSON should contain only one
           set of data, meaning arrays nested in this JSON should be size 1.
           Otherwise, a "ValueError" will be raised.

        field (str, optional): A string representing the structure of the
           nested JSON. Keys for each level are separated by ".". If None,  the
           nested JSON must contain only one data (rather than one set of
           data),  i.e. only one key at each level of nesting. Otherwise, a
           "ValueError" will be raised. Defaults to None.

     Returns:
        str or else: One data specified by the "field" or the only data
        contained in the nested JSON. It should be str most of the time but  can
        be any types except list and dict.

* **xena_dataset module**

  This module mainly provides a XenaDataset class representing one  Xena  matrix
  in a Xena cohort.

  The XenaDataset class contains 3 methods, ``download_gdc``, ``transform`` and
  ``metadata``, which can be used for quickly assembling an ETL pipeline
  importing GDC data into Xena.

  **class xena_dataset.XenaDataset(projects, xena_dtype, root_dir='.',
  raw_data_dir=None, matrix_dir=None)**

     Bases: ``object``

     XenaDataset represents for one Xena matrix in a Xena cohort.

     This class provides a set of method for downloading and transforming GDC
     data, as well as generating associated metadata for a transformed Xena
     matrix.

     Attributes:
        projects (str or list): One (string) or a list of GDC's
           "cases.project.project_id". All corresponding projects will be
           included in this dataset.

        xena_dtype (str): A dataset type supported by this class. To get a
           list of supported types, use ``XenaDataset.get_supported_dtype()``.

        root_dir (str, optional): Defines the root directory for this dataset.
           By default, all files related to this class, such as raw data,  Xena
           matrix, metadata, should and highly recommended to be  organized and
           saved under this directory. The default directory  structure is:

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

        raw_data_list (list): A list of file paths for all GDC raw data
           related to this dataset. It will be automatically set by the
           ``download_gdc`` method; or it can be assigned directly as a  public
           attribute. This ``raw_data_list`` attribute will be used by
           ``transform`` method for making a Xena matrix from GDC raw data.

        matrix_dir (str, optional): A path for saving the Xena matrix for this
           dataset. Defaults to None. By default, it will be set by the
           ``set_default_dir_tree`` method, based on the ``root_dir``.

        matrix (str, optional): A path for the Xena matrix of this dataset.
           This attribute will be used but not validated (i.e. it can be any
           desired directory and filename) by the ``transform`` method for
           saving newly generated Xena matrix. This attribute will also be  used
           and validated (i.e. it has to point to a valid file) by the
           ``metadata`` method for making metadata assciated with the Xena
           matrix and with this dataset.

           By default, when used by the ``transform`` method, this ``matrix``
           attribute will adapte a pattern as "projects.xena_type.tsv". There
           are two ways to pass custom matrix name: 1) you can set this
           ``matrix`` attribute before calling ``transform``. This allows you
           to customize both directory and filename for the Xena matrix.  2) you
           can pass a ``matrix_name`` argument to the ``transform``  method. You
           can only customize the name (not its directory) of the  Xena matrix.

     **download_gdc()**

        Download GDC's open access data for this dataset.

        Data is selected based on "projects" and "xena_dtype" properties. A  GDC
        query will be built from "xena_dtype" using "__XENA_GDC_DTYPE".  All
        open access data matches the criteria will be put together under  one
        directory as a single dataset. Though this method won't check,  putting
        different types of data into one dataset is NOT recommended.

        By default, data files will be renamed to
        "<cases.samples.submitter_id>.<UUID>.<file extension>" and saved under
        the directory defined by the "raw_data_dir" property. For SNV  datasets,
        the file will be renamed to "<UUID>.<file extension>"  because GDC's
        mutation data (MAF) is one single aggreated data for one  project, not
        per sample based. Check the "_get_gdc_data_dict" method  for details
        about file naming. Check the "root_dir" property for  details about the
        default directory structure.

        A list of paths for downloaded files will be assigned to the
        "raw_data_list" property which can be used for Xena matrix "transform"
        processing. Check the "transform" method for details.

        Returns:
           self: allow method chaining.

     ``classmethod get_supported_dtype()``

        Return a list of dataset type codes supported by this class.

     **metadata()**

        Make "metadata.json" for Xena data loading

        One metadata will be created for one Xena matrix defined by the
        "matrix" property. The metadata JSON file will be saved under the same
        directory as the matrix file and named with a ".json" postfix appended
        to the filename of Xena matrix. The "transform" method will generate a
        Xena matrix from raw data and assign it to "matrix"; or "matrix" can  be
        assigned directly. JSON templates for metatdata are defined by
        "__METADATA_TEMPLATE_DIR" and "__METADATA_TEMPLATE" constants;
        pre-defined Xena cohort names are defined by the "__XENA_COHORT"
        constant; data type labels in the metadata are defined by the
        "__METADATA_GDC_TYPE" constant. Specific raw data (MAF) urls for
        "Masked Somatic Mutation" data will be queried according to "projects"
        and "xena_dtype" properties.

        Returns:
           self: allow method chaining.

     **set_default_dir_tree(root_dir, reset_default=False)**

        Set the default directory structure for this dataset.

        For default directory structure, please check the ``root_dir``
        property.

        Args:
           root_dir (str): The root directory for keep the new default
              directory structure of this dataset.

           reset_default (bool): Whether to overide current settings for
              "raw_data_dir" and "matrix_dir" properties if they have  already
              been set. Defaults to False.

        Returns:
           self: allow method chaining.

     **transform(raw_data_list=None, matrix_name=None)**

        Transform raw data in a dataset into Xena matrix

        Args:
           raw_data_list (str or list, optional): One (str) or a list of raw
              data file(s) used for building Xena matrix. This  "raw_data_list"
              argument has a higher priority than the  "raw_data_list" property,
              i.e. if "raw_data_list" is not None,  it will be used to define or
              overwrite the "raw_data_list"  property. If this "raw_data_list"
              argument is None, the  "raw_data_list" property will be checked
              first. If  "raw_data_list" property is defined and is not None, it
              will  be used. The "raw_data_list" property can be assigned
              directly  or be set by the "download" method. Check the "download"
              method for details. If the "raw_data_list" property is not
              usable, the "raw_data_dir" property will be checked. All files
              under this directory will be treated as data and used for
              building Xena matrix. Defaults to None.

           matrix_name (str, optional): File name of the Xena matrix. The
              matrix will be save under the directory defined by the
              "matrix_dir" property. If None, default filename will be used,
              which has a pattern as "projects.xena_type.tsv". Full path of
              this matrix will be assigned to the "matrix" property which  can
              be used for making metadata. Check the "metadata" method  for
              details. Defaults to None.

        Returns:
           self: allow method chaining.

  **xena_dataset.mkdir_p(dir_name)**

     Make the directory as needed: no error if existing.

     Args:
        dir_name (str): Directory name or path.

     Returns:
        str: The absolute path for the directory.

  **xena_dataset.process_average_log(df)**

     Process Xena data matrix by first averaging columns having the same  name
     and then transform it by log(x + 1).

     Args:
        df (pandas.core.frame.DataFrame): Input raw data matrix.

     Returns:
        pandas.core.frame.DataFrame: Transformed pandas DataFrame.

  **xena_dataset.process_maf(df)**

     Transform pre-sliced GDC's MAF data into Xena data matrix.

     A new column of DNA variant allele frequencies named "dna_vaf" will
     calculated by division "t_alt_count" / "t_depth". Columns "t_alt_count"
     and "t_depth" will then be dropped. In the end, columns will be renamed
     accordingly and row index will be set as sample ID.

     Args:
        df (pandas.core.frame.DataFrame): Input raw matrix for mutation data
           (from MAF file).

     Returns:
        pandas.core.frame.DataFrame: Transformed pandas DataFrame.

  **xena_dataset.read_by_ext(filename, mode='r')**

     Automatically decide how to open a file which might be compressed.

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
        file object: A filehandle to be used with *with*.
