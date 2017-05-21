# xena-GDC-ETL (To be finished...)
Extract, transform and load GDC data onto UCSC Xena.
- GDC module

   Handle data search and downloading on GDC.
   
   |Function|Description|
   |---|---|
   |get_all_project_ids|Get project ids (as a list) for all projects on GDC.|
   |get_files_uuids|Get UUIDs (as a list) for files matching query conditions.|
   |and_eq_filter_constructor|A simple constructor converting a query dictionary into GDC API specific [`filters`](https://docs.gdc.cancer.gov/API/Users_Guide/Search_and_Retrieval/#filters-specifying-the-query).|
   |download_data|Download open access data from GDC according to input UUIDs.|
   |label_files|Label a list of file UUIDs with their corresponding field of selection, such as their corresponding "aliquots.submitter_id".|
   |get_all_case_info|Get some basic information for all cases on GDC.|
   
- Xena module

   Has the main pipeline for importing GDC data into Xena. It will use the GDC module to retrieve data. Within the module, it has functions for data transformation and assembly for individual data types.
   
   |Function|Description|
   |---|---|
   |xena_matrix_update_RNA|Transformation (log2(x+1)) and data matrix assembly for RNA-seq data (HTSeq - Counts, HTSeq - FPKM and HTSeq - FPKM-UQ).|
   |main|The main pipeline for importing GDC data into Xena.|

### To be finished...
