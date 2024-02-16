import argparse
import os
import pandas as pd

from xena_gdc_etl import gdc


def create_parser():
    """
    Construct the program options.
    """

    parser = argparse.ArgumentParser(
        description='Postprocess CPTAC-3 data',
    )
    parser.add_argument(
        '-p',
        '--project',
        type=str,
        required=True,
        help='The project name.',
    )
    parser.add_argument(
        '-t',
        '--datatype',
        type=str,
        required=True,
        help='The Xena data type of the file.',
    )
    parser.add_argument(
        '-f',
        '--file',
        type=str,
        required=True,
        help='Path to CPTAC-3 data file. File will be read by pandas.read_csv with sep="\t".',
    )

    return parser


def get_gdc_data(project):
    """Get cases and samples through the GDC API for CPTAC-3. 
    
    Args:
        project (str): The project name.

    Returns:
        dict: A dict whose keys are cases' submitter_ids and values are samples'
        submitter_id.
    """

    gdc_data = {}
    cases = gdc.search('cases', 
                        in_filter={'project.project_id': project}, 
                        fields=['submitter_id', 'samples.submitter_id'],
                        typ='json',
                    )
    for id in cases: 
        samples = [s['submitter_id'] for s in id['samples']]
        gdc_data[id['submitter_id']] = samples

    return gdc_data


def postprocess(data_type, df, gdc_data):
    """Postprocess the Xena data by mapping the samples to cases and removing duplicated data.

    Args: 
        data_type (str): The Xena data type of the file.
        df (pandas.core.frame.DataFrame): The Xena dataframe to be processed.
        gdc_data (dict): A dict whose keys are cases' submitter_ids and values are samples'
            submitter_id.

    Returns:
        pandas.core.frame.DataFrame: Transformed pandas DataFrame.
    """

    sample_columns = ['star_counts', 'star_tpm', 'star_fpkm', 'star_fpkm-uq', 'gene-level_ascat-ngs', 'methylation_epic']
    if data_type in sample_columns:
        for column in df.columns:
            for key, value in gdc_data.items():
                if column in value: 
                    df.rename(columns={column: key}, inplace=True)
        processed_df = df.T.groupby(level=0).first().T 
    else:
        df.set_index('sample', inplace=True)
        for index, row in df.iterrows():
            for key, value in gdc_data.items():
                if index in value: 
                    df.rename(index={index: key}, inplace=True) 
        processed_df = df.drop_duplicates()
        
    return processed_df


def main():
    parser = create_parser()
    options = parser.parse_args()
    file_df = pd.read_csv(options.file, sep='\t')
    filename = os.path.basename(options.file)
    new_dir = os.path.join(os.path.dirname(os.path.dirname(options.file)), 'Postprocessed_Matrices')
    if not os.path.exists(new_dir):
        os.mkdir(new_dir)
    gdc_data = get_gdc_data(options.project)
    processed_df = postprocess(options.datatype, file_df, gdc_data)
    processed_df.to_csv(os.path.join(new_dir, filename), sep='\t')

if __name__ == '__main__':
    main()
