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

    sample_columns = ['star_counts', 'star_tpm', 'star_fpkm', 'star_fpkm-uq', 'mirna', 'gene-level_ascat-ngs', 'methylation_epic', 'methylation_epic_v2']
    if data_type in sample_columns:
        for column in df.columns:
            for key, value in gdc_data.items():
                if column in value: 
                    msg = '{} has been changed to {}.'
                    print(msg.format(column, key))
                    df.rename(columns={column: key}, inplace=True)
        duplicated_cols = df.columns.duplicated()
        duplicated = {list(df.columns)[i] for i, val in enumerate(duplicated_cols) if val}
        processed_df = df.loc[:, ~duplicated_cols]
        msg = '\r{} samples had duplicated data. The samples are: {}.'
        print(msg.format(len(duplicated), ', '.join(duplicated)))
    else:
        df.set_index('sample', inplace=True)
        for index, row in df.iterrows():
            for key, value in gdc_data.items():
                if index in value: 
                    msg = '{} has been changed to {}.'
                    print(msg.format(index, key))
                    df.rename(index={index: key}, inplace=True) 
        df.reset_index(inplace=True)
        duplicated_rows = df.duplicated()
        duplicated = {val for i, val in enumerate(df[duplicated_rows]['sample'])}
        processed_df = df.drop_duplicates()
        msg = '\r{} samples had duplicated data. The samples are: {}.'
        print(msg.format(len(duplicated), ', '.join(duplicated)))
    processed_df.set_index(df.columns[0], inplace=True)

    return processed_df


def main():
    parser = create_parser()
    options = parser.parse_args()
    file_df = pd.read_csv(options.file, sep='\t')
    filename = os.path.basename(options.file)
    project = options.project
    datatype = options.datatype
    new_dir = os.path.join(os.path.dirname(os.path.dirname(options.file)), 'Postprocessed_Matrices')
    print('{} data will be postprocessed for the following datatype: {}.'.format(project, datatype))
    if not os.path.exists(new_dir):
        print('Creating new directory for postprocessed matrices ...')
        os.mkdir(new_dir)
    if project == 'BEATAML1.0-COHORT':
        if datatype.startswith('star_'):
            processed_df = file_df.set_index('Ensembl_ID')
            file_df.set_index('Ensembl_ID', inplace=True)
            for column in file_df.columns:
                msg = '{} has been changed to {}.'
                print(msg.format(column, column[:-1]))
                processed_df.rename(columns={column: column[:-1]}, inplace=True)
        else:
            print('Old sample column:\n')
            print(file_df['sample'])
            file_df['sample'] = file_df['sample'].str[:-1]
            print('New sample column:\n')
            print(file_df['sample'])
            processed_df = file_df.set_index('sample')
    elif project == 'CMI-MPC' or project == 'CMI-MBC':
        if datatype.startswith('star_'):
            processed_df = file_df.set_index('Ensembl_ID')
            file_df.set_index('Ensembl_ID', inplace=True)
            for column in file_df.columns:
                msg = '{} has been changed to {}.'
                print(msg.format(column, column[:-4]))
                processed_df.rename(columns={column: column[:-4]}, inplace=True)
        else: 
            print('Old sample column:\n')
            print(file_df['sample'])
            file_df['sample'] = file_df['sample'].str[:-4]
            print('New sample column:\n')
            print(file_df['sample'])
            processed_df = file_df.set_index('sample')
    else: 
        gdc_data = get_gdc_data(project)
        processed_df = postprocess(datatype, file_df, gdc_data)
    status = 'Postprocessed {} data is ready for {}.'
    print(status.format(datatype, project))
    processed_df.to_csv(os.path.join(new_dir, filename), sep='\t')
    print('Postprocessed matrix is saved at {}.'.format(os.path.join(new_dir, filename)))


if __name__ == '__main__':
    main()
