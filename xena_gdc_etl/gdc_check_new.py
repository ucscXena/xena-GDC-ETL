import sys

import pandas as pd

import xena_gdc_etl.gdc as gdc


def gdc_check_new(url):
    """
    This function provides a command line tool for checking GDC's list of
    updated files and summarize impacted project(s), data_type(s) and
    analysis.workflow_type(s).
    """
    new_uuids = pd.read_table(url)['New File UUID'].tolist()
    df_list = []
    for uuids in (new_uuids[i:i + 20000]
                  for i in range(0, len(new_uuids), 20000)):
        df = gdc.search(
                'files',
                in_filter={'access': 'open', 'file_id': uuids},
                fields=['cases.project.project_id', 'data_type',
                        'analysis.workflow_type'],
                method='POST'
            )
        try:
            df['cases'] = df['cases'].map(
                    lambda c: ', '.join({p['project']['project_id']
                                         for p in c})
                )
        except:
            pass
        df_list.append(df)
    df = pd.concat(df_list, axis=0)
    try:
        df = df.drop('id', axis=1)
    except:
        pass
    try:
        df = df.drop_duplicates()
    except:
        pass
    df.to_csv(sys.stdout, sep='\t', index=False)
