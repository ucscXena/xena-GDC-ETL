#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This script provides a command line tool for checking GDC's list of updated
files and summarize impacted project(s), data_type(s) and
analysis.workflow_type(s). Example usage is::

    python gdc_check_new.py [URL]

Please check help message with "-h" option for details.
"""

# Ensure Python 2 and 3 compatibility
from __future__ import division
from __future__ import print_function

import argparse
import sys

import pandas as pd

import gdc

def main():
    parser = argparse.ArgumentParser(
            description="Check GDC's list of updated files and summarize "
                        "impacted project(s), data_type(s) and "
                        "analysis.workflow_type(s)."
        )
    parser.add_argument('url', type=str, metavar='URL',
                        help='URL for GDC\'s list of updated files. It can be '
                             'a compressed file with a supported extension, '
                             'which includes ".gz", ".bz2", ".zip", or "xz". '
                             'New files should be listed under a column named '
                             'by "New File UUID".')
    args = parser.parse_args()
    new_uuids = pd.read_table(args.url)['New File UUID'].tolist()
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


if __name__ == '__main__':
    main()
