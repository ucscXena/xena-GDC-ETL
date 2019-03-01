#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This script provides a command line tool for testing the equality of 2 Xena
matrices.
"""

# Ensure Python 2 and 3 compatibility
from __future__ import print_function

import argparse

import pandas as pd
from pandas.util.testing import assert_frame_equal


def main():
    parser = argparse.ArgumentParser(
            description='Test the equality of 2 Xena matrices.'
        )
    parser.add_argument('df1', type=str,
                        help='Directory for the first matrix.')
    parser.add_argument('df2', type=str,
                        help='Directory for the second matrix.')
    args = parser.parse_args()
    df1 = pd.read_table(args.df1, header=0,
                        index_col=0).sort_index(axis=0).sort_index(axis=1)
    df2 = pd.read_table(args.df2, header=0,
                        index_col=0).sort_index(axis=0).sort_index(axis=1)
    try:
        assert_frame_equal(df1, df2)
        print('Equal.')
    except:  # appeantly AssertionError doesn't catch all
        print('Not equal.')

if __name__ == '__main__':
    main()
