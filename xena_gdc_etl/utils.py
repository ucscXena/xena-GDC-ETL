from __future__ import print_function
import os

import pandas as pd
from pandas.util.testing import assert_frame_equal


def mkdir_p(dir_name):
    """Make the directory as needed: no error if existing.

    Args:
        dir_name (str): Directory name or path.

    Returns:
        str: The absolute path for the directory.
    """

    dir_path = os.path.abspath(dir_name)
    try:
        os.makedirs(dir_path)
    except OSError:
        if not os.path.isdir(dir_path):
            raise
    return dir_path


def equal_matrices(df1, df2):
    """This function is used for testing the equality of 2 Xena matrices

    Args:
        xena matrices df1 and df2 whose equality is to be checked
    """
    df1 = pd.read_table(df1, header=0,
                        index_col=0).sort_index(axis=0).sort_index(axis=1)
    df2 = pd.read_table(df2, header=0,
                        index_col=0).sort_index(axis=0).sort_index(axis=1)
    try:
        assert_frame_equal(df1, df2)
        print('Equal.')
    except:  # appeantly AssertionError doesn't catch all
        print('Not equal.')
