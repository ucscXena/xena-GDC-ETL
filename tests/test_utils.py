import os
from mock import patch, call
try:
    from StringIO import StringIO 
except ImportError:
    from io import StringIO
from unittest import TestCase

import pandas as pd

from xena_gdc_etl import utils


def test_mkdir_p():
    input_ = "tests/tmp_dir"
    result = utils.mkdir_p(input_)
    assert os.path.isdir(input_) is True
    assert input_ in result
    os.rmdir(input_)


@patch("sys.stdout", new_callable=StringIO)
def test_equal_matrices(mocked_print):
        utils.equal_matrices(
            "tests/fixtures/df1.csv",
            "tests/fixtures/df2.csv"
        )
        utils.equal_matrices(
            "tests/fixtures/df1.csv",
            "tests/fixtures/df3.csv"
        )
        assert mocked_print.getvalue() == "Not equal.\nEqual.\n"
