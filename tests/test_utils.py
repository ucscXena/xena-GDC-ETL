import json
import os
from mock import Mock, patch
import time
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

import pandas as pd
import pytest

from xena_gdc_etl import utils


def test_mkdir_p():
    input_ = "tests/tmp_dir"
    result = utils.mkdir_p(input_)
    assert os.path.isdir(input_) is True
    assert input_ in result
    os.rmdir(input_)


@pytest.fixture
@patch("sys.stdout", new_callable=StringIO)
def test_equal_matrices(mocked_print):
    utils.equal_matrices("df1.csv", "df2.csv")
    utils.equal_matrices("df1.csv", "df3.csv")
    assert mocked_print.getvalue() == "Not equal.\nEqual.\n"


@pytest.fixture
@patch.object(utils, 'time', Mock(wraps=time))
def test_metadata():
    utils.time.gmtime.return_value = time.struct_time(
        (2019, 4, 28, 0, 0, 0, 0, 0, 0)
    )
    path_to_df = "HTSeq-FPKM-UQ.csv"
    utils.metadata(path_to_df, "htseq_fpkm-uq")
    with open("HTSeq-FPKM-UQ.csv.json", "r") as actual:
        actual = json.load(actual)
    os.unlink("HTSeq-FPKM-UQ.csv.json")
    with open("HTSeq-FPKM-UQ.json", "r") as expected:
        expected = json.load(expected)
    assert actual == expected


@pytest.fixture
def test_merge_xena():
    name = "merged_GDC_TARGET-CCSK.tsv"
    filelists = ["merge-xena1.csv", "merge-xena2.csv"]
    cohort = "GDC TARGET-CCSK"
    datatype = "GDC_phenotype"
    path = "tests/fixtures"
    outdir = path
    utils.handle_merge_xena(name, filelists, cohort, datatype, outdir)
    with open("MergedCohort04292019.GDC_phenotype.tsv", "r") as actual:
        actual = pd.read_csv(actual)
    with open(path + name, "r") as expected:
        expected = pd.read_csv(expected)
    os.unlink(path + name)
    actual.equals(expected)
