import json
import os
from mock import Mock, patch
import time
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

import pandas as pd

from xena_gdc_etl import utils

PATH = "tests/fixtures/"  # path to static files


def test_mkdir_p():
    input_ = "tests/tmp_dir"
    result = utils.mkdir_p(input_)
    assert os.path.isdir(input_) is True
    assert input_ in result
    os.rmdir(input_)


@patch("sys.stdout", new_callable=StringIO)
def test_equal_matrices(mocked_print):
    utils.equal_matrices(PATH + "df1.csv", PATH + "df2.csv")
    utils.equal_matrices(PATH + "df1.csv", PATH + "df3.csv")
    assert mocked_print.getvalue() == "Not equal.\nEqual.\n"


@patch.object(utils, 'time', Mock(wraps=time))
def test_metadata():
    utils.time.gmtime.return_value = time.struct_time(
        (2019, 4, 28, 0, 0, 0, 0, 0, 0)
    )
    path_to_df = PATH + "HTSeq-FPKM-UQ.csv"
    utils.metadata(path_to_df, "htseq_fpkm-uq")
    with open(PATH + "HTSeq-FPKM-UQ.csv.json", "r") as actual:
        actual = json.load(actual)
    os.unlink(PATH + "HTSeq-FPKM-UQ.csv.json")
    with open(PATH + "HTSeq-FPKM-UQ.json", "r") as expected:
        expected = json.load(expected)
    assert actual == expected


def test_merge_xena():
    name = "merged_GDC_TARGET-CCSK.tsv"
    filelists = ["merge-xena1.csv", "merge-xena2.csv"]
    cohort = "GDC TARGET-CCSK"
    datatype = "GDC_phenotype"
    outdir = PATH
    utils.handle_merge_xena(name, filelists, cohort, datatype, outdir)
    with open(PATH + "MergedCohort04292019.GDC_phenotype.tsv", "r") as actual:
        actual = pd.read_csv(actual)
    with open(PATH + name, "r") as expected:
        expected = pd.read_csv(expected)
    os.unlink(PATH + name)
    actual.equals(expected)
