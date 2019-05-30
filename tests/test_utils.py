import os

import pandas as pd

from xena_gdc_etl import utils

PATH = "tests/fixtures/"  # path to static files


def test_mkdir_p():
    input_ = "tests/tmp_dir"
    result = utils.mkdir_p(input_)
    assert os.path.isdir(input_) is True
    assert input_ in result
    os.rmdir(input_)


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
