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


def test_get_intersection():
    list_1 = ["hello", "world", "earth"]
    list_2 = ["hello", "world", "world", "earth"]
    actual = utils.get_intersection(list_1, list_2)
    expected = ["hello", "world"]
    assert actual == expected


def test_extract_leaves():
    data = {
        "analytes": [
            {
                "aliquots": [
                    {"aliquot_id": "35f8d837-f78c-4f88-a4cd-50ea2f9f9437"},
                    {"aliquot_id": "48de3d5a-35d8-456d-a23e-e2dc25a840ac"},
                    {"aliquot_id": "e6443c75-5f1d-45c6-8796-992dd51a0496"},
                ]
            }
        ]
    }
    actual = []
    utils.extract_leaves(
        data=data,
        field_name="analytes.aliquots.aliquot_id",
        leaves=actual,
    )
    expected = [
        "35f8d837-f78c-4f88-a4cd-50ea2f9f9437",
        "48de3d5a-35d8-456d-a23e-e2dc25a840ac",
        "e6443c75-5f1d-45c6-8796-992dd51a0496",
    ]
    assert actual == expected
