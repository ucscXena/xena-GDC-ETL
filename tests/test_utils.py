import os

import pandas as pd
import pytest

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


@pytest.mark.parametrize(
    'raw,path,expect',
    [
        ({'a': 'b'}, '', []),
        ([{'a': 'b'}, {'c': [{'d': 'e'}]}], 'c.d', ['e']),
        (
            {
                "submitter_id": "TCGA-AX-A064",
                "samples": [{"portions": [{"analytes": [{"aliquots": [
                    {
                        "submitter_id": "TCGA-AX-A064-10A-01W-A027-09",
                        "aliquot_id": "14163707-5628-4c1b-9efd-87f7dd11d300",
                    },
                    {
                        "submitter_id": "TCGA-AX-A064-10A-01W-A028-08",
                    }
                ]}]}]}]
            },
            'samples.portions.analytes.aliquots.aliquot_id',
            [
                '14163707-5628-4c1b-9efd-87f7dd11d300',
            ]
        ),
    ]
)
def test_get_json_objects(raw, path, expect):
    assert utils.get_json_objects(raw, path) == expect
