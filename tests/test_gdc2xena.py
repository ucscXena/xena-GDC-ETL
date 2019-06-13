import json
import os

import pytest

from xena_gdc_etl import gdc2xena


PATH = "tests/fixtures/"


@pytest.mark.CI
def test_logging():
    gdc2xena.gdc2xena(
        root_dir=".",
        projects=["SOME_PROJECT1", "SOME_PROJECT2"],
        xena_dtypes=["htseq_counts", "methylation450"],
    )
    with open(PATH + "unfinished_test.json", "r") as infile:
        expected = json.load(infile)
    with open("unfinished.json", "r") as infile:
        actual = json.load(infile)
    assert actual == expected
    os.unlink("unfinished.json")
