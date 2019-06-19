import json
import os
import shutil

import pytest

from xena_gdc_etl import gdc2xena


PATH = "tests/fixtures/"


@pytest.mark.CI
def test_unfinished_generation():
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


@pytest.mark.CI
def test_deletion_of_raw_data():
    os.mkdir("test_deletion")
    gdc2xena.gdc2xena(
        root_dir="test_deletion",
        projects=["TCGA-UVM"],
        xena_dtypes=["somaticsniper_snv"],
        delete_raw_data=True,
    )
    path = "/test_deletion/TCGA-UVM/Raw-Data/"
    all_files = []
    for root, _, files in os.walk(path):
        for file in files:
            all_files.append(os.path.join(root, file))
    assert all_files == []
    shutil.rmtree("test_deletion")
