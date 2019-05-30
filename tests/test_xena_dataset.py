import pandas as pd

from xena_gdc_etl import xena_dataset


PATH = "tests/fixtures/xena_dataset/"


def test_read_biospecimen():
    file_name = PATH + "nationwidechildrens.org_biospecimen.TCGA-AR-A0TQ"
    expected = xena_dataset.read_biospecimen(file_name + ".xml")
    actual = pd.read_csv(file_name + ".csv", sep="\t")
    expected.equals(actual)


def test_read_clinical():
    file_name = PATH + "nationwidechildrens.org_omf.TCGA-RW-A68A"
    expected = xena_dataset.read_clinical(file_name + ".xml")
    actual = pd.read_csv(file_name + ".csv", sep="\t")
    expected.equals(actual)
