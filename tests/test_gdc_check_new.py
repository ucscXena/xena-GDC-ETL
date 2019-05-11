try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

import pandas as pd
import pytest

from xena_gdc_etl import gdc_check_new


@pytest.mark.CI
def test_gdc_check_new(capfd):
    url = "https://docs.gdc.cancer.gov/Data/Release_Notes/DR9.0_files_swap.txt.gz"  # noqa
    gdc_check_new.gdc_check_new(url)
    out, err = capfd.readouterr()
    actual = pd.read_csv(StringIO(out), sep='\t')
    expected = pd.read_csv("tests/fixtures/gdc_check_new_DR9.0_files_swap.csv",
                           sep='\t')
    expected = expected.head()
    actual.equals(expected)
