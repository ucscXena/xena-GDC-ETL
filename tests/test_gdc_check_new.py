try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

import pandas as pd
import pytest
from mock import patch

from xena_gdc_etl import gdc_check_new


@pytest.mark.CI
@patch("sys.stdout", new_callable=StringIO)
def test_gdc_check_new(mocked_stdout):
    url = "https://docs.gdc.cancer.gov/Data/Release_Notes/DR9.0_files_swap.txt.gz"  # noqa
    gdc_check_new.gdc_check_new(url)
    actual = StringIO(mocked_stdout.getvalue())
    actual = pd.read_csv(actual, sep='\t')
    expected = pd.read_csv("tests/fixtures/gdc_check_new_DR9.0_files_swap.csv",
                           sep='\t')
    expected = expected.head()
    actual.equals(expected)
