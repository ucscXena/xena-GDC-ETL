import os
import shutil

import pytest

from xena_gdc_etl import xena_dataset


@pytest.mark.CI
def test_download():
    dataset = xena_dataset.GDCOmicset(
        projects="TCGA-BRCA",
        root_dir=r".",
        xena_dtype="htseq_counts",
    )
    dataset.download(no_of_files=1)
    files = next(os.walk("TCGA-BRCA/Raw_Data/htseq_counts"))[2]
    assert len(files) == 1
    shutil.rmtree("./TCGA-BRCA", ignore_errors=True)
