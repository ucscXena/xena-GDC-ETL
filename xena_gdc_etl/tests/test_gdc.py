from gdc import get_project_info


def test_get_project_info():
    assert 'TCGA-BRCA' in get_project_info().index
    assert get_project_info(['TCGA-BRCA']).index.tolist() == ['TCGA-BRCA']
