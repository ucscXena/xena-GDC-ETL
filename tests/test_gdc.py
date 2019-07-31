try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
import os

import pandas as pd
import pytest

from xena_gdc_etl import gdc


def test_simple_and_filter():
    in_dict_1 = {}
    exclude_dict_1 = {}
    output_1 = gdc.simple_and_filter(in_dict_1, exclude_dict_1)
    assert output_1 == in_dict_1
    in_dict_2 = {'a': 'b'}
    exclude_dict_2 = {'c': 'd'}
    expected = {
        "content": [
            {"content": {"field": "a", "value": ["b"]}, "op": "in"},
            {"content": {"field": "c", "value": ["d"]}, "op": "exclude"},
        ],
        "op": "and",
    }
    actual = gdc.simple_and_filter(in_dict_2, exclude_dict_2)
    assert actual == expected


def test_reduce_json_array():
    input_1 = [{'a': 'hello', 'b': [1, 2, 3], 'c': [10]}]
    input_2 = [{'a': 'b'}]
    actual_1 = gdc.reduce_json_array(input_1)
    expected_1 = {"a": "hello", "b": [1, 2, 3], "c": 10}
    assert actual_1 == expected_1
    actual_2 = gdc.reduce_json_array(input_2)
    expected_2 = {'a': 'b'}
    assert actual_2 == expected_2


def test_get_ext():
    input_1 = "txt.vcf.xls"
    actual_1 = gdc.get_ext(input_1)
    expected_1 = "txt.vcf.xls"
    assert actual_1 == expected_1
    input_2 = "abc.xyz.pqr"
    expected_2 = "pqr"
    actual_2 = gdc.get_ext(input_2)
    assert actual_2 == expected_2
    input_3 = "name.txt.vcf.xls"
    actual_3 = gdc.get_ext(input_3)
    expected_3 = "txt.vcf.xls"
    assert actual_3 == expected_3


@pytest.mark.CI
def test_download():
    uuid = "53a637ce-8aaf-4cec-b02d-89202bbb0890"
    gdc.download(uuid, download_dir="./tests")
    file_path = "./tests/" + uuid + ".svs"
    assert os.path.isfile(file_path) is True
    os.unlink(file_path)


@pytest.mark.CI
def test_get_project_info():
    project_name = "TCGA-THCA"
    assert 'TCGA-BRCA' in gdc.get_project_info().index
    assert gdc.get_project_info(['TCGA-BRCA']).index.tolist() == ['TCGA-BRCA']
    actual = gdc.get_project_info([project_name]).head()
    expected = {
        "id": ["TCGA-OV", "Ovarian", "Serous"],
        "name": "Cystadenocarcinoma",
        "primary_site": "Ovary",
        "program.name": "TCGA",
        "project_id": "TCGA-OV",
    }
    actual.equals(expected)


@pytest.mark.CI
def test_get_samples_clinical():
    project_id = "TCGA-OV"
    actual = gdc.get_samples_clinical(project_id)
    in_filter = {"project.project_id": project_id}
    expected = gdc.search(
        in_filter=in_filter,
        fields="case_id",
        endpoint="cases"
    )
    assert expected["case_id"][0] == actual["case_id"][0]


@pytest.mark.CI
def test_search():
    endpoint = "cases"
    in_filter = {"project.project_id": "TARGET-CCSK"}
    fields = ["submitter_id"]
    actual = gdc.search(endpoint=endpoint, in_filter=in_filter, fields=fields)
    expected = {
        "id": "d1a15919-f5e2-5e60-aed9-cb52a8b4a7a1",
        "target": "TARGET-51-PAKWMM",
    }
    actual.equals(expected)
    with pytest.raises(ValueError) as exception_info:
        gdc.search(
            endpoint=endpoint, in_filter=in_filter, fields=fields, method="PUT"
        )
    error_str = 'Invalid method: PUT\n method must be either "GET" or "POST".'
    assert exception_info.value.args[0] == error_str


@pytest.mark.CI
def test_gdc_check_new(capfd):
    url = "https://docs.gdc.cancer.gov/Data/Release_Notes/DR9.0_files_swap.txt.gz"  # noqa
    new_file_uuids = pd.read_csv(url, sep='\t')['New File UUID'].tolist()
    gdc.gdc_check_new(new_file_uuids)
    out, err = capfd.readouterr()
    actual = pd.read_csv(StringIO(out), sep='\t')
    expected = pd.read_csv(
        "tests/fixtures/gdc_check_new_DR9.0_files_swap.csv", sep='\t'
    )
    expected = expected.head()
    actual.equals(expected)


@pytest.mark.CI
@pytest.mark.parametrize(
    'endpoint,input_field,output_field,input_values,expect',
    [
        (
            'cases',
            'samples.portions.analytes.aliquots.aliquot_id',
            'samples.submitter_id',
            [
                '35f8d837-f78c-4f88-a4cd-50ea2f9f9437',
                '48de3d5a-35d8-456d-a23e-e2dc25a840ac',
            ],
            {
                '35f8d837-f78c-4f88-a4cd-50ea2f9f9437': ['TCGA-4G-AAZO-11A'],
                '48de3d5a-35d8-456d-a23e-e2dc25a840ac': ['TCGA-4G-AAZO-11A'],
            }
        ),
        (
            'cases',
            'samples.portions.analytes.aliquots.aliquot_id',
            'invalid.path',
            ['35f8d837-f78c-4f88-a4cd-50ea2f9f9437'],
            {
                '35f8d837-f78c-4f88-a4cd-50ea2f9f9437': [],
            }
        ),
        (
            'files',
            'analysis.input_files.file_id',
            'analysis.input_files.file_name',
            [
                '4b0b0fb4-99c1-4b86-be77-0bce51ced90d',
                'cd01b316-357d-4535-b7ae-a5645e619b35',
            ],
            {
                '4b0b0fb4-99c1-4b86-be77-0bce51ced90d': ['COTES_p_TCGAaffxB8_9a_S_GenomeWideSNP_6_F08_293030.nocnv_grch38.seg.v2.txt'],  # noqa: E501
                'cd01b316-357d-4535-b7ae-a5645e619b35': ['QUANT_p_TCGA_Batch12_AFFX_GenomeWideSNP_6_H02_437614.nocnv_grch38.seg.v2.txt'],  # noqa: E501
            }
        ),
        (
            'projects',
            'project_id',
            'disease_type',
            [
                'NCICCR-DLBCL',
                'CPTAC-3',
                'BEATAML1.0-CRENOLANIB',
                'CTSP-DLBCL1',
            ],
            {
                'NCICCR-DLBCL': [
                    'Lymphoid Neoplasm Diffuse Large B-cell Lymphoma'
                ],
                'CPTAC-3': ['Adenomas and Adenocarcinomas'],
                'BEATAML1.0-CRENOLANIB': ['Myeloid Leukemias'],
                'CTSP-DLBCL1': [
                    'Lymphoid Neoplasm Diffuse Large B-cell Lymphoma'
                ],
            }
        ),
        (
            'annotations',
            'annotation_id',
            'project.program.name',
            [
                '6e5a8aeb-9c98-4f82-bd5b-08ce6f764485',
                'invalid-value',
            ],
            {
                '6e5a8aeb-9c98-4f82-bd5b-08ce6f764485': ['TARGET'],
                'invalid-value': [],
            }
        ),
    ]
)
def test_map_two_fields(
    endpoint,
    input_field,
    output_field,
    input_values,
    expect
):
    actual = gdc.map_two_fields(
        endpoint, input_field, output_field, input_values
    )
    assert actual == expect
