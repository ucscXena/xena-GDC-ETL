import json
import os

from xena_gdc_etl import gdc
from tests.utils import compare_dict


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
                {"content": {"field": "c", "value": ["d"]}, "op": "exclude"}
            ],
            "op": "and"
        }
    actual = gdc.simple_and_filter(in_dict_2, exclude_dict_2)
    compare_dict(expected, actual)


def test_reduce_json_array():
    input_1 = [{
        'a': 'hello',
        'b': [1, 2, 3],
        'c': [10]
    }]
    input_2 = [{'a': 'b'}]
    actual_1 = json.dumps(gdc.reduce_json_array(input_1), sort_keys=True)
    expected_1 = json.dumps({"a": "hello", "b": [1, 2, 3], "c": 10},
                            sort_keys=True)
    assert actual_1 == expected_1
    actual_2 = gdc.reduce_json_array(input_2)
    expected_2 = {'a': 'b'}
    compare_dict(actual_2, expected_2)


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


def test_mkdir_p():
    input_ = "tests/tmp_dir"
    gdc.mkdir_p(input_)
    assert os.path.isdir(input_) is True
    os.rmdir(input_)


def test_download():
    uuid = "53a637ce-8aaf-4cec-b02d-89202bbb0890"
    gdc.download(uuid, download_dir="./tests")
    file_path = "./tests/" + uuid + ".svs"
    assert os.path.isfile(file_path) is True
    os.unlink(file_path)


def test_download_error():
    uuid = "something-invalid-0101"
    gdc.download(uuid)


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
        "project_id": "TCGA-OV"
    }
    actual.equals(expected)


def test_get_samples_clinical():
    project_id = "TCGA-OV"
    actual = gdc.get_samples_clinical(project_id)
    assert actual['case_id'][0] == "71faa2c1-0d5b-4dcc-bdf9-f2405f29907c"


test_download_error()