import unittest

from xena_gdc_etl import gdc2xena


class ParserTest(unittest.TestCase):
    def setUp(self):
        self.parser = gdc2xena.create_parser()

    def test_etl(self):
        parsed = self.parser.parse_args([
            "etl",
            "-r", 
            "path/to/dir",
            "-p",
            "project_name",
            "-t",
            "datatype",
        ])
        assert parsed.subcomm == "etl"
        assert parsed.root == "path/to/dir"
        assert parsed.projects == ["project_name"]
        assert parsed.datatype == ["datatype"]
        # for mutually exclusive groups
        parsed = self.parser.parse_args([
            "etl",
            "-r", 
            "path/to/dir",
            "-P",
            "not_this_project_name",
            "-T",
            "not_this_datatype",
        ])
        assert parsed.subcomm == "etl"
        assert parsed.root == "path/to/dir"
        assert parsed.not_projects == ["not_this_project_name"]
        assert parsed.not_datatype == ["not_this_datatype"]

    def test_metaparser(self):
        parsed = self.parser.parse_args([
            "metadata",
            "-p",
            "project_name",
            "-t",
            "datatype",
            "-m",
            "path/to/matrix",
            "-r",
            "10",
        ])
        assert parsed.subcomm == "metadata"
        assert parsed.project == "project_name"
        assert parsed.datatype == "datatype"
        assert parsed.matrix == "path/to/matrix"
        assert parsed.release == 10.0
