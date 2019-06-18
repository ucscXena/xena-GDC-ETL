import unittest
import pkg_resources

import pytest

from xena_gdc_etl import main


class ParserTest(unittest.TestCase):
    def setUp(self):
        self.parser = main.create_parser()

    @pytest.fixture(autouse=True)
    def _capture(self, capfd):
        self.capfd = capfd

    def test_xena_eql(self):
        parsed = self.parser.parse_args(["xena-eql", "df1", "df2"])
        assert parsed.subcomm == "xena-eql"
        assert parsed.df1 == "df1"
        assert parsed.df2 == "df2"

    def test_gdc_check_new(self):
        parsed = self.parser.parse_args(
            ["gdc-check-new", "https://example.com/data.gz"]
        )
        assert parsed.subcomm == "gdc-check-new"
        assert parsed.url == "https://example.com/data.gz"

    def test_merge_xena(self):
        parsed = self.parser.parse_args(
            [
                "merge-xena",
                "-f",
                "path/to/matrix1",
                "path/to/matrix2",
                "-t",
                "datatype",
                "-o",
                "path/to/dir",
                "-n",
                "new_name",
                "-c",
                "cohort_name",
            ]
        )
        assert parsed.subcomm == "merge-xena"
        assert parsed.files == ["path/to/matrix1", "path/to/matrix2"]
        assert parsed.datatype == "datatype"
        assert parsed.outdir == "path/to/dir"
        assert parsed.name == "new_name"
        assert parsed.cohort == "cohort_name"

    def test_etl(self):
        parsed = self.parser.parse_args(
            [
                "etl",
                "-r",
                "path/to/dir",
                "-p",
                "project_name",
                "-t",
                "datatype",
                "-D",
            ]
        )
        assert parsed.subcomm == "etl"
        assert parsed.root == "path/to/dir"
        assert parsed.projects == ["project_name"]
        assert parsed.datatype == ["datatype"]
        assert parsed.delete is True
        # for mutually exclusive groups
        parsed = self.parser.parse_args(
            [
                "etl",
                "-r",
                "path/to/dir",
                "-P",
                "not_this_project_name",
                "-T",
                "not_this_datatype",
            ]
        )
        assert parsed.subcomm == "etl"
        assert parsed.root == "path/to/dir"
        assert parsed.not_projects == ["not_this_project_name"]
        assert parsed.not_datatype == ["not_this_datatype"]

    def test_metaparser(self):
        parsed = self.parser.parse_args(
            [
                "metadata",
                "-p",
                "project_name",
                "-t",
                "datatype",
                "-m",
                "path/to/matrix",
                "-r",
                "10",
            ]
        )
        assert parsed.subcomm == "metadata"
        assert parsed.project == "project_name"
        assert parsed.datatype == "datatype"
        assert parsed.matrix == "path/to/matrix"
        assert parsed.release == 10.0

    def test_version(self):
        with pytest.raises(SystemExit):
            self.parser.parse_args(['--version'])
            out, _ = self.capfd.readouterr()
            __version__ = pkg_resources.get_distribution(
                "xena_gdc_etl"
            ).version
            assert out == "xge " + __version__ + "\n"
