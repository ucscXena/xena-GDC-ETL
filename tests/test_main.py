import unittest

from xena_gdc_etl import main


class ParserTest(unittest.TestCase):
    def setUp(self):
        self.parser = main.create_parser()

    def test_xena_eql(self):
        parsed = self.parser.parse_args(["xena-eql", "df1", "df2"])
        assert parsed.df1 == "df1"
        assert parsed.df2 == "df2"

    def test_make_metadata(self):
        parsed = self.parser.parse_args(["make-metadata", "-m",
                                        "path/to/matrix", "-t", "datatype"])
        assert parsed.matrix == "path/to/matrix"
        assert parsed.datatype == "datatype"

    def test_gdc_check_new(self):
        parsed = self.parser.parse_args(["gdc_check_new",
                                        "https://example.com/data.gz"])
        assert parsed.url == "https://example.com/data.gz"
