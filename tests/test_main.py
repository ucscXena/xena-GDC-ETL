import unittest

from xena_gdc_etl import main


class ParserTest(unittest.TestCase):
    def setUp(self):
        self.parser = main.create_parser()

    def test_xena_eql(self):
        parsed = self.parser.parse_args(["xena-eql", "df1", "df2"])
        self.assertEqual(parsed.df1, "df1")
        self.assertEqual(parsed.df2, "df2")
