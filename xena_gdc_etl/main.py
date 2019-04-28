import argparse

from .utils import equal_matrices, metadata
from .gdc_check_new import gdc_check_new
from .constants import valid_dtype


def main():
    """
    Program entry point
    """
    parser = create_parser()
    options = vars(parser.parse_args())
    # handle check xena equality matrices
    if 'df1' in options and 'df2' in options:
        equal_matrices(options['df1'], options['df2'])
    # handle make metadata
    elif "matrix" in options and "datatype" in options:
        metadata(options["matrix"], options["datatype"])
    # handle gdc_check_new
    elif "url" in options:
        gdc_check_new(options["url"])


def create_parser():
    """
    Construct the program options
    """
    parser = argparse.ArgumentParser(
        prog="xge",
        description="Extract, transform and load GDC data onto UCSC Xena"
    )
    subparsers = parser.add_subparsers(help="Sub-parsers for xena-gdc-ETL")
    # equal_matrices subparser
    equality_parser = subparsers.add_parser(
        "xena-eql",
        help="Test the equality of 2 Xena matrices."
    )
    equality_parser.add_argument(
        "df1", type=str,
        help='Directory for the first matrix.'
    )
    equality_parser.add_argument(
        "df2", type=str,
        help='Directory for the second matrix.'
    )
    # make_metadata subparser
    make_metadata_parser = subparsers.add_parser(
        "make-metadata",
        help="Generating Xena metadata for a Xena matrix."
    )
    make_metadata_parser.add_argument(
        '-m', '--matrix', type=str, required=True,
        help='The path, including the filename, of the Xena matrix.'
    )
    make_metadata_parser.add_argument(
        '-t', '--datatype', type=str, required=True,
        help='One data type code indication the data type in matrices to be '
        'merged. Supported data type codes include: {}'.format(
            str(valid_dtype)
        )
    )
    # gdc_check_new subparser
    gdc_check_new_parser = subparsers.add_parser(
        "gdc_check_new",
        description="Check GDC's list of updated files and summarize "
                    "impacted project(s), data_type(s) and "
                    "analysis.workflow_type(s)."
    )
    gdc_check_new_parser.add_argument(
        'url', type=str, metavar='URL',
        help='URL for GDC\'s list of updated files. It can be a compressed '
             'file with a supported extension, which includes ".gz", ".bz2", '
             '".zip", or "xz". New files should be listed under a column named'
             ' by "New File UUID".')
    return parser
