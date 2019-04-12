import argparse

from .utils import equal_matrices


def main():
    """
    Program entry point
    """
    parser = create_parser()
    options = vars(parser.parse_args())
    # handle check xena equality matrices
    if 'df1' in options and 'df2' in options:
        equal_matrices(options['df1'], options['df2'])


def create_parser():
    """
    Construct the program options
    """
    parser = argparse.ArgumentParser(
        prog="xge",
        description="Extract, transform and load GDC data onto UCSC Xena"
    )
    subparsers = parser.add_subparsers(help="Sub-parsers for xena-gdc-ETL")
    equality_parser = subparsers.add_parser(
        "xena-eql",
        help="Test the equality of 2 Xena matrices."
    )
    equality_parser.add_argument(
        'df1', type=str,
        help='Directory for the first matrix.'
    )
    equality_parser.add_argument(
        'df2', type=str,
        help='Directory for the second matrix.'
    )
    return parser
