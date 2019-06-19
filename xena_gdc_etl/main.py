from __future__ import print_function
import argparse
from datetime import date
import os
import pkg_resources

import pandas as pd
from pandas.util.testing import assert_frame_equal

from .utils import handle_merge_xena
from .gdc import gdc_check_new, get_project_info
from .constants import valid_dtype, GDC_RELEASE_URL
from .xena_dataset import GDCOmicset, GDCPhenoset, GDCSurvivalset
from .gdc2xena import gdc2xena


__version__ = pkg_resources.get_distribution("xena_gdc_etl").version


def main():
    """
    Program entry point
    """
    parser = create_parser()
    options = parser.parse_args()
    # handle check xena equality matrices
    if options.subcomm == "xena-eql":
        df1 = (
            pd.read_csv(options.df1, sep='\t', header=0, index_col=0)
            .sort_index(axis=0)
            .sort_index(axis=1)
        )
        df2 = (
            pd.read_csv(options.df2, sep='\t', header=0, index_col=0)
            .sort_index(axis=0)
            .sort_index(axis=1)
        )
        try:
            assert_frame_equal(df1, df2)
            print('Equal.')
        # apparently AssertionError doesn't catch all
        except:  # noqa: E722
            print('Not equal.')
    # handle gdc_check_new
    elif options.subcomm == "gdc-check-new":
        new_file_uuids = pd.read_csv(options.url, sep='\t')[
            'New File UUID'
        ].tolist()
        gdc_check_new(new_file_uuids)
    # handle merge_xena
    elif options.subcomm == "merge-xena":
        handle_merge_xena(
            options.name,
            options.files,
            options.cohort,
            options.datatype,
            options.outdir,
        )
    # handle etl
    elif options.subcomm == 'etl':
        delete_raw_data = options.delete
        root_dir = os.path.abspath(options.root)
        projects = options.projects
        if 'all' in [p.lower() for p in projects]:
            projects = [str(x) for x in get_project_info().index]
        for p in options.not_projects:
            projects.remove(p)
        xena_dtypes = options.datatype
        if 'all' in [t.lower() for t in xena_dtypes]:
            xena_dtypes = valid_dtype
        for t in options.not_datatype:
            xena_dtypes.remove(t)
        print('#### GDC to Xena Importing Settings ####')
        total_projects = len(projects)
        print('Import the following {} projects:'.format(total_projects))
        print(repr(projects), end='\n\n')
        print('for the following {} types of data:'.format(len(xena_dtypes)))
        print(str(xena_dtypes), end='\n\n')
        print('into this directory: {}'.format(root_dir))
        print('########################################', end='\n\n')
        gdc2xena(root_dir, projects, xena_dtypes, delete_raw_data)
    # handle metadata
    elif options.subcomm == 'metadata':
        root_dir = os.path.dirname(options.matrix)
        if options.datatype == 'survival':
            dataset = GDCSurvivalset(options.project, root_dir)
        elif options.datatype == 'raw_phenotype':
            if options.project.startswith('TCGA'):
                dataset = GDCPhenoset(
                    options.project, 'raw_phenotype', root_dir
                )
            if options.project.startswith('TARGET'):
                dataset = GDCPhenoset(options.project, 'clinical', root_dir)
        elif options.datatype == 'GDC_phenotype':
            dataset = GDCPhenoset(options.project, 'GDC_phenotype', root_dir)
        else:
            dataset = GDCOmicset(options.project, options.datatype, root_dir)
        dataset.matrix = options.matrix
        dataset.gdc_release = (
            GDC_RELEASE_URL
            + '#data-release-'
            + str(options.release).replace('.', '')
        )
        dataset.metadata()


def create_parser():
    """
    Construct the program options
    """
    parser = argparse.ArgumentParser(
        prog="xge",
        description="Extract, transform and load GDC data onto UCSC Xena",
    )
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s {v}".format(v=__version__),
    )
    subparsers = parser.add_subparsers(
        help="Sub-parsers for xena-gdc-ETL", dest="subcomm"
    )
    # equal_matrices subparser
    equality_parser = subparsers.add_parser(
        "xena-eql", help="Test the equality of 2 Xena matrices."
    )
    equality_parser.add_argument(
        "df1", type=str, help='Directory for the first matrix.'
    )
    equality_parser.add_argument(
        "df2", type=str, help='Directory for the second matrix.'
    )
    # gdc_check_new subparser
    gdc_check_new_parser = subparsers.add_parser(
        "gdc-check-new",
        description="Check GDC's list of updated files and summarize "
        "impacted project(s), data_type(s) and "
        "analysis.workflow_type(s).",
    )
    gdc_check_new_parser.add_argument(
        'url',
        type=str,
        metavar='URL',
        help='URL for GDC\'s list of updated files. It can be a compressed '
        'file with a supported extension, which includes ".gz", ".bz2", '
        '".zip", or "xz". New files should be listed under a column named'
        ' by "New File UUID".',
    )
    # merge-xena subparser
    merge_xena_subparser = subparsers.add_parser(
        "merge-xena",
        description='Pipeline for merging Xena matrices of the same data'
        'type.',
    )
    merge_xena_subparser.add_argument(
        '-f',
        '--files',
        type=str,
        nargs='+',
        required=True,
        help='A list of paths for Xena matrices files to be merged. All paths '
        'in this list support UNIX style pathname pattern expansion with '
        '"glob". Files will be read by pandas.read_csv with sep="\t".',
    )
    merge_xena_subparser.add_argument(
        '-t',
        '--datatype',
        type=str,
        required=True,
        help='One data type code indication the data type in matrices to be '
        'merged. Supported data type codes include: {}'.format(
            str(valid_dtype)
        ),
    )
    merge_xena_subparser.add_argument(
        '-o',
        '--outdir',
        type=str,
        default='.',
        help='A directory to put the merged matrix. Defaults to the current '
        'working directory of python.',
    )
    merge_xena_subparser.add_argument(
        '-n',
        '--name',
        type=str,
        default=None,
        help='Filename for the merged matrix. Defaults to None. If None, the '
        'filename will be derived from the cohort name and the data type. '
        'Check "-t" and "-c" options for details.',
    )
    merge_xena_subparser.add_argument(
        '-c',
        '--cohort',
        type=str,
        default=None,
        help='A cohort name for the merged matrix. Defaults to None. If '
        'None, it will be set to a format of "MergedCohort<date>" by default. '
        'For example, "MergedCohort{}".'.format(
            date.today().strftime('%m%d%Y')
        ),
    )
    # Subcommand for full ETL (download, transform, and metadata)
    etlparser = subparsers.add_parser(
        'etl',
        help='Download and transform GDC data into Xena matrix, '
        'and generate corresponding metadata.',
        epilog='Supported data types are: {}'.format(str(valid_dtype)),
    )
    etlparser.add_argument(
        '-r',
        '--root',
        type=str,
        default='.',
        help='Root directory for imported data.',
    )
    etlparser.add_argument(
        '-D',
        '--delete',
        action='store_true',
        help='Deletes raw data upon generation of Xena_matrix.',
    )
    projects_group = etlparser.add_mutually_exclusive_group()
    projects_group.add_argument(
        '-p',
        '--projects',
        type=str,
        nargs='+',
        help='GDC project ID(s) to be imported; or "all" if all projects on'
        'GDC are going to be imported. Defaults to "all".',
        default=['all'],
    )
    projects_group.add_argument(
        '-P',
        '--not-projects',
        type=str,
        nargs='+',
        help='Import all projects on GDC except projects specified by this'
        'option. This option and the "-p" option are mutually exclusive.',
        default=[],
    )
    datatype_group = etlparser.add_mutually_exclusive_group()
    datatype_group.add_argument(
        '-t',
        '--datatype',
        type=str,
        nargs='+',
        help='Data type code(s) to be imported; or "all" if all supported'
        'types are going to be imported. Defaults to "all".',
        default=['all'],
    )
    datatype_group.add_argument(
        '-T',
        '--not-datatype',
        type=str,
        nargs='+',
        help='Import all supported types except projects specified by this'
        'option. This option and the "-t" option are mutually exclusive.',
        default=[],
    )
    # Subcommand for making metadata
    metaparser = subparsers.add_parser(
        'metadata',
        help='Generate metadata for a Xena matrix',
        epilog='Supported data types are: {}'.format(str(valid_dtype)),
    )
    metaparser.add_argument(
        '-p',
        '--project',
        type=str,
        required=True,
        help='The project of the matrix.',
    )
    metaparser.add_argument(
        '-t',
        '--datatype',
        type=str,
        required=True,
        help='One data type code for the matrix.',
    )
    metaparser.add_argument(
        '-m', '--matrix', type=str, required=True, help='Path to a Xena matrix'
    )
    metaparser.add_argument(
        '-r',
        '--release',
        type=float,
        required=True,
        help='GDC data release number.',
    )
    return parser
