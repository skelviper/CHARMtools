#!/usr/bin/env python3
"""
CHARMtools main entry point
This script serves as the main entry point for CHARMtools package
to maintain compatibility with existing pipelines.
"""

import sys
import os
import argparse

# Add the current directory to Python path to ensure imports work
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Import only the preprocessing modules that are actually used in the pipeline
from charm_preprocess import clean_leg
from charm_preprocess import clean_splicing
from charm_preprocess import clean_isolated
from charm_preprocess import sep_clean
from charm_preprocess import clean3

def cli():
    parser = argparse.ArgumentParser(prog="CHARMtools", description="Functions for data-analysis in HiRES etc. projects")
    subcommands = parser.add_subparsers(title="These are sub-commands", metavar="command")
    
    # clean_leg sub command
    clean_leg_arg = subcommands.add_parser(
        "clean_leg",
        help="clean promiscuous legs that contacts with multiple legs")
    clean_leg_arg.set_defaults(handle=clean_leg.cli)
    clean_leg_arg.add_argument(
        dest="filename",
        metavar="INPUT_FILE",
        nargs=1)
    clean_leg_arg.add_argument(
        "-t", "--thread",
        type=int,
        dest="thread",
        action="store",
        default=4,
        help="set thread number")
    clean_leg_arg.add_argument(
        "-d", "--distance",
        dest="max_distance",
        metavar="MAX_DISTANCE",
        type=int,
        action="store",
        default=1000,
        help="max distance to calculate adjacent legs")
    clean_leg_arg.add_argument(
        "-n", "--count",
        metavar="MAX_COUNT",
        dest="max_count",
        type=int,
        action="store",
        default=15,
        help="number threshold of adjacent legs")
    clean_leg_arg.add_argument(
        "-o", "--output",
        dest="out_name",
        action="store",
        metavar="OUTPUT_FILE",
        help="output file name")
    
    # clean_isolated sub command
    clean_isolated_arg = subcommands.add_parser(
        "clean_isolated",
        help="clean isolated contact points")
    clean_isolated_arg.set_defaults(handle=clean_isolated.cli)
    clean_isolated_arg.add_argument(
        dest="filename",
        metavar="INPUT_FILE",
        nargs=1)
    clean_isolated_arg.add_argument(
        "-t", "--thread",
        type=int,
        dest="thread",
        action="store",
        default=4,
        help="set thread number")
    clean_isolated_arg.add_argument(
        "-o", "--output",
        dest="output_file",
        action="store",
        metavar="OUTPUT_FILE",
        help="output file name")
    clean_isolated_arg.add_argument(
        "--dense",
        type=int,
        dest="dense",
        action="store",
        default=5,
        help="dense threshold")
    clean_isolated_arg.add_argument(
        "--distance",
        type=int,
        dest="distance",
        action="store",
        default=10000000,
        help="distance threshold")
    
    # clean_splicing sub command
    clean_splicing_arg = subcommands.add_parser(
        "clean_splicing",
        help="clean exon splicing from mRNA in contact files")
    clean_splicing_arg.set_defaults(handle=clean_splicing.cli)
    clean_splicing_arg.add_argument(
        dest="filename",
        metavar="INPUT_FILE",
        nargs=1)
    clean_splicing_arg.add_argument(
        "-r", "--reference",
        dest="gtf_filename",
        action="store",
        metavar="ANNOTATION_FILE",
        help="annotation file")
    clean_splicing_arg.add_argument(
        "-t", "--thread",
        type=int,
        dest="num_thread",
        action="store",
        default=4,
        help="set thread number")
    clean_splicing_arg.add_argument(
        "-o", "--output",
        dest="out_name",
        action="store",
        metavar="OUTPUT_FILE",
        help="output file name")
    
    # sep_clean sub command
    sep_clean_arg = subcommands.add_parser(
        "sep_clean",
        help="separation cleaning tools")
    sep_clean_arg.set_defaults(handle=sep_clean.cli)
    sep_clean_arg.add_argument(
        dest="filename",
        metavar="INPUT_FILE",
        nargs=1)
    sep_clean_arg.add_argument(
        "-n", "--thread",
        type=int,
        dest="num_thread",
        action="store",
        default=4,
        help="set thread number")
    sep_clean_arg.add_argument(
        "-o1", "--output1",
        dest="output_file1",
        action="store",
        metavar="OUTPUT_FILE1",
        help="output file name 1")
    sep_clean_arg.add_argument(
        "-o2", "--output2",
        dest="output_file2",
        action="store",
        metavar="OUTPUT_FILE2",
        help="output file name 2")
    sep_clean_arg.add_argument(
        "-d", "--dense",
        dest="dense",
        type=int,
        action="store",
        default=5,
        help="dense threshold")
    sep_clean_arg.add_argument(
        "--distance",
        dest="distance",
        type=int,
        action="store",
        default=1000000,
        help="distance threshold")
    
    # clean3 sub command
    clean3_arg = subcommands.add_parser(
        "clean3",
        help="3D structure cleaning")
    clean3_arg.set_defaults(handle=clean3.cli)
    clean3_arg.add_argument(
        "-i", "--input",
        dest="filename",
        action="store",
        metavar="INPUT_FILE",
        help="input 3dg file")
    clean3_arg.add_argument(
        "-r", "--reference",
        dest="ref_filename",
        action="store",
        metavar="REFERENCE_FILE",
        help="reference pairs file")
    clean3_arg.add_argument(
        "-o", "--output",
        dest="output",
        action="store",
        metavar="OUTPUT_FILE",
        help="output file name")
    clean3_arg.add_argument(
        "-q", "--quantile",
        dest="quantile",
        type=float,
        action="store",
        default=0.1,
        help="clean quantile threshold")
    clean3_arg.add_argument(
        "-d", "--distance",
        dest="distance",
        type=float,
        action="store",
        default=2000000.0,
        help="max clean distance")
    
    # Parse arguments and call appropriate handler
    args = parser.parse_args()
    if hasattr(args, 'handle'):
        args.handle(args)
    else:
        parser.print_help()

if __name__ == "__main__":
    cli()