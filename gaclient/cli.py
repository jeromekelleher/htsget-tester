"""
Command line interface for the GA4GH read data API.
"""
from __future__ import division
from __future__ import print_function

import argparse

import gaclient


def run(args):
    client = gaclient.Client(
        base_url=args.url,
        id_=args.id,
        reference_name=args.reference_name,
        start=args.start,
        end=args.end,
        output_file=args.output)

    client.download()
    client.close()


def get_gaclient_parser():
    parser = argparse.ArgumentParser(
        description="Command line interface for the GA4GH Streaming API.")
    parser.add_argument(
        "-V", "--version", action='version',
        version='%(prog)s {}'.format(gaclient.__version__))

    parser.add_argument(
        "url", type=str, help="The URL prefix of the server")
    parser.add_argument(
        "id", type=str, help="The ID of the ReadGroupSet")
    parser.add_argument(
        "--format", "-F", type=str, default=None,
        help="The format of data to request.")
    parser.add_argument(
        "--reference-name", "-r", type=str, default=None,
        help="The reference name. If not specified return all reads")
    parser.add_argument(
        "--start", "-s", type=int, default=None,
        help="The starting coordinate.")
    parser.add_argument(
        "--end", "-e", type=int, default=None,
        help="The end coordinate")
    parser.add_argument(
        "--output", "-O", type=str, default=None,
        help="The output file path. Defaults to stdout")
    return parser


def main():
    parser = get_gaclient_parser()
    args = parser.parse_args()
    run(args)