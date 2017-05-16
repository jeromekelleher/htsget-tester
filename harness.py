"""
Runs tests specified in a input file in parallel.
"""
from __future__ import print_function
from __future__ import division

import os.path
import argparse
import random
import sys
import logging
import multiprocessing

import tqdm
import hjson
import humanize
import pandas as pd

import tester


class TestCase(object):
    """
    A single client/server combination for a data file.
    """
    def __init__(self, input_file, server_name, url, client):
        self.input_file = input_file
        self.server_name = server_name
        self.url = url
        self.client = client
        self.random_seed = 1
        self.max_references = 100
        self.num_random_queries = 5
        self.tmpdir = None
        self.filter_unmapped = False
        self.log_file_prefix = None

    @property
    def log_file_name(self):
        return os.path.join(self.log_file_prefix, "{}.log".format(self.test_name))

    @property
    def test_name(self):
        data_name = os.path.split(self.input_file)[-1]
        return "{}_{}_{}".format(data_name, self.server_name, self.client)


def parse_input(json_cases, data_file_prefix, clients):
    """
    Parses the input JSON into test cases.
    """
    cases = []
    for json_case in json_cases["cases"]:
        input_file = os.path.join(data_file_prefix, json_case["file"])
        for json_server in json_case["servers"]:
            for client in clients:
                cases.append(TestCase(
                    input_file, json_server["name"], json_server["url"], client))
                if "filter_unmapped" in json_server:
                    cases[-1].filter_unmapped = json_server["filter_unmapped"]
    return cases


def run_test(case):
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    log_handler = logging.FileHandler(case.log_file_name, mode="w")
    formatter = logging.Formatter('%(asctime)s %(message)s')
    log_handler.setFormatter(formatter)
    logger.addHandler(log_handler)

    random.seed(case.random_seed)
    server_tester = tester.ServerTester(
        case.input_file, case.url, filter_unmapped=case.filter_unmapped,
        tmpdir=case.tmpdir, max_references=case.max_references,
        client=tester.client_map[case.client])
    completed = False
    try:
        server_tester.initialise()
        server_tester.run_full_contig_fetch()
        server_tester.run_start_reads()
        server_tester.run_end_reads()
        server_tester.run_random_reads(case.num_random_queries)
        completed = True
    except tester.DownloadFailedException as dfe:
        print("Download failed: ", dfe, "see ", case.log_file_name, file=sys.stderr)
    except tester.TestFailedException as tfe:
        print(
            "Test failed:", tfe, ". Is the input file correct?", case.log_file_name,
            file=sys.stderr)
    finally:
        server_tester.cleanup()
        logger.removeHandler(log_handler)
        log_handler.close()
    ret = None
    if completed:
        ret = server_tester.report()
    return case, ret


def main():

    parser = argparse.ArgumentParser(description=(
        "A test harness for the GA4GH streaming API to test multiple client/server "
        "combinations in parallel."))
    parser.add_argument(
        "-V", "--version", action='version',
        version=tester.__version__)
    parser.add_argument('--verbose', '-v', action='count', default=0)
    parser.add_argument(
        "input_file", type=str, help=(
            "A file describing the input files and corresponding remote URLs"))
    parser.add_argument(
        "clients", nargs="*", type=str, help=(
            "The clients to run this against. If not specified use all clients. "
            "Choices = {}".format(sorted(tester.client_map.keys()))))
    parser.add_argument(
        "-d", "--data-file-prefix", default="./",
        help="The path prefix to attach to all data files in the input JSON.")
    parser.add_argument(
        "-l", "--log-file-prefix", default="./",
        help="The path prefix for log files.")
    parser.add_argument(
        "-O", "--output-file", default=None,
        help="The file to write the CSV report summary.")
    parser.add_argument(
        "-P", "--processes", default=1, type=int,
        help="The number of worker processes to use.")
    parser.add_argument(
        "--tmpdir", type=str,
        help=(
            "Directory in which to store temporary files. "
            "Defaults to platform default"))
    parser.add_argument(
        "--max-references", type=int, default=100,
        help="The maximum number of references to consider")
    parser.add_argument(
        "--num-random-queries", type=int, default=20,
        help="The number of random queries to send")
    parser.add_argument(
        "--random-seed", type=int, default=1,
        help=(
            "The random seed. Use this to ensure that the same set of queries is "
            "run across different client/server combinations."))

    args = parser.parse_args()
    clients = list(tester.client_map.keys())
    if len(args.clients) > 0:
        for client in args.clients:
            if client not in tester.client_map:
                raise KeyError("Unknown client:{}".format(client))
        clients = args.clients
    with open(args.input_file) as infile:
        cases = parse_input(hjson.load(infile), args.data_file_prefix, clients)
    for case in cases:
        case.log_file_prefix = args.log_file_prefix
        case.max_references = args.max_references
        case.random_seed = args.random_seed
        case.num_random_queries = args.num_random_queries
        case.tmpdir = args.tmpdir

    progress = tqdm.tqdm(total=len(cases))
    failures = []
    results = []

    def store_result(case, result):
        progress.update()
        if result is None:
            failures.append(case.log_file_name)
        else:
            results.append((case, result))

    pool = multiprocessing.Pool(args.processes)
    try:
        for case, result in pool.imap_unordered(run_test, cases):
            store_result(case, result)
    finally:
        pool.terminate()

    results.sort(key=lambda x: x[0].test_name)
    rows = []
    for case, report in results:
        d = {}
        d["server"] = case.server_name
        d["file"] = os.path.split(case.input_file)[-1]
        d["client"] = case.client
        d["total_downloaded_data"] = report.total_downloaded_data
        d["total_queries"] = report.total_queries
        d["num_references"] = report.num_references
        d["total_mismatches"] = sum(count for count in report.mismatch_counts.values())
        rows.append(d)
    df = pd.DataFrame(rows)
    print(df)
    if args.output_file is not None:
        df.to_csv(args.output_file)

    if len(failures) > 0:
        print(
            "WARNING!! {} cases failed. See following log files:".format(len(failures)))
        for log_file_name in failures:
            print("\t", log_file_name)
    if sum(len(report.mismatch_counts) for _, report in results) > 0:
        print("WARNING!! Mismatches present")
        for case, report in results:
            if len(report.mismatch_counts) > 0:
                print(case.test_name)
                for k, v in report.mismatch_counts.items():
                    print("\t", k, v, sep="\t")



if __name__ == "__main__":
    main()
