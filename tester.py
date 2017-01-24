"""
Data tester for the GA4GH streaming API.
"""
from __future__ import print_function
from __future__ import division

import argparse
import collections
import logging
import os
import os.path
import subprocess
import random
import tempfile
import time

import htsget
import pysam
import humanize

from six.moves import zip


__version__ = "0.1.0"


def htsget_api(url, filename, reference_name=None, start=None, end=None):
    logging.info("htsget-api: request ('{}', {}, {})".format(reference_name, start, end))
    with open(filename, "wb") as tmp:
        htsget.get(url, tmp, reference_name=reference_name, start=start, end=end)


def htsget_cli(url, filename, reference_name=None, start=None, end=None):
    cmd = ["htsget", url, "-O", filename]
    if reference_name is not None:
        cmd.extend(["-r", str(reference_name)])
    if start is not None:
        cmd.extend(["-s", str(start)])
    if end is not None:
        cmd.extend(["-e", str(end)])
    logging.info("htsget-cli: run {}".format(" ".join(cmd)))
    subprocess.check_call(cmd)


class TestFailedException(Exception):
    """
    Exception raised when we know we've failed.
    """


def hashkey(read):
    """
    Returns a unique hash key for the specified read. We cannot use the qname
    as mates and secondary alignments have the same qname.
    # (first of pair, second of pair, not primary alignment, supplementary alignment).
    """
    key = "{}_{}_{}_{}_{}_{}_{}".format(
        read.query_name, read.pos, read.cigar,
        int(read.is_read1), int(read.is_read2),
        int(read.is_secondary), int(read.is_supplementary))
    return key


def filter_unmapped(iterator):
    """
    Returns an iterator over the specified iterator of reads with
    all unmapped reads removed.
    """
    filtered = 0
    for read in iterator:
        if read.is_unmapped:
            filtered += 1
        else:
            yield read
    logging.info("Filtered {} reads".format(filtered))


class Contig(object):
    """
    Represents a single contig within the BAM file.
    """
    def __init__(
            self, reference_name, length, start_positions, end_positions):
        self.reference_name = reference_name
        self.length = length
        self.start_positions = start_positions
        self.end_positions = end_positions


class Field(object):
    """
    A class to simplify comparing a single field in two different Pysam Read
    records.
    """
    def __init__(self, name, equals):
        self.name = name
        self.equals = equals


def simple_equality(v1, v2):
    """
    Equality for simple fields where equality can be tested directly.
    """
    return v1 == v2


def sorted_equality(v1, v2):
    """
    Equality for fields where the values must be sorted before equality tested.
    """
    return sorted(v1) == sorted(v2)


class ServerTester(object):
    """
    Class to run systematic tests on a server based on a local indexed
    BAM file.
    """
    fields = [
        Field("pos", simple_equality),
        Field("query_name", simple_equality),
        Field("reference_name", simple_equality),
        Field("cigarstring", simple_equality),
        Field("query_alignment_sequence", simple_equality),
        Field("query_alignment_qualities", simple_equality),
        Field("template_length", simple_equality),
        Field("next_reference_id", simple_equality),
        Field("next_reference_start", simple_equality),
        Field("flag", simple_equality),
        Field("mapping_quality", simple_equality),
        Field("tags", sorted_equality),
        # TODO fill in remaining BAM fields.
    ]

    def __init__(
            self, source_file_name, server_url, filter_unmapped=False, tmpdir=None,
            num_boundary_reads=10, max_references=100, max_random_query_length=10**6,
            client=None):
        self.source_file_name = source_file_name
        self.server_url = server_url
        self.filter_unmapped = filter_unmapped
        self.tmpdir = tmpdir
        self.max_random_query_length = max_random_query_length
        self.alignment_file = pysam.AlignmentFile(self.source_file_name)
        self.temp_file_name = None
        self.num_boundary_reads = num_boundary_reads
        self.max_references = max_references
        self.client = client
        self.contigs = []
        self.mismatch_counts = collections.Counter()
        self.total_queries = 0
        self.total_downloaded_data = 0
        self.total_download_time = 1e-8  # Avoid zero division problems

    def get_start_positions(self, reference_name):
        """
        Returns the positions of the first num_boundary_reads on the specified reference.
        """
        start_positions = []
        for read in self.alignment_file.fetch(reference_name):
            start_positions.append(read.pos)
            if len(start_positions) == self.num_boundary_reads:
                break
        return start_positions

    def get_end_positions(self, reference_name, length):
        """
        Returns the positions of the last num_boundary_reads on the specified reference.
        """
        found = False
        x = length
        while not found:
            for read in self.alignment_file.fetch(reference_name, x):
                found = True
                break
            if not found:
                # Skip back 1% of the length of the chromosome
                x = x - length / 100
        positions = collections.deque([], maxlen=self.num_boundary_reads)
        for read in self.alignment_file.fetch(reference_name, x):
            positions.append(read.pos)
        return list(positions)

    def initialise(self):
        """
        Scans the input BAM file and initialises the data structures we need for the
        tests.
        """
        fd, self.temp_file_name = tempfile.mkstemp(
            prefix="gastream_test_", dir=self.tmpdir, suffix=".bam")
        os.close(fd)
        # Determine the bounds of the individual contigs.
        total_references = len(self.alignment_file.lengths)
        num_references = min(self.max_references, total_references)
        logging.info("Reading file {}".format(self.source_file_name))
        for j in range(num_references):
            reference_name = self.alignment_file.references[j]
            length = self.alignment_file.lengths[j]
            start_positions = self.get_start_positions(reference_name)
            if len(start_positions) > 0:
                end_positions = self.get_end_positions(reference_name, length)
                contig = Contig(reference_name, length, start_positions, end_positions)
                self.contigs.append(contig)
                msg = (
                    "Read contig {}: got {} start positions and {} end "
                    "positions".format(
                        contig.reference_name, len(start_positions), len(end_positions)))
                logging.info(msg)
            else:
                logging.info("Skipping empty contig {}".format(reference_name))

    def verify_reads_equal(self, r1, r2):
        for field in self.fields:
            v1 = getattr(r1, field.name)
            v2 = getattr(r2, field.name)
            if not field.equals(v1, v2):
                logging.info("Mismatch at {}:{}.{}:: {} != {}".format(
                    r1.reference_name, r1.pos, field.name, v1, v2))
                self.mismatch_counts[field.name] += 1

    def verify_reads(self, iter1, iter2):
        """
        Verifies that the specified iterators contain the same set of
        reads. Returns the number of reads in the iterators.
        """
        d1 = {}
        d2 = {}
        last_pos = -1
        num_reads = 0
        num_checks = 0
        for r1, r2 in zip(iter1, iter2):
            num_reads += 1
            if r1.pos != last_pos:
                assert len(d1) == len(d2)
                for k in d1.keys():
                    # TODO change these asserts to raise a TestFailedException
                    assert k in d2
                    self.verify_reads_equal(d1[k], d2[k])
                    num_checks += 1
                d1 = {}
                d2 = {}
                last_pos = r1.pos
            k = hashkey(r1)
            assert k not in d1
            d1[k] = r1
            k = hashkey(r2)
            assert k not in d2
            d2[k] = r2
        # Check the last set of reads in the dictionaries.
        assert len(d1) == len(d2)
        for k in d1.keys():
            assert k in d2
            self.verify_reads_equal(d1[k], d2[k])
            num_checks += 1
        assert num_checks == num_reads
        return num_reads

    def verify_query(self, reference_name, start=None, end=None):
        """
        Runs the specified query and verifies the result.
        """
        self.total_queries += 1
        # We use the wall-clock time here
        before = time.time()
        self.client(
            self.server_url, self.temp_file_name, reference_name=reference_name,
            start=start, end=end)
        duration = time.time() - before
        size = os.path.getsize(self.temp_file_name)
        self.total_downloaded_data += size
        self.total_download_time += duration
        logging.info("Downloaded {} in {:.3f} seconds".format(
            humanize.naturalsize(size, binary=True), duration))
        # Index the downloaded file and compare the reads.
        before = time.clock()
        pysam.index(self.temp_file_name)
        duration = time.clock() - before
        logging.debug("Indexed BAM file in {:.3f} CPU seconds".format(duration))

        # Analyse the reads
        before = time.clock()
        subset_reads = pysam.AlignmentFile(self.temp_file_name)
        iter1 = self.alignment_file.fetch(reference_name, start, end)
        try:
            iter2 = subset_reads.fetch(reference_name, start, end)
        except ValueError as ve:
            raise TestFailedException("Reading downloaded data: {}".format(ve))
        if self.filter_unmapped:
            iter1 = filter_unmapped(iter1)
        num_reads = self.verify_reads(iter1, iter2)
        # Now count the reads outside the original region.
        subset_reads.reset()
        all_reads = 0
        for read in subset_reads:
            all_reads += 1
        extra = all_reads - num_reads
        duration = time.clock() - before
        logging.debug("Checked reads in {:.3f} CPU seconds".format(duration))
        logging.info("Downloaded {} reads with {} extra".format(num_reads, extra))

    def run_random_reads(self, num_tests):
        """
        Runs tests for valid ranges chosen randomly across the available
        regions.
        """
        logging.info("Starting random queries")
        for _ in range(num_tests):
            contig = random.choice(self.contigs)
            start = random.randint(0, contig.length)
            end = random.randint(
                start, min(start + self.max_random_query_length, contig.length))
            self.verify_query(contig.reference_name, start, end)

    def run_full_contig_fetch(self):
        """
        Gets all reads for contigs < 1Mb
        """
        logging.info("Starting full contig queries")
        for contig in self.contigs:
            if contig.length < 10**6:
                self.verify_query(contig.reference_name)
            else:
                logging.info("Skipping full fetch for {}; too long".format(
                    contig.reference_name))

    def run_start_reads(self):
        """
        Gets the first few reads from each contig.
        """
        logging.info("Starting contig start queries")
        for contig in self.contigs:
            self.verify_query(
                contig.reference_name,
                contig.start_positions[0],
                contig.start_positions[-1] + 1)
            self.verify_query(
                contig.reference_name,
                None,
                contig.start_positions[-1] + 1)
            self.verify_query(
                contig.reference_name,
                max(0, contig.start_positions[0] - 100),
                contig.start_positions[0] + 1)
            self.verify_query(
                contig.reference_name,
                contig.start_positions[-1],
                contig.start_positions[-1] + 1)

    def run_end_reads(self):
        """
        Gets the first few reads from each contig.
        """
        logging.info("Starting contig end queries")
        for contig in self.contigs:
            if len(contig.end_positions) > 0:
                self.verify_query(
                    contig.reference_name,
                    contig.end_positions[0],
                    contig.end_positions[-1] + 1)
                self.verify_query(
                    contig.reference_name,
                    contig.end_positions[0],
                    None)
                self.verify_query(
                    contig.reference_name,
                    max(0, contig.end_positions[0] - 100),
                    contig.end_positions[0] + 1)
            else:
                logging.info("Skipping end reads for {}".format(contig.reference_name))

    def cleanup(self):
        self.alignment_file.close()
        if os.path.exists(self.temp_file_name):
            os.unlink(self.temp_file_name)
            index_file = self.temp_file_name + ".bai"
            if os.path.exists(index_file):
                os.unlink(index_file)

    def report(self):
        # TODO provide options to make this machine readable.
        megabytes = self.total_downloaded_data / (1024**2)
        average_bandwidth = megabytes / self.total_download_time
        print("{} queries downloaded {} of data at an average of {:.2f} Mib/s".format(
            humanize.intword(self.total_queries),
            humanize.naturalsize(self.total_downloaded_data, binary=True),
            average_bandwidth))
        print("Total mismatches = ", sum(self.mismatch_counts.values()))
        for name in sorted(self.mismatch_counts.keys()):
            print("", name, self.mismatch_counts[name], sep="\t")


if __name__ == "__main__":

    client_map = {
        "htsget-api": htsget_api,
        "htsget-cli": htsget_cli,
    }

    parser = argparse.ArgumentParser(
        description="A simple tester application for the GA4GH streaming API.")
    parser.add_argument(
        "-V", "--version", action='version',
        version='%(prog)s {}'.format(__version__))
    parser.add_argument('--verbose', '-v', action='count', default=0)
    parser.add_argument(
        "bam_file", type=str, help="The local BAM file to compare to")
    parser.add_argument(
        "url", type=str, help="The server url corresponding to the input BAM file.")
    parser.add_argument(
        "--filter-unmapped", action='store_true',
        help="Filter out unmapped reads before comparing reads.")
    parser.add_argument(
        "--tmpdir", type=str,
        help=(
            "Directory in which to store temporary files. "
            "Defaults to platform default"))
    parser.add_argument(
        "--max-references", type=int, default=100,
        help="The maximum number of references to consider")
    parser.add_argument(
        "--num-random-reads", type=int, default=20,
        help="The number of random queries to send")
    parser.add_argument(
        "--random-seed", type=int, default=1,
        help=(
            "The random seed. Use this to ensure that the same set of queries is "
            "run across different client/server combinations."))
    parser.add_argument(
        "--max-random-query-length", type=int, default=10**7,
        help="The maximum length of a random query in bases")
    parser.add_argument(
        "--client", choices=list(client_map.keys()), default="htsget-api",
        help="The client to use for running the transfer")

    args = parser.parse_args()
    log_level = logging.WARNING
    if args.verbose == 1:
        log_level = logging.INFO
    if args.verbose >= 2:
        log_level = logging.DEBUG
    logging.basicConfig(format='%(asctime)s %(message)s', level=log_level)
    tester = ServerTester(
        args.bam_file, args.url, filter_unmapped=args.filter_unmapped,
        tmpdir=args.tmpdir, max_references=args.max_references,
        client=client_map[args.client])
    try:
        tester.initialise()
        tester.run_full_contig_fetch()
        tester.run_start_reads()
        tester.run_end_reads()
        tester.run_random_reads(args.num_random_reads)
    except KeyboardInterrupt:
        logging.warn("Interrupted!")
    finally:
        tester.cleanup()
    tester.report()
