from __future__ import print_function
from __future__ import division

import pysam
import os
import itertools
import time
import tempfile
import logging
import random
import os.path
import argparse

import gaclient


def decode_virtual_offset(virtual_offset):
    return virtual_offset >> 16, virtual_offset & 0xFFFF


class Contig(object):
    """
    Represents a single contig within the BAM file.
    """
    def __init__(
            self, reference_name, length, initial_positions, final_positions):
        self.reference_name = reference_name
        self.length = length
        self.initial_positions = initial_positions
        self.final_positions = final_positions


class ServerTester(object):
    """
    Class to run systematic tests on a server based on a local indexed
    BAM file.
    """
    def __init__(self, source_file_name, server_url, read_group_set_id):
        self.source_file_name = source_file_name
        self.server_url = server_url
        self.read_group_set_id = read_group_set_id
        self.alignment_file = pysam.AlignmentFile(self.source_file_name)
        fd, self.temp_file_name = tempfile.mkstemp(prefix="gastream_")
        os.close(fd)
        self.num_initial_reads = 10
        # Determine the bounds of the individual contigs.
        self.contigs = []
        num_references = len(self.alignment_file.lengths)
        for j in range(num_references):
            reference_name = self.alignment_file.references[j]
            length = self.alignment_file.lengths[j]
            # Find the offset of the first k reads
            initial_positions = []
            for read in self.alignment_file.fetch(reference_name):
                initial_positions.append(read.pos)
                if len(initial_positions) == self.num_initial_reads:
                    break
            last_positions = []
            if j < num_references - 1:
                iterator = self.alignment_file.fetch(
                    self.alignment_file.references[j + 1])
                ret = next(iterator, None)
                if ret is not None:
                    # Sometimes this works, somtimes it doesn't. Don't know why.
                    position = self.alignment_file.tell()
                    file_offset, block_offset = decode_virtual_offset(position)
                    if block_offset != 0:
                        block_start = file_offset << 16
                        self.alignment_file.seek(block_start)
                        for read in self.alignment_file:
                            if read.reference_name == reference_name:
                                last_positions.append(read.pos)
                            else:
                                break
            else:
                print("GET EOF??")
            contig = Contig(
                reference_name, length, initial_positions, last_positions)
            self.contigs.append(contig)
            print(
                "READ", contig.reference_name, len(initial_positions),
                "initial reads", len(last_positions), "final reads")

    def verify_reads(self, iter1, iter2):
        """
        Verifies that the specified iterators contain the same set of
        reads. Returns the number of reads in the iterators.
        """
        num_reads = 0
        for r1, r2 in itertools.izip(iter1, iter2):
            num_reads += 1
            if r1 != r2:
                print("ERROR")
                print(r1)
                print(r2)
            assert r1 == r2
        assert next(iter1, None) is None
        assert next(iter2, None) is None
        return num_reads

    def verify_query(self, reference_name, start=None, end=None):
        """
        Runs the specified query and verifies the result.
        """
        print("reference_name, start, end = '{}', {}, {}".format(
            reference_name, start, end))
        # TODO refactor client to be a persistent object
        client = gaclient.Client(
            self.server_url, self.read_group_set_id,
            reference_name, start, end, self.temp_file_name)
        before = time.clock()
        client.download()
        client.close()
        duration = time.clock() - before
        print("Downloaded in {} seconds".format(duration))
        # Index the downloaded file and compare the reads.
        pysam.index(self.temp_file_name)
        subset_reads = pysam.AlignmentFile(self.temp_file_name)
        iter1 = self.alignment_file.fetch(reference_name, start, end)
        iter2 = subset_reads.fetch(reference_name, start, end)
        num_reads = self.verify_reads(iter1, iter2)
        # Now count the reads outside the original region.
        subset_reads.reset()
        all_reads = 0
        for read in subset_reads:
            all_reads += 1
        extra = all_reads - num_reads
        print(num_reads, "read with", extra, "extra reads")

    def run_random_reads(self, num_tests):
        """
        Runs tests for valid ranges chosen randomly across the available
        regions.
        """
        for _ in range(num_tests):
            contig = random.choice(self.contigs)
            start = random.randint(0, contig.length)
            end = random.randint(start, contig.length)
            self.verify_query(contig.reference_name, start, end)

    def run_full_contig_fetch(self):
        """
        Gets all reads for contigs < 1Mb
        """
        for contig in self.contigs:
            if contig.length < 10**6:
                self.verify_query(contig.reference_name)

    def run_initial_reads(self):
        """
        Gets the first few reads from each contig.
        """
        for contig in self.contigs:
            self.verify_query(
                contig.reference_name,
                contig.initial_positions[0],
                contig.initial_positions[-1] + 1)
            self.verify_query(
                contig.reference_name,
                None,
                contig.initial_positions[-1] + 1)
            self.verify_query(
                contig.reference_name,
                max(0, contig.initial_positions[0] - 100),
                contig.initial_positions[0] + 1)
            self.verify_query(
                contig.reference_name,
                contig.initial_positions[-1],
                contig.initial_positions[-1] + 1)

    def run_final_reads(self):
        """
        Gets the first few reads from each contig.
        """
        for contig in self.contigs:
            if len(contig.final_positions) > 0:
                self.verify_query(
                    contig.reference_name,
                    contig.final_positions[0],
                    contig.final_positions[-1] + 1)
                self.verify_query(
                    contig.reference_name,
                    contig.final_positions[0],
                    None)
                self.verify_query(
                    contig.reference_name,
                    max(0, contig.final_positions[0] - 100),
                    contig.final_positions[0] + 1)
            else:
                print("Skipping final for ", contig.reference_name)

    def cleanup(self):
        self.alignment_file.close()
        os.unlink(self.temp_file_name)
        index_file = self.temp_file_name + ".bai"
        if os.path.exists(index_file):
            os.unlink(index_file)


if __name__ == "__main__":
    logging.basicConfig(format='%(asctime)s %(message)s', level=logging.DEBUG)

    parser = argparse.ArgumentParser(
        description="A simple tester application for the GA4GH streaming API.")
    parser.add_argument(
        "-V", "--version", action='version',
        version='%(prog)s {}'.format(gaclient.__version__))
    parser.add_argument('--verbose', '-v', action='count', default=0)
    parser.add_argument(
        "bam_file", type=str, help="The local BAM file to compare to")
    parser.add_argument(
        "url", type=str, help="The URL prefix of the server")
    parser.add_argument(
        "id", type=str, help="The ID of the ReadGroupSet")

    args = parser.parse_args()
    tester = ServerTester(args.bam_file, args.url, args.id)
    try:
        # TODO These don't work for the 8660_5#17.bam; why???
        tester.run_full_contig_fetch()
        tester.run_initial_reads()
        tester.run_final_reads()
        tester.run_random_reads(10**6)
    finally:
        tester.cleanup()
