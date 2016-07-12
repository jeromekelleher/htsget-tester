from __future__ import print_function
from __future__ import division

import pysam
import os
import itertools
import time
import tempfile
import random
import os.path

import gaclient


class ServerTester(object):

    def __init__(source_bam, url, id_):
        self.source_bam = source_bam
        self.outfile = "tmp__NOBACKUP__/test_{}.bam".format(os.getpid())
        boundary_block_size = 10
        in_af = pysam.AlignmentFile(source_bam)
        self.reference_names = in_af.references
        self.reference_lengths = in_af.lengths

        # TODO find the coordinates of the first and last k records
        # so that we can test behaviour at the edges.
        self.initial_block_coordinates = []
        for ref_name in self.references_names:
            num_reads = 0
            for read in in_af.fetch(ref_name):
                if num_reads == boundary_block_size:
                    break
        in_af.close()


def run_random_tests(source_bam, url, id_):

    in_af = pysam.AlignmentFile(source_bam)
    ref_lengths = []
    for r, l in zip(in_af.references, in_af.lengths):
        ref_lengths.append((r, l))
    in_af.close()
    outfile = "tmp__NOBACKUP__/test_{}.bam".format(os.getpid())
    while True:
        reference_name, length = random.choice(ref_lengths)
        start = random.randint(0, length)
        end = random.randint(start, length)
        # reference_name, start, end = "GL000246.1", 12866, 38114
        print("reference_name, start, end = {}, {}, {}".format(
            reference_name, start, end))
        client = gaclient.Client(
            url, id_, reference_name, start, end, outfile)
        before = time.clock()
        client.download()
        client.close()
        duration = time.clock() - before
        print("Downloaded in {} seconds".format(duration))

        all_reads = pysam.AlignmentFile(source_bam)
        pysam.index(outfile)
        subset_reads = pysam.AlignmentFile(outfile)
        iter1 = all_reads.fetch(reference_name, start, end)
        iter2 = subset_reads.fetch(reference_name, start, end)
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
        # Now count the reads outside the original region.
        subset_reads.reset()
        all_reads = 0
        for read in subset_reads:
            all_reads += 1
        extra = all_reads - num_reads
        print(num_reads, "read with", extra, "trailing reads")

if __name__ == "__main__":
    na12878 = "/home/jk/public_html/1kg-data/NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.bam"
    server_url = "http://holly:8000/reads/fetch"
    server_id = "1046405015"
    # tester = ServerTester(na12878, server_url, server_id)
    # run_edge_tests(na12878, server_url, server_id)
    run_random_tests(na12878, server_url, server_id)

