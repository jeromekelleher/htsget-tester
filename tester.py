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

import gaclient


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
    logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)

    # jk
    bam_file = "/home/jk/public_html/1kg-data/NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.bam"
    server_url = "http://holly:8000/reads/fetch"
    server_id = "1046405015"

    # htsnexus
    server_url = "http://htsnexus.rnd.dnanex.us/v1/reads/"
    server_id = "1000genomes_low_coverage/NA12878"

    # EBI
    bam_file = "tmp__NOBACKUP__/input_data/8660_5#17.bam"
    server_url = "http://ga4gh.ebi.ac.uk/ticket"
    server_id = "ERR217910"

    # tester = ServerTester(na12878, server_url, server_id)
    # run_edge_tests(na12878, server_url, server_id)
    run_random_tests(bam_file, server_url, server_id)


