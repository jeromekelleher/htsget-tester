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
import random
import subprocess
import sys
import tempfile
import time
import requests
import json

import htsget
import pysam
import humanize

from six.moves import zip
from six.moves.urllib.parse import urlencode


__version__ = "0.1.1"

FORMAT_BAM = "BAM"
FORMAT_CRAM = "CRAM"


def retry_command(cmd, filename=None, max_retries=5, retry_wait=5):
    """
    Runs the specified command and writes the stdout to the specified file (if
    supplied). If an exception occurs, retry a maximumum of max_retries times,
    waiting retry_wait seconds after every failure. If the maximum number of
    retries has been exceeded, raise the original exception.
    """
    done = False
    num_retries = 0
    while not done:
        try:
            outfile = None
            if filename is not None:
                outfile = open(filename, "w")
            subprocess.check_call(cmd, stdout=outfile)
            done = True
        except subprocess.CalledProcessError as cse:
            if num_retries == max_retries:
                raise cse
            logging.warn("Command failed, retrying: '{}'".format(" ".join(cmd)))
            time.sleep(retry_wait)
            num_retries += 1
        finally:
            if outfile is not None:
                outfile.close()


def htsget_api(
        url, filename, reference_name=None, reference_md5=None, start=None, end=None,
        data_format=None):
    logging.info("htsget-api: request (name='{}', md5='{}', start={}, end={})".format(
        reference_name, reference_md5, start, end))
    with open(filename, "wb") as tmp:
        htsget.get(
            url, tmp, reference_name=reference_name, reference_md5=reference_md5,
            start=start, end=end, data_format=data_format, timeout=120)


def htsget_cli(
        url, filename, reference_name=None, reference_md5=None, start=None, end=None,
        data_format=None):
    """
    Runs the htsget CLI program. Assumes that htsget has been installed and the
    CLI program is in PATH.
    """
    cmd = ["htsget", url, "-O", filename, "--timeout", "120"]

    if args.bearer_token:   cmd.extend(["--bearer-token", args.bearer_token])
    if args.ca_bundle:      os.environ["CURL_CA_BUNDLE"] = args.ca_bundle
    
    if reference_name is not None:
        cmd.extend(["-r", str(reference_name)])
    if start is not None:
        cmd.extend(["-s", str(start)])
    if end is not None:
        cmd.extend(["-e", str(end)])
    if data_format is not None:
        cmd.extend(["-f", str(data_format)])
    logging.info("htsget-cli run: {}".format(" ".join(cmd)))
    # We don't need to retry here because htsget automatically retries on errors.
    subprocess.check_call( cmd )


def dnanexus_cli(
        url, filename, reference_name=None, start=None, end=None, data_format=None):
    """
    Runs the htsnexus CLI program. Assumes that the script has been downloaded
    into the current working directory, and is executable.

    See https://github.com/dnanexus-rnd/htsnexus for details on how to download
    the script and use it.
    """
    url_segments = url.split('/')
    namespace = url_segments[-2]
    accession = url_segments[-1]
    server = "/".join(url_segments[:-2])
    cmd = ["./htsnexus.py", "-s", server]
   
    if reference_name is not None:
        ref = str(reference_name)
        # htsnexus doesn't support specifying start and end on their own,
        # so we have to work around this with boundary values.
        s = 0 if start is None else start
        # This should be a safe upper bound according to the SAM spec.
        e = 2**31 - 1 if end is None else end
        if start is not None or end is not None:
            ref += ":{}-{}".format(s, e)
        cmd.extend(["-r", ref])

    if args.bearer_token: cmd.extend(["--token", args.bearer_token])
    if args.ca_bundle:    cmd.extend(["--insecure"])

    cmd.extend([namespace, accession])    
    
    if data_format is not None:
        cmd.append(data_format)
        
    logging.info("htsnexus run: {}".format(" ".join(cmd)))
    retry_command(cmd, filename)


def sanger_cli(
        url, filename, reference_name=None, start=None, end=None, data_format=None):
    """
    Runs the Sanger Javascript client. Make sure you run these commands in the
    same directory as the current script, i.e.

    $ npm install npg_ranger
    $ python tester.py --client=sanger-cli <other args>
    """
    params = {}
    if reference_name is not None:
        params["referenceName"] = reference_name
    if start is not None:
        params["start"] = str(start)
    if end is not None:
        params["end"] = str(end)
    if data_format is not None:
        params["format"] = data_format
    if len(params) > 0:
        url += "?{}".format(urlencode(params))

    cmd = ["node_modules/.bin/npg_ranger_client", url, filename] 

    if args.ca_bundle: cmd.extend(["--with_ca", args.ca_bundle])  

    if args.bearer_token: 
        with open('./sanger_token.json', 'w') as outfile:
            json.dump({'token': ""+str(args.bearer_token)}, outfile)        
        cmd.extend(["--token_config",'./sanger_token.json'])

    logging.info("sanger client run: {}".format(" ".join(cmd)))
    retry_command(cmd)


def ena_cli(url, filename, reference_name=None, start=None, end=None, data_format=None):
    """
    Runs the ENA client ( https://github.com/enasequence/ena-ga4gh-read-client )
    """
    url_segments = url.split('/')
    dataset_id = url_segments[-1]
    endpoint_url = "/".join(url_segments[:-1]) + "/"

    cmd = ["ga4gh-dapi-client"]
    cmd.extend(["--endpoint-url", endpoint_url])
    cmd.extend(["--dataset-id", dataset_id])
    cmd.extend(["--output-file", filename])
    if data_format is not None:
        cmd.extend(["--format", data_format])
    if reference_name is not None:
        cmd.extend(["--reference-name", reference_name])
    if start is not None:
        cmd.extend(["--alignment-start", str(start)])
    if end is not None:
        cmd.extend(["--alignment-stop", str(end)])

    logging.info("ENA run: {}".format(" ".join(cmd)))
    retry_command(cmd)

def ega_cli(url, filename, reference_name=None, start=None, end=None, data_format=None):
    """
    Runs the EGA client by Alexander Senf ( https://github.com/EbiEga/ega-htsget-client )
    """
    url_segments = url.split('/')
    dataset_id = url_segments[-1]
    endpoint_url = "/".join(url_segments[:-1]) + "/"


    cmd = ["java"]
    cmd.extend(["-jar", "EgaHtsgetClient.jar"])

    if args.bearer_token: cmd.extend(["--oauth-token", args.bearer_token])
    
    cmd.extend(["--endpoint-url", endpoint_url])
    cmd.extend(["--dataset-id", dataset_id])
    cmd.extend(["--output-file", filename])
    cmd.extend(["--debug"])
    
    if data_format is not None:
        cmd.extend(["--format", data_format])
    if reference_name is not None:
        cmd.extend(["--reference-name", reference_name])
    if start is not None:
        cmd.extend(["--alignment-start", str(start)])
    if end is not None:
        cmd.extend(["--alignment-stop", str(end)])

    logging.info("EGA run: {}".format(" ".join(cmd)))
    retry_command(cmd)
    


def samtools_cli( url, filename, reference_name=None, start=None, end=None, data_format=None):
     """
     Runs the samtools CLI program.
     https://github.com/samtools/samtools
     """

     url += '?format=' + ('BAM' if data_format is None else str(data_format))
     format_flag = '-b' if data_format is None or data_format == FORMAT_BAM else '-C'

     if reference_name is not None:
         url += '&referenceName='+str(reference_name)
     if start is not None:
         url += '&start='+ str(start)
     if end is not None:
         url += '&end='+ str(end)        
     
     cmd = ["samtools", "view", format_flag, url]

     #--- uncomment these 2 lines to debug samtools output
     #format_flag = "-c"
     #cmd = ["htsfile", "-vvvvvvvvv", format_flag, url]

     if args.bearer_token:         
         with open('samtools_token', 'w') as outfile:
             outfile.write( args.bearer_token )          
         os.environ["HTS_AUTH_LOCATION"] = "samtools_token"           

     if args.ca_bundle: os.environ["CURL_CA_BUNDLE"] = args.ca_bundle;

     logging.info("samtools run: {}".format(" ".join(cmd)))
     logging.info("tmp file : "+str(filename) )
     retry_command(cmd, filename)


# All clients must be registered in this map to be available for use.
client_map = {
    "htsget-api": htsget_api,
    "htsget-cli": htsget_cli,
    "dnanexus-cli": dnanexus_cli,
    "sanger-cli": sanger_cli,
    "ena-cli": ena_cli,
    "ega-cli": ega_cli,
    "samtools-cli": samtools_cli
}

class TestFailedException(Exception):
    """
    Exception raised when we know we've failed.
    """


class DownloadFailedException(Exception):
    """
    Exception raised when we downloading the data failed for some reason.
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


class TestReport(object):
    """
    A report on a single test run.
    """
    def __init__(
            self, num_references=None, total_queries=None, total_downloaded_data=None,
            total_download_time=None, mismatch_counts=None):
        self.num_references = num_references
        self.total_queries = total_queries
        self.total_downloaded_data = total_downloaded_data
        self.total_download_time = total_download_time
        self.mismatch_counts = mismatch_counts


class Contig(object):
    """
    Represents a single contig within the BAM file.
    """
    def __init__(
            self, reference_name, reference_md5, length, start_positions, end_positions):
        self.reference_name = reference_name
        self.reference_md5 = reference_md5
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


def simple_equality(v1, v2, read):
    """
    Equality for simple fields where equality can be tested directly.
    """
    return v1 == v2


def sorted_equality(v1, v2, read):
    """
    Equality for fields where the values must be sorted before equality tested.
    """
    return sorted(v1) == sorted(v2)


class ServerTester(object):
    """
    Class to run systematic tests on a server based on a local indexed
    BAM file.
    """

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
        extension = os.path.splitext(self.source_file_name)[1].lower()
        if extension == ".cram":
            self.data_format = FORMAT_CRAM
        elif extension == ".bam":
            self.data_format = FORMAT_BAM
        else:
            raise ValueError("Unknown file format: {}. Please use .bam or .cram".format(
                extension))
        self.temp_file_name = None
        self.subset_alignment_file = None
        self.num_boundary_reads = num_boundary_reads
        self.max_references = max_references
        self.client = client
        self.contigs = []
        self.mismatch_counts = collections.Counter()
        self.mismatch_report_threshold = 20  # Avoid flooding log.
        self.total_queries = 0
        self.total_downloaded_data = 0
        self.total_download_time = 1e-8  # Avoid zero division problems
        self.fields = [
            Field("pos", simple_equality),
            Field("query_name", simple_equality),
            Field("reference_name", simple_equality),
            Field("cigarstring", simple_equality),
            Field("query_alignment_sequence", simple_equality),
            Field("query_alignment_qualities", simple_equality),
            Field("template_length", simple_equality),
            Field("mapping_quality", simple_equality),
            Field("tags", sorted_equality),
            Field("next_reference_id", self.next_reference_equality),
            Field("next_reference_start", self.next_reference_start_equality),
            Field("flag", self.flags_equality),
            # TODO fill in remaining BAM fields.
        ]

    def next_reference_equality(self, local_rid, remote_rid, local_read):
        """
        Compares the two specified references.
        """
        # We need to convert the reference IDs back into names so that
        # we can compare them.
        r1 = self.alignment_file.references[local_rid]
        r2 = self.subset_alignment_file.references[remote_rid]
        if self.filter_unmapped and local_read.mate_is_unmapped:
            ret = True
            # We allow the remote reference ID to be unset if we are filtering
            # out unmapped reads and the mate is unmapped
            if remote_rid != -1:
                ret = r1 == r2
        else:
            # If both are unset, they are equal. We need to check this as the
            # order of the references is abitrary, and references[-1] is just
            # the last element in the list.
            if local_rid == -1 and remote_rid == -1:
                ret = True
            else:
                ret = r1 == r2
        return ret

    def next_reference_start_equality(self, v1, v2, local_read):
        """
        Compares the two specified values.
        """
        if self.filter_unmapped and local_read.mate_is_unmapped:
            ret = True
            # We allow the remote value to be unset if we are filtering out unmapped
            # reads.
            if v2 != -1:
                ret = v1 == v2
        else:
            ret = v1 == v2
        return ret

    def flags_equality(self, local_flags, remote_flags, local_read):
        """
        Compares the specified flags values.
        """
        condition = (
            self.filter_unmapped and
            local_read.mate_is_unmapped and
            local_read.mate_is_reverse)
        if condition:
            # If we are filtering out unmapped reads and the mate is unmapped, then
            # the mate_is_reverse flag can't be set on the remote flags. This is
            # what the offset by 32 effectively means.
            ret = local_flags == remote_flags + 32
        else:
            ret = local_flags == remote_flags
        return ret

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
        Scans the input file and initialises the data structures we need for the
        tests.
        """
        fd, self.temp_file_name = tempfile.mkstemp(
            prefix="gastream_test_", dir=self.tmpdir, suffix="." + self.data_format)
        os.close(fd)
        # Determine the bounds of the individual contigs.
        total_references = len(self.alignment_file.lengths)
        num_references = min(self.max_references, total_references)
        logging.info("Reading file {}".format(self.source_file_name))
        for j in range(num_references):
            reference_name = self.alignment_file.references[j]
            sq = self.alignment_file.header['SQ'][j]
            reference_md5 = sq.get('M5', None)
            length = self.alignment_file.lengths[j]
            assert sq['LN'] == length
            assert sq['SN'] == reference_name
            start_positions = self.get_start_positions(reference_name)
            if len(start_positions) > 0:
                end_positions = self.get_end_positions(reference_name, length)
                contig = Contig(
                    reference_name, reference_md5, length, start_positions,
                    end_positions)
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
            if not field.equals(v1, v2, r1):
                self.mismatch_counts[field.name] += 1
                if self.mismatch_counts[field.name] < self.mismatch_report_threshold:
                    logging.warning("Mismatch at {}:{}.{}:: {} != {}".format(
                        r1.reference_name, r1.pos, field.name, v1, v2))
                    if type(v1) is list and type(v2) is list:
                        v1 = set(v1)
                        v2 = set(v2)
                        logging.warning(
                            "Extra elements in 1st set (v1-v2) : {}".format(
                                sorted(v1 - v2)))
                        logging.warning(
                            "Extra elements in 2nd set (v2-v1) : {}".format(
                                sorted(v2 - v1)))
                elif self.mismatch_counts[field.name] == self.mismatch_report_threshold:
                    logging.warning(
                        "More than {} mismatches logged. Not reporting any more".format(
                            self.mismatch_report_threshold))

    def verify_reads(self, iter1, iter2):
        """
        Verifies that the specified iterators contain the same set of
        reads. Returns the number of reads in the iterators.
        """

        def check_reads_for_position(d1, d2):
            num_checks = 0
            if len(d1) != len(d2):
                raise TestFailedException("different numbers of reads returned")
            for k in d1.keys():
                if k not in d2:
                    raise TestFailedException("{} not found".format(k))
                self.verify_reads_equal(d1[k], d2[k])
                num_checks += 1
            return num_checks

        # Because the order of reads for a given position is unspecified,
        # we gather the reads for each iterator into dictionaries indexed
        # by a (hopefully) unique combination of fields.
        d1 = {}
        d2 = {}
        last_pos = -1
        num_reads = 0
        total_checks = 0
        for r1, r2 in zip(iter1, iter2):
            num_reads += 1
            if r1.pos != r2.pos:
                raise TestFailedException("Non-matching read positions(local != remote) : ({}){}!={}({})".format( getattr(r1,"query_name"), r1.pos, r2.pos, getattr(r2,"query_name")) )
            if r1.pos != last_pos:
                total_checks += check_reads_for_position(d1, d2)
                d1.clear()
                d2.clear()
                last_pos = r1.pos
            k = hashkey(r1)
            if k in d1:
                raise TestFailedException("Duplicate read: {}".format(k))
            d1[k] = r1
            k = hashkey(r2)
            if k in d2:
                raise TestFailedException("Duplicate read: {}".format(k))
            d2[k] = r2

        # Check the last set of reads in the dictionaries.
        total_checks += check_reads_for_position(d1, d2)
        # make sure both iterators are empty
        r1 = next(iter1, None)
        r2 = next(iter2, None)
        if r1 is not None or r2 is not None:
            extra_reads = 0;            
            while r1 is not None: 
                extra_reads += 1
                r1 = next(iter1, None)
            while r2 is not None: 
                extra_reads += 1                
                r2 = next(iter2, None)                
            raise TestFailedException("Total number of reads not matching {} vs {}".format(num_reads, num_reads+extra_reads+1))
        assert total_checks == num_reads
        return num_reads

    def verify_query(self, reference_name=None, start=None, end=None):
        """
        Runs the specified query and verifies the result.
        """
        self.total_queries += 1
        # We use the wall-clock time here
        before = time.time()
        self.client(
            self.server_url, self.temp_file_name, reference_name=reference_name,
            start=start, end=end, data_format=self.data_format)
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
        logging.debug("Indexed {} file in {:.3f} CPU seconds".format(
            self.data_format, duration))

        # Analyse the reads
        before = time.clock()
        iter1 = self.alignment_file.fetch(reference_name, start, end)
        if self.data_format == FORMAT_CRAM:
            # Due to a bug in htslib, we cannot use the standard code path below
            # for CRAM files that have no records. The workaround is to open the
            # CRAM file without an index and see if it is emtpy. If it is, we
            # skip the rest of the tests. Once the upstream bug in htslib has
            # been fixed and incorporated into pysam, we should remove this
            # special case.
            # See https://github.com/pysam-developers/pysam/issues/483
            tmp = pysam.AlignmentFile(
                self.temp_file_name, filepath_index="/no/such/index/exists.crai")
            empty = True
            for _ in tmp.head(1):
                empty = False
            tmp.close()
            if empty:
                # Make sure the original iterator is also empty...
                count = 0
                for _ in iter1:
                    count += 1
                if count != 0:
                    raise TestFailedException(
                        "Downloaded CRAM empty, but original contains {} reads".format(
                            count))
                return
        try:
            self.subset_alignment_file = pysam.AlignmentFile(self.temp_file_name)
            iter2 = self.subset_alignment_file.fetch(reference_name, start, end)
        except ValueError as ve:
            raise DownloadFailedException("Reading downloaded data: {}".format(ve))
        if self.filter_unmapped:
            iter1 = filter_unmapped(iter1)
        num_reads = self.verify_reads(iter1, iter2)
        duration = time.clock() - before
        logging.debug("Checked reads in {:.3f} CPU seconds".format(duration))
        if self.data_format == FORMAT_BAM:
            # Count the reads outside the original region. This only works for
            # BAM files, so we don't bother for CRAM.
            self.subset_alignment_file.reset()
            all_reads = 0
            for read in self.subset_alignment_file:
                all_reads += 1
            extra = all_reads - num_reads
            logging.info("Downloaded {} reads with {} extra".format(num_reads, extra))
        else:
            logging.info("Downloaded {} reads".format(num_reads))

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
            self.verify_query(contig.reference_name, start=start, end=end)

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
            values = [
                (contig.start_positions[0], contig.start_positions[-1] + 1),
                (None, contig.start_positions[-1] + 1),
                (max(0, contig.start_positions[0] - 100), contig.start_positions[0] + 1),
                (contig.start_positions[-1], contig.start_positions[-1] + 1)
            ]
            for start, end in values:
                self.verify_query(
                    reference_name=contig.reference_name, start=start, end=end)
                # self.verify_query(
                #     reference_md5=contig.reference_md5, start=start, end=end)

    def run_end_reads(self):
        """
        Gets the first few reads from each contig.
        """
        logging.info("Starting contig end queries")
        for contig in self.contigs:
            if len(contig.end_positions) > 0:
                values = [
                    (contig.end_positions[0], contig.end_positions[-1] + 1),
                    (contig.end_positions[0], None),
                    (max(0, contig.end_positions[0] - 100), contig.end_positions[0] + 1),
                ]
                for start, end in values:
                    self.verify_query(
                        reference_name=contig.reference_name, start=start, end=end)
            else:
                logging.info("Skipping end reads for {}".format(contig.reference_name))

    def cleanup(self):
        self.alignment_file.close()
        if os.path.exists(self.temp_file_name):
            os.unlink(self.temp_file_name)
            index_file = self.temp_file_name + ".bai"
            if os.path.exists(index_file):
                os.unlink(index_file)
            index_file = self.temp_file_name + ".crai"
            if os.path.exists(index_file):
                os.unlink(index_file)

    def report(self):
        return TestReport(
            num_references=len(self.contigs),
            total_queries=self.total_queries,
            total_downloaded_data=self.total_downloaded_data,
            total_download_time=self.total_download_time,
            mismatch_counts=self.mismatch_counts)

    def print_report(self):
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

    version_report = "%(prog)s {} (htsget {}; Python {})".format(
            __version__, htsget.__version__, ".".join(map(str, sys.version_info[:3])))

    parser = argparse.ArgumentParser(
        description="A simple tester application for the GA4GH streaming API.")
    parser.add_argument(
        "-V", "--version", action='version',
        version=version_report)
    parser.add_argument('--verbose', '-v', action='count', default=0)
    parser.add_argument(
        "source_file", type=str, help="The local BAM/CRAM file to compare to")
    parser.add_argument(
        "url", type=str, help="The server url corresponding to the input source file.")
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

    parser.add_argument(
        "--bearer-token", type=str, help="Bearer token")
        
    parser.add_argument(
                "--ca-bundle", type=str, help="CA bundle")
        

    args = parser.parse_args()
    log_level = logging.WARNING
    if args.verbose == 1:
        log_level = logging.INFO
    if args.verbose >= 2:
        log_level = logging.DEBUG
    logging.basicConfig(format='%(asctime)s %(message)s', level=log_level)
    random.seed(args.random_seed)
    tester = ServerTester(
        args.source_file, args.url, filter_unmapped=args.filter_unmapped,
        tmpdir=args.tmpdir, max_references=args.max_references,
        client=client_map[args.client])
    exit_status = 1
    try:
        tester.initialise()
        tester.run_full_contig_fetch()
        tester.run_start_reads()
        tester.run_end_reads()
        tester.run_random_reads(args.num_random_reads)
        exit_status = 0
    except DownloadFailedException as dfe:
        print("Download failed: ", dfe, file=sys.stderr)
        exit_status = 1
    except TestFailedException as tfe:
        print("Test failed:", tfe, ". Is the file correct?", file=sys.stderr)
        exit_status = 1
    except KeyboardInterrupt:
        logging.warn("Interrupted!")
        exit_status = 0
    finally:
        tester.cleanup()
    # If everything went OK, write out a report.
    if exit_status == 0:
        tester.print_report()
    sys.exit(exit_status)
