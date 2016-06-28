"""
Client API implementation of the draft GA4GH specification for bulk transfer
of read data.
"""
from __future__ import division
from __future__ import print_function

import base64
import ftplib
import os.path
import sys

import requests

IS_PY2 = sys.version_info[0] < 3

if IS_PY2:
    from urlparse import urlparse
    from urlparse import urlparse
else:
    from urllib.parse import urlparse
    from urllib.parse import unquote

class Client(object):
    """
    Client for the GA4GH streaming API.
    """
    def __init__(
            self, base_url=None, id_=None, reference_name=None, start=None,
            end=None, output_file=None):
        self.__base_url = base_url
        self.__id = id_
        self.__reference_name = reference_name
        self.__start = start
        self.__end = end
        self.__chunk_size = 64 * 1024  # TODO add as CLI parameter
        if output_file is None:
            if IS_PY2:
                self.__output = sys.stdout
            else:
                self.__output = sys.stdout.buffer
        else:
            self.__output = open(output_file, "wb")

    def __store_chunk(self, chunk):
        """
        Writes the specified chunk of data to the output.
        """
        self.__output.write(chunk)

    def download(self):
        """
        Runs the query.
        """
        params = {}
        if self.__reference_name is not None:
            params["referenceName"] = self.__reference_name
        if self.__start is not None:
            params["start"] = self.__start
        if self.__end is not None:
            params["end"] = self.__end
        params["format"] = "BAM"
        url = os.path.join(self.__base_url, self.__id)
        response = requests.get(url, params=params)
        response.raise_for_status()
        ticket = response.json()
        if "prefix" in ticket:
            prefix = base64.b64decode(ticket["prefix"])
            self.__output.write(prefix)
        byte_ranges = None
        if "byteRanges" in ticket:
            if ticket["byteRanges"] is not None:
                byte_ranges = ticket["byteRanges"]
                # TODO raise a proper exception here.
                assert len(byte_ranges) == len(ticket["urls"])
        extra_headers = {}
        if "httpRequestHeaders" in ticket:
            if ticket["httpRequestHeaders"] is not None:
                extra_headers.update(ticket["httpRequestHeaders"])
        for j, url in enumerate(ticket["urls"]):
            parsed_url = urlparse(url)
            if parsed_url.scheme.startswith("http"):
                headers = dict(extra_headers)
                if byte_ranges is not None:
                    start = byte_ranges[j]["start"]
                    end = byte_ranges[j]["end"]
                    headers["Range"] = "bytes={}-{}".format(start, end - 1)
                response = requests.get(url, stream=True, headers=headers)
                response.raise_for_status()
                length = 0
                for chunk in response.iter_content(self.__chunk_size):
                    self.__store_chunk(chunk)
                    length += len(chunk)
                if "Content-Length" in response.headers:
                    content_length = int(response.headers["Content-Length"])
                    if content_length != length:
                        raise ValueError(
                            "Mismatch in downloaded length:{} != {}".format(
                                content_length, length))
            elif parsed_url.scheme == "ftp":
                print(url)
                print(byte_ranges)
                with ftplib.FTP(parsed_url.netloc) as ftp:
                    ftp.set_debuglevel(1)
                    ftp.login()
                    # TODO: what does the protocol say about quoting of URLs??
                    # Should we really have to unquote these URLs for FTP?
                    directory, filename = map(
                        unquote, os.path.split(parsed_url.path))
                    ftp.cwd(directory)
                    ftp.retrbinary(
                        "RETR {}".format(filename), self.__store_chunk)
            else:
                raise ValueError("Unsupported URL:{}".format(url))
        if "suffix" in ticket:
            prefix = base64.b64decode(ticket["suffix"])
            self.__output.write(prefix)

    def close(self):
        """
        Closes any open connections and files.
        """
        self.__output.close()
