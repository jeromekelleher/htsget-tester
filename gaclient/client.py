"""
Client API implementation of the draft GA4GH specification for bulk transfer
of read data.
"""
from __future__ import division
from __future__ import print_function

import base64
import os.path
import logging
import sys

import requests

IS_PY2 = sys.version_info[0] < 3

if IS_PY2:
    from urlparse import urlparse
else:
    from urllib.parse import urlparse


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
        self.__http_timeout = 5  # TODO add as a CLI parameter
        self.__max_attempts = 5
        self.__checkpoint_offset = None
        self.__is_stdout = False
        if output_file is None:
            self.__is_stdout = True
            # We cannot do retries when writing to stdout.
            self.__max_attempts = 1
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

    def __set_checkpoint(self):
        """
        Sets the checkpoint offset where we restart the read to.
        """
        if not self.__is_stdout:
            self.__checkpoint_offset = self.__output.tell()

    def __resume_checkpoint(self):
        """
        Resets the file pointer to the checkpoint offset.
        """
        if not self.__is_stdout:
            self.__output.seek(self.__checkpoint_offset)

    def __handle_http(self, http_url):
        """
        Handles a single HTTP request.
        """
        url = http_url["url"]
        method = http_url.get("method", "GET")
        # The dict may have had a null for method
        if method is None:
            method = "GET"
        headers = http_url.get("headers", None)
        body = http_url.get("body", None)
        logging.info("HTTP: {}: {}: Header = {}".format(method, url, headers))
        self.__set_checkpoint()
        successful = False
        num_attempts = 0
        while not successful:
            if num_attempts == self.__max_attempts:
                # TODO proper exception.
                raise ValueError("Too many retries")
            try:
                self.__resume_checkpoint()
                response = requests.request(
                    method, url, headers=headers, data=body, stream=True,
                    timeout=self.__http_timeout)
                response.raise_for_status()
                length = 0
                for chunk in response.iter_content(self.__chunk_size):
                    self.__store_chunk(chunk)
                    length += len(chunk)
                if "Content-Length" in response.headers:
                    content_length = int(response.headers["Content-Length"])
                    if content_length != length:
                        # TODO proper exception.
                        raise ValueError(
                            "Mismatch in downloaded length:{} != {}".format(
                                content_length, length))
                successful = True
            except requests.Timeout as te:
                logging.info("Timeout:{}".format(te))
            except requests.ConnectionError as ce:
                logging.info("Connection error: {}".format(ce))
            num_attempts += 1

    def __handle_data(self, data_uri):
        """
        Handles a single data URI.
        """
        parsed_url = urlparse(data_uri["url"])
        # TODO parse out the encoding properly.
        split = parsed_url.path.split(",", 1)
        data = base64.b64decode(split[1])
        logging.info("DATA: {}: length = {}".format(split[0], len(split[1])))
        self.__store_chunk(data)

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
        # TODO add format and other query parameters.
        url = os.path.join(self.__base_url, self.__id)
        response = requests.get(url, params=params)
        response.raise_for_status()
        ticket = response.json()
        for url_object in ticket["urls"]:
            url = url_object["url"]
            if url.startswith("http"):
                self.__handle_http(url_object)
            elif url.startswith("data"):
                self.__handle_data(url_object)
            else:
                raise ValueError("Unsupported URL:{}".format(url))

    def close(self):
        """
        Closes any open connections and files.
        """
        self.__output.close()
