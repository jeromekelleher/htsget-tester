"""
Client API implementation of the draft GA4GH specification for bulk transfer
of read data.
"""
from __future__ import division
from __future__ import print_function

import base64
import os.path
import sys

import requests

IS_PY2 = sys.version_info[0] < 3

if IS_PY2:
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

        # print("url = ", url)
        # print("headers = ", headers)
        response = requests.request(
            method, url, headers=headers, data=body, stream=True)
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


    def __handle_data(self, data_uri):
        """
        Handles a single data URI.
        """
        parsed_url = urlparse(data_uri["url"])
        # TODO parse out the encoding properly.
        data = base64.b64decode(parsed_url.path.split(",", 1)[1])
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
