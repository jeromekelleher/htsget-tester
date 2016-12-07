# GA4GH Streaming API tester.

A tester for the draft API spec at 

https://docs.google.com/document/d/1OSPfxdJ3uPoCfUVzMaekCOPF5sNEwqkJEUj-SjlECy0

## Usage

```bash
$ git clone https://github.com/jeromekelleher/ga4gh-streaming-tester.git
$ pip install -r requirements.txt
$ python tester.py -vv ENCFF000VWO.bam http://htsnexus.rnd.dnanex.us/v1/reads/ENCODE/ENCFF000VWO
```

This should work with Python 2.7 or 3.x. For full command line options, 
please use ``python tester.py --help``.

# Supporting another client

To add support for another client:

1. Add a function following the example of ``htsget_cli``.
2. Update the ``client_map`` dictionary in ``main()``.
 

