# Htsget server tester.

A test application for the GA4GH [htsget](http://samtools.github.io/hts-specs/htsget.html) protocol.

## Usage

```bash
$ git clone https://github.com/jeromekelleher/htsget-tester.git
$ pip install -r requirements.txt
$ python tester.py -vv ENCFF000VWO.bam http://htsnexus.rnd.dnanex.us/v1/reads/ENCODE/ENCFF000VWO
```

This should work with Python 2.7 or 3.x. For full command line options, 
please use ``python tester.py --help``.

# Supporting another client

To add support for another client:

1. Add a function following the example of ``htsget_cli``.
2. Update the ``client_map`` dictionary in ``main()``.
 

