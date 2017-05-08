# do not exit on session disconnect
trap "" HUP  

#do not stop on errors
set +e 


set +x
echo ------------------------------------------------------------
echo Testing 2 small files ENCFF000VWO.bam and ENCFF284YOU.bam...
echo ------------------------------------------------------------
set -x

python tester.py -vv ENCFF000VWO.bam http://52.213.144.237:9090/npg_ranger/ga4gh/v.0.1/get/sample/ENCFF000VWO "$@"
python tester.py -vv ENCFF284YOU.bam http://52.213.144.237:9090/npg_ranger/ga4gh/v.0.1/get/sample/ENCFF284YOU_GRCH38 "$@"


if [[ "$@" == *fulltest ]]
then
    set +x
	echo ------------------------------------------
    echo NOW TESTING BIG FILES...might take hours..
    echo ------------------------------------------
    set -x

	python tester.py -vv NA12878.bam http://ga4gh.ebi.ac.uk/ticket/NA12878.bam  "$@"
	python tester.py -vv NA12878.cram http://ga4gh.ebi.ac.uk/ticket/NA12878.cram "$@"

	python tester.py -vv NA12891.bam http://ga4gh.ebi.ac.uk/ticket/NA12891.bam "$@"
	python tester.py -vv NA12891.cram http://ga4gh.ebi.ac.uk/ticket/NA12891.cram "$@"

	python tester.py -vv NA12892.bam http://ga4gh.ebi.ac.uk/ticket/NA12892.bam "$@"
	python tester.py -vv NA12892.cram http://ga4gh.ebi.ac.uk/ticket/NA12892.cram "$@"
fi

set -e