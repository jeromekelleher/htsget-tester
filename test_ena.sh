# do not exit on session disconnect
trap "" HUP  

#do not stop on errors
set +e 


if [[ "$@" == *fulltest ]]
then
    # remove last arg ( 'fulltest' )
    set -- "${@:1:$(($#-1))}"

    set +x
	echo ------------------------------------------
    echo TESTING BIG FILES.....might take hours....
    echo ------------------------------------------
    set -x

	python tester.py -vv NA12878.bam http://ga4gh.ebi.ac.uk/ticket/NA12878.bam  "$@"
	python tester.py -vv NA12878.cram http://ga4gh.ebi.ac.uk/ticket/NA12878.cram "$@"

	python tester.py -vv NA12891.bam http://ga4gh.ebi.ac.uk/ticket/NA12891.bam "$@"
	python tester.py -vv NA12891.cram http://ga4gh.ebi.ac.uk/ticket/NA12891.cram "$@"

	python tester.py -vv NA12892.bam http://ga4gh.ebi.ac.uk/ticket/NA12892.bam "$@"
	python tester.py -vv NA12892.cram http://ga4gh.ebi.ac.uk/ticket/NA12892.cram "$@"
fi

set +x
echo ------------------------------------------------------------
echo Testing 2 small files ENCFF000VWO.bam and ENCFF284YOU.bam...
echo ------------------------------------------------------------
set -x

python tester.py -vv ENCFF000VWO.bam http://ga4gh.ebi.ac.uk/ticket/ENCFF000VWO "$@"
python tester.py -vv ENCFF284YOU.bam http://ga4gh.ebi.ac.uk/ticket/ENCFF284YOU "$@"

set -e