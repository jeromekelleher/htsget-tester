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
    echo NOW TESTING BIG FILES...might take hours..
    echo ------------------------------------------
    set -x

	python tester.py -vv --filter-unmapped NA12878.bam http://104.196.18.135/reads/CO_Ph7XUCRDW0-riiPr48fgB "$@"
	#python tester.py -vv --filter-unmapped NA12878.cram http://104.196.18.135/reads/CO_Ph7XUCRCX3-no2-WqtqcB "$@"

	python tester.py -vv --filter-unmapped NA12891.bam http://104.196.18.135/reads/CO_Ph7XUCRDfpdqpjO3Tim0  "$@"
	#python tester.py -vv --filter-unmapped NA12891.cram http://104.196.18.135/reads/CO_Ph7XUCRDngbisk52h2KEB "$@"

	python tester.py -vv --filter-unmapped NA12892.bam http://104.196.18.135/reads/CO_Ph7XUCRCMtrWtwPmA4dgB "$@"
	#python tester.py -vv --filter-unmapped NA12892.cram http://104.196.18.135/reads/CO_Ph7XUCRDnpcT9pKrW-tUB "$@"
fi


set +x
echo ------------------------------------------------------------
echo Testing 2 small files ENCFF000VWO.bam and ENCFF284YOU.bam...
echo ------------------------------------------------------------
set -x

python tester.py -vv --filter-unmapped ENCFF000VWO.bam http://104.196.18.135/reads/CO_Ph7XUCRCSn-mRrofkhXY "$@"
python tester.py -vv --filter-unmapped ENCFF284YOU.bam http://104.196.18.135/reads/CO_Ph7XUCRCi9Lasl7yY5Z0B "$@"

set -e