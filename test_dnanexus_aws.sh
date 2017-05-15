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

	python tester.py -vv NA12878.bam http://htsnexus.rnd.dnanex.us/v1/reads/BroadHiSeqX_b37/NA12878 "$@"
	python tester.py -vv NA12878.cram http://htsnexus.rnd.dnanex.us/v1/reads/BroadHiSeqX_hg38/NA12878  "$@"

	python tester.py -vv NA12891.bam http://htsnexus.rnd.dnanex.us/v1/reads/BroadHiSeqX_b37/NA12891  "$@"
	python tester.py -vv NA12891.cram http://htsnexus.rnd.dnanex.us/v1/reads/BroadHiSeqX_hg38/NA12891  "$@"

	python tester.py -vv NA12892.bam http://htsnexus.rnd.dnanex.us/v1/reads/BroadHiSeqX_b37/NA12892  "$@"
	python tester.py -vv NA12892.cram http://htsnexus.rnd.dnanex.us/v1/reads/BroadHiSeqX_hg38/NA12892  "$@"
fi

set +x
echo ------------------------------------------------------------
echo Testing 2 small files ENCFF000VWO.bam and ENCFF284YOU.bam...
echo ------------------------------------------------------------
set -x

python tester.py -vv ENCFF000VWO.bam http://htsnexus.rnd.dnanex.us/v1/reads/ENCODE/ENCFF284YOU  "$@"
python tester.py -vv ENCFF284YOU.bam http://htsnexus.rnd.dnanex.us/v1/reads/ENCODE/ENCFF000VWO  "$@"

set -e