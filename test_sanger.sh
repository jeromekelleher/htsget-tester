# do not exit on session disconnect
trap "" HUP  

#do not stop on errors
set +e 


set +x
echo ------------------------------------------------------------
echo Testing 2 small files ENCFF000VWO.bam and ENCFF284YOU.bam...
echo ------------------------------------------------------------
set -x

python tester.py -vv ENCFF000VWO.bam http://52.213.144.237:9090/npg_ranger/ga4gh/v.0.1/get/sample/ENCFF000VWO --client sanger-cli --max-references 20
python tester.py -vv ENCFF284YOU.bam http://52.213.144.237:9090/npg_ranger/ga4gh/v.0.1/get/sample/ENCFF284YOU_GRCH38 --client sanger-cli --max-references 20


if [[ "$@" == "fulltest" ]]
then
    set +x
	echo ------------------------------------------
    echo NOW TESTING BIG FILES...might take hours..
    echo ------------------------------------------
    set -x

	python tester.py -vv NA12878.bam http://52.213.144.237:9090/npg_ranger/ga4gh/v.0.1/get/sample/NA12878 --client sanger-cli --max-references 20
	python tester.py -vv NA12878.cram http://52.213.144.237:9090/npg_ranger/ga4gh/v.0.1/get/sample/NA12878_GRCH38 --client sanger-cli --max-references 20

	python tester.py -vv NA12891.bam http://52.213.144.237:9090/npg_ranger/ga4gh/v.0.1/get/sample/NA12891 --client sanger-cli --max-references 20
	python tester.py -vv NA12891.cram http://52.213.144.237:9090/npg_ranger/ga4gh/v.0.1/get/sample/NA12891_GRCH38 --client sanger-cli --max-references 20

	python tester.py -vv NA12892.bam http://52.213.144.237:9090/npg_ranger/ga4gh/v.0.1/get/sample/NA12892 --client sanger-cli --max-references 20
	python tester.py -vv NA12892.cram http://52.213.144.237:9090/npg_ranger/ga4gh/v.0.1/get/sample/NA12892_GRCH38 --client sanger-cli --max-references 20
fi

set -e