# do not exit on session disconnect
trap "" HUP

get_token()
{
        set +x
        echo "Obtaining access token..."
        read token < <(echo "GPqHSr-TfwC0nIkwR71SRhYTrrmMsdar")
        #curl -d "grant_type=password&client_id=f20cd2d3-682a-4568-a53e-4262ef54c8f4&client_secret=AMenuDLjVdVo4BSwi0QD54LL6NeVDEZRzEQUJ7hJOM3g4imDZBHHX0hNfKHPeQIGkskhtCmqAJtt_jm7EKq-rWw&username=ega-test-data@ebi.ac.uk&password=egarocks&scope=openid" -H "Content-Type: application/x-www-form-urlencoded" -k https://ega.ebi.ac.uk:8443/ega-openid-connect-server/token | python -c "import sys, json; print json.load(sys.stdin)[\"access_token\"]")
        set -x
}




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

	get_token
	python tester.py -vv NA12878.bam https://htsget.wtsi-npg-test.co.uk:10090/npg_ranger/authtoken/ga4gh/sample/NA12878 --bearer-token $token "$@"
	get_token
	python tester.py -vv NA12878.cram https://htsget.wtsi-npg-test.co.uk:10090/npg_ranger/authtoken/ga4gh/sample/NA12878_GRCH38 --bearer-token $token "$@"

	get_token
	python tester.py -vv NA12891.bam https://htsget.wtsi-npg-test.co.uk:10090/npg_ranger/authtoken/ga4gh/sample/NA12891 --bearer-token $token "$@"
	get_token
	python tester.py -vv NA12891.cram https://htsget.wtsi-npg-test.co.uk:10090/npg_ranger/authtoken/ga4gh/sample/NA12891_GRCH38 --bearer-token $token "$@"

	get_token
	python tester.py -vv NA12892.bam https://htsget.wtsi-npg-test.co.uk:10090/npg_ranger/authtoken/ga4gh/sample/NA12892 --bearer-token $token "$@"
	get_token
	python tester.py -vv NA12892.cram https://htsget.wtsi-npg-test.co.uk:10090/npg_ranger/authtoken/ga4gh/sample/NA12892_GRCH38 --bearer-token $token "$@"
fi

set +x
echo ------------------------------------------------------------
echo Testing 2 small files ENCFF000VWO.bam and ENCFF284YOU.bam...
echo ------------------------------------------------------------
set -x

get_token
python tester.py -vv ENCFF000VWO.bam https://htsget.wtsi-npg-test.co.uk:10090/npg_ranger/authtoken/ga4gh/sample/ENCFF000VWO --bearer-token $token "$@"
get_token
python tester.py -vv ENCFF284YOU.bam https://htsget.wtsi-npg-test.co.uk:10090/npg_ranger/authtoken/ga4gh/sample/ENCFF284YOU_GRCH38 --bearer-token $token "$@"

set -e
