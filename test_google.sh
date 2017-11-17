# do not exit on session disconnect
trap "" HUP

get_token()
{
        set +x
        echo "Obtaining access token..."
        read token < <(curl --request POST --data '{"token_uri":"https://www.googleapis.com/oauth2/v4/token","refresh_token":"1/A3OHTTywlfJAqCU95uxNYGx9a8PUhbZjhGC4O8A5_sM","client_secret":"T-hgpqETtzseRHYn6egZSPXK","client_id":"819467526948-ochlfup6bvc444455c7cirqovkikl7ik.apps.googleusercontent.com"}' https://developers.google.com/oauthplayground/refreshAccessToken | python -c "import sys, json; print json.load(sys.stdin)[\"access_token\"]")
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
	python tester.py -vv --filter-unmapped NA12878.bam  https://35.196.212.220/reads/ga4gh-demo/NA12878.bam --ca-bundle htsget-demo.pem  --bearer-token $token "$@"
	get_token
	python tester.py -vv --filter-unmapped NA12878.cram https://35.196.212.220/reads/ga4gh-demo/NA12878.cram --ca-bundle htsget-demo.pem  --bearer-token $token "$@"

    get_token
	python tester.py -vv --filter-unmapped NA12891.bam https://35.196.212.220/reads/ga4gh-demo/NA12891.bam --ca-bundle htsget-demo.pem   --bearer-token $token "$@"
	get_token
	python tester.py -vv --filter-unmapped NA12891.cram https://35.196.212.220/reads/ga4gh-demo/NA12891.cram --ca-bundle htsget-demo.pem  --bearer-token $token "$@"

	get_token
	python tester.py -vv --filter-unmapped NA12892.bam https://35.196.212.220/reads/ga4gh-demo/NA12892.bam --ca-bundle htsget-demo.pem  --bearer-token $token "$@"
	get_token
	python tester.py -vv --filter-unmapped NA12892.cram https://35.196.212.220/reads/ga4gh-demo/NA12892.cram --ca-bundle htsget-demo.pem  --bearer-token $token "$@"
fi


set +x
echo ------------------------------------------------------------
echo Testing 2 small files ENCFF000VWO.bam and ENCFF284YOU.bam...
echo ------------------------------------------------------------
set -x

get_token
python tester.py -vv --filter-unmapped ENCFF000VWO.bam https://35.196.212.220/reads/ga4gh-demo/ENCFF000VWO.bam --ca-bundle htsget-demo.pem  --bearer-token $token "$@"
get_token
python tester.py -vv --filter-unmapped ENCFF284YOU.bam https://35.196.212.220/reads/ga4gh-demo/ENCFF284YOU.bam --ca-bundle htsget-demo.pem  --bearer-token $token "$@"

set -e