#!/bin/bash
set -e

DIR=$(dirname $(realpath "$0")) 	# locate folder where this sh-script is located in
SCRIPT="./run_tests.inp"
PACKAGE_NAME="logging"

cd $DIR
echo "Switched to ${DIR}"

gretlcli -b -e -q ${SCRIPT}

if [ $? -eq 0 ]
then
  echo "Success: All tests passed for 'logging'."
  exit 0
else
  echo "Failure: Tests not passed for 'logging'." >&2
  exit 1
fi

