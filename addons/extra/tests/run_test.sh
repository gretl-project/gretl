#!/bin/bash
set -e

DIR=$(dirname $(realpath "$0")) 	# locate folder where this sh-script is located in
SCRIPT_1="test_extra.inp"
SCRIPT_2="test_extra_matrix.inp"
SCRIPT_3="test_add_outofsample.inp"
PROJECT="extra"

cd $DIR
echo "Switched to ${DIR}"

gretlcli -b -e -q ${SCRIPT_1}
if [ $? -eq 0 ]
then
  echo "Success: All tests passed for '${SCRIPT_1}'."
else
  echo "Failure: Tests not passed for '${SCRIPT_1}'." >&2
  exit 1
fi

gretlcli -b -e -q ${SCRIPT_2}
if [ $? -eq 0 ]
then
  echo "Success: All tests passed for '${SCRIPT_2}'."
else
  echo "Failure: Tests not passed for '${SCRIPT_2}'." >&2
  exit 1
fi

gretlcli -b -e -q ${SCRIPT_3}
if [ $? -eq 0 ]
then
  echo "Success: All tests passed for '${SCRIPT_3}'."
else
  echo "Failure: Tests not passed for '${SCRIPT_3}'." >&2
  exit 1
fi


echo "Success: All tests passed."
exit 0
