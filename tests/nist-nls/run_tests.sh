#!/bin/sh

./nist-nls-test > numeric.out
./nist-nls-test -a > analytic.out
GRETL_USE_QR=1 ./nist-nls-test > qr-numeric.out
GRETL_USE_QR=1 ./nist-nls-test -a > qr-analytic.out

