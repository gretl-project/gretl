#!/bin/sh

DSTR=`date +%Y-%m-%d`
OUT="out.$DSTR"
ERRS="errlog.$DSTR"

./nist-nls-test "$@" >$OUT 2>$ERRS

grep ERROR $ERRS
grep 'failed to conv' $OUT
