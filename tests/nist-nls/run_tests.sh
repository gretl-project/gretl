#!/bin/sh

./nist-nls-test "$@" 2>errlog

grep ERROR errlog
