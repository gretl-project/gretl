#!/bin/sh

GRETLCLI="$1"

export GRETL_INCLUDE_PATH=`pwd`

rm -f test.out

for f in examples/*.inp ; do
  if $GRETLCLI -b $f >> test.out ; then 
     echo "$f: OK"
  elif [ `basename $f` = "example7.inp" ] ; then
     # requires R's fGarch package
     echo "$f: skipped"
  else
     echo "$f: failed"
     break
  fi
done
