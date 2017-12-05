#!/bin/sh

TARG=$1

if [ "x$TARG" = "x" ] ; then
   echo "Please give the name of a file to be adjusted"
   exit 1
fi

sed -i -e 's/\r*$/\r/' $TARG
