#!/bin/sh

if [ "x$1" = "x" ] ; then
   echo "Please specify a chapter number for NBER macrohistory data."
   exit 1
fi  

CH=$1
if [ $CH -lt 10 ] && ! echo $CH | grep '0' >/dev/null ; then
   CH="0$CH"
fi

echo "getting data for chapter $CH"

wget -r -l2 -A db -X icons http://www.nber.org/databases/macrohistory/data/$1/

mv www.nber.org/databases/macrohistory/data/$1 www.nber.org
rm -rf www.nber.org/databases
rm -f www.nber.org/robots.txt


