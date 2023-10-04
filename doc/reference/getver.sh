#!/bin/sh

TOPSRC="$1"

# the program version
GRETL_VERSION=`grep GRETL_VERSION ${TOPSRC}/lib/src/version.h | \
 awk '{print $NF}' | sed s/-git// | sed s/\"//g`
# the library major version
LIB_VERSION=`grep LIBGRETL_CURRENT ${TOPSRC}/lib/src/version.h | \
 awk '{print $NF}'`
 
if [ "x$GRETL_VERSION" = "x" ] ; then
   exit 1
fi

cat ${TOPSRC}/doc/reference/libgretl-docs.xml.in | \
 sed s/LIB_VERSION/${LIB_VERSION}/ | \
 sed s/GRETL_VERSION/${GRETL_VERSION}/ > libgretl-docs.xml

