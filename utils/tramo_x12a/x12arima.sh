#!/bin/sh

# grab and build x12a for Linux

# You can supply an installation root, e.g. /usr/local or /opt,
# otherwise it goes right here.  In any case it goes into a subdir
# named "x12arima".

X12ROOT=$1
if [ "x$X12ROOT" = "x" ] ; then
   X12ROOT=.
fi

SITE="http://www.census.gov/ts/x12a/final/unix"
SRC="omegasrc.tar.Z"
EXAMPLES="finexam.tar.Z"
PSDOCS="omegaps.tar.Z"

THISDIR=`pwd`
DEST=$X12ROOT/x12arima

# create a directory to hold the distribution and cd to it
mkdir -p $DEST
mkdir -p $DEST/src
cd $DEST

# grab the relevant archives and extract
wget $SITE/$SRC
cd $DEST/src && tar xvfz $DEST/$SRC && rm $DEST/$SRC && cd $DEST
wget $SITE/$EXAMPLES 
tar xvfz $EXAMPLES && rm $EXAMPLES
wget $SITE/$PSDOCS
tar xvfz $PSDOCS && rm $PSDOCS

# move the aux files up from src
mv ./src/x12a.mdl ./src/test.spc $DEST

# patch the Makefile for Linux
patch -p0 < $THISDIR/x12a_make.diff

# remove DOS file termination chars
echo "Fixing DOS-terminated source files..."
for f in ./src/*.f ./src/*.var ./src/*.prm ./src/*.cmn ; do
   cat $f | tr -d '\032' > $f.tmp && mv $f.tmp $f
done

# build the program, then compress the sources
make -C src 
strip ./x12a
make -C src clean
tar cvfz src.tgz src && rm -rf src




