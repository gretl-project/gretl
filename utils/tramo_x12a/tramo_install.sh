#!/bin/sh

# grab and install TRAMO on linux

# You can supply an installation root, e.g. /usr/local or /opt,
# otherwise it goes right here.  In any case it goes into a subdir
# named "tramo".

TRAMOROOT=$1
if [ "x$TRAMOROOT" = "x" ] ; then
   TRAMOROOT=.
fi

SITE="http://www.bde.es/servicio/software/tramo"
PKG="linux.exe"

DEST=$TRAMOROOT/tramo

mkdir -p $DEST
cd $DEST
wget $SITE/$PKG
unzip $PKG && rm $PKG
tar xvfz TSU.TGZ && rm TSU.TGZ
