#!/bin/sh

# grab and install TRAMO on linux

# You can supply an installation root, e.g. /usr/local or /opt,
# on the command line.  In any case it goes into a subdir
# named "tramo".

TRAMOROOT=$1
if [ "x$TRAMOROOT" = "x" ] ; then
   echo -n "Path for installation? (default is /opt): " 
   read TRAMOROOT
   if [ "x$TRAMOROOT" = "x" ] ; then
      TRAMOROOT="/opt"
   fi  
fi

echo "Installing in $TRAMOROOT"

SITE="http://www.bde.es/servicio/software/tramo"
PKG="linux.exe"

DEST=$TRAMOROOT/tramo

mkdir -p $DEST
cd $DEST
wget $SITE/$PKG
unzip $PKG && rm $PKG
tar xvfz TSU.TGZ && rm TSU.TGZ
mkdir docs
cd docs
wget $SITE/manualdos.pdf
wget $SITE/updatedos.pdf
wget $SITE/guide.pdf
