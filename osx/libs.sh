#!/bin/sh

VERSION="1.2.6"
PKGDIR="/Users/allin/gretl-$VERSION"

otool -L ./bin/gretl_x11 | awk '{ print $1 }' | grep /sw/lib \
 > liblist
 
for f in `cat liblist` ; do
  cp $f ./lib
done

cd /sw && tar cvfz ~/lib.tgz lib/gtk-2.0 lib/pango && \
cd $PKGDIR && tar xvfz ~/lib.tgz && \
rm ~/lib.tgz

cd /sw && tar cvfz ~/etc.tgz etc/gtk-2.0 etc/pango && \
cd $PKGDIR && tar xvfz ~/etc.tgz && \
rm ~/etc.tgz
