#!/bin/sh

TOPDIR="/Users/allin/dist"
PKGDIR="$TOPDIR/Gretl.app/Contents/Resources"

cd $PKGDIR

otool -L ./bin/gretl_x11 | awk '{ print $1 }' | grep /sw/lib > liblist
otool -L ./bin/gretlcli | awk '{ print $1 }' | grep /sw/lib >> liblist
otool -L ./bin/gnuplot | awk '{ print $1 }' | grep /sw/lib >> liblist
 
for f in `cat liblist | uniq` ; do
  cp $f ./lib
done

rm liblist

exit 0

# The following is only a vague approximation to what needs
# to be done: these actions will grab a lot of files that are
# not really needed; and various configuration files need to be
# edited.

cd /sw && tar cvfz ~/gtklib.tgz lib/gtk-2.0 lib/pango && \
cd $PKGDIR && tar xvfz ~/gtklib.tgz && \
rm ~/gtklib.tgz

cd /sw && tar cvfz ~/gtketc.tgz etc/gtk-2.0 etc/pango && \
cd $PKGDIR && tar xvfz ~/gtketc.tgz && \
rm ~/gtketc.tgz
