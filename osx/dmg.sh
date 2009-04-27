#!/bin/bash

# The directory above Gretl.app
TOPDIR="/Users/allin/dist"

# this has to be GNU cp -- may be in /sw/bin
CP=/usr/local/bin/cp

rm -f gretl.dmg tmp.dmg

HERE=`pwd`
KB=`du -ks ${TOPDIR} | awk '{ print $1 }'`
KB=$((KB+1024))
hdiutil create -size ${KB}k tmp.dmg -layout NONE
MYDEV=`hdid -nomount tmp.dmg`
sudo newfs_hfs -v gretl $MYDEV
hdiutil eject $MYDEV
hdid tmp.dmg
cd ${TOPDIR} && \
$CP -a Gretl.app /Volumes/gretl && \
$CP -a README.pdf /Volumes/gretl
cd $HERE
hdiutil eject $MYDEV
hdiutil convert -format UDZO tmp.dmg -o gretl.dmg && rm tmp.dmg

