#!/bin/bash

# the directory above Gretl.app
TOPDIR=/Users/allin/dist

HERE=`pwd`
KB=`du -ks ~/dist | awk '{ print $1 }'`
KB=$((KB+640))
hdiutil create -size ${KB}k tmp.dmg -layout NONE
MYDEV=`hdid -nomount tmp.dmg`
sudo newfs_hfs -v gretl $MYDEV
hdiutil eject $MYDEV
hdid tmp.dmg
cd $TOPDIR && \
cp -a Gretl.app /Volumes/gretl && \
cp -a README.pdf /Volumes/gretl
cd $HERE
hdiutil eject $MYDEV
hdiutil convert -format UDZO tmp.dmg -o gretl.dmg && rm tmp.dmg
