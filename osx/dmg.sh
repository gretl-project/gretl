#!/bin/bash

KB=`du -ks ~/dist | awk '{ print $1 }'`
KB=$((KB+640))
hdiutil create -size ${KB}k tmp.dmg -layout NONE
MYDEV=`hdid -nomount tmp.dmg`
sudo newfs_hfs -v gretl $MYDEV
hdiutil eject $MYDEV
hdid tmp.dmg
# volume "gretl" should now appear in Finder
# copy files into it, then eject, then...

# hdiutil convert -format UDZO tmp.dmg -o gretl.dmg && rm tmp.dmg
