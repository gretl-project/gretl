#!/bin/sh

MANUAL=$GTK_EXE_PREFIX/manual.pdf
ASF=`echo $MANUAL | sed s+/++ | sed s+/+:+g`
osascript -e "tell app \"Finder\" 
   open the file \"$ASF\"
end tell"

