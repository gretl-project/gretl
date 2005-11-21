#!/bin/sh

if [ "x$1" = "xuguide" ] ; then
   MANUAL=$GTK_EXE_PREFIX/gretl-guide.pdf 
else
   MANUAL=$GTK_EXE_PREFIX/gretl-ref.pdf  
fi

ASF=`echo $MANUAL | sed s+/++ | sed s+/+:+g`
osascript -e "tell app \"Finder\" 
   open the file \"$ASF\"
end tell"

