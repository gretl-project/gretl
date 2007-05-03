#!/bin/sh

# if no fonts are available, do
#   sudo fc-cache
# before running this program


CWD="`(cd \"\`dirname \\\"$0\\\"\`\"; echo $PWD)`"
TOP="`dirname \"$CWD\"`"

echo "TOP=$TOP" > ~/where

if [ -f ~/.profile ] ; then
  . ~/.profile
fi  

export "GRETL_HOME=$TOP/share/gretl/"
export "GRETL_LANGS_DIR=$TOP/share/gretl/gtksourceview"
export "GTK_EXE_TOP=$TOP"
export "GDK_PIXBUF_MODULE_FILE=$TOP/etc/gtk-2.0/gdk-pixbuf.loaders"
export "PANGO_RC_FILE=$TOP/etc/pangorc"
export "GTK_IM_MODULE_FILE=$TOP/etc/gtk-2.0/gtk.immodules"
export "DYLD_LIBRARY_PATH=$TOP/lib"
export "XDG_DATA_DIRS=$TOP/share"
export "XDG_DATA_HOME=$TOP/share"

# location of gnuplot help file
export "GNUHELP=$TOP/share/gnuplot/4.0/gnuplot.gih"
# location of gnuplot X11 driver
export "GNUPLOT_DRIVER_DIR=$TOP/libexec/gnuplot/4.0"
# path for fonts for GD fonts (check this)
export "GDFONTPATH=$TOP/fonts:/usr/X11R6/lib/X11/fonts/TTF"
# default font for gnuplot PNG
export "GNUPLOT_DEFAULT_GDFONT=$TOP/fonts/Vera.ttf"

export "PATH=$CWD:$PATH"

# echo "pwd is `pwd`" >>~/where

cd $CWD
exec "$CWD/gretl_x11" "$@"
