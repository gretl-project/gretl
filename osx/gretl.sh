#!/bin/sh

# if no fonts are available, do
#   sudo fc-cache
# before running this program

PREFIX=`pwd | sed s+/bin$++`

# echo "PREFIX=$PREFIX" > ~/where

if [ -f ~/.profile ] ; then
  . ~/.profile
fi  

GRETL_HOME="$PREFIX/share/gretl/"
export GRETL_HOME
GTK_EXE_PREFIX="$PREFIX"
export GTK_EXE_PREFIX
GDK_PIXBUF_MODULE_FILE="$PREFIX/etc/gtk-2.0/gdk-pixbuf.loaders"
export GDK_PIXBUF_MODULE_FILE
PANGO_RC_FILE="$PREFIX/etc/pangorc"
export PANGO_RC_FILE
GTK_IM_MODULE_FILE="$PREFIX/etc/gtk-2.0/gtk.immodules"
export GTK_IM_MODULE_FILE
DYLD_LIBRARY_PATH="$PREFIX/lib"
export DYLD_LIBRARY_PATH

# location of gnuplot help file
GNUHELP=$PREFIX/share/gnuplot/4.0/gnuplot.gih
export GNUHELP
# location of gnuplot X11 driver
GNUPLOT_DRIVER_DIR=$PREFIX/libexec/gnuplot/4.0
export GNUPLOT_DRIVER_DIR
# path for fonts for GD fonts (check this)
GDFONTPATH=$PREFIX/fonts:/usr/X11R6/lib/X11/fonts/TTF
export GDFONTPATH
# default font for gnuplot PNG
GNUPLOT_DEFAULT_GDFONT=$PREFIX/fonts/Vera.ttf
export GNUPLOT_DEFAULT_GDFONT

export PATH=$PREFIX/bin:$PATH
exec "$PREFIX/bin/gretl_x11" "$@"
