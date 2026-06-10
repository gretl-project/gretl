#!/bin/bash

# get this script's directory
BINDIR=`dirname "$0"`

# determine the Resources directory
RES="$BINDIR"/..

export GRETL_HOME="$RES"/share/gretl/

# GTK2-specific info (harmless for GTK3?)
export PANGO_SYSCONFDIR="$RES"/etc
export PANGO_RC_FILE="$RES"/etc/pango/pangorc
export PANGO_LIBDIR="$RES"/lib

export G_FILENAME_ENCODING=UTF-8

# location of gnuplot help file
export GNUHELP="$RES"/share/gnuplot/6.0/gnuplot.gih
# location of PostScript resources
export GNUPLOT_PS_DIR="$RES"/share/gnuplot/6.0/PostScript

export "PATH=$BINDIR:$PATH"

# Strip out the argument added by the OS, if any
if /bin/expr "x$1" : '^x-psn_' > /dev/null; then
    shift 1
fi

exec "$BINDIR"/gretlcli "$@"
