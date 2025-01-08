#!/bin/bash

# get this script's directory
BINDIR=`dirname "$0"`

# determine the Resources directory
RES="$BINDIR"/..

export G_FILENAME_ENCODING=UTF-8

# location of gnuplot help file
export GNUHELP="$RES"/share/gnuplot/6.0/gnuplot.gih
# location of PostScript resources
export GNUPLOT_PS_DIR="$RES"/share/gnuplot/6.0/PostScript
# preferred terminal
export GNUTERM=wxt

# Strip out the argument added by the OS, if any
if /bin/expr "x$1" : '^x-psn_' > /dev/null; then
    shift 1
fi

exec "$BINDIR"/gnuplot "$@"
