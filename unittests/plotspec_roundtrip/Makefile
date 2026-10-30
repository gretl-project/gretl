# Round-trip test for the plotspec reader/writer: generate a set of
# gnuplot scripts using gretlcli, then check that each survives a
# read -> write -> read -> write cycle (see roundtrip.c).
#
#   make check
#
# Assumes an in-tree layout with libgretl built in $(BUILD).

TOP = ../..
BUILD = $(TOP)/build
GRETLCLI = $(BUILD)/cli/gretlcli
LIBDIR = $(BUILD)/lib/.libs
INC = -I$(TOP)/lib/src -I$(BUILD) $(shell pkg-config --cflags glib-2.0 libxml-2.0)
LIBS = -L$(LIBDIR) -Wl,-rpath,$(abspath $(LIBDIR)) -lgretl-1.0 \
	$(shell pkg-config --libs glib-2.0)
CFLAGS = -g -O1 -Wall -W -Wno-unused-parameter

PLOTS = plots

# Known round-trip failures. Remove an entry once the underlying
# problem is fixed (the test then reports XPASS).
#
# freq_normal.plt: the "set label ... at graph X, graph Y front" lines
# carrying the normality test result are dropped. In
# parse_gp_set_line() the result of parse_label_line() is ignored, so
# when the label position is not a pair of plain numbers the line is
# treated as handled by maybe_promote_literal_lines() and removed
# from the literal lines, although nothing was stored.
XFAIL = -x freq_normal.plt

all: roundtrip

roundtrip: roundtrip.c $(TOP)/lib/src/plotspec.h
	$(CC) $(CFLAGS) $(INC) roundtrip.c -o $@ $(LIBS)

$(PLOTS)/.stamp: gen_plots.inp
	rm -rf $(PLOTS) && mkdir $(PLOTS)
	cd $(PLOTS) && DISPLAY= $(abspath $(GRETLCLI)) -b -q -e $(abspath gen_plots.inp) > gen.log 2>&1 \
	  || { tail gen.log; exit 1; }
	touch $@

check: roundtrip $(PLOTS)/.stamp
	./roundtrip $(XFAIL) $(PLOTS)/*.plt

clean:
	rm -rf $(PLOTS) roundtrip

.PHONY: all check clean
