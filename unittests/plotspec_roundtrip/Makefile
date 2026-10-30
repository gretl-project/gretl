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
	$(shell pkg-config --libs glib-2.0) -lm
CFLAGS = -g -O1 -Wall -W -Wno-unused-parameter

PLOTS = plots

SAN = -g -O1 -fsanitize=address,undefined -fno-omit-frame-pointer

all: roundtrip

roundtrip: roundtrip.c $(TOP)/lib/src/plotspec.h
	$(CC) $(CFLAGS) $(INC) roundtrip.c -o $@ $(LIBS)

$(PLOTS)/.stamp: gen_plots.inp
	rm -rf $(PLOTS) && mkdir $(PLOTS)
	cd $(PLOTS) && DISPLAY= $(abspath $(GRETLCLI)) -b -q -e $(abspath gen_plots.inp) > gen.log 2>&1 \
	  || { tail gen.log; exit 1; }
	touch $@

# The reader and plotspec.c are compiled into the test binary itself
# (so that they are instrumented, unlike the copies in libgretl); the
# reader is fed damaged files.
roundtrip-asan: roundtrip.c $(TOP)/lib/src/plotspec_read.c $(TOP)/lib/src/plotspec.c $(TOP)/lib/src/plotspec.h
	$(CC) $(SAN) -Wall -W -Wno-unused-parameter $(INC) roundtrip.c \
	  $(TOP)/lib/src/plotspec_read.c $(TOP)/lib/src/plotspec.c -o $@ $(LIBS)

check-asan: roundtrip-asan $(PLOTS)/.stamp
	./roundtrip-asan $(PLOTS)/*.plt
	./roundtrip-asan -f $(PLOTS)/*.plt

check: roundtrip $(PLOTS)/.stamp
	./roundtrip $(PLOTS)/*.plt

clean:
	rm -rf $(PLOTS) roundtrip roundtrip-asan

.PHONY: all check check-asan clean
