top_builddir = ..

ifeq ($(prefix),)
  prefix = /opt/esl
endif
ifeq ($(bindir),)
  bindir = $(prefix)/bin
endif
ifeq ($(CC),)
  CC = gcc
endif
ifeq ($(CFLAGS),)
  CFLAGS = -g -Wall -pedantic
endif
ifeq ($(MAKE),)
  MAKE = make
endif
ifeq ($(INSTALL),)
  INSTALL = /bin/install -c
endif
ifeq ($(INSTALL_PROGRAM),)
  INSTALL_PROGRAM = $(INSTALL) -m 755
endif
ifeq ($(INSTALL_DATA),)
  INSTALL_DATA = $(INSTALL) -m 644
endif
ifeq ($(LN),)
  LN = ln -sf
endif

# Directories
srcdir  = ./src
VPATH = ./src
libdir = $(prefix)/lib
includedir = $(prefix)/include/gretl
share = $(prefix)/share/gretl/

# General compile options and libs:
LIBGRETL  = libgretl.la
LIBS = -lm -L/usr/local/lib -lz

#### End of system configuration section. ####

SRCS = estimate.c compare.c utils.c errors.c pvalues.c dataio.c \
	generate.c strutils.c describe.c printout.c texprint.c var.c \
	interact.c graphing.c monte_carlo.c nonparam.c discrete.c \
	subsample.c
HDRS = libgretl.h commands.h estimate.h utils.h errors.h pvalues.h \
	compare.h dataio.h generate.h strutils.h describe.h printout.h \
	version.h texprint.h var.h interact.h graphing.h monte_carlo.h \
	nonparam.h discrete.h subsample.h internal.h
OBJS = estimate.o compare.o utils.o errors.o pvalues.o dataio.o \
	generate.o strutils.o describe.o printout.o texprint.o var.o \
	interact.o graphing.o monte_carlo.o nonparam.o discrete.o \
	subsample.o
LOBJS = estimate.lo compare.lo utils.lo errors.lo pvalues.lo dataio.lo \
	generate.lo strutils.lo describe.lo printout.lo texprint.lo var.lo \
	interact.lo graphing.lo monte_carlo.lo nonparam.lo discrete.lo \
	subsample.lo

COMPILE = ../libtool --mode=compile $(CC) -c -I.. -I./src $(CFLAGS)

.c.o:        
	$(COMPILE) $<

$(LIBGRETL): $(OBJS) $(LOBJS) ../Makefile
	../libtool --mode=link $(CC) -o $(LIBGRETL) $(LOBJS) \
	-rpath $(prefix)/lib $(LIBS) -version-info 6:0:0

$(LOBJS): Makefile

$(OBJS): Makefile

.PHONY: 

install: $(LIBGRETL) installdirs
	../libtool $(INSTALL_PROGRAM) $(LIBGRETL) $(libdir)/$(LIBGRETL)
	$(INSTALL_DATA) $(srcdir)/*.h $(includedir)
	$(INSTALL_PROGRAM) ../gretl-config $(bindir)

install-strip: $(LIBGRETL) installdirs
	../libtool $(INSTALL_PROGRAM) -s $(LIBGRETL) $(libdir)/$(LIBGRETL)
	$(INSTALL_DATA) $(srcdir)/*.h $(includedir)
	$(INSTALL_PROGRAM) ../gretl-config $(bindir)

installdirs:
	../tools/mkinstalldirs $(libdir) $(includedir)

clean:
	rm -f *.o *.lo $(LIBGRETL) core 
	rm -rf ./.libs

tags:
	cd $(srcdir)/../.. && etags lib/src/*.c cli/*.c gui/*.c


# DO NOT DELETE

src/compare.o: src/libgretl.h src/commands.h src/errors.h src/estimate.h
src/compare.o: src/generate.h src/utils.h src/compare.h src/pvalues.h
src/compare.o: src/dataio.h src/strutils.h src/describe.h src/printout.h
src/compare.o: src/texprint.h src/var.h src/interact.h src/graphing.h
src/compare.o: src/monte_carlo.h src/nonparam.h src/discrete.h
src/compare.o: src/subsample.h src/internal.h
src/dataio.o: src/libgretl.h src/commands.h src/errors.h src/estimate.h
src/dataio.o: src/generate.h src/utils.h src/compare.h src/pvalues.h
src/dataio.o: src/dataio.h src/strutils.h src/describe.h src/printout.h
src/dataio.o: src/texprint.h src/var.h src/interact.h src/graphing.h
src/dataio.o: src/monte_carlo.h src/nonparam.h src/discrete.h src/subsample.h
src/dataio.o: src/internal.h
src/describe.o: src/libgretl.h src/commands.h src/errors.h src/estimate.h
src/describe.o: src/generate.h src/utils.h src/compare.h src/pvalues.h
src/describe.o: src/dataio.h src/strutils.h src/describe.h src/printout.h
src/describe.o: src/texprint.h src/var.h src/interact.h src/graphing.h
src/describe.o: src/monte_carlo.h src/nonparam.h src/discrete.h
src/describe.o: src/subsample.h src/internal.h
src/discrete.o: src/libgretl.h src/commands.h src/errors.h src/estimate.h
src/discrete.o: src/generate.h src/utils.h src/compare.h src/pvalues.h
src/discrete.o: src/dataio.h src/strutils.h src/describe.h src/printout.h
src/discrete.o: src/texprint.h src/var.h src/interact.h src/graphing.h
src/discrete.o: src/monte_carlo.h src/nonparam.h src/discrete.h
src/discrete.o: src/subsample.h src/internal.h
src/errors.o: src/libgretl.h src/commands.h src/errors.h src/estimate.h
src/errors.o: src/generate.h src/utils.h src/compare.h src/pvalues.h
src/errors.o: src/dataio.h src/strutils.h src/describe.h src/printout.h
src/errors.o: src/texprint.h src/var.h src/interact.h src/graphing.h
src/errors.o: src/monte_carlo.h src/nonparam.h src/discrete.h src/subsample.h
src/estimate.o: src/libgretl.h src/commands.h src/errors.h src/estimate.h
src/estimate.o: src/generate.h src/utils.h src/compare.h src/pvalues.h
src/estimate.o: src/dataio.h src/strutils.h src/describe.h src/printout.h
src/estimate.o: src/texprint.h src/var.h src/interact.h src/graphing.h
src/estimate.o: src/monte_carlo.h src/nonparam.h src/discrete.h
src/estimate.o: src/subsample.h src/internal.h
src/generate.o: src/libgretl.h src/commands.h src/errors.h src/estimate.h
src/generate.o: src/generate.h src/utils.h src/compare.h src/pvalues.h
src/generate.o: src/dataio.h src/strutils.h src/describe.h src/printout.h
src/generate.o: src/texprint.h src/var.h src/interact.h src/graphing.h
src/generate.o: src/monte_carlo.h src/nonparam.h src/discrete.h
src/generate.o: src/subsample.h src/internal.h
src/graphing.o: src/libgretl.h src/commands.h src/errors.h src/estimate.h
src/graphing.o: src/generate.h src/utils.h src/compare.h src/pvalues.h
src/graphing.o: src/dataio.h src/strutils.h src/describe.h src/printout.h
src/graphing.o: src/texprint.h src/var.h src/interact.h src/graphing.h
src/graphing.o: src/monte_carlo.h src/nonparam.h src/discrete.h
src/graphing.o: src/subsample.h src/internal.h src/graphing.h
src/interact.o: src/libgretl.h src/commands.h src/errors.h src/estimate.h
src/interact.o: src/generate.h src/utils.h src/compare.h src/pvalues.h
src/interact.o: src/dataio.h src/strutils.h src/describe.h src/printout.h
src/interact.o: src/texprint.h src/var.h src/interact.h src/graphing.h
src/interact.o: src/monte_carlo.h src/nonparam.h src/discrete.h
src/interact.o: src/subsample.h src/internal.h src/cmdlist.h
src/monte_carlo.o: src/libgretl.h src/commands.h src/errors.h src/estimate.h
src/monte_carlo.o: src/generate.h src/utils.h src/compare.h src/pvalues.h
src/monte_carlo.o: src/dataio.h src/strutils.h src/describe.h src/printout.h
src/monte_carlo.o: src/texprint.h src/var.h src/interact.h src/graphing.h
src/monte_carlo.o: src/monte_carlo.h src/nonparam.h src/discrete.h
src/monte_carlo.o: src/subsample.h src/internal.h
src/nonparam.o: src/libgretl.h src/commands.h src/errors.h src/estimate.h
src/nonparam.o: src/generate.h src/utils.h src/compare.h src/pvalues.h
src/nonparam.o: src/dataio.h src/strutils.h src/describe.h src/printout.h
src/nonparam.o: src/texprint.h src/var.h src/interact.h src/graphing.h
src/nonparam.o: src/monte_carlo.h src/nonparam.h src/discrete.h
src/nonparam.o: src/subsample.h src/internal.h
src/printout.o: src/libgretl.h src/commands.h src/errors.h src/estimate.h
src/printout.o: src/generate.h src/utils.h src/compare.h src/pvalues.h
src/printout.o: src/dataio.h src/strutils.h src/describe.h src/printout.h
src/printout.o: src/texprint.h src/var.h src/interact.h src/graphing.h
src/printout.o: src/monte_carlo.h src/nonparam.h src/discrete.h
src/printout.o: src/subsample.h src/internal.h src/version.h
src/pvalues.o: src/libgretl.h src/commands.h src/errors.h src/estimate.h
src/pvalues.o: src/generate.h src/utils.h src/compare.h src/pvalues.h
src/pvalues.o: src/dataio.h src/strutils.h src/describe.h src/printout.h
src/pvalues.o: src/texprint.h src/var.h src/interact.h src/graphing.h
src/pvalues.o: src/monte_carlo.h src/nonparam.h src/discrete.h
src/pvalues.o: src/subsample.h
src/strutils.o: src/libgretl.h src/commands.h src/errors.h src/estimate.h
src/strutils.o: src/generate.h src/utils.h src/compare.h src/pvalues.h
src/strutils.o: src/dataio.h src/strutils.h src/describe.h src/printout.h
src/strutils.o: src/texprint.h src/var.h src/interact.h src/graphing.h
src/strutils.o: src/monte_carlo.h src/nonparam.h src/discrete.h
src/strutils.o: src/subsample.h
src/subsample.o: src/libgretl.h src/commands.h src/errors.h src/estimate.h
src/subsample.o: src/generate.h src/utils.h src/compare.h src/pvalues.h
src/subsample.o: src/dataio.h src/strutils.h src/describe.h src/printout.h
src/subsample.o: src/texprint.h src/var.h src/interact.h src/graphing.h
src/subsample.o: src/monte_carlo.h src/nonparam.h src/discrete.h
src/subsample.o: src/subsample.h src/internal.h
src/texprint.o: src/libgretl.h src/commands.h src/errors.h src/estimate.h
src/texprint.o: src/generate.h src/utils.h src/compare.h src/pvalues.h
src/texprint.o: src/dataio.h src/strutils.h src/describe.h src/printout.h
src/texprint.o: src/texprint.h src/var.h src/interact.h src/graphing.h
src/texprint.o: src/monte_carlo.h src/nonparam.h src/discrete.h
src/texprint.o: src/subsample.h
src/utils.o: ../config.h src/libgretl.h src/commands.h src/errors.h
src/utils.o: src/estimate.h src/generate.h src/utils.h src/compare.h
src/utils.o: src/pvalues.h src/dataio.h src/strutils.h src/describe.h
src/utils.o: src/printout.h src/texprint.h src/var.h src/interact.h
src/utils.o: src/graphing.h src/monte_carlo.h src/nonparam.h src/discrete.h
src/utils.o: src/subsample.h src/internal.h
src/var.o: src/libgretl.h src/commands.h src/errors.h src/estimate.h
src/var.o: src/generate.h src/utils.h src/compare.h src/pvalues.h
src/var.o: src/dataio.h src/strutils.h src/describe.h src/printout.h
src/var.o: src/texprint.h src/var.h src/interact.h src/graphing.h
src/var.o: src/monte_carlo.h src/nonparam.h src/discrete.h src/subsample.h
src/var.o: src/internal.h
src/win32.o: src/libgretl.h src/commands.h src/errors.h src/estimate.h
src/win32.o: src/generate.h src/utils.h src/compare.h src/pvalues.h
src/win32.o: src/dataio.h src/strutils.h src/describe.h src/printout.h
src/win32.o: src/texprint.h src/var.h src/interact.h src/graphing.h
src/win32.o: src/monte_carlo.h src/nonparam.h src/discrete.h src/subsample.h
