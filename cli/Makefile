top_builddir = ..

ifeq ($(prefix),)
  prefix = /opt/esl
  CLI_LIBS = -lreadline -ltermcap
  RL_OBJ = complete.o
endif
ifeq ($(CC),)
  CC = gcc
endif
ifeq ($(CFLAGS),)
  CFLAGS = -g -Wall
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
ifeq ($(LN),)
  LN = ln -sf
endif

# General compile options and libs:
LIBS = $(CLI_LIBS)
INCLUDES = -I.

# Directories
srcdir  = .
VPATH   = .
tooldir = ../tools
bindir = $(prefix)/bin

#### End of system configuration section. ####

CLI = gretlcli

COMPILE = $(CC) -c $(CFLAGS) $(INCLUDES) $(DEFINES) 

.c.o:
	$(COMPILE) $<

$(CLI): gretlcli.o $(RL_OBJ) Makefile ../Makefile 
	../libtool --mode=link $(CC) -o $(CLI) gretlcli.o \
	$(RL_OBJ) ../lib/libgretl.la $(CLI_LIBS) 

install: $(CLI) installdirs
	../libtool $(INSTALL_PROGRAM) $(CLI) $(bindir)/$(CLI)

install-strip: $(CLI) installdirs
	../libtool $(INSTALL_PROGRAM) -s $(CLI) $(bindir)/$(CLI)

installdirs: $(tooldir)/mkinstalldirs
	$(tooldir)/mkinstalldirs $(bindir)

clean:
	rm -f *.o $(CLI) core 
	rm -rf ./.libs

# DO NOT DELETE

complete.o: ../config.h ../lib/src/libgretl.h ../lib/src/commands.h
complete.o: ../lib/src/errors.h ../lib/src/estimate.h ../lib/src/generate.h
complete.o: ../lib/src/utils.h ../lib/src/compare.h ../lib/src/pvalues.h
complete.o: ../lib/src/dataio.h ../lib/src/strutils.h ../lib/src/describe.h
complete.o: ../lib/src/printout.h ../lib/src/texprint.h ../lib/src/var.h
complete.o: ../lib/src/interact.h ../lib/src/graphing.h
complete.o: ../lib/src/monte_carlo.h ../lib/src/nonparam.h
complete.o: ../lib/src/discrete.h ../lib/src/subsample.h
gretlcli.o: ../config.h ../lib/src/libgretl.h ../lib/src/commands.h
gretlcli.o: ../lib/src/errors.h ../lib/src/estimate.h ../lib/src/generate.h
gretlcli.o: ../lib/src/utils.h ../lib/src/compare.h ../lib/src/pvalues.h
gretlcli.o: ../lib/src/dataio.h ../lib/src/strutils.h ../lib/src/describe.h
gretlcli.o: ../lib/src/printout.h ../lib/src/texprint.h ../lib/src/var.h
gretlcli.o: ../lib/src/interact.h ../lib/src/graphing.h
gretlcli.o: ../lib/src/monte_carlo.h ../lib/src/nonparam.h
gretlcli.o: ../lib/src/discrete.h ../lib/src/subsample.h common.c
