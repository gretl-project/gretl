# Makefile for cross-compiling gretl for win32 on Linux,
# using mingw compiler

include config.mk

# directories
topsrc = ..
libdir  = $(topsrc)/lib/src
extradir = ../gui2/gtkextra-lite
plugindir = ../plugin
imports = ./import-libs

# tools
CC = $(MGW_PREFIX)-gcc -Wall -O2 -mms-bitfields
AS = $(MGW_PREFIX)-as
DLLWRAP = $(MGW_PREFIX)-dllwrap

# libraries
GTKLIBS := -L$(imports) $(shell PKG_CONFIG_PATH=$(PKG_CONFIG_PATH) \
                          pkg-config --libs gtk+-win32-2.0)
GLIBLIB := -L$(imports) $(shell PKG_CONFIG_PATH=$(PKG_CONFIG_PATH) \
                          pkg-config --libs glib-2.0)
ZLIB_LIB = -lz.dll
CLI_LIBS = -L$(imports) -lreadline.dll
GUI_LIBS = $(GTKLIBS) -lgtkextra-lite -lwsock32 #-lhhctrl
GRETLLIB = -L$(imports) -lgretl -lintl.dll
LIBS = -lm $(GRETLLIB) -ladvapi32

# flags/defines
LDFLAGS = -s 
CFLAGS := -DOS_WIN32 -I. -I$(libdir) -I$(extradir) -I$(plugindir) \
$(shell PKG_CONFIG_PATH=$(PKG_CONFIG_PATH) pkg-config --cflags gtk+-win32-2.0)

# source paths
vpath %.c $(libdir)
vpath %.c $(topsrc)/cephes
vpath %.c $(topsrc)/cli
vpath %.c $(topsrc)/gui2
vpath %.c $(extradir)
vpath %.c $(topsrc)/plugin
vpath %.c $(topsrc)/plugin/libole2
vpath %.c $(topsrc)/tests
vpath %.h $(libdir):$(topsrc)/cli:$(topsrc)/gui2:$(extradir):$(topsrc)/plugin:.
vpath %.h $(topsrc)/cephes

#### End of system configuration section. ####

CLIPROG = gretlcli.exe
GUIPROG = gretlw32.exe
UPDATER = gretl_updater.exe
CHECKER = nisttest.exe

PLUGINS = plugins/stats_tables.dll plugins/panel_data.dll plugins/longname.dll \
	plugins/progress_bar.dll plugins/gnumeric_import.dll \
	plugins/excel_import.dll plugins/mp_ols.dll plugins/sur.dll plugins/johansen.dll \
	plugins/range-mean.dll plugins/lad.dll plugins/tramo-x12a.dll

WRES    = gretlw32_res.o
RC      = $(MGW_PREFIX)-windres
RCFLAGS = --define __WIN32__ --define __WIN95__ --define MSRC \
	--define __GNUWIN32__
RCOUT   = gretlw32_res.o
RES2COFF = echo gretlw32_res.o

LIBSRC = estimate.c compare.c gretl_utils.c errors.c pvalues.c dataio.c \
       generate.c strutils.c describe.c printout.c modelprint.c texprint.c var.c \
       interact.c graphing.c gretl_win32.c monte_carlo.c nonparam.c discrete.c \
       subsample.c calendar.c system.c
LIBHDR = libgretl.h commands.h estimate.h gretl_utils.h errors.h pvalues.h \
       compare.h dataio.h generate.h strutils.h describe.h printout.h modelprint.h \
       version.h texprint.h var.h interact.h graphing.h monte_carlo.h \
       nonparam.h discrete.h subsample.h calendar.h system.h gretl_win32.h
LIBOBJ = $(LIBSRC:.c=.o)

CLISRC = gretlcli.c complete.c 
CLIOBJ = $(CLISRC:.c=.o)

CHECKSRC = nistcheck.c

GUISRC = gretl.c library.c gui_utils.c callbacks.c dialogs.c ssheet.c \
	datafiles.c database.c gpt_control.c session.c webget.c \
	calculator.c fileselect.c guiprint.c boxplots.c series_view.c \
	selector.c settings.c helpfiles.c treeutils.c console.c 
GUIHDR = gretl.h library.h gui_utils.h callbacks.h dialogs.h session.h \
	webget.h guiprint.h selector.h settings.h helpfiles.h \
	treeutils.h build.h
GUIOBJ = $(GUISRC:.c=.o)

HACK_SRC = gtkfontselhack.c font_filter.c
HACK_HDR = gtkfontselhack.h
HACK_OBJ = gtkfontselhack.o

EXTRASRC = gtkpsfont.c gtkplotpc.c gtkplotps.c
EXTRAHDR = gtkplot-lite.h gtkplotps.h gtkplotpc.h gtkpsfont.h
EXTRAOBJ = $(EXTRASRC:.c=.o)

PROBSRC = bdtr.c btdtr.c chdtr.c drand.c expx2.c fdtr.c gamma.c gdtr.c \
	igam.c igami.c incbet.c incbi.c mtherr.c nbdtr.c ndtr.c ndtri.c pdtr.c \
	stdtr.c unity.c polevl.c const.c
PROBHDR = mconf.h
PROBOBJ = $(PROBSRC:.c=.o)

MO = mo/es.mo mo/fr.mo

# targets begin

all: buildstamp dirs prob-dll dll extra-dll plugin $(CLIPROG) $(GUIPROG) \
	$(CHECKER) $(MO)

$(CLIPROG): $(CLIOBJ) 
	$(CC) -s -o $@ $^ $(LIBS) $(CLI_LIBS) $(ZLIB_LIB) -lmingwex

$(GUIPROG): $(GUIOBJ) $(HACK_OBJ) $(WRES) 
	$(CC) -mwindows -s -o $@ $^ $(LIBS) $(GUI_LIBS) -lpng $(ZLIB_LIB) -lmingwex

$(CHECKER): $(CHECKSRC)
	$(CC) -s -I. -I../lib/src -DOS_WIN32 -DUSE_GMP -o $@ $^ $(LIBS) $(ZLIB_LIB)

$(CLIOBJ): $(CLI_SRC) 

$(HACK_OBJ): $(HACK_SRC) $(HACK_HDR)

buildstamp:
	./builddate.pl

$(EXTRAOBJ): $(EXTRAHDR) $(EXTRASRC) 

$(LIBOBJ): $(LIBSRC) $(LIBHDR)

$(PROBOBJ): $(PROBSRC) $(PROBHDR) 

$(WRES): gretlw32.rc 
	$(RC) $(RCFLAGS) gretlw32.rc $(RCOUT)
	$(RES2COFF)

$(UPDATER):
	make -C updater

$(MO):
	make -C mo

dirs:
	mkdir -p dlls 
	mkdir -p plugins 

dist: dirs prob-dll dll extra-dll plugin $(CLIPROG) $(GUIPROG) $(MO) $(UPDATER)
	./mkwindist 

topclean:
	rm -f *.o *.exe

clean:
	rm -f *.o *.def $(CLIPROG) $(GUIPROG) $(CHECKER) dlls/libgretl.dll
	rm -f plugins/* $(PLUGINS)
	rm -f dlls/libgtkextra-lite.dll $(imports)/libgtkextra-lite.a $(imports)/libgtkextra-lite.def
	rm -f gretl.stamp build.h
	rm -f windist/MANIFEST windist/gretl.iss 
	rm -rf windist/gretl windist/gnuplot windist/Output
	make -C updater clean
	make -C mo clean

depend:
	make -f Make.windep depend 2>/dev/null && \
	../tools/win_depend.pl

dll: dlls/libgretl.dll
extra-dll: dlls/libgtkextra-lite.dll
prob-dll: dlls/libprob.dll
plugin: $(PLUGINS)

DLLWRAP_FLAGS = --as=$(AS) --export-all --driver-name $(CC) -s

dlls/libgretl.dll: $(LIBOBJ) dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) \
	--output-def libgretl.def --implib $(imports)/libgretl.a \
	-o $@ $^ -lm -L$(imports) -lxml -lz.dll -lintl.dll -lprob -lgmp -lmingwex $(GLIBLIB)

dlls/libprob.dll: $(PROBOBJ) dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) \
	--output-def libprob.def --implib $(imports)/libprob.a \
	-o $@ $^ 

dlls/libgtkextra-lite.dll: $(EXTRAOBJ) dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) \
	--output-def libgtkextra-lite.def --implib $(imports)/libgtkextra-lite.a \
	-o $@ $^ $(GTKLIBS)

plugins/stats_tables.dll: stats_tables.o dllinit.o 
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB)

plugins/range-mean.dll: range-mean.o dllinit.o 
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB)

plugins/lad.dll: lad.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) 

plugins/panel_data.dll: panel_data.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) $(GTKLIBS)

plugins/tramo-x12a.dll: tramo-x12a.o tramo_options.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) $(GTKLIBS)

plugins/progress_bar.dll: progress_bar.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GTKLIBS) -L$(imports) -lintl.dll

plugins/gnumeric_import.dll: gnumeric_import.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ -L$(imports) -lxml \
	$(GRETLLIB) $(GTKLIBS)

plugins/excel_import.dll: excel_import.o workbook.o ms-ole.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) $(GTKLIBS) -lmoldname

plugins/longname.dll: longname.o 
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ -lkernel32

plugins/mp_ols.dll: mp_ols.o dllinit.o 
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) -lgmp.dll

plugins/sur.dll: sur.o gretl_matrix.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) -llapack -lblas -lf2c

plugins/johansen.dll: johansen.o gretl_matrix.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) -llapack -lblas -lf2c

ms-ole.o: ms-ole.c $(topsrc)/plugin/libole2/ms-ole.h
	$(CC) -c $(CFLAGS) $(GLIBINC) $<

workbook.o: workbook.c $(topsrc)/plugin/libole2/ms-ole.h
	$(CC) -c $(CFLAGS) $(GLIBINC) $<

# DO NOT DELETE

boxplots.o: ../gui2/gretl.h winconfig.h ../lib/src/libgretl.h
boxplots.o: ../lib/src/commands.h ../lib/src/errors.h
boxplots.o: ../lib/src/estimate.h ../lib/src/generate.h
boxplots.o: ../lib/src/gretl_utils.h ../lib/src/compare.h
boxplots.o: ../lib/src/pvalues.h ../lib/src/dataio.h
boxplots.o: ../lib/src/strutils.h ../lib/src/describe.h
boxplots.o: ../lib/src/printout.h ../lib/src/modelprint.h
boxplots.o: ../lib/src/texprint.h ../lib/src/var.h
boxplots.o: ../lib/src/interact.h ../lib/src/graphing.h
boxplots.o: ../lib/src/monte_carlo.h ../lib/src/nonparam.h
boxplots.o: ../lib/src/discrete.h ../lib/src/subsample.h
boxplots.o: ../lib/src/calendar.h ../lib/src/system.h
boxplots.o: ../lib/src/gretl_win32.h ../gui2/gretltypes.h
boxplots.o: ../gui2/gui_utils.h ../gui2/callbacks.h ../gui2/dialogs.h
boxplots.o: ../gui2/library.h ../gui2/session.h ../gui2/settings.h
boxplots.o: ../gui2/helpfiles.h ../gui2/boxplots.h
calculator.o: ../gui2/gretl.h winconfig.h ../lib/src/libgretl.h
calculator.o: ../lib/src/commands.h ../lib/src/errors.h
calculator.o: ../lib/src/estimate.h ../lib/src/generate.h
calculator.o: ../lib/src/gretl_utils.h ../lib/src/compare.h
calculator.o: ../lib/src/pvalues.h ../lib/src/dataio.h
calculator.o: ../lib/src/strutils.h ../lib/src/describe.h
calculator.o: ../lib/src/printout.h ../lib/src/modelprint.h
calculator.o: ../lib/src/texprint.h ../lib/src/var.h
calculator.o: ../lib/src/interact.h ../lib/src/graphing.h
calculator.o: ../lib/src/monte_carlo.h ../lib/src/nonparam.h
calculator.o: ../lib/src/discrete.h ../lib/src/subsample.h
calculator.o: ../lib/src/calendar.h ../lib/src/system.h
calculator.o: ../lib/src/gretl_win32.h ../gui2/gretltypes.h
calculator.o: ../gui2/gui_utils.h ../gui2/callbacks.h
calculator.o: ../gui2/dialogs.h ../gui2/library.h ../gui2/session.h
calculator.o: ../gui2/settings.h ../gui2/helpfiles.h
callbacks.o: ../gui2/gretl.h winconfig.h ../lib/src/libgretl.h
callbacks.o: ../lib/src/commands.h ../lib/src/errors.h
callbacks.o: ../lib/src/estimate.h ../lib/src/generate.h
callbacks.o: ../lib/src/gretl_utils.h ../lib/src/compare.h
callbacks.o: ../lib/src/pvalues.h ../lib/src/dataio.h
callbacks.o: ../lib/src/strutils.h ../lib/src/describe.h
callbacks.o: ../lib/src/printout.h ../lib/src/modelprint.h
callbacks.o: ../lib/src/texprint.h ../lib/src/var.h
callbacks.o: ../lib/src/interact.h ../lib/src/graphing.h
callbacks.o: ../lib/src/monte_carlo.h ../lib/src/nonparam.h
callbacks.o: ../lib/src/discrete.h ../lib/src/subsample.h
callbacks.o: ../lib/src/calendar.h ../lib/src/system.h
callbacks.o: ../lib/src/gretl_win32.h ../gui2/gretltypes.h
callbacks.o: ../gui2/gui_utils.h ../gui2/callbacks.h
callbacks.o: ../gui2/dialogs.h ../gui2/library.h ../gui2/session.h
callbacks.o: ../gui2/settings.h ../gui2/helpfiles.h
callbacks.o: ../gui2/treeutils.h ../gui2/selector.h
console.o: ../gui2/gretl.h winconfig.h ../lib/src/libgretl.h
console.o: ../lib/src/commands.h ../lib/src/errors.h
console.o: ../lib/src/estimate.h ../lib/src/generate.h
console.o: ../lib/src/gretl_utils.h ../lib/src/compare.h
console.o: ../lib/src/pvalues.h ../lib/src/dataio.h
console.o: ../lib/src/strutils.h ../lib/src/describe.h
console.o: ../lib/src/printout.h ../lib/src/modelprint.h
console.o: ../lib/src/texprint.h ../lib/src/var.h
console.o: ../lib/src/interact.h ../lib/src/graphing.h
console.o: ../lib/src/monte_carlo.h ../lib/src/nonparam.h
console.o: ../lib/src/discrete.h ../lib/src/subsample.h
console.o: ../lib/src/calendar.h ../lib/src/system.h
console.o: ../lib/src/gretl_win32.h ../gui2/gretltypes.h
console.o: ../gui2/gui_utils.h ../gui2/callbacks.h ../gui2/dialogs.h
console.o: ../gui2/library.h ../gui2/session.h ../gui2/settings.h
console.o: ../gui2/helpfiles.h
database.o: ../gui2/gretl.h winconfig.h ../lib/src/libgretl.h
database.o: ../lib/src/commands.h ../lib/src/errors.h
database.o: ../lib/src/estimate.h ../lib/src/generate.h
database.o: ../lib/src/gretl_utils.h ../lib/src/compare.h
database.o: ../lib/src/pvalues.h ../lib/src/dataio.h
database.o: ../lib/src/strutils.h ../lib/src/describe.h
database.o: ../lib/src/printout.h ../lib/src/modelprint.h
database.o: ../lib/src/texprint.h ../lib/src/var.h
database.o: ../lib/src/interact.h ../lib/src/graphing.h
database.o: ../lib/src/monte_carlo.h ../lib/src/nonparam.h
database.o: ../lib/src/discrete.h ../lib/src/subsample.h
database.o: ../lib/src/calendar.h ../lib/src/system.h
database.o: ../lib/src/gretl_win32.h ../gui2/gretltypes.h
database.o: ../gui2/gui_utils.h ../gui2/callbacks.h ../gui2/dialogs.h
database.o: ../gui2/library.h ../gui2/session.h ../gui2/settings.h
database.o: ../gui2/helpfiles.h ../gui2/treeutils.h
database.o: ../gui2/boxplots.h
datafiles.o: ../gui2/gretl.h winconfig.h ../lib/src/libgretl.h
datafiles.o: ../lib/src/commands.h ../lib/src/errors.h
datafiles.o: ../lib/src/estimate.h ../lib/src/generate.h
datafiles.o: ../lib/src/gretl_utils.h ../lib/src/compare.h
datafiles.o: ../lib/src/pvalues.h ../lib/src/dataio.h
datafiles.o: ../lib/src/strutils.h ../lib/src/describe.h
datafiles.o: ../lib/src/printout.h ../lib/src/modelprint.h
datafiles.o: ../lib/src/texprint.h ../lib/src/var.h
datafiles.o: ../lib/src/interact.h ../lib/src/graphing.h
datafiles.o: ../lib/src/monte_carlo.h ../lib/src/nonparam.h
datafiles.o: ../lib/src/discrete.h ../lib/src/subsample.h
datafiles.o: ../lib/src/calendar.h ../lib/src/system.h
datafiles.o: ../lib/src/gretl_win32.h ../gui2/gretltypes.h
datafiles.o: ../gui2/gui_utils.h ../gui2/callbacks.h
datafiles.o: ../gui2/dialogs.h ../gui2/library.h ../gui2/session.h
datafiles.o: ../gui2/settings.h ../gui2/helpfiles.h
datafiles.o: ../gui2/treeutils.h
dialogs.o: ../gui2/gretl.h winconfig.h ../lib/src/libgretl.h
dialogs.o: ../lib/src/commands.h ../lib/src/errors.h
dialogs.o: ../lib/src/estimate.h ../lib/src/generate.h
dialogs.o: ../lib/src/gretl_utils.h ../lib/src/compare.h
dialogs.o: ../lib/src/pvalues.h ../lib/src/dataio.h
dialogs.o: ../lib/src/strutils.h ../lib/src/describe.h
dialogs.o: ../lib/src/printout.h ../lib/src/modelprint.h
dialogs.o: ../lib/src/texprint.h ../lib/src/var.h
dialogs.o: ../lib/src/interact.h ../lib/src/graphing.h
dialogs.o: ../lib/src/monte_carlo.h ../lib/src/nonparam.h
dialogs.o: ../lib/src/discrete.h ../lib/src/subsample.h
dialogs.o: ../lib/src/calendar.h ../lib/src/system.h
dialogs.o: ../lib/src/gretl_win32.h ../gui2/gretltypes.h
dialogs.o: ../gui2/gui_utils.h ../gui2/callbacks.h ../gui2/dialogs.h
dialogs.o: ../gui2/library.h ../gui2/session.h ../gui2/settings.h
dialogs.o: ../gui2/helpfiles.h ../gui2/ssheet.h ../lib/src/version.h
dialogs.o: build.h ../gui2/selector.h
fileselect.o: ../gui2/gretl.h winconfig.h ../lib/src/libgretl.h
fileselect.o: ../lib/src/commands.h ../lib/src/errors.h
fileselect.o: ../lib/src/estimate.h ../lib/src/generate.h
fileselect.o: ../lib/src/gretl_utils.h ../lib/src/compare.h
fileselect.o: ../lib/src/pvalues.h ../lib/src/dataio.h
fileselect.o: ../lib/src/strutils.h ../lib/src/describe.h
fileselect.o: ../lib/src/printout.h ../lib/src/modelprint.h
fileselect.o: ../lib/src/texprint.h ../lib/src/var.h
fileselect.o: ../lib/src/interact.h ../lib/src/graphing.h
fileselect.o: ../lib/src/monte_carlo.h ../lib/src/nonparam.h
fileselect.o: ../lib/src/discrete.h ../lib/src/subsample.h
fileselect.o: ../lib/src/calendar.h ../lib/src/system.h
fileselect.o: ../lib/src/gretl_win32.h ../gui2/gretltypes.h
fileselect.o: ../gui2/gui_utils.h ../gui2/callbacks.h
fileselect.o: ../gui2/dialogs.h ../gui2/library.h ../gui2/session.h
fileselect.o: ../gui2/settings.h ../gui2/helpfiles.h
fileselect.o: ../gui2/fileselect.h ../gui2/boxplots.h
fileselect.o: ../gui2/gpt_control.h
gpt_control.o: ../gui2/gretl.h winconfig.h ../lib/src/libgretl.h
gpt_control.o: ../lib/src/commands.h ../lib/src/errors.h
gpt_control.o: ../lib/src/estimate.h ../lib/src/generate.h
gpt_control.o: ../lib/src/gretl_utils.h ../lib/src/compare.h
gpt_control.o: ../lib/src/pvalues.h ../lib/src/dataio.h
gpt_control.o: ../lib/src/strutils.h ../lib/src/describe.h
gpt_control.o: ../lib/src/printout.h ../lib/src/modelprint.h
gpt_control.o: ../lib/src/texprint.h ../lib/src/var.h
gpt_control.o: ../lib/src/interact.h ../lib/src/graphing.h
gpt_control.o: ../lib/src/monte_carlo.h ../lib/src/nonparam.h
gpt_control.o: ../lib/src/discrete.h ../lib/src/subsample.h
gpt_control.o: ../lib/src/calendar.h ../lib/src/system.h
gpt_control.o: ../lib/src/gretl_win32.h ../gui2/gretltypes.h
gpt_control.o: ../gui2/gui_utils.h ../gui2/callbacks.h
gpt_control.o: ../gui2/dialogs.h ../gui2/library.h ../gui2/session.h
gpt_control.o: ../gui2/settings.h ../gui2/helpfiles.h
gpt_control.o: ../gui2/gpt_control.h
gretl.o: ../gui2/gretl.h winconfig.h ../lib/src/libgretl.h
gretl.o: ../lib/src/commands.h ../lib/src/errors.h
gretl.o: ../lib/src/estimate.h ../lib/src/generate.h
gretl.o: ../lib/src/gretl_utils.h ../lib/src/compare.h
gretl.o: ../lib/src/pvalues.h ../lib/src/dataio.h
gretl.o: ../lib/src/strutils.h ../lib/src/describe.h
gretl.o: ../lib/src/printout.h ../lib/src/modelprint.h
gretl.o: ../lib/src/texprint.h ../lib/src/var.h ../lib/src/interact.h
gretl.o: ../lib/src/graphing.h ../lib/src/monte_carlo.h
gretl.o: ../lib/src/nonparam.h ../lib/src/discrete.h
gretl.o: ../lib/src/subsample.h ../lib/src/calendar.h
gretl.o: ../lib/src/system.h ../lib/src/gretl_win32.h
gretl.o: ../gui2/gretltypes.h ../gui2/gui_utils.h ../gui2/callbacks.h
gretl.o: ../gui2/dialogs.h ../gui2/library.h ../gui2/session.h
gretl.o: ../gui2/settings.h ../gui2/helpfiles.h ../gui2/treeutils.h
gretl.o: ../gui2/ssheet.h ../gui2/gpt_control.h ../gui2/guiprint.h
gretl.o: ../pixmaps/mini.calc.xpm ../pixmaps/mini.edit.xpm
gretl.o: ../pixmaps/mini.sh.xpm ../pixmaps/mini.session.xpm
gretl.o: ../pixmaps/mini.manual.xpm ../pixmaps/mini.netscape.xpm
gretl.o: ../pixmaps/mini.pdf.xpm ../pixmaps/mini.plot.xpm
gretl.o: ../pixmaps/mini.ofolder.xpm
gtkfontselhack.o: ../gui2/gretl.h winconfig.h ../lib/src/libgretl.h
gtkfontselhack.o: ../lib/src/commands.h ../lib/src/errors.h
gtkfontselhack.o: ../lib/src/estimate.h ../lib/src/generate.h
gtkfontselhack.o: ../lib/src/gretl_utils.h ../lib/src/compare.h
gtkfontselhack.o: ../lib/src/pvalues.h ../lib/src/dataio.h
gtkfontselhack.o: ../lib/src/strutils.h ../lib/src/describe.h
gtkfontselhack.o: ../lib/src/printout.h ../lib/src/modelprint.h
gtkfontselhack.o: ../lib/src/texprint.h ../lib/src/var.h
gtkfontselhack.o: ../lib/src/interact.h ../lib/src/graphing.h
gtkfontselhack.o: ../lib/src/monte_carlo.h ../lib/src/nonparam.h
gtkfontselhack.o: ../lib/src/discrete.h ../lib/src/subsample.h
gtkfontselhack.o: ../lib/src/calendar.h ../lib/src/system.h
gtkfontselhack.o: ../lib/src/gretl_win32.h ../gui2/gretltypes.h
gtkfontselhack.o: ../gui2/gui_utils.h ../gui2/callbacks.h
gtkfontselhack.o: ../gui2/dialogs.h ../gui2/library.h
gtkfontselhack.o: ../gui2/session.h ../gui2/settings.h
gtkfontselhack.o: ../gui2/helpfiles.h ../gui2/gtkfontselhack.h
gtkfontselhack.o: ../gui2/font_filter.c
guiprint.o: ../gui2/gretl.h winconfig.h ../lib/src/libgretl.h
guiprint.o: ../lib/src/commands.h ../lib/src/errors.h
guiprint.o: ../lib/src/estimate.h ../lib/src/generate.h
guiprint.o: ../lib/src/gretl_utils.h ../lib/src/compare.h
guiprint.o: ../lib/src/pvalues.h ../lib/src/dataio.h
guiprint.o: ../lib/src/strutils.h ../lib/src/describe.h
guiprint.o: ../lib/src/printout.h ../lib/src/modelprint.h
guiprint.o: ../lib/src/texprint.h ../lib/src/var.h
guiprint.o: ../lib/src/interact.h ../lib/src/graphing.h
guiprint.o: ../lib/src/monte_carlo.h ../lib/src/nonparam.h
guiprint.o: ../lib/src/discrete.h ../lib/src/subsample.h
guiprint.o: ../lib/src/calendar.h ../lib/src/system.h
guiprint.o: ../lib/src/gretl_win32.h ../gui2/gretltypes.h
guiprint.o: ../gui2/gui_utils.h ../gui2/callbacks.h ../gui2/dialogs.h
guiprint.o: ../gui2/library.h ../gui2/session.h ../gui2/settings.h
guiprint.o: ../gui2/helpfiles.h ../gui2/selector.h
gui_utils.o: ../gui2/gretl.h winconfig.h ../lib/src/libgretl.h
gui_utils.o: ../lib/src/commands.h ../lib/src/errors.h
gui_utils.o: ../lib/src/estimate.h ../lib/src/generate.h
gui_utils.o: ../lib/src/gretl_utils.h ../lib/src/compare.h
gui_utils.o: ../lib/src/pvalues.h ../lib/src/dataio.h
gui_utils.o: ../lib/src/strutils.h ../lib/src/describe.h
gui_utils.o: ../lib/src/printout.h ../lib/src/modelprint.h
gui_utils.o: ../lib/src/texprint.h ../lib/src/var.h
gui_utils.o: ../lib/src/interact.h ../lib/src/graphing.h
gui_utils.o: ../lib/src/monte_carlo.h ../lib/src/nonparam.h
gui_utils.o: ../lib/src/discrete.h ../lib/src/subsample.h
gui_utils.o: ../lib/src/calendar.h ../lib/src/system.h
gui_utils.o: ../lib/src/gretl_win32.h ../gui2/gretltypes.h
gui_utils.o: ../gui2/gui_utils.h ../gui2/callbacks.h
gui_utils.o: ../gui2/dialogs.h ../gui2/library.h ../gui2/session.h
gui_utils.o: ../gui2/settings.h ../gui2/helpfiles.h
gui_utils.o: ../gui2/guiprint.h ../gui2/series_view.h
helpfiles.o: ../gui2/gretl.h winconfig.h ../lib/src/libgretl.h
helpfiles.o: ../lib/src/commands.h ../lib/src/errors.h
helpfiles.o: ../lib/src/estimate.h ../lib/src/generate.h
helpfiles.o: ../lib/src/gretl_utils.h ../lib/src/compare.h
helpfiles.o: ../lib/src/pvalues.h ../lib/src/dataio.h
helpfiles.o: ../lib/src/strutils.h ../lib/src/describe.h
helpfiles.o: ../lib/src/printout.h ../lib/src/modelprint.h
helpfiles.o: ../lib/src/texprint.h ../lib/src/var.h
helpfiles.o: ../lib/src/interact.h ../lib/src/graphing.h
helpfiles.o: ../lib/src/monte_carlo.h ../lib/src/nonparam.h
helpfiles.o: ../lib/src/discrete.h ../lib/src/subsample.h
helpfiles.o: ../lib/src/calendar.h ../lib/src/system.h
helpfiles.o: ../lib/src/gretl_win32.h ../gui2/gretltypes.h
helpfiles.o: ../gui2/gui_utils.h ../gui2/callbacks.h
helpfiles.o: ../gui2/dialogs.h ../gui2/library.h ../gui2/session.h
helpfiles.o: ../gui2/settings.h ../gui2/helpfiles.h
helpfiles.o: ../gui2/treeutils.h
library.o: ../gui2/gretl.h winconfig.h ../lib/src/libgretl.h
library.o: ../lib/src/commands.h ../lib/src/errors.h
library.o: ../lib/src/estimate.h ../lib/src/generate.h
library.o: ../lib/src/gretl_utils.h ../lib/src/compare.h
library.o: ../lib/src/pvalues.h ../lib/src/dataio.h
library.o: ../lib/src/strutils.h ../lib/src/describe.h
library.o: ../lib/src/printout.h ../lib/src/modelprint.h
library.o: ../lib/src/texprint.h ../lib/src/var.h
library.o: ../lib/src/interact.h ../lib/src/graphing.h
library.o: ../lib/src/monte_carlo.h ../lib/src/nonparam.h
library.o: ../lib/src/discrete.h ../lib/src/subsample.h
library.o: ../lib/src/calendar.h ../lib/src/system.h
library.o: ../lib/src/gretl_win32.h ../gui2/gretltypes.h
library.o: ../gui2/gui_utils.h ../gui2/callbacks.h ../gui2/dialogs.h
library.o: ../gui2/library.h ../gui2/session.h ../gui2/settings.h
library.o: ../gui2/helpfiles.h ../gui2/gpt_control.h
library.o: ../lib/src/cmdlist.h ../gui2/selector.h ../gui2/boxplots.h
library.o: ../gui2/series_view.h ../cli/common.c
selector.o: ../gui2/gretl.h winconfig.h ../lib/src/libgretl.h
selector.o: ../lib/src/commands.h ../lib/src/errors.h
selector.o: ../lib/src/estimate.h ../lib/src/generate.h
selector.o: ../lib/src/gretl_utils.h ../lib/src/compare.h
selector.o: ../lib/src/pvalues.h ../lib/src/dataio.h
selector.o: ../lib/src/strutils.h ../lib/src/describe.h
selector.o: ../lib/src/printout.h ../lib/src/modelprint.h
selector.o: ../lib/src/texprint.h ../lib/src/var.h
selector.o: ../lib/src/interact.h ../lib/src/graphing.h
selector.o: ../lib/src/monte_carlo.h ../lib/src/nonparam.h
selector.o: ../lib/src/discrete.h ../lib/src/subsample.h
selector.o: ../lib/src/calendar.h ../lib/src/system.h
selector.o: ../lib/src/gretl_win32.h ../gui2/gretltypes.h
selector.o: ../gui2/gui_utils.h ../gui2/callbacks.h ../gui2/dialogs.h
selector.o: ../gui2/library.h ../gui2/session.h ../gui2/settings.h
selector.o: ../gui2/helpfiles.h ../gui2/selector.h
selector.o: ../gui2/treeutils.h
series_view.o: ../gui2/gretl.h winconfig.h ../lib/src/libgretl.h
series_view.o: ../lib/src/commands.h ../lib/src/errors.h
series_view.o: ../lib/src/estimate.h ../lib/src/generate.h
series_view.o: ../lib/src/gretl_utils.h ../lib/src/compare.h
series_view.o: ../lib/src/pvalues.h ../lib/src/dataio.h
series_view.o: ../lib/src/strutils.h ../lib/src/describe.h
series_view.o: ../lib/src/printout.h ../lib/src/modelprint.h
series_view.o: ../lib/src/texprint.h ../lib/src/var.h
series_view.o: ../lib/src/interact.h ../lib/src/graphing.h
series_view.o: ../lib/src/monte_carlo.h ../lib/src/nonparam.h
series_view.o: ../lib/src/discrete.h ../lib/src/subsample.h
series_view.o: ../lib/src/calendar.h ../lib/src/system.h
series_view.o: ../lib/src/gretl_win32.h ../gui2/gretltypes.h
series_view.o: ../gui2/gui_utils.h ../gui2/callbacks.h
series_view.o: ../gui2/dialogs.h ../gui2/library.h ../gui2/session.h
series_view.o: ../gui2/settings.h ../gui2/helpfiles.h
session.o: ../gui2/gretl.h winconfig.h ../lib/src/libgretl.h
session.o: ../lib/src/commands.h ../lib/src/errors.h
session.o: ../lib/src/estimate.h ../lib/src/generate.h
session.o: ../lib/src/gretl_utils.h ../lib/src/compare.h
session.o: ../lib/src/pvalues.h ../lib/src/dataio.h
session.o: ../lib/src/strutils.h ../lib/src/describe.h
session.o: ../lib/src/printout.h ../lib/src/modelprint.h
session.o: ../lib/src/texprint.h ../lib/src/var.h
session.o: ../lib/src/interact.h ../lib/src/graphing.h
session.o: ../lib/src/monte_carlo.h ../lib/src/nonparam.h
session.o: ../lib/src/discrete.h ../lib/src/subsample.h
session.o: ../lib/src/calendar.h ../lib/src/system.h
session.o: ../lib/src/gretl_win32.h ../gui2/gretltypes.h
session.o: ../gui2/gui_utils.h ../gui2/callbacks.h ../gui2/dialogs.h
session.o: ../gui2/library.h ../gui2/session.h ../gui2/settings.h
session.o: ../gui2/helpfiles.h ../gui2/selector.h ../gui2/boxplots.h
session.o: ../gui2/ssheet.h ../gui2/gpt_control.h ../gui2/guiprint.h
session.o: ../pixmaps/model.xpm ../pixmaps/boxplot.xpm
session.o: ../pixmaps/gnuplot.xpm ../pixmaps/xfm_sc.xpm
session.o: ../pixmaps/xfm_info.xpm ../pixmaps/xfm_text.xpm
session.o: ../pixmaps/xfm_make.xpm ../pixmaps/rhohat.xpm
session.o: ../pixmaps/summary.xpm
settings.o: ../gui2/gretl.h winconfig.h ../lib/src/libgretl.h
settings.o: ../lib/src/commands.h ../lib/src/errors.h
settings.o: ../lib/src/estimate.h ../lib/src/generate.h
settings.o: ../lib/src/gretl_utils.h ../lib/src/compare.h
settings.o: ../lib/src/pvalues.h ../lib/src/dataio.h
settings.o: ../lib/src/strutils.h ../lib/src/describe.h
settings.o: ../lib/src/printout.h ../lib/src/modelprint.h
settings.o: ../lib/src/texprint.h ../lib/src/var.h
settings.o: ../lib/src/interact.h ../lib/src/graphing.h
settings.o: ../lib/src/monte_carlo.h ../lib/src/nonparam.h
settings.o: ../lib/src/discrete.h ../lib/src/subsample.h
settings.o: ../lib/src/calendar.h ../lib/src/system.h
settings.o: ../lib/src/gretl_win32.h ../gui2/gretltypes.h
settings.o: ../gui2/gui_utils.h ../gui2/callbacks.h ../gui2/dialogs.h
settings.o: ../gui2/library.h ../gui2/session.h ../gui2/settings.h
settings.o: ../gui2/helpfiles.h ../gui2/gtkfontselhack.h
ssheet.o: ../gui2/gretl.h winconfig.h ../lib/src/libgretl.h
ssheet.o: ../lib/src/commands.h ../lib/src/errors.h
ssheet.o: ../lib/src/estimate.h ../lib/src/generate.h
ssheet.o: ../lib/src/gretl_utils.h ../lib/src/compare.h
ssheet.o: ../lib/src/pvalues.h ../lib/src/dataio.h
ssheet.o: ../lib/src/strutils.h ../lib/src/describe.h
ssheet.o: ../lib/src/printout.h ../lib/src/modelprint.h
ssheet.o: ../lib/src/texprint.h ../lib/src/var.h
ssheet.o: ../lib/src/interact.h ../lib/src/graphing.h
ssheet.o: ../lib/src/monte_carlo.h ../lib/src/nonparam.h
ssheet.o: ../lib/src/discrete.h ../lib/src/subsample.h
ssheet.o: ../lib/src/calendar.h ../lib/src/system.h
ssheet.o: ../lib/src/gretl_win32.h ../gui2/gretltypes.h
ssheet.o: ../gui2/gui_utils.h ../gui2/callbacks.h ../gui2/dialogs.h
ssheet.o: ../gui2/library.h ../gui2/session.h ../gui2/settings.h
ssheet.o: ../gui2/helpfiles.h ../gui2/treeutils.h ../gui2/ssheet.h
treeutils.o: ../gui2/gretl.h winconfig.h ../lib/src/libgretl.h
treeutils.o: ../lib/src/commands.h ../lib/src/errors.h
treeutils.o: ../lib/src/estimate.h ../lib/src/generate.h
treeutils.o: ../lib/src/gretl_utils.h ../lib/src/compare.h
treeutils.o: ../lib/src/pvalues.h ../lib/src/dataio.h
treeutils.o: ../lib/src/strutils.h ../lib/src/describe.h
treeutils.o: ../lib/src/printout.h ../lib/src/modelprint.h
treeutils.o: ../lib/src/texprint.h ../lib/src/var.h
treeutils.o: ../lib/src/interact.h ../lib/src/graphing.h
treeutils.o: ../lib/src/monte_carlo.h ../lib/src/nonparam.h
treeutils.o: ../lib/src/discrete.h ../lib/src/subsample.h
treeutils.o: ../lib/src/calendar.h ../lib/src/system.h
treeutils.o: ../lib/src/gretl_win32.h ../gui2/gretltypes.h
treeutils.o: ../gui2/gui_utils.h ../gui2/callbacks.h
treeutils.o: ../gui2/dialogs.h ../gui2/library.h ../gui2/session.h
treeutils.o: ../gui2/settings.h ../gui2/helpfiles.h
treeutils.o: ../gui2/treeutils.h
webget.o: ../gui2/gretl.h winconfig.h ../lib/src/libgretl.h
webget.o: ../lib/src/commands.h ../lib/src/errors.h
webget.o: ../lib/src/estimate.h ../lib/src/generate.h
webget.o: ../lib/src/gretl_utils.h ../lib/src/compare.h
webget.o: ../lib/src/pvalues.h ../lib/src/dataio.h
webget.o: ../lib/src/strutils.h ../lib/src/describe.h
webget.o: ../lib/src/printout.h ../lib/src/modelprint.h
webget.o: ../lib/src/texprint.h ../lib/src/var.h
webget.o: ../lib/src/interact.h ../lib/src/graphing.h
webget.o: ../lib/src/monte_carlo.h ../lib/src/nonparam.h
webget.o: ../lib/src/discrete.h ../lib/src/subsample.h
webget.o: ../lib/src/calendar.h ../lib/src/system.h
webget.o: ../lib/src/gretl_win32.h ../gui2/gretltypes.h
webget.o: ../gui2/gui_utils.h ../gui2/callbacks.h ../gui2/dialogs.h
webget.o: ../gui2/library.h ../gui2/session.h ../gui2/settings.h
webget.o: ../gui2/helpfiles.h ../gui2/webget.h
