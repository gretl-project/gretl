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
CFLAGS = -DOS_WIN32 -I. -I$(libdir) -I$(extradir) -I$(plugindir)
GTK_CFLAGS := $(shell PKG_CONFIG_PATH=$(PKG_CONFIG_PATH) pkg-config --cflags gtk+-win32-2.0)

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
vpath %.h $(topsrc)/cephes:..
vpath %.rc $(topsrc)/win32

#### End of system configuration section. ####

CLIPROG = gretlcli.exe
GUIPROG = gretlw32.exe
UPDATER = gretl_updater.exe
CHECKER = nisttest.exe

PLUGINSRC = \
	excel_import.c \
	gnumeric_import.c \
	gretl_matrix.c \
	johansen.c \
	lad.c \
	leverage.c \
	longname.c \
	mp_ols.c \
	panel_data.c \
	progress_bar.c \
	range-mean.c \
	stats_tables.c \
	sur.c \
	tramo_options.c \
	tramo-x12a.c \
	workbook.c

PLUGINOBJ = $(PLUGINSRC:.c=.o)

PLUGINS = plugins/stats_tables.dll plugins/panel_data.dll plugins/longname.dll \
	plugins/progress_bar.dll plugins/gnumeric_import.dll \
	plugins/excel_import.dll plugins/mp_ols.dll plugins/sur.dll plugins/johansen.dll \
	plugins/range-mean.dll plugins/lad.dll plugins/tramo-x12a.dll \
	plugins/leverage.dll

WRES    = gretlw32_res.o
RC      = $(MGW_PREFIX)-windres
RCFLAGS = --define __WIN32__ --define __WIN95__ --define MSRC \
	--define __GNUWIN32__
RCOUT   = gretlw32_res.o
RES2COFF = echo gretlw32_res.o

LIBSRC = estimate.c compare.c gretl_utils.c gretl_errors.c pvalues.c dataio.c \
       generate.c strutils.c describe.c printout.c modelprint.c texprint.c var.c \
       interact.c graphing.c gretl_win32.c monte_carlo.c nonparam.c discrete.c \
       subsample.c calendar.c system.c random.c
LIBOBJ = $(LIBSRC:.c=.o)

CLISRC = gretlcli.c complete.c 
CLIOBJ = $(CLISRC:.c=.o)

CHECKSRC = nistcheck.c

GUISRC = gretl.c library.c gui_utils.c callbacks.c dialogs.c ssheet.c \
	datafiles.c database.c gpt_control.c session.c webget.c \
	calculator.c fileselect.c guiprint.c boxplots.c series_view.c \
	selector.c settings.c helpfiles.c treeutils.c console.c model_table.c
GUIOBJ = $(GUISRC:.c=.o)

HACK_SRC = gtkfontselhack.c font_filter.c
HACK_OBJ = gtkfontselhack.o

EXTRASRC = gtkpsfont.c gtkplotpc.c gtkplotps.c
EXTRAOBJ = $(EXTRASRC:.c=.o)

PROBSRC = bdtr.c btdtr.c chdtr.c drand.c expx2.c fdtr.c gamma.c gdtr.c \
	igam.c igami.c incbet.c incbi.c mtherr.c nbdtr.c ndtr.c ndtri.c pdtr.c \
	stdtr.c unity.c polevl.c const.c
PROBOBJ = $(PROBSRC:.c=.o)

OBJS = $(LIBOBJ) $(CLIOBJ) $(GUIOBJ) $(PLUGINOBJ)

MO = mo/es.mo mo/fr.mo

%.o: %.c
	$(CC) -c $(CFLAGS) $(GTK_CFLAGS) $<
	$(CC) -MM -MG $(CFLAGS) $< > $*.d 2>/dev/null

# targets begin

all: buildstamp dirs prob-dll dll extra-dll plugin $(CLIPROG) $(GUIPROG) \
	$(CHECKER) $(MO)

-include $(OBJS:.o=.d)

$(CLIPROG): $(CLIOBJ) 
	$(CC) -s -o $@ $^ $(LIBS) $(CLI_LIBS) $(ZLIB_LIB) -lmingwex

$(GUIPROG): $(GUIOBJ) $(HACK_OBJ) $(WRES) 
	$(CC) -mwindows -s -o $@ $^ $(LIBS) $(GUI_LIBS) -lpng $(ZLIB_LIB) -lmingwex

$(CHECKER): $(CHECKSRC)
	$(CC) -s -I. -I../lib/src -DOS_WIN32 -DUSE_GMP -o $@ $^ $(LIBS) $(ZLIB_LIB)

buildstamp:
	./builddate.pl

$(WRES): gretlw32.rc 
	$(RC) $(RCFLAGS) $< $(RCOUT)
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
	rm -f *.o *.d *.def $(CLIPROG) $(GUIPROG) $(CHECKER) dlls/*
	rm -f plugins/* $(PLUGINS)
	rm -f dlls/libgtkextra-lite.dll $(imports)/libgtkextra-lite.a $(imports)/libgtkextra-lite.def
	rm -f gretl.stamp build.h
	rm -f windist/MANIFEST windist/gretl.iss 
	rm -rf windist/gretl windist/gnuplot windist/Output
	make -C updater clean
	make -C mo clean

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

plugins/leverage.dll: leverage.o gretl_matrix.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) -llapack -lblas -lf2c

ms-ole.o: ms-ole.c $(topsrc)/plugin/libole2/ms-ole.h
	$(CC) -c $(CFLAGS) $(GTK_CFLAGS) $<

workbook.o: workbook.c $(topsrc)/plugin/libole2/ms-ole.h
	$(CC) -c $(CFLAGS) $(GTK_CFLAGS) $<

