# Makefile for cross-compiling gretl for win32 on Linux,
# using mingw compiler

include config.mk

# directories
topsrc = ..
libdir  = $(topsrc)/lib/src
extradir = ../gui2/gtkextra-lite
imports = ./import-libs

# tools
CC = $(MGW_PREFIX)-gcc -fnative-struct -Wall -O2
AS = $(MGW_PREFIX)-as
DLLWRAP = $(MGW_PREFIX)-dllwrap
DLLTOOL = $(MGW_PREFIX)-dlltool

# libraries
GTKLIBS := -L$(imports) $(shell PKG_CONFIG_PATH=$(PKG_CONFIG_PATH) \
                          pkg-config --libs gtk+-win32-2.0)
ZLIB_LIB = -lz.dll
CLI_LIBS = -L$(imports) -lreadline.dll
GUI_LIBS = $(GTKLIBS) -lgtkextra-lite -lwsock32 #-lhhctrl
GRETLLIB = -L$(imports) -lgretl -lintl.dll
PIXBUFLIB = -L$(imports) -lgdk-pixbuf -lpng -lz.dll
LIBS = -lm $(GRETLLIB) -ladvapi32

# flags/defines
LDFLAGS = -s 
CFLAGS := -DHAVE_READLINE -DOS_WIN32 -I. -I$(libdir) -I$(extradir) \
$(shell PKG_CONFIG_PATH=$(PKG_CONFIG_PATH) pkg-config --cflags gtk+-win32-2.0)

# source paths
vpath %.c $(libdir)
vpath %.c $(topsrc)/cephes
vpath %.c $(topsrc)/cli
vpath %.c $(topsrc)/gui2
vpath %.c $(extradir)
vpath %.c $(topsrc)/plugin
vpath %.c $(topsrc)/plugin/libole2
vpath %.h $(libdir):$(topsrc)/cli:$(topsrc)/gui2:$(extradir):$(topsrc)/plugin:.
vpath %.h $(topsrc)/cephes

#### End of system configuration section. ####

CLIPROG = gretlcli.exe
GUIPROG = gretlw32.exe
UPDATER = gretl_updater.exe
PLUGINS = plugins/stats_tables.dll plugins/panel_data.dll plugins/longname.dll \
	plugins/progress_bar-2.dll plugins/gnumeric_import-2.dll \
	plugins/excel_import-2.dll plugins/des_import.dll plugins/mp_ols.dll \
	plugins/range-mean.dll

WRES    = gretlw32_res.o
RC      = $(MGW_PREFIX)-windres
RCFLAGS = --include-dir $(MGW_INC) \
	--define __WIN32__ --define __WIN95__ --define MSRC \
	--define __GNUWIN32__
RCOUT   = gretlw32_res.o
RES2COFF = echo gretlw32_res.o

LIBSRC = estimate.c compare.c gretl_utils.c errors.c pvalues.c dataio.c \
       generate.c strutils.c describe.c printout.c modelprint.c texprint.c var.c \
       interact.c graphing.c win32.c monte_carlo.c nonparam.c discrete.c \
       subsample.c calendar.c
LIBHDR = libgretl.h commands.h estimate.h gretl_utils.h errors.h pvalues.h \
       compare.h dataio.h generate.h strutils.h describe.h printout.h modelprint.h \
       version.h texprint.h var.h interact.h graphing.h monte_carlo.h \
       nonparam.h discrete.h subsample.h calendar.h win32.h
LIBOBJ = $(LIBSRC:.c=.o)

CLISRC = gretlcli.c complete.c 
CLIOBJ = $(CLISRC:.c=.o)

GUISRC = gretl.c library.c gui_utils.c callbacks.c dialogs.c ssheet.c \
	datafiles.c database.c gpt_control.c session.c webget.c \
	calculator.c fileselect.c guiprint.c boxplots.c \
	selector.c settings.c helpfiles.c treeutils.c console.c 
GUIHDR = gretl.h library.h gui_utils.h callbacks.h dialogs.h session.h \
	webget.h guiprint.h selector.h settings.h helpfiles.h \
	treeutils.h
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

MO = mo/es.mo

# targets begin
all: dirs prob-dll dll extra-dll plugin $(CLIPROG) $(GUIPROG) $(MO)

$(CLIPROG): $(CLIOBJ) 
	$(CC) -s -o $@ $^ $(LIBS) $(CLI_LIBS) $(ZLIB_LIB)

$(GUIPROG): $(GUIOBJ) $(HACK_OBJ) $(WRES) 
	$(CC) -mwindows -s -o $@ $^ $(LIBS) $(GUI_LIBS) $(ZLIB_LIB)

$(CLIOBJ): $(CLI_SRC) 

$(GUIOBJ): $(GUISRC) $(GUIHDR) build.h

$(HACK_OBJ): $(HACK_SRC) $(HACK_HDR)

build.h:
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
	rm -f *.o *.def $(CLIPROG) $(GUIPROG) dlls/libgretl.dll
	rm -f plugins/* $(PLUGINS)
	rm -f dlls/libgtkextra-lite.dll $(imports)/libgtkextra-lite.a $(imports)/libgtkextra-lite.def
	rm -f gretl.stamp build.h
	rm -f windist/MANIFEST windist/gretl2.iss 
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
	-o $@ $^ -lm -L$(imports) -lxml -lz.dll -lintl.dll -lprob

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

plugins/panel_data.dll: panel_data.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) $(GTKLIBS)

plugins/progress_bar-2.dll: progress_bar-2.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GTKLIBS) -L$(imports) -lintl.dll

plugins/gnumeric_import-2.dll: gnumeric_import.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ -L$(imports) -lxml \
	$(GRETLLIB) $(GTKLIBS)

plugins/excel_import-2.dll: excel_import.o workbook.o ms-ole.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) $(GTKLIBS) -lmoldname

plugins/des_import.dll: des_import.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) -lz.dll

plugins/longname.dll: longname.o 
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ -lkernel32

plugins/mp_ols.dll: mp_ols.o dllinit.o 
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) -lgmp

panel_data.o: panel_data.c 
	$(CC) -c $(CFLAGS) $(GLIBINC) $<

OLEFLAGS = -I$(topsrc)/plugin

ms-ole.o: ms-ole.c $(topsrc)/plugin/libole2/ms-ole.h
	$(CC) -c $(CFLAGS) $(OLEFLAGS) $(GLIBINC) $<

workbook.o: workbook.c $(topsrc)/plugin/libole2/ms-ole.h
	$(CC) -c $(CFLAGS) $(OLEFLAGS) $(GLIBINC) $<


