# Makefile for cross-compiling gretl for win32 on Linux,
# using mingw compiler

include config.mk

# directories
topsrc = ..
libdir  = $(topsrc)/lib/src
plugindir = $(topsrc)/plugin
guidir = $(topsrc)/gui2
extradir = $(guidir)/gtkextra-lite
imports = ./import-libs

# tools
CC = $(MGW_PREFIX)-gcc -Wall -O2 -mms-bitfields -DWIN32
AS = $(MGW_PREFIX)-as
DLLWRAP = $(MGW_PREFIX)-dllwrap
DLLWRAP_FLAGS = --as=$(AS) --export-all --driver-name $(CC) -s
RC = $(MGW_PREFIX)-windres
RCFLAGS = --define __WIN32__ --define __WIN95__ --define MSRC \
	--define __GNUWIN32__

# libraries
GTKLIBS := -L$(imports) $(shell PKG_CONFIG_PATH=$(PKG_CONFIG_PATH) \
                          pkg-config --libs gtk+-win32-2.0)
GLIBLIB := -L$(imports) $(shell PKG_CONFIG_PATH=$(PKG_CONFIG_PATH) \
                          pkg-config --libs glib-2.0)
CLI_LIBS = -L$(imports) -lreadline.dll
GUI_LIBS = $(GTKLIBS) -lgtkextra-lite -lgtksourceview -lwsock32 #-lhhctrl
GRETLLIB = -L$(imports) -lgretl -lintl
LIBS = -lm $(GRETLLIB) -ladvapi32
LAPACK_LIBS = -llapack -lblas -lf2c

# flags/defines
CFLAGS = -I. -I$(libdir) -I$(extradir) -I$(plugindir) -I$(guidir) -DUSE_GMP
GTK_CFLAGS := $(shell PKG_CONFIG_PATH=$(PKG_CONFIG_PATH) pkg-config --cflags gtk+-win32-2.0)

# source paths
vpath %.c $(libdir)
vpath %.c $(topsrc)/cephes
vpath %.c $(topsrc)/minpack
vpath %.c $(topsrc)/cli
vpath %.c $(guidir)
vpath %.c $(extradir)
vpath %.c $(plugindir)
vpath %.c $(plugindir)/libole2
vpath %.c $(topsrc)/tests
vpath %.h $(libdir):$(topsrc)/cli:$(guidir):$(extradir):$(plugindir):.
vpath %.h $(topsrc)/cephes:..
vpath %.rc $(topsrc)/win32

#### End of system configuration section. ####

PLUGINSRC = \
	excel_import.c \
	gnumeric_import.c \
	johansen.c \
	kernel.c \
	lad.c \
	vif.c \
	leverage.c \
	mp_ols.c \
	panel_data.c \
	pca.c \
	progress_bar.c \
	range-mean.c \
	stats_tables.c \
	sysest.c \
	fiml.c \
	liml.c \
	tramo_options.c \
	tramo-x12a.c \
	nistcheck.c \
	arma.c \
        arma_x12.c \
	workbook.c \
	logistic.c \
	tobit.c \
	fcp.c \
	garch.c \
        audio.c \
	urcdist.c 

PLUGINOBJ = $(PLUGINSRC:.c=.o)

PLUGINS = plugins/stats_tables.dll \
	plugins/panel_data.dll \
	plugins/progress_bar.dll \
	plugins/gnumeric_import.dll \
	plugins/excel_import.dll \
	plugins/mp_ols.dll \
	plugins/sysest.dll \
	plugins/johansen.dll \
	plugins/kernel.dll \
	plugins/arma.dll \
	plugins/arma_x12.dll \
	plugins/range-mean.dll \
	plugins/lad.dll \
	plugins/tramo-x12a.dll \
	plugins/leverage.dll \
	plugins/pca.dll \
	plugins/nistcheck.dll \
	plugins/logistic.dll \
	plugins/tobit.dll \
	plugins/garch.dll \
        plugins/audio.dll \
	plugins/vif.dll \
	plugins/urcdist.dll

LIBSRC = bhhh_max.c \
	calendar.c \
	compare.c \
	compat.c \
	dataio.c \
	dbread.c \
	describe.c \
        discrete.c \
	estimate.c \
	generate.c \
	genstack.c \
	graphing.c \
	gretl_commands.c \
	gretl_errors.c \
	gretl_func.c \
	gretl_list.c \
	gretl_matrix.c \
	gretl_model.c \
	gretl_paths.c \
	gretl_restrict.c \
	gretl_string_table.c \
	gretl_utils.c \
	gretl_win32.c \
	interact.c \
	libset.c \
	missing.c \
	modelprint.c \
	modelspec.c \
	monte_carlo.c \
	nls.c \
	nonparam.c \
	options.c \
	plugins.c \
	printout.c \
	pvalues.c \
	qr_estimate.c \
	random.c \
	strutils.c \
	subsample.c \
	system.c \
	texprint.c \
	transforms.c \
	var.c

LIBOBJ = $(LIBSRC:.c=.o)

CLISRC = gretlcli.c complete.c 

CLIOBJ = $(CLISRC:.c=.o)

GUISRC = about.c \
	boxplots.c \
	calculator.c \
	callbacks.c \
	cmdstack.c \
	console.c \
	database.c \
	datafiles.c \
	dialogs.c \
	dlgutils.c \
	filelists.c \
	fileselect.c \
	gpt_control.c \
	gpt_dialog.c \
	graph_page.c \
	gretl.c \
	gretlwin32.c \
	guiprint.c \
	gui_utils.c \
	helpfiles.c \
	library.c \
	menustate.c \
	model_table.c \
	objectsave.c \
	obsbutton.c \
	selector.c \
	series_view.c \
	session.c \
	settings.c \
	ssheet.c \
	textbuf.c \
	textutil.c \
	toolbar.c \
	treeutils.c \
	webget.c

GUIOBJ = $(GUISRC:.c=.o)

EXTRASRC = gtkpsfont.c gtkplotpc.c gtkplotps.c
EXTRAOBJ = $(EXTRASRC:.c=.o)

PROBSRC = bdtr.c btdtr.c chdtr.c drand.c expx2.c fdtr.c gamma.c gdtr.c \
	igam.c igami.c incbet.c incbi.c mtherr.c nbdtr.c ndtr.c ndtri.c pdtr.c \
	stdtr.c unity.c polevl.c const.c
PROBOBJ = $(PROBSRC:.c=.o)

MINSRC = chkder.c dpmpar.c enorm.c fdjac2.c lmder1.c lmder.c lmdif.c lmpar.c \
	qrfac.c qrsolv.c
MINOBJ = $(MINSRC:.c=.o)

OBJS = $(LIBOBJ) $(CLIOBJ) $(GUIOBJ) $(PLUGINOBJ)

MOFILES = mo/es.mo mo/fr.mo mo/it.mo

%.o: %.c
	$(CC) -c $(CFLAGS) $(GTK_CFLAGS) $<
	$(CC) -MM $(CFLAGS) $(GTK_CFLAGS) $< > $*.d 

# targets begin

all: buildstamp dirs dlls/libprob.dll dlls/libgretl.dll dlls/libgtkextra-lite.dll \
	$(PLUGINS) gretlcli.exe gretlw32.exe $(MOFILES)

gretlcli.exe: $(CLIOBJ)
	$(CC) -o $@ $^ $(LIBS) $(CLI_LIBS) -lz -lmingwex -s

gretlw32.exe: $(GUIOBJ) gretlw32_res.o
	$(CC) -mwindows -o $@ $^ $(LIBS) $(GUI_LIBS) -lpng -lz -lmingwex -s

gretlw32_res.o: gretlw32.rc 
	$(RC) $(RCFLAGS) $< $@

gretl_updater.exe:
	make -C updater

$(MOFILES):
	make -C mo

dlls/libgretl.dll: $(LIBOBJ) $(MINOBJ) dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) \
	--output-def libgretl.def --implib $(imports)/libgretl.a \
	-o $@ $^ -lf2c -lm -L$(imports) -lxml2 -lz -lintl -lprob -lgmp \
	-lmingwex $(GLIBLIB) $(LAPACK_LIBS)

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

plugins/kernel.dll: kernel.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) 

plugins/lad.dll: lad.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) 

plugins/vif.dll: vif.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) 

plugins/urcdist.dll: urcdist.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ -lz -lm -lintl 

plugins/panel_data.dll: panel_data.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(LAPACK_LIBS) $(GRETLLIB) $(GTKLIBS) 

plugins/tramo-x12a.dll: tramo-x12a.o tramo_options.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) $(GTKLIBS)

plugins/progress_bar.dll: progress_bar.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GTKLIBS) -L$(imports) -lintl

plugins/gnumeric_import.dll: gnumeric_import.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ -L$(imports) -lxml2 \
	$(GRETLLIB) $(GTKLIBS)

plugins/excel_import.dll: excel_import.o workbook.o ms-ole.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) $(GTKLIBS) -lmoldname

plugins/mp_ols.dll: mp_ols.o dllinit.o 
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) -lgmp.dll

plugins/sysest.dll: sysest.o fiml.o liml.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) $(LAPACK_LIBS)

plugins/johansen.dll: johansen.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) $(LAPACK_LIBS)

plugins/leverage.dll: leverage.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) $(GTKLIBS) $(LAPACK_LIBS)

plugins/pca.dll: pca.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) $(GTKLIBS) $(LAPACK_LIBS)

plugins/nistcheck.dll: nistcheck.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) $(LAPACK_LIBS)

plugins/arma.dll: arma.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) $(LAPACK_LIBS)

plugins/arma_x12.dll: arma_x12.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) $(GTKLIBS) $(LAPACK_LIBS)

plugins/logistic.dll: logistic.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) $(GTKLIBS)

plugins/tobit.dll: tobit.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB)

plugins/garch.dll: garch.o fcp.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) $(LAPACK_LIBS)

plugins/audio.dll: audio.o midi_utils.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) $(GTKLIBS) -lole32 -lsapi

ms-ole.o: ms-ole.c $(topsrc)/plugin/libole2/ms-ole.h
	$(CC) -c $(CFLAGS) $(GTK_CFLAGS) $<

workbook.o: workbook.c $(topsrc)/plugin/libole2/ms-ole.h
	$(CC) -c $(CFLAGS) $(GTK_CFLAGS) $<

midi_utils.o: midi_utils.c $(topsrc)/plugin/midi_utils.h
	$(CC) -c $(CFLAGS) $(GTK_CFLAGS) $<

-include $(OBJS:.o=.d)

# Distribution executable

dist.stamp:
	sh ./mkwindist

.PHONY:

dirs:
	mkdir -p dlls 
	mkdir -p plugins 

buildstamp:
	./builddate.pl

clean:
	rm -f *.o *.d *.def *.exe $(CLIPROG) $(GUIPROG) $(CHECKER) dlls/*
	rm -f plugins/* $(PLUGINS)
	rm -f dlls/libgtkextra-lite.dll $(imports)/*.a $(imports)/libgtkextra-lite.def
	rm -f gretl.stamp build.h
	rm -f windist/MANIFEST windist/gretl.iss 
	rm -rf windist/gretl windist/gnuplot windist/Output
	make -C updater clean
	make -C mo clean

cleandist:
	rm -f windist/MANIFEST windist/gretl.iss
	rm -rf windist/gretl windist/gnuplot windist/Output

dist: all dist.stamp


