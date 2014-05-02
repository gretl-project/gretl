# This makefile can be used for cross-compiling gretl for Windows
# on Linux, or for a native build on MS Windows.  In both cases we
# assume that the mingw compiler is used.  The build type is
# controlled by the BUILDTYPE variable in config.mk, included
# below.

include config.mk

# directories
topsrc = ..
libsrcdir  = $(topsrc)/lib/src
plugindir = $(topsrc)/plugin
zipdir = $(topsrc)/plugin/zipunzip
guidir = $(topsrc)/gui2
imports = ./import-libs
host_builddir = $(topsrc)/build

# tools (can add -fopenmp to CC if wanted)
CC = $(MGW_PREFIX)gcc -msse2 -fopenmp -Wall -O2 -mms-bitfields -DWIN32
AS = $(MGW_PREFIX)as
DLLWRAP = $(MGW_PREFIX)dllwrap
DLLWRAP_FLAGS = --as=$(AS) --export-all --driver-name $(CC) -s
RC = $(MGW_PREFIX)windres
RCFLAGS = --define __WIN32__ --define __WIN95__ --define MSRC \
	--define __GNUWIN32__

# libraries
GTKLIBS := $(shell PKG_CONFIG_PATH=$(PKG_CONFIG_PATH) pkg-config --libs gtk+-win32-2.0)
GLIBLIB := $(shell PKG_CONFIG_PATH=$(PKG_CONFIG_PATH) pkg-config --libs glib-2.0)
RL_LIBS = -lreadline -lhistory
GUI_LIBS = -lwsock32 -lole32 $(GTKLIBS) -lgtksourceview-2.0 -lxml2.dll -L$(imports)
GRETLLIB = -L$(imports) -lgretl -lintl
LIBS = -lm $(GRETLLIB) -ladvapi32
LAPACK_LIBS = -llapack.dll -lblas.dll
FFTW_LIBS = -lfftw3
CURL_LIBS = -lcurl.dll

# flags/defines
CFLAGS = -I. -I$(libsrcdir) -I$(plugindir) -I$(guidir) -I$(zipdir) -I$(XML2_INC) $(RLIB_CFLAGS) 
GTK_CFLAGS := $(shell PKG_CONFIG_PATH=$(PKG_CONFIG_PATH) pkg-config --cflags gtk+-win32-2.0 gtksourceview-2.0)

ifeq ($(BUILDTYPE),cross)
  GRETL_LANG = cross_lang
  GRETL_NEWS = cross_news
  BUILDDATE = cross_builddate
else
  GRETL_LANG = native_lang
  GRETL_NEWS = native_news
  BUILDDATE = native_builddate
endif

ifeq ($(have_mpi),yes)
  MPICLI = gretlcli-mpi.exe
  MPI_LIB = -lmsmpi
endif

# source paths
vpath %.c $(libsrcdir)
vpath %.c $(topsrc)/cephes
vpath %.c $(topsrc)/minpack
vpath %.c $(topsrc)/dcmt
vpath %.c $(topsrc)/cli
vpath %.c $(guidir)
vpath %.c $(plugindir)
vpath %.c $(plugindir)/libole2
vpath %.c $(plugindir)/zipunzip
vpath %.c $(plugindir)/rq
vpath %.c $(topsrc)/tests
vpath %.h $(libsrcdir):$(topsrc)/cli:$(guidir):$(plugindir):.
vpath %.h $(topsrc)/cephes:..
vpath %.rc $(topsrc)/win32

#### End of system configuration section. ####

PLUGINSRC = \
	arbond.c \
	excel_import.c \
	xlsx_import.c \
	gnumeric_import.c \
	eviews_import.c \
	stata_import.c \
	spss_import.c \
	sas_import.c \
	jmulti_import.c \
	ods_import.c \
	johansen.c \
	jrestrict.c \
	jalpha.c \
	kernel.c \
	vif.c \
	leverage.c \
	mp_ols.c \
	pca.c \
	progress_bar.c \
	range-mean.c \
	fractals.c \
	poisson.c \
	duration.c \
	stats_tables.c \
	sysest.c \
	fiml.c \
	liml.c \
	tramo_options.c \
	tramo-x12a.c \
	nistcheck.c \
	arma_cml.c \
	arma_init.c \
	arma.c \
        arma_x12.c \
	workbook.c \
	mod_garch.c \
	fcp.c \
	garch.c \
        audio.c \
	urcdist.c \
	gretlzip.c \
	heckit.c \
	odbc_import.c \
	quantreg.c \
	interval.c \
	anova.c \
	biprobit.c \
	reprobit.c

PLUGINOBJ = $(PLUGINSRC:.c=.o)

PLUGINS = plugins/stats_tables.dll \
	plugins/progress_bar.dll \
	plugins/gnumeric_import.dll \
	plugins/eviews_import.dll \
	plugins/excel_import.dll \
	plugins/xlsx_import.dll \
	plugins/stata_import.dll \
	plugins/spss_import.dll \
	plugins/sas_import.dll \
	plugins/jmulti_import.dll \
	plugins/ods_import.dll \
	plugins/mp_ols.dll \
	plugins/sysest.dll \
	plugins/johansen.dll \
	plugins/kernel.dll \
	plugins/arbond.dll \
	plugins/arma.dll \
	plugins/arma_x12.dll \
	plugins/range-mean.dll \
	plugins/fractals.dll \
	plugins/poisson.dll \
	plugins/duration.dll \
	plugins/tramo-x12a.dll \
	plugins/leverage.dll \
	plugins/pca.dll \
	plugins/nistcheck.dll \
	plugins/garch.dll \
        plugins/audio.dll \
	plugins/vif.dll \
	plugins/urcdist.dll \
	plugins/gretlzip.dll \
	plugins/heckit.dll \
	plugins/odbc_import.dll \
	plugins/quantreg.dll \
	plugins/interval.dll \
	plugins/anova.dll \
	plugins/interpolate.dll \
	plugins/biprobit.dll \
	plugins/reprobit.dll \
	plugins/panurc.dll

LIBSRC = adf_kpss.c \
	bhhh_max.c \
	bootstrap.c \
	boxplots.c \
	calendar.c \
	compare.c \
	compat.c \
	csvdata.c \
	dataio.c \
	dataset.c \
	dbread.c \
	dbwrite.c \
	describe.c \
        discrete.c \
	estimate.c \
	flow_control.c \
	forecast.c \
	geneval.c \
	genfuncs.c \
	genlex.c \
	genmain.c \
	gensyntax.c \
        gmm.c \
	graphing.c \
	gretl_bfgs.c \
	gretl_bundle.c \
	gretl_commands.c \
	gretl_data_io.c \
	gretl_errors.c \
	gretl_fft.c \
	gretl_foreign.c \
	gretl_func.c \
	gretl_intl.c \
	gretl_list.c \
	gretl_matrix.c \
	gretl_model.c \
	gretl_panel.c \
	gretl_paths.c \
	gretl_prn.c \
	gretl_restrict.c \
	gretl_string_table.c \
	gretl_untar.c \
	gretl_utils.c \
	gretl_win32.c \
	gretl_www.c \
	gretl_xml.c \
	gretl_zip.c \
	interact.c \
	kalman.c \
	libglue.c \
	libset.c \
	matrix_extra.c \
	missing.c \
	modelprint.c \
	monte_carlo.c \
	nls.c \
	nonparam.c \
	objstack.c \
	options.c \
	plotspec.c \
	plugins.c \
	printout.c \
	printscan.c \
	pvalues.c \
	qr_estimate.c \
	random.c \
	strutils.c \
	subsample.c \
	system.c \
	texprint.c \
	transforms.c \
	tsls.c \
        usermat.c \
	uservar.c \
	var.c \
	varprint.c \
	vartest.c \
	irfboot.c
	
ifeq ($(have_mpi),yes)
  LIBSRC += gretl_mpi.c
endif

LIBOBJ = $(LIBSRC:.c=.o)

CLISRC = gretlcli.c complete.c 

CLIOBJ = $(CLISRC:.c=.o)

GUISRC = about.c \
	calculator.c \
	callbacks.c \
	cmdstack.c \
	console.c \
	database.c \
	datafiles.c \
	datawiz.c \
	dialogs.c \
	dlgutils.c \
	filelists.c \
	fileselect.c \
	filters.c \
	fncall.c \
	fnsave.c \
	gpt_control.c \
	gpt_dialog.c \
	graphics.c \
	graph_page.c \
	gretl.c \
	gretlwin32.c \
	gretl_ipc.c \
	guiprint.c \
	gui_recode.c \
	gui_utils.c \
	helpfiles.c \
	lagpref.c \
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
	tabwin.c \
	textbuf.c \
	textutil.c \
	toolbar.c \
	treeutils.c \
	update.c \
	varinfo.c \
	winstack.c

ifeq ($(HAVE_GTK_SPINNER),no)
  GUISRC += spinner.c
endif

GUIOBJ = $(GUISRC:.c=.o)

PROBSRC = airy.c const.c hyp2f1.c igami.c j1.c kn.c pdtr.c struve.c \
	bdtr.c expx2.c hyperg.c incbet.c jn.c mtherr.c polevl.c unity.c \
	btdtr.c fdtr.c i0.c incbi.c jv.c  nbdtr.c polrt.c yn.c \
	chbevl.c gamma.c i1.c iv.c k0.c kv.c ndtr.c psi.c \
	chdtr.c gdtr.c igam.c j0.c k1.c ndtri.c stdtr.c

PROBOBJ = $(PROBSRC:.c=.o)

MINSRC = chkder.c enorm.c fdjac2.c lmder1.c lmder.c lmdif.c lmpar.c \
	qrfac.c qrsolv.c lbfgsb.c
MINOBJ = $(MINSRC:.c=.o)

DCMTSRC = check32.c eqdeg.c genmtrand.c mt19937.c prescr.c seive.c
DCMTOBJ = $(DCMTSRC:.c=.o)

OBJS = $(LIBOBJ) $(CLIOBJ) $(GUIOBJ) $(PLUGINOBJ)

MOFILES = mo/de.po mo/el.mo mo/es.mo mo/eu.mo mo/fr.mo mo/it.mo mo/pl.mo \
	mo/pt.mo mo/ru.mo mo/rw.mo mo/tr.mo mo/cs.mo mo/sq.mo

%.o: %.c
	$(CC) -c $(GTK_CFLAGS) $(CFLAGS) $<
	$(CC) -MM $(GTK_CFLAGS) $(CFLAGS) $< > $*.d 

# targets begin

all: buildstamp dirs dlls/libprob.dll dlls/libgretl.dll \
	$(PLUGINS) gretlcli.exe gretl.exe $(MPICLI) $(MOFILES) $(GRETL_LANG) $(GRETL_NEWS)

gretlcli.exe: $(CLIOBJ)
	$(CC) -o $@ $^ $(LIBS) $(RL_LIBS) -lz -lmingwex -s

gretlcli-mpi.exe: gretlcli-mpi.o
	$(CC) -o $@ $< $(MPI_LIB) $(LIBS) -lz -lmingwex -s

gretl.exe: $(GUIOBJ) gretlres.o
	$(CC) -mwindows -o $@ $^ $(LIBS) $(GUI_LIBS) -lz -lmingwex -s

gretlres.o: gretl.rc 
	$(RC) $(RCFLAGS) $< $@

cross_lang: $(libsrcdir)/options.c $(libsrcdir)/genlex.c $(libsrcdir)/gretl_commands.c
	HERE=$$PWD ; cd $(host_builddir)/gui2 && make mklang && ./mklang > $$HERE/gretl.lang ; cd $$HERE

cross_news: $(topsrc)/ChangeLog
	HERE=$$PWD ; cd $(host_builddir)/gui2 && make mknews && cd $$HERE && $(topsrc)/build/gui2/mknews $< > NEWS

mklang.o: mklang.c $(libsrcdir)/options.c $(libsrcdir)/genlex.c $(libsrcdir)/gretl_commands.c
	$(CC) -c -I. -I$(libsrcdir) -I$(XML2_INC) $<

mklang.exe: mklang.o
	$(CC) -o $@ $< $(GRETLLIB) -lxml2.dll -lz

mknews.exe: mknews.c
	$(CC) -o $@ $<

native_lang: mklang.exe
	PATH=/bin:./dlls:/mingw/winbuild/misc-dll ./mklang.exe > gretl.lang

native_news: mknews.exe
	./mknews.exe $(topsrc)/ChangeLog > NEWS

$(MOFILES):
	make -C mo

dlls/libgretl.dll: $(LIBOBJ) $(MINOBJ) $(DCMTOBJ) dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) \
	--output-def libgretl.def --implib $(imports)/libgretl.a \
	-o $@ $^ -lm -L$(imports) -lxml2.dll -lz -lintl -lprob -lgmp.dll \
	-lmingwex -lwsock32 $(GLIBLIB) $(LAPACK_LIBS) $(FFTW_LIBS) $(CURL_LIBS)

dlls/libprob.dll: $(PROBOBJ) dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) \
	--output-def libprob.def --implib $(imports)/libprob.a \
	-o $@ $^ 

plugins/stats_tables.dll: stats_tables.o dllinit.o 
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ -lz $(GRETLLIB)

plugins/range-mean.dll: range-mean.o dllinit.o 
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB)

plugins/fractals.dll: fractals.o dllinit.o 
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB)

plugins/poisson.dll: poisson.o dllinit.o 
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) -lprob

plugins/duration.dll: duration.o dllinit.o 
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB)

plugins/kernel.dll: kernel.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) 

plugins/vif.dll: vif.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) $(LAPACK_LIBS)

plugins/urcdist.dll: urcdist.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ -lz -lm -lintl $(GRETLLIB)

plugins/tramo-x12a.dll: tramo-x12a.o tramo_options.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) $(GTKLIBS)

plugins/progress_bar.dll: progress_bar.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GTKLIBS) -L$(imports) -lintl

plugins/gnumeric_import.dll: gnumeric_import.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ -L$(imports) -lxml2.dll \
	$(GRETLLIB) $(GTKLIBS)

plugins/ods_import.dll: ods_import.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ -L$(imports) -lxml2.dll \
	$(GRETLLIB) $(GTKLIBS)

plugins/xlsx_import.dll: xlsx_import.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ -L$(imports) -lxml2.dll \
	$(GRETLLIB) $(GTKLIBS)

plugins/excel_import.dll: excel_import.o workbook.o ms-ole.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) $(GTKLIBS) -lmoldname

plugins/eviews_import.dll: eviews_import.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB)

plugins/stata_import.dll: stata_import.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GLIBLIB) $(GRETLLIB)

plugins/spss_import.dll: spss_import.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GLIBLIB) $(GRETLLIB)

plugins/sas_import.dll: sas_import.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB)

plugins/jmulti_import.dll: jmulti_import.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB)

plugins/mp_ols.dll: mp_ols.o dllinit.o 
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) -lmpfr.dll -lgmp.dll

plugins/sysest.dll: sysest.o fiml.o liml.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) $(LAPACK_LIBS)

plugins/johansen.dll: johansen.o jrestrict.o jalpha.o dllinit.o
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

plugins/arbond.dll: arbond.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) $(LAPACK_LIBS)

plugins/arma.dll: arma.o arma_cml.o arma_init.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) $(LAPACK_LIBS) -lprob

plugins/arma_x12.dll: arma_x12.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) $(GTKLIBS) $(LAPACK_LIBS)

plugins/garch.dll: garch.o fcp.o mod_garch.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) $(LAPACK_LIBS)

plugins/audio.dll: audio.o midi_utils.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) $(GTKLIBS) -lole32 -lsapi

plugins/gretlzip.dll: gretlzip.o zfileio.o zipfile.o zipmain.o zipwork.o \
	zsystem.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GLIBLIB) -lz $(GRETLLIB)

plugins/heckit.dll: heckit.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB)

plugins/odbc_import.dll: odbc_import.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) -lodbc32

plugins/quantreg.dll: quantreg.o rqfnb.o rqbr.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) $(LAPACK_LIBS)

plugins/interval.dll: interval.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) 

plugins/anova.dll: anova.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) $(GLIBLIB) 

plugins/interpolate.dll: interpolate.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB)

plugins/biprobit.dll: biprobit.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB)

plugins/reprobit.dll: reprobit.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB)

plugins/panurc.dll: panurc.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) $(GLIBLIB)

ms-ole.o: ms-ole.c $(topsrc)/plugin/libole2/ms-ole.h
	$(CC) -c $(CFLAGS) $(GTK_CFLAGS) $<

workbook.o: workbook.c $(topsrc)/plugin/libole2/ms-ole.h
	$(CC) -c $(CFLAGS) $(GTK_CFLAGS) $<

midi_utils.o: midi_utils.c $(topsrc)/plugin/midi_utils.h
	$(CC) -c $(CFLAGS) $(GTK_CFLAGS) $<

-include $(OBJS:.o=.d)

# Distribution executable

dist.stamp:
	bash ./mkwindist $(BUILDTYPE)

.PHONY:

dirs:
	mkdir -p dlls 
	mkdir -p plugins 
	mkdir -p import-libs

cross_builddate: $(topsrc)/builddate.c
	gcc -Wall -o builddate $<

native_builddate: $(topsrc)/builddate.c
	$(CC) -Wall -o builddate $<

buildstamp: $(BUILDDATE)
	./builddate

clean:
	rm -f *.o *.d *.def *.exe $(CLIPROG) $(GUIPROG) $(MPICLI) $(CHECKER) dlls/*
	rm -f plugins/* $(PLUGINS)
	rm -f $(imports)/*.a 
	rm -f gretl.stamp gretl.lang builddate build.h
	rm -f windist/MANIFEST windist/gretl.iss 
	rm -rf windist/gretl windist/Output
	make -C mo clean
	make -C db clean

cleandist:
	rm -f windist/MANIFEST windist/gretl.iss
	rm -rf windist/gretl windist/gnuplot windist/Output

dist: all dist.stamp

