# tools
CC = mingw32-gcc -fnative-struct -Wall
AS = mingw32-as
DLLWRAP = mingw32-dllwrap
DLLTOOL = mingw32-dlltool

# outputs
CLIPROG = gretlcli.exe
GUIPROG = gretlw32.exe

# libraries
LIBS = -lm -L./libgretl -lgretl -ladvapi32
GTKLIBS = -L./a_libs -lglib-2.0.dll -lgtk.dll -lgdk.dll
ZLIB_LIB = -lz.dll
CLI_LIBS = -L./a_libs -lreadline.dll
GUI_LIBS = $(GTKLIBS) -lgtkextra-lite -lwsock32 #-lhhctrl
GRETLLIB = -L./libgretl -lgretl
PIXBUFLIB = -L./a_libs -lgdk-pixbuf -lpng -lz.dll

# flags/defines
LDFLAGS = -s 
GLIBINC = -I./include/glib-2.0 -I$(extradir) 
CFLAGS = -I. -I$(libdir) -I./include -DHAVE_READLINE -DOS_WIN32 $(GLIBINC)

# Directories
gretldir = /home/cottrell/stats/esl/gretl
libdir  = $(gretldir)/lib/src
extradir = ./gtkextra-lib
vpath %.c $(libdir)
vpath %.c $(gretldir)/cli
vpath %.c $(gretldir)/gui
vpath %.c $(extradir)
vpath %.c $(gretldir)/plugin
vpath %.c $(gretldir)/plugin/libole2
vpath %.h $(libdir):$(gretldir)/cli:$(gretldir)/gui:$(extradir):$(gretldir)/plugin

#### End of system configuration section. ####

WRES    = gretlw32_res.o
RC      = mingw32-windres
RCFLAGS = --include-dir /usr/local/cross-tools/mingw32/include \
	--define __WIN32__ --define __WIN95__ --define MSRC \
	--define __GNUWIN32__
RCOUT   = gretlw32_res.o
RES2COFF = echo gretlw32_res.o

LIBSRC = estimate.c compare.c utils.c errors.c pvalues.c dataio.c \
       generate.c strutils.c describe.c printout.c texprint.c var.c \
       interact.c graphing.c win32.c monte_carlo.c nonparam.c discrete.c \
       subsample.c calendar.c
LIBHDR = libgretl.h commands.h estimate.h utils.h errors.h pvalues.h \
       compare.h dataio.h generate.h strutils.h describe.h printout.h \
       version.h texprint.h var.h interact.h graphing.h monte_carlo.h \
       nonparam.h discrete.h subsample.h calendar.h win32.h
LIBOBJ = $(LIBSRC:.c=.o)

CLISRC = gretlcli.c complete.c 
CLIOBJ = $(CLISRC:.c=.o)

GUISRC = gretl.c library.c gui_utils.c callbacks.c dialogs.c ssheet.c \
	datafiles.c database.c gpt_control.c session.c webget.c \
	calculator.c fileselect.c htmlprint.c guiprint.c boxplots.c
GUIHDR = gretl.h library.h gui_utils.h callbacks.h dialogs.h session.h \
	webget.h htmlprint.h guiprint.h
GUIOBJ = $(GUISRC:.c=.o)

EXTRASRC = gtkitementry.c gtksheet.c gtkiconlist.c \
	gtkplotpc.c gtkplotps.c gtkpsfont.c gtkplot.c
EXTRAHDR = $(EXTRASRC:.c=.h)
EXTRAOBJ = $(EXTRASRC:.c=.o)

# targets begin
all: dll extra-dll plugin $(CLIPROG) $(GUIPROG) 

$(CLIPROG): $(CLIOBJ) 
	$(CC) -s -o $@ $^ $(LIBS) $(CLI_LIBS) $(ZLIB_LIB)

$(GUIPROG): $(GUIOBJ) $(WRES) 
	$(CC) -mwindows -s -o $@ $^ $(LIBS) $(GUI_LIBS) $(ZLIB_LIB)

$(CLIOBJ): $(CLI_SRC) 

$(GUIOBJ): $(GUISRC) $(GUIHDR) # guistamp

guistamp:
	./builddate.pl > build.h

$(EXTRAOBJ): $(EXTRAHDR) $(EXTRASRC) 

$(LIBOBJ): $(LIBSRC) $(LIBHDR) 

$(WRES): gretlw32.rc 
	$(RC) $(RCFLAGS) gretlw32.rc $(RCOUT)
	$(RES2COFF)

windist: $(CLIPROG) $(GUIPROG) dll extra-dll plugin
	./mkwindist 

topclean:
	rm -f *.o *.exe

clean:
	rm -f *.o $(CLIPROG) $(GUIPROG) dlls/libgretl.dll libgretl/* 
	rm -f dlls/stats_tables.dll dlls/panel_data.dll plugins/*
	rm -f dlls/libgtkextra-lite.dll

# make a dll of libgretl...

dll: dlls/libgretl.dll
extra-dll: dlls/libgtkextra-lite.dll
plugin: dlls/stats_tables.dll dlls/panel_data.dll dlls/longname.dll \
	dlls/progress_bar.dll dlls/gnumeric_import.dll dlls/excel_import.dll

DLLWRAP_FLAGS = --as=$(AS) --export-all --driver-name $(CC) -s

dlls/libgretl.dll: $(LIBOBJ) dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) \
	--output-def libgretl/libgretl.def \
	--implib libgretl/libgretl.a \
	-o $@ $^ -lm -L./a_libs -lxml -lz

dlls/libgtkextra-lite.dll: $(EXTRAOBJ) dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) \
	--output-def a_libs/libgtkextra-lite.def \
	--implib a_libs/libgtkextra-lite.a \
	-o $@ $^ $(GTKLIBS)

dlls/stats_tables.dll: stats_tables.o dllinit.o 
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB)

dlls/panel_data.dll: panel_data.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) $(GTKLIBS)

dlls/progress_bar.dll: progress_bar.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GTKLIBS)

dlls/gnumeric_import.dll: gnumeric_import.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ -L./a_libs -lxml \
	$(GRETLLIB) $(GTKLIBS)

dlls/excel_import.dll: excel_import.o workbook.o ms-ole.o dllinit.o
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ $(GRETLLIB) $(GTKLIBS) -lmoldname

dlls/longname.dll: longname.o 
	$(DLLWRAP) $(DLLWRAP_FLAGS) --output-def $(<:.o=.def) -o $@ \
	--implib plugins/$(<:.o=.a) $^ -lkernel32

panel_data.o: panel_data.c 
	$(CC) -c $(CFLAGS) $(GLIBINC) $<

OLEFLAGS = -I$(gretldir)/plugin

ms-ole.o: ms-ole.c $(gretldir)/plugin/libole2/ms-ole.h
	$(CC) -c $(CFLAGS) $(OLEFLAGS) $(GLIBINC) $<

workbook.o: workbook.c $(gretldir)/plugin/libole2/ms-ole.h
	$(CC) -c $(CFLAGS) $(OLEFLAGS) $(GLIBINC) $<
