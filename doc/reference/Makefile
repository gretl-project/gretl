# Copy this to your project's API docs directory and modify the variables to
# suit your project. See the GTK+ Makefiles in gtk+/docs/reference for examples
# of using the various options.

# modify for local use
HTML_DIR = /home/cottrell/stats/esl/website/API/new/gretl

srcdir = .
LIBTOOL = ../../libtool

# The name of the module.
DOC_MODULE = gretl

# Uncomment for versioned docs and specify the version of the module, e.g. '2'.
#DOC_MODULE_VERSION=2

# The top-level SGML file. Change it if you want.
DOC_MAIN_SGML_FILE = $(DOC_MODULE)-docs.sgml

# Directories containing the source code, relative to $(srcdir).
# gtk-doc will search all .c and .h files beneath these paths
# for inline comments documenting functions and macros.
# e.g. DOC_SOURCE_DIR=../../../gtk ../../../gdk
DOC_SOURCE_DIR = ../../lib/src

# Extra options to supply to gtkdoc-scan.
# e.g. SCAN_OPTIONS=--deprecated-guards="GTK_DISABLE_DEPRECATED"
SCAN_OPTIONS =

# Header files to ignore when scanning.
IGNORE_SRCS = \
        gretl_f2c.h \
	bhhh_max.c \
	bhhh_max.h \
	clapack_double.h \
	cmd_private.h \
	csvdata.c \
	csvdata.h \
        estim_private.h \
	flow_control.c \
	flow_control.h \
	gretl_bundle.c \
	gretl_bundle.h \
	gretl_intl.c \
	gretl_intl.h \
	gretl_func.c \
	gretl_func.h \
	gretl_win32.c \
	gretl_win32.h \
	gretl_www.c \
	gretl_www.h \
	interact.c \
	interact.h \
	libglue.c \
	libglue.h \
	missing_private.h \
	monte_carlo.c \
	monte_carlo.h \
	objstack.c \
	objstack.h \
	options.c \
	options.h \
	printscan.c \
	printscan.h \
	tsls.c \
	tsls.h \
	usermat.c \
	usermat.h \
	version.h

# Extra options to supply to gtkdoc-mkdb.
MKDB_OPTIONS = --output-format=xml --ignore-files="$(IGNORE_SRCS)" \
               --sgml-mode --expand-content-files=""

# Extra options to supply to gtkdoc-fixref.
FIXXREF_OPTIONS =

# Used for dependencies. The docs will be rebuilt if any of these change.
HFILE_GLOB=$(DOC_SOURCE_DIR)/*.h
CFILE_GLOB=$(DOC_SOURCE_DIR)/*.c

# Extra header to include when scanning, which are not under DOC_SOURCE_DIR
# e.g. EXTRA_HFILES=$(top_srcdir}/contrib/extra.h
EXTRA_HFILES=

# Extra SGML files that are included by $(DOC_MAIN_SGML_FILE).
content_files =

# Other files to distribute.
extra_files =

# CFLAGS and LDFLAGS for compiling scan program. Only needed if your app/lib
# contains GtkObjects/GObjects and you want to document signals and properties.
GTKDOC_CFLAGS =
GTKDOC_LIBS =

GTKDOC_CC = $(LIBTOOL) --mode=compile $(CC)
GTKDOC_LD = $(LIBTOOL) --mode=link $(CC)

# If you need to override some of the declarations, place them in this file
# and uncomment this line.
#DOC_OVERRIDES = $(DOC_MODULE)-overrides.txt

###########################################################################
# Everything below here is generic
###########################################################################

DOC_STAMPS=scan-build.stamp tmpl-build.stamp sgml-build.stamp html-build.stamp \
	pdf-build.stamp \
	$(srcdir)/tmpl.stamp $(srcdir)/sgml.stamp $(srcdir)/html.stamp \
	$(srcdir)/pdf.stamp

all-local: html-build.stamp

#### scan ####

scan-build.stamp: $(HFILE_GLOB)
	@echo '*** Scanning header files ***'
	cd $(srcdir) && \
	  gtkdoc-scan --module=$(DOC_MODULE) --source-dir=$(DOC_SOURCE_DIR) \
	--ignore-headers="$(IGNORE_SRCS)" $(SCAN_OPTIONS) $(EXTRA_HFILES)
	touch scan-build.stamp

$(DOC_MODULE)-decl.txt: scan-build.stamp
	@true

#### templates ####

tmpl-build.stamp: $(DOC_MODULE)-decl.txt $(SCANOBJ_FILES) $(DOC_MODULE)-sections.txt $(DOC_OVERRIDES)
	@echo '*** Rebuilding template files ***'
	cd $(srcdir) && gtkdoc-mktmpl --module=$(DOC_MODULE)
	touch tmpl-build.stamp

tmpl.stamp: tmpl-build.stamp
	@true

#### sgml ####

sgml-build.stamp: tmpl.stamp $(CFILE_GLOB) # $(srcdir)/tmpl/*.sgml
	@echo '*** Building SGML ***'
	cd $(srcdir) && \
	gtkdoc-mkdb --module=$(DOC_MODULE) --source-dir=$(DOC_SOURCE_DIR) \
	--main-sgml-file=$(DOC_MAIN_SGML_FILE) $(MKDB_OPTIONS)
	touch sgml-build.stamp

sgml.stamp: sgml-build.stamp
	@true

#### html ####

html-build.stamp: sgml.stamp $(DOC_MAIN_SGML_FILE) $(content_files)
	@echo '*** Building HTML ***'
	test -d $(srcdir)/html || mkdir $(srcdir)/html
	cd $(srcdir)/html && gtkdoc-mkhtml $(DOC_MODULE) ../$(DOC_MAIN_SGML_FILE)
	@echo '-- Fixing Crossreferences' 
	cd $(srcdir) && gtkdoc-fixxref --module-dir=html $(FIXXREF_OPTIONS)
	touch html-build.stamp

##############

man: sgml.stamp $(DOC_MAIN_SGML_FILE) $(content_files)
	test -d $(srcdir)/man || mkdir $(srcdir)/man
	gtkdoc-mkman 

clean:
	rm -f *~ *.bak $(SCANOBJ_FILES) *-unused.txt $(DOC_STAMPS)
	rm -f $(DOC_MODULE)-sections.txt

realclean: clean
	cd $(srcdir) && rm -rf sgml html $(DOC_MODULE)-decl-list.txt \
	$(DOC_MODULE)-decl.txt

install:
	mkdir -p $(HTML_DIR)
	rm -f $(HTML_DIR)/*.html
	(installfiles=`echo $(srcdir)/html/*.html`; \
	if test "$$installfiles" = '$(srcdir)/html/*.html'; \
	then echo '-- Nothing to install' ; \
	else \
	  for i in $$installfiles; do \
	    echo '-- Installing '$$i ; \
	    install -m 644 $$i $(HTML_DIR); \
	  done \
	fi)

