# Copy this to your project's API docs directory and modify the variables to
# suit your project. See the GTK+ Makefiles in gtk+/docs/reference for examples
# of using the various options.

# modify for local use
HTML_DIR = /home/cottrell/stats/esl/website/API

srcdir = .
LIBTOOL = ../../libtool

# The name of the module.
DOC_MODULE = gretl

# The top-level SGML file. Change it if you want.
DOC_MAIN_SGML_FILE = $(DOC_MODULE)-docs.sgml

# The directory containing the source code. Relative to $(srcdir).
# gtk-doc will search all .c & .h files beneath here for inline comments
# documenting functions and macros.
DOC_SOURCE_DIR = ../../lib/src

# Extra options to supply to gtkdoc-scan.
SCAN_OPTIONS =

# Header files to ignore when scanning.
IGNORE_HFILES = \
        gretl_f2c.h \
	clapack_double.h \
	cmd_private.h \
        estim_private.h \
	missing_private.h \
	version.h

# Extra options to supply to gtkdoc-mkdb.
MKDB_OPTIONS = --output-format=xml --ignore-files="$(IGNORE_HFILES)" \
               --expand-content-files=""

# Extra options to supply to gtkdoc-fixref.
FIXXREF_OPTIONS =

# Used for dependencies.
HFILE_GLOB=
CFILE_GLOB=

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
	   $(srcdir)/tmpl.stamp $(srcdir)/sgml.stamp $(srcdir)/html.stamp

all-local: html-build.stamp

#### scan ####

scan-build.stamp: $(HFILE_GLOB)
	@echo '*** Scanning header files ***'
	cd $(srcdir) && \
	  gtkdoc-scan --module=$(DOC_MODULE) --source-dir=$(DOC_SOURCE_DIR) \
	--ignore-headers="$(IGNORE_HFILES)" $(SCAN_OPTIONS) $(EXTRA_HFILES)
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
	(installfiles=`echo $(srcdir)/html/*.html`; \
	if test "$$installfiles" = '$(srcdir)/html/*.html'; \
	then echo '-- Nothing to install' ; \
	else \
	  for i in $$installfiles; do \
	    echo '-- Installing '$$i ; \
	    install -m 644 $$i $(HTML_DIR); \
	  done \
	fi)

