top_builddir = ..

ifeq ($(prefix),)
  prefix = /opt/esl
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

# Directories
srcdir  = ./src
VPATH = ./src
share = $(prefix)/share/gretl

GTK_CFLAGS = `gtk-config --cflags`

# General compile options and libs:
PLUGINS  = stats_tables.so panel_data.so

#### End of system configuration section. ####

COMPILE = $(CC) -c $(CFLAGS) $(GTK_CFLAGS) -I../lib/src -fPIC -DPIC

all: $(PLUGINS)

.c.o:        
	$(COMPILE) $< 

stats_tables.so: stats_tables.o
	$(CC) -shared -o stats_tables.so stats_tables.o 

panel_data.so: panel_data.o
	$(CC) -shared -o panel_data.so panel_data.o 

.PHONY: 

install: $(PLUGINS) installdirs
	$(INSTALL_PROGRAM) $(PLUGINS) $(share)/plugins

install-strip: $(PLUGINS) installdirs
	$(INSTALL_PROGRAM) -s $(PLUGINS) $(share)/plugins

installdirs:
	../tools/mkinstalldirs $(share)/plugins

clean:
	rm -f *.o $(PLUGINS) core 





