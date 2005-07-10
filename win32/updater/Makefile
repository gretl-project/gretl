# Makefile for cross-compiling gretl_updater for win32 on Linux,
# using mingw compiler

include ../config.mk

CC = $(MGW_PREFIX)gcc -O2 -Wall -mms-bitfields
CFLAGS = -DUPDATER -I../../lib/src -I../../gui2 -I../../win32/updater
SRC = update_client.c webget.c untgz.c progress_bar.c
OBJ = update_client.o webget.o untgz.o progress_bar.o

GTKLIBS = $(shell PKG_CONFIG_PATH=$(PKG_CONFIG_PATH) \
	pkg-config --libs gtk+-win32-2.0)
LIBS = -lz -lwsock32 $(GTKLIBS)

CFLAGS += $(shell PKG_CONFIG_PATH=$(PKG_CONFIG_PATH) \
	pkg-config --cflags gtk+-win32-2.0)

vpath %.h ../../lib/src ../../win32/updater
vpath %.c ../../win32/updater ../../gui2 ../../plugin

gretl_updater.exe: $(OBJ) version.h
	 $(CC) -mwindows -s -o $@ $(OBJ) $(LIBS)

.PHONY :

clean:
	rm -f *.o gretl_updater.exe


