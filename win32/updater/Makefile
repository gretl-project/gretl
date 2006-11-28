# Makefile for cross-compiling gretl_updater for win32 on Linux,
# using mingw compiler

include ../config.mk

CC = $(MGW_PREFIX)gcc -O2 -Wall -mms-bitfields
CFLAGS = -DUPDATER -I../../lib/src -I../../win32/updater
SRC = update_client.c webget.c untgz.c
OBJ = update_client.o webget.o untgz.o

LIBS = -lz -lwsock32

CFLAGS += $(shell PKG_CONFIG_PATH=$(PKG_CONFIG_PATH) \
	pkg-config --cflags gtk+-win32-2.0)

vpath %.h ../../lib/src ../../win32/updater
vpath %.c ../../win32/updater ../../plugin

gretl_updater.exe: $(OBJ) version.h
	 $(CC) -mwindows -s -o $@ $^ $(LIBS)

.PHONY :

clean:
	rm -f *.o gretl_updater.exe


