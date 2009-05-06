# Makefile for cross-compiling gretl_updater for win32 on Linux,
# using mingw compiler

include ../config.mk

CC = $(MGW_PREFIX)gcc -O2 -Wall -mms-bitfields
CFLAGS = -I../../lib/src -I../../win32/updater
SRC = update_client.c gretl_www.c gretl_untar.c
OBJ = update_client.o 

LIBS = -lz -lwsock32

CFLAGS += $(shell PKG_CONFIG_PATH=$(PKG_CONFIG_PATH) \
	pkg-config --cflags gtk+-win32-2.0)

vpath %.h ../../win32/updater ../../lib/src
vpath %.c ../../win32/updater ../../lib/src

gretl_updater.exe: $(OBJ) version.h
	 $(CC) -mwindows -s -o $@ $^ $(LIBS)

update_client.o: update_client.c gretl_www.c gretl_untar.c

.PHONY :

clean:
	rm -f *.o gretl_updater.exe


