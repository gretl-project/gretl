# Makefile for partial build of libole2 for gretl's excel importer

CC = gcc

GLIB_CFLAGS = `glib-config --cflags`

INCLUDES = -I.. $(GLIB_CFLAGS)
CFLAGS = -g -O2 -Wall -DHAVE_UNISTD_H
COMPILE = $(CC) $(CFLAGS) $(INCLUDES)

.c.o:
	$(COMPILE) -c $<

all: ms-ole.o 

clean:
	rm -f *.o 
