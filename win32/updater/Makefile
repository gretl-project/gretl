CC = mingw32-gcc -W -Wall -DOS_WIN32 -mms-bitfields
SRC = update_client.c webget.c untgz.c
OBJ = update_client.o webget.o untgz.o
CFLAGS = -I../include
LIBS = -L../import-libs -lz.dll -lwsock32
libsrc = ../../lib/src

PROG = gretl_updater.exe

vpath %.h $(libsrc)

$(PROG): $(OBJ) version.h
	 $(CC) -mwindows -s -o $@ $(OBJ) $(LIBS)

.PHONY :

clean:
	rm -f *.o $(PROG)

tags:
	etags *.c
