CC = mingw32-gcc -W -Wall -fnative-struct -DOS_WIN32
SRC = update_client.c webget.c untgz.c
OBJ = update_client.o webget.o untgz.o
CFLAGS = -I../include
LIBS = -L../import-libs -lz.dll -lwsock32
libsrc = ../../lib/src

PROG = gretl_updater.exe

$(PROG): $(OBJ) 
	 $(CC) -mwindows -s -o $@ $^ $(LIBS)

$(OBJ): version.h

version.h:
	cp $(libsrc)/$@ .

.PHONY :

clean:
	rm -f *.o $(PROG) rm -f version.h

tags:
	etags *.c
