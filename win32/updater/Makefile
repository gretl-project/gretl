CC = mingw32-gcc -W -Wall -DOS_WIN32 -mms-bitfields
# CC = gcc -Wall -DWDEBUG
SRC = update_client.c webget.c untgz.c
OBJ = update_client.o webget.o untgz.o
CFLAGS = -I../include -I../../lib/src
TESTLIBS = -L/usr/local/lib -lz
WINLIBS = -L../import-libs -lz.dll -lwsock32

PROG = gretl_updater.exe

vpath %.h ../../lib/src

$(PROG): $(OBJ) version.h
	 $(CC) -mwindows -s -o $@ $(OBJ) $(WINLIBS)

testprog: $(OBJ) version.h
	$(CC) -o $@ $(OBJ) $(TESTLIBS)

.PHONY :

clean:
	rm -f *.o $(PROG) testprog

tags:
	etags *.c
