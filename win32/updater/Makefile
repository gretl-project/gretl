CC = mingw32-gcc -W -Wall -DUPDATER -mms-bitfields
SRC = update_client.c webget.c untgz.c
OBJ = update_client.o webget.o untgz.o
CFLAGS = -I../../lib/src -I../../gui2 -I../../win32/updater
TESTLIBS = -L/usr/local/lib -lz
WINLIBS = -L../import-libs -lz.dll -lwsock32

PROG = gretl_updater.exe

vpath %.h ../../lib/src ../../win32/updater
vpath %.c ../../win32/updater ../../gui2

$(PROG): $(OBJ) version.h
	 $(CC) -mwindows -s -o $@ $(OBJ) $(WINLIBS)

testprog: $(OBJ) version.h
	$(CC) -o $@ $(OBJ) $(TESTLIBS)

.PHONY :

clean:
	rm -f *.o $(PROG) testprog

tags:
	etags *.c ../../gui2/webget.c
