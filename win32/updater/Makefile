CC = mingw32-gcc -O2 -W -Wall -DUPDATER -DNEWICC -mms-bitfields
SRC = update_client.c webget.c untgz.c
OBJ = update_client.o webget.o untgz.o
CFLAGS = -I../../lib/src -I../../gui2 -I../../win32/updater
WINLIBS = -lz.dll -lwsock32 -lcomctl32

vpath %.h ../../lib/src ../../win32/updater
vpath %.c ../../win32/updater ../../gui2

gretl_updater.exe: $(OBJ) version.h
	 $(CC) -mwindows -s -o $@ $(OBJ) $(WINLIBS)

.PHONY :

clean:
	rm -f *.o gretl_updater.exe testprog

tags:
	etags *.c ../../gui2/webget.c
