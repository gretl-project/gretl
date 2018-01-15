PROG = gretlzip
ZLIB = -lz
CC = gcc -Wall -W -O -ansi -pedantic

$(PROG): gretlzip.c
	$(CC) -o $(PROG) gretlzip.c $(ZLIB)

clean:
	rm -f $(PROG)
