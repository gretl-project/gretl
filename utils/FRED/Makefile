CC = gcc -Wall -O2

CURL_CFLAGS := $(shell curl-config --cflags)
CURL_LIBS := $(shell curl-config --libs)
XML_CFLAGS := $(shell pkg-config --cflags libxml-2.0)
XML_LIBS := $(shell pkg-config --libs libxml-2.0)

CFLAGS = $(CURL_CFLAGS) $(XML_CFLAGS)
LIBS = $(CURL_LIBS) $(XML_LIBS) -lm

fred2db: fred2db.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)

fedstl.idx: fred2db
	./fred2db

clean:
	rm -f fred2db *.idx *.bin
