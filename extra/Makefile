ifeq ($(CC),)
  CC = gcc
endif
ifeq ($(CFLAGS),)
  CFLAGS = -g -O2 -Wall
endif
ifeq ($(MAKE),)
  MAKE = make
endif

all: client

COMPILE = $(CC) -c $(CFLAGS) `gretl-config --cflags`

.c.o:
	$(COMPILE) $<

client: simple_client.o Makefile
	$(CC) -o client simple_client.o `gretl-config --libs`

clean:
	rm -f *.o ./client core 
