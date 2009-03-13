topsrc = ../..

HOSTCC = gcc -Wall
BCIHDIR = $(topsrc)/share/bcih

vpath %.c $(BCIHDIR)

fedstl.bin: mkbin $(BCIHDIR)/fedstl.dat
	./mkbin $(BCIHDIR) 

mkbin: mkbin.c
	$(HOSTCC) -o $@ $<

clean:
	rm -f mkbin fedstl.bin
