topsrc = ../..

HOSTCC = gcc -Wall
BCIHDIR = $(topsrc)/share/bcih

fedstl.bin: mkbin $(BCIHDIR)/fedstl.dat
	./mkbin $(BCIHDIR) 

mkbin: mkbin.c
	$(HOSTCC) -o $@ $<

clean:
	rm -f mkbin fedstl.bin
