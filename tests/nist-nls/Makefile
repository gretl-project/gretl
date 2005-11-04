CC = gcc -g -Wall
CFLAGS = `pkg-config --cflags gretl`
LIBS = `pkg-config --libs gretl`

nist-nls-test: nist-nls-test.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)

clean:
	rm -f nist-nls-test

distclean: clean
	rm -f *.dat *.out
