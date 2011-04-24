CC = gcc -g -Wall -O2
CFLAGS = `pkg-config --cflags gretl`
LIBS = `pkg-config --libs gretl`

test: nist-nls-test run_tests.sh
	-./run_tests.sh
	
testa: nist-nls-test run_tests.sh
	-./run_tests.sh -a

nist-nls-test: nist-nls-test.o
	$(CC) -o $@ $< $(LIBS)
	
nist-nls-test.o: nist-nls-test.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f nist-nls-test *.o
	
distclean: clean
	rm -f *.dat errlog.* out.*
	
depend:
	makedepend $(CFLAGS) -Y nist-nls-test.c	
	
# DO NOT DELETE
