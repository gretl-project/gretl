CC = gcc -O2 -Wall
XML_CFLAGS = -I/usr/include/libxml2  
XML_LIBS = -lxml2 -lpthread -lz -lm 

formatter: formatter.c
	$(CC) $(XML_CFLAGS) -o $@ $^ $(XML_LIBS)

clean:
	rm -f formatter
