CC = gcc -Wall -O2

XML_CFLAGS := $(shell pkg-config --cflags libxml-2.0)
XML_LIBS := $(shell pkg-config --libs libxml-2.0)

frb2db: frb2db.c
	$(CC) $(XML_CFLAGS) -o $@ $< $(XML_LIBS)

idx: frb2db 
	./frb2db H15_data.xml

clean:
	rm -f frb2db *.idx *.bin
