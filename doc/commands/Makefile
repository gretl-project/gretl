CC = gcc -O2 -Wall

CFLAGS = `pkg-config --cflags libxml-2.0`
XML_LIBS = `pkg-config --libs libxml-2.0`

GRETL_CFLAGS = `pkg-config --cflags gretl`
GRETL_LIBS = `pkg-config --libs gretl`

XMLDCL = /opt/jade/share/xml.dcl

all: formatter skeleton

formatter: formatter.o xml_print.o
	$(CC) -o $@ $^ $(XML_LIBS)

skeleton: skeleton.c
	$(CC) $(GRETL_CFLAGS) -o $@ $^ $(GRETL_LIBS)

# test generation of XML skeleton from libgretl info
# the skeleton could be used as a template for translation
skeleton.xml: skeleton gretl_commands.dtd
	./skeleton > $@

# test validity of skeleton.xml, and ability of formatter
# to generate dummy 'formatted' output (text, and XML to trial.xml)
testprint.txt: skeleton.xml formatter
	nsgmls -wxml -s $(XMLDCL) $<
	./formatter skeleton.xml > $@

clean:
	rm -f formatter skeleton skeleton.xml testprint.txt \
	trial.xml *~ *.o
