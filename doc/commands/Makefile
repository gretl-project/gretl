CC = gcc -O2 -Wall

# These flags may be have to changed depending on the system
XML_CFLAGS = `pkg-config --cflags libxml-2.0 libxslt`
XML_LIBS = `pkg-config --libs libxml-2.0 libxslt`

# and these too
GRETL_CFLAGS = `gretl-config --cflags`
GRETL_LIBS = `gretl-config --libs`

progs: formatter skeleton reflow

formatter: formatter.c 
	$(CC) $(XML_CFLAGS) -o $@ $^ $(XML_LIBS)

skeleton: skeleton.c
	$(CC) $(GRETL_CFLAGS) -o $@ $^ $(GRETL_LIBS)

reflow: dolines.c
	$(CC) -o $@ $^ 

# test generation of XML skeleton from libgretl info:
# the skeleton could be used as a template for translation
skeleton.xml: skeleton gretl_commands.dtd
	./skeleton > $@

# test validity of skeleton.xml, and ability of formatter
# to generate (incomplete, as yet) formatted output (text and XML)
test: skeleton.xml formatter gretlman.xsl gretltxt.xsl
	xmllint --valid --noout $< && ./formatter skeleton.xml

# work in progress: format and reflow the English command reference
sample: formatter reflow inprogress.xml gretlman.xsl gretltxt.xsl
	./formatter inprogress.xml
	./reflow < commands.txt > commands.hlp

clean:
	rm -f formatter skeleton reflow skeleton.xml *~
