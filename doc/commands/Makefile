CC = gcc -O2 -Wall

XML_CFLAGS = `pkg-config --cflags libxml-2.0 libxslt`
XML_LIBS = `pkg-config --libs libxml-2.0 libxslt`

GRETL_CFLAGS = `gretl-config --cflags`
GRETL_LIBS = `gretl-config --libs`

XMLDCL = /opt/jade/share/xml.dcl

progs: formatter skeleton reflow

formatter: formatter.c 
	$(CC) $(XML_CFLAGS) -o $@ $^ $(XML_LIBS)

skeleton: skeleton.c
	$(CC) $(GRETL_CFLAGS) -o $@ $^ $(GRETL_LIBS)

reflow: dolines.c
	$(CC) -o $@ $^ 

# test generation of XML skeleton from libgretl info
# the skeleton could be used as a template for translation
skeleton.xml: skeleton gretl_commands.dtd
	./skeleton > $@

# test validity of skeleton.xml, and ability of formatter
# to generate dummy 'formatted' output (text and XML)
test: skeleton.xml formatter gretlman.xsl gretltxt.xsl
	xmllint --valid --noout $< && ./formatter skeleton.xml

# a sample: the ols help
sample: formatter reflow test_command.xml gretlman.xsl gretltxt.xsl
	./formatter inprogress.xml # test_command.xml 
	./reflow < commands.txt > commands.hlp

clean:
	rm -f formatter skeleton reflow skeleton.xml *~
