CC = gcc -O2 -Wall

# These flags may be have to changed depending on the system
XML_CFLAGS = `pkg-config --cflags libxml-2.0 libxslt`
XML_LIBS = `pkg-config --libs libxml-2.0 libxslt`

# and these too
GRETL_CFLAGS = `pkg-config --cflags gretl`
GRETL_LIBS = `pkg-config --libs gretl`

progs: xsltrans skeleton reflow validate

xsltrans: xsltrans.c 
	$(CC) $(XML_CFLAGS) -o $@ $^ $(XML_LIBS)

skeleton: skeleton.c
	$(CC) $(GRETL_CFLAGS) -o $@ $^ $(GRETL_LIBS)

validate: validate.c
	$(CC) $(XML_CFLAGS) $(GRETL_CFLAGS) -o $@ $^ $(GRETL_LIBS) $(XML_LIBS)

reflow: reflow.c
	$(CC) -o $@ $^ 

# test generation of XML skeleton from libgretl info:
# the skeleton could be used as a template for translation
skeleton.xml: skeleton gretl_commands.dtd
	./skeleton > $@

# test validity of skeleton.xml, and ability of xsltrans
# to generate formatted output (text and XML)
test: skeleton.xml xsltrans gretlman.xsl gretltxt.xsl
	xmllint --valid --noout $< && ./xsltrans skeleton.xml

# Generate docbook XML for the command-reference section of the 
# full manual, cmdlist.xml
docbook: xsltrans gretl_commands.xml gretlman.xsl gretltxt.xsl
	./xsltrans gretl_commands.xml && \
	cp cmdlist.xml ../chapters/cmdlist.xml

# Generate "online" help file, gretlcli.hlp and gretl.hlp
hlp: xsltrans reflow gretl_commands.xml gretlman.xsl gretltxt.xsl
	./xsltrans gretl_commands.xml && \
	./reflow < cmdlist.txt > ../../share/gretlcli.hlp
	./reflow < guilist.txt > ../../share/gretl.hlp

docclean:
	 rm -f skeleton.xml cmdlist.*

clean: docclean
	rm -f xsltrans skeleton reflow validate
