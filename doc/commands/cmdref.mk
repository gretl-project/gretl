CC = gcc -O2 -Wall

# These flags may be have to changed depending on the system
XML_CFLAGS = `pkg-config --cflags libxml-2.0 libxslt`
XML_LIBS = `pkg-config --libs libxml-2.0 libxslt`

docbook: xsltrans gretl_commands.xml gretlman.xsl
	./xsltrans --docbook gretl_commands.xml && \
	cp cmdlist.xml ../chapters/cmdlist.xml

docbook-standalone: xsltrans gretl_commands_it.xml gretlman.xsl
	./xsltrans --docbook-standalone gretl_commands_it.xml && \
	cp cmdlist_standalone.xml ../chapters/cmdlist_standalone_it.xml

xsltrans: xsltrans.c 
	$(CC) $(XML_CFLAGS) -o $@ $^ $(XML_LIBS)

clean: 
	rm -f cmdlist.xml
