CC = gcc -O2 -Wall

# These flags may be have to changed depending on the system
XML_CFLAGS = `pkg-config --cflags libxml-2.0 libxslt`
XML_LIBS = `pkg-config --libs libxml-2.0 libxslt`

dbtex: xsltrans gretl_commands.xml gretlman.xsl
	cat tex.entities > tmp.xml
	cat gretl_commands.xml | grep -v DOCTYPE | grep -v 'xml version' >> tmp.xml
	./xsltrans --docbook tmp.xml && \
	cp cmdlist.xml ../chapters/cmdlist.xml && \
	rm tmp.xml

dbtex_it: xsltrans gretl_commands_it.xml gretlman.xsl
	cat tex.entities > tmp.xml
	cat gretl_commands_it.xml | grep -v DOCTYPE | grep -v 'xml version' >> tmp.xml
	./xsltrans --docbook tmp.xml && \
	echo "<?xml version=\"1.0\" encoding=\"ISO-8859-1\" ?>" > ../chapters_it/cmdlist.xml && \
	cat cmdlist.xml >> ../chapters_it/cmdlist.xml && \
	rm tmp.xml

docbook: xsltrans gretl_commands.xml gretlman.xsl
	./xsltrans --docbook gretl_commands.xml && \
	cp cmdlist.xml ../chapters/cmdlist.xml

docbook_it: xsltrans gretl_commands_it.xml gretlman.xsl
	./xsltrans --docbook gretl_commands_it.xml && \
	cp cmdlist.xml ../chapters_it/cmdlist.xml

xsltrans: xsltrans.c 
	$(CC) $(XML_CFLAGS) -o $@ $^ $(XML_LIBS)

clean: 
	rm -f cmdlist.xml
