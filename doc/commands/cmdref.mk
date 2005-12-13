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

dbtex_es: xsltrans gretl_commands_es.xml gretlman.xsl
	cat tex.entities > tmp.xml
	cat gretl_commands_es.xml | grep -v DOCTYPE | grep -v 'xml version' >> tmp.xml
	./xsltrans --docbook tmp.xml && \
	echo "<?xml version=\"1.0\" encoding=\"ISO-8859-1\" ?>" > ../chapters_es/cmdlist.xml && \
	cat cmdlist.xml >> ../chapters_es/cmdlist.xml && \
	rm tmp.xml

xsltrans: xsltrans.c 
	$(CC) $(XML_CFLAGS) -o $@ $^ $(XML_LIBS)

reflow: reflow.c 
	$(CC) -o $@ $^

olddoc: xsltrans reflow gretl_commands.xml gretl_commands_it.xml gretl_commands_es.xml gretltxt.xsl
	mkdir -p ../../share/latin1
	./xsltrans --hlp gretl_commands.xml && \
	./reflow < guilist.txt > ../../share/latin1/gretlgui.hlp
	./xsltrans --hlp gretl_commands_it.xml && \
	./reflow < guilist.txt > ../../share/latin1/gretlgui.hlp.it
	./xsltrans --hlp gretl_commands_es.xml && \
	./reflow < guilist.txt > ../../share/latin1/gretlgui.hlp.es

clean: 
	rm -f cmdlist.xml *list.txt xsltrans reflow
