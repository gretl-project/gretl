XMLDCL = /opt/jade/share/xml.dcl

pdf: manual.tex jadetex.cfg lucidabr.tex
	pdfjadetex manual
	pdfjadetex manual
	pdfjadetex manual

manual.tex: manual.xml gretl.dsl MathMLrev.dsl titlepage.dsl
	jade -t tex -d ./gretl.dsl $(XMLDCL) manual.xml

html: manual.xml 
	./db2html manual.xml


