XMLDCL = /opt/jade/share/xml.dcl
MANUAL = manual
DSL = gretl_texmath.dsl
MATHSS = TeXMath.dsl

pdf: $(MANUAL).tex jadetex.cfg lucidabr.tex
	pdfjadetex $(MANUAL)
	pdfjadetex $(MANUAL)
	pdfjadetex $(MANUAL)

$(MANUAL).tex: $(MANUAL).xml $(DSL) $(MATHSS) titlepage.dsl
	jade -t tex -d $(DSL) $(XMLDCL) $(MANUAL).xml
	./unescape_math.pl

html: manual.xml 
	./db2html $(MANUAL).xml


