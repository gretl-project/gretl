XMLDCL = /opt/jade/share/xml.dcl
MANUAL = manual
DSL = gretl_texmath.dsl
HDSL = gretl_htmlmath.dsl
MATHSS = TeXMath.dsl
HMATHSS = HTMLMath.dsl

pdf: $(MANUAL).tex jadetex.cfg lucidabr.tex
	pdfjadetex $(MANUAL)
	pdfjadetex $(MANUAL)
	pdfjadetex $(MANUAL)

$(MANUAL).tex: $(MANUAL).xml $(DSL) $(MATHSS) titlepage.dsl
	jade -t tex -d $(DSL) $(XMLDCL) $(MANUAL).xml
	./unescape_math.pl

html: manual.xml db2html $(HDSL) $(HMATHSS)
	./db2html $(MANUAL).xml


