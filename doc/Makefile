MANUAL = manual
DSL = gretl-both.dsl
MATHSS = TeXMath.dsl
HMATHSS = HTMLMath.dsl

pdf: $(MANUAL).tex jadetex.cfg lucidabr.tex
	pdfjadetex $(MANUAL)
	pdfjadetex $(MANUAL)
	pdfjadetex $(MANUAL)

$(MANUAL).tex: $(MANUAL).sgml $(DSL) $(MATHSS) titlepage.dsl
	jade -t tex -d $(DSL)\#print $(MANUAL).sgml
	./unescape_math.pl

html: manual.sgml db2html $(HDSL) $(HMATHSS) pix
	./db2html $(MANUAL).sgml

clean:
	rm -rf manual manual.junk manual.tex *.log *.out *.aux

pix: manual.sgml
	./grab_equations.pl


