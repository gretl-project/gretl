DOC = manual
DSL = gretl-both.dsl
MATHSS = TeXMath.dsl
HMATHSS = HTMLMath.dsl

pdf: $(DOC).tex jadetex.cfg lucidabr.tex
	pdfjadetex $(DOC)
	pdfjadetex $(DOC)
	pdfjadetex $(DOC)

$(DOC).tex: $(DOC).sgml $(DSL) $(MATHSS) titlepage.dsl unescape_math.pl
	jade -t tex -d $(DSL)\#print $(DOC).sgml
	./unescape_math.pl $(DOC).tex

html: $(DOC).sgml db2html $(DSL) $(HMATHSS) pix
	./db2html $(DOC).sgml

clean:
	rm -rf $(DOC) $(DOC).junk $(DOC).tex *.log *.out *.aux

pix: $(DOC).sgml grab_equations.pl
	spam -x -x $(DOC).sgml > tmp.sgml && \
	./grab_equations.pl tmp.sgml && rm tmp.sgml


