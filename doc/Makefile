DOC = manual
SRC = manual.xml
DSL = gretl-both.dsl
MATHSS = TeXMath.dsl
HMATHSS = HTMLMath.dsl
XSLDIR = /usr/share/sgml/docbook-xsl-1.48

CHAPTERS = chapters/cmdref.xml chapters/intro.xml chapters/modes.xml \
	chapters/panel.xml chapters/appendices.xml chapters/datafiles.xml \
	chapters/looping.xml chapters/nist.xml chapters/starting.xml \
	chapters/cli.xml chapters/graphs.xml chapters/menus.xml \
	chapters/optarg.xml chapters/trouble.xml

pdf: $(DOC).tex jadetex.cfg lucidabr.tex
	-pdfjadetex $(DOC)
	-pdfjadetex $(DOC)
	-pdfjadetex $(DOC)

thumbs: $(SRC) $(CHAPTERS)
	thumbpdf $(DOC)

$(DOC).tex: $(SRC) $(CHAPTERS) $(DSL) $(MATHSS)
	jade -t tex -d $(DSL)\#print /opt/jade/share/xml.dcl $(SRC)
	./unescape_math.pl $(DOC).tex

html: $(SRC) db2html $(DSL) $(HMATHSS) 
	./db2html $(SRC)
	./texmath2png.pl manual-html/equation-list.sgml

htmlhelp:
	cd htmlhelp && \
	saxon ../$(SRC) $(XSLDIR)/htmlhelp/htmlhelp.xsl tex.math.in.alt=latex

clean:
	rm -f *.aux *.log *.out *.pdf $(DOC).tex
	rm -rf manual-html manual-html.junk