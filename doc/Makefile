DOC = manual
SRC = $(DOC).xml
DSL = $(DOC)-both.dsl
ALTDSL = gretl-a4.dsl
MATHSS = TeXMath.dsl
HMATHSS = HTMLMath.dsl
XSLDIR = /usr/share/sgml/docbook-xsl-1.48
HTMLDIR = $(DOC)-html

CHAPTERS = chapters/cmdref.xml chapters/intro.xml chapters/modes.xml \
	chapters/panel.xml chapters/appendices.xml chapters/datafiles.xml \
	chapters/looping.xml chapters/starting.xml \
	chapters/cli.xml chapters/graphs.xml chapters/menus.xml \
	chapters/optarg.xml chapters/trouble.xml

pdf: $(DOC).tex jadetex.cfg lucidabr.tex
	-pdfjadetex $(DOC)
	-pdfjadetex $(DOC)
	-pdfjadetex $(DOC)

altpdf: $(DOC)_a4.tex jadetex.cfg lucidabr.tex
	-pdfjadetex $(DOC)_a4
	-pdfjadetex $(DOC)_a4
	-pdfjadetex $(DOC)_a4

thumbs: $(SRC) $(CHAPTERS)
	thumbpdf $(DOC)

$(DOC).tex: $(SRC) $(CHAPTERS) $(DSL) $(MATHSS)
	jade -t tex -d $(DSL)\#print /opt/jade/share/xml.dcl $(SRC)
	./unescape_math.pl $(DOC).tex

$(DOC)_a4.tex: $(SRC) $(CHAPTERS) $(ALTDSL) $(MATHSS)
	jade -t tex -d $(ALTDSL)\#print /opt/jade/share/xml.dcl $(SRC)
	mv $(DOC).tex $(DOC)_a4.tex
	./unescape_math.pl $(DOC)_a4.tex

html: $(SRC) db2html $(DSL) $(HMATHSS) 
	./db2html $(SRC)
	./texmath2png.pl $(DOC)-html/equation-list.sgml

htmlhelp:
	cd gretlhlp && \
	saxon ../$(SRC) ../gretlhlp.xsl tex.math.in.alt=latex
	mkdir -p gretlhlp/figures
	cp figures/*.png gretlhlp/figures
	cd gretlhlp/figures && ./mkgif

clean:
	rm -f *.aux *.log *.out *.pdf $(DOC).tex
	rm -rf $(DOC)-html 
