DOC = manual
SRC = $(DOC).xml
DSL = $(DOC)-both.dsl
ALTDOC = $(DOC)-a4
ALTSRC = $(ALTDOC).xml
ALTDSL = manual-a4.dsl
MATHSS = TeXMath.dsl
HMATHSS = HTMLMath.dsl
XSLDIR = /usr/share/sgml/docbook-xsl-1.48
HTMLDIR = $(DOC)-html

CHAPTERS = chapters/cmdref.xml chapters/intro.xml chapters/modes.xml \
	chapters/panel.xml chapters/appendices.xml chapters/datafiles.xml \
	chapters/looping.xml chapters/starting.xml \
	chapters/cli.xml chapters/graphs.xml chapters/menus.xml \
	chapters/optarg.xml chapters/trouble.xml

all: thumbs pdf altthumbs altpdf html htmlhelp

pdf: $(DOC).tex jadetex.cfg lucidabr.tex
	-pdfjadetex $(DOC)
	-pdfjadetex $(DOC)
	-pdfjadetex $(DOC)

altpdf: $(ALTDOC).tex jadetex.cfg lucidabr.tex
	echo "Running pdfjadetex on $(ALTDOC)"
	-pdfjadetex $(ALTDOC)
	-pdfjadetex $(ALTDOC)
	-pdfjadetex $(ALTDOC)

thumbs: $(SRC) $(CHAPTERS)
	thumbpdf $(DOC)
        
altthumbs: $(ALTSRC) $(CHAPTERS)
	thumbpdf $(ALTDOC)

$(DOC).tex: $(SRC) $(CHAPTERS) $(DSL) $(MATHSS)
	jade -t tex -d $(DSL)\#print /opt/jade/share/xml.dcl $(SRC)
	./unescape_math.pl $(DOC).tex

$(ALTDOC).tex: $(ALTSRC) $(CHAPTERS) $(ALTDSL) $(MATHSS)
	echo "Running jade on $(ALTSRC)"
	jade -t tex -d $(ALTDSL)\#print /opt/jade/share/xml.dcl $(ALTSRC)
	echo "Unescaping math in $(ALTDOC).tex"
	./unescape_math.pl $(ALTDOC).tex

$(ALTSRC): $(SRC)
	echo "Copying $(SRC) to $(ALTSRC)"
	cp $(SRC) $(ALTSRC)       

html: $(SRC) db2html $(DSL) $(HMATHSS) 
	./db2html $(SRC)
	./texmath2png.pl $(DOC)-html/equation-list.sgml

htmlhelp:
	mkdir -p gretlhlp
	cd gretlhlp && \
	saxon ../$(SRC) ../gretlhlp.xsl tex.math.in.alt=latex
	mkdir -p gretlhlp/figures
	cp figures/*.png gretlhlp/figures
	for f in gretlhlp/figures/*.png ; do
	convert $f ${f%.png}.gif ; done
	zip -r gretlhlp.zip gretlhlp/*

clean:
	rm -f *.aux *.log *.out *.pdf $(DOC).tex $(ALTDOC).tex
	rm -rf $(DOC)-html 
