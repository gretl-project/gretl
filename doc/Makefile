# location of various things

XSLDIR = /usr/share/sgml/docbook-xsl
XMLDCL = /opt/jade/share/xml.dcl
JADE = jade
SAXON = saxon
CONVERT = convert # from ImageMagick suite
HHC = wine "c:/Program Files/HTML Help Workshop/hhc.exe"
WEBDIR = /home/cottrell/stats/esl/website

DOC = manual
SRC = $(DOC).xml
DSL = $(DOC)-both.dsl
ALTDOC = $(DOC)-a4
ALTSRC = $(ALTDOC).xml
ALTDSL = manual-a4.dsl
TESTDOC = $(DOC)-times
TESTSRC = $(TESTDOC).xml
TESTDSL = manual-times.dsl
MATHSS = TeXMath.dsl
HMATHSS = HTMLMath.dsl
HTMLDIR = $(DOC)-html

CHAPTERS = chapters/cmdref.xml chapters/intro.xml chapters/modes.xml \
	chapters/panel.xml chapters/appendices.xml chapters/datafiles.xml \
	chapters/nls.xml chapters/looping.xml chapters/starting.xml \
	chapters/cli.xml chapters/graphs.xml chapters/optarg.xml \
	chapters/trouble.xml chapters/biblio.xml

all: $(DOC).pdf $(ALTDOC).pdf html chm # thumbs altthumbs

$(DOC).pdf: $(DOC).tex jadetex.cfg lucidabr.tex
	pdfjadetex $(DOC)
	-grep undefined $(DOC).log >/dev/null && pdfjadetex $(DOC)
	-grep undefined $(DOC).log >/dev/null && pdfjadetex $(DOC)

$(ALTDOC).pdf: $(ALTDOC).tex jadetex.cfg lucidabr.tex
	pdfjadetex $(ALTDOC)
	-grep undefined $(ALTDOC).log >/dev/null && pdfjadetex $(ALTDOC)
	-grep undefined $(ALTDOC).log >/dev/null && pdfjadetex $(ALTDOC)

pdftest: $(TESTDOC).tex jadetex.cfg 
	-pdfjadetex $(TESTDOC)
	-pdfjadetex $(TESTDOC)
	-pdfjadetex $(TESTDOC)

thumbs: $(SRC) $(CHAPTERS)
	thumbpdf $(DOC)
	pdfjadetex $(DOC)
        
altthumbs: $(ALTSRC) $(CHAPTERS)
	thumbpdf $(ALTDOC)
	pdfjadetex $(ALTDOC)

$(DOC).tex: $(SRC) $(CHAPTERS) $(DSL) $(MATHSS)
	$(JADE) -t tex -d $(DSL)\#print $(XMLDCL) $(SRC)
	./unescape_math.pl $(DOC).tex

$(ALTDOC).tex: $(ALTSRC) $(CHAPTERS) $(ALTDSL) $(MATHSS)
	echo "Running jade on $(ALTSRC)"
	$(JADE) -t tex -d $(ALTDSL)\#print $(XMLDCL) $(ALTSRC)
	echo "Unescaping math in $(ALTDOC).tex"
	./unescape_math.pl $(ALTDOC).tex

$(TESTDOC).tex: $(TESTSRC) $(CHAPTERS) $(TESTDSL) $(MATHSS)
	$(JADE) -t tex -d $(TESTDSL)\#print $(XMLDCL) $(TESTSRC)
	./unescape_math.pl $(TESTDOC).tex

$(ALTSRC): $(SRC)
	echo "Copying $(SRC) to $(ALTSRC)"
	cp $(SRC) $(ALTSRC) 

$(TESTSRC): $(SRC)
	echo "Copying $(SRC) to $(TESTSRC)"
	cp $(SRC) $(TESTSRC) 

html: $(SRC) db2html $(DSL) $(HMATHSS) 
	./db2html $(SRC)
	./texmath2png.pl $(DOC)-html/equation-list.sgml

# make windows HTML help file for gretl

topdoc = .

chm: $(topdoc)/manual.xml $(topdoc)/chapters/*.xml 
	mkdir -p htmlhelp/figures
	cp $(topdoc)/figures/*.png htmlhelp/figures
	for f in htmlhelp/figures/*.png ; do $(CONVERT) $$f $${f%.png}.gif ; done
	cd htmlhelp && \
	$(SAXON) ../$(topdoc)/manual.xml ../$(topdoc)/gretlhlp.xsl tex.math.in.alt=latex
	$(HHC) htmlhelp/gretl.hhp

clean:
	rm -f *.aux *.log *.out *.pdf $(DOC).tex $(ALTDOC).tex $(TESTDOC).tex
	rm -rf $(DOC)-html htmlhelp

texclean:
	rm -f *.aux *.log *.out *.pdf $(DOC).tex $(ALTDOC).tex $(TESTDOC).tex

install: 
	cp manual.pdf manual-a4.pdf $(WEBDIR)
