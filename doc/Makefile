# location of various things

XSLDIR = /usr/share/sgml/docbook-xsl
XMLDCL = /opt/jade/share/xml.dcl
JADE = jade
SAXON = saxon
CONVERT = convert # from ImageMagick suite
HHC = wine "c:/Program Files/HTML Help Workshop/hhc.exe"
WEBDIR = /home/cottrell/stats/esl/website

DSL = manual-both.dsl
ALTDSL = manual-a4.dsl
MATHSS = TeXMath.dsl
HMATHSS = HTMLMath.dsl

CHAPTERS = chapters/cmdref.xml chapters/intro.xml chapters/modes.xml \
	chapters/panel.xml chapters/appendices.xml chapters/datafiles.xml \
	chapters/nls.xml chapters/looping.xml chapters/starting.xml \
	chapters/cli.xml chapters/graphs.xml chapters/optarg.xml \
	chapters/trouble.xml chapters/biblio.xml 

CMDLIST = chapters/cmdlist.xml

CMDREF = commands/gretl_commands.xml
MANXSL = commands/gretlman.xsl

all: manual.pdf manual-a4.pdf html chm # thumbs altthumbs

manual.pdf: manual.tex jadetex.cfg lucidabr.tex
	pdfjadetex manual
	-grep undefined manual.log >/dev/null && pdfjadetex manual
	-grep undefined manual.log >/dev/null && pdfjadetex manual
	-grep undefined manual.log >/dev/null && pdfjadetex manual

manual-a4.pdf: manual-a4.tex jadetex.cfg lucidabr.tex
	pdfjadetex manual-a4
	-grep undefined manual-a4.log >/dev/null && pdfjadetex manual-a4
	-grep undefined manual-a4.log >/dev/null && pdfjadetex manual-a4
	-grep undefined manual-a4.log >/dev/null && pdfjadetex manual-a4

thumbs: manual.pdf
	thumbpdf manual
	pdfjadetex manual
        
altthumbs: manual-a4.pdf
	thumbpdf manual-a4
	pdfjadetex manual-a4

manual.tex: manual.xml $(CHAPTERS) $(DSL) $(MATHSS) $(CMDLIST)
	$(JADE) -t tex -d $(DSL)\#print $(XMLDCL) $<
	./unescape_math.pl $@

manual-a4.tex: manual-a4.xml $(CHAPTERS) $(ALTDSL) $(MATHSS) $(CMDLIST)
	$(JADE) -t tex -d $(ALTDSL)\#print $(XMLDCL) $<
	./unescape_math.pl $@

manual-a4.xml: manual.xml
	cp $< $@

$(CMDLIST): $(CMDREF) $(MANXSL) 
	make -C commands docbook

html: manual.xml db2html $(DSL) $(HMATHSS) 
	./db2html manual.xml
	./texmath2png.pl manual.html/equation-list.sgml

# make windows HTML help file for gretl

topdoc = .

chm: $(topdoc)/manual.xml $(topdoc)/chapters/*.xml 
	mkdir -p manual.chm/figures
	cp $(topdoc)/figures/*.png manual.chm/figures
	for f in manual.chm/figures/*.png ; do $(CONVERT) $$f $${f%.png}.gif ; done
	cd manual.chm && \
	$(SAXON) ../$(topdoc)/manual.xml ../$(topdoc)/gretlhlp.xsl \
	tex.math.in.alt=latex && $(HHC) gretl.hhp

clean:
	rm -f *.aux *.log *.out *.pdf manual.tex manual-a4.tex
	rm -rf manual.html manual.chm

texclean:
	rm -f *.aux *.log *.out *.pdf manual.tex manual-a4.tex

install: 
	cp manual.pdf manual-a4.pdf $(WEBDIR)
