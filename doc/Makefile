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
CMDLIST_STANDALONE = chapters/cmdlist_standalone_it.xml
CMDREF_WRAPPER = chapters/cmdref_standalone_it.xml

CMDREF = commands/gretl_commands.xml
CMDREF_IT = commands/gretl_commands_it.xml
MANXSL = commands/gretlman.xsl

%.pdf: %.tex
	pdfjadetex $<
	-grep undefined $(<:.tex=.log) >/dev/null && pdfjadetex $<
	-grep undefined $(<:.tex=.log) >/dev/null && pdfjadetex $<
	-grep undefined $(<:.tex=.log) >/dev/null && pdfjadetex $<

all: manual.pdf manual-a4.pdf html chm # thumbs altthumbs

manual.pdf: manual.tex jadetex.cfg lucidabr.tex

manual-a4.pdf: manual-a4.tex jadetex.cfg lucidabr.tex

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

# Italian standalone PDF command reference

cmdbook_it.pdf: cmdbook_it.tex jadetex.cfg lucidabr.tex

cmdbook_it.tex: cmdbook_it.xml $(CMDREF_WRAPPER) $(ALTDSL) $(MATHSS) $(CMDLIST_STANDALONE)
	$(JADE) -t tex -d $(ALTDSL)\#print $(XMLDCL) $<
	./unescape_math.pl $@

$(CMDLIST): $(CMDREF) $(MANXSL) 
	make -C commands -f cmdref.mk docbook

$(CMDLIST_STANDALONE): $(CMDREF) $(MANXSL) 
	make -C commands -f cmdref.mk docbook-standalone

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
	rm -f *.aux *.log *.out *.pdf manual.tex manual-a4.tex cmdbook.tex

install: 
	cp manual.pdf manual-a4.pdf $(WEBDIR)
