# location of various things

XSLDIR = /usr/share/sgml/docbook-xsl
XMLDCL = /opt/jade/share/xml.dcl
JADE = jade
XSLTPROC = xsltproc
CONVERT = convert # from ImageMagick suite
HHC = wine "c:/Program Files/HTML Help Workshop/hhc.exe"
WEBDIR = /home/cottrell/stats/esl/website

DSL = manual-both.dsl
ALTDSL = manual-a4.dsl
MATHSS = TeXMath.dsl
HMATHSS = HTMLMath.dsl
TITLESS = titlepage.dsl

CHAPTERS = chapters/cmdref.xml chapters/intro.xml chapters/modes.xml \
	chapters/panel.xml chapters/appendices.xml chapters/datafiles.xml \
	chapters/nls.xml chapters/looping.xml chapters/starting.xml \
	chapters/cli.xml chapters/graphs.xml chapters/optarg.xml \
	chapters/trouble.xml chapters/biblio.xml chapters/sampling.xml \
	chapters/functions.xml chapters/genr.xml chapters/vecm.xml

CHAPTERS_IT = chapters_it/cmdref.xml chapters_it/intro.xml chapters_it/modes.xml \
	chapters_it/panel.xml chapters_it/appendices.xml chapters_it/datafiles.xml \
	chapters_it/nls.xml chapters_it/looping.xml chapters_it/starting.xml \
	chapters_it/cli.xml chapters_it/graphs.xml chapters_it/optarg.xml \
	chapters_it/trouble.xml chapters_it/biblio.xml chapters_it/sampling.xml \
        chapters_it/functions.xml chapters_it/genr.xml chapters_it/vecm.xml

CMDLIST = chapters/cmdlist.xml
CMDLIST_IT = chapters_it/cmdlist.xml

CMDREF = commands/gretl_commands.xml
CMDREF_IT = commands/gretl_commands_it.xml

MANXSL = commands/gretlman.xsl

%.pdf: %.tex
	pdfjadetex $<
	-grep undefined $(<:.tex=.log) >/dev/null && pdfjadetex $<
	-grep undefined $(<:.tex=.log) >/dev/null && pdfjadetex $<
	-grep undefined $(<:.tex=.log) >/dev/null && pdfjadetex $<
	-grep 'Rerun LaTeX' $(<:.tex=.log) >/dev/null && pdfjadetex $<

all: thumbs altthumbs thumbs_it html chm html_it chm_it

# thumbs altthumbs

manual.pdf: manual.tex

manual-a4.pdf: manual-a4.tex

manual_it.pdf: manual_it.tex

thumbs: manual.pdf
	thumbpdf manual
	pdfjadetex manual

altthumbs: manual-a4.pdf
	thumbpdf manual-a4
	pdfjadetex manual-a4

thumbs_it: manual_it.pdf
	thumbpdf manual_it
	pdfjadetex manual_it

manual.tex: manual.xml $(CHAPTERS) $(CMDLIST) mathextra.inc $(DSL) $(MATHSS) $(TITLESS) \
	$(CMDLIST) jadetex.cfg.en lucidabr.tex
	rm -f jadetex.cfg && ln -s jadetex.cfg.en jadetex.cfg
	$(JADE) -t tex -d $(DSL)\#print $(XMLDCL) $<
	./unescape_math.pl $@

manual-a4.tex: manual-a4.xml $(CHAPTERS) mathextra.inc $(ALTDSL) $(MATHSS) \
	$(TITLESS) $(CMDLIST) jadetex.cfg.en lucidabr.tex
	rm -f jadetex.cfg && ln -s jadetex.cfg.en jadetex.cfg
	$(JADE) -t tex -d $(ALTDSL)\#print $(XMLDCL) $<
	./unescape_math.pl $@

manual_it.tex: manual_it.xml $(CHAPTERS_IT) $(CMDLIST_IT) mathextra.inc $(ALTDSL) $(MATHSS) \
	$(TITLESS) $(CMDLIST_IT) jadetex.cfg.it lucidabr.tex
	rm -f jadetex.cfg && ln -s jadetex.cfg.it jadetex.cfg
	$(JADE) -t tex -d $(ALTDSL)\#print $(XMLDCL) $<
	./unescape_math.pl $@

manual-a4.xml: manual.xml
	cp $< $@

# command reference section

$(CMDLIST): $(CMDREF) $(MANXSL) 
	make -C commands -f cmdref.mk docbook

$(CMDLIST_IT): $(CMDREF) $(MANXSL) 
	make -C commands -f cmdref.mk docbook_it

html: manual.xml db2html $(DSL) $(HMATHSS) $(CMDLIST)
	./db2html manual.xml
	./texmath2png.pl manual.html/equation-list.sgml

html_it: manual_it.xml db2html $(DSL) $(HMATHSS) $(CMDLIST_IT)
	./db2html manual_it.xml
	./texmath2png.pl manual_it.html/equation-list.sgml

# make windows HTML help file for gretl

topdoc = .

chmfigs: $(topdoc)/figures/*.png
	mkdir -p manual.chm/figures
	cp $(topdoc)/figures/*.png manual.chm/figures
	for f in manual.chm/figures/*.png ; do $(CONVERT) $$f $${f%.png}.gif ; done
	rm -f manual.chm/figures/*.png

chm: $(topdoc)/manual.xml $(topdoc)/chapters/*.xml $(CMDLIST) chmfigs
	cd manual.chm && \
	$(XSLTPROC) --param tex.math.in.alt latex \
	../$(topdoc)/gretlhlp.xsl ../$(topdoc)/manual.xml && \
	../fixtoc.sh && $(HHC) gretl.hhp

chmfigs_it: $(topdoc)/figures/*.png $(topdoc)/figures_it/*.png
	mkdir -p manual_it.chm/figures
	cp $(topdoc)/figures/*.png manual_it.chm/figures
	cp $(topdoc)/figures_it/*.png manual_it.chm/figures
	for f in manual_it.chm/figures/*.png ; do $(CONVERT) $$f $${f%.png}.gif ; done
	rm -f manual_it.chm/figures/*.png

chm_it: $(topdoc)/manual_it.xml $(topdoc)/chapters_it/*.xml $(CMDLIST_IT) chmfigs_it
	cd manual_it.chm && \
	saxon ../$(topdoc)/manual_it.xml ../$(topdoc)/gretlhlp.xsl \
	"tex.math.in.alt=latex" "lang=it" "htmlhelp.hhp=gretl_it.hhp" \
	"htmlhelp.chm=gretl_it.chm" && \
	../fixfigs.sh && ../fixtoc.sh && $(HHC) gretl_it.hhp

clean:
	rm -f *.aux *.log *.out *.pdf manual.tex manual-a4.tex manual_it.tex
	rm -rf manual.html manual.chm manual_it.html manual_it.chm

texclean:
	rm -f *.aux *.log *.out *.pdf manual.tex manual-a4.tex manual_it.tex

install: 
	cp manual.pdf manual-a4.pdf manual_it.pdf $(WEBDIR)
