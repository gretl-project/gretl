GUIDE = gretl-guide.pdf
REF = gretl-ref.pdf
ALTGUIDE = gretl-guide-a4.pdf
ALTREF = gretl-ref-a4.pdf

WEBDIR = /home/cottrell/stats/esl/website

XMLFILES = ../chapters/cmdlist.xml \
	../chapters/cmdref.xml \
	../chapters/optarg.xml \
	../chapters/cmdtopics.xml \
	../shell.xml

TEXXSLT = ../gretltex.xsl
MANXSLT = ../commands/gretlman.xsl

$(GUIDE): gretl-guide.tex links textables
	pdflatex $<
	pdflatex $<
	pdflatex $<

$(REF): gretl-ref.tex refbody.tex reftables.tex textables links
	pdflatex $<
	pdflatex $<
	pdflatex $<

$(ALTGUIDE): gretl-guide-a4.tex links
	pdflatex $<
	pdflatex $<
	pdflatex $<

$(ALTREF): gretl-ref-a4.tex links
	pdflatex $<
	pdflatex $<
	pdflatex $<

letter: $(GUIDE) $(REF)

a4: $(ALTGUIDE) $(ALTREF)

gretl-guide-a4.tex: gretl-guide.tex
	cat $< | sed s/letterpaper/a4paper/ > $@

gretl-ref-a4.tex: gretl-ref.tex
	cat $< | sed s/letterpaper/a4paper/ > $@

refbody.tex: $(XMLFILES) $(TEXXSLT)
	xsltproc -o refbody.tex $(TEXXSLT) ../shell.xml

../chapters/cmdlist.xml: ../commands/gretl_commands.xml $(MANXSLT)
	make -C ../commands -f cmdref.mk dbtex

../chapters/cmdtopics.xml: ../commands/gretl_commands.xml $(MANXSLT)
	make -C ../commands -f cmdref.mk dbtex

.PHONY : links textables alt onepass refpass thumbs install clean

links:
	rm -f gretl.sty && ln -s ../gretl.sty ./gretl.sty
	rm -f figures && ln -s ../figures ./figures

textables:
	make -C ../../build/doc/commands textables

alt: gretl-guide-a4.pdf gretl-ref-a4.pdf

onepass:
	pdflatex gretl-guide

refpass: refbody.tex textables
	pdflatex gretl-ref

thumbs: $(DOC)
	thumbpdf $@ && pdflatex $< 	

install: $(GUIDE) $(REF)
	cp $(GUIDE) $(REF) $(WEBDIR)

mildclean:
	rm -f *.aux *.log *.out *.toc *.tpt *~

clean:
	rm -f *.aux *.log *.out *.toc *.tpt $(ALTGUIDE) $(ALTREF) *~
	rm -f refbody.tex gretl-guide-a4.tex gretl-ref-a4.tex
	rm -rf auto

distclean: clean
	rm -rf CVS

