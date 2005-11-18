GUIDE = gretl-guide.pdf
REF = gretl-ref.pdf
ALTGUIDE = gretl-guide-a4.pdf
ALTREF = gretl-ref-a4.pdf

WEBDIR = /home/cottrell/stats/esl/website

XMLFILES = ../chapters/cmdlist.xml \
	../chapters/cmdref.xml \
	../chapters/optarg.xml \
	../shell.xml

TEXXSLT = ../gretltex.xsl
MANXSLT = ../commands/gretlman.xsl

$(GUIDE): gretl-guide.tex links
	pdflatex $<
	pdflatex $<
	pdflatex $<

$(REF): gretl-ref.tex refbody.tex links
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

gretl-guide-a4.tex: gretl-guide.tex
	cat $< | sed s/letterpaper/a4paper/ > $@

gretl-ref-a4.tex: gretl-ref.tex
	cat $< | sed s/letterpaper/a4paper/ > $@

refbody.tex: $(XMLFILES) $(TEXXSLT)
	xsltproc -o refbody.tex $(TEXXSLT) ../shell.xml

../chapters/cmdlist.xml: ../commands/gretl_commands.xml $(MANXSLT)
	make -C ../commands -f cmdref.mk dbtex

.PHONY : links thumbs install clean

links:
	rm -f gretl.sty && ln -s ../gretl.sty ./gretl.sty
	rm -f figures && ln -s ../figures ./figures

thumbs: $(DOC)
	thumbpdf $@ && pdflatex $< 	

install: $(GUIDE) $(REF)
	cp $(GUIDE) $(REF) $(WEBDIR)

clean:
	rm -f *.aux *.log *.out *.toc *.tpt $(GUIDE) $(REF) *~
