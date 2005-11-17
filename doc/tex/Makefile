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

$(GUIDE): links gretl-guide.tex
	pdflatex $<
	pdflatex $<

$(REF): links gretl-ref.tex refbody.tex
	pdflatex $<
	pdflatex $<

$(ALTGUIDE): links gretl-guide-a4.tex
	pdflatex $<
	pdflatex $<

$(ALTREF): links gretl-ref-a4.tex
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

links:
	rm -f gretl.sty && ln -s ../gretl.sty ./gretl.sty
	rm -f figures && ln -s ../figures ./figures

thumbs: $(DOC)
	thumbpdf $@ && pdflatex $< 	

install: $(GUIDE) $(REF)
	cp $(GUIDE) $(REF) $(WEBDIR)

clean:
	rm -f *.aux *.log *.out *.toc *.tpt $(GUIDE) $(REF) *~
