GUIDE = gretl-guide.pdf
REF = gretl-ref.pdf
ALTGUIDE = gretl-guide-a4.pdf
ALTREF = gretl-ref-a4.pdf

WEBDIR = /home/cottrell/stats/esl/website

XMLFILES = ../chapters/cmdlist.xml \
	../chapters/cmdref.xml \
	../chapters/optarg.xml \
	../shell.xml

XSLT = ../gretltex.xsl

$(GUIDE): gretl-guide.tex
	pdflatex $<
	pdflatex $<

$(REF): gretl-ref.tex refbody.tex
	pdflatex $<
	pdflatex $<

$(ALTGUIDE): gretl-guide-a4.tex
	pdflatex $<
	pdflatex $<

$(ALTREF): gretl-ref-a4.tex
	pdflatex $<
	pdflatex $<

gretl-guide-a4.tex: gretl-guide.tex
	cat $< | sed s/letterpaper/a4paper/ > $@

gretl-ref-a4.tex: gretl-ref.tex
	cat $< | sed s/letterpaper/a4paper/ > $@

refbody.tex: $(XMLFILES) $(XSLT)
	xsltproc -o refbody.tex $(XSLT) ../shell.xml

thumbs: $(DOC)
	thumbpdf $@ && pdflatex $< 	

install: $(GUIDE) $(REF)
	cp $(GUIDE) $(REF) $(WEBDIR)

clean:
	rm -f *.aux *.log *.out *.toc *.tpt $(GUIDE) $(REF) *~
