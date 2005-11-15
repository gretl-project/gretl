GUIDE = gretl-guide-it.pdf
REF = gretl-ref-it.pdf

WEBDIR = /home/cottrell/stats/esl/website

XMLFILES = ../chapters_it/cmdlist.xml \
	../chapters_it/cmdref.xml \
	../chapters_it/optarg.xml \
	../shell_it.xml

XSLT = ../gretltex.xsl

$(GUIDE): gretl-guide-it.tex
	pdflatex $<
	pdflatex $<

$(REF): gretl-ref-it.tex refbody.tex
	pdflatex $<
	pdflatex $<

refbody.tex: $(XMLFILES) $(XSLT) 
	xsltproc -o refbody.tex $(XSLT) ../shell_it.xml

thumbs: $(DOC)
	thumbpdf $@ && pdflatex $< 	

install: $(GUIDE) $(REF)
	cp $(GUIDE) $(REF) $(WEBDIR)

clean:
	rm -f *.aux *.log *.out *.toc *.tpt $(GUIDE) $(REF) *~
