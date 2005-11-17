GUIDE = gretl-guide-es.pdf
REF = gretl-ref-es.pdf

WEBDIR = /home/cottrell/stats/esl/website

XMLFILES = ../chapters_es/cmdlist.xml \
	../chapters_es/cmdref.xml \
	../chapters_es/optarg.xml \
	../shell_es.xml

TEXXSLT = ../gretltex.xsl
MANXSLT = ../commands/gretlman.xsl

$(GUIDE): links gretl-guide-es.tex
	pdflatex $<
	pdflatex $<

$(REF): links gretl-ref-es.tex refbody.tex
	pdflatex $<
	pdflatex $<

refbody.tex: $(XMLFILES) $(TEXXSLT) 
	xsltproc -o refbody.tex $(TEXXSLT) ../shell_es.xml

../chapters_es/cmdlist.xml: ../commands/gretl_commands_es.xml $(MANXSLT)
	make -C ../commands -f cmdref.mk dbtex_es

links:
	rm -f gretl.sty && ln -s ../gretl.sty ./gretl.sty
	rm -f figures && ln -s ../figures_es ./figures

thumbs: $(DOC)
	thumbpdf $@ && pdflatex $< 	

install: $(GUIDE) $(REF)
	cp $(GUIDE) $(REF) $(WEBDIR)

clean:
	rm -f *.aux *.log *.out *.toc *.tpt $(GUIDE) $(REF) *~
