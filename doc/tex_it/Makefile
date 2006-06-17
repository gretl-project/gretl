GUIDE = gretl-guide-it.pdf
REF = gretl-ref-it.pdf

WEBDIR = /home/cottrell/stats/esl/website

XMLFILES = ../chapters_it/cmdlist.xml \
	../chapters_it/cmdref.xml \
	../chapters_it/optarg.xml \
	../shell_it.xml

TEXXSLT = ../gretltex.xsl
MANXSLT = ../commands/gretlman.xsl

$(GUIDE): gretl-guide-it.tex links textables
	pdflatex $<
	pdflatex $<
	pdflatex $<

$(REF): gretl-ref-it.tex refbody.tex reftables.tex textables links
	pdflatex $<
	pdflatex $<
	pdflatex $<

refbody.tex: $(XMLFILES) $(TEXXSLT) 
	xsltproc -o refbody.tex $(TEXXSLT) ../shell_it.xml

../chapters_it/cmdlist.xml: ../commands/gretl_commands_it.xml $(MANXSLT)
	make -C ../commands -f cmdref.mk dbtex_it

.PHONY : textables links onepass refpass thumbs install clean

textables:
	make -C ../commands textables

links:
	rm -f gretl.sty && ln -s ../gretl.sty ./gretl.sty
	rm -f figures && ln -s ../figures_it ./figures

onepass:
	pdflatex gretl-guide-it

refpass: refbody.tex
	pdflatex gretl-ref-it

thumbs: $(DOC)
	thumbpdf $@ && pdflatex $< 	

install: $(GUIDE) $(REF)
	cp $(GUIDE) $(REF) $(WEBDIR)

clean:
	rm -f *.aux *.log *.out *.toc *.tpt $(GUIDE) $(REF) *~
	rm -f refbody.tex

distclean: clean
	rm -rf CVS
