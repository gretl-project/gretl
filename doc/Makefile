DOC = manual
DSL = gretl-both.dsl
MATHSS = TeXMath.dsl
HMATHSS = HTMLMath.dsl

CHAPTERS = chapters/cmdref.sgml chapters/intro.sgml chapters/modes.sgml \
	chapters/panel.sgml chapters/appendices.sgml chapters/datafiles.sgml \
	chapters/looping.sgml chapters/nist.sgml chapters/starting.sgml \
	chapters/cli.sgml chapters/graphs.sgml chapters/menus.sgml \
	chapters/optarg.sgml chapters/trouble.sgml

pdf: $(DOC).tex jadetex.cfg lucidabr.tex
	pdfjadetex $(DOC)
	pdfjadetex $(DOC)
	pdfjadetex $(DOC)

$(DOC).tex: $(DOC).sgml $(CHAPTERS) $(DSL) $(MATHSS) titlepage.dsl unescape_math.pl
	jade -t tex -d $(DSL)\#print $(DOC).sgml
	./unescape_math.pl $(DOC).tex

html: $(DOC).sgml db2html $(DSL) $(HMATHSS) 
	./db2html $(DOC).sgml
	./texmath2png.pl manual/equation-list.sgml

clean:
	rm -rf $(DOC) $(DOC).junk $(DOC).tex *.log *.out *.aux
