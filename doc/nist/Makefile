DOC = nist
SRC = $(DOC).xml
DSL = $(DOC)-both.dsl
ALTDSL = ../gretl-a4.dsl
MATHSS = ../TeXMath.dsl
HMATHSS = ../HTMLMath.dsl

pdf: $(DOC).tex jadetex.cfg lucidabr.tex
	-pdfjadetex $(DOC)
	-pdfjadetex $(DOC)
	-pdfjadetex $(DOC)

altpdf: $(DOC)_a4.tex jadetex.cfg lucidabr.tex
	-pdfjadetex $(DOC)_a4
	-pdfjadetex $(DOC)_a4
	-pdfjadetex $(DOC)_a4

$(DOC).tex: $(SRC) $(CHAPTERS) $(DSL) $(MATHSS)
	jade -t tex -d $(DSL)\#print /opt/jade/share/xml.dcl $(SRC)
	../unescape_math.pl $(DOC).tex

html: $(SRC) db2html $(DSL) $(HMATHSS) 
	mkdir -p $(DOC)-html
	cd $(DOC)-html && jade -t sgml -ihtml -V nochunks -d ../$(DSL)\#html \
	/opt/jade/share/xml.dcl ../$(SRC) > index.html
	cd $(DOC)-html && ../../texmath2png.pl equation-list.sgml
        
clean:
	rm -f *.aux *.log *.out *.pdf $(DOC).tex
	rm -rf $(DOC)-html 
        
