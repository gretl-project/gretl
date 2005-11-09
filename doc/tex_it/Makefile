DOC = gretl-guide-it.pdf

$(DOC): gretl-guide-it.tex
	pdflatex $<
	pdflatex $<
	pdflatex $<
	thumbpdf $@ && pdflatex $< 

clean:
	rm -f *.aux *.log *.out *.toc *.tpt $(DOC) *~
