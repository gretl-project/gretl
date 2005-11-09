DOC = gretl-guide.pdf

$(DOC): gretl-guide.tex
	pdflatex $<
	pdflatex $<
	pdflatex $<
	thumbpdf $@ && pdflatex $< 

clean:
	rm -f *.aux *.log *.out *.toc *.tpt $(DOC) *~
