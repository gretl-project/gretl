SVAR.pdf : bib
	pdflatex SVAR; pdflatex SVAR

SVAR.aux : SVAR.tex
	pdflatex SVAR

bib : SVAR.aux
	bibtex $<

clean : 
	rm -rf *~ *.aux *.log *.bbl *.blg *.toc auto/ SVAR.pdf
