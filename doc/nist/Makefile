DOC = nist
SRC = $(DOC).xml
DSL = $(DOC)-both.dsl
HMATHSS = ../HTMLMath.dsl
HTMLDIR = $(DOC)-html

all: html

html: $(SRC) ../db2html $(DSL) $(HMATHSS) 
	../db2html $(SRC)
	../texmath2png.pl $(DOC)-html/equation-list.sgml

clean:
	rm -f *.html
