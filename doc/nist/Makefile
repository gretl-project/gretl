DOC = nist
SRC = $(DOC).xml
DSL = $(DOC)-both.dsl
XMLDCL = /opt/jade/share/xml.dcl
HMATHSS = ../HTMLMath.dsl

all: nist.html

nist.html: $(SRC) $(DSL) $(HMATHSS) figures
	jade -t sgml -i html -V nochunks -d $(DSL)\#html $(XMLDCL) $< > $@
	../texmath2png.pl equation-list.sgml

.PHONY:

figures:
	mkdir -p figures

clean:
	rm -f *.html equation-list.sgml
	rm -rf figures
