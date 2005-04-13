include ../config.mk

vpath %.po ../../po

%.mo: %.po
	$(WIN32_MSGFMT) -o $@ $<

all: de.mo es.mo eu.mo fr.mo it.mo ja.mo pl.mo rw.mo tr.mo

clean:
	rm -f *.mo
