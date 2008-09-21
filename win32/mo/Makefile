include ../config.mk

vpath %.po ../../po

%.mo: %.po
	$(WIN32_MSGFMT) -o $@ $<

all: de.mo es.mo eu.mo fr.mo it.mo ja.mo pl.mo ru.mo rw.mo tr.mo pt.mo pt_BR.mo

clean:
	rm -f *.mo
