include ../config.mk

vpath %.po ../../po

%.mo: %.po
	$(WIN32_MSGFMT) $<
	mv messages.mo $@

all: es.mo fr.mo it.mo ja.mo pl.mo

clean:
	rm -f *.mo
