include ../config.mk

PODIR = ../../po

all: es.mo fr.mo it.mo

es.mo: $(PODIR)/es.po
	$(WIN32_MSGFMT) $<
	mv messages $@

fr.mo: $(PODIR)/fr.po
	$(WIN32_MSGFMT) $<
	mv messages $@

it.mo: $(PODIR)/it.po
	$(WIN32_MSGFMT) $<
	mv messages $@

clean:
	rm -f es.mo fr.mo it.mo
