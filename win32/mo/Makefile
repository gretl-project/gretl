include ../config.mk

PODIR = ../../po

all: es.mo fr.mo it.mo ja.mo

es.mo: $(PODIR)/es.po
	$(WIN32_MSGFMT) $<
	mv messages.mo $@

fr.mo: $(PODIR)/fr.po
	$(WIN32_MSGFMT) $<
	mv messages.mo $@

it.mo: $(PODIR)/it.po
	$(WIN32_MSGFMT) $<
	mv messages.mo $@

ja.mo: $(PODIR)/ja.po
	$(WIN32_MSGFMT) $<
	mv messages.mo $@

clean:
	rm -f es.mo fr.mo it.mo ja.mo
