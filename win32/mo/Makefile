include ../config.mk

PODIR = ../../po

es.mo: $(PODIR)/es.po
	$(WIN32_MSGFMT) $<
	mv messages $@

fr.mo: $(PODIR)/fr.po
	$(WIN32_MSGFMT) $<
	mv messages $@

clean:
	rm -f es.mo fr.mo
