include ../config.mk

PODIR = ../../po

es.mo: $(PODIR)/es.po
	cp $< .
	msgconv -t utf-8 -o es-utf8.po es.po
	$(WIN32_MSGFMT) es-utf8.po && \
	rm es-utf8.po
	mv messages $@ && rm es.po

clean:
	rm -f es.po es-utf8.po es.mo
