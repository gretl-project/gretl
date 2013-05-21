include ../config.mk

LANGS = $(shell cat ../../po/LINGUAS | grep -v \\\#)

all:
	for f in $(LANGS) ; do $(WIN32_MSGFMT) -o $$f.mo ../../po/$$f.po ; done

clean:
	rm -f *.mo
