include ../config.mk

LANGS = $(shell cat ../../po/LINGUAS | grep -v \\\#)
MO = $(LANGS:=.mo)

vpath %.po ../../po

%.mo: %.po
	$(WIN32_MSGFMT) -o $@ $<

all: $(MO)

clean:
	rm -f *.mo
