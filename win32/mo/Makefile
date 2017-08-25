include ../config.mk

LANGS = $(shell cat $(podir)/LINGUAS | grep -v \\\#)
MO = $(LANGS:=.mo)

vpath %.po $(podir)

%.mo: %.po
	$(WIN32_MSGFMT) -o $@ $<

all: $(MO)

clean:
	rm -f *.mo
