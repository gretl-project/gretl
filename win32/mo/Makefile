PODIR = ../../po

es.mo: $(PODIR)/es.po
	cp $(PODIR)/es.po .
	msgconv -t utf-8 -o es-utf8.po es.po
	wine c:/Userdata/gettext/msgfmt.exe es-utf8.po
