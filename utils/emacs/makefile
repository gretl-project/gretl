# adjust as needed
MKLANG = ~/src/build/gui2/mklang

gretl.el : gretl.el.top gretl.el.bottom gretl.el.middle
	cat gretl.el.top > gretl.el
	cat gretl.el.middle >> gretl.el
	cat gretl.el.bottom >> gretl.el

gretl.el.middle :
	$(MKLANG) --emacs > $@

