CSRC := $(wildcard ../lib/src/*.c ../cli/*.c ../gui/*.c ../plugin/*.c)

version := $(shell grep version_ ../lib/src/version.h | sed s/[^0-9\.]*//g)

gretl.pot: gretl.po pottop
	cat pottop | sed s/VERSION/$(version) > $@
	tail +10 $< >> $@

gretl.po: $(CSRC)
	xgettext -d gretl -k_ $^

