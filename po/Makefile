CSRC := $(wildcard ../lib/src/*.c ../cli/*.c ../gui/*.c ../plugin/*.c)
VERSION = $(cat ../lib/src/version.h | grep version_ | \
	awk '{ print $NF}' | sed s/[\"\;]//g)

gretl.pot: gretl.po pottop
	./maketop && tail +10 $< >> $@

gretl.po: $(CSRC)
	xgettext -d gretl -k_ $^

