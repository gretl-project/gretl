HOSTCC = gcc -W -Wall

include ../../config.mk

ifeq ($(BUILDTYPE),cross)
  TODOS = todos
endif

gretl.iss: make_iss
	./make_iss < MANIFEST | $(TODOS) > gretl.iss 

make_iss: make_iss.c
	$(HOSTCC) -o $@ $< 

versiondiff: versiondiff.c
	$(HOSTCC) -o $@ $< 

clean:
	rm -rf ./gretl 
	rm -f errlog Output/*
