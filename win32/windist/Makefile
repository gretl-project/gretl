HOSTCC = gcc -Wall

include ../config.mk

ifeq ($(BUILDTYPE),cross)
  DOIT = ./make_iss < MANIFEST | todos
else
  DOIT = ./make_iss < MANIFEST
endif

gretl.iss: make_iss MANIFEST
	$(DOIT) > gretl.iss

make_iss: make_iss.c
	$(HOSTCC) -o $@ $< 

versiondiff: versiondiff.c
	$(HOSTCC) -o $@ $< 

clean:
	rm -rf ./gretl 
	rm -f errlog Output/*

