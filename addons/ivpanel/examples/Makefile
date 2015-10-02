CHECKFILES = nc_crime.inp 

check: $(CHECKFILES) 
	for f in $(CHECKFILES) ; do \
	if gretlcli -b $$f > $$f.out ; then echo "$$f: OK" ; \
	else echo "*** Error executing $$f" ; break ; fi \
	done
