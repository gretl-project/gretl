CHECKFILES = example1.inp \
	example2.inp \
	example4.inp \
	example5.inp \
	example6.inp \
	test2.inp \
	test3.inp \
	test4.inp \
	test.inp

check: $(CHECKFILES) 
	for f in $(CHECKFILES) ; do \
	if gretlcli -b $$f ; then echo "$$f: OK" ; \
	else echo "*** Error executing $$f" ; break ; fi \
	done
