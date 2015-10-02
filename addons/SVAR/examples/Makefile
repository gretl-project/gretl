CHECKFILES = awm.inp \
	BlQuah.inp \
	BlQuah_w_sr.inp \
	IS-LM.inp \
	plain.inp \
	simple_C.inp

check: $(CHECKFILES) 
	for f in $(CHECKFILES) ; do \
	if gretlcli -b $$f ; then echo "$$f: OK" ; \
	else echo "*** Error executing $$f" ; break ; fi \
	done
