# make windows HTML help file for gretl

topdoc = ../../doc

htmlhelp: $(topdoc)/manual.xml $(topdoc)/chapters/*.xml 
	mkdir -p gretlhlp/figures
	cp $(topdoc)/figures/*.png gretlhlp/figures
	for f in gretlhlp/figures/*.png ; do convert $$f $${f%.png}.gif ; done
	cd gretlhlp && \
	saxon ../$(topdoc)/manual.xml ../$(topdoc)/gretlhlp.xsl tex.math.in.alt=latex
	wine "c:/Program Files/HTML Help Workshop/hhc.exe" gretlhlp/gretl.hhp
	mv gretlhlp/gretl.chm ./gretl.chm && rm -rf gretlhlp

clean:
	rm -rf gretlhlp
	rm -f gretl.chm
