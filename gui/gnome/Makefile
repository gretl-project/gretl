helpdir = `gnome-config --prefix`/share/gnome/help/gretl/C
prog = gretl

html: $(prog).sgml
	db2html $(prog).sgml
	cp -ap figures $(prog)

pdf: $(prog).sgml
	db2pdf $(prog).sgml

ps:  $(prog).sgml
	db2ps $(prog).sgml       

install: 
	echo "Install routine not properly ready yet"
	cp -ap gretl/* $(helpdir)

clean:
	rm -rf $(prog).junk
	rm -rf $(prog)
	rm -f $(prog).ps $(prog).pdf $(prog).tex 


	

