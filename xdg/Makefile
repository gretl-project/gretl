# gretl's mime types:
#  application/gretl.data (*.gdt, gretl-data-icon.png)
#  application/gretl.session (*.gretl, gretl-session-icon.png)
#  text/gretl.script (*.inp, gretl-script-icon.png)

install: gretl.desktop gretl.xml gretl-data-icon.png gretl-session-icon.png gretl-script-icon.png
	xdg-desktop-menu install gretl.desktop --novendor
	xdg-mime install gretl.xml
	xdg-icon-resource install --context mimetypes --size 32 gretl-data-icon.png application-gretl.data
	xdg-icon-resource install --context mimetypes --size 32 gretl-session-icon.png application-gretl.session
	xdg-icon-resource install --context mimetypes --size 32 gretl-script-icon.png text-gretl.script
