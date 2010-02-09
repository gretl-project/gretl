# make OS X icns files based on the win32 gretl icons

ICODIR=../win32/icons

cp $ICODIR/icon1.ico .
cp $ICODIR/icon2.ico .
cp $ICODIR/icon3.ico .

# convert with ImageMagick
convert icon1.ico gdt.png
convert icon2.ico script.png
convert icon3.ico session.png

# see libicns at http://sourceforge.net/projects/icns/ 
png2icns gdt.icns gdt.png
png2icns script.icns script.png
png2icns session.icns session.png


