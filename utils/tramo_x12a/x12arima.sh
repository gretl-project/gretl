#!/bin/sh

# grab and build x12a for Linux

# You can supply an installation root, e.g. /usr/local or /opt,
# on the command line.  In any case it goes into a subdir
# named "x12arima".

X12ROOT=$1
if [ "x$X12ROOT" = "x" ] ; then
   echo -n "Path for installation? (default is /opt): " 
   read X12ROOT
   if [ "x$X12ROOT" = "x" ] ; then
      X12ROOT="/opt"
   fi  
fi

echo "Installing in $X12ROOT"

SITE="http://www.census.gov/ts/x12a/final/unix"
SRC="omegasrc.tar.Z"
EXAMPLES="finexam.tar.Z"
PSDOCS="omegaps.tar.Z"

THISDIR=`pwd`
DEST=$X12ROOT/x12arima

# create a directory to hold the distribution and cd to it
mkdir -p $DEST
mkdir -p $DEST/src
cd $DEST

# grab the relevant archives and extract
wget $SITE/$SRC
cd $DEST/src && tar xvfz $DEST/$SRC && rm $DEST/$SRC && cd $DEST
wget $SITE/$EXAMPLES 
tar xvfz $EXAMPLES && rm $EXAMPLES
wget $SITE/$PSDOCS
tar xvfz $PSDOCS && rm $PSDOCS

# move the aux files up from src
mv ./src/x12a.mdl ./src/test.spc $DEST

# patch the Makefile for Linux
cat > $THISDIR/x12a_make.diff << EOF
--- src/Makefile.orig	Mon May  6 11:22:42 2002
+++ src/Makefile	Thu Dec  5 09:47:12 2002
@@ -1,12 +1,11 @@
-#
-FC        = f77
-LINKER    = f77
-PROGRAM         = x12a
+FC        = g77
+LINKER    = g77
+PROGRAM   = ../x12a
 DEST      = .
 EXTHDRS         =
-FFLAGS    = -c -C -g -Nl99
+FFLAGS    = -O2
 HDRS            =
-LDFLAGS   = -g -dn -Bstatic
+LDFLAGS   = -s
 LDMAP     = 
 LIBS      =
 MAKEFILE  = Makefile
@@ -208,11 +207,11 @@
 \$(PROGRAM):     \$(OBJS) \$(LIBS)
 	\$(LINKER) -o \$@ \$(OBJS) \$(LDMAP) \$(LIBS) \$(LDFLAGS)
 
-clean:;         @del -f \$(OBJS)
+clean:;         @rm -f \$(OBJS)
 
 install:   \$(PROGRAM)
 	@echo Installing \$(PROGRAM) in \$(DEST)
-	@if not \$(DEST)x==.x copy \$(PROGRAM) \$(DEST)
+	@install \$(PROGRAM) \$(DEST)
 ### OPUS MKMF:  Do not remove this line!  Automatic dependencies follow.
 
 aaamain.o:  cchars.i error.cmn hiddn.cmn lex.i ssap.prm stdio.i title.cmn \\
EOF

patch -p0 < $THISDIR/x12a_make.diff

# remove DOS file termination chars
echo "Fixing DOS-terminated source files..."
for f in ./src/*.f ./src/*.var ./src/*.prm ./src/*.cmn ; do
   cat $f | tr -d '\032' > $f.tmp && mv $f.tmp $f
done

# build the program, then compress the sources
make -C src 
strip ./x12a
make -C src clean
tar cvfz src.tgz src && rm -rf src




