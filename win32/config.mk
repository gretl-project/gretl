# Configuration variables governing the build of gretl for win32

# directory containing the cross tools
MGWDIR = /opt/cross-tools/mingw32
# MGWDIR - c:/mingw

# prefix for mingw tools (e.g. mingw32-gcc)
MGW_PREFIX = mingw32-
# MGW_PREFIX = 

# mingw include dir
MGW_INC = $(MGWDIR)/include

# libxml2 includes: adjust to match your system
XML2_INC = $(MGW_INC)/libxml2

# msgfmt command for producing win32 messages file
WIN32_MSGFMT = wine c:/bin/msgfmt.exe
# WIN32_MSGFMT = msgfmt.exe

# pkgconfig path
PKG_CONFIG_PATH = $(MGWDIR)/lib/pkgconfig
