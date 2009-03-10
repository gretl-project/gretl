# Configuration variables governing the build of gretl for win32

# change to 'native' if building on MS Windows
BUILDTYPE = cross

ifeq ($(BUILDTYPE),cross)
  # directory containing the cross tools
  MGWDIR = /opt/cross-tools/mingw32
  # prefix for mingw tools (e.g. mingw32-gcc)
  MGW_PREFIX = mingw32-
  # msgfmt command for producing win32 messages file
  WIN32_MSGFMT = wine c:/bin/msgfmt.exe 
  # location of pkgconfig files
  PKG_CONFIG_PATH = $(MGWDIR)/lib/pkgconfig
else
  # directory containing mingw
  MGWDIR = c:/mingw
  # no prefix for gcc, etc.
  MGW_PREFIX =
  # msgfmt command for producing win32 messages file
  WIN32_MSGFMT = msgfmt.exe
  # location of pkgconfig files
  PKG_CONFIG_PATH = c:/msys/1.0/lib/pkgconfig:$(MGWDIR)/lib/pkgconfig
endif

# libxml2 includes: adjust to match your system
XML2_INC = $(MGWDIR)/libxml2
