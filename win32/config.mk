# Configuration variables governing the build of gretl for win32

topsrc = ..
podir = ../$(topsrc)/po

# change to 'native' if building on MS Windows
BUILDTYPE = cross

ifeq ($(BUILDTYPE),cross)
  # directory containing the cross tools
  MGWDIR = /opt/win32/i686-w64-mingw32
  # prefix for mingw tools
  MGW_PREFIX = i686-w64-mingw32-
  # msgfmt command for producing win32 messages file
  WIN32_MSGFMT = wine c:/bin/msgfmt.exe 
  # location of pkgconfig files
  PKG_CONFIG_PATH = $(MGWDIR)/lib/pkgconfig
else
  # directory containing mingw
  MGWDIR = /mingw32
  # no prefix for gcc, etc.
  MGW_PREFIX =
  # msgfmt command for producing win32 messages file
  WIN32_MSGFMT = msgfmt.exe
  # location of pkgconfig files
  PKG_CONFIG_PATH = $(MGWDIR)/lib/pkgconfig:/usr/lib/pkgconfig
endif

# GTK version switch: set HAVE_GTK_SPINNER to "yes" if you
# are building (and running) with GTK version 2.20.0 or higher,
# otherwise set to "no".
HAVE_GTK_SPINNER = yes

# set this for an MPI-enabled build
have_mpi = yes

# lapack/blas: adjust to match your system
LAPACK_LIBS = -lopenblas.dll

# libxml2 includes: adjust to match your system
XML2_INC = -I$(MGWDIR)/include/libxml2

# libsvm includes: adjust/add if required
# SVM_INC = -I/usr/local/include

# libR includes: likewise
RLIB_CFLAGS = -I/opt/R/lib64/R/include

