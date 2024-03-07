dnl
dnl Process this file with autoconf to produce a configure script.
dnl

AC_INIT(lib/src/version.h)
AC_PREREQ(2.12)
AC_CONFIG_HEADER(config.h)
AC_CONFIG_AUX_DIR($srcdir/tools)

PACKAGE=gretl

VERSION=`grep GRETL_VERSION $srcdir/lib/src/version.h | awk '{ print $NF }' | sed -e 's/\"//g'`
LIBGRETL_CURRENT=`grep LIBGRETL_CURRENT $srcdir/lib/src/version.h | awk '{ print $NF }'`
LIBGRETL_REVISION=`grep LIBGRETL_REVISION $srcdir/lib/src/version.h | awk '{ print $NF }'`
LIBGRETL_AGE=`grep LIBGRETL_AGE $srcdir/lib/src/version.h | awk '{ print $NF }'`
LIBVERSION="$LIBGRETL_CURRENT.$LIBGRETL_AGE.$LIBGRETL_REVISION"

AC_DEFINE_UNQUOTED(PACKAGE, "$PACKAGE")
AC_DEFINE_UNQUOTED(VERSION, "$VERSION")
AC_DEFINE_UNQUOTED(LIBVERSION, "$LIBVERSION")
AC_SUBST(PACKAGE)
AC_SUBST(VERSION)
AC_SUBST(LIBVERSION)

AC_SUBST(LIBGRETL_CURRENT)
AC_SUBST(LIBGRETL_REVISION)
AC_SUBST(LIBGRETL_AGE)

dnl PWD=`pwd`
dnl ABS_SRCDIR=`(cd $srcdir && pwd)`
dnl if test $PWD = $ABS_SRCDIR ; then
dnl  echo "Building of gretl in the main source tree is not supported. Please"
dnl  echo "create a build directory either parallel to the gretl source or"
dnl  echo "under it, then configure and make from there."
dnl  exit 1
dnl fi

echo "configuring for $PACKAGE $VERSION (library version $LIBVERSION)"

dnl
dnl Variables controlling the gretl build
dnl
build_gui="yes"
prefer_gtk2="no"
build_docs="no"
build_addons="yes"
build_editor="no"
use_lucida="no"
build_po="no"
have_readline="no"
new_readline="no"
have_zlib="no"
gtk_version="none"
use_xdg="yes"
use_xdg_utils="no"
try_xdg_utils="yes"
use_gsf="no"
try_gmp="yes"
have_gmp="no"
have_mpfr="no"
gmp_fail="no"
mpfr_fail="no"
have_lapack="no"
have_fftw3="no"
seek_inet_addr="no"
have_tramo="yes"
have_x12a="yes"
have_sourceview="no"
gtksv_completion="no"
have_gnu_regex="no"
have_libxml2="no"
have_xslt="no"
try_json_glib="yes"
have_json_glib="no"
try_odbc="no"
have_odbc="no"
try_mpi="yes"
have_mpi="no"
try_libR="yes"
have_libR="no"
try_bc="no"
check_gnuplot="yes"
use_curl="yes"
darwin_build="no"
win32_build="no"
win32pkg="no"
win32bit="no"
macpkg="no"
quiet_build="no"
swap_ends="no"
pkg_build="no"
mac_native="no"
mac_integration="no"
gui_disabled="no"
gp_3d="yes"
MAKE="make"
XML_LIBS=
FONTREQ="\\RequirePackage{ae}"
RT_LIB=

AC_CHECK_HEADERS(libintl.h,
	try_nls="yes",
	try_nls="no"
)

if test [ "$try_nls" = "yes" ]
then
   AM_GNU_GETTEXT([external],,)
fi

dnl
dnl Process options supplied to the configure script
dnl
AC_ARG_ENABLE(gui,
[  --enable-gui            Enable gui [[default=yes]]],
[if test "$enableval" = "yes"
then
  build_gui=yes
else
  build_gui=no
  gui_disabled=yes
fi])

AC_ARG_ENABLE(gtk2,
[  --enable-gtk2           Use GTK 2.0 in preference to GTK 3.0 [[default=no]]],
[if test "$enableval" = "yes"
then
  prefer_gtk2=yes
else
  prefer_gtk2=no
fi])

AC_ARG_ENABLE(build-doc,
[  --enable-build-doc      Enable building of gretl docs [[git only]]],
[if test "$enableval" = "yes"
then
  if ! test -d $srcdir/.git ; then
    echo "The --enable-build-doc option is valid only when building from git"
    exit 1
  else
    build_docs=yes
  fi
else
  build_docs=no
fi])

AC_ARG_ENABLE(build-addons,
[  --enable-build-addons   Enable building of gretl addons [[default=yes]]],
[if test "$enableval" = "yes"
then
  build_addons=yes
else
  build_addons=no
fi])

AC_ARG_ENABLE(build-editor,
[  --enable-build-editor   Enable building of gretl_edit [[default=no]]],
[if test "$enableval" = "yes"
then
  build_editor=yes
else
  build_editor=no
fi])

AC_ARG_ENABLE(quiet-build,
[  --enable-quiet-build    Enable quiet building],
[if test "$enableval" = "yes"
then
  quiet_build=yes
fi])

AC_ARG_ENABLE(swap-ends,
[  --enable-swap-ends      Build FRED database with swapped endianness ],
[if test "$enableval" = "yes"
then
  swap_ends=yes
fi])

AC_ARG_ENABLE(xdg,
[  --disable-xdg           Don't Install XDG files],
[if test "$enableval" = "no"
then
  use_xdg=no
  try_xdg_utils=no
fi])

AC_ARG_ENABLE(json,
[  --disable-json          Don't try to include JSON support],
[if test "$enableval" = "no"
then
  try_json_glib=no
fi])

AC_ARG_ENABLE(gmp,
[  --disable-gmp           Don't try to include GMP support],
[if test "$enableval" = "no"
then
  try_gmp=no
fi])

AC_ARG_ENABLE(xdg-utils,
[  --enable-xdg-utils      Use xdg-utils for installing XDG files [[default=auto]]],
[if test "$enableval" = "no"
then
  try_xdg_utils=no
fi])

AC_ARG_ENABLE(pkgbuild,
[  --enable-pkgbuild       For building relocatable package (e.g. OS X) ],
[if test "$enableval" = "yes"
then
  AC_DEFINE(PKGBUILD)
  pkg_build="yes"
fi])

AC_ARG_ENABLE(gnuplot-checks,
[  --disable-gnuplot-checks do not run checks for gnuplot],
[if test "$enableval" = "yes"
then
  check_gnuplot=yes
  echo "check_gnuplot = yes"
else
  check_gnuplot=no
  echo "check_gnuplot = no"
fi])

AC_ARG_ENABLE(gnuplot-3d,
[  --disable-gnuplot-3d     assume gnuplot can't do interactive 3D plots],
[if test "$enableval" = "yes"
then
  gp_3d=yes
  echo "test gp3d = yes"
else
  gp_3d=no
  echo "test gp3d = no"
fi])

AC_ARG_ENABLE(www,
[  --disable-www           do not link to libcurl for http],
[if test "$enableval" = "yes"
then
  use_curl=yes
else
  if test "$gui_disabled" = "no" ; then
    echo ""
    echo "The --disable-www option is available only if you also"
    echo "give the option --disable-gui"
    exit 1
  else
    use_curl=no
  fi
fi])

AC_ARG_WITH(lucida,
[  --with-lucida           Use Lucida fonts for manual [[git only]]],
[ case "$withval" in
  yes)
        FONTREQ="\\RequirePackage[[T1]]{fontenc,lucidabr}"
	use_lucida=yes
        ;;
  *)
        ;;
  esac ],
)
AC_SUBST(FONTREQ)

AC_ARG_WITH(gmake,
[  --with-gmake            Use GNU make explicitly [[default=no]]],
[ case "$withval" in
  yes)
        MAKE="gmake"
        ;;
  *)
        ;;
  esac ],
)
AC_SUBST(MAKE)

AC_ARG_WITH(readline,
[  --with-readline        Use readline library [[default=auto]]],
if test "x${withval}" = "xno"
then
  try_readline=no
else
  try_readline=yes
fi,
try_readline=yes)

AC_ARG_WITH(gsf,
[  --with-gsf             Use libgsf [[default=auto]]],
if test "x${withval}" = "xno"
then
  try_gsf=no
else
  try_gsf=yes
fi,
try_gsf=yes)

AC_ARG_WITH(x-12-arima,
[  --with-x-12-arima      include X-12-ARIMA support [[default=yes]]],
if test "x${withval}" = "xno"
then
  have_x12a="no"
else
  AC_DEFINE(HAVE_X12A)
  HAVE_X12A=1
fi,
AC_DEFINE(HAVE_X12A)
HAVE_X12A=1)
AC_SUBST(HAVE_X12A)

AC_ARG_WITH(tramo-seats,
[  --with-tramo-seats     include TRAMO/SEATS support [[default=yes]]],
if test "$withval" = "no"
then
  have_tramo="no"
else
  AC_DEFINE(HAVE_TRAMO)
  HAVE_TRAMO=1
fi,
AC_DEFINE(HAVE_TRAMO)
HAVE_TRAMO=1)
AC_SUBST(HAVE_TRAMO)

AC_ARG_WITH(libR,
[  --with-libR            include libR support [[default=auto]]],
if test "$withval" = "no"
then
  try_libR=no
fi,
try_libR=yes)

AC_ARG_WITH(mpi,
[  --with-mpi             include MPI support [[default=auto]]],
if test "$withval" = "no"
then
  try_mpi=no
fi,
try_mpi=yes)

AC_ARG_WITH(odbc,
[  --with-odbc            include ODBC support],
if test "$withval" = "yes"
then
  try_odbc=yes
fi,
try_odbc=no)

AC_PROG_CC

AC_CANONICAL_HOST
AC_PROG_INSTALL

if test "$CC" = "icx" ; then
  AX_APPEND_FLAG([-fp-model precise])
fi

dnl check for Windows build
case $host_os in
  mingw*)
    win32_build=yes
    if test $pkg_build = yes ; then
      win32pkg=yes
    fi
    if test "x$WINDRES" = x ; then
      WINDRES=windres
    fi
    AC_SUBST(WINDRES)
  ;;
esac

case $host_cpu in
  i?86*)
    win32bit=yes
  ;;
esac

AC_DISABLE_STATIC
AC_LIBTOOL_DLOPEN
AC_PROG_LIBTOOL
AC_PROG_LN_S
AC_PROG_GREP

AC_CHECK_LIB(m, sin)
AC_CHECK_LIB(c, fopen)

dnl Checks for header files, etc.
AC_HEADER_SYS_WAIT
AC_CHECK_HEADERS(stdint.h unistd.h dirent.h fnmatch.h fenv.h byteswap.h sys/times.h)
AC_CHECK_HEADERS(libproc.h sys/proc_info.h)
AC_C_CONST
AC_CHECK_FUNCS(vasprintf)
AC_CHECK_FUNCS(posix_memalign)
AC_CHECK_SIZEOF(int)
AC_C_BIGENDIAN

AC_ARG_VAR(HOSTCC, host compiler (for cross-build))
AC_ARG_VAR(MKLANG, builder for gtksourceview files (cross-build))
AC_ARG_VAR(MKNEWS, builder for NEWS file (cross-build))
AC_ARG_VAR(MKPNGLIST, builder for icons index (cross-build))
AC_ARG_VAR(COMPRES, GLib resource compiler (cross-build))
AC_ARG_VAR(MPICC, path to MPI compiler wrapper)
AC_ARG_VAR(MPILINK, linker commands for MPI (cross-build))
AC_ARG_VAR(OMP_LIB, linker command for OpenMP (auto))

AC_C_SSE2
AC_C_AVX
AC_C_OPENMP

dnl check for SIMD intrinsics headers
if test "x$sse2_result" = "xyes" || test "x$avx_result" = "xyes" ; then
  AC_CHECK_HEADERS(immintrin.h)
fi

dnl ensure that the source will get rebuilt if openmp
dnl is turned on or off, by flagging this in config.h
if test "$ac_openmp_result" = "yes" ; then
  AC_DEFINE(OPENMP_BUILD)
fi

dnl Check for MPI apparatus?
if test "$try_mpi" = "yes" ; then
  if test "x${MPICC}" = x ; then
    AC_DEFUN([AC_PROG_MPICC],
    [AC_CHECK_PROG(HAVE_MPICC,mpicc,yes)])

    AC_PROG_MPICC
  else
    HAVE_MPICC="yes"
  fi
  if test "x${HAVE_MPICC}" = "xyes" ; then
    if test $win32_build = "yes" ; then
      MPI_LIBS="-lmsmpi"
      have_mpi="yes"
      AC_DEFINE(HAVE_MPI)
    else
      AC_C_MPI
    fi
  fi
fi

dnl
dnl Check for libreadline for use with gretlcli
dnl
if test "x$CROSS_READLINE_LIBS" != x ; then
   READLINE_LIBS="$CROSS_READLINE_LIBS"
   have_readline=yes
   AC_DEFINE(HAVE_READLINE)
   AC_DEFINE(NEW_READLINE)
   AC_SUBST(READLINE_CFLAGS)
   AC_SUBST(READLINE_LIBS)
elif test "$try_readline" = "yes" ; then
   AM_PATH_READLINE(have_readline="yes", have_readline="no")
   if test "$have_readline" = "no" ; then
     READLINE_LIBS=""
     READLINE_CFLAGS=""
   fi
fi

dnl
dnl Check for C++ compiler (for svm plugin)
dnl
AC_PROG_CXX

dnl
dnl Check for XDG for handling of gretl mime-types
dnl
if test "$try_xdg_utils" = "yes" ; then
  AC_DEFUN([AC_PROG_XDG_MIME],
  [AC_CHECK_PROG(XDG_MIME,xdg-mime,yes)])

  AC_PROG_XDG_MIME
  if test x"${XDG_MIME}" = xyes ; then
    use_xdg_utils=yes
  fi
fi

dnl
dnl Check for zlib
dnl
AC_CHECK_LIB(z, gzopen,have_zlib="yes" ; \
  ZLIB="-lz" ; AC_DEFINE(HAVE_ZLIB),,)
AC_SUBST(have_zlib)

dnl
dnl Check for GMP and MPFR
dnl
if test "$try_gmp" = "yes" ; then
  AM_PATH_GMP(4.0.1, have_gmp="yes", gmp_fail="yes")
  AM_PATH_MPFR(2.2.0, have_mpfr="yes", mpfr_fail="yes")
  if test "$have_gmp" = "yes" -a "$have_mpfr" = "yes" ; then
    AC_DEFINE(HAVE_GMP)
  fi
fi

dnl
dnl Check for LAPACK
dnl
case $host_os in
  *arwin*)
    if test "x$LAPACK_LIBS" = "x" ; then
       LAPACK_LIBS="-Wl,-framework -Wl,Accelerate"
    fi
    have_lapack="yes"
    darwin_build="yes"
    if test $pkg_build = yes ; then
       macpkg="yes"
    fi
    AC_DEFINE(OS_OSX)
  ;;
  *mingw*)
    have_lapack="yes"
  ;;
  *cygwin*)
    AM_PATH_LAPACK(, have_lapack="yes")
    if test "$have_lapack" = "yes" ; then
      LAPACK_LIBS="${LAPACK_LIBS} -XCClinker,-lgfortran"
    fi
  ;;
  *)
    AM_PATH_LAPACK(, have_lapack="yes")
  ;;
esac

dnl More detail for macOS build
if test "$darwin_build" = "yes" ; then
  echo "darwin host_os: $host_os"
  case $host_os in
    *darwin19* | *darwin2*)
      MACLIB="-Wl,-framework -Wl,CoreServices"
      AC_SUBST(MACLIB)
      AC_DEFINE(USE_SSE2)
      sse2_result="yes"
      echo "macOS: using CoreServices, not Carbon"
      ;;
    *)
      MACLIB="-Wl,-framework -Wl,Carbon"
      AC_SUBST(MACLIB)
      AC_DEFINE(USE_CARBON)
      echo "macOS: using Carbon"
      ;;    
  esac
  if test "$pkg_build" = "yes" ; then
    DARWIN_RPATH="@executable_path/../lib"
    AC_SUBST(DARWIN_RPATH)
  fi
fi

dnl
dnl Is MMAP available?
dnl
AC_MSG_CHECKING([for mmap])
AC_CACHE_VAL(ac_cv_func_mmap_ok,
 [AC_TRY_LINK(
 changequote(<<, >>)dnl
 <<
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
 >>,
 changequote([, ])dnl
 [mmap((void *)0, 0, PROT_READ, 0, 0, 0);],
 ac_cv_func_mmap_ok=yes,
 ac_cv_func_mmap_ok=no)] )
AC_MSG_RESULT($ac_cv_func_mmap_ok)
if test $ac_cv_func_mmap_ok = yes
then
  AC_DEFINE(HAVE_MMAP)
  if test $darwin_build = no ; then
    RT_LIB="-lrt"
  fi
fi

dnl
dnl Check for unixODBC?
dnl
if test "$try_odbc" = "yes" ; then
  case $host_os in
    *arwin*)
      ODBC_LIBS="-liodbc"
      have_odbc="yes"
    ;;
    mingw*)
      ODBC_LIBS="-lodbc32"
      have_odbc="yes"
    ;;
    *)
      AM_PATH_ODBC(,
        have_odbc="yes",
        have_odbc="no"
      )
    ;;
  esac
fi

dnl
dnl Check for gnuplot and its PNG capacity, unless the
dnl user has said not to
dnl
if test "$check_gnuplot" = "no" ; then
  have_gnuplot=yes
else
  AC_DEFUN([AC_PROG_GNUPLOT],
  [AC_CHECK_PROG(GNUPLOT,gnuplot,yes)])

  AC_PROG_GNUPLOT
  if test x"${GNUPLOT}" = xyes ; then
    have_gnuplot=yes
  fi

  if test "$build_gui" = "yes" ; then
    have_gnuplot=no
    AC_MSG_CHECKING(for gnuplot >= 5.2 with cairo support)
    echo "set term pngcairo" | GNUTERM=dumb gnuplot 2>/dev/null && have_gnuplot=yes
    if test "$have_gnuplot" = yes ; then
      GPVERSION=`gnuplot --version | awk '{print $2}'`
      if test "x$GPVERSION" = x ; then
        have_gnuplot=no
      else
        GPMAJOR=`echo $GPVERSION | cut -c1 2>/dev/null`
	GPMINOR=`echo $GPVERSION | cut -c3 2>/dev/null`
        if test $GPMAJOR -eq 5 && test $GPMINOR -ge 2; then
          AC_MSG_RESULT(yes)
	elif test $GPMAJOR -gt 5; then
	  AC_MSG_RESULT(yes)
        else
          have_gnuplot=no
        fi
      fi
    fi
    if test "$have_gnuplot" = no ; then
      AC_MSG_RESULT(no)
      echo
      echo "* gretl needs gnuplot >= 5.2, with cairo support."
      echo "* The current version of gnuplot is available from www.gnuplot.info"
      exit 1
    fi
  fi
fi

dnl
dnl Check for pdflatex executable
dnl
AC_DEFUN([AC_PROG_PDFLATEX],
[AC_CHECK_PROG(PDFLATEX,pdflatex,yes)])

AC_PROG_PDFLATEX
test x"${PDFLATEX}" = xyes && AC_DEFINE(HAVE_PDFLATEX)
if test x"${PDFLATEX}" = xyes ; then
  have_pdflatex="yes"
else
  have_pdflatex="no"
fi

dnl
dnl Preliminary check for pkg-config
dnl
if test -z "$PKG_CONFIG"; then
   AC_PATH_PROG(PKG_CONFIG, pkg-config, no)
fi
if test "$PKG_CONFIG" = "no" ; then
   echo "
Please install pkg-config and then reconfigure gretl.
pkg-config is available from http://www.freedesktop.org/
"
   exit 1
fi

orig_CFLAGS=$CFLAGS

dnl
dnl fftw3
dnl
PKG_CHECK_MODULES(FFTW, fftw3, have_fftw3="yes", have_fftw3="no")

dnl
dnl libcurl
dnl
if test "$use_curl" = "yes" ; then
  if test "x$CURL_CFLAGS" != "x" || test "x$CURL_LIBS" != "x" ; then
    have_curl="yes"
  else
    PKG_CHECK_MODULES(CURL, libcurl >= 7.13.0, have_curl="yes", have_curl="no")
  fi
  if test "$have_curl" = "yes" ; then
    AC_DEFINE(USE_CURL)
  fi
fi

dnl
dnl glib 2.0 (needed even if not building GUI)
dnl
PKG_CHECK_MODULES(GLIB, glib-2.0 >= 2.46.0, have_glib2="yes", have_glib2="no")
if test "${have_glib2}" = "no" ; then
   echo "
Please install glib-2.0 >= 2.46.0 and then reconfigure gretl.
glib-2.0 is available from http://www.gtk.org/
"
   exit 1
fi

dnl
dnl see if we have json-glib-1.0 (not a hard dependency, yet)
dnl
if test "$try_json_glib" = "yes" ; then
   PKG_CHECK_MODULES(JSON_GLIB, json-glib-1.0 >= 0.14,
   have_json_glib="yes",
   have_json_glib="no")
   if test $have_json_glib = "no" ; then
echo "
The json-glib-1.0 library was not found, or the version found
was not new enough. While we recommend including JSON support
you can disable this by using the option --disable-json when
running configure.
"
      exit 1
   fi
fi

dnl
dnl Initial check for GTK availability
dnl
if test "$prefer_gtk2" = "yes" ; then
    PKG_CHECK_MODULES(GTK2, gtk+-2.0 >= 2.20.0,
      GTK_CFLAGS=$GTK2_CFLAGS
      GTK_LIBS=$GTK2_LIBS
      CFLAGS="$CFLAGS $GTK_CFLAGS"
      LIBS="$LIBS $GTK_LIBS"
      gtk_version="2.0",
      echo "GTK 2 >= 2.20.0 not found"
    )
elif test "$build_gui" = "yes" ; then
    PKG_CHECK_MODULES(GTK3, gtk+-3.0 >= 3.20.0,
      GTK_CFLAGS=$GTK3_CFLAGS
      GTK_LIBS=$GTK3_LIBS
      CFLAGS="$CFLAGS $GTK_CFLAGS"
      LIBS="$LIBS $GTK_LIBS"
      AC_DEFINE(USE_GTK3)
      gtk_version="3.0",
      echo "GTK 3 >= 3.20.0 not found"
    )
fi

dnl GTK 2.0 fallback in case preference not given
if test "$prefer_gtk2" = "no" && test "$gtk_version" = "none" ; then
   if test "$build_gui" = "yes" ; then
      PKG_CHECK_MODULES(GTK2, gtk+-2.0 >= 2.20.0,
        GTK_CFLAGS=$GTK2_CFLAGS
        GTK_LIBS=$GTK2_LIBS
        CFLAGS="$CFLAGS $GTK_CFLAGS"
        LIBS="$LIBS $GTK_LIBS"
        gtk_version="2.0",
        build_gui="no"
        echo "GTK 2 >= 2.20.0 not found"
    )
   fi
fi

dnl
dnl Mac: check for Quartz version of GTK
dnl
if test "$darwin_build" = "yes" ; then
  if test "$gtk_version" = "3.0" ; then
    PKG_CHECK_EXISTS(gtk+-quartz-3.0,
     echo "Found quartz-based GTK"
     AC_DEFINE(MAC_NATIVE)
     mac_native="yes",
     mac_native="no")
  elif test "$gtk_version" = "2.0" ; then
    PKG_CHECK_EXISTS(gtk+-quartz-2.0,
     echo "Found quartz-based GTK"
     AC_DEFINE(MAC_NATIVE)
     mac_native="yes",
     mac_native="no")
     PKG_CHECK_EXISTS(gretl-mac-themes,
     AC_DEFINE(HAVE_MAC_THEMES),)
  fi
fi

dnl
dnl Mac: if using Quartz GTK, require mac-integration
dnl
if test "$mac_native" = "yes" ; then
    PKG_CHECK_MODULES(GTKMAC, gtk-mac-integration >= 2.0.5,
      mac_integration="yes",
      mac_integration="no",
    )
fi

dnl
dnl Determine the gtk prefix
dnl
if test "$gtk_version" = "2.0" ; then
    GTK_PREFIX=`$PKG_CONFIG --variable=prefix gtk+-2.0`
else
    GTK_PREFIX=`$PKG_CONFIG --variable=prefix gtk+-3.0`
fi

dnl
dnl Check for installed gtksourceview (GTK 2 variant)
dnl
if test "$gtk_version" = "2.0" ; then
    PKG_CHECK_MODULES(GTKSOURCEVIEW, gtksourceview-2.0 >= 2.10.5,
      GTKSOURCEVIEW_PREFIX=`$PKG_CONFIG --variable=prefix gtksourceview-2.0`
      have_sourceview="2.0",
      echo "version 2.0 >= 2.10.5 not found"
    )
fi

dnl
dnl Check for installed gtksourceview (GTK 3 variant, current)
dnl
if test "$gtk_version" = "3.0" ; then
    PKG_CHECK_MODULES(GTKSOURCEVIEW, gtksourceview-4,
      GTKSOURCEVIEW_PREFIX=`$PKG_CONFIG --variable=prefix gtksourceview-4`
      have_sourceview="4",
        PKG_CHECK_MODULES(GTKSOURCEVIEW, gtksourceview-3.0 >= 3.20.0,
        GTKSOURCEVIEW_PREFIX=`$PKG_CONFIG --variable=prefix gtksourceview-3.0`
        have_sourceview="3.0",
        echo "version 3.0 >= 3.20.0 not found"
      )
   )
fi

if test $have_sourceview = "3.0" || test $have_sourceview = "4" || test $pkg_build = "yes" ; then
   AC_DEFINE(HAVE_GTKSV_COMPLETION)
   gtksv_completion="yes"
fi

if test x"${MSGFMT}" != x ; then
   build_po="yes"
else
   echo "*** msgfmt not found, can't build message catalogs"
fi

dnl countermand any positive results for xdg if we're not
dnl building the gretl GUI
if test "$gui_disabled" = "yes" ; then
  use_xdg=no
  use_xdg_utils=no
fi

dnl record the results of the above tests
AC_SUBST(build_gui)
AC_SUBST(gtk_version)
AC_SUBST(mac_native)
AC_SUBST(win32_build)
AC_SUBST(win32pkg)
AC_SUBST(win32bit)
AC_SUBST(macpkg)
AC_SUBST(have_sourceview)
AC_SUBST(gtksv_completion)
AC_SUBST(have_gnu_regex)
AC_SUBST(have_odbc)
AC_SUBST(have_mpi)
AC_SUBST(build_po)
AC_SUBST(use_curl)
AC_SUBST(use_xdg)
AC_SUBST(use_xdg_utils)
AC_SUBST(GLIB_CFLAGS)
AC_SUBST(GLIB_LIBS)
AC_SUBST(GTK_CFLAGS)
AC_SUBST(GTK_LIBS)
AC_SUBST(CURL_CFLAGS)
AC_SUBST(CURL_LIBS)
AC_SUBST(GTKSOURCEVIEW_PREFIX)
AC_SUBST(GTK_PREFIX)
AC_SUBST(RT_LIB)

dnl Check for gdk-pixbuf support
if test ${build_gui} = "yes" ; then
    PKG_CHECK_MODULES(GDK_PIXBUF, gdk-pixbuf-2.0 >= 0.7.0,
        pixbuf_failed="no",
	pixbuf_failed="yes"
   )
fi
if test "${pixbuf_failed}" = "yes" ; then
    echo "* gretl needs but did not find the gdk-pixbuf library."
    exit 1
fi

dnl Check for libxml2
PKG_CHECK_MODULES(XML, libxml-2.0 >= 2.5.0)

dnl Check for libxslt (for the helpfiles)
if test "${build_docs}" = "yes" ; then
  PKG_CHECK_MODULES(XSLT, libxslt >= 1.0.15,
    have_xslt="yes",
    have_xslt="no"
  )
fi

dnl Check for libR (optional)
if test "$try_libR" = "yes" ; then
  if test $win32_build = "yes" ; then
     have_libR="yes"
     AC_DEFINE(USE_RLIB)
  else
    PKG_CHECK_MODULES(RLIB, libR,
      AC_DEFINE(USE_RLIB)
      if test "$darwin_build" = "yes" ; then
        RLIBDIR="/Library/Frameworks/R.framework/Resources/lib"
      elif echo "$RLIB_LIBS" | grep ^-L > /dev/null ; then
        RLIBDIR=`echo $RLIB_LIBS | cut -b 3- | awk '{ print $1 }'`
      else
        RLIBDIR="/usr/lib"
      fi
      if test "$darwin_build" = "yes" ; then
        RLIBNAME="libR.dylib"
      else
        RLIBNAME="libR.so"
      fi
      RLIB="${RLIBDIR}/${RLIBNAME}"
      AC_DEFINE_UNQUOTED(RLIBPATH, "$RLIB")
      have_libR="yes",
      have_libR="no"
    )
  fi
fi

dnl
dnl Check for libgsf-1, if requested
dnl
if test "$try_gsf" = "yes" ; then
    PKG_CHECK_MODULES(GSF, libgsf-1 >= 1.14.31,
      CFLAGS="$CFLAGS $GSF_CFLAGS"
      LIBS="$LIBS $GSF_LIBS"
      AC_DEFINE(USE_GSF)
      use_gsf="yes",
      echo "libgsf >= 1.14.31 not found"
    )
fi

dnl If building docs was requested, see if we can do it
if test "${build_docs}" = "yes" ; then
   if test "$have_xslt" = "no" || test "$have_pdflatex" = "no" ; then
      build_docs="no"
      echo "*** Can't build gretl docs: XSLT and/or pdflatex not found"
   fi
fi

dnl If building addons was requested, see if we can do it
if test "${build_addons}" = "yes" ; then
   if test "$have_pdflatex" = "no" ; then
      echo "*** Can't build documentation for addons: pdflatex not found."
      build_addons="no"
   fi
fi

AC_SUBST(build_docs)
AC_SUBST(build_addons)
AC_SUBST(build_editor)
AC_SUBST(quiet_build)
AC_SUBST(pkg_build)
AC_SUBST(swap_ends)
AC_SUBST(use_gsf)

GTK_HUSH=""
CFLAGS=-Wno-deprecated-declarations
AC_MSG_CHECKING([whether CC supports -Wno-deprecated-declarations])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([])],
   [AC_MSG_RESULT([yes])]
   [GTK_HUSH=-Wno-deprecated-declarations],
   [AC_MSG_RESULT([no])]
)
AC_SUBST(GTK_HUSH)

if test "$gtk_version" = "3.0" ; then
   plugsub="gretl-gtk3"
else
   plugsub="gretl-gtk2"
fi
AC_SUBST(plugsub)

dnl record results on various libraries
AC_SUBST(XML_CFLAGS)
AC_SUBST(XML_LIBS)
AC_SUBST(XSLT_CFLAGS)
AC_SUBST(XSLT_LIBS)
AC_SUBST(GDK_PIXBUF_CFLAGS)
AC_SUBST(GDK_PIXBUF_LIBS)
AC_SUBST(RLIB_CFLAGS)
AC_SUBST(RLIB_LIBS)
AC_SUBST(RLIB)
AC_SUBST(USE_SSE2)
AC_SUBST(USE_AVX)

dnl Determine installation prefix
if test "${prefix}" = "NONE" ; then
   prefix="/usr/local"
fi
AC_DEFINE_UNQUOTED(GRETL_PREFIX, "${prefix}")

if test "${gp_3d}" = "yes" ; then
  AC_DEFINE(GNUPLOT3D)
fi

CFLAGS="$orig_CFLAGS"

AC_CONFIG_FILES([
Makefile
gretl_sh
gretl.pc
gui/Makefile
gui/debug
xdg/gretl.desktop
xdg/gretl_edit.desktop
xdg/Makefile
cli/Makefile
cli/debug
lib/Makefile
plugin/Makefile
po/Makefile.in
share/Makefile
share/bcih/Makefile
tests/Makefile
unittests/Makefile
apidemo/Makefile
redhat/gretl.spec
utils/emacs/Makefile
utils/dbzip/Makefile
doc/gretl.sty
doc/gretl-lite.sty
])

if test "${darwin_build}" = "yes" ; then
  AC_CONFIG_FILES([
    osx/Info.plist
  ])
fi

if test "${win32_build}" = "yes" ; then
  if test "${pkg_build}" = "yes" ; then
    mingw_prefix="$MINGW_PREFIX"
    AC_SUBST(prefix)
    AC_SUBST(mingw_prefix)
    AC_CONFIG_FILES([
      win32/windist/mkiss.sh
      win32/windist/post-install.sh
      win32/windist/install-runtime.sh
      win32/windist/install-runtime-gtk3.sh
    ])
  else
    GHOME="${prefix}/share/gretl"
    GEXEPATH="${prefix}/bin"
    GEXE="${prefix}/bin/gretl.exe"
    MINGWBIN="$MINGW_PREFIX/bin"
    MSYSBIN="/bin"
    AC_PROG_SED
    AC_CHECK_PROG(CYGPATH_CHECK,cygpath,yes)
    if test x"$CYGPATH_CHECK" = x"yes" ; then
      GHOME=`cygpath -w $GHOME | $SED 's+\\\\+\\\\\\\\+g'`
      GEXEPATH=`cygpath -w $GEXEPATH | $SED 's+\\\\+\\\\\\\\+g'`
      GEXE=`cygpath -m $GEXE`
      MINGWBIN=`cygpath -m $MINGWBIN`
      MSYSBIN=`cygpath -m $MSYSBIN`
      MINGWPRE=`cygpath -m $MINGW_PREFIX`
    fi
    AC_SUBST(GHOME)
    AC_SUBST(GEXE)
    AC_SUBST(GEXEPATH)
    AC_SUBST(MINGWBIN)
    AC_SUBST(MSYSBIN)
    AC_SUBST(MINGWPRE)
    AC_CONFIG_FILES([
      win32/gretlrun.cmd
      win32/gretlmime.reg
      win32/mkshortcut.sh
    ])
  fi
fi

if test "${build_docs}" = "yes" ; then
  AC_CONFIG_FILES([
    doc/Makefile
    doc/commands/Makefile
    doc/tex/Makefile
    doc/tex_it/Makefile
    doc/tex_es/Makefile
    doc/tex_pt/Makefile
    doc/tex_gl/Makefile
    doc/tex_ru/Makefile
    doc/tex/extract.mk
    doc/reference/Makefile
  ])
fi

if test "${build_addons}" = "yes" ; then
  AC_CONFIG_FILES([
    addons/Makefile
    addons/gig/Makefile
    addons/gig/gig.spec
    addons/gig/pkg.inp
    addons/gig/doc/Makefile
    addons/HIP/Makefile
    addons/HIP/HIP.spec
    addons/HIP/pkg.inp
    addons/HIP/doc/Makefile
    addons/ivpanel/Makefile
    addons/ivpanel/ivpanel.spec
    addons/ivpanel/pkg.inp
    addons/ivpanel/doc/Makefile
    addons/SVAR/Makefile
    addons/SVAR/SVAR.spec
    addons/SVAR/pkg.inp
    addons/SVAR/doc/Makefile
    addons/dbnomics/Makefile
    addons/dbnomics/dbnomics.spec
    addons/dbnomics/pkg.inp
    addons/dbnomics/doc/Makefile
    addons/extra/extra.spec
    addons/extra/pkg.inp
    addons/extra/Makefile
    addons/extra/doc/Makefile
    addons/geoplot/Makefile
    addons/geoplot/geoplot.spec
    addons/geoplot/pkg.inp
    addons/geoplot/doc/Makefile
    addons/geoplot/examples/Makefile
    addons/regls/Makefile
    addons/regls/regls.spec
    addons/regls/pkg.inp
    addons/regls/doc/Makefile
    addons/logging/Makefile
    addons/logging/logging.spec
    addons/logging/pkg.inp
    addons/logging/doc/Makefile
    addons/KFgui/Makefile
    addons/KFgui/KFgui.spec
    addons/KFgui/pkg.inp
  ])
fi

if test "${build_editor}" = "yes" ; then
  AC_CONFIG_FILES([
    editor/Makefile
  ])
fi

AC_OUTPUT

if test "${have_zlib}" = "no" ; then
echo "
Please install zlib (compression library) and then reconfigure gretl.
zlib is available via http://www.info-zip.org/pub/infozip/zlib/
"

elif test x"${XML_LIBS}" = x ; then
echo "
Please install libxml2 and then reconfigure gretl.
libxml2 is available from http://xmlsoft.org/
"

elif test "${have_lapack}" = "no" ; then
echo "
Please install lapack and then reconfigure gretl.
Lapack is available from http://www.netlib.org/lapack/
"

elif test "${have_fftw3}" = "no" ; then
echo "
Please install fftw3 and then reconfigure gretl.
fftw3 is available from http://www.fftw.org/
"

elif test "${gmp_fail}" = "yes" ; then
echo "
Please install GMP and then reconfigure gretl.
GMP is available from http://gmplib.org/
"

elif test "${mpfr_fail}" = "yes" ; then
echo "
Please install MPFR and then reconfigure gretl.
MPFR is available from http://www.mpfr.org/
"

elif test "${have_curl}" = "no" ; then
echo "
Please install libcurl and then reconfigure gretl.
libcurl is available from http://curl.haxx.se/libcurl/
"

elif test "${have_glib2}" = "no" ; then
echo "
Please install glib-2.0 >= 2.28.0 and then reconfigure gretl.
glib-2.0 is available from http://www.gtk.org/
"

elif test "${gtk_version}" = "2.0" && test "${have_sourceview}" = "no" ; then
echo "
Please install gtksourceview-2.0 and then reconfigure gretl.
gtksourceview-2.0 is available from http://www.gnome.org/
"

elif test "${gtk_version}" = "3.0" && test "${have_sourceview}" = "no" ; then
echo "
Please install gtksourceview 3.0 or 4 and then reconfigure gretl.
gtksourceview is available from http://www.gnome.org/
"

elif test "${darwin_build}" = "yes" && test "${mac_native}" = "no" ; then
echo "
Building for macOS requires gtk+-quartz-2.0 or gtk+-quartz-3.0
"

elif test "${darwin_build}" = "yes" && test "${mac_integration}" = "no" ; then
echo "
Building for macOS requires gtk-mac-integration version 2.0.5 or higher
"

else
  if test "${use_xdg_utils}" = "yes" ; then
    xdg_utils_msg="if DESTDIR not set"
  else
    xdg_utils_msg="no"
  fi
echo "
Configuration:

  Installation path:                      ${prefix}
  Use readline library:                   ${have_readline}
  Use gnuplot for graphs:                 ${have_gnuplot}
  Use pdflatex for typesetting:           ${have_pdflatex}
  Use libgsf for zip/unzip:               ${use_gsf}
  SIMD support for RNG:                   ${sse2_result}
  OpenMP support:                         ${ac_openmp_result}
  MPI support:                            ${have_mpi}
  AVX support for arithmetic:             ${avx_result}
  Build with GTK version:                 ${gtk_version}
  Use gtksourceview version:              ${have_sourceview}
  Build gretl documentation:              ${build_docs}
  Use Lucida fonts:                       ${use_lucida}
  Build message catalogs:                 ${build_po}
  Build gretl addons:                     ${build_addons}
  Build gretl_edit:                       ${build_editor}
  X-12-ARIMA support:                     ${have_x12a}
  TRAMO/SEATS support:                    ${have_tramo}
  libR support:                           ${have_libR}
  ODBC support:                           ${have_odbc}
  GMP support:                            ${have_gmp}
  JSON parsing support:                   ${have_json_glib}
  Use xdg-utils in installation:          ${xdg_utils_msg}

  LAPACK libraries:
    ${LAPACK_LIBS}

Now type '$MAKE' to build gretl."

  if test "${build_docs}" = "yes" ; then
    echo "You can also do 'make pdfdocs' to build the PDF documentation."
  fi
  echo

  if test "${gui_disabled}" = "no" && test "${gtk_version}" = "none" ; then
    echo "*** Warning: since GTK is not available the gretl GUI will"
    echo "*** not be built, only the library and command-line program."
    echo
  fi

  if test "$ac_openmp_result" = "yes" && test "x${blas_pthreads}" != "x"; then
    echo "*** Warning: you are building gretl with OpenMP, but linking"
    echo "*** against an OpenBLAS library that uses pthreads. This mixture"
    echo "*** is not recommended. Ideally, OpenBLAS should use OpenMP."
    echo
  fi
fi
