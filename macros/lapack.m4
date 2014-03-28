# Configure paths for lapack
# Allin Cottrell <cottrell@wfu.edu>, last updated September 2013

dnl AM_PATH_LAPACK([MINIMUM-VERSION, [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl Test for LAPACK, and define LAPACK_CFLAGS and LAPACK_LIBS.
dnl
AC_DEFUN([AM_PATH_LAPACK],
[dnl 
AC_ARG_WITH(lapack-prefix,[  --with-lapack-prefix=PFX   Prefix where LAPACK is installed (optional)],
            lapack_config_prefix="$withval", lapack_config_prefix="")
AC_ARG_VAR([LAPACK_LIBS],[linker flags for lapack, overriding auto-detection])

  if test x"${LAPACK_LIBS}" = x ; then 
     AC_MSG_CHECKING(for libgfortran, libg2c or libf2c)
     AC_CHECK_LIB(gfortran,_gfortran_abort,FLIB="-lgfortran",FLIB="none")
     if test $FLIB = "none" ; then
        AC_CHECK_LIB(g2c,c_sqrt,FLIB="-lg2c",FLIB="none")
     fi
     if test $FLIB = "none" ; then
        AC_CHECK_LIB(f2c,c_sqrt,FLIB="-lf2c",FLIB="none")
     fi
     if test $FLIB = "none" ; then
        echo "*** Couldn't find libgfortran, libg2c or libf2c"
        FLIB=""
     fi
  fi

  AC_MSG_CHECKING(for LAPACK)
  if test x"${LAPACK_LIBS}" = x ; then
     if test "x$lapack_config_prefix" = x ; then
        LAPACK_LIBS="-llapack -lblas $FLIB"
     else
        LAPACK_LIBS="-L$lapack_config_prefix -llapack -lblas $FLIB"
     fi
  fi

  ac_save_CFLAGS="$CFLAGS"
  ac_save_LIBS="$LIBS"
  CFLAGS="$LAPACK_CFLAGS $CFLAGS"
  LIBS="$LAPACK_LIBS $LIBS"

dnl
dnl Check the installed LAPACK.
dnl
  rm -f conf.lapacktest
  AC_TRY_RUN([
#include <stdlib.h>
#include "gretl_f2c.h"

int 
main ()
{
  integer ispec;
  real zero = 0.0;
  real one = 1.0;

  ieeeck_(&ispec, &zero, &one);
  system ("touch conf.lapacktest");
  return 0;
}
],, no_lapack=yes,[echo $ac_n "cross compiling; assumed OK... $ac_c"])
  CFLAGS="$ac_save_CFLAGS"
  LIBS="$ac_save_LIBS"

  if test "x$no_lapack" = x ; then
     AC_MSG_RESULT(yes)
     ifelse([$2], , :, [$2])     
  else
     AC_MSG_RESULT(no)
     if test -f conf.lapacktest ; then
       :
     else
       echo "*** Could not run LAPACK test program, checking why..."
       CFLAGS="$LAPACK_CFLAGS $CFLAGS $OPENMP_CFLAGS"
       LIBS="$LIBS $LAPACK_LIBS"
       AC_TRY_LINK([
#include <stdio.h>
],     [ return (1); ],
       [ echo "*** The test program compiled, but did not run. This usually means"
         echo "*** that the run-time linker is not finding LAPACK. If it is not"
         echo "*** finding LAPACK, you'll need to set your LD_LIBRARY_PATH "
         echo "*** environment variable, or edit /etc/ld.so.conf to point"
         echo "*** to the installed location.  Also, make sure you have run"
         echo "*** ldconfig if that is required on your system."
         echo "***" ],
       [ echo "*** The test program failed to compile or link. See config.log for the"
         echo "*** exact error that occured. This may mean LAPACK was incorrectly installed"
         echo "*** or that you have moved LAPACK since it was installed." ])
         CFLAGS="$ac_save_CFLAGS"
         LIBS="$ac_save_LIBS"
     fi
     LAPACK_CFLAGS=""
     LAPACK_LIBS=""
     ifelse([$3], , :, [$3])
  fi
  AC_SUBST(LAPACK_CFLAGS)
  AC_SUBST(LAPACK_LIBS)
  AC_SUBST(FLIB)
  rm -f conf.lapacktest
])
