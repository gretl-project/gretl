# Configure paths for lapack
# Allin Cottrell <cottrell@wfu.edu>, last updated January 2020
#
# This is designed to handle openblas (in a standard location) or
# liblapack plus libblas (again, in a standard location). For
# fancier setups, supply LAPACK_LIBS.

dnl AM_PATH_LAPACK([MINIMUM-VERSION, [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl Test for LAPACK, and define LAPACK_LIBS.
dnl
AC_DEFUN([AM_PATH_LAPACK],
[dnl
AC_ARG_VAR([LAPACK_LIBS],[linker flags for lapack, overriding auto-detection])

if test x"${LAPACK_LIBS}" = x ; then
   # check for OpenMP openblas using OpenBLAS naming
   AC_CHECK_LIB(openblas_omp,ddot_,LAPACK_LIBS="-lopenblas_omp",LAPACK_LIBS="none")
   if test $LAPACK_LIBS = "none" ; then
      # OpenMP OpenBLAS, Fedora naming
      AC_CHECK_LIB(openblaso,ddot_,LAPACK_LIBS="-lopenblaso",LAPACK_LIBS="none")
   fi
   if test $LAPACK_LIBS = "none" ; then
      # OpenBLAS at all?
      AC_CHECK_LIB(openblas,ddot_,LAPACK_LIBS="-lopenblas",LAPACK_LIBS="none")
   fi
   if test $LAPACK_LIBS = "none" ; then
      # fall back to liblapack plus libblas
      AC_CHECK_LIB(blas,ddot_,LAPACK_LIBS="-llapack -lblas",LAPACK_LIBS="none")
   fi
   if test $LAPACK_LIBS = "none" ; then
      echo "*** Couldn't find libopenblas or libblas"
      LAPACK_LIBS=""
   fi
else
   AC_MSG_CHECKING(for lapack)
fi

ac_save_LIBS="$LIBS"
LIBS="$LAPACK_LIBS $LIBS"

dnl Check a LAPACK function in the specified or detected setup
dnl
  rm -f conf.lapacktest
  AC_TRY_RUN([
#include <stdlib.h>
int dpotrf_(char *, int *, double *, int *, int *);
int main (void) {
  char uplo = 'L';
  int one = 1;
  int info = 0;
  double x = 1;
  dpotrf_(&uplo, &one, &x, &one, &info);
  system("touch conf.lapacktest");
  return 0;
}
],,no_lapack=yes,[echo $ac_n "cross compiling; assumed OK... $ac_c"])

if test "x$no_lapack" = x ; then
  AC_MSG_RESULT(yes)
  ifelse([$2], , :, [$2])
else
  AC_MSG_RESULT(no)
  if test -f conf.lapacktest ; then
    :
  else
    echo "*** Could not run LAPACK test program, checking why..."
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
         LIBS="$ac_save_LIBS"
  fi
     LAPACK_LIBS=""
     ifelse([$3], , :, [$3])
fi
dnl finalize
LIBS="$ac_save_LIBS"
AC_SUBST(LAPACK_LIBS)
rm -f conf.lapacktest
])
