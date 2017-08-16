# Configure paths for libsvm
# Allin Cottrell, August 2017

# Based on macros by Owen Taylor.

dnl AM_PATH_LIBSVM([MINIMUM-VERSION, [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl Test for LIBSVM, and define LIBSVM_CFLAGS and LIBSVM_LIBS.
dnl
AC_DEFUN([AM_PATH_LIBSVM],
[dnl 
AC_ARG_WITH(libsvm-prefix,[  --with-libsvm-prefix=PFX   Prefix where LIBSVM is installed (optional)],
            libsvm_config_prefix="$withval", libsvm_config_prefix="")

  if test x$libsvm_config_prefix != x ; then
     libsvm_config_args="$libsvm_config_args --prefix=$libsvm_config_prefix"
  fi

  min_libsvm_version=ifelse([$1], ,321,$1)

  AC_MSG_CHECKING(for LIBSVM - version >= $min_libsvm_version)

  LIBSVM_CFLAGS="-I$libsvm_config_prefix/include"
  LIBSVM_LIBS="-L$libsvm_config_prefix/lib -lsvm"

  ac_save_CFLAGS="$CFLAGS"
  ac_save_LIBS="$LIBS"
  CFLAGS="$CFLAGS $LIBSVM_CFLAGS"
  LIBS="$LIBSVM_LIBS $LIBS"

dnl
dnl Now check if the installed LIBSVM is sufficiently new.
dnl
  rm -f conf.libsvmtest
  AC_TRY_RUN([
#include <svm.h>
#include <stdio.h>
#include <stdlib.h>

int 
main ()
{
  int my_libsvm_version = 0;
  const char *tmp = "$min_libsvm_version";
  int req = atoi(tmp);
  
  system ("touch conf.libsvmtest");

#ifdef LIBSVM_VERSION
  my_libsvm_version = LIBSVM_VERSION;
#endif

  if (my_libsvm_version >= req)
  {
    return 0;
  }
  else
  {
    printf("*** You need LIBSVM >= %d. The latest version is available\n", req);
    printf("*** from https://www.csie.ntu.edu.tw/~cjlin/libsvm/.\n");
    printf("***\n");
  }

  return 1;
}
],, no_libsvm=yes,[echo $ac_n "cross compiling; assumed OK... $ac_c"])
  CFLAGS="$ac_save_CFLAGS"
  LIBS="$ac_save_LIBS"

  if test "x$no_libsvm" = x ; then
     AC_MSG_RESULT(yes)
     ifelse([$2], , :, [$2])     
  else
     AC_MSG_RESULT(no)
     if test -f conf.libsvmtest ; then
       :
     else
       echo "*** Could not run LIBSVM test program, checking why..."
       CFLAGS="$CFLAGS $LIBSVM_CFLAGS"
       LIBS="$LIBS $LIBSVM_LIBS"
       AC_TRY_LINK([
#include <svm.h>
#include <stdio.h>
],     [ return (1); ],
       [ echo "*** The test program compiled, but did not run. This usually means"
         echo "*** that the run-time linker is not finding LIBSVM or finding the wrong"
         echo "*** version of LIBSVM. If it is not finding LIBSVM, you'll need to set your"
         echo "*** LD_LIBRARY_PATH environment variable, or edit /etc/ld.so.conf to point"
         echo "*** to the installed location  Also, make sure you have run ldconfig if that"
         echo "*** is required on your system"
         echo "***"
         echo "*** If you have an old version installed, it is best to remove it, although"
         echo "*** you may also be able to get things to work by modifying LD_LIBRARY_PATH"
         echo "***" ],
       [ echo "*** The test program failed to compile or link. See the file config.log for the"
         echo "*** exact error that occured. This usually means LIBSVM was incorrectly installed"
         echo "*** or that you have moved LIBSVM since it was installed." ])
         CFLAGS="$ac_save_CFLAGS"
         LIBS="$ac_save_LIBS"
     fi
     LIBSVM_CFLAGS=""
     LIBSVM_LIBS=""
     ifelse([$3], , :, [$3])
  fi
  if test "$LIBSVM_CFLAGS" = "-I/include" ; then 
     LIBSVM_CFLAGS=""
  fi
  if test "$LIBSVM_LIBS" = "-L/lib -lsvm" ; then 
     LIBSVM_LIBS="-lsvm"
  fi
  AC_SUBST(LIBSVM_CFLAGS)
  AC_SUBST(LIBSVM_LIBS)
  rm -f conf.libsvmtest
])
