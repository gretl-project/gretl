# Configure paths for MPFR
# Allin Cottrell, June 2006

# Based on macros by Owen Taylor.

dnl AM_PATH_MPFR([MINIMUM-VERSION, [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl Test for MPFR, and define MPFR_CFLAGS and MPFR_LIBS.
dnl
AC_DEFUN([AM_PATH_MPFR],
[dnl 
AC_ARG_WITH(mpfr-prefix,[  --with-mpfr-prefix=PFX   Prefix where MPFR is installed (optional)],
            mpfr_config_prefix="$withval", mpfr_config_prefix="")

  if test x$mpfr_config_prefix != x ; then
     mpfr_config_args="$mpfr_config_args --prefix=$mpfr_config_prefix"
  fi

  min_mpfr_version=ifelse([$1], ,1.0.0,$1)

  AC_MSG_CHECKING(for MPFR - version >= $min_mpfr_version)

  MPFR_CFLAGS="-I$mpfr_config_prefix/include"
  MPFR_LIBS="-L$mpfr_config_prefix/lib -lmpfr"

  ac_save_CFLAGS="$CFLAGS"
  ac_save_LIBS="$LIBS"
  CFLAGS="$CFLAGS $MPFR_CFLAGS"
  LIBS="$MPFR_LIBS $LIBS"

dnl
dnl Now check if the installed MPFR is sufficiently new.
dnl
  rm -f conf.mpfrtest
  AC_TRY_RUN([
#include <mpfr.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int 
main ()
{
  int mpfr_major_version = 0, mpfr_minor_version = 0, mpfr_micro_version = 0;
  int major, minor, micro;
  char *tmp_version;
  mpfr_t a, b, c;

  mpfr_init (a);
  mpfr_init (b);
  mpfr_init (c);
  mpfr_mul (c, a, b, GMP_RNDN);

  system ("touch conf.mpfrtest");

#ifdef MPFR_VERSION_MAJOR
  mpfr_major_version = MPFR_VERSION_MAJOR;
#endif

#ifdef MPFR_VERSION_MINOR
  mpfr_minor_version = MPFR_VERSION_MINOR;
#endif

#ifdef MPFR_VERSION_PATCHLEVEL
  mpfr_micro_version = MPFR_VERSION_PATCHLEVEL;
#endif

  /* HP/UX 9 (%@#!) writes to sscanf strings */
  tmp_version = strdup("$min_mpfr_version");
  if (sscanf(tmp_version, "%d.%d.%d", &major, &minor, &micro) != 3) {
     printf("%s, bad version string\n", "$min_mpfr_version");
     exit(1);
  }

  if ((mpfr_major_version > major) ||
     ((mpfr_major_version == major) && (mpfr_minor_version > minor)) ||
     ((mpfr_major_version == major) && (mpfr_minor_version == minor) && (mpfr_micro_version >= micro)))
  {
    return 0;
  }
  else
  {
    printf("\n*** An old version of MPFR (%d.%d.%d) was found.\n",
           mpfr_major_version, mpfr_minor_version, mpfr_micro_version);
    printf("*** You need a version of MPFR newer than %d.%d.%d. The latest version of\n",
           major, minor, micro);

    printf("*** MPFR is available from http://www.mpfr.org/.\n");
    printf("***\n");
  }

  return 1;
}
],, no_mpfr=yes,[echo $ac_n "cross compiling; assumed OK... $ac_c"])
  CFLAGS="$ac_save_CFLAGS"
  LIBS="$ac_save_LIBS"

  if test "x$no_mpfr" = x ; then
     AC_MSG_RESULT(yes)
     ifelse([$2], , :, [$2])     
  else
     AC_MSG_RESULT(no)
     if test -f conf.mpfrtest ; then
       :
     else
       echo "*** Could not run MPFR test program, checking why..."
       CFLAGS="$CFLAGS $MPFR_CFLAGS"
       LIBS="$LIBS $MPFR_LIBS"
       AC_TRY_LINK([
#include <mpfr.h>
#include <stdio.h>
],     [ return (1); ],
       [ echo "*** The test program compiled, but did not run. This usually means"
         echo "*** that the run-time linker is not finding MPFR or finding the wrong"
         echo "*** version of MPFR. If it is not finding MPFR, you'll need to set your"
         echo "*** LD_LIBRARY_PATH environment variable, or edit /etc/ld.so.conf to point"
         echo "*** to the installed location  Also, make sure you have run ldconfig if that"
         echo "*** is required on your system"
         echo "***"
         echo "*** If you have an old version installed, it is best to remove it, although"
         echo "*** you may also be able to get things to work by modifying LD_LIBRARY_PATH"
         echo "***" ],
       [ echo "*** The test program failed to compile or link. See the file config.log for the"
         echo "*** exact error that occured. This usually means MPFR was incorrectly installed"
         echo "*** or that you have moved MPFR since it was installed." ])
         CFLAGS="$ac_save_CFLAGS"
         LIBS="$ac_save_LIBS"
     fi
     MPFR_CFLAGS=""
     MPFR_LIBS=""
     ifelse([$3], , :, [$3])
  fi
  if test "$MPFR_CFLAGS" = "-I/include" ; then 
     MPFR_CFLAGS=""
  fi
  if test "$MPFR_LIBS" = "-L/lib -lmpfr" ; then 
     MPFR_LIBS="-lmpfr"
  fi
  AC_SUBST(MPFR_CFLAGS)
  AC_SUBST(MPFR_LIBS)
  rm -f conf.mpfrtest
])
