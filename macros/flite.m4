# Configure paths for flite
# Allin Cottrell <cottrell@wfu.edu>, May 2004

dnl AM_PATH_FLITE([MINIMUM-VERSION, [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl Test for FLITE, and define FLITE_CFLAGS and FLITE_LIBS.
dnl
AC_DEFUN([AM_PATH_FLITE],
[dnl 
AC_ARG_WITH(flite-prefix,[  --with-flite-prefix=PFX   Prefix where FLITE is installed (optional)],
            flite_config_prefix="$withval", flite_config_prefix="")

  if test x$flite_config_prefix != x ; then
     flite_config_args="$flite_config_args --prefix=$flite_config_prefix"
  fi

  AC_MSG_CHECKING(for FLITE)
  if test x"${FLITE_LIBS}" = x ; then  
     FLITE_CFLAGS="-I$flite_config_prefix/include"
     FLITE_LIBS="-L$flite_config_prefix/lib -lflite_cmu_us_kal16 -lflite_usenglish -lflite_cmulex -lflite -lm"
  fi

  ac_save_CFLAGS="$CFLAGS"
  ac_save_LIBS="$LIBS"
  CFLAGS="$FLITE_CFLAGS $CFLAGS"
  LIBS="$FLITE_LIBS $LIBS"

dnl
dnl Check the installed FLITE.
dnl
  rm -f conf.flitetest
  AC_TRY_RUN([
#include <flite/flite.h>
#include <stdio.h>
#include <stdlib.h>
extern cst_voice *register_cmu_us_kal (void);
int main (void)
{
  flite_init();
  register_cmu_us_kal();
  system ("touch conf.flitetest");
  return 0;
}
],, no_flite=yes,[echo $ac_n "cross compiling; assumed OK... $ac_c"])
  CFLAGS="$ac_save_CFLAGS"
  LIBS="$ac_save_LIBS"

  if test "x$no_flite" = x ; then
     AC_MSG_RESULT(yes)
     ifelse([$2], , :, [$2])     
  else
     AC_MSG_RESULT(no)
     if test -f conf.flitetest ; then
       :
     else
       echo "*** Could not run FLITE test program, checking why..."
       CFLAGS="$FLITE_CFLAGS $CFLAGS"
       LIBS="$LIBS $FLITE_LIBS"
       AC_TRY_LINK([
#include <flite/flite.h>
#include <stdio.h>
#include <stdlib.h>
],     [ return (1); ],
       [ echo "*** The test program compiled, but did not run. This usually means"
         echo "*** that the run-time linker is not finding FLITE. If it is not"
         echo "*** finding FLITE, you'll need to set your LD_LIBRARY_PATH "
         echo "*** environment variable, or edit /etc/ld.so.conf to point"
         echo "*** to the installed location.  Also, make sure you have run"
         echo "*** ldconfig if that is required on your system."
         echo "***" ],
       [ echo "*** The test program failed to compile or link. See config.log for the"
         echo "*** exact error that occured. This may mean FLITE was incorrectly installed"
         echo "*** or that you have moved FLITE since it was installed." ])
         CFLAGS="$ac_save_CFLAGS"
         LIBS="$ac_save_LIBS"
     fi
     FLITE_CFLAGS=""
     FLITE_LIBS=""
     ifelse([$3], , :, [$3])
  fi
  AC_SUBST(FLITE_CFLAGS)
  AC_SUBST(FLITE_LIBS)
  rm -f conf.flitetest
])
