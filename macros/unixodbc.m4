# Configure paths for unixODBC
# Allin Cottrell <cottrell@wfu.edu>

dnl AM_PATH_ODBC([MINIMUM-VERSION, [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl Test for unixODBC, and define ODBC_CFLAGS and ODBC_LIBS.
dnl
AC_DEFUN([AM_PATH_ODBC],
[dnl 
AC_ARG_WITH(ODBC-prefix,[  --with-ODBC-prefix=PFX   Prefix where unixODBC is installed (optional)],
            ODBC_config_prefix="$withval", ODBC_config_prefix="")

  if test x$ODBC_config_prefix != x ; then
     ODBC_config_args="$ODBC_config_args --prefix=$ODBC_config_prefix"
  fi

  AC_MSG_CHECKING(for ODBC)
  if test x"${ODBC_LIBS}" = x ; then  
     ODBC_CFLAGS="-I$ODBC_config_prefix/include"
     ODBC_LIBS="-L$ODBC_config_prefix/lib -lodbc"
  fi

  ac_save_CFLAGS="$CFLAGS"
  ac_save_LIBS="$LIBS"
  CFLAGS="$ODBC_CFLAGS $CFLAGS"
  LIBS="$ODBC_LIBS $LIBS"

dnl
dnl Check the installed ODBC.
dnl
  rm -f conf.ODBCtest
  AC_TRY_RUN([
#include <stdlib.h>
#include <sql.h>
#include <sqlext.h>
#include <sqltypes.h>
int main (void)
{
  SQLFreeConnect(NULL);
  system ("touch conf.ODBCtest");
  return 0;
}
],, no_ODBC=yes,[echo $ac_n "cross compiling; assumed OK... $ac_c"])
  CFLAGS="$ac_save_CFLAGS"
  LIBS="$ac_save_LIBS"

  if test "x$no_ODBC" = x ; then
     AC_MSG_RESULT(yes)
     ifelse([$2], , :, [$2])     
  else
     AC_MSG_RESULT(no)
     if test -f conf.ODBCtest ; then
       :
     else
       echo "*** Could not run ODBC test program, checking why..."
       CFLAGS="$ODBC_CFLAGS $CFLAGS"
       LIBS="$LIBS $ODBC_LIBS"
       AC_TRY_LINK([
#include <stdlib.h>
#include <sql.h>
#include <sqlext.h>
#include <sqltypes.h>
],     [ return (1); ],
       [ echo "*** The test program compiled, but did not run. This usually means"
         echo "*** that the run-time linker is not finding ODBC. If it is not"
         echo "*** finding ODBC, you'll need to set your LD_LIBRARY_PATH "
         echo "*** environment variable, or edit /etc/ld.so.conf to point"
         echo "*** to the installed location.  Also, make sure you have run"
         echo "*** ldconfig if that is required on your system."
         echo "***" ],
       [ echo "*** The test program failed to compile or link. See config.log for the"
         echo "*** exact error that occured. This may mean ODBC was incorrectly installed"
         echo "*** or that you have moved ODBC since it was installed." ])
         CFLAGS="$ac_save_CFLAGS"
         LIBS="$ac_save_LIBS"
     fi
     ODBC_CFLAGS=""
     ODBC_LIBS=""
     ifelse([$3], , :, [$3])
  fi
  AC_SUBST(ODBC_CFLAGS)
  AC_SUBST(ODBC_LIBS)
  rm -f conf.ODBCtest
])
