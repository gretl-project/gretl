# Configure paths for GRETL2 (broken!!)
# Allin Cottrell July, 2002
# Based on macros by Owen Taylor

dnl AM_PATH_GRETL([MINIMUM-VERSION, [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl Test for GRETL, and define GRETL_CFLAGS and GRETL_LIBS
dnl
AC_DEFUN(AM_PATH_GRETL,
[dnl 
dnl Get the cflags and libraries from pkgconfig
dnl
AC_ARG_ENABLE(gretltest, [  --disable-gretltest       Do not try to compile and run a test GRETL program],
		    , enable_gretltest=yes)

  AC_PATH_PROG(GRETL_CONFIG, gretl-config, no)
  min_gretl_version=ifelse([$1], ,15.0.0,$1)
  AC_MSG_CHECKING(for libgretl - version >= $min_gretl_version)
  no_gretl=""
  if test "$GRETL_CONFIG" = "no" ; then
    no_gretl=yes
  else
    GRETL_CFLAGS=`pkg-config --cflags gretl2`
    GRETL_LIBS=`pkg-config --libs gretl2`

    gretl_major_version=`pkg-config --modversion gretl2 | \
           sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\1/'`
    gretl_minor_version=`pkg-config --modversion gretl2 | \
           sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\2/'`
    gretl_micro_version=`pkg-config --modversion gretl2 | \
           sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\3/'`
    if test "x$enable_gretltest" = "xyes" ; then
      ac_save_CFLAGS="$CFLAGS"
      ac_save_LIBS="$LIBS"
      CFLAGS="$CFLAGS $GRETL_CFLAGS"
      LIBS="$LIBS $GRETL_LIBS"
dnl
dnl Now check if the installed GRETL is sufficiently new. (Also sanity
dnl checks the results of gretl-config to some extent)
dnl
      rm -f conf.gretltest
      AC_TRY_RUN([
#include <gretl/libgretl.h>

char*
my_strdup (char *str)
{
  char *new_str;
  
  if (str)
    {
      new_str = malloc ((strlen (str) + 1));
      strcpy (new_str, str);
    }
  else
    new_str = NULL;
  
  return new_str;
}

int main ()
{
  int major, minor, micro;
  char *tmp_version;

  system ("touch conf.gretltest");

  /* HP/UX 9 (%@#!) writes to sscanf strings */
  tmp_version = my_strdup("$min_gretl_version");
  if (sscanf(tmp_version, "%d.%d.%d", &major, &minor, &micro) != 3) {
     printf("%s, bad version string\n", "$min_gretl_version");
     exit(1);
   }

   if (($gretl_major_version > major) ||
      (($gretl_major_version == major) && ($gretl_minor_version > minor)) ||
      (($gretl_major_version == major) && ($gretl_minor_version == minor) && ($gretl_micro_version >= micro)))
    {
      return 0;
    }
  else
    {
      printf("\n*** 'pkg-config --modversion gretl2' returned %d.%d.%d, but the\n", $gretl_major_version, $gretl_minor_version, $gretl_micro_version);
      printf("*** minimum version of GRETL required is %d.%d.%d.\n", major, minor, micro);
      return 1;
    }
}

],, no_gretl=yes,[echo $ac_n "cross compiling; assumed OK... $ac_c"])
       CFLAGS="$ac_save_CFLAGS"
       LIBS="$ac_save_LIBS"
     fi
  fi
  if test "x$no_gretl" = x ; then
     AC_MSG_RESULT(yes)
     ifelse([$2], , :, [$2])     
  else
     AC_MSG_RESULT(no)
     if test "$GRETL_CONFIG" = "no" ; then
       echo "*** The gretl-config script installed by GRETL could not be found"
       echo "*** If GRETL was installed in PREFIX, make sure PREFIX/bin is in"
       echo "*** your path, or set the GRETL_CONFIG environment variable to the"
       echo "*** full path to gretl-config."
     else
       if test -f conf.gretltest ; then
        :
       else
          echo "*** Could not run GRETL test program, checking why..."
          CFLAGS="$CFLAGS $GRETL_CFLAGS"
          LIBS="$LIBS $GRETL_LIBS"
          AC_TRY_LINK([
#include <stdio.h>
#include <gretl/libgretl.h>
],      [ return 0; ],
        [ echo "*** The test program compiled, but did not run. This usually means"
          echo "*** that the run-time linker is not finding GRETL or finding the wrong"
          echo "*** version of GRETL. If it is not finding GRETL, you'll need to set your"
          echo "*** LD_LIBRARY_PATH environment variable, or edit /etc/ld.so.conf to point"
          echo "*** to the installed location.  Also, make sure you have run ldconfig if that"
          echo "*** is required on your system"
	  echo "***"
          echo "*** If you have an old version installed, it is best to remove it, although"
          echo "*** you may also be able to get things to work by modifying LD_LIBRARY_PATH"],
        [ echo "*** The test program failed to compile or link. See the file config.log for the"
          echo "*** exact error that occured. This usually means GRETL was incorrectly installed"
          echo "*** or that you have moved GRETL since it was installed. In the latter case, you"
          echo "*** may want to edit the gretl-config script: $GRETL_CONFIG" ])
          CFLAGS="$ac_save_CFLAGS"
          LIBS="$ac_save_LIBS"
       fi
     fi
     GRETL_CFLAGS=""
     GRETL_LIBS=""
     ifelse([$3], , :, [$3])
  fi
  AC_SUBST(GRETL_CFLAGS)
  AC_SUBST(GRETL_LIBS)
  rm -f conf.gretltest
])

