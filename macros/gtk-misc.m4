# Configure paths for GTK+
# Owen Taylor     97-11-3

dnl AM_PATH_GTK([MINIMUM-VERSION, [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND [, MODULES]]]])
dnl Test for GTK, and define GTK_CFLAGS and GTK_LIBS
dnl
AC_DEFUN(AM_PATH_GTK,
[dnl 
dnl Get the cflags and libraries from the gtk-config script
dnl
AC_ARG_WITH(gtk-prefix,[  --with-gtk-prefix=PFX   Prefix where GTK is installed (optional)],
            gtk_config_prefix="$withval", gtk_config_prefix="")
AC_ARG_WITH(gtk-exec-prefix,[  --with-gtk-exec-prefix=PFX Exec prefix where GTK is installed (optional)],
            gtk_config_exec_prefix="$withval", gtk_config_exec_prefix="")
AC_ARG_ENABLE(gtktest, [  --disable-gtktest       Do not try to compile and run a test GTK program],
		    , enable_gtktest=yes)

  for module in . $4
  do
      case "$module" in
         gthread) 
             gtk_config_args="$gtk_config_args gthread"
         ;;
      esac
  done

  if test x$gtk_config_exec_prefix != x ; then
     gtk_config_args="$gtk_config_args --exec-prefix=$gtk_config_exec_prefix"
     if test x${GTK_CONFIG+set} != xset ; then
        GTK_CONFIG=$gtk_config_exec_prefix/bin/gtk-config
     fi
  fi
  if test x$gtk_config_prefix != x ; then
     gtk_config_args="$gtk_config_args --prefix=$gtk_config_prefix"
     if test x${GTK_CONFIG+set} != xset ; then
        GTK_CONFIG=$gtk_config_prefix/bin/gtk-config
     fi
  fi

  AC_PATH_PROG(GTK_CONFIG, gtk-config, no)
  min_gtk_version=ifelse([$1], ,0.99.7,$1)
  AC_MSG_CHECKING(for GTK - version >= $min_gtk_version)
  no_gtk=""
  if test "$GTK_CONFIG" = "no" ; then
    no_gtk=yes
  else
    GTK_CFLAGS=`$GTK_CONFIG $gtk_config_args --cflags`
    GTK_LIBS=`$GTK_CONFIG $gtk_config_args --libs`
    gtk_config_major_version=`$GTK_CONFIG $gtk_config_args --version | \
           sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\1/'`
    gtk_config_minor_version=`$GTK_CONFIG $gtk_config_args --version | \
           sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\2/'`
    gtk_config_micro_version=`$GTK_CONFIG $gtk_config_args --version | \
           sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\3/'`
    if test "x$enable_gtktest" = "xyes" ; then
      ac_save_CFLAGS="$CFLAGS"
      ac_save_LIBS="$LIBS"
      CFLAGS="$CFLAGS $GTK_CFLAGS"
      LIBS="$GTK_LIBS $LIBS"
dnl
dnl Now check if the installed GTK is sufficiently new. (Also sanity
dnl checks the results of gtk-config to some extent
dnl
      rm -f conf.gtktest
      AC_TRY_RUN([
#include <gtk/gtk.h>
#include <stdio.h>
#include <stdlib.h>

int 
main ()
{
  int major, minor, micro;
  char *tmp_version;

  system ("touch conf.gtktest");

  /* HP/UX 9 (%@#!) writes to sscanf strings */
  tmp_version = g_strdup("$min_gtk_version");
  if (sscanf(tmp_version, "%d.%d.%d", &major, &minor, &micro) != 3) {
     printf("%s, bad version string\n", "$min_gtk_version");
     exit(1);
   }

  if ((gtk_major_version != $gtk_config_major_version) ||
      (gtk_minor_version != $gtk_config_minor_version) ||
      (gtk_micro_version != $gtk_config_micro_version))
    {
      printf("\n*** 'gtk-config --version' returned %d.%d.%d, but GTK+ (%d.%d.%d)\n", 
             $gtk_config_major_version, $gtk_config_minor_version, $gtk_config_micro_version,
             gtk_major_version, gtk_minor_version, gtk_micro_version);
      printf ("*** was found! If gtk-config was correct, then it is best\n");
      printf ("*** to remove the old version of GTK+. You may also be able to fix the error\n");
      printf("*** by modifying your LD_LIBRARY_PATH enviroment variable, or by editing\n");
      printf("*** /etc/ld.so.conf. Make sure you have run ldconfig if that is\n");
      printf("*** required on your system.\n");
      printf("*** If gtk-config was wrong, set the environment variable GTK_CONFIG\n");
      printf("*** to point to the correct copy of gtk-config, and remove the file config.cache\n");
      printf("*** before re-running configure\n");
    } 
#if defined (GTK_MAJOR_VERSION) && defined (GTK_MINOR_VERSION) && defined (GTK_MICRO_VERSION)
  else if ((gtk_major_version != GTK_MAJOR_VERSION) ||
	   (gtk_minor_version != GTK_MINOR_VERSION) ||
           (gtk_micro_version != GTK_MICRO_VERSION))
    {
      printf("*** GTK+ header files (version %d.%d.%d) do not match\n",
	     GTK_MAJOR_VERSION, GTK_MINOR_VERSION, GTK_MICRO_VERSION);
      printf("*** library (version %d.%d.%d)\n",
	     gtk_major_version, gtk_minor_version, gtk_micro_version);
    }
#endif /* defined (GTK_MAJOR_VERSION) ... */
  else
    {
      if ((gtk_major_version > major) ||
        ((gtk_major_version == major) && (gtk_minor_version > minor)) ||
        ((gtk_major_version == major) && (gtk_minor_version == minor) && (gtk_micro_version >= micro)))
      {
        return 0;
       }
     else
      {
        printf("\n*** An old version of GTK+ (%d.%d.%d) was found.\n",
               gtk_major_version, gtk_minor_version, gtk_micro_version);
        printf("*** You need a version of GTK+ newer than %d.%d.%d. The latest version of\n",
	       major, minor, micro);
        printf("*** GTK+ is always available from ftp://ftp.gtk.org.\n");
        printf("***\n");
        printf("*** If you have already installed a sufficiently new version, this error\n");
        printf("*** probably means that the wrong copy of the gtk-config shell script is\n");
        printf("*** being found. The easiest way to fix this is to remove the old version\n");
        printf("*** of GTK+, but you can also set the GTK_CONFIG environment to point to the\n");
        printf("*** correct copy of gtk-config. (In this case, you will have to\n");
        printf("*** modify your LD_LIBRARY_PATH enviroment variable, or edit /etc/ld.so.conf\n");
        printf("*** so that the correct libraries are found at run-time))\n");
      }
    }
  return 1;
}
],, no_gtk=yes,[echo $ac_n "cross compiling; assumed OK... $ac_c"])
       CFLAGS="$ac_save_CFLAGS"
       LIBS="$ac_save_LIBS"
     fi
  fi
  if test "x$no_gtk" = x ; then
     AC_MSG_RESULT(yes)
     ifelse([$2], , :, [$2])     
  else
     AC_MSG_RESULT(no)
     if test "$GTK_CONFIG" = "no" ; then
       echo "*** The gtk-config script installed by GTK could not be found"
       echo "*** If GTK was installed in PREFIX, make sure PREFIX/bin is in"
       echo "*** your path, or set the GTK_CONFIG environment variable to the"
       echo "*** full path to gtk-config."
     else
       if test -f conf.gtktest ; then
        :
       else
          echo "*** Could not run GTK test program, checking why..."
          CFLAGS="$CFLAGS $GTK_CFLAGS"
          LIBS="$LIBS $GTK_LIBS"
          AC_TRY_LINK([
#include <gtk/gtk.h>
#include <stdio.h>
],      [ return ((gtk_major_version) || (gtk_minor_version) || (gtk_micro_version)); ],
        [ echo "*** The test program compiled, but did not run. This usually means"
          echo "*** that the run-time linker is not finding GTK or finding the wrong"
          echo "*** version of GTK. If it is not finding GTK, you'll need to set your"
          echo "*** LD_LIBRARY_PATH environment variable, or edit /etc/ld.so.conf to point"
          echo "*** to the installed location  Also, make sure you have run ldconfig if that"
          echo "*** is required on your system"
	  echo "***"
          echo "*** If you have an old version installed, it is best to remove it, although"
          echo "*** you may also be able to get things to work by modifying LD_LIBRARY_PATH"
          echo "***"
          echo "*** If you have a RedHat 5.0 system, you should remove the GTK package that"
          echo "*** came with the system with the command"
          echo "***"
          echo "***    rpm --erase --nodeps gtk gtk-devel" ],
        [ echo "*** The test program failed to compile or link. See the file config.log for the"
          echo "*** exact error that occured. This usually means GTK was incorrectly installed"
          echo "*** or that you have moved GTK since it was installed. In the latter case, you"
          echo "*** may want to edit the gtk-config script: $GTK_CONFIG" ])
          CFLAGS="$ac_save_CFLAGS"
          LIBS="$ac_save_LIBS"
       fi
     fi
     GTK_CFLAGS=""
     GTK_LIBS=""
     ifelse([$3], , :, [$3])
  fi
  AC_SUBST(GTK_CFLAGS)
  AC_SUBST(GTK_LIBS)
  rm -f conf.gtktest
])

# Configure paths for GTK+EXTRA
# Owen Taylor     97-11-3
# Adrian Feiguin  01-04-03 

dnl AM_PATH_GTK_EXTRA([MINIMUM-VERSION, [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND [, MODULES]]]])
dnl Test for GTK_EXTRA, and define GTK_EXTRA_CFLAGS and GTK_EXTRA_LIBS
dnl
AC_DEFUN(AM_PATH_GTK_EXTRA,
[dnl 
dnl Get the cflags and libraries from the gtkextra-config script
dnl
AC_ARG_WITH(gtkextra-prefix,[  --with-gtkextra-prefix=PFX   Prefix where GTK_EXTRA is installed (optional)],
            gtkextra_config_prefix="$withval", gtkextra_config_prefix="")
AC_ARG_WITH(gtkextra-exec-prefix,[  --with-gtkextra-exec-prefix=PFX Exec prefix where GTK_EXTRA is installed (optional)],
            gtkextra_config_exec_prefix="$withval", gtkextra_config_exec_prefix="")
AC_ARG_ENABLE(gtkextratest, [  --disable-gtkextratest       Do not try to compile and run a test GTK_EXTRA program],
		    , enable_gtkextratest=yes)

  for module in . $4
  do
      case "$module" in
         gthread) 
             gtkextra_config_args="$gtkextra_config_args gthread"
         ;;
      esac
  done

  if test x$gtkextra_config_exec_prefix != x ; then
     gtkextra_config_args="$gtkextra_config_args --exec-prefix=$gtkextra_config_exec_prefix"
     if test x${GTK_EXTRA_CONFIG+set} != xset ; then
        GTK_EXTRA_CONFIG=$gtkextra_config_exec_prefix/bin/gtkextra-config
     fi
  fi
  if test x$gtkextra_config_prefix != x ; then
     gtkextra_config_args="$gtkextra_config_args --prefix=$gtkextra_config_prefix"
     if test x${GTK_EXTRA_CONFIG+set} != xset ; then
        GTK_EXTRA_CONFIG=$gtkextra_config_prefix/bin/gtkextra-config
     fi
  fi

  AC_PATH_PROG(GTK_EXTRA_CONFIG, gtkextra-config, no)
  min_gtkextra_version=ifelse([$1], ,0.99.13,$1)
  AC_MSG_CHECKING(for GTK_EXTRA - version >= $min_gtkextra_version)
  no_gtkextra=""
  if test "$GTK_EXTRA_CONFIG" = "no" ; then
    no_gtkextra=yes
  else
    GTK_EXTRA_CFLAGS=`$GTK_EXTRA_CONFIG $gtkextra_config_args --cflags`
    GTK_EXTRA_LIBS=`$GTK_EXTRA_CONFIG $gtkextra_config_args --libs`
    gtkextra_config_major_version=`$GTK_EXTRA_CONFIG $gtkextra_config_args --version | \
           sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\1/'`
    gtkextra_config_minor_version=`$GTK_EXTRA_CONFIG $gtkextra_config_args --version | \
           sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\2/'`
    gtkextra_config_micro_version=`$GTK_EXTRA_CONFIG $gtkextra_config_args --version | \
           sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\3/'`
    if test "x$enable_gtkextratest" = "xyes" ; then
      ac_save_CFLAGS="$CFLAGS"
      ac_save_LIBS="$LIBS"
      CFLAGS="$CFLAGS $GTK_EXTRA_CFLAGS"
      LIBS="$GTK_EXTRA_LIBS $LIBS"
dnl
dnl Now check if the installed GTK_EXTRA is sufficiently new. (Also sanity
dnl checks the results of gtkextra-config to some extent
dnl
      rm -f conf.gtkextratest
      AC_TRY_RUN([
#include <gtkextra/gtkextra.h>
#include <stdio.h>
#include <stdlib.h>

int 
main ()
{
  int major, minor, micro;
  char *tmp_version;

  system ("touch conf.gtkextratest");

  /* HP/UX 9 (%@#!) writes to sscanf strings */
  tmp_version = g_strdup("$min_gtkextra_version");
  if (sscanf(tmp_version, "%d.%d.%d", &major, &minor, &micro) != 3) {
     printf("%s, bad version string\n", "$min_gtkextra_version");
     exit(1);
   }

  if ((gtkextra_major_version != $gtkextra_config_major_version) ||
      (gtkextra_minor_version != $gtkextra_config_minor_version) ||
      (gtkextra_micro_version != $gtkextra_config_micro_version))
    {
      printf("\n*** 'gtkextra-config --version' returned %d.%d.%d, but GTK_EXTRA+ (%d.%d.%d)\n", 
             $gtkextra_config_major_version, $gtkextra_config_minor_version, $gtkextra_config_micro_version,
             gtkextra_major_version, gtkextra_minor_version, gtkextra_micro_version);
      printf ("*** was found! If gtkextra-config was correct, then it is best\n");
      printf ("*** to remove the old version of GTK_EXTRA+. You may also be able to fix the error\n");
      printf("*** by modifying your LD_LIBRARY_PATH enviroment variable, or by editing\n");
      printf("*** /etc/ld.so.conf. Make sure you have run ldconfig if that is\n");
      printf("*** required on your system.\n");
      printf("*** If gtkextra-config was wrong, set the environment variable GTK_EXTRA_CONFIG\n");
      printf("*** to point to the correct copy of gtkextra-config, and remove the file config.cache\n");
      printf("*** before re-running configure\n");
    } 
#if defined (GTK_EXTRA_MAJOR_VERSION) && defined (GTK_EXTRA_MINOR_VERSION) && defined (GTK_EXTRA_MICRO_VERSION)
  else if ((gtkextra_major_version != GTK_EXTRA_MAJOR_VERSION) ||
	   (gtkextra_minor_version != GTK_EXTRA_MINOR_VERSION) ||
           (gtkextra_micro_version != GTK_EXTRA_MICRO_VERSION))
    {
      printf("*** GTK_EXTRA+ header files (version %d.%d.%d) do not match\n",
	     GTK_EXTRA_MAJOR_VERSION, GTK_EXTRA_MINOR_VERSION, GTK_EXTRA_MICRO_VERSION);
      printf("*** library (version %d.%d.%d)\n",
	     gtkextra_major_version, gtkextra_minor_version, gtkextra_micro_version);
    }
#endif /* defined (GTK_EXTRA_MAJOR_VERSION) ... */
  else
    {
      if ((gtkextra_major_version > major) ||
        ((gtkextra_major_version == major) && (gtkextra_minor_version > minor)) ||
        ((gtkextra_major_version == major) && (gtkextra_minor_version == minor) && (gtkextra_micro_version >= micro)))
      {
        return 0;
       }
     else
      {
        printf("\n*** An old version of GTK_EXTRA+ (%d.%d.%d) was found.\n",
               gtkextra_major_version, gtkextra_minor_version, gtkextra_micro_version);
        printf("*** You need a version of GTK_EXTRA+ newer than %d.%d.%d. The latest version of\n",
	       major, minor, micro);
        printf("*** GTK_EXTRA+ is always available from ftp://ftp.gtkextra.org.\n");
        printf("***\n");
        printf("*** If you have already installed a sufficiently new version, this error\n");
        printf("*** probably means that the wrong copy of the gtkextra-config shell script is\n");
        printf("*** being found. The easiest way to fix this is to remove the old version\n");
        printf("*** of GTK_EXTRA+, but you can also set the GTK_EXTRA_CONFIG environment to point to the\n");
        printf("*** correct copy of gtkextra-config. (In this case, you will have to\n");
        printf("*** modify your LD_LIBRARY_PATH enviroment variable, or edit /etc/ld.so.conf\n");
        printf("*** so that the correct libraries are found at run-time))\n");
      }
    }
  return 1;
}
],, no_gtkextra=yes,[echo $ac_n "cross compiling; assumed OK... $ac_c"])
       CFLAGS="$ac_save_CFLAGS"
       LIBS="$ac_save_LIBS"
     fi
  fi
  if test "x$no_gtkextra" = x ; then
     AC_MSG_RESULT(yes)
     ifelse([$2], , :, [$2])     
  else
     AC_MSG_RESULT(no)
     if test "$GTK_EXTRA_CONFIG" = "no" ; then
       echo "*** The gtkextra-config script installed by GTK_EXTRA could not be found"
       echo "*** If GTK_EXTRA was installed in PREFIX, make sure PREFIX/bin is in"
       echo "*** your path, or set the GTK_EXTRA_CONFIG environment variable to the"
       echo "*** full path to gtkextra-config."
     else
       if test -f conf.gtkextratest ; then
        :
       else
          echo "*** Could not run GTK_EXTRA test program, checking why..."
          CFLAGS="$CFLAGS $GTK_EXTRA_CFLAGS"
          LIBS="$LIBS $GTK_EXTRA_LIBS"
          AC_TRY_LINK([
#include <gtkextra/gtkextra.h>
#include <stdio.h>
],      [ return ((gtkextra_major_version) || (gtkextra_minor_version) || (gtkextra_micro_version)); ],
        [ echo "*** The test program compiled, but did not run. This usually means"
          echo "*** that the run-time linker is not finding GTK_EXTRA or finding the wrong"
          echo "*** version of GTK_EXTRA. If it is not finding GTK_EXTRA, you'll need to set your"
          echo "*** LD_LIBRARY_PATH environment variable, or edit /etc/ld.so.conf to point"
          echo "*** to the installed location  Also, make sure you have run ldconfig if that"
          echo "*** is required on your system"
	  echo "***"
          echo "*** If you have an old version installed, it is best to remove it, although"
          echo "*** you may also be able to get things to work by modifying LD_LIBRARY_PATH"
          echo "***"
          echo "*** If you have a RedHat 5.0 system, you should remove the GTK_EXTRA package that"
          echo "*** came with the system with the command"
          echo "***"
          echo "***    rpm --erase --nodeps gtkextra gtkextra-devel" ],
        [ echo "*** The test program failed to compile or link. See the file config.log for the"
          echo "*** exact error that occured. This usually means GTK_EXTRA was incorrectly installed"
          echo "*** or that you have moved GTK_EXTRA since it was installed. In the latter case, you"
          echo "*** may want to edit the gtkextra-config script: $GTK_EXTRA_CONFIG" ])
          CFLAGS="$ac_save_CFLAGS"
          LIBS="$ac_save_LIBS"
       fi
     fi
     GTK_EXTRA_CFLAGS=""
     GTK_EXTRA_LIBS=""
     ifelse([$3], , :, [$3])
  fi
  AC_SUBST(GTK_EXTRA_CFLAGS)
  AC_SUBST(GTK_EXTRA_LIBS)
  rm -f conf.gtkextratest
])

# Configure paths for libole2
# Arturo Tena (tenix)	  1999-09-24
#
# This file is a modified version of the
# file glib.m4 that came with glib 1.2.3:
# Configure paths for GLIB
# Owen Taylor     97-11-3

# This m4 macro don't depend on GNOME, just glib.

dnl AM_PATH_LIBOLE2([MINIMUM-VERSION, [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl Test for libole2, and define LIBOLE2_CFLAGS and LIBOLE2_LIBS
dnl
dnl Example of use (write the following in your 'configure.in'):
dnl AM_PATH_LIBOLE2(0.0.1, [LIBS="$LIBS $LIBOLE2_LIBS" CFLAGS="$CFLAGS $LIBOLE2_CFLAGS"], AC_MSG_ERROR([Can't find libole2.]))
dnl
dnl
AC_DEFUN(AM_PATH_LIBOLE2,
[dnl 
dnl Get the cflags and libraries from the libole2-config script
dnl
AC_ARG_WITH(libole2-prefix,[  --with-libole2-prefix=PFX   Prefix where libole2 is installed (optional)],
           libole2_config_prefix="$withval", libole2_config_prefix="")
AC_ARG_WITH(libole2-exec-prefix,[  --with-libole2-exec-prefix=PFX Exec prefix where libole2 is installed (optional)],
           libole2_config_exec_prefix="$withval", libole2_config_exec_prefix="")
AC_ARG_ENABLE(libole2test, [  --disable-libole2test       Do not try to compile and run a test libole2 program],
		    , enable_libole2test=yes)

  if test x$libole2_config_exec_prefix != x ; then
     libole2_config_args="$libole2_config_args --exec-prefix=$libole2_config_exec_prefix"
     if test x${LIBOLE2_CONFIG+set} != xset ; then
        LIBOLE2_CONFIG=$libole2_config_exec_prefix/bin/libole2-config
     fi
  fi
  if test x$libole2_config_prefix != x ; then
     libole2_config_args="$libole2_config_args --prefix=$libole2_config_prefix"
     if test x${LIBOLE2_CONFIG+set} != xset ; then
        LIBOLE2_CONFIG=$libole2_config_prefix/bin/libole2-config
     fi
  fi

dnl  for module in . $4
dnl  do
dnl      case "$module" in
dnl         gmodule) 
dnl             glib_config_args="$glib_config_args gmodule"
dnl         ;;
dnl         gthread) 
dnl             glib_config_args="$glib_config_args gthread"
dnl         ;;
dnl      esac
dnl  done

  AC_PATH_PROG(LIBOLE2_CONFIG, libole2-config, no)
  min_libole2_version=ifelse([$1], ,0.0.1,$1)
  AC_MSG_CHECKING(for libole2 - version >= $min_libole2_version)
  no_libole2=""
  if test "$LIBOLE2_CONFIG" = "no" ; then
    no_libole2=yes
  else
    LIBOLE2_CFLAGS="`$LIBOLE2_CONFIG $libole2_config_args --cflags` `glib-config --cflags`"
    LIBOLE2_LIBS=`$LIBOLE2_CONFIG $libole2_config_args --libs`
    libole2_config_major_version=`$LIBOLE2_CONFIG $libole2_config_args --version | \
           sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\1/'`
    libole2_config_minor_version=`$LIBOLE2_CONFIG $libole2_config_args --version | \
           sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\2/'`
    libole2_config_micro_version=`$LIBOLE2_CONFIG $libole2_config_args --version | \
           sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\3/'`

    if test "x$enable_libole2test" = "xyes" ; then
      ac_save_CFLAGS="$CFLAGS"
      ac_save_LIBS="$LIBS"
      CFLAGS="$CFLAGS $LIBOLE2_CFLAGS"
      LIBS="$LIBOLE2_LIBS $LIBS"
dnl
dnl Now check if the installed libole2 is sufficiently new. (Also sanity
dnl checks the results of libole2-config to some extent)
dnl
      rm -f conf.libole2test
      AC_TRY_RUN([
#include <libole2/libole2.h>
#include <stdio.h>
#include <stdlib.h>

int
mystrlen (char * s)
{
	char *p;
	for (p = s; *p; p++) ;
	return p - s;
}

char *
mystrdup (char * s)
{
	char *ret, *p1, *p2;

	ret = malloc (mystrlen (s));
	if (ret == NULL) return NULL;

	p1 = s; p2 = ret;
	while (*p1) {
		*p2 = *p1;
		p1++;
		p2++;
	}

	return ret;
}

int 
main ()
{
  int major, minor, micro;
  char *tmp_version;

  system ("touch conf.libole2test");

  /* HP/UX 9 (%@#!) writes to sscanf strings */
  tmp_version = mystrdup("$min_libole2_version");
  if (sscanf(tmp_version, "%d.%d.%d", &major, &minor, &micro) != 3) {
     printf("%s, bad version string\n", "$min_libole2_version");
     exit(1);
   }

  if ((libole2_major_version != $libole2_config_major_version) ||
      (libole2_minor_version != $libole2_config_minor_version) ||
      (libole2_micro_version != $libole2_config_micro_version))
    {
      printf("\n*** 'libole2-config --version' returned %d.%d.%d, but libole2 (%d.%d.%d)\n", 
             $libole2_config_major_version, $libole2_config_minor_version, $libole2_config_micro_version,
             libole2_major_version, libole2_minor_version, libole2_micro_version);
      printf ("*** was found! If libole2-config was correct, then it is best\n");
      printf ("*** to remove the old version of libole2. You may also be able to fix the error\n");
      printf("*** by modifying your LD_LIBRARY_PATH enviroment variable, or by editing\n");
      printf("*** /etc/ld.so.conf. Make sure you have run ldconfig if that is\n");
      printf("*** required on your system.\n");
      printf("*** If libole2-config was wrong, set the environment variable LIBOLE2_CONFIG\n");
      printf("*** to point to the correct copy of libole2-config, and remove the file\n");
      printf("*** config.cache before re-running configure\n");
    } 
  else if ((libole2_major_version != LIBOLE2_MAJOR_VERSION) ||
	   (libole2_minor_version != LIBOLE2_MINOR_VERSION) ||
           (libole2_micro_version != LIBOLE2_MICRO_VERSION))
    {
      printf("*** LIBOLE2 header files (version %d.%d.%d) do not match\n",
	   LIBOLE2_MAJOR_VERSION, LIBOLE2_MINOR_VERSION, LIBOLE2_MICRO_VERSION);
      printf("*** library (version %d.%d.%d)\n",
	   libole2_major_version, libole2_minor_version, libole2_micro_version);
    }
  else
    {
      if ((libole2_major_version > major) ||
        ((libole2_major_version == major) && (libole2_minor_version > minor)) ||
        ((libole2_major_version == major) && (libole2_minor_version == minor) && (libole2_micro_version >= micro)))
      {
        return 0;
       }
     else
      {
        printf("\n*** An old version of libole2 (%d.%d.%d) was found.\n",
               libole2_major_version, libole2_minor_version, libole2_micro_version);
        printf("*** You need a version of libole2 newer than %d.%d.%d. The latest version of\n", major, minor, micro);
        printf("*** libole2 is always available from ftp://ftp.gnome.org.\n");
        printf("***\n");
        printf("*** If you have already installed a sufficiently new version, this error\n");
        printf("*** probably means that the wrong copy of the libole2-config shell script is\n");
        printf("*** being found. The easiest way to fix this is to remove the old version\n");
        printf("*** of libole2, but you can also set the LIBOLE2_CONFIG environment to point to\n");
        printf("*** the correct copy of libole2-config. (In this case, you will have to\n");
        printf("*** modify your LD_LIBRARY_PATH enviroment variable, or edit /etc/ld.so.conf\n");
        printf("*** so that the correct libraries are found at run-time))\n");
      }
    }
  return 1;
}
],, no_libole2=yes,[echo $ac_n "cross compiling; assumed OK... $ac_c"])
       CFLAGS="$ac_save_CFLAGS"
       LIBS="$ac_save_LIBS"
     fi
  fi
  if test "x$no_libole2" = x ; then
     AC_MSG_RESULT(yes)
     ifelse([$2], , :, [$2])     
  else
     AC_MSG_RESULT(no)
     if test "$LIBOLE2_CONFIG" = "no" ; then
       echo "*** The libole2-config script installed by libole2 could not be found"
       echo "*** If libole2 was installed in PREFIX, make sure PREFIX/bin is in"
       echo "*** your path, or set the LIBOLE2_CONFIG environment variable to the"
       echo "*** full path to libole2-config."
     else
       if test -f conf.libole2test ; then
        :
       else
          echo "*** Could not run libole2 test program, checking why..."
          CFLAGS="$CFLAGS $LIBOLE2_CFLAGS"
          LIBS="$LIBS $LIBOLE2_LIBS"
          AC_TRY_LINK([
#include <libole2/libole2.h>
#include <stdio.h>
],      [ return ((libole2_major_version) || (libole2_minor_version) || (libole2_micro_version)); ],
        [ echo "*** The test program compiled, but did not run. This usually means"
          echo "*** that the run-time linker is not finding libole2 or finding the wrong"
          echo "*** version of libole2. If it is not finding libole2, you'll need to set your"
          echo "*** LD_LIBRARY_PATH environment variable, or edit /etc/ld.so.conf to point"
          echo "*** to the installed location. Also, make sure you have run ldconfig if that"
          echo "*** is required on your system"
	  echo "***"
          echo "*** If you have an old version installed, it is best to remove it, although"
          echo "*** you may also be able to get things to work by modifying LD_LIBRARY_PATH" ],
        [ echo "*** The test program failed to compile or link. See the file config.log for the"
          echo "*** exact error that occured. This usually means libole2 was incorrectly"
          echo "*** installed or that you have moved libole2 since it was installed. In the"
          echo "*** latter case, you may want to edit the libole2-config script:"
          echo "*** $LIBOLE2_CONFIG" ])
          CFLAGS="$ac_save_CFLAGS"
          LIBS="$ac_save_LIBS"
       fi
     fi
     LIBOLE2_CFLAGS=""
     LIBOLE2_LIBS=""
     ifelse([$3], , :, [$3])
  fi
  AC_SUBST(LIBOLE2_CFLAGS)
  AC_SUBST(LIBOLE2_LIBS)
  rm -f conf.libole2test
])

# Configure paths for gdk-pixbuf
# Elliot Lee 2000-01-10
# stolen from Raph Levien 98-11-18
# stolen from Manish Singh    98-9-30
# stolen back from Frank Belew
# stolen from Manish Singh
# Shamelessly stolen from Owen Taylor

dnl AM_PATH_GDK_PIXBUF([MINIMUM-VERSION, [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl Test for GDK_PIXBUF, and define GDK_PIXBUF_CFLAGS and GDK_PIXBUF_LIBS
dnl
AC_DEFUN(AM_PATH_GDK_PIXBUF,
[dnl 
dnl Get the cflags and libraries from the gdk-pixbuf-config script
dnl
AC_ARG_WITH(gdk-pixbuf-prefix,[  --with-gdk-pixbuf-prefix=PFX   Prefix where GDK_PIXBUF is installed (optional)],
            gdk_pixbuf_prefix="$withval", gdk_pixbuf_prefix="")
AC_ARG_WITH(gdk-pixbuf-exec-prefix,[  --with-gdk-pixbuf-exec-prefix=PFX Exec prefix where GDK_PIXBUF is installed (optional)],
            gdk_pixbuf_exec_prefix="$withval", gdk_pixbuf_exec_prefix="")
AC_ARG_ENABLE(gdk_pixbuftest, [  --disable-gdk_pixbuftest       Do not try to compile and run a test GDK_PIXBUF program],
		    , enable_gdk_pixbuftest=yes)

  if test x$gdk_pixbuf_exec_prefix != x ; then
     gdk_pixbuf_args="$gdk_pixbuf_args --exec-prefix=$gdk_pixbuf_exec_prefix"
     if test x${GDK_PIXBUF_CONFIG+set} = xset ; then
        GDK_PIXBUF_CONFIG=$gdk_pixbuf_exec_prefix/gdk-pixbuf-config
     fi
  fi
  if test x$gdk_pixbuf_prefix != x ; then
     gdk_pixbuf_args="$gdk_pixbuf_args --prefix=$gdk_pixbuf_prefix"
     if test x${GDK_PIXBUF_CONFIG+set} = xset ; then
        GDK_PIXBUF_CONFIG=$gdk_pixbuf_prefix/bin/gdk-pixbuf-config
     fi
  fi

  AC_PATH_PROG(GDK_PIXBUF_CONFIG, gdk-pixbuf-config, no)
  min_gdk_pixbuf_version=ifelse([$1], ,0.2.5,$1)
  AC_MSG_CHECKING(for GDK_PIXBUF - version >= $min_gdk_pixbuf_version)
  no_gdk_pixbuf=""
  if test "$GDK_PIXBUF_CONFIG" = "no" ; then
    no_gdk_pixbuf=yes
  else
    GDK_PIXBUF_CFLAGS=`$GDK_PIXBUF_CONFIG $gdk_pixbufconf_args --cflags`
    GDK_PIXBUF_LIBS=`$GDK_PIXBUF_CONFIG $gdk_pixbufconf_args --libs`

    gdk_pixbuf_major_version=`$GDK_PIXBUF_CONFIG $gdk_pixbuf_args --version | \
           sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\1/'`
    gdk_pixbuf_minor_version=`$GDK_PIXBUF_CONFIG $gdk_pixbuf_args --version | \
           sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\2/'`
    gdk_pixbuf_micro_version=`$GDK_PIXBUF_CONFIG $gdk_pixbuf_config_args --version | \
           sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\3/'`
    if test "x$enable_gdk_pixbuftest" = "xyes" ; then
      ac_save_CFLAGS="$CFLAGS"
      ac_save_LIBS="$LIBS"
      CFLAGS="$CFLAGS $GDK_PIXBUF_CFLAGS"
      LIBS="$LIBS $GDK_PIXBUF_LIBS"
dnl
dnl Now check if the installed GDK_PIXBUF is sufficiently new. (Also sanity
dnl checks the results of gdk-pixbuf-config to some extent
dnl
      rm -f conf.gdk_pixbuftest
      AC_TRY_RUN([
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gdk-pixbuf/gdk-pixbuf.h>

char*
my_strdup (char *str)
{
  char *new_str;
  
  if (str)
    {
      new_str = malloc ((strlen (str) + 1) * sizeof(char));
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

  system ("touch conf.gdk_pixbuftest");

  /* HP/UX 9 (%@#!) writes to sscanf strings */
  tmp_version = my_strdup("$min_gdk_pixbuf_version");
  if (sscanf(tmp_version, "%d.%d.%d", &major, &minor, &micro) != 3) {
     printf("%s, bad version string\n", "$min_gdk_pixbuf_version");
     exit(1);
   }

   if (($gdk_pixbuf_major_version > major) ||
      (($gdk_pixbuf_major_version == major) && ($gdk_pixbuf_minor_version > minor)) ||
      (($gdk_pixbuf_major_version == major) && ($gdk_pixbuf_minor_version == minor) && ($gdk_pixbuf_micro_version >= micro)))
    {
      return 0;
    }
  else
    {
      printf("\n*** 'gdk-pixbuf-config --version' returned %d.%d.%d, but the minimum version\n", $gdk_pixbuf_major_version, $gdk_pixbuf_minor_version, $gdk_pixbuf_micro_version);
      printf("*** of GDK_PIXBUF required is %d.%d.%d. If gdk-pixbuf-config is correct, then it is\n", major, minor, micro);
      printf("*** best to upgrade to the required version.\n");
      printf("*** If gdk-pixbuf-config was wrong, set the environment variable GDK_PIXBUF_CONFIG\n");
      printf("*** to point to the correct copy of gdk-pixbuf-config, and remove the file\n");
      printf("*** config.cache before re-running configure\n");
      return 1;
    }
}

],, no_gdk_pixbuf=yes,[echo $ac_n "cross compiling; assumed OK... $ac_c"])
       CFLAGS="$ac_save_CFLAGS"
       LIBS="$ac_save_LIBS"
     fi
  fi
  if test "x$no_gdk_pixbuf" = x ; then
     AC_MSG_RESULT(yes)
     ifelse([$2], , :, [$2])     
  else
     AC_MSG_RESULT(no)
     if test "$GDK_PIXBUF_CONFIG" = "no" ; then
       echo "*** The gdk-pixbuf-config script installed by GDK_PIXBUF could not be found"
       echo "*** If GDK_PIXBUF was installed in PREFIX, make sure PREFIX/bin is in"
       echo "*** your path, or set the GDK_PIXBUF_CONFIG environment variable to the"
       echo "*** full path to gdk-pixbuf-config."
     else
       if test -f conf.gdk_pixbuftest ; then
        :
       else
          echo "*** Could not run GDK_PIXBUF test program, checking why..."
          CFLAGS="$CFLAGS $GDK_PIXBUF_CFLAGS"
          LIBS="$LIBS $GDK_PIXBUF_LIBS"
          AC_TRY_LINK([
#include <stdio.h>
#include <gdk-pixbuf/gdk-pixbuf.h>
],      [ return 0; ],
        [ echo "*** The test program compiled, but did not run. This usually means"
          echo "*** that the run-time linker is not finding GDK_PIXBUF or finding the wrong"
          echo "*** version of GDK_PIXBUF. If it is not finding GDK_PIXBUF, you'll need to set your"
          echo "*** LD_LIBRARY_PATH environment variable, or edit /etc/ld.so.conf to point"
          echo "*** to the installed location  Also, make sure you have run ldconfig if that"
          echo "*** is required on your system"
	  echo "***"
          echo "*** If you have an old version installed, it is best to remove it, although"
          echo "*** you may also be able to get things to work by modifying LD_LIBRARY_PATH"],
        [ echo "*** The test program failed to compile or link. See the file config.log for the"
          echo "*** exact error that occured. This usually means GDK_PIXBUF was incorrectly installed"
          echo "*** or that you have moved GDK_PIXBUF since it was installed. In the latter case, you"
          echo "*** may want to edit the gdk-pixbuf-config script: $GDK_PIXBUF_CONFIG" ])
          CFLAGS="$ac_save_CFLAGS"
          LIBS="$ac_save_LIBS"
       fi
     fi
     GDK_PIXBUF_CFLAGS=""
     GDK_PIXBUF_LIBS=""
     ifelse([$3], , :, [$3])
  fi
  AC_SUBST(GDK_PIXBUF_CFLAGS)
  AC_SUBST(GDK_PIXBUF_LIBS)
  rm -f conf.gdk_pixbuftest
])
