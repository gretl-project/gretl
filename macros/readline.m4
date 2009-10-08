
# Configure paths for readline, for gretl
# Allin Cottrell, April 2003

dnl AM_PATH_READLINE([ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl Test for READLINE, and define READLINE_CFLAGS and READLINE_LIBS.
dnl
AC_DEFUN([AM_PATH_READLINE],
[dnl 
AC_ARG_WITH(readline-prefix,[  --with-readline-prefix=PFX   Prefix where readline is installed (optional)],
            readline_config_prefix="$withval", readline_config_prefix="")

  READLINE_CFLAGS="-I$readline_config_prefix/include"
  READLINE_LIBS="-L$readline_config_prefix/lib -lreadline"

  ac_save_CFLAGS="$CFLAGS"
  ac_save_LIBS="$LIBS"
  CFLAGS="$CFLAGS $READLINE_CFLAGS"
  LIBS="$READLINE_LIBS $LIBS"

  dnl check for readline header
  AC_CHECK_HEADER(readline/readline.h, have_rl_header="yes", have_rl_header="no")
  
  if test "$have_rl_header" = "yes" ; then
    dnl check for the libraries that readline depends on
    AC_CHECK_LIB(ncurses, tgetent, termcap_lib=-lncurses,
    [AC_CHECK_LIB(termcap, tgetent, termcap_lib=-ltermcap,
    [AC_CHECK_LIB(curses, tgetent, termcap_lib=-lcurses,
    termcap_lib='')])])

    AC_CHECK_LIB(readline, readline, have_readline="yes" ; \
       AC_DEFINE(HAVE_READLINE),,$termcap_lib)
    AC_SUBST(have_readline)

    if test "$have_readline" = "yes" ; then 
      AC_CHECK_LIB(readline, rl_completion_matches, AC_DEFINE(NEW_READLINE),,$termcap_lib)
      AC_SUBST(new_readline)

      dnl remove any extraneous stuff from the flags and libs lines
      if test "$READLINE_CFLAGS" = "-I/include" ; then 
         READLINE_CFLAGS=""
      fi
      if test "$READLINE_LIBS" = "-L/lib -lreadline" ; then 
         READLINE_LIBS="-lreadline"
      fi

      dnl append the termcap lib portion
      READLINE_LIBS="$READLINE_LIBS $termcap_lib"

      dnl we're now ready to write stuff out
      AC_SUBST(READLINE_CFLAGS)
      AC_SUBST(READLINE_LIBS)
    fi
  fi

  CFLAGS="$ac_save_CFLAGS"
  LIBS="$ac_save_LIBS"
],[])
