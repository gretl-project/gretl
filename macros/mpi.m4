# Configure path for MPI header
# Note: we're assuming that mpicc will take care of finding the MPI library
# (which is only required for linking of gretlmpi), so here we're mostly
# just checking that mpi.h can be found, and recording the include path if
# necessary
# Allin Cottrell, February 2014, revised November 2017

AC_DEFUN([AC_C_MPI],
[
# mpi-lib not used, just retained for backward compatibility
AC_ARG_WITH(mpi-lib,
    [  --with-mpi-lib=PFX path to MPI library],,)

AC_ARG_WITH(mpi-include,
    [  --with-mpi-include=PFX path to MPI header],
    [MPI_inc_check="$with_mpi_include $with_mpi_include/include $with_mpi_include/include/openmpi"],
    [MPI_inc_check="/usr/include /usr/local/include /usr/include/openmpi /opt/openmpi/include /usr/lib/openmpi/include"])

# can we find the header just using the default search paths?
AC_CHECK_HEADER(mpi.h,have_mpi="yes",,)

# if not, try our own search
if test $have_mpi = no ; then
  ARCH=`uname -m`
  ARCHPATHS="/usr/include/openmpi-${ARCH} /usr/lib/${ARCH}-linux-gnu/openmpi/include /usr/local/include/openmpi-${ARCH}"
  MPI_incdir=
  for m in $MPI_inc_check $ARCHPATHS ; do
    if test -d "$m" && test -f "$m/mpi.h" ; then
      MPI_incdir=$m
      break
    fi
  done
  if test -z "$MPI_incdir" ; then
    AC_MSG_RESULT([couldn't find mpi.h])
  else
    AC_MSG_RESULT([found mpi.h in $MPI_incdir])
    have_mpi=yes
    if test "$MPI_incdir" != "/usr/include" ; then
      MPI_CFLAGS="-I${MPI_incdir}"
    fi
  fi
fi

if test $have_mpi = yes ; then
  AC_DEFINE(HAVE_MPI)
  if test x"${MPICC}" = x ; then
    MPICC=mpicc
  fi 
  AC_SUBST(MPICC)
  AC_SUBST(MPI_CFLAGS)
fi
])
