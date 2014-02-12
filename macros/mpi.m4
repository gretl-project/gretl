# Configure path for MPI header
# Allin Cottrell, February 2014

AC_DEFUN([AC_MPI_HEADER],
[
dnl
dnl Check the architecture
dnl
ARCH=`uname -m`
#
# Set up configure script macro
#
AC_ARG_WITH(mpi-include,
	[  --with-mpi-include=PFX path to MPI header],
	[MPI_inc_check="$with_mpi_include $with_mpi_include/include $with_mpi_include/include/openmpi"],
	[MPI_inc_check="/usr/include /usr/include/openmpi /usr/include/openmpi-${ARCH} /usr/local/include /usr/local/include/openmpi-${ARCH} /opt/openmpi/include /usr/include/openmpi"])
#
# Look for MPI headers
#
AC_MSG_CHECKING([for MPI include directory])
MPI_incdir=
for m in $MPI_inc_check ; do
  if test -d "$m" && test -f "$m/mpi.h" ; then
    MPI_incdir=$m
    break
  fi
done

if test -z "$MPI_incdir" ; then
  AC_MSG_RESULT([couldn't find mpi.h])
  have_mpi=no  
else 
  AC_MSG_RESULT([$MPI_incdir])
  have_mpi=yes
  AC_DEFINE(HAVE_MPI)
  if test "$MPI_incdir" = "/usr/include" ; then
    MPI_CFLAGS=""
  else
    MPI_CFLAGS="-I${MPI_incdir}"
  fi
  AC_SUBST(MPI_CFLAGS)
fi
])
