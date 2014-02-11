# Configure path for MPI header
# Allin Cottrell, February 2014

AC_DEFUN([AC_OPENMPI_HEADER],
[
#
# Set up configure script macro
#
AC_ARG_WITH(openmpi-include,
	[  --with-openmpi-include=PFX path to Openmpi header],
	[OPENMPI_inc_check="$with_openmpi_include $with_openmpi_include/include $with_openmpi_include/include/openmpi"],
	[OPENMPI_inc_check="/usr/include /usr/local/include /opt/openmpi/include /usr/include/openmpi"])
#
# Look for OpenMPI headers
#
AC_MSG_CHECKING([for OpenMPI include directory])
OPENMPI_incdir=
for m in $OPENMPI_inc_check ; do
  if test -d "$m" && test -f "$m/mpi.h" ; then
    OPENMPI_incdir=$m
    break
  fi
done

if test -z "$OPENMPI_incdir" ; then
  AC_MSG_RESULT([couldn't find mpi.h])
  have_openmpi=no  
else 
  AC_MSG_RESULT([$OPENMPI_incdir])
  have_openmpi=yes
  AC_DEFINE(HAVE_OPENMPI)
  if test "$OPENMPI_incdir" = "/usr/include" ; then
    OPENMPI_CFLAGS=""
  else
    OPENMPI_CFLAGS="-I${OPENMPI_incdir}"
  fi
  AC_SUBST(OPENMPI_CFLAGS)
fi
])
