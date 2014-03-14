# Configure paths for MPI library and header
# Allin Cottrell, February 2014

AC_DEFUN([AC_C_MPI],
[
dnl
dnl Check the architecture
dnl
ARCH=`uname -m`
if test "${ARCH}" = "i386" -o "${ARCH}" = "i686" -o "${ARCH}" = "ppc" -o "${ARCH}" = "sparcv9" ; then
   MODE=32
elif test "${ARCH}" = "ia64" ; then
   MODE=64
elif test "${ARCH}" = "s390" ; then
   MODE=31
elif test "${ARCH}" = "s390x" ; then
   MODE=64
elif test "${ARCH}" = "x86_64" -o "${ARCH}" = "ppc64" -o "${ARCH}" = "sparc64" ; then
   MODE=64
fi

#
# Set up configure script macros
#
AC_ARG_WITH(mpi-lib,
    [  --with-mpi-lib=PFX path to MPI library],
    [MPI_lib_check="$with_mpi_lib $with_mpi_lib/lib${MODE} $with_mpi_lib/lib $with_mpi_lib/lib${MODE}/openmpi $with_mpi_lib/lib/openmpi"],
    [MPI_lib_check="/usr/lib /usr/lib/openmpi /usr/lib64/openmpi/lib"])

AC_ARG_WITH(mpi-include,
    [  --with-mpi-include=PFX path to MPI header],
    [MPI_inc_check="$with_mpi_include $with_mpi_include/include $with_mpi_include/include/openmpi"],
    [MPI_inc_check="/usr/include /usr/include/openmpi /usr/include/openmpi-${ARCH} /usr/local/include /usr/local/include/openmpi-${ARCH} /opt/openmpi/include /usr/include/openmpi /usr/lib/openmpi/include"])

#
# Look for MPI library
#
AC_MSG_CHECKING([for MPI library])
MPI_libdir=
for m in ${MPI_lib_check} ; do
  if test -d "$m" ; then
    if (test -f "$m/libmpi.so" || test -f "$m/libmpi.a" || test -f "$m/libmpi.dylib"); then
      MPI_libdir=$m
      break
    fi
  fi
done

#
# if not found, try for the MPICH library instead
#
if test -z "$MPI_libdir" ; then
  for m in ${MPI_lib_check} ; do
    if test -d "$m" ; then
      if (test -f "$m/libmpich.so" || test -f "$m/libmpich.a" || test -f "$m/libmpich.dylib"); then
        MPI_libdir=$m
        break
      fi
    fi
  done
fi

if test -z "$MPI_libdir"
then
  AC_MSG_RESULT([Couldn't find MPI library])
  have_mpi=no
else
  AC_MSG_RESULT([$MPI_libdir])
  if test "MPI_libdir" = "/usr/lib" ; then
    MPI_LIBS="-lmpi"
  else
    MPI_LIBS="-L$MPI_libdir -lmpi"
  fi
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
    if test x"${MPICC}" = x ; then
      MPICC=mpicc
    fi 
    AC_SUBST(MPICC)
    AC_SUBST(MPI_CFLAGS)
    AC_SUBST(MPI_LIBS)
  fi
fi
])
