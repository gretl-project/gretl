# Configure paths for openmpi
# Allin Cottrell <cottrell@wfu.edu>, last updated February 2014
# based on AX_OPENMPI by Denis Arnaud

AC_DEFUN([AM_PATH_OPENMPI],
[

dnl
dnl Check the CPU (and OS) architecture: 32 or 64 bits
dnl
ARCH=`uname -m`
if test "${ARCH}" = "i386" -o "${ARCH}" = "i686" -o "${ARCH}" = "ppc" -o "${ARCH}" = "sparcv9"
then
   MODE=32
   MODEFLAG="-m32"
elif test "${ARCH}" = "ia64"
then
   MODE=64
elif test "${ARCH}" = "s390"
then
   MODE=31
elif test "${ARCH}" = "s390x"
then
   MODE=64
elif test "${ARCH}" = "x86_64" -o "${ARCH}" = "ppc64" -o "${ARCH}" = "sparc64"
then
   MODE=64
   MODEFLAG="-m64"
fi

GENARCH="${ARCH}"
if test "${ARCH}" = "i686"
then
  GENARCH="i386" 
fi

OPENMPI_DEF_VERSION="1.4"
OPENMPI_DEF_CC="gcc"
MPIDIR="${OPENMPI_DEF_VERSION}-${OPENMPI_DEF_CC}"

	#
	# Set up configure script macros
	#
	AC_ARG_WITH(openmpi,
		[  --with-openmpi=<path>    root directory path of Openmpi installation],
		[OPENMPI_lib_check="$with_openmpi/lib${MODE}/openmpi $with_openmpi/lib/openmpi $with_openmpi/lib${MODE} $with_openmpi/lib"
		OPENMPI_inc_check="$with_openmpi/include $with_openmpi/include/openmpi $with_openmpi/include/openmpi/${MPIDIR}"
		OPENMPI_bin_check="$with_openmpi/bin $with_openmpi/lib${MODE}/openmpi/bin $with_openmpi/lib/openmpi/bin $with_openmpi/lib${MODE}/bin $with_openmpi/lib/bin"],
		[OPENMPI_lib_check="/usr/lib${MODE} /usr/lib /usr/lib${MODE}/openmpi /usr/lib/openmpi /usr/lib${MODE}/openmpi/lib /usr/lib/openmpi/lib /usr/lib${MODE}/openmpi/${MPIDIR} /usr/lib/openmpi/${MPIDIR} /usr/local/lib${MODE} /usr/local/lib /usr/local/lib${MODE}/openmpi /usr/local/lib/openmpi /usr/local/lib${MODE}/openmpi/${MPIDIR} /usr/local/lib/openmpi/${MPIDIR} /usr/local/openmpi/lib${MODE} /usr/local/openmpi/lib /usr/local/openmpi/${MPIDIR}/lib${MODE} /usr/local/openmpi/${MPIDIR}/lib /opt/openmpi/lib /opt/openmpi/${MPIDIR}/lib /opt/openmpi/lib/openmpi"
		OPENMPI_inc_check="/usr/include /usr/include/openmpi /usr/include/openmpi-${ARCH} /usr/include/openmpi-${GENARCH} /usr/include/openmpi/${MPIDIR} /usr/local/include /usr/local/include/openmpi /usr/local/include/openmpi/${MPIDIR} /opt/openmpi/include /opt/openmpi/${MPIDIR}/include"
		OPENMPI_bin_check="/usr/bin /usr/lib${MODE}/bin /usr/lib/bin /usr/lib${MODE}/openmpi/bin /usr/lib/openmpi/bin /usr/lib${MODE}/openmpi/lib/bin /usr/lib/openmpi/lib/bin /usr/lib${MODE}/openmpi/${MPIDIR}/bin /usr/lib/openmpi/${MPIDIR}/bin /usr/local/bin /usr/local/lib${MODE}/bin /usr/local/lib/bin /usr/local/lib${MODE}/openmpi/bin /usr/local/lib/openmpi/bin /usr/local/lib${MODE}/openmpi/${MPIDIR}/bin /usr/local/lib/openmpi/${MPIDIR}/bin /usr/local/openmpi/lib${MODE}/bin /usr/local/openmpi/lib/bin /usr/local/openmpi/${MPIDIR}/lib${MODE}/bin /usr/local/openmpi/${MPIDIR}/lib/bin /opt/bin /opt/openmpi/bin /opt/openmpi/lib/bin /opt/openmpi/${MPIDIR}/lib/bin /opt/openmpi/lib/openmpi/bin"])

	AC_ARG_WITH(openmpi-lib,
		[  --with-openmpi-lib=<path> directory path of Openmpi library installation],
		[OPENMPI_lib_check="$with_openmpi_lib $with_openmpi_lib/lib${MODE} $with_openmpi_lib/lib $with_openmpi_lib/lib${MODE}/openmpi $with_openmpi_lib/lib/openmpi $with_openmpi_lib/lib${MODE}/openmpi/${MPIDIR} $with_openmpi_lib/lib/openmpi/${MPIDIR}"])

	AC_ARG_WITH(openmpi-include,
		[  --with-openmpi-include=<path> directory path of Openmpi header installation],
		[OPENMPI_inc_check="$with_openmpi_include $with_openmpi_include/include $with_openmpi_include/include/openmpi $with_openmpi_include/include/openmpi/${MPIDIR}"])

	# OpenMPI library
	OPENMPI_LIB=mpi

	#
	# Look for OpenMPI library
	#
	AC_MSG_CHECKING([for OpenMPI library directory])
	OPENMPI_libdir=
	for m in ${OPENMPI_lib_check}
	do
		if test -d "$m"
		then
			if (test -f "$m/lib$OPENMPI_LIB.so" || test -f "$m/lib$OPENMPI_LIB.a")
			then
				OPENMPI_libdir=$m
				break
			fi
		fi
	done

	if test -z "$OPENMPI_libdir"
	then
		AC_MSG_ERROR([Didn't find $OPENMPI_LIB library in '$OPENMPI_lib_check'])
	fi

	case "$OPENMPI_libdir" in
		/* ) ;;
		* )  AC_MSG_ERROR([The Openmpi library directory ($OPENMPI_libdir) must be an absolute path.]) ;;
	esac

	AC_MSG_RESULT([$OPENMPI_libdir])

	#
	# Look for OpenMPI headers
	#
	AC_MSG_CHECKING([for OpenMPI include directory])
	OPENMPI_incdir=
	for m in $OPENMPI_inc_check
	do
		if test -d "$m" && test -f "$m/mpi.h"
		then
			OPENMPI_incdir=$m
			break
		fi
	done

	if test -z "$OPENMPI_incdir"
	then
		AC_MSG_ERROR([Didn't find the Openmpi include dir in '$OPENMPI_inc_check'])
	fi

	case "$OPENMPI_incdir" in
		/* ) ;;
		* )  AC_MSG_ERROR([The Openmpi include directory ($OPENMPI_incdir) must be an absolute path.]) ;;
	esac

	AC_MSG_RESULT([$OPENMPI_incdir])

	#
	# Look for OpenMPI configuration
	#
	# OpenMPI configuration binary
	OPENMPI_VERBIN=ompi_info

	AC_MSG_CHECKING([for OpenMPI configuration])
	OPENMPI_verdir=
	for m in ${OPENMPI_bin_check}
	do
		if test -d "$m"
		then
			if (test -x "$m/$OPENMPI_VERBIN")
			then
				OPENMPI_verdir=$m
				break
			fi
		fi
	done

	if test -z "$OPENMPI_verdir"
	then
		AC_MSG_RESULT([Did not find $OPENMPI_VERBIN tool, but will try to determine the version with other tricks])
	else
		case "$OPENMPI_verdir" in
			/* ) ;;
			* )  AC_MSG_ERROR([The Openmpi tool directory ($OPENMPI_verdir) must be an absolute path.]) ;;
		esac
	fi

	AC_MSG_RESULT([$OPENMPI_verdir])

	#
	# Look for OpenMPI version
	#
	min_openmpi_version=ifelse([$1], , 1.4.1, $1)
	AC_MSG_CHECKING([for OpenMPI version >= $min_openmpi_version])

	# Derive the full version, e.g., 1.3.3
	if test ! -z "$OPENMPI_verdir"
	then
		OPENMPI_VERSION=`${OPENMPI_verdir}/${OPENMPI_VERBIN} -v ompi full | grep "Open MPI:" | cut -d: -f2 | cut -d\  -f2`
	fi
	if test -z "${OPENMPI_VERSION}" ; then
	   OPENMPI_VERSION=${OPENMPI_DEF_VERSION}
	fi

	# Calculation of major, minor, micro
    openmpi_major_version=`echo $OPENMPI_VERSION | sed 's/^\([[0-9]]*\).*/\1/'`
    if test "x${openmpi_major_version}" = "x" ; then
       openmpi_major_version=0
    fi

    openmpi_minor_version=`echo $OPENMPI_VERSION | \
						sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\2/'`
    if test "x${openmpi_minor_version}" = "x" ; then
       openmpi_minor_version=0
    fi

    openmpi_micro_version=`echo $OPENMPI_VERSION | \
          sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\3/'`
    if test "x${openmpi_micro_version}" = "x" ; then
       openmpi_micro_version=0
    fi
	
	AC_MSG_RESULT([$OPENMPI_VERSION])

	#
	# Derive the flags for compilation and linking
	#
	OPENMPI_CFLAGS="-I${OPENMPI_incdir} -I${OPENMPI_incdir}/openmpi -I${OPENMPI_incdir}/${MODE} ${MODEFLAG}"
	OPENMPI_LIBS="-L${OPENMPI_libdir} -lmpi_cxx -lmpi"

	#
	# Make the changes permanent
	#
	HAVE_OPENMPI=1
	MPIGEN_VERSION=${OPENMPI_VERSION}
	MPIGEN_CFLAGS=${OPENMPI_CFLAGS}
	MPIGEN_LIBS=${OPENMPI_LIBS}
	AC_SUBST(OPENMPI_VERSION)
	AC_SUBST(OPENMPI_CFLAGS)
	AC_SUBST(OPENMPI_LIBS)

	# Test linking with openmpi
	save_LIBS="$LIBS"
	if test -z "$OPENMPI_LIBS"
	then
		OPENMPI_LIBS="-L/usr/lib${MODE}/openmpi/${MPIDIR} -lmpi_cxx -lmpi"
	fi
	LIBS="$LIBS $OPENMPI_LIBS"
	AC_CHECK_LIB($OPENMPI_LIB, ompi_mpi_char,
                 [], 
                 [AC_MSG_ERROR([Could not find working Openmpi library!])]
                )
	LIBS="$save_LIBS"
])
