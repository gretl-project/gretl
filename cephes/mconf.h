/*							mconf.h
 *
 *	Common include file for math routines
 *
 *
 *
 * SYNOPSIS:
 *
 * #include "mconf.h"
 *
 *
 *
 * DESCRIPTION:
 *
 * This file contains definitions for error codes that are
 * passed to the common error handling routine mtherr()
 * (which see).
 *
 * The file also includes a conditional assembly definition
 * for the type of computer arithmetic (IEEE, DEC, Motorola
 * IEEE, or UNKnown).
 * 
 * For Digital Equipment PDP-11 and VAX computers, certain
 * IBM systems, and others that use numbers with a 56-bit
 * significand, the symbol DEC should be defined.  In this
 * mode, most floating point constants are given as arrays
 * of octal integers to eliminate decimal to binary conversion
 * errors that might be introduced by the compiler.
 *
 * For little-endian computers, such as IBM PC, that follow the
 * IEEE Standard for Binary Floating Point Arithmetic (ANSI/IEEE
 * Std 754-1985), the symbol IBMPC should be defined.  These
 * numbers have 53-bit significands.  In this mode, constants
 * are provided as arrays of hexadecimal 16 bit integers.
 *
 * Big-endian IEEE format is denoted MIEEE.  On some RISC
 * systems such as Sun SPARC, double precision constants
 * must be stored on 8-byte address boundaries.  Since integer
 * arrays may be aligned differently, the MIEEE configuration
 * may fail on such machines.
 *
 * To accommodate other types of computer arithmetic, all
 * constants are also provided in a normal decimal radix
 * which one can hope are correctly converted to a suitable
 * format by the available C language compiler.  To invoke
 * this mode, define the symbol UNK.
 *
 * An important difference among these modes is a predefined
 * set of machine arithmetic constants for each.  The numbers
 * MACHEP (the machine roundoff error), MAXNUM (largest number
 * represented), and several other parameters are preset by
 * the configuration symbol.  Check the file const.c to
 * ensure that these values are correct for your computer.
 *
 * Configurations NANS, INFINITIES, MINUSZERO, and DENORMAL
 * may fail on many systems.  Verify that they are supposed
 * to work on your computer.
 */
/*
Cephes Math Library Release 2.3:  June, 1995
Copyright 1984, 1987, 1989, 1995 by Stephen L. Moshier
*/

#ifndef CEPHES_MCONF_H
#define CEPHES_MCONF_H

#include <math.h>

#ifdef _WIN32
# include "winconfig.h"
#else
# include "config.h"
#endif

/* Name of package */
#define PROB_PACKAGE "cephes"

/* Version number of package */
#define VERSION "2.7"

/* math error conditions */
enum {
    CEPHES_DOMAIN = 1,  /* argument domain error */
    CEPHES_SING,        /* argument singularity */
    CEPHES_OVERFLOW,    /* overflow range error */
    CEPHES_UNDERFLOW,   /* underflow range error */
    CEPHES_TLOSS,       /* total loss of precision */
    CEPHES_PLOSS,       /* partial loss of precision */
    CEPHES_UNKNOWN      /* unspecified error */
};

/* complex number */
typedef struct {
    double r;
    double i;
} cmplx;

#define CMPLX

#ifdef HAVE_LONG_DOUBLE
/* Long double complex number */
typedef struct {
    long double r;
    long double i;
} cmplxl;
#endif

/* Type of computer arithmetic */

/* UNKnown arithmetic, invokes coefficients given in
 * normal decimal format.  Beware of range boundary
 * problems (MACHEP, MAXLOG, etc. in const.c) and
 * roundoff problems in pow.c:
 * (Sun SPARCstation)
 */
#define UNK 1

#ifndef isnan
# define isnan(x) ((x) != (x))
#endif

#ifndef isfinite
# define isfinite(x) (!isnan(x) && !isinf(x))
#endif

/* Define to support tiny denormal numbers, else undefine. */
#undef DENORMAL

/* Define to distinguish between -0.0 and +0.0 */
#define MINUSZERO 

/* Mechanism for error reporting.  See mtherr.c. */
int mtherr (char *, int);
int mtherr_with_arg (char *, int, double);

/* includes shared private functions */
#include "cephes.h"

#endif /* CEPHES_MCONF_H */
