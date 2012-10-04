/*							mtherr.c
 *
 *	Library common error handling routine
 *
 *
 *
 * SYNOPSIS:
 *
 * char *fctnam;
 * int code;
 * int mtherr();
 *
 * mtherr( fctnam, code );
 *
 *
 *
 * DESCRIPTION:
 *
 * This routine may be called to report one of the following
 * error conditions (in the include file mconf.h).
 *  
 *   Mnemonic        Value          Significance
 *
 *    DOMAIN            1       argument domain error
 *    SING              2       function singularity
 *    OVERFLOW          3       overflow range error
 *    UNDERFLOW         4       underflow range error
 *    TLOSS             5       total loss of precision
 *    PLOSS             6       partial loss of precision
 *    EDOM             33       Unix domain error code
 *    ERANGE           34       Unix range error code
 *
 * The default version of the file prints the function name,
 * passed to it by the pointer fctnam, followed by the
 * error condition.  The display is directed to the standard
 * output device.  The routine then returns to the calling
 * program.  Users may wish to modify the program to abort by
 * calling exit() under severe error conditions such as domain
 * errors.
 *
 * Since all error conditions pass control to this function,
 * the display may be easily changed, eliminated, or directed
 * to an error logging device.
 *
 * SEE ALSO:
 *
 * mconf.h
 *
 */

/*
Cephes Math Library Release 2.0:  April, 1987
Copyright 1984, 1987 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#include <stdio.h>
#include "mconf.h"

static int cephes_errno = 0;

/* Notice: the order of appearance of the following
 * messages is bound to the error codes defined
 * in mconf.h.
 */
static const char *ermsg[] = {
    "no error",
    "domain",       /* error code 1 */
    "singularity",  /* et seq.      */
    "overflow",
    "underflow",
    "total loss of precision",
    "partial loss of precision",
    "unknown"
};

/* @name is supposed to be the name of the function in
   which the error occurred; @code is an index into
   the array of error messages above; @arg is the offending
   argument, if @have_arg is non-zero, otherwise it is
   ignored.
*/

static int real_mtherr (char *name, int code, double arg,
			int have_arg, int quiet)
{
    static int hush;

    if (name == NULL) {
	hush = quiet;
	return 0;
    }

    if (hush) {
	return 0;
    }

    fprintf(stderr, "%s ", name);

    if (code <= 0 || code > CEPHES_UNKNOWN)
	code = CEPHES_UNKNOWN;

    cephes_errno = code;

    if (have_arg) {
	fprintf(stderr, "%s error (arg = %g)\n", ermsg[code], arg);
    } else {
	fprintf(stderr, "%s error\n", ermsg[code]);
    }

    return 0;
}

int mtherr_with_arg (char *name, int code, double arg)
{
    return real_mtherr(name, code, arg, 1, 0);
}

int mtherr (char *name, int code)
{
    return real_mtherr(name, code, 0, 0, 0);
}

int get_cephes_errno (void)
{
    int ret = cephes_errno;

    cephes_errno = 0; /* clear the code */

    return ret;
}

void set_cephes_hush (int s)
{
    real_mtherr(NULL, 0, 0, 0, s);
}
