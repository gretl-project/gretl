/* 
 *  gretl -- Gnu Regression, Econometrics and Time-series Library
 *  Copyright (C) 2001 Allin Cottrell and Riccardo "Jack" Lucchetti
 * 
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 * 
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */

#include "libgretl.h"

#include <errno.h>

char gretl_errmsg[ERRLEN];

static const char *gretl_error_messages[] = {
    NULL,
    NULL,
    N_("Data error"),                                            /* E_DATA = 2 */
    N_("Exact or near collinearity encountered"),                /* E_SINGULAR */
    N_("Insufficient degrees of freedom for regression"),        /* E_DF */
    N_("Dependent variable is all zeros, aborting regression"),  /* E_ZERO */
    N_("Total sum of squares was not positive"),                 /* E_TSS */
    N_("Sum of squared residuals negative!"),                    /* E_ESS */
    N_("Sorry, command not available for this estimator"),       /* E_NOTIMP */
    N_("Unspecified error -- FIXME"),                            /* E_UNSPEC */
    N_("Syntax error in genr formula"),                          /* E_SYNTAX */
    N_("This command won't work with the current periodicity"),  /* E_PDWRONG */
    N_("Error attempting to open file"),                         /* E_FOPEN */
    N_("Out of memory!"),                                        /* E_ALLOC */
    N_("No formula supplied in genr"),                           /* E_EQN */
    N_("Unknown variable name in command"),                      /* E_UNKVAR */
    N_("Command has insufficient arguments"),                    /* E_ARGS */
    N_("This command is implemented only for OLS models"),       /* E_OLSONLY */
    N_("Invalid argument for function"),                         /* E_INVARG */
    N_("Syntax error in command line"),                          /* E_PARSE */
    N_("No independent variables left after omissions"),         /* E_NOVARS */
    N_("No independent variables were omitted"),                 /* E_NOOMIT */
    N_("No new independent variables were added"),               /* E_NOADD */
    N_("One or more \"added\" vars were already present"),       /* E_ADDDUP */
    N_("Error generating logarithms"),                           /* E_LOGS */
    N_("Error generating squares"),                              /* E_SQUARES */
    N_("Error generating lagged variables"),                     /* E_LAGS */
    N_("Attempting to take square root of negative number"),     /* E_SQRT */
    N_("Excessive exponent in genr formula"),                    /* E_HIGH */
    N_("Need valid starting and ending observations"),           /* E_OBS */
    N_("You must include a constant in this sort of model"),     /* E_NOCONST */
    N_("The statistic you requested is not available"),          /* E_BADSTAT */
    N_("Missing sub-sample information; can't merge data"),      /* E_NOMERGE */
    N_("The convergence criterion was not met"),                 /* E_NOCONV */
    N_("The operation was canceled"),                            /* E_CANCEL */
    N_("Missing values encountered"),                            /* E_MISSDATA */
    N_("Not a Number in calculation"),                           /* E_NAN */
    N_("Matrices not conformable for operation"),                /* E_NONCONF */
    N_("Data types not conformable for operation"),              /* E_TYPES */
    N_("Wrong data type"),                                       /* E_DATATYPE */
    N_("Incompatible options"),                                  /* E_BADOPT */
    N_("The restrictions do not identify the parameters"),       /* E_NOIDENT */
    N_("External command failed"),                               /* E_EXTERNAL */
    N_("Maximum length of command line (8192 bytes) exceeded"),  /* E_TOOLONG */ 
    N_("Argument is not a scalar"),                              /* E_NOTSCALAR */ 
    NULL,                                                        /* E_DB_DUP */
    NULL,                                                        /* E_OK */
    NULL                                                         /* E_MAX */
};

static const char *look_up_errmsg (int err)
{
    const char *ret = NULL;

    if (err > 0 && err < E_MAX) {
	ret = gretl_error_messages[err];
    } else {
	fprintf(stderr, "look_up_errmsg: out of bounds errcode %d\n", 
		err);
    }

    return ret;
}

static int error_printed;

/**
 * errmsg_get_with_default: 
 * @err: gretl error code (see #error_codes).
 *
 * Returns: a specific error message if available,
 * otherwise a generic error message corresponding to the 
 * given @err.
 */

const char *errmsg_get_with_default (int err)
{
    const char *msg = "";

    if (*gretl_errmsg != '\0') {
	msg = gretl_errmsg;
    } else {
	const char *deflt = look_up_errmsg(err);

	if (deflt != NULL) {
	    msg = _(deflt);
	}
    }

    return msg;
}

/**
 * errmsg:
 * @err: gretl error code (see #error_codes).
 * @prn: gretl printing struct.
 *
 * Prints to @prn a specific error message if available, 
 * otherwise a generic error message corresponding to the 
 * given @err.
 */

void errmsg (int err, PRN *prn)
{
    if (!error_printed && prn != NULL) {
	const char *msg = errmsg_get_with_default(err);

	pprintf(prn, "%s\n", msg);
	error_printed = 1;
    } 
}

/**
 * gretl_errmsg_get:
 *
 * Returns: a specific error message if available,
 * otherwise an empty string.
 */

const char *gretl_errmsg_get (void)
{
    return gretl_errmsg;
}

/**
 * gretl_errmsg_set:
 * @str: an error message.
 *
 * If %gretl_errmsg is currently blank, copy the given string into
 * the message space; or if the error message is not blank but
 * sufficient space remains, append @str to the message.
 */

void gretl_errmsg_set (const char *str)
{
    if (*gretl_errmsg == '\0') {
	strncat(gretl_errmsg, str, ERRLEN - 1);
    } else {
	/* should we do the following? */
	int n = strlen(gretl_errmsg);
	int m = strlen(str);

	if (n + m + 1 < ERRLEN) {
	    strcat(gretl_errmsg, "\n");
	    strcat(gretl_errmsg, str);
	}
    }
}

/**
 * gretl_errmsg_sprintf:
 * @fmt: format string.
 *
 * If %gretl_errmsg is currently blank, print a formatted
 * message into place.
 */

void gretl_errmsg_sprintf (const char *fmt, ...)
{
    if (*gretl_errmsg == '\0' && fmt != NULL) {
	va_list args;

	va_start(args, fmt);
	vsnprintf(gretl_errmsg, ERRLEN, fmt, args);
	va_end(args);
    }
}

/**
 * gretl_errmsg_set_from_errno:
 *
 * If %gretl_errmsg is currently blank, copy the string 
 * returned by %strerror into the message space; or if the 
 * error message is not blank but sufficient space remains, 
 * append the new error info to the message.
 */

void gretl_errmsg_set_from_errno (const char *s)
{
    char *msg = NULL;

    if (errno) {
	msg = strerror(errno);
	errno = 0;
    }

    if (msg != NULL) {
	if (s != NULL) {
	    gretl_errmsg_sprintf("%s: %s", s, msg);
	} else {
	    gretl_errmsg_set(msg);
	}
    }
}

/**
 * gretl_error_clear:
 *
 * Blank out any previously recorded error message.
 */

void gretl_error_clear (void)
{
    *gretl_errmsg = '\0';
    error_printed = 0;
    errno = 0;
}

#if 0 /* not yet */

static char *warnbuf;

void gretl_warning_add (const char *s)
{
    free(warnbuf);
    warnbuf = gretl_strdup(s);
}

const char *gretl_warning_get (void)
{
    if (warnbuf != NULL) {
	return warnbuf;
    } else {
	return "";
    }
}

#endif


