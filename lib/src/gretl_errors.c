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
#include "gretl_func.h"
#include "libset.h"
#include "monte_carlo.h"

#include <errno.h>

#define EDEBUG 0

#define ERRLEN 2048

static char gretl_errmsg[ERRLEN];
static char gretl_warnmsg[ERRLEN];

static int gretl_errno;
static int gretl_warnnum;

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
    N_("This command won't work with the current periodicity"),  /* E_PDWRONG */
    N_("Error attempting to open file"),                         /* E_FOPEN */
    N_("Out of memory!"),                                        /* E_ALLOC */
    N_("No formula supplied in genr"),                           /* E_EQN */
    N_("Unknown variable name in command"),                      /* E_UNKVAR */
    N_("Command has insufficient arguments"),                    /* E_ARGS */
    N_("This command is implemented only for OLS models"),       /* E_OLSONLY */
    N_("Invalid argument"),                                      /* E_INVARG */
    N_("Syntax error"),                                          /* E_PARSE */
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
    N_("Incompatible options"),                                  /* E_BADOPT */
    N_("The restrictions do not identify the parameters"),       /* E_NOIDENT */
    N_("External command failed"),                               /* E_EXTERNAL */
    N_("Maximum length of command line (65536 bytes) exceeded"), /* E_TOOLONG */
    N_("No dataset is in place"),                                /* E_NODATA */
    N_("Matrix is not positive definite"),                       /* E_NOTPD */
    N_("Failed to calculate Jacobian"),                          /* E_JACOBIAN */
    N_("Insufficient observations for this operation"),          /* E_TOOFEW */
    N_("You cannot define a function within a function"),        /* E_FNEST */
    N_("Out of bounds error"),                                   /* E_BOUNDS */
    N_("Error executing function"),                              /* E_FUNCERR */
    N_("Execution aborted by request"),                          /* E_STOP */
    N_("'catch' cannot be used in this context"),                /* E_BADCATCH */
    N_("complex arguments/operands not supported"),              /* E_CMPLX */
    N_("mixed complex/real arguments not supported"),            /* E_MIXED */
    N_("Dependencies not met"),                                  /* E_DEPENDS */
    NULL,                                                        /* E_DB_DUP */
    NULL,                                                        /* E_OK */
    NULL                                                         /* E_MAX */
};

static const char *gretl_warning_messages[] = {
    NULL,
    N_("gradient is not close to zero"),                 /* W_GRADIENT */
    N_("generated missing values"),                      /* W_GENMISS */
    N_("generated non-finite values"),                   /* W_GENNAN */
    NULL                                                 /* W_MAX */
};

static const char *look_up_errmsg (int err)
{
    if (err > 0 && err < E_MAX) {
	return gretl_error_messages[err];
    } else if (err == 0) {
	return "";
    } else {
	fprintf(stderr, "look_up_errmsg: out of bounds code %d\n", err);
	return "missing error message!";
    }
}

static const char *look_up_warnmsg (int w)
{
    if (w > 0 && w < W_MAX) {
	return gretl_warning_messages[w];
    } else {
	fprintf(stderr, "look_up_warnmsg: out of bounds code %d\n", w);
	return "missing warning message!";
    }
}

static int error_printed;
static int alarm_set;

/**
 * errmsg_get_with_default:
 * @err: gretl error code (see #gretl_error_codes).
 *
 * Returns: a specific error message if available,
 * otherwise a generic error message corresponding to the
 * given @err.
 */

const char *errmsg_get_with_default (int err)
{
    const char *ret = "";

#if EDEBUG
    fprintf(stderr, "errmsg_get_with_default: msg='%s'\n",
	    gretl_errmsg);
#endif

    if ((err > 0 && err < E_MAX) || err == -1) {
	if (*gretl_errmsg != '\0') {
	    ret = gretl_errmsg;
	} else if (err > 0) {
	    const char *deflt = look_up_errmsg(err);

	    if (deflt != NULL) {
		ret = _(deflt);
	    }
	}
    }

    return ret;
}

/**
 * gretl_warnmsg_get:
 *
 * Returns: the current gretl warning message, or %NULL if no
 * warning is currently in place.
 */

const char *gretl_warnmsg_get (void)
{
    const char *ret = NULL;

    if (gretl_warnnum) {
	if (*gretl_warnmsg != '\0') {
	    ret = gretl_warnmsg;
	    /* note; can't zero the message here! */
	} else {
	    const char *deflt = look_up_warnmsg(gretl_warnnum);

	    if (deflt != NULL) {
		ret = _(deflt);
	    }
	}
	gretl_warnnum = 0;
    }

    return ret;
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
#if EDEBUG
    fprintf(stderr, "errmsg: err=%d, error_printed=%d\n",
	    err, error_printed);
#endif

    if (!error_printed && prn != NULL) {
	const char *msg = errmsg_get_with_default(err);

	if (print_redirection_level(prn) > 0) {
	    /* FIXME can we get this message to appear at the
	       "top level" of @prn? */
	    const char *fname = print_redirection_filename(prn);

	    if (fname != NULL) {
		fprintf(stderr, "error when 'outfile' (%s) active\n %s\n",
			fname, msg);
	    } else {
		fprintf(stderr, "error when 'outfile' active\n %s\n", msg);
	    }
	}
	pprintf(prn, "%s\n", msg);
	error_printed = 1;
    }
}

static void print_function_info (PRN *prn)
{
    if (gretl_function_depth() > 0) {
	const char *fname = NULL;
	const char *pname = NULL;

	current_function_info(&fname, &pname);
	if (fname != NULL && pname != NULL) {
	    pprintf(prn, "%s %s (%s %s):\n", _("In regard to function"),
		    fname, _("package"), pname);
	} else if (fname != NULL) {
	    pprintf(prn, "%s %s:\n", _("In regard to function"), fname);
	}
    }
}

/**
 * warnmsg:
 * @prn: gretl printing struct.
 *
 * If a gretl warning is set, prints a message to @prn
 * and zeros the warning signal.
 */

void warnmsg (PRN *prn)
{
    if (prn == NULL || gretl_warnnum == 0) {
	return;
    }

    if (!gretl_warnings_on()) {
	*gretl_warnmsg = '\0';
	gretl_warnnum = 0;
	return;
    }

    if (*gretl_warnmsg != '\0') {
	print_function_info(prn);
	pprintf(prn, "%s: %s\n", _("Warning"), gretl_warnmsg);
	*gretl_warnmsg = '\0';
    } else {
	const char *s = look_up_warnmsg(gretl_warnnum);

	print_function_info(prn);
	pprintf(prn, "%s: %s\n", _("Warning"), _(s));
    }

    gretl_warnnum = 0;
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
#if EDEBUG
    fprintf(stderr, "gretl_errmsg_set: '%s'\n", str);
#endif

    if (alarm_set && *gretl_errmsg != '\0') {
	/* leave the current error message in place */
	return;
    }

    if (*gretl_errmsg == '\0') {
	strncat(gretl_errmsg, str, ERRLEN - 1);
    } else if (strcmp(gretl_errmsg, str)) {
	/* should we do the following? */
	int n = strlen(gretl_errmsg);
	int m = strlen(str);

	if (n + m + 2 < ERRLEN) {
	    strcat(gretl_errmsg, "\n");
	    strcat(gretl_errmsg, str);
	}
    }

#if EDEBUG
    fprintf(stderr, "gretl_errmsg now: '%s'\n", gretl_errmsg);
#endif
}

/**
 * gretl_errmsg_append:
 * @str: an error message.
 * @err: the current error state, if any.
 *
 * Add @str to the current gretl error message, starting a
 * new line, if space permits.
 */

void gretl_errmsg_append (const char *str, int err)
{
    const char *s = NULL;
    gchar *tmp = NULL;

    if (*gretl_errmsg) {
	s = tmp = g_strdup(gretl_errmsg);
    } else if (err > 0 && err < E_MAX) {
	s = look_up_errmsg(err);
    }
    if (s == NULL) {
	g_snprintf(gretl_errmsg, ERRLEN, "%s", str);
    } else {
	g_snprintf(gretl_errmsg, ERRLEN, "%s\n%s", s, str);
    }

    g_free(tmp);
}

/**
 * gretl_errmsg_prepend:
 * @str: an error message.
 * @err: the current error state, if any.
 *
 * Prepend @str to the current gretl error message, ending with a
 * newline, if space permits.
 */

void gretl_errmsg_prepend (const char *str, int err)
{
    const char *s = NULL;
    gchar *tmp = NULL;

    if (*gretl_errmsg) {
	s = tmp = g_strdup(gretl_errmsg);
    } else if (err > 0 && err < E_MAX) {
	s = look_up_errmsg(err);
    }
    if (s == NULL) {
	g_snprintf(gretl_errmsg, ERRLEN, "%s", str);
    } else if (str[strlen(str)-1] != '\n') {
	g_snprintf(gretl_errmsg, ERRLEN, "%s\n%s", str, s);
    } else {
	g_snprintf(gretl_errmsg, ERRLEN, "%s%s", str, s);
    }

    g_free(tmp);
}

/**
 * gretl_errmsg_ensure:
 * @str: an error message.
 *
 * If %gretl_errmsg is currently blank, copy the given string into
 * the message space.
 */

void gretl_errmsg_ensure (const char *str)
{
    if (*gretl_errmsg == '\0') {
	strncat(gretl_errmsg, str, ERRLEN - 1);
    }
}

/**
 * gretl_warnmsg_set:
 * @str: a warning message.
 *
 * Copy the given string into the warning message space.
 */

void gretl_warnmsg_set (const char *str)
{
    *gretl_warnmsg = '\0';
    strncat(gretl_warnmsg, str, ERRLEN - 1);
    gretl_warnnum = W_MAX;
}

/**
 * gretl_errmsg_sprintf:
 * @fmt: format string.
 * @...: arguments, as to sprintf.
 *
 * Append a formatted message to the current gretl
 * error message.
 */

void gretl_errmsg_sprintf (const char *fmt, ...)
{
#if EDEBUG
    fprintf(stderr, "gretl_errmsg_sprintf: fmt='%s'\n", fmt);
#endif

    if (*gretl_errmsg == '\0') {
	va_list ap;

	va_start(ap, fmt);
	vsnprintf(gretl_errmsg, ERRLEN, fmt, ap);
	va_end(ap);
    } else if (strstr(gretl_errmsg, "*** error in fun") &&
	       strstr(fmt, "*** error in fun")) {
	/* don't print more than one "error in function"
	   message, as this gets confusing
	*/
	;
    } else {
	/* find the number of characters left */
	int len0 = strlen(gretl_errmsg);
	int n = ERRLEN - len0 - 2;

	if (n > 31) {
	    char tmp[ERRLEN];
	    va_list ap;

	    *tmp = '\0';
	    va_start(ap, fmt);
	    vsnprintf(tmp, n, fmt, ap);
	    va_end(ap);

	    if (gretl_errmsg[len0 - 1] != '\n') {
		strcat(gretl_errmsg, "\n");
	    }
	    strcat(gretl_errmsg, tmp);
	}
    }
}

void gretl_errmsg_sprintf_replace (const char *fmt, ...)
{
    va_list ap;

    gretl_errmsg[0] = '\0';
    va_start(ap, fmt);
    vsnprintf(gretl_errmsg, ERRLEN, fmt, ap);
    va_end(ap);
}

/**
 * gretl_warnmsg_sprintf:
 * @fmt: format string.
 * @...: arguments, as to sprintf.
 *
 * Write a formatted message to the current gretl
 * warning message space.
 */

void gretl_warnmsg_sprintf (const char *fmt, ...)
{
    va_list ap;

    *gretl_warnmsg = '\0';

    va_start(ap, fmt);
    vsnprintf(gretl_warnmsg, ERRLEN, fmt, ap);
    va_end(ap);

    gretl_warnnum = W_MAX;
}

char *gretl_strerror (int errnum)
{
#if 0 /* doesn't work, AND not cross-platform */
    static locale_t loc = (locale_t) 0;

    if (loc == (locale_t) 0) {
	loc = newlocale(LC_ALL_MASK, "", loc);
    }

    if (loc != (locale_t) 0) {
	uselocale(loc);
	return strerror_l(errnum, loc);
    } else {
	return strerror(errnum);
    }
#else
    return strerror(errnum);
#endif
}

/**
 * gretl_errmsg_set_from_errno:
 * @s: string to prepend to error message, or %NULL.
 *
 * If %gretl_errmsg is currently blank, copy the string
 * returned by %strerror into the message space; or if the
 * error message is not blank but sufficient space remains,
 * append the new error info to the message.
 */

void gretl_errmsg_set_from_errno (const char *s, int errnum)
{
    char *msg = NULL;

    if (errnum) {
	msg = gretl_strerror(errnum);
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

int gretl_error_clear (void)
{
#if EDEBUG
    fprintf(stderr, "gretl_error_clear\n");
#endif
    if (!alarm_set) {
	*gretl_errmsg = '\0';
    }
    error_printed = 0;
    errno = 0;

    return 0;
}

/**
 * gretl_errmsg_is_set:
 *
 * Returns: 1 if the gretl error message is currently
 * set (not blank), otherwise 0.
 */

int gretl_errmsg_is_set (void)
{
    return (*gretl_errmsg != '\0');
}

/**
 * maybe_save_gretl_errmsg:
 *
 * Returns: An allocated copy of the current gretl error
 * message, or NULL if @err is zero or the current message
 * is blank.
 */

char *maybe_save_gretl_errmsg (int err)
{
    if (err && *gretl_errmsg != '\0') {
	return gretl_strdup(gretl_errmsg);
    } else {
	return NULL;
    }
}

/* setting the "alarm" prevents gretl_errmsg from being
   overwritten
*/

void set_gretl_alarm (int val)
{
    alarm_set = val;
}

static void real_set_gretl_errno (int err)
{
    gretl_errno = err;
}

void set_gretl_errno (int err)
{
    real_set_gretl_errno(err);
}

void set_gretl_warning (int w)
{
    gretl_warnnum = w;
}

int get_gretl_errno (void)
{
    int err = gretl_errno;

    real_set_gretl_errno(0);
    return err;
}

int check_gretl_errno (void)
{
    return gretl_errno;
}

int check_gretl_warning (void)
{
    return gretl_warnnum;
}

int gretl_error_is_fatal (void)
{
    if (gretl_compiling_function()) {
	return 1;
    } else if (gretl_compiling_loop()) {
	return 1;
    } else {
	return gretl_in_batch_mode();
    }
}

int invalid_field_error (const char *s)
{
    gretl_errmsg_sprintf(_("field '%s' in command is invalid"), s);
    return E_DATA;
}
