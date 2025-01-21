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

/* support for C-like (s)printf and sscanf */

#include "libgretl.h"
#include "libset.h"
#include "gretl_func.h"
#include "matrix_extra.h"
#include "gretl_string_table.h"
#include "uservar.h"

#include <errno.h>

#define PSDEBUG 0

/* For the point of "x *= (1 + 0x1p-52)" below, see
   http://sourceware.org/bugzilla/show_bug.cgi?id=4943 ,
   contribution from Eric Postpischil, 2007-10-05.
   It is designed to produce the generally expected
   rounding behavior on printing a double.
*/

static void printf_series (int v, const DATASET *dset,
			   const char *fmt,
			   int wid, int prec,
			   int wstar, int pstar,
			   PRN *prn)
{
    int obslen = max_obs_marker_length(dset);
    double x;
    int t;

    for (t=dset->t1; t<=dset->t2; t++) {
	print_obs_marker(t, dset, obslen, prn);
	x = dset->Z[v][t];
	if (!na(x)) {
	    x *= (1 + 0x1p-52);
	    if (wstar && pstar) {
		pprintf(prn, fmt, wid, prec, x);
	    } else if (wstar || pstar) {
		wid = (wstar)? wid : prec;
		pprintf(prn, fmt, wid, x);
	    } else {
		pprintf(prn, fmt, x);
	    }
	}
	pputc(prn, '\n');
    }
}

static int printf_escape (int c, PRN *prn)
{
    switch (c) {
    case 'n':
	pputc(prn, '\n');
	break;
    case 'r':
        pputc(prn, '\r');
	break;
    case 't':
	pputc(prn, '\t');
	break;
    case 'v':
	pputc(prn, '\v');
	break;
    case '\\':
	pputc(prn, '\\');
	break;
    case '"':
	pputc(prn, '"');
	break;
    default:
	/* treat as literal backslash */
	pputc(prn, '\\');
	pputc(prn, c);
    }

    return 0;
}

static char *printf_get_string (char *s, DATASET *dset,
				int t, int *err)
{
    char *ret = NULL;

    if (!strncmp(s, "marker", 6) && dset != NULL && dset->S != NULL) {
	/* special, not in regular genr */
	const char *p = dset->S[t];
	const char *q = s + 6;
	int len = strlen(p);
	int offset = 0;

	while (isspace(*q)) q++;

	if (*q == '+') {
	    offset = atoi(++q);
	    len -= offset;
	}
	if (len >= 0) {
	    ret = gretl_strndup(p + offset, len);
	    if (ret == NULL) {
		*err = E_ALLOC;
	    }
	}
    } else {
	ret = generate_string(s, dset, err);
	if (*err && dset != NULL && t >= 0) {
	    /* using a string-valued series? */
	    int v = current_series_index(dset, s);
	    const char *str;

	    if (v > 0 && is_string_valued(dset, v)) {
		str = series_get_string_for_obs(dset, v, t);
		if (str != NULL) {
		    ret = gretl_strdup(str);
		    *err = 0;
		}
	    }
	}
    }

    if (ret == NULL && !*err) {
	*err = E_DATA;
    }

    return ret;
}

/* here we're trying to get a field width or precision
   specifier */

static int printf_get_int (const char *s, DATASET *dset,
			   int *err)
{
    int ret = 0;

#if PSDEBUG
    fprintf(stderr, "printf_get_int: looking at '%s'\n", s);
#endif

    ret = gretl_int_from_string(s, err);

    if (*err) {
	*err = 0;
	gretl_error_clear();
	/* we'll allow a bit of slop with regard to
	   integerhood here */
	ret = (int) generate_scalar(s, dset, err);
    }

    if (!*err && abs(ret) > 255) {
	*err = E_DATA;
    }

#if PSDEBUG
    fprintf(stderr, "printf_get_int: returning %d\n", ret);
#endif

    return ret;
}

/* gen_arg_val: when this function is called we have already
   determined that the argument string @s is not associated with a
   string format specifier, is not a literal numerical value, and is
   not the name of an existing scalar, matrix or series.  So now we
   attempt to generate a value (of type unknown).
*/

static int gen_arg_val (const char *s, DATASET *dset,
			double *px, gretl_matrix **pm,
			int *pv, int *scalar)
{
    const char *name = "$ptmp";
    char *genstr;
    int err = 0;

    genstr = gretl_strdup_printf("%s=%s", name, s);
    if (genstr == NULL) {
	return E_ALLOC;
    }

    err = generate(genstr, dset, GRETL_TYPE_ANY, OPT_P, NULL);

    if (!err) {
	int type = genr_get_last_output_type();

	if (type == GRETL_TYPE_DOUBLE) {
	    *scalar = 1;
	    *px = gretl_scalar_get_value(name, NULL);
	    user_var_delete_by_name(name, NULL);
	} else if (type == GRETL_TYPE_MATRIX) {
	    gretl_matrix *m = get_matrix_by_name(name);

	    *pm = gretl_matrix_copy(m);
	    if (*pm == NULL) {
		err = E_ALLOC;
	    }
	    user_var_delete_by_name(name, NULL);
	} else if (type == GRETL_TYPE_SERIES) {
	    *pv = current_series_index(dset, name);
	} else {
	    err = E_TYPES; /* FIXME cleanup? */
	}
    }

    free(genstr);

    return err;
}

/* dup argv (@s) up to the next free comma */

static char *get_next_arg (const char *s, int *len, int *err)
{
    const char *p;
    char *arg = NULL;
    int par = 0, br = 0;
    int quoted = 0;
    int n = 0;

    if (s == NULL) {
	*err = E_PARSE;
	return NULL;
    }

    *len = strspn(s, ", ");
    s += *len;
    p = s;

    while (*p) {
	if (*p == '"' && (p == s || *(p-1) != '\\')) {
	    quoted = !quoted;
	}
	if (!quoted) {
	    if (*p == '(') par++;
	    else if (*p == ')') par--;
	    else if (*p == '[') br++;
	    else if (*p == ']') br--;
	}
	if (!quoted && !par && !br && *p == ',') {
	    break;
	}
	p++;
    }

    n = p - s;
    *len += n;

    if (n > 0) {
	arg = gretl_strndup(s, n);
	if (arg == NULL) {
	    *err = E_ALLOC;
	}
    } else {
	*err = E_ARGS;
    }

#if PSDEBUG
    fprintf(stderr, "get_next_arg: got '%s'\n", arg);
#endif

    return arg;
}

/* conversion characters that may be preceded by 'l' */
#define l_ok(c) (c == 'd' || c == 'u' || c == 'x')

/* conversion characters that may be preceded by 'v' */
#define v_ok(c) (c == 'f' || c == 'g' || c == 'e' || c == 'E')

/* dup format string up to the end of the current conversion:
   e.g. "%10.4f", "%6g", "%.8g", %3s", "%.*s", "%#12.4g"
*/

static char *
get_printf_format_chunk (const char *s, char *conv,
			 int *len, int *wstar, int *pstar,
			 int *err)
{
    const char *cnvchars = "eEfgGduxs";
    const char *numchars = "0123456789";
    char *chunk = NULL;
    const char *p = s;
    int n;

    p++; /* move past % */

    /* '#', if present, must come first */
    if (*p == '#') {
	p++;
    }

    /* zero padding? */
    if (*p == '0') {
	p++;
    }

    /* left justification? */
    if (*p == '-') {
	p++;
    }

    /* leading space before positive number? */
    if (*p == ' ') {
	p++;
    }

    /* always print sign? */
    if (*p == '+') {
	p++;
    }

    /* optional width? */
    n = strspn(p, numchars);
    if (n == 0 && *p == '*') {
	/* variable version */
	*wstar = 1;
	p++;
    } else {
	p += n;
    }

    /* optional dot plus precision? */
    if (*p == '.') {
	p++;
	n = strspn(p, numchars);
	if (n == 0 && *p == '*') {
	    *pstar = 1;
	    p++;
	} else {
	    p += n;
	}
    }

    /* now we should have a conversion character */
    if (*p == '\0' || strchr(cnvchars, *p) == NULL) {
	if (*p == 'l' && l_ok(*(p+1))) {
	    conv[0] = *p;
	    conv[1] = *(p+1);
	    p++;
	} else if (*p == 'v') {
	    if (v_ok(*(p+1))) {
		/* allow "vf", "vg" for matrices */
		conv[0] = *p;
		conv[1] = *(p+1);
		p++;
	    } else {
		/* allow plain 'v' -> "vg" */
		conv[0] = *p;
	    }
	} else {
	    fprintf(stderr, "bad conversion '%c'\n", *p);
	    *err = E_PARSE;
	    return NULL;
	}
    }

    if (!conv[0]) {
	conv[0] = *p;
    }
    p++;
    *len = n = p - s;

    if (n > 0) {
	chunk = gretl_strndup(s, n);
	if (chunk == NULL) {
	    *err = E_ALLOC;
	}
    } else {
	*err = E_PARSE;
    }

#if PSDEBUG
    fprintf(stderr, "get_printf_format_chunk: got '%s'\n", chunk);
#endif

    return chunk;
}

static int printf_as_int (double x, int fc,
			  int lflag, char *fmt,
			  int wstar, int pstar,
			  int wid, int prec,
			  PRN *prn)
{
    if ((fc == 'u' || fc == 'x') && lflag) {
	if (x < 0 || x > (double) ULONG_MAX) {
	    return E_INVARG;
	}
	if (wstar && pstar) {
	    pprintf(prn, fmt, wid, prec, (unsigned long) x);
	} else if (wstar || pstar) {
	    wid = (wstar)? wid : prec;
	    pprintf(prn, fmt, wid, (unsigned long) x);
	} else {
	    pprintf(prn, fmt, (unsigned long) x);
	}
    } else if (fc == 'u' || fc == 'x') {
	if (x < 0 || x > (double) UINT_MAX) {
	    return E_INVARG;
	}
	if (wstar && pstar) {
	    pprintf(prn, fmt, wid, prec, (unsigned) x);
	} else if (wstar || pstar) {
	    wid = (wstar)? wid : prec;
	    pprintf(prn, fmt, wid, (unsigned) x);
	} else {
	    pprintf(prn, fmt, (unsigned) x);
	}
    } else if (fc == 'd' && lflag) {
	if (x < (double) LONG_MIN || x > (double) LONG_MAX) {
	    return E_INVARG;
	}
	if (wstar && pstar) {
	    pprintf(prn, fmt, wid, prec, (long) x);
	} else if (wstar || pstar) {
	    wid = (wstar)? wid : prec;
	    pprintf(prn, fmt, wid, (long) x);
	} else {
	    pprintf(prn, fmt, (long) x);
	}
    } else {
	if (x < (double) INT_MIN || x > (double) INT_MAX) {
	    return E_INVARG;
	}
	if (wstar && pstar) {
	    pprintf(prn, fmt, wid, prec, (int) x);
	} else if (wstar || pstar) {
	    wid = (wstar)? wid : prec;
	    pprintf(prn, fmt, wid, (int) x);
	} else {
	    pprintf(prn, fmt, (int) x);
	}
    }

    return 0;
}

/* Allow for UTF-8 strings where the number of characters
   is not equal to the number of bytes: we'll adjust a
   printf specification such as "%10s" so that the result
   comes out 10 characters wide, rather than 10 bytes.
*/

static void maybe_adjust_string_format (char *str,
					int *wstar,
					int *pstar,
					int *wid,
					int *prec,
					char **fmt)
{
    int adj = strlen(str) - g_utf8_strlen(str, -1);

    if (adj > 0) {
	const char *s = *fmt;
	char *p, *tmp = NULL;
	int uwid = 0, uprec = 0;

	if (*wstar) {
	    /* got width in variable form */
	    uwid = *wid;
	} else {
	    /* maybe got inline width? */
	    if (isdigit(s[1]) || s[1] == '-') {
		uwid = atoi(s + 1);
	    }
	}

	if (*pstar) {
	    /* got "precision" in variable form */
	    p = g_utf8_offset_to_pointer(str, *prec);
	    uprec = p - str;
	} else if ((p = strchr(s, '.')) != NULL) {
	    /* got inline precision? */
	    uprec = atoi(p + 1);
	    p = g_utf8_offset_to_pointer(str, uprec);
	    uprec = p - str;
	}

	if (uwid != 0) {
	    if (uprec > 0 && uprec < strlen(str)) {
		gchar *trunc = g_strndup(str, uprec);

		adj = strlen(trunc) - g_utf8_strlen(trunc, -1);
		g_free(trunc);
	    }
	    uwid = (uwid < 0)? uwid - adj : uwid + adj;
	}

	if (uwid != 0 && uprec > 0) {
	    tmp = g_strdup_printf("%%%d.%ds", uwid, uprec);
	} else if (uwid != 0) {
	    tmp = g_strdup_printf("%%%ds", uwid);
	} else if (uprec > 0) {
	    tmp = g_strdup_printf("%%.%ds", uprec);
	}

	if (tmp != NULL) {
	    /* replace the incoming format */
	    free(*fmt);
	    *fmt = gretl_strdup(tmp);
	    g_free(tmp);
	    *wstar = *pstar = 0;
	}
    }
}

/* Extract the next conversion spec from *pfmt, find and evaluate the
   corresponding elements in *pargs, and print the result. Note
   that @t should be negative (and therefore ignored) unless
   we're doing the "genr markers" special thing.
*/

static int print_arg (const char **pfmt, const char **pargs,
		      DATASET *dset, int t, PRN *prn)
{
    const char *intconv = "dxu";
    char *fmt = NULL;
    char *arg = NULL;
    char *str = NULL;
    char conv[3] = {0};
    gretl_matrix *m = NULL;
    int free_m = 0;
    int series_v = -1;
    int free_v = 0;
    double x = NADBL;
    int flen = 0, alen = 0;
    int wstar = 0, pstar = 0;
    int wid = 0, prec = 0;
    int got_scalar = 0;
    int fc = 0;
    int err = 0;

#if PSDEBUG
    fprintf(stderr, "print_arg: *pfmt='%s', *pargs='%s'\n",
	    *pfmt, *pargs);
#endif

    /* select current conversion format */
    fmt = get_printf_format_chunk(*pfmt, conv, &flen, &wstar, &pstar, &err);
    if (err) {
	return err;
    }

    *pfmt += flen;

    if (conv[0] == 'l' && l_ok(conv[1])) {
	/* e.g. "ld", "lx" */
	fc = conv[1];
    } else {
	fc = conv[0];
    }

    if (wstar) {
	/* evaluate field width specifier */
	arg = get_next_arg(*pargs, &alen, &err);
	if (!err) {
	    x = printf_get_int(arg, dset, &err);
	}
	free(arg);
	if (err) {
	    goto bailout;
	}
	*pargs += alen;
	wid = x;
    }

    if (pstar) {
	/* evaluate precision specifier */
	arg = get_next_arg(*pargs, &alen, &err);
	if (!err) {
	    x = printf_get_int(arg, dset, &err);
	}
	free(arg);
	if (err) {
	    goto bailout;
	}
	*pargs += alen;
	prec = x;
    }

    /* get next substantive arg */
    arg = get_next_arg(*pargs, &alen, &err);

    if (!err && !strcmp(arg, "time") &&
	current_series_index(dset, arg) < 0) {
	/* pre-generate "time" if needed */
	err = gen_time(dset, 1, NULL);
    }

    if (!err) {
	int v = -1;

	if (fc == 's') {
	    /* must be printing a string */
	    str = printf_get_string(arg, dset, t, &err);
	} else if (numeric_string(arg)) {
	    /* printing a scalar value */
	    x = dot_atof(arg);
	    got_scalar = 1;
	} else if (gretl_is_scalar(arg)) {
	    /* printing a named scalar */
	    x = gretl_scalar_get_value(arg, NULL);
	    got_scalar = 1;
	} else if ((m = get_matrix_by_name(arg)) != NULL) {
	    ; /* printing a named matrix */
	} else if ((v = current_series_index(dset, arg)) >= 0) {
	    /* printing a named series, unless t >= 0 */
	    if (t < 0) {
		series_v = v;
	    } else {
		x = dset->Z[v][t];
		got_scalar = 1;
	    }
	} else {
	    /* try treating as generated value */
	    err = gen_arg_val(arg, dset, &x, &m, &v, &got_scalar);
	    if (m != NULL) {
		free_m = 1;
	    } else if (v > 0) {
		series_v = v;
		free_v = 1;
	    }
	}
	*pargs += alen;
    }

    free(arg);

    if (!err && fc == 'v' && m == NULL) {
	/* 'v' for variant is for matrices only */
	err = E_PARSE;
    }

    if (err) {
	goto bailout;
    }

    if (m != NULL && m->rows == 1 && m->cols == 1) {
	/* print 1 x 1 matrix as scalar? */
	if (!isnan(m->val[0])) {
	    x = m->val[0];
	    got_scalar = 1;
	}
    }

    if (got_scalar && na(x)) {
	fc = fmt[flen - 1] = 's';
	str = gretl_strdup("NA");
    }

    /* do the actual printing */

    if (m != NULL && !got_scalar) {
	/* printing a matrix */
	if (!wstar) wid = -1;
	if (!pstar) prec = -1;
	gretl_matrix_print_with_format(m, fmt, wid, prec, prn);
    } else if (series_v >= 0) {
	/* printing a series */
	printf_series(series_v, dset, fmt, wid, prec,
		      wstar, pstar, prn);
    } else if (fc == 's') {
	/* printing a string */
	if (!strcmp(fmt+1, "s")) {
	    pputs(prn, str);
	} else {
	    maybe_adjust_string_format(str, &wstar, &pstar,
				       &wid, &prec, &fmt);
	    if (wstar && pstar) {
		pprintf(prn, fmt, wid, prec, str);
	    } else if (wstar || pstar) {
		wid = (wstar)? wid : prec;
		pprintf(prn, fmt, wid, str);
	    } else {
		pprintf(prn, fmt, str);
	    }
	}
    } else if (strchr(intconv, fc)) {
	/* printing a scalar as int */
	int lflag = (conv[0] == 'l');

	err = printf_as_int(x, fc, lflag, fmt, wstar, pstar,
			    wid, prec, prn);
    } else if (got_scalar) {
	/* printing a (non-missing) scalar value */
	if (!isnan(x) && !isinf(x)) {
	    char *s;

	    if (pstar && prec >= 15) {
		; /* don't mess with it */
	    } else if ((s = strchr(fmt, '.')) && isdigit(*(s+1))) {
		if (atoi(s+1) < 15) {
		    x *= (1 + 0x1p-52);
		}
	    } else {
		x *= (1 + 0x1p-52);
	    }
	}
	if (wstar && pstar) {
	    pprintf(prn, fmt, wid, prec, x);
	} else if (wstar || pstar) {
	    wid = (wstar)? wid : prec;
	    pprintf(prn, fmt, wid, x);
	} else {
	    pprintf(prn, fmt, x);
	}
    }

 bailout:

    if (err) {
	pputc(prn, '\n');
    }

    free(fmt);
    free(str);

    if (free_m) {
	gretl_matrix_free(m);
    } else if (free_v) {
	dataset_drop_last_variables(dset, 1);
    }

#if PSDEBUG
    fprintf(stderr, "print_arg: returning %d\n", err);
#endif

    return err;
}

/* supports both regular printf and generation of obs markers */

static int real_do_printf (const char *format,
			   const char *args,
			   DATASET *dset, PRN *inprn,
			   int *nchars, int t)
{
    PRN *prn = NULL;
    int err = 0;

    if (inprn == NULL) {
	return 0;
    }

    gretl_error_clear();

#if PSDEBUG
    fprintf(stderr, "real_do_printf:\n");
    fprintf(stderr, " format = '%s'\n", format);
    fprintf(stderr, " args = '%s'\n", args);
#endif

    /* We'll buffer the output locally in case there's an
       error part-way through the printing.
    */
    prn = gretl_print_new(GRETL_PRINT_BUFFER, &err);
    if (err) {
        return err;
    }

    if (format != NULL) {
	const char *p = format;
	const char *q = args;

	while (*p && !err) {
	    if (*p == '%' && *(p+1) == '%') {
		pputc(prn, '%');
		p += 2;
	    } else if (*p == '%') {
		err = print_arg(&p, &q, dset, t, prn);
	    } else if (*p == '\\') {
		err = printf_escape(*(p+1), prn);
		p += 2;
	    } else {
		pputc(prn, *p);
		p++;
	    }
	}

	if (!err && q != NULL && *q != '\0') {
	    gchar *s = g_strdup_printf(_("unprocessed argument(s): '%s'"), q);

	    gretl_errmsg_sprintf("printf: %s", s);
	    g_free(s);
	    err = E_PARSE;
	}
    }

    if (err) {
	pputc(inprn, '\n');
    } else {
	const char *buf = gretl_print_get_buffer(prn);

	if (buf != NULL) {
	    if (nchars != NULL) {
		*nchars = strlen(buf);
	    }
	    pputs(inprn, buf);
	}
    }

    gretl_print_destroy(prn);

    return err;
}

/**
 * do_printf:
 * @format: format string.
 * @args: list of arguments as string.
 * @dset: dataset struct.
 * @prn: printing struct.
 * @nchars: location to receive number of characters printed.
 *
 * Implements a somewhat limited version of C's printf()
 * for use in gretl scripts.
 *
 * Returns: 0 on success, non-zero on error.
 */

int do_printf (const char *format, const char *args,
	       DATASET *dset, PRN *prn, int *nchars)
{
    return real_do_printf(format, args, dset, prn, nchars, -1);
}

/**
 * do_sprintf_function:
 * @format: format string.
 * @args: list of arguments as string.
 * @dset: dataset struct.
 * @err: location to receive error code.
 *
 * Implements the hansl-function form of sprintf.
 *
 * Returns: constructed string, or NULL on error.
 */

char *do_sprintf_function (const char *format, const char *args,
			   DATASET *dset, int *err)
{
    const char *p = format;
    const char *q = args;
    char *buf = NULL;
    PRN *prn;

    if (format == NULL || *format == '\0') {
	gretl_errmsg_set("sprintf: format string is missing");
	*err = E_DATA;
	return NULL;
    }

    prn = gretl_print_new(GRETL_PRINT_BUFFER, err);
    if (*err) {
	fprintf(stderr, "sprintf: creation of buffer failed\n");
	return NULL;
    }

    while (*p && !*err) {
	if (*p == '%' && *(p+1) == '%') {
	    pputc(prn, '%');
	    p += 2;
	} else if (*p == '%') {
	    *err = print_arg(&p, &q, dset, -1, prn);
	} else if (*p == '\\') {
	    printf_escape(*(p+1), prn);
	    p += 2;
	} else {
	    pputc(prn, *p);
	    p++;
	}
    }

    if (!*err && q != NULL && *q != '\0') {
	gchar *s = g_strdup_printf(_("unprocessed argument(s): '%s'"), q);

	gretl_errmsg_sprintf("sprintf: %s", s);
	g_free(s);
	*err = E_PARSE;
    }

    if (!*err) {
	buf = gretl_print_steal_buffer(prn);
    }

    gretl_print_destroy(prn);

    return buf;
}

/* below: sscanf apparatus */

struct bracket_scan {
    int reverse;
    char *chrs;
};

static int bscan_escape (char *s)
{
    int i, n = strlen(s);

    for (i=0; i<n; i++) {
	if (s[i] == '\\' && (i == 0 || s[i-1] != '\\')) {
	    if (s[i+1] == 'n') {
		s[i] = '\n';
		shift_string_left(s + i + 1, 1);
		i++;
	    } else if (s[i+1] == 't') {
		s[i] = '\t';
		shift_string_left(s + i + 1, 1);
		i++;
	    }
	}
    }

    return 0;
}

static int
parse_bracket_scanner (const char **pfmt, struct bracket_scan *bscan)
{
    const char *p, *s = *pfmt;
    int k = 0, n = 0, closed = 0;

    s++;
    k++;

    if (*s == '^') {
	bscan->reverse = 1;
	s++;
	k++;
    }
    if (*s == ']') {
	s++;
	k++;
    }

    p = s;
    while (*p) {
	if (*p == ']') {
	    closed = 1;
	    break;
	}
	p++;
	n++;
    }

    if (!closed || n == 0) {
	return E_PARSE;
    }

    bscan->chrs = gretl_strndup(s, n);
    if (bscan->chrs == NULL) {
	return E_ALLOC;
    }

    bscan_escape(bscan->chrs);

#if PSDEBUG
    fprintf(stderr, "bracket_scan: reverse = %d, chars = '%s'\n",
	    bscan->reverse, bscan->chrs);
#endif

    *pfmt += (k + n + 1);

    return 0;
}

#define WIDTH_UNSPEC -1

static char *
get_scanf_format_chunk (const char *s, int *fc, int *len,
			int *width, int *supp,
			struct bracket_scan *bscan,
			int *err)
{
    const char *cnvchars = "dfgsmx";
    const char *numchars = "0123456789";
    char *chunk = NULL;
    const char *p = s;
    int n;

    p++; /* % */

    if (*p == '*') {
	/* suppress assignment */
	*supp = 1;
	p++;
    }

    /* optional max width? */
    n = strspn(p, numchars);
    if (n > 0) {
	*width = atoi(p);
	p += n;
    } else {
	*width = WIDTH_UNSPEC;
    }

    /* now we should have a conversion spec */
    if (*p == '[') {
	*err = parse_bracket_scanner(&p, bscan);
	if (!*err) {
	    *fc = 's';
	}
    } else {
	if (*p == 'l' && *(p+1) == 'f') {
	    p++;
	}
	if (*p == '\0' || strchr(cnvchars, *p) == NULL) {
	    fprintf(stderr, "bad conversion '%c'\n", *p);
	    *err = E_PARSE;
	    return NULL;
	}
	*fc = *p++;
    }

    *len = n = p - s;

    if (n > 0) {
	chunk = gretl_strndup(s, n);
	if (chunk == NULL) {
	    *err = E_ALLOC;
	}
    } else {
	*err = E_PARSE;
    }

#if PSDEBUG
    fprintf(stderr, "get_scanf_format_chunk: got '%s'\n", chunk);
#endif

    return chunk;
}

static int
scan_string (char *targ, const char **psrc, int width,
	     struct bracket_scan *bscan, int *ns)
{
    int n, err = 0;

    *psrc += strspn(*psrc, " \t\n");

    if (bscan->chrs != NULL) {
	if (bscan->reverse) {
	    n = strcspn(*psrc, bscan->chrs);
	} else {
	    n = strspn(*psrc, bscan->chrs);
	}
    } else {
	n = strcspn(*psrc, " \t\n");
    }

#if PSDEBUG
    fprintf(stderr, "scan_string: *psrc='%s', n=%d\n", *psrc, n);
    int m = strcspn(*psrc, "\n");
    fprintf(stderr, "m=%d (chrs='%s')\n", m, bscan->chrs);
#endif

    if (width != WIDTH_UNSPEC && n > width) {
	n = width;
    }

    if (targ != NULL) {
	user_var *uvar;

	uvar = get_user_var_of_type_by_name(targ, GRETL_TYPE_STRING);

	if (uvar == NULL) {
	    gretl_errmsg_sprintf(_("%s: not a string variable"), targ);
	    err = E_DATA;
	} else if (n > 0) {
	    char *conv = gretl_strndup(*psrc, n);

	    if (conv == NULL) {
		err = E_ALLOC;
	    } else {
		user_var_replace_value(uvar, conv, GRETL_TYPE_STRING);
		*psrc += n;
	    }
	    if (!err) {
		*ns += 1;
	    }
	}
    } else {
	*psrc += n;
    }

    return err;
}

static int scan_scalar (char *targ, const char **psrc,
			int fc, int width, DATASET *dset,
			int *ns)
{
    char *endp = NULL;
    long k = 0;
    unsigned long u = 0;
    double x = 0;
    int err = 0;

    errno = 0;

    if (fc == 'd') {
	k = strtol(*psrc, &endp, 10);
    } else if (fc == 'x') {
	u = strtoul(*psrc, &endp, 16);
    } else {
	x = strtod(*psrc, &endp);
    }

    if (width != WIDTH_UNSPEC) {
	int eaten = endp - *psrc;

	if (eaten > width) {
	    char *tmp = gretl_strndup(*psrc, width);

	    if (tmp == NULL) {
		return E_ALLOC;
	    }
	    if (fc == 'd') {
		k = strtol(tmp, &endp, 10);
	    } else if (fc == 'x') {
		u = strtoul(tmp, &endp, 16);
	    } else {
		x = strtod(tmp, &endp);
	    }
	    endp = (char *) *psrc + width;
	    free(tmp);
	}
    }

    /* "genr" operates in the C locale regardless */
    gretl_push_c_numeric_locale();

    if (!errno && endp != *psrc) {
	char *genline;

	if (fc == 'd') {
	    genline = gretl_strdup_printf("%s=%ld", targ, k);
	} else if (fc == 'x') {
	    genline = gretl_strdup_printf("%s=%lu", targ, u);
	} else {
	    genline = gretl_strdup_printf("%s=%.16g", targ, x);
	}
	/* we use GRETL_TYPE_ANY below to allow for the possibility
	   that we're assigning to, e.g., an element of a matrix
	*/
	err = generate(genline, dset, GRETL_TYPE_ANY, OPT_P, NULL);
#if PSDEBUG
	fprintf(stderr, "genline '%s', err=%d\n", genline, err);
#endif
	if (!err) {
	    *ns += 1;
	}
	free(genline);
    }

    /* match the push above */
    gretl_pop_c_numeric_locale();

    *psrc = endp;

    return err;
}

static int scan_matrix (char *targ, const char **psrc,
			int rows, int *ns)
{
    const char *src = *psrc;
    char *endp = NULL;
    gretl_matrix *m = NULL;
    double x;
    int r = 0, c = 0, c0 = 0;
    int i, j, err = 0;

    errno = 0;

    while (*src) {
	strtod(src, &endp);
	if (endp == src || errno == ERANGE) {
	    break;
	}
	c++;
	src = endp;
	src += strspn(src, " \t");
	if (*src == '\r' && *(src+1) == '\n') {
	    src++;
	}
	if (*src == '\n' || *src == '\r' || *src == '\0') {
	    if (r == 0) {
		c0 = c;
	    } else if (c != c0) {
		break;
	    }
	    r++;
	    if (r == rows) {
		break;
	    }
	    c = 0;
	}
    }

    if (targ == NULL) {
	*psrc = src;
	return 0;
    }

    if (get_matrix_by_name(targ) == NULL) {
	return E_UNKVAR;
    }

    if (r > 0 && c0 > 0) {
	m = gretl_matrix_alloc(r, c0);
	if (m == NULL) {
	    return E_ALLOC;
	}

	src = *psrc;
	for (i=0; i<r; i++) {
	    for (j=0; j<c0; j++) {
		x = strtod(src, &endp);
		gretl_matrix_set(m, i, j, x);
		src = endp;
		src += strspn(src, " \t\n\r");
	    }
	}

	err = user_matrix_replace_matrix_by_name(targ, m);
	if (err) {
	    gretl_matrix_free(m);
	} else {
	    /* number of elements scanned */
	    *ns += m->rows * m->cols;
	}
    }

    *psrc = src;

    return err;
}

static int scan_arg (const char **psrc, const char **pfmt, const char **pargs,
		     DATASET *dset, int *ns)
{
    char *fmt = NULL;
    char *arg = NULL;
    char *str = NULL;
    char *targ = NULL;
    struct bracket_scan bscan;
    int flen = 0, alen = 0;
    int wid = 0, supp = 0;
    int fc = 0;
    int err = 0;

#if PSDEBUG
    fprintf(stderr, "scan_arg: *psrc='%s', *pfmt='%s', *pargs='%s'\n",
	    *psrc, *pfmt, *pargs);
#endif

    bscan.reverse = 0;
    bscan.chrs = NULL;

    /* select current conversion format */
    fmt = get_scanf_format_chunk(*pfmt, &fc, &flen, &wid, &supp,
				 &bscan, &err);
    if (err) {
	return err;
    }

    *pfmt += flen;

    if (!supp) {
	/* get next substantive arg */
	arg = get_next_arg(*pargs, &alen, &err);
	if (!err) {
	    targ = arg;
	    if (*targ == '&') {
		targ++;
	    }
	    if (*targ == '\0') {
		err = E_PARSE;
	    } else {
		*pargs += alen;
	    }
	}
    }

#if PSDEBUG
    fprintf(stderr, "scan_arg: targ = '%s'\n", targ);
#endif

    /* do the actual scanning */

    if (!err && fc == 's') {
	err = scan_string(targ, psrc, wid, &bscan, ns);
    } else if (!err && fc == 'm') {
	err = scan_matrix(targ, psrc, wid, ns);
    } else if (!err) {
	err = scan_scalar(targ, psrc, fc, wid, dset, ns);
    }

    free(fmt);
    free(arg);
    free(str);
    free(bscan.chrs);

    return err;
}

/**
 * do_sscanf:
 * @src: the source string.
 * @format: the format string.
 * @args: the arguments, as composite string.
 * @dset: dataset struct.
 * @n_items: location to receive the number of scanned items.
 *
 * Implements a somewhat simplified version of C's sscanf()
 * for use in hansl scripts.
 *
 * Returns: 0 on success, non-zero on error.
 */

int do_sscanf (const char *src, const char *format, const char *args,
	       DATASET *dset, int *n_items)
{
    const char *r, *p, *q;
    int nscan = 0;
    int err = 0;

    gretl_error_clear();

#if PSDEBUG
    fprintf(stderr, "do_sscanf: src = '%s'\n", src);
    fprintf(stderr, "do_sscanf: format = '%s'\n", format);
    fprintf(stderr, "do_sscanf: args = '%s'\n", args);
#endif

    r = src;
    p = format;
    q = args;

    while (*r && *p && !err) {
	if (isspace(*p)) {
	    while (isspace(*r)) r++;
	    while (isspace(*p)) p++;
	}
	if (*r == '\t' && *p == '\\' && *(p+1) == 't') {
	    r++;
	    p += 2;
	} else if (*r == '\n' && *p == '\\' && *(p+1) == 'n') {
	    r++;
	    p += 2;
	} else if (*r == *p) {
	    r++;
	    p++;
	} else if (*p == '%') {
	    err = scan_arg(&r, &p, &q, dset, &nscan);
	} else {
	    break;
	}
    }

    if (!err) {
	*n_items = nscan;
    }

    return err;
}

/* Apparatus to support the command form of printf: the @quoted
   argument indicates whether @parm was in quotes on input.
*/

int do_printf_command (const char *parm, const char *args,
		       DATASET *dset, PRN *prn, int quoted)
{
    const char *format = NULL;

#if PSDEBUG
    fprintf(stderr, "do_printf_command:\n"
	    "  parm='%s', args='%s', quoted=%d\n", parm, args, quoted);
#endif

    if (quoted) {
	format = parm;
    } else if (!strncmp(parm, "T_(", 3)) {
        gretl_errmsg_set("A 'T_()' format-string requires the function "
                         "form of printf()");
    } else {
	/* @parm should be a string variable */
	format = get_string_by_name(parm);
    }

    if (format == NULL) {
	return E_INVARG;
    } else {
	return do_printf(format, args, dset, prn, NULL);
    }
}

/* The originating command is of form:

     genr markers = format, arg1,...

   We're assuming that we're just getting the part after
   the equals sign here, in the variable @s.
*/

int generate_obs_markers (const char *s, DATASET *dset)
{
    char format[16] = {0};
    char args[32] = {0};
    PRN *prn;
    int n, t, err = 0;

    n = sscanf(s, "\"%15[^\"]\", %31[^\r\n]", format, args);
    if (n != 2) {
	return E_PARSE;
    }

    prn = gretl_print_new(GRETL_PRINT_BUFFER, &err);
    if (err) {
	return err;
    }

    if (dset->S == NULL) {
	err = dataset_allocate_obs_markers(dset);
    }

    if (!err) {
	const char *buf;

	for (t=0; t<dset->n && !err; t++) {
	    gretl_print_reset_buffer(prn);
	    err = real_do_printf(format, args, dset, prn, NULL, t);
	    if (!err) {
		buf = gretl_print_get_buffer(prn);
		dset->S[t][0] = '\0';
		strncat(dset->S[t], buf, OBSLEN - 1);
	    }
	}
    }

    gretl_print_destroy(prn);

    return err;
}
