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
#include "uservar.h"
#include "dbwrite.h"
#include "libset.h"
#include "gretl_func.h"
#include "gretl_xml.h"
#include "gretl_panel.h"
#include "gretl_string_table.h"
#include "matrix_extra.h"
#include "csvdata.h"
#include "usermat.h"

#include <time.h>
#include <errno.h>
#include <unistd.h>

#define MERGE_DEBUG 0
#define DATES_DEBUG 0

/**
 * SECTION:dataio
 * @short_description: data handling (internal)
 * @title: Data support
 * @include: gretl/libgretl.h
 *
 * The following data handling functions are basically internal to
 * gretl and not in a state where they can be readily
 * documented as public APIs.
 *
 */

typedef enum {
    GRETL_FMT_GDT,       /* standard gretl XML data */
    GRETL_FMT_BINARY,    /* native gretl binary data */
    GRETL_FMT_OCTAVE,    /* data in Gnu Octave format */
    GRETL_FMT_CSV,       /* data in Comma Separated Values format */
    GRETL_FMT_R,         /* data in Gnu R format */
    GRETL_FMT_DAT,       /* data in PcGive format */
    GRETL_FMT_DB,        /* gretl native database format */
    GRETL_FMT_JM,        /* JMulti ascii data */
    GRETL_FMT_DTA,       /* Stata .dta format */
    GRETL_FMT_JSON       /* geojson (maps) */
} GretlDataFormat;

#define IS_DATE_SEP(c) (c == '.' || c == ':' || c == ',')

#define PROGRESS_BAR "progress_bar"

/**
 * get_date_x:
 * @pd: frequency of data.
 * @obs: observation string.
 *
 * Returns: the floating-point representation of @obs.
 */

double get_date_x (int pd, const char *obs)
{
    double x = 1.0;

    if ((pd == 5 || pd == 6 || pd == 7 || pd == 52) && strlen(obs) > 4) {
	/* calendar data */
	guint32 ed = get_epoch_day(obs);

	if (ed > 0) {
	    x = ed;
	}
    } else {
	x = obs_str_to_double(obs);
    }

    return x;
}

static int real_check_varname (const char *vname,
			       int is_series)
{
    int testchar = 'a';
    int firstbad = 0;
    int err = 0;

    gretl_error_clear();

    if (vname == NULL || *vname == '\0') {
	gretl_errmsg_set(_("Expected an identifier"));
	return E_PARSE;
    }

    if (strlen(vname) >= VNAMELEN) {
	gretl_errmsg_set(_("Varname exceeds the maximum of 31 characters"));
	err = E_DATA;
    } else if (strcmp(vname, "return") && gretl_reserved_word(vname)) {
	err = E_DATA;
    } else if (!(isalpha((unsigned char) *vname))) {
	firstbad = 1;
	testchar = *vname;
        err = E_DATA;
    } else if (is_series && (function_lookup(vname) ||
			     get_user_function_by_name(vname) ||
			     is_function_alias(vname))) {
	gretl_warnmsg_sprintf(_("'%s' shadows a function of the same name"),
			      vname);
    } else {
	const char *p = vname;

	while (*p && testchar == 'a') {
	    if (!(isalpha((unsigned char) *p))
		&& !(isdigit((unsigned char) *p))
		&& *p != '_') {
		testchar = *p;
		err = E_DATA;
	    }
	    p++;
	}
    }

    if (err && !is_series && strlen(vname) == 2 && is_greek_letter(vname)) {
	return 0;
    }

    if (testchar != 'a') {
	if (isprint((unsigned char) testchar)) {
	    if (firstbad) {
		gretl_errmsg_sprintf(_("First char of varname '%s' is bad\n"
				       "(first must be alphabetical)"),
				     vname);
	    } else {
		gretl_errmsg_sprintf(_("Varname '%s' contains illegal character '%c'\n"
				       "Use only letters, digits and underscore"),
				     vname, (unsigned char) testchar);
	    }
	} else {
	    if (firstbad) {
		gretl_errmsg_sprintf(_("First char of varname (0x%x) is bad\n"
				       "(first must be alphabetical)"),
				     (unsigned) testchar);
	    } else {
		gretl_errmsg_sprintf(_("Varname contains illegal character 0x%x\n"
				       "Use only letters, digits and underscore"),
				     (unsigned) testchar);
	    }
	}
    }

    return err;
}

/**
 * check_varname:
 * @varname: putative name for variable (or object).
 *
 * Check a variable/object name for legality: the name
 * must start with a letter, and be composed of letters,
 * numbers or the underscore character, and nothing else.
 *
 * Returns: 0 if name is OK, non-zero if not.
 */

int check_varname (const char *varname)
{
    return real_check_varname(varname, 1);
}

int check_identifier (const char *varname)
{
    /* FIXME series case? */
    return real_check_varname(varname, 0);
}

static int bad_date_string (const char *s)
{
    int err = 0;

    gretl_error_clear();

    while (*s && !err) {
	if (!isdigit((unsigned char) *s) && !IS_DATE_SEP(*s)) {
	    if (isprint((unsigned char) *s)) {
		gretl_errmsg_sprintf(_("Bad character '%c' in date string"), *s);
	    } else {
		gretl_errmsg_sprintf(_("Bad character %d in date string"), *s);
	    }
	    err = 1;
	}
	s++;
    }

    return err;
}

static void maybe_unquote_label (char *targ, const char *src)
{
    if (*src == '"' || *src == '\'') {
	int n;

	strcpy(targ, src + 1);
	n = strlen(targ);
	if (n > 0 && (targ[n-1] == '"' || targ[n-1] == '\'')) {
	    targ[n-1] = '\0';
	}
    } else {
	strcpy(targ, src);
    }
}

static int get_dot_pos (const char *s)
{
    int i, pos = 0;

    for (i=0; *s != '\0'; i++, s++) {
	if (IS_DATE_SEP(*s)) {
	    pos = i;
	    break;
	}
    }

    return pos;
}

static int match_obs_marker (const char *s, const DATASET *dset)
{
    char test[OBSLEN];
    int t;

#if DATES_DEBUG
    fprintf(stderr, "dateton: checking '%s' against marker strings\n", s);
#endif

    maybe_unquote_label(test, s);

    for (t=0; t<dset->n; t++) {
	if (!strcmp(test, dset->S[t])) {
	    /* handled */
	    return t;
	}
    }

    return -1;
}

static int caldate_precheck (const char *s, int *y2, int *slashed)
{
    int n = strlen(s);
    int ok = 0;

    if (n < 8 || !isdigit(s[0]) || !isdigit(s[1])) {
	/* can't be date */
	return 0;
    }

    if (n == 10) {
	if (s[4] == '-' && s[7] == '-') {
	    ok = 1; /* could be YYYY-MM-DD */
	} else if (s[4] == '/' && s[7] == '/') {
	    *slashed = 1; /* could be YYYY/MM/DD */
	    ok = 1;
	}
    } else if (n == 8) {
	if (s[2] == '-' && s[5] == '-') {
	    ok = 1; /* could be YY-MM-DD */
	} else if (s[2] == '/' && s[5] == '/') {
	    *slashed = 1; /* could be YY/MM/DD */
	    ok = 1;
	}
	if (ok) *y2 = 1;
    }

    return ok;
}

static int datecmp (const char *s1, int y21, int slash1,
		    const char *s2, int y22, int slash2)
{
    if (y21 && !y22) {
	s2 += 2;
    } else if (!y21 && y22) {
	s1 += 2;
    }
    if (slash1 && !slash2) {
	char revs1[OBSLEN];

	strcpy(revs1, s1);
	gretl_charsub(revs1, '/', '-');
	return strcmp(revs1, s2);
    } else if (slash2 && !slash1) {
	char revs2[OBSLEN];

	strcpy(revs2, s2);
	gretl_charsub(revs2, '/', '-');
	return strcmp(s1, revs2);
    } else {
	return strcmp(s1, s2);
    }
}

static int
real_dateton (const char *date, const DATASET *dset, int nolimit)
{
    int handled = 0;
    int t, n = -1;

    /* first check if this is calendar data and if so,
       treat accordingly */

    if (calendar_data(dset)) {
#if DATES_DEBUG
	fprintf(stderr, "dateton: treating as calendar data\n");
#endif
	if (dataset_has_markers(dset)) {
	    /* "hard-wired" calendar dates as strings */
	    int tryit, y21 = 0, slash1 = 0;
	    int y22 = 0, slash2 = 0;

	    tryit = caldate_precheck(date, &y21, &slash1);
	    if (!tryit) {
		return -1;
	    }
	    tryit = caldate_precheck(dset->S[0], &y22, &slash2);
	    if (!tryit) {
		return -1;
	    }
	    for (t=0; t<dset->n; t++) {
		if (!datecmp(date, y21, slash1, dset->S[t], y22, slash2)) {
		    /* handled */
		    return t;
		}
	    }
	    /* out of options: abort */
	    return -1;
	} else {
	    /* automatic calendar dates */
	    n = calendar_obs_number(date, dset, nolimit);
	    handled = 1;
	}
    } else if (dataset_is_daily(dset) ||
	       dataset_is_weekly(dset)) {
#if DATES_DEBUG
	fprintf(stderr, "dateton: trying undated time series\n");
#endif
	t = positive_int_from_string(date);
	if (t > 0) {
	    n = t - 1;
	    handled = 1;
	}
    } else if (dataset_is_decennial(dset)) {
	t = positive_int_from_string(date);
	if (t > 0) {
	    n = (t - dset->sd0) / 10;
	    handled = 1;
	}
    } else if (dataset_has_markers(dset)) {
	t = match_obs_marker(date, dset);
	if (t >= 0) {
	    return t;
	}
	/* else maybe just a straight obs number */
	t = positive_int_from_string(date);
	if (t > 0) {
	    n = t - 1;
	    handled = 1;
	}
    }

    if (!handled) {
	int pos1, pos2;

#if DATES_DEBUG
	fprintf(stderr, "dateton: treating %s as regular numeric obs\n",
		date);
#endif
	if (bad_date_string(date)) {
	    return -1;
	}

	pos1 = get_dot_pos(date);
	pos2 = get_dot_pos(dset->stobs);

	if ((pos1 && !pos2) || (pos2 && !pos1)) {
	    gretl_errmsg_sprintf(_("'%s': invalid observation index"),
				 date);
	} else if (!pos1 && !pos2) {
	    n = atoi(date) - atoi(dset->stobs);
	} else if (pos1 > OBSLEN - 2) {
	    gretl_errmsg_sprintf(_("'%s': invalid observation index"),
				 date);
	} else {
	    char tmp[OBSLEN];
	    int maj, min;

	    *tmp = '\0';
	    strncat(tmp, date, OBSLEN-1);
	    tmp[pos1] = '\0';
	    maj = positive_int_from_string(tmp);
	    min = positive_int_from_string(tmp + pos1 + 1);

	    if (maj <= 0 || min <= 0 || min > dset->pd) {
		gretl_errmsg_sprintf(_("'%s': invalid observation index"),
				     date);
		n = -1;
	    } else {
		int maj0, min0;

		*tmp = '\0';
		strncat(tmp, dset->stobs, OBSLEN-1);
		tmp[pos2] = '\0';
		maj0 = atoi(tmp);
		min0 = atoi(tmp + pos2 + 1);

		n = dset->pd * (maj - maj0) + (min - min0);
	    }
	}
    }

    if (!nolimit && dset->n > 0 && n >= dset->n) {
	fprintf(stderr, "n = %d, dset->n = %d: out of bounds\n", n, dset->n);
	gretl_errmsg_set(_("Observation number out of bounds"));
	n = -1;
    }

#if DATES_DEBUG
    fprintf(stderr, "dateton: returning %d\n", n);
#endif

    return n;
}

/**
 * dateton:
 * @date: string representation of date for processing.
 * @dset: pointer to data information struct.
 *
 * Determines the observation number corresponding to @date,
 * relative to @dset. It is an error if @date represents an
 * observation that lies outside of the full data range
 * specified in @dset.
 *
 * Returns: zero-based observation number, or -1 on error.
 */

int dateton (const char *date, const DATASET *dset)
{
    return real_dateton(date, dset, 0);
}

/**
 * merge_dateton:
 * @date: string representation of date for processing.
 * @dset: pointer to data information struct.
 *
 * Works just as dateton(), except that for this function it
 * is not an error if @date represents an observation that
 * lies beyond the data range specified in @dset. This is
 * intended for use when merging data, or when creating a new
 * dataset.
 *
 * Returns: zero-based observation number, or -1 on error.
 */

int merge_dateton (const char *date, const DATASET *dset)
{
    return real_dateton(date, dset, 1);
}

static char *panel_obs (char *s, int t, const DATASET *dset)
{
    int i = t / dset->pd + 1;
    int j = (t + 1) % dset->pd;
    int d = 1 + floor(log10(dset->pd));

    if (j == 0) {
	j = dset->pd;
    }

    sprintf(s, "%d:%0*d", i, d, j);

    return s;
}

/**
 * ntolabel:
 * @datestr: char array to which date is to be printed.
 * @t: zero-based observation number.
 * @dset: data information struct.
 *
 * Prints to @datestr (which must be at least #OBSLEN bytes)
 * the calendar representation of observation number @t.
 *
 * Returns: the observation string.
 */

char *ntolabel (char *datestr, int t, const DATASET *dset)
{
    double x;

#if 0
    fprintf(stderr, "ntolabel: t=%d, pd=%d, sd0=%g, incoming stobs='%s'\n",
	    t, dset->pd, dset->sd0, dset->stobs);
    fprintf(stderr, " calendar_data(dset) %d\n", calendar_data(dset));
#endif

    if (calendar_data(dset)) {
	/* handles both daily and dated weekly data */
	if (dataset_has_markers(dset)) {
	    strcpy(datestr, dset->S[t]);
	    if (strchr(datestr, '/')) {
		gretl_charsub(datestr, '/', '-');
	    }
	} else {
	    calendar_date_string(datestr, t, dset);
	}
	return datestr;
    } else if (dataset_is_daily(dset) ||
	       dataset_is_weekly(dset)) {
	/* undated time series */
	x = date_as_double(t, 1, dset->sd0);
	sprintf(datestr, "%d", (int) x);
	return datestr;
    } else if (dataset_is_decennial(dset)) {
	x = dset->sd0 + 10 * t;
	sprintf(datestr, "%d", (int) x);
	return datestr;
    } else if (dataset_is_panel(dset)) {
	panel_obs(datestr, t, dset);
	return datestr;
    }

    x = date_as_double(t, dset->pd, dset->sd0);

    if (dset->pd == 1) {
        sprintf(datestr, "%d", (int) x);
    } else {
	int pdp = dset->pd;
	short len = 1;
	char fmt[10];

	while ((pdp = pdp / 10)) len++;
	sprintf(fmt, "%%.%df", len);
	sprintf(datestr, fmt, x);
	colonize_obs(datestr);
    }

    return datestr;
}

/* print observation date in ISO 8601 extended format */

char *ntolabel_8601 (char *datestr, int t, const DATASET *dset)
{
    *datestr = '\0';

    if (calendar_data(dset)) {
	if (dataset_has_markers(dset)) {
	    strcpy(datestr, dset->S[t]);
	} else {
	    calendar_date_string(datestr, t, dset);
	}
    } else if (dataset_is_decennial(dset)) {
	double x = dset->sd0 + 10 * t;
	int yr = lrint(x);

	sprintf(datestr, "%d-01-01", yr);
    } else {
	double x = date_as_double(t, dset->pd, dset->sd0);
	int maj = lrint(floor(x));

	if (dset->pd == 1) {
	    sprintf(datestr, "%d-01-01", maj);
	} else if (dset->pd == 12) {
	    int min = lrint(100 * (x - floor(x)));

	    sprintf(datestr, "%d-%02d-01", maj, min);
	} else if (dset->pd == 4) {
	    int min = lrint(10 * (x - floor(x)));
	    int mo = min==2 ? 4 : min==3? 7 : min==4? 10 : min;

	    sprintf(datestr, "%d-%02d-01", maj, mo);
	}
    }

    return datestr;
}

#define xround(x) (((x-floor(x))>.5)? ceil(x) : floor(x))

/**
 * get_subperiod:
 * @t: zero-based observation number.
 * @dset: data information struct.
 * @err: location to receive error code, or NULL.
 *
 * For "seasonal" time series data (in a broad sense),
 * determines the sub-period at observation @t. The "sub-period"
 * might be a quarter, month, hour or whatever.  The value
 * returned is zero-based (e.g. first quarter = 0).
 * If the data are not "seasonal", 0 is returned and if
 * @err is non-NULL it receives a non-zero error code.
 *
 * Returns: the sub-period.
 */

int get_subperiod (int t, const DATASET *dset, int *err)
{
    int ret = 0;

    if (!dataset_is_seasonal(dset)) {
	if (err != NULL) {
	    *err = E_PDWRONG;
	}
	return 0;
    }

    if (dataset_is_weekly(dset)) {
	/* bodge -- what else to do? */
	ret = t % dset->pd;
    } else if (dated_daily_data(dset)) {
	char obs[OBSLEN];

	ntolabel(obs, t, dset);
	ret = weekday_from_date(obs);
        /* 2023-08-06: for now, ensure backward compatibility,
           though I'm not sure that's actually right.
        */
        ret = ret == 7 ? 0 : ret;
    } else if (dataset_is_daily(dset)) {
	/* bodge, again */
	ret = t % dset->pd;
    } else {
	/* quarterly, monthly, hourly... */
	double x = date_as_double(t, dset->pd, dset->sd0);
	int i, d = ceil(log10(dset->pd));

	x -= floor(x);
	for (i=0; i<d; i++) {
	    x *= 10;
	}
	ret = xround(x) - 1;
    }

    return ret;
}

/**
 * get_precision:
 * @x: data array.
 * @n: length of @x.
 * @placemax: the maximum number of decimal places to try.
 *
 * Find the number of decimal places required to represent a given
 * data series uniformly and accurately, if possible.
 *
 * Returns: the required number of decimal places or
 * #PMAX_NOT_AVAILABLE if it can't be done.
 */

int get_precision (const double *x, int n, int placemax)
{
    int t, p, pmax = 0;
    char *s, numstr[64];
    int len, n_ok = 0;
    double z;

    for (t=0; t<n; t++) {
	if (!na(x[t])) {
	    z = fabs(x[t]);
	    /* escape clause: numbers are too big or too small for
	       this treatment */
	    if (z > 0 && (z < 1.0e-6 || z > 1.0e+8)) {
		return PMAX_NOT_AVAILABLE;
	    }
	    n_ok++;
	}
    }

    if (n_ok == 0) {
	return PMAX_NOT_AVAILABLE;
    }

    for (t=0; t<n; t++) {
	if (!na(x[t])) {
	    p = placemax;
	    sprintf(numstr, "%.*f", p, fabs(x[t]));
	    /* go to the end and drop trailing zeros */
	    len = strlen(numstr);
	    s = numstr + len - 1;
	    while (*s-- == '0') {
		p--;
		len--;
	    }
	    if (len > 10) {
		/* this is going to be too big */
		return PMAX_NOT_AVAILABLE;
	    }
	    if (p > pmax) {
		pmax = p;
	    }
	}
    }

    return pmax;
}

struct extmap {
    GretlFileType ftype;
    const char *ext;
};

static struct extmap data_ftype_map[] = {
    { GRETL_XML_DATA,     ".gdt" },
    { GRETL_BINARY_DATA,  ".gdtb" },
    { GRETL_CSV,          ".csv" },
    { GRETL_OCTAVE,       ".m" },
    { GRETL_GNUMERIC,     ".gnumeric" },
    { GRETL_XLS,          ".xls" },
    { GRETL_XLSX,         ".xlsx" },
    { GRETL_ODS,          ".ods" },
    { GRETL_WF1,          ".wf1" },
    { GRETL_DTA,          ".dta" },
    { GRETL_SAV,          ".sav" },
    { GRETL_SAS,          ".xpt" }
};

static const char *map_suffixes[] = {
    ".json",
    ".geojson",
    ".dbf",
    ".shp"
};

static const char *get_filename_extension (const char *fname)
{
    const char *ext = strrchr(fname, '.');

    if (ext != NULL && strchr(ext, '/')) {
	/* the rightmost dot is not in the basename */
	ext = NULL;
    }

#ifdef WIN32
    if (ext != NULL && strchr(ext, '\\')) {
	ext = NULL;
    }
#endif

    return ext;
}

static GretlFileType dat_filetype (const char *fname)
{
    FILE *fp = gretl_fopen(fname, "rb");
    GretlFileType ret = GRETL_CSV;

    if (fp != NULL) {
	char test[3] = {0};

	if (fread(test, 1, 2, fp) == 2 && !strcmp(test, "/*")) {
	    ret = GRETL_JMULTI;
	}
	fclose(fp);
    }

    return ret;
}

static GretlFileType data_file_type_from_extension (const char *ext,
						    const char *fname)
{
    int i, n = G_N_ELEMENTS(data_ftype_map);

    for (i=0; i<n; i++) {
	if (!g_ascii_strcasecmp(ext, data_ftype_map[i].ext)) {
	    return data_ftype_map[i].ftype;
	}
    }

    /* try to handle JMulTi */
    if (!g_ascii_strcasecmp(ext, ".dat")) {
	return dat_filetype(fname);
    }

    /* a few extras */
    if (!g_ascii_strcasecmp(ext, ".txt") ||
	!g_ascii_strcasecmp(ext, ".asc") ||
	!g_ascii_strcasecmp(ext, ".gz")) {
	return GRETL_CSV;
    }

    /* map metadata */
    n = G_N_ELEMENTS(map_suffixes);
    for (i=0; i<n; i++) {
	if (!g_ascii_strcasecmp(ext, map_suffixes[i])) {
	    return GRETL_MAP;
	}
    }

    return GRETL_UNRECOGNIZED;
}

GretlFileType data_file_type_from_name (const char *fname)
{
    const char *ext = strrchr(fname, '.');

    if (ext != NULL && strchr(ext, '/')) {
	/* the rightmost dot is not in the basename */
	ext = NULL;
    }

#ifdef WIN32
    if (ext != NULL && strchr(ext, '\\')) {
	ext = NULL;
    }
#endif

    if (ext != NULL) {
	return data_file_type_from_extension(ext, fname);
    }

    return GRETL_UNRECOGNIZED;
}

#define non_native(o) (o & (OPT_M | OPT_R | OPT_D | OPT_G | OPT_J))

static GretlDataFormat
format_from_opt_or_name (gretlopt opt, const char *fname,
			 char *delim, int *add_ext,
			 int *gzip, int *err)
{
    GretlDataFormat fmt = GRETL_FMT_GDT;

    if (has_suffix(fname, ".gdt")) {
	if (non_native(opt)) {
	    *err = E_BADOPT;
	}
	return GRETL_FMT_GDT;
    } else if (has_suffix(fname, ".gdtb")) {
	if (non_native(opt)) {
	    *err = E_BADOPT;
	}
	return GRETL_FMT_BINARY;
    } else if (has_suffix(fname, ".geojson")) {
	return GRETL_FMT_JSON;
    }

    if (opt & OPT_M) {
	fmt = GRETL_FMT_OCTAVE;
    } else if (opt & OPT_R) {
	fmt = GRETL_FMT_R;
    } else if (opt & OPT_D) {
	fmt = GRETL_FMT_DB;
    } else if (opt & OPT_G) {
	fmt = GRETL_FMT_DAT;
    } else if (opt & OPT_J) {
	fmt = GRETL_FMT_JM;
    }

    if (fmt == GRETL_FMT_GDT) {
	if (has_suffix(fname, ".R")) {
	    fmt = GRETL_FMT_R;
	} else if (has_suffix(fname, ".csv")) {
	    fmt = GRETL_FMT_CSV;
	} else if (has_suffix(fname, ".csv.gz")) {
	    fmt = GRETL_FMT_CSV;
	    *gzip = 1;
	} else if (has_suffix(fname, ".m")) {
	    fmt = GRETL_FMT_OCTAVE;
	} else if (has_suffix(fname, ".txt") ||
		   has_suffix(fname, ".asc")) {
	    fmt = GRETL_FMT_CSV;
	    *delim = ' ';
	} else if (has_suffix(fname, ".dta")) {
	    fmt = GRETL_FMT_DTA;
	} else if (has_suffix(fname, ".bin")) {
	    fmt = GRETL_FMT_DB;
	}
    }

    if (fmt == GRETL_FMT_GDT) {
	*add_ext = 1;
    }

    return fmt;
}

void date_maj_min (int t, const DATASET *dset, int *maj, int *min)
{
    char obs[OBSLEN];

    ntolabel(obs, t, dset);

    if (maj != NULL) {
	*maj = atoi(obs);
    }

    if (min != NULL) {
	char *s, sep = ':';

	if (strchr(obs, sep) == NULL) {
	    if (dset->pd == 4 && strchr(obs, 'Q')) {
		sep = 'Q';
	    } else if (dset->pd == 12 && strchr(obs, 'M')) {
		sep = 'M';
	    }
	}

	s = strchr(obs, sep);
	if (s != NULL && strlen(s) > 1) {
	    *min = atoi(s + 1);
	} else {
	    *min = 1;
	}
    }
}

#define NO_PMAX(p,k) (p == NULL || p[k-1] == PMAX_NOT_AVAILABLE)

#define TMPLEN 64

static void csv_data_out (const DATASET *dset, const int *list,
			  int print_obs, int digits, char decpoint,
			  char delim, FILE *fp, gzFile fz)
{
    const char *NA = get_csv_na_write_string();
    char tmp[TMPLEN];
    double xt;
    int popit = 0, dotsub = 0;
    int t, i, vi;

    if (decpoint == '.' && get_local_decpoint() == ',') {
	gretl_push_c_numeric_locale();
	popit = 1;
    } else if (decpoint == ',' && get_local_decpoint() == '.') {
	dotsub = 1;
    }

    for (t=dset->t1; t<=dset->t2; t++) {
	if (print_obs) {
	    if (dset->S != NULL) {
		if (fz != NULL) {
		    gzprintf(fz, "\"%s\"%c", dset->S[t], delim);
		} else {
		    fprintf(fp, "\"%s\"%c", dset->S[t], delim);
		}
	    } else {
		ntolabel(tmp, t, dset);
		if (quarterly_or_monthly(dset)) {
		    modify_date_for_csv(tmp, dset->pd);
		}
		if (fz != NULL) {
		    gzprintf(fz, "%s%c", tmp, delim);
		} else {
		    fprintf(fp, "%s%c", tmp, delim);
		}
	    }
	}

	for (i=1; i<=list[0]; i++) {
	    vi = list[i];
	    xt = dset->Z[vi][t];
	    if (na(xt)) {
		fputs(NA, fp);
	    } else {
		if (is_string_valued(dset, vi)) {
		    const char *st;

		    st = series_get_string_for_obs(dset, vi, t);
		    if (st != NULL) {
			*tmp = '\0';
			strcat(tmp, "\"");
			strncat(tmp, st, TMPLEN - 3);
			strcat(tmp, "\"");
		    } else {
			fprintf(stderr, "missing string at t=%d, vi=%d, xt=%g\n",
				t, vi, xt);
			strcpy(tmp, "\"NA\"");
		    }
		} else if (series_is_coded(dset, vi)) {
		    sprintf(tmp, "\"%d\"", (int) xt);
		} else {
		    sprintf(tmp, "%.*g", digits, xt);
		}
		if (dotsub) {
		    gretl_charsub(tmp, '.', ',');
		}
		if (fz != NULL) {
		    gzputs(fz, tmp);
		} else {
		    fputs(tmp, fp);
		}
	    }
	    if (fz != NULL) {
		gzputc(fz, i < list[0] ? delim : '\n');
	    } else {
		fputc(i < list[0] ? delim : '\n', fp);
	    }
	}
    }

    if (popit) {
	gretl_pop_c_numeric_locale();
    }
}

static int markers_are_unique (const DATASET *dset)
{
    int t, s;

    for (t=dset->t1; t<dset->t2; t++) {
	for (s=t+1; s<=dset->t2; s++) {
	    if (strcmp(dset->S[t], dset->S[s]) == 0) {
		return 0;
	    }
	}
    }

    return 1;
}

static void R_data_out (const DATASET *dset, const int *list,
			int digits, FILE *fp)
{
    int print_markers = 0;
    double xt;
    int t, i, vi;

    if (dset->S != NULL) {
	print_markers = markers_are_unique(dset);
    }

    for (t=dset->t1; t<=dset->t2; t++) {
	if (print_markers) {
	    fprintf(fp, "\"%s\" ", dset->S[t]);
	}
	for (i=1; i<=list[0]; i++) {
	    vi = list[i];
	    xt = dset->Z[vi][t];
	    if (na(xt)) {
		fputs("NA", fp);
	    } else if (is_string_valued(dset, vi)) {
		fprintf(fp, "\"%s\"", series_get_string_for_obs(dset, vi, t));
	    } else if (series_is_coded(dset, vi)) {
		fprintf(fp, "\"%d\"", (int) xt);
	    } else {
		fprintf(fp, "%.*g", digits, xt);
	    }
	    fputc(i < list[0] ? ' ' : '\n', fp);
	}
    }
}

static int write_dta_data (const char *fname, const int *list,
			   gretlopt opt, const DATASET *dset)
{
    int (*exporter) (const char *, const int *, gretlopt,
		     const DATASET *);
    int err = 0;

    exporter = get_plugin_function("stata_export");

    if (exporter == NULL) {
        err = 1;
    } else {
	err = (*exporter)(fname, list, opt, dset);
    }

    return err;
}

static int write_map_data (const char *fname,
			   const int *list,
			   const DATASET *dset)
{
    gretl_bundle *b = NULL;
    int err = 0;

    b = get_current_map(dset, list, &err);

    if (!err) {
	err = gretl_bundle_write_to_file(b, fname, 0);
    }

    gretl_bundle_destroy(b);

    return err;
}

#define DEFAULT_CSV_DIGITS 15

static int real_write_data (const char *fname, int *list,
			    const DATASET *dset, gretlopt opt,
			    int progress, PRN *prn)
{
    int i, t, v, l0;
    GretlDataFormat fmt;
    char datfile[MAXLEN];
    int n = dset->n;
    int pop_locale = 0;
    char delim = 0;
    FILE *fp = NULL;
    gzFile fz = NULL;
    int freelist = 0;
    int csv_digits = 0;
    int add_ext = 0;
    int gzip = 0;
    double xx;
    int err = 0;

    gretl_error_clear();

    if (list != NULL && list[0] == 0) {
	return E_ARGS;
    }

    fmt = format_from_opt_or_name(opt, fname, &delim, &add_ext,
				  &gzip, &err);
    if (err) {
	return err;
    }

    if (list == NULL) {
	list = full_var_list(dset, &l0);
	if (l0 == 0) {
	    return E_ARGS;
	} else if (list == NULL) {
	    return E_ALLOC;
	} else {
	    freelist = 1;
	}
    }

    l0 = list[0];
    fname = gretl_maybe_switch_dir(fname);

    if (fmt == GRETL_FMT_GDT || fmt == GRETL_FMT_BINARY) {
	/* write native data file (.gdt or .gdtb) */
	err = gretl_write_gdt(fname, list, dset, opt, progress);
	goto write_exit;
    }

    if (fmt == GRETL_FMT_DB) {
	/* native type database file */
	err = write_db_data(fname, list, opt, dset);
	goto write_exit;
    }

    if (fmt == GRETL_FMT_DTA) {
	/* Stata */
	err = write_dta_data(fname, list, opt, dset);
	goto write_exit;
    }

    if (fmt == GRETL_FMT_JSON) {
	/* writing map as geojson */
	err = write_map_data(fname, list, dset);
	goto write_exit;
    }

    strcpy(datfile, fname);

    /* open file for output */
    if (gzip) {
	fz = gretl_gzopen(datfile, "wb1");
    } else {
	fp = gretl_fopen(datfile, "wb");
    }
    if (fp == NULL && fz == NULL) {
	err = E_FOPEN;
	goto write_exit;
    }

    csv_digits = libset_get_int(CSV_DIGITS);

    if (csv_digits <= 0) {
	csv_digits = DEFAULT_CSV_DIGITS;
    }

    if (fmt != GRETL_FMT_CSV) {
	/* ensure C locale for data output */
	gretl_push_c_numeric_locale();
	pop_locale = 1;
    }

    if (fmt == GRETL_FMT_CSV) {
	const char *msg = get_optval_string(STORE, OPT_E);
	char decpoint = get_data_export_decpoint();
	int print_obs = 0;

	if (opt & OPT_I) {
	    /* the CSV --decimal-comma option */
	    decpoint = ',';
	    delim = ';';
	} else if (delim == 0) {
	    delim = get_data_export_delimiter();
	}

	if (msg != NULL && *msg != '\0') {
	    if (gzip) {
		gzprintf(fz, "# %s\n", msg);
	    } else {
		fprintf(fp, "# %s\n", msg);
	    }
	}

	if (!(opt & OPT_X)) {
	    /* OPT_X prohibits printing of observation strings */
	    print_obs = dataset_is_time_series(dset) || dset->S != NULL;
	}

	if (!(opt & OPT_N)) {
	    /* header: variable names */
	    if (print_obs && (dset->S != NULL || dset->structure != CROSS_SECTION)) {
		if (gzip) {
		    gzprintf(fz, "obs%c", delim);
		} else {
		    fprintf(fp, "obs%c", delim);
		}
	    }
	    for (i=1; i<l0; i++) {
		if (gzip) {
		    gzprintf(fz, "%s%c", dset->varname[list[i]], delim);
		} else {
		    fprintf(fp, "%s%c", dset->varname[list[i]], delim);
		}
	    }
	    if (gzip) {
		gzprintf(fz, "%s\n", dset->varname[list[l0]]);
	    } else {
		fprintf(fp, "%s\n", dset->varname[list[l0]]);
	    }
	}

	csv_data_out(dset, list, print_obs, csv_digits,
		     decpoint, delim, fp, fz);
    } else if (fmt == GRETL_FMT_R) {
	/* friendly to GNU R */
	if (dataset_is_time_series(dset)) {
	    char datestr[OBSLEN];

	    ntolabel(datestr, dset->t1, dset);
	    fprintf(fp, "# time-series data: start = %s, frequency = %d\n",
		    datestr, dset->pd);
	}

	for (i=1; i<l0; i++) {
	    fprintf(fp, "%s ", dset->varname[list[i]]);
	}
	fprintf(fp, "%s\n", dset->varname[list[l0]]);

	R_data_out(dset, list, csv_digits, fp);
    } else if (fmt == GRETL_FMT_OCTAVE) {
	/* GNU Octave: write out data as several matrices (one per
	   series) in the same file */

	for (i=1; i<=list[0]; i++) {
	    v = list[i];
	    fprintf(fp, "# name: %s\n# type: matrix\n# rows: %d\n# columns: 1\n",
		    dset->varname[v], n);
	    for (t=dset->t1; t<=dset->t2; t++) {
		xx = dset->Z[v][t];
		if (na(xx)) {
		    fputs("NaN ", fp);
		} else {
		    fprintf(fp, "%.*g\n", csv_digits, xx);
		}
	    }
	}
    } else if (fmt == GRETL_FMT_DAT) {
	/* PcGive: data file with load info */
	int pd = dset->pd;

	for (i=1; i<=list[0]; i++) {
	    fprintf(fp, ">%s ", dset->varname[list[i]]);
	    if (dset->structure == TIME_SERIES &&
		(pd == 1 || pd == 4 || pd == 12)) {
		int maj, min;

		date_maj_min(dset->t1, dset, &maj, &min);
		fprintf(fp, "%d %d ", maj, min);
		date_maj_min(dset->t2, dset, &maj, &min);
		fprintf(fp, "%d %d %d", maj, min, pd);
	    } else {
		fprintf(fp, "%d 1 %d 1 1", dset->t1, dset->t2);
	    }

	    fputc('\n', fp);

	    for (t=dset->t1; t<=dset->t2; t++) {
		v = list[i];
		xx = dset->Z[v][t];
		if (na(xx)) {
		    fprintf(fp, "-9999.99");
		} else {
		    fprintf(fp, "%.*g", csv_digits, xx);
		}
		fputc('\n', fp);
	    }
	    fputc('\n', fp);
	}
    } else if (fmt == GRETL_FMT_JM) {
	/* JMulti: ascii with comments and date info */
	const char *vlabel;
	int maj, min;

	fputs("/*\n", fp);
	for (i=1; i<=list[0]; i++) {
	    v = list[i];
	    vlabel = series_get_label(dset, v);
	    fprintf(fp, " %s: %s\n", dset->varname[v],
		    vlabel == NULL ? "" : vlabel);
	}
	fputs("*/\n", fp);
	date_maj_min(dset->t1, dset, &maj, &min);
	if (dset->pd == 4 || dset->pd == 12) {
	    fprintf(fp, "<%d %c%d>\n", maj, (dset->pd == 4)? 'Q' : 'M', min);
	} else if (dset->pd == 1) {
	    fprintf(fp, "<%d>\n", maj);
	} else {
	    fputs("<1>\n", fp);
	}
	for (i=1; i<=list[0]; i++) {
	    v = list[i];
	    fprintf(fp, " %s", dset->varname[v]);
	}
	fputc('\n', fp);
	for (t=dset->t1; t<=dset->t2; t++) {
	    for (i=1; i<=list[0]; i++) {
		v = list[i];
		if (na(dset->Z[v][t])) {
		    fputs("NaN ", fp);
		} else {
		    fprintf(fp, "%.*g ", csv_digits, dset->Z[v][t]);
		}
	    }
	    fputc('\n', fp);
	}
    }

    if (pop_locale) {
	gretl_pop_c_numeric_locale();
    }

    if (fz != NULL) {
	gzclose(fz);
    } else if (fp != NULL) {
	fclose(fp);
    }

 write_exit:

    if (!err && prn != NULL && gretl_messages_on()) {
	if (add_ext) {
	    pprintf(prn, _("wrote %s.gdt\n"), fname);
	} else {
	    pprintf(prn, _("wrote %s\n"), fname);
	}
    }

    if (freelist) {
	free(list);
    }

    return err;
}

/**
 * write_data:
 * @fname: name of file to write.
 * @list: list of variables to write (or %NULL to write all series).
 * @dset: dataset struct.
 * @opt: option flag indicating format in which to write the data.
 * @prn: gretl printer or NULL.
 *
 * Write out a data file containing the values of the given set
 * of variables.
 *
 * Returns: 0 on successful completion, non-zero on error.
 */

int write_data (const char *fname, int *list, const DATASET *dset,
		gretlopt opt, PRN *prn)
{
    return real_write_data(fname, list, dset, opt, 0, prn);
}

int gui_write_data (const char *fname, int *list, const DATASET *dset,
		    gretlopt opt)
{
    return real_write_data(fname, list, dset, opt, 1, NULL);
}

/**
 * is_gzipped:
 * @fname: filename to examine.
 *
 * Determine if the given file is gzipped.
 *
 * Returns: 1 in case of a gzipped file, 0 if not gzipped or
 * inaccessible.
 *
 */

int is_gzipped (const char *fname)
{
    FILE *fp;
    int gz = 0;

    if (fname == NULL || *fname == '\0') {
	return 0;
    }

    fp = gretl_fopen(fname, "rb");
    if (fp == NULL) {
	return 0;
    }

    if (fgetc(fp) == 037 && fgetc(fp) == 0213) {
	gz = 1;
    }

    fclose(fp);

    return gz;
}

/**
 * gretl_seek_data:
 * @fname: name of file to try.
 * @dset: dataset struct.
 * @opt: option flags.
 * @prn: where messages should be written.
 *
 * Reads data from file (either "native" or "CSV") into
 * gretl's workspace, allocating storage as required.
 *
 * Path-searching is performed on @fname, and the function
 * will try adding the .gdt or .gdtb extension to @fname
 * if this is not given.
 *
 * The only applicable option is that @opt may contain
 * OPT_T when appending data to a panel dataset: in
 * that case we try to interpret the new data as time
 * series, in common across all panel units. In most
 * cases, just give OPT_NONE.
 *
 * Returns: 0 on successful completion, non-zero otherwise.
 */

int gretl_seek_data (char *fname, DATASET *dset,
                     gretlopt opt, PRN *prn)
{
    gretlopt myopt = OPT_NONE;
    char *test;
    int err = 0;

    gretl_error_clear();

#if 0
    fprintf(stderr, "gretl_seek_data: calling addpath\n");
#endif

    test = gretl_addpath(fname, 0);
    if (test == NULL) {
	return E_FOPEN;
    }

    if (opt & OPT_T) {
	myopt = OPT_T;
    }

    if (has_native_data_suffix(fname)) {
	/* specific processing for gretl datafiles  */
	err = gretl_read_gdt(fname, dset, myopt, prn);
    } else {
	/* try fallback to a "csv"-type import */
	err = import_csv(fname, dset, myopt, prn);
    }

    return err;
}

/**
 * gretl_get_data:
 * @fname: name of file to try.
 * @dset: dataset struct.
 * @opt: option flags.
 * @prn: where messages should be written.
 *
 * Works like gretl_seek_data() except that no path searching
 * is performed; @fname is always taken "as is" and is never
 * modified.
 *
 * Returns: 0 on successful completion, non-zero otherwise.
 */

int gretl_get_data (const char *fname, DATASET *dset,
		    gretlopt opt, PRN *prn)
{
    gretlopt myopt = OPT_NONE;
    int err = 0;

    gretl_error_clear();

    if (opt & OPT_T) {
	myopt = OPT_T;
    }

    if (has_native_data_suffix(fname)) {
	/* specific processing for gretl datafiles  */
	err = gretl_read_gdt(fname, dset, myopt, prn);
    } else {
        /* try for an "import" of some kind */
        GretlFileType ft = gretl_detect_filetype(fname);

        if (ft == GRETL_UNRECOGNIZED) {
            gretl_errmsg_set(_("Unknown data import type"));
            err = E_DATA;
        } else {
            err = gretl_read_foreign_data(fname, ft, dset, prn);
        }
    }

    return err;
}

/**
 * open_nulldata:
 * @dset: dataset struct.
 * @data_status: indicator for whether a data file is currently open
 * in gretl's work space (1) or not (0).
 * @length: desired length of data series.
 * @opt: may contain OPT_N to suppress addition of an index series.
 * @prn: gretl printing struct.
 *
 * Create an empty "dummy" data set, suitable for simulations.
 *
 * Returns: 0 on successful completion, non-zero otherwise.
 *
 */

int open_nulldata (DATASET *dset, int data_status, int length,
		   gretlopt opt, PRN *prn)
{
    int t;

    /* clear any existing data info */
    if (data_status) {
	clear_datainfo(dset, CLEAR_FULL);
    }

    /* dummy up the data info */
    dset->n = length;
    dset->v = (opt & OPT_N)? 1 : 2;
    dataset_obs_info_default(dset);

    if (dataset_allocate_varnames(dset)) {
	return E_ALLOC;
    }

    /* allocate data storage */
    if (allocate_Z(dset, 0)) {
	return E_ALLOC;
    }

    if (dset->v > 1) {
	/* add an index var */
	strcpy(dset->varname[1], "index");
	series_set_label(dset, 1, _("index variable"));
	for (t=0; t<dset->n; t++) {
	    dset->Z[1][t] = (double) (t + 1);
	}
    }

    if (prn != NULL && gretl_messages_on()) {
	/* print basic info */
	pprintf(prn, _("periodicity: %d, maxobs: %d\n"
		       "observations range: %s to %s\n"),
		dset->pd, dset->n, dset->stobs, dset->endobs);
    }

    /* Set sample range to entire length of data-set by default */
    dset->t1 = 0;
    dset->t2 = dset->n - 1;

    set_current_dataset(dset);

    return 0;
}

static int extend_markers (DATASET *dset, int old_n, int new_n)
{
    char **S = realloc(dset->S, new_n * sizeof *S);
    int t, err = 0;

    if (S == NULL) {
	err = 1;
    } else {
	dset->S = S;
	for (t=old_n; t<new_n && !err; t++) {
	    S[t] = malloc(OBSLEN);
	    if (S[t] == NULL) {
		err = 1;
	    }
	}
    }

    return err;
}

static void merge_error (const char *msg, PRN *prn)
{
    pputs(prn, msg);
    if (!printing_to_standard_stream(prn)) {
        /* hack to avoid revising translations */
        gchar *s = g_strndup(msg, strlen(msg) - 1);

	gretl_errmsg_set(s);
        g_free(s);
    }
}

static void merge_warning (const char *msg, PRN *prn)
{
    pprintf(prn, "%s\n", msg);
    if (!printing_to_standard_stream(prn)) {
	gretl_warnmsg_set(msg);
    }
}

static void merge_name_error (const char *objname, PRN *prn)
{
    gchar *msg;

    msg = g_strdup_printf(_("Can't replace %s with a series"), objname);
    pprintf(prn, "%s\n", msg);
    if (!printing_to_standard_stream(prn)) {
	gretl_errmsg_set(msg);
    }
    g_free(msg);
}

static int count_new_vars (const DATASET *d1, const DATASET *d2,
			   PRN *prn)
{
    const char *vname;
    int addvars = d2->v - 1;
    int i, j;

    /* We start by assuming that all the series in @d2 are new,
       then subtract those we find to be already present. We also
       check for collision between the names of series to be added and
       the names of existing objects other than series.
    */

    for (i=1; i<d2->v && addvars >= 0; i++) {
	vname = d2->varname[i];
	if (gretl_is_user_var(vname)) {
	    merge_name_error(vname, prn);
	    addvars = -1;
	} else if (gretl_function_depth() > 0) {
	    if (current_series_index(d1, vname) > 0) {
		addvars--;
	    }
	} else {
	    for (j=1; j<d1->v; j++) {
		if (!strcmp(vname, d1->varname[j])) {
		    addvars--;
		    break;
		}
	    }
	}
    }

#if MERGE_DEBUG
    if (gretl_function_depth() == 0) {
	int found;

	for (i=1; i<d2->v; i++) {
	    found = 0;
	    for (j=1; j<d1->v && !found; j++) {
		if (!strcmp(d2->varname[i], d1->varname[j])) {
		    found = 1;
		}
	    }
	    if (!found) {
		fprintf(stderr, "'%s' in import but not current dataset\n",
			d2->varname[i]);
	    }
	}
    }
#endif

    return addvars;
}

static int year_special_markers (const DATASET *dset,
				 const DATASET *addset)
{
    char *test;
    int overlap = 0;
    int i, t, err = 0;

    /* See if we can match obs markers in @addset
       against years in @dset: we'll try this if all
       the markers in addset are integer strings, at
       least some of them are within the obs range of
       @dset, and none of them are outside of the
       "sanity" range of 1 to 2500.
    */

    if (!dataset_is_time_series(dset) || dset->pd != 1) {
	return 0;
    }

    if (dset->markers || !addset->markers) {
	return 0;
    }

    errno = 0;

    for (i=0; i<addset->n; i++) {
	t = strtol(addset->S[i], &test, 10);
	if (*test || errno) {
	    errno = 0;
	    err = 1;
	    break;
	}
	if (t < 1 || t > 2500) {
	    err = 1;
	    break;
	}
	if (!overlap) {
	    t = dateton(addset->S[i], dset);
	    if (t >= 0 && t < dset->n) {
		overlap = 1;
	    }
	}
    }

    return !err && overlap;
}

static int compare_ranges (const DATASET *targ,
			   const DATASET *src,
			   int newvars,
			   int *offset,
			   int *yrspecial,
			   int *err)
{
    int ed0 = dateton(targ->endobs, targ);
    int sd1, ed1, addobs = -1;
    int range_err = 0;

    if (dataset_is_cross_section(targ) &&
	dataset_is_cross_section(src) &&
	!(targ->markers && src->markers)) {
	if (newvars == 0) {
	    if (src->markers) {
		/* pass the problem on to just_append_rows */
		return 0;
	    } else {
		/* assume the new data should be appended length-wise */
		*offset = ed0 + 1;
		return src->n;
	    }
	} else {
	    /* we've already determined that the series length in
	       @src doesn't match either the full series length or
	       the current sample range in @targ; we therefore have
	       no information with which to match rows for new
	       series
	    */
	    gretl_errmsg_set(_("append: don't know how to align the new series!"));
	    *err = E_DATA;
	    return -1;
	}
    }

    sd1 = merge_dateton(src->stobs, targ);
    ed1 = merge_dateton(src->endobs, targ);

#if DATES_DEBUG
    fprintf(stderr, "compare_ranges:\n"
	    " targ->n = %d, src->n = %d, targ->stobs = '%s'\n"
	    " src->stobs =  '%s' -> sd1 %d\n"
	    " src->endobs = '%s' -> ed1 = %d\n",
	    targ->n, src->n, targ->stobs, src->stobs,
	    sd1, src->endobs, ed1);
#endif

    if (sd1 < 0) {
	/* case: new data start earlier than old */
	if (ed1 < 0) {
	    range_err = 1;
	} else if (ed1 > ed0) {
	    range_err = 2;
	} else {
	    *offset = sd1;
	    addobs = 0;
	}
    } else if (sd1 == 0 && ed1 == ed0) {
	/* case: exact match of ranges */
	*offset = 0;
	addobs = 0;
    } else if (sd1 == 0) {
	/* case: starting obs the same */
	*offset = 0;
	if (ed1 > ed0) {
	    addobs = ed1 - ed0;
	} else {
	    addobs = 0;
	}
    } else if (sd1 == ed0 + 1) {
	/* case: new data start right after end of old */
	*offset = sd1;
	addobs = src->n;
    } else if (sd1 > 0) {
	/* case: new data start later than old */
	if (sd1 <= ed0) {
	    /* but there's some overlap */
	    *offset = sd1;
	    if (ed1 > ed0) {
		addobs = ed1 - ed0;
	    } else {
		addobs = 0;
	    }
	}
    }

    if (range_err) {
	/* try another approach? */
	*yrspecial = year_special_markers(targ, src);
	if (*yrspecial) {
	    addobs = 0;
	}
    }

    if (addobs < 0) {
	if (range_err == 1) {
	    fputs("compare_ranges: no overlap, can't merge\n", stderr);
	} else if (range_err == 2) {
	    fputs("compare ranges: new data start earlier, end later\n", stderr);
	} else {
	    fputs("compare_ranges: flagging error\n", stderr);
	}
    }

    return addobs;
}

/* Determine whether there's any overlap between the calendar
   in @addset and that in @dset. Return 0 on success (there is
   an overlap), non-zero otherwise.
*/

static int check_for_overlap (const DATASET *dset,
			      const DATASET *addset,
			      int *offset)
{
    int at1 = merge_dateton(addset->stobs, dset);
    int at2 = merge_dateton(addset->endobs, dset);

    if (!(at1 >= dset->n) && !(at2 < 0)) {
	/* OK, there must be some overlap */
	*offset = at1;
	return 0;
    } else {
	/* either the "add" data start after the original data end,
	   or they end before the originals start, no there's no
	   overlap
	*/
	gretl_errmsg_set(_("No overlap in data ranges"));
	return E_DATA;
    }
}

/* When appending data to a current panel dataset, and the length of
   the series in the new data is less than the full panel size
   (n * T), try to determine if it's OK to expand the incoming data to
   match.

   We'll say it's OK if the new series length equals the panel T: in
   that case we'll take the new data to be time-series, which should
   be replicated for each panel unit.

   A second possibility arises if the length of the new series
   equals the panel n: in that case we could treat it as a time-
   invariant characteristic of the panel unit, which should be
   replicated for each time period.  But note that if OPT_T is
   given, this second expansion is forbidden: the user has
   stipulated that the new data are time-varying.
*/

static int panel_expand_ok (DATASET *dset, DATASET *addinfo,
			    gretlopt opt)
{
    int n = dset->n / dset->pd;
    int T = dset->pd;
    int ok = 0;

    if (addinfo->n == T) {
	ok = 1;
    } else if (!(opt & OPT_T) &&
	       addinfo->n == n &&
	       addinfo->pd == 1) {
	ok = 1;
    }

    return ok;
}

static int panel_append_special (int addvars,
				 DATASET *dset,
				 DATASET *addset,
				 gretlopt opt,
				 PRN *prn)
{
    int n = dset->n / dset->pd;
    int T = dset->pd;
    int k = dset->v;
    int tsdata;
    int i, j, s, p, t;
    int err = 0;

    if (addvars > 0 && dataset_add_series(dset, addvars)) {
	merge_error(_("Out of memory!\n"), prn);
	err = E_ALLOC;
    }

    tsdata = ((opt & OPT_T) || addset->n != n);

    for (i=1; i<addset->v && !err; i++) {
	int v = series_index(dset, addset->varname[i]);

	if (v >= k) {
	    /* a new variable */
	    v = k++;
	    strcpy(dset->varname[v], addset->varname[i]);
	    copy_varinfo(dset->varinfo[v], addset->varinfo[i]);
	}

	s = 0;
	for (j=0; j<n; j++) {
	    /* loop across units */
	    for (t=0; t<T; t++) {
		/* loop across periods */
		p = (tsdata)? t : j;
		dset->Z[v][s++] = addset->Z[i][p];
	    }
	}
    }

    return err;
}

static int markers_compatible (const DATASET *d1, DATASET *d2,
			       int *offset)
{
    int ret = 0;

    if (d1->markers == 0 && d2->markers == 0) {
	*offset = d1->n;
	ret = 1;
    } else if (d1->markers == 0) {
	/* markers "on the right only": are they consecutive
	   integers starting between 1 and d1->n + 1?
	*/
	if (integer_string(d2->S[0])) {
	    int k0 = atoi(d2->S[0]);

	    if (k0 >= 1 && k0 <= d1->n + 1) {
		int i, k1;

		ret = 1;
		for (i=1; i<d2->n && ret; i++) {
		    if (!integer_string(d2->S[i])) {
			ret = 0;
		    } else if ((k1 = atoi(d2->S[i])) != k0 + 1) {
			ret = 0;
		    } else {
			k0 = k1;
		    }
		}
	    }
	    if (ret) {
		*offset = atoi(d2->S[0]) - 1;
		/* the @d2 markers have done their job -- yielding
		   an @offset value -- and they can now be trashed
		*/
		dataset_destroy_obs_markers(d2);
	    }
	}
    } else {
	/* markers on both sides: are they totally distinct? */
	int i, j;

	ret = 1;
	for (i=0; i<d2->n && ret; i++) {
	    for (j=0; j<d1->n && ret; j++) {
		if (!strcmp(d2->S[i], d1->S[j])) {
		    /* no, not totally distinct */
		    ret = 0;
		}
	    }
	}
    }

#if MERGE_DEBUG
    fprintf(stderr, " markers_compatible: ret=%d, offset=%d\n", ret, *offset);
#endif

    return ret;
}

static int
just_append_rows (const DATASET *targ, DATASET *src, int *offset)
{
    int ret = 0;

    if (targ->structure == CROSS_SECTION &&
	src->structure == CROSS_SECTION &&
	targ->sd0 == 1 && src->sd0 == 1) {
	int ok, test_offset = -1;

	ok = markers_compatible(targ, src, &test_offset);
	if (ok) {
	    /* note: we do this only if we're not adding any new
	       series: we'll append to existing series lengthwise
	       (or perhaps write data into existing rows)
	    */
	    *offset = test_offset;
	    ret = src->n - (targ->n - *offset);
	    if (ret < 0) {
		ret = 0;
	    }
	}
    }

    return ret;
}

static int simple_range_match (const DATASET *targ, const DATASET *src,
			       int *offset)
{
    int ret = 0;

    if (src->pd == 1 && src->structure == CROSS_SECTION) {
	if (src->n == targ->n) {
	    ret = 1;
	} else if (src->n == targ->t2 - targ->t1 + 1) {
	    ret = 1;
	    *offset = targ->t1;
	}
    }

    return ret;
}

static int merge_lengthen_series (DATASET *dset,
				  const DATASET *addset,
				  int addobs,
				  int offset)
{
    int i, t, new_n = dset->n + addobs;
    int err = 0;

    if (dset->markers) {
	err = extend_markers(dset, dset->n, new_n);
	if (!err) {
	    for (t=dset->n; t<new_n; t++) {
		strcpy(dset->S[t], addset->S[t-offset]);
	    }
	}
    }

    for (i=0; i<dset->v && !err; i++) {
	double *x;

	x = realloc(dset->Z[i], new_n * sizeof *x);
	if (x == NULL) {
	    err = E_ALLOC;
	    break;
	}

	for (t=dset->n; t<new_n; t++) {
	    if (i == 0) {
		x[t] = 1.0;
	    } else {
		x[t] = NADBL;
	    }
	}
	dset->Z[i] = x;
    }

    if (!err) {
	dset->n = new_n;
	ntolabel(dset->endobs, new_n - 1, dset);
	dset->t2 = dset->n - 1;
    }

    return err;
}

#if 0 /* not yet (maybe usable with DND) */

int basic_data_merge_check (const DATASET *dset,
			    DATASET *addset)
{
    int dayspecial = 0;
    int yrspecial = 0;
    int addsimple = 0;
    int addvars = 0;
    int addobs = 0;
    int offset = 0;
    int err = 0;

    /* first see how many new vars we have */
    addvars = count_new_vars(dset, addset, NULL);
    if (addvars < 0) {
	return 1;
    }

    if (dated_daily_data(dset) && dated_daily_data(addset)) {
	dayspecial = 1;
    }

    if (simple_range_match(dset, addset, &offset)) {
	addsimple = 1;
    } else if (dset->pd != addset->pd) {
	err = 1;
    }

    if (!err) {
	if (!addsimple) {
	    addobs = compare_ranges(dset, addset, addvars, &offset,
				    &yrspecial, &err);
	}
	if (!err && addobs <= 0 && addvars == 0) {
	    addobs = just_append_rows(dset, addset, &offset);
	}
    }

    if (!err && (addobs < 0 || addvars < 0)) {
	err = E_DATA;
    }

    if (!err && dset->markers != addset->markers) {
	if (addobs == 0 && addvars == 0) {
	    err = E_DATA;
	} else if (addset->n != dset->n && !yrspecial && !dayspecial) {
	    err = E_DATA;
	}
    }

    return err;
}

#endif

/* When we're merging dated weekly data, the day-of-week on which
   weeks are supposed to start should be consistent across the
   two datasets.
*/

static int check_week_start_dates (const DATASET *ds0,
				   const DATASET *ds1)
{
    int wd0 = weekday_from_epoch_day((guint32) ds0->sd0);
    int wd1 = weekday_from_epoch_day((guint32) ds1->sd0);

    if (wd1 != wd0) {
	gretl_errmsg_set("Error appending data: inconsistent dating of weeks");
	return E_DATA;
    } else {
	return 0;
    }
}

#define simple_structure(p) (p->structure == TIME_SERIES ||		\
			     p->structure == SPECIAL_TIME_SERIES ||	\
			     (p->structure == CROSS_SECTION &&		\
			      p->S == NULL))

/**
 * merge_data:
 * @dset: dataset struct.
 * @addset: dataset to be merged in.
 * @opt: may include OPT_T to force a time-series interpretation
 * when appending to a panel dataset; may include OPT_U to update
 * values of overlapping observations.
 * @prn: print struct to accept messages.
 *
 * Attempt to merge the content of a newly opened data file into
 * gretl's current working data set.
 *
 * Returns: 0 on successful completion, non-zero otherwise.
 */

static int merge_data (DATASET *dset, DATASET *addset,
		       gretlopt opt, PRN *prn)
{
    int update_overlap = (opt & OPT_U);
    int orig_n = dset->n;
    int dayspecial = 0;
    int yrspecial = 0;
    int fixsample = 0;
    int addsimple = 0;
    int addpanel = 0;
    int addvars = 0;
    int addobs = 0;
    int offset = 0;
    int err = 0;

#if MERGE_DEBUG
    debug_print_option_flags("merge_data", opt);
#endif

    /* first see how many new vars we have */
    addvars = count_new_vars(dset, addset, prn);
    if (addvars < 0) {
	return 1;
    }

#if MERGE_DEBUG
    fprintf(stderr, " new series count = %d\n", addvars);
#endif

    if (dated_daily_data(dset) && dated_daily_data(addset)) {
#if MERGE_DEBUG
	fprintf(stderr, " special: merging daily data\n");
#endif
	dayspecial = 1;
    } else if (dated_weekly_data(dset) && dated_weekly_data(addset)) {
	err = check_week_start_dates(dset, addset);
	if (err) {
	    return err;
	}
    }

    if (opt & OPT_X) {
	fixsample = 1;
    } else if (simple_range_match(dset, addset, &offset)) {
	/* we'll allow undated data to be merged with the existing
	   dateset, sideways, provided the number of observations
	   matches OK */
	addsimple = 1;
    } else if (dataset_is_panel(dset) &&
	       panel_expand_ok(dset, addset, opt)) {
	/* allow appending to panel when the number of obs matches
	   either the cross-section size or the time-series length
	*/
	addpanel = 1;
    } else if (dset->pd != addset->pd) {
	merge_error(_("Data frequency does not match\n"), prn);
	err = 1;
    }

    if (!err && fixsample) {
	err = check_for_overlap(dset, addset, &offset);
    } else if (!err && gretl_function_depth() > 0) {
	/* we won't add observations within a function, but
	   we should still check for an error from compare_ranges()
	*/
	if (!addsimple && !addpanel) {
	    addobs = compare_ranges(dset, addset, addvars, &offset,
				    &yrspecial, &err);
	    if (!err && addobs > 0) {
		addobs = 0;
	    }
	}
    } else if (!err) {
	if (!addsimple && !addpanel) {
	    addobs = compare_ranges(dset, addset, addvars, &offset,
				    &yrspecial, &err);
#if MERGE_DEBUG
	    fprintf(stderr, " added obs, from compare_ranges: %d (offset %d)\n",
		    addobs, offset);
#endif
	}
	if (!err && addobs <= 0 && addvars == 0) {
	    addobs = just_append_rows(dset, addset, &offset);
#if MERGE_DEBUG
	    fprintf(stderr, " added obs, from just_append_rows: %d\n", addobs);
#endif
	}
    }

    if (!err && (addobs < 0 || addvars < 0)) {
	merge_error(_("New data not conformable for appending\n"), prn);
	err = E_DATA;
    } else if (!err && addobs == 0 && addvars == 0) {
	merge_warning(_("No data were found for appending"), prn);
	return 0;
    }

    if (!err && !addpanel && dset->markers != addset->markers) {
	if (addobs == 0 && addvars == 0) {
	    if (update_overlap) {
		; /* might be OK? */
	    } else {
		gretl_errmsg_set(_("Found no data conformable for appending"));
		err = E_DATA;
	    }
	} else if (addset->n != dset->n && !yrspecial && !dayspecial) {
	    merge_error(_("Inconsistency in observation markers\n"), prn);
	    err = E_DATA;
	} else if (addset->markers && !dset->markers &&
		   !yrspecial && !dayspecial) {
	    dataset_destroy_obs_markers(addset);
	}
    }

#if MERGE_DEBUG
    if (!err) {
	fprintf(stderr, " after preliminaries: addvars = %d, addobs = %d\n",
		addvars, addobs);
    } else {
	fprintf(stderr, " after preliminaries: err = %d\n", err);
    }
#endif

    /* if checks are passed, try merging the data */

    if (!err && addobs > 0) {
	err = merge_lengthen_series(dset, addset, addobs, offset);
	if (err) {
	    merge_error(_("Out of memory!\n"), prn);
	}
    }

    if (!err && addpanel) {
	err = panel_append_special(addvars, dset, addset,
				   opt, prn);
    } else if (!err) {
	int k = dset->v;
	int i, t;

	if (addvars > 0 && dataset_add_series(dset, addvars)) {
	    merge_error(_("Out of memory!\n"), prn);
	    err = E_ALLOC;
	}

	for (i=1; i<addset->v && !err; i++) {
	    int v = series_index(dset, addset->varname[i]);
	    int tmin, newvar = v >= k;

	    if (!newvar && !update_overlap) {
		tmin = orig_n;
	    } else {
		tmin = 0;
	    }

	    if (newvar) {
		v = k++;
		strcpy(dset->varname[v], addset->varname[i]);
		copy_varinfo(dset->varinfo[v], addset->varinfo[i]);
		if (is_string_valued(addset, i) &&
		    addset->n == dset->n && offset == 0 &&
		    addobs == 0) {
		    /* attach the string table to the target
		       series and detach it from @addset
		    */
		    series_table *st;

		    st = series_get_string_table(addset, i);
		    series_attach_string_table(dset, v, st);
		    series_attach_string_table(addset, i, NULL);
		}
	    } else {
		/* not a new series */
		int lsval = is_string_valued(dset, v);
		int rsval = is_string_valued(addset, i);

		if (lsval + rsval == 1) {
		    gretl_errmsg_set(_("Can't concatenate string-valued and numeric series"));
		    err = E_DATA;
		} else if (lsval) {
		    err = merge_string_tables(dset, v, addset, i);
		}
	    }

	    if (dayspecial) {
		char obs[OBSLEN];
		int s;

		for (t=tmin; t<dset->n; t++) {
		    ntolabel(obs, t, dset);
		    s = dateton(obs, addset);
		    if (s >= 0 && s < addset->n) {
			dset->Z[v][t] = addset->Z[i][s];
		    } else {
			dset->Z[v][t] = NADBL;
		    }
		}
	    } else if (yrspecial) {
		int s;

		if (newvar) {
		    for (t=0; t<dset->n; t++) {
			dset->Z[v][t] = NADBL;
		    }
		}
		for (s=0; s<addset->n; s++) {
		    t = dateton(addset->S[s], dset);
		    if (t >= tmin && t < dset->n) {
			dset->Z[v][t] = addset->Z[i][s];
		    }
		}
	    } else {
		for (t=tmin; t<dset->n; t++) {
		    if (t >= offset && t - offset < addset->n) {
			dset->Z[v][t] = addset->Z[i][t - offset];
		    } else if (newvar) {
			dset->Z[v][t] = NADBL;
		    }
		}
	    }
	}
    }

    if (!err && (addvars || addobs) && gretl_messages_on()) {
	pputs(prn, _("Data appended OK\n"));
    }

    return err;
}

/* We want to ensure that calendar dates are recorded as per
   ISO 8601 -- that is, YYYY-MM-DD; here we remedy dates
   recorded in the form YYYY/MM/DD.
*/

static void maybe_fix_calendar_dates (DATASET *dset)
{
    if (strchr(dset->stobs, '/') != NULL) {
	gretl_charsub(dset->stobs, '/', '-');
	gretl_charsub(dset->endobs, '/', '-');
	if (dset->S != NULL && dset->markers == DAILY_DATE_STRINGS) {
	    int t;

	    for (t=0; t<dset->n; t++) {
		gretl_charsub(dset->S[t], '/', '-');
	    }
	}
    }
}

/**
 * get_merge_opts:
 * @opt: gretl options flags.
 *
 * Returns: just those components of @opt (if any) that
 * can be passed to merge_or_replace_data(); may be
 * useful when calling that function in the context
 * of a command only some of whose options should be
 * forwarded.
 */

gretlopt get_merge_opts (gretlopt opt)
{
    gretlopt merge_opt = OPT_NONE;

    if (opt & OPT_T) {
	/* panel, common time-series */
	merge_opt |= OPT_T;
    }
    if (opt & OPT_U) {
	/* update overlapping observations */
	merge_opt |= OPT_U;
    }
    if (opt & OPT_X) {
	/* fixed sample range */
	merge_opt |= OPT_X;
    }

    return merge_opt;
}

/* Apparatus for converting a dataset read from file
   into a gretl matrix, as opposed to replacing an
   existing dataset or merging with it.
*/

static gretl_matrix **dset_matrix;

void set_dset_matrix_target (gretl_matrix **pm)
{
    dset_matrix = pm;
}

void *get_dset_matrix_target (void)
{
    return dset_matrix;
}

/**
 * merge_or_replace_data:
 * @dset0: original dataset struct.
 * @pdset1: new dataset struct.
 * @opt: zero or more option flags (OPT_K presrves @pdset1,
 * otherwise it is destroyed).
 * @prn: print struct to accept messages.
 *
 * Given a newly-created dataset, pointed to by @pdset1, either
 * attempt to merge it with @dset0, if the original data array
 * is non-NULL, or replace the content of the original pointer
 * with the new dataset.
 *
 * In case merging is not successful, the new dataset is
 * destroyed.
 *
 * Returns: 0 on successful completion, non-zero otherwise.
 */

int merge_or_replace_data (DATASET *dset0, DATASET **pdset1,
			   gretlopt opt, PRN *prn)
{
    int keep = (opt & OPT_K);
    int err = 0;

    if (dset_matrix != NULL) {
	/* Convert the new dataset to matrix; don't touch
	   the existing dataset, if any.
	*/
	const DATASET *dset = *pdset1;

	*dset_matrix = gretl_matrix_data_subset(NULL, dset,
						0, dset->n - 1,
						M_MISSING_OK,
						&err);
	destroy_dataset(*pdset1);
	*pdset1 = NULL;
	return err;
    }

    if (dset0->Z != NULL) {
	/* we have an existing dataset into which the new data
	   should be merged */
	gretlopt merge_opt = OPT_NONE;

	if (opt & OPT_T) {
	    /* panel, common time-series */
	    merge_opt |= OPT_T;
	}
	if (opt & OPT_U) {
	    /* update overlapping observations */
	    merge_opt |= OPT_U;
	}
	if (opt & OPT_X) {
	    /* fixed sample range */
	    merge_opt |= OPT_X;
	}
	err = merge_data(dset0, *pdset1, merge_opt, prn);
	if (!keep) {
	    destroy_dataset(*pdset1);
	}
    } else {
	/* starting from scratch */
	*dset0 = **pdset1;
	free(*pdset1);
	if (calendar_data(dset0)) {
	    maybe_fix_calendar_dates(dset0);
	}
    }

    if (!keep) {
	*pdset1 = NULL;
    }

    set_current_dataset(dset0);

    return err;
}

static int check_imported_string (char *src, int i, size_t len)
{
    int err = 0;

    if (!g_utf8_validate(src, -1, NULL)) {
	gchar *trstr = NULL;
	gsize bytes;

	trstr = g_locale_to_utf8(src, -1, NULL, &bytes, NULL);

	if (trstr == NULL) {
	    gretl_errmsg_sprintf(_("Invalid characters in imported string, line %d"), i);
	    err = E_DATA;
	} else {
	    *src = '\0';
	    strncat(src, trstr, len - 1);
	    g_free(trstr);
	}
    }

    return err;
}

static int count_markers (FILE *fp, char *line, int linelen,
			  char *marker)
{
    int n = 0;

    while (fgets(line, linelen, fp)) {
	if (sscanf(line, "%31[^\n\r]", marker) == 1) {
	    g_strstrip(marker);
	    if (*marker != '\0') {
		n++;
	    }
	}
    }

    rewind(fp);

    return n;
}

static void finalize_add_markers (DATASET *dset, char **S, int err)
{
    if (err) {
	strings_array_free(S, dset->n);
    } else {
	if (dset->S != NULL) {
	    strings_array_free(dset->S, dset->n);
	}
	dset->markers = REGULAR_MARKERS;
	dset->S = S;
    }
}

/**
 * add_obs_markers_from_file:
 * @dset: data information struct.
 * @fname: name of file containing case markers.
 *
 * Read case markers (strings of %OBSLEN - 1 characters or less that identify
 * the observations) from a file, and associate them with the
 * current data set.  The file should contain one marker per line,
 * with a number of lines equal to the number of observations in
 * the current data set.
 *
 * Returns: 0 on successful completion, non-zero otherwise.
 */

int add_obs_markers_from_file (DATASET *dset, const char *fname)
{
    char **S = NULL;
    FILE *fp;
    char line[128], marker[32];
    int ns, panel_units = 0;
    int t, err = 0;

    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	return E_FOPEN;
    }

    ns = count_markers(fp, line, sizeof line, marker);
    if (dataset_is_panel(dset) && ns == dset->n / dset->pd) {
	panel_units = 1; /* OK */
    } else if (ns < dset->n) {
	return E_INVARG;
    }

    S = strings_array_new_with_length(dset->n, OBSLEN);
    if (S == NULL) {
	fclose(fp);
	return E_ALLOC;
    }

    if (panel_units) {
	int T = dset->pd;
	int t, i = 0;

	while (fgets(line, sizeof line, fp) && !err) {
	    *marker = '\0';
	    if (sscanf(line, "%31[^\n\r]", marker) == 1) {
		g_strstrip(marker);
		strncat(S[i], marker, OBSLEN - 1);
		err = check_imported_string(S[i], i+1, OBSLEN);
		if (!err) {
		    /* copy to remaining observations */
		    for (t=1; t<T; t++) {
			strcpy(S[i+t], S[i]);
		    }
		}
		i += T;
	    }
	}
    } else {
	for (t=0; t<dset->n && !err; t++) {
	    if (fgets(line, sizeof line, fp) == NULL) {
		gretl_errmsg_sprintf(_("Expected %d markers; found %d\n"),
				     dset->n, t);
		err = E_DATA;
	    } else if (sscanf(line, "%31[^\n\r]", marker) != 1) {
		gretl_errmsg_sprintf(_("Couldn't read marker on line %d"), t+1);
		err = E_DATA;
	    } else {
		g_strstrip(marker);
		strncat(S[t], marker, OBSLEN - 1);
		err = check_imported_string(S[t], t+1, OBSLEN);
	    }
	}
    }

    finalize_add_markers(dset, S, err);

    return err;
}

static int add_obs_markers_from_array (DATASET *dset, const char *aname)
{
    gretl_array *a;
    char **S = NULL;
    const char *s;
    const char *white = " \t\r\n";
    char marker[32];
    int ns, panel_units = 0;
    int t, err = 0;

    a = get_strings_array_by_name(aname);
    if (a == NULL) {
	return E_INVARG;
    }

    ns = gretl_array_get_length(a);
    if (dataset_is_panel(dset) && ns == dset->n / dset->pd) {
	panel_units = 1; /* OK */
    } else if (ns < dset->n) {
	return E_INVARG;
    }

    S = strings_array_new_with_length(dset->n, OBSLEN);
    if (S == NULL) {
	return E_ALLOC;
    }

    if (panel_units) {
	int N = dset->n / dset->pd;
	int T = dset->pd;
	int i, j = 0;

	for (i=0; i<N; i++) {
	    s = gretl_array_get_data(a, i);
	    s += strspn(s, white);
	    *marker = '\0';
	    strncat(marker, s, 31);
	    g_strstrip(marker);
	    gretl_utf8_truncate_b(marker, OBSLEN-1);
	    strcpy(S[j++], marker);
	    for (t=1; t<T; t++) {
		strcpy(S[j++], marker);
	    }
	}
    } else {
	for (t=0; t<dset->n; t++) {
	    s = gretl_array_get_data(a, t);
	    s += strspn(s, white);
	    *marker = '\0';
	    strncat(marker, s, 31);
	    g_strstrip(marker);
	    gretl_utf8_truncate_b(marker, OBSLEN-1);
	    strcpy(S[t], marker);
	}
    }

    finalize_add_markers(dset, S, err);

    return err;
}

static int add_obs_markers_from_series (DATASET *dset, const char *vname)
{
    char **S = NULL;
    const char *s;
    const char *white = " \t\r\n";
    char marker[32];
    int v, t, err = 0;

    v = current_series_index(dset, vname);
    if (v <= 0) {
	return E_UNKVAR;
    } else if (!is_string_valued(dset, v)) {
	return E_INVARG;
    }

    S = strings_array_new_with_length(dset->n, OBSLEN);
    if (S == NULL) {
	return E_ALLOC;
    }

    for (t=0; t<dset->n; t++) {
	s = series_get_string_for_obs(dset, v, t);
	if (s == NULL) {
	    gretl_errmsg_sprintf(_("Missing marker at observation %d"), t+1);
	    err = E_DATA;
	    break;
	} else {
	    s += strspn(s, white);
	    *marker = '\0';
	    strncat(marker, s, 31);
	    g_strstrip(marker);
	    gretl_utf8_truncate_b(marker, OBSLEN-1);
	    strcpy(S[t], marker);
	}
    }

    finalize_add_markers(dset, S, err);

    return err;
}

/**
 * dataset_has_var_labels:
 * @dset: data information struct.
 *
 * Returns: 1 if at least one variable in the current dataset
 * has a descriptive label, otherwise 0.
 */

int dataset_has_var_labels (const DATASET *dset)
{
    const char *vlabel;
    int i, imin = 1;

    if (dset->v > 1) {
	if (!strcmp(dset->varname[1], "index")) {
	    vlabel = series_get_label(dset, 1);
	    if (vlabel != NULL && !strcmp(vlabel, _("index variable"))) {
		imin = 2;
	    }
	}
    }

    for (i=imin; i<dset->v; i++) {
	vlabel = series_get_label(dset, i);
	if (vlabel != NULL && *vlabel != '\0') {
	    return 1;
	}
    }

    return 0;
}

/**
 * save_var_labels_to_file:
 * @dset: data information struct.
 * @fname: name of file containing labels.
 *
 * Writes to @fname the descriptive labels for the series in
 * the current dataset.
 *
 * Returns: 0 on successful completion, non-zero otherwise.
 */

int save_var_labels_to_file (const DATASET *dset,
			     const char *fname)
{
    const char *vlabel;
    FILE *fp;
    int i, err = 0;

    fp = gretl_fopen(fname, "w");

    if (fp == NULL) {
	err = E_FOPEN;
    } else {
	for (i=1; i<dset->v; i++) {
	    vlabel = series_get_label(dset, i);
	    fprintf(fp, "%s\n", vlabel == NULL ? "" : vlabel);
	}
	fclose(fp);
    }

    return err;
}

static int save_dset_description_to_file (const DATASET *dset,
					  const char *fname)
{
    FILE *fp;
    int err = 0;

    fp = gretl_fopen(fname, "w");

    if (fp == NULL) {
	err = E_FOPEN;
    } else {
	fputs(dset->descrip, fp);
	fclose(fp);
    }

    return err;
}

static int save_var_labels_to_array (const DATASET *dset,
				     const char *aname)
{
    gretl_array *a = NULL;
    int err = 0;

    if (gretl_is_series(aname, dset)) {
	err = E_TYPES;
    } else {
	err = check_identifier(aname);
    }

    if (!err) {
	a = gretl_array_new(GRETL_TYPE_STRINGS, dset->v - 1, &err);
    }

    if (!err) {
	err = user_var_add_or_replace(aname, GRETL_TYPE_STRINGS, a);
    }

    if (!err) {
	char *vlabel;
	int i;

	for (i=1; i<dset->v; i++) {
	    vlabel = (char *) series_get_label(dset, i);
	    gretl_array_set_element(a, i-1, vlabel != NULL ? vlabel : "",
				    GRETL_TYPE_STRING, 1);
	}
    }

    if (err && a != NULL) {
	gretl_array_destroy(a);
	a = NULL;
    }

    return err;
}

static int save_obs_markers_to_array (const DATASET *dset,
				      const char *aname)
{
    gretl_array *a = NULL;
    int err = 0;

    if (gretl_is_series(aname, dset)) {
	err = E_TYPES;
    } else {
	err = check_identifier(aname);
    }

    if (!err) {
	a = gretl_array_new(GRETL_TYPE_STRINGS, dset->n, &err);
    }

    if (!err) {
	err = user_var_add_or_replace(aname, GRETL_TYPE_STRINGS, a);
    }

    if (!err) {
	char *marker;
	int i;

	for (i=0; i<dset->n; i++) {
	    marker = dset->S[i];
	    gretl_array_set_element(a, i, marker != NULL ? marker : "",
				    GRETL_TYPE_STRING, 1);
	}
    }

    if (err && a != NULL) {
	gretl_array_destroy(a);
	a = NULL;
    }

    return err;
}

/**
 * add_var_labels_from_file:
 * @dset: data information struct.
 * @fname: name of file containing labels.
 *
 * Read descriptive variables for labels (strings of %MAXLABEL - 1
 * characters or less) from a file, and associate them with the
 * current data set.  The file should contain one label per line,
 * with a number of lines equal to the number of variables in
 * the current data set, excluding the constant.
 *
 * Returns: 0 on successful completion, non-zero otherwise.
 */

int add_var_labels_from_file (DATASET *dset, const char *fname)
{
    FILE *fp;
    char line[1024];
    gchar *label;
    int nlabels = 0;
    int i, err = 0;

    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	return E_FOPEN;
    }

    for (i=1; i<dset->v && !err; i++) {
	if (fgets(line, sizeof line, fp) == NULL) {
	    break;
	} else {
	    label = g_strstrip(g_strdup(line));
	    if (strlen(label) > 0) {
		if (!g_utf8_validate(label, -1, NULL)) {
		    gchar *trstr = NULL;
		    gsize bytes;

		    trstr = g_locale_to_utf8(label, -1, NULL,
					     &bytes, NULL);
		    if (trstr != NULL) {
			series_set_label(dset, i, trstr);
			nlabels++;
			g_free(trstr);
		    }
		} else {
		    series_set_label(dset, i, label);
		    nlabels++;
		}
	    }
	    g_free(label);
	}
    }

    if (!err && nlabels == 0) {
	gretl_errmsg_set(_("No labels found"));
	err = E_DATA;
    }

    return err;
}

static int add_dset_description_from_file (DATASET *dset,
					   const char *fname)
{
    gchar *buf = NULL;
    GError *gerr = NULL;
    gsize size;
    int err = 0;

    g_file_get_contents(fname, &buf, &size, &gerr);

    if (gerr != NULL) {
	gretl_errmsg_set(gerr->message);
	g_error_free(gerr);
	return E_DATA;
    }

    if (size == 0) {
	gretl_errmsg_set("No content found");
	err = E_DATA;
    } else if (g_utf8_validate(buf, -1, NULL)) {
	free(dset->descrip);
	dset->descrip = gretl_strdup(buf);
    } else {
	gchar *trstr = NULL;

	trstr = g_locale_to_utf8(buf, -1, NULL,
				 &size, &gerr);
	if (gerr != NULL) {
	    gretl_errmsg_set(gerr->message);
	    g_error_free(gerr);
	    err = E_DATA;
	} else {
	    free(dset->descrip);
	    dset->descrip = gretl_strdup(trstr);
	    g_free(trstr);
	}
    }

    g_free(buf);

    return err;
}

static int add_var_labels_from_array (DATASET *dset, const char *aname)
{
    gretl_array *a = get_array_by_name(aname);
    int i, err = 0;

    if (a == NULL) {
	gretl_errmsg_sprintf(_("%s: no such array"), aname);
	err = E_DATA;
    } else if (gretl_array_get_type(a) != GRETL_TYPE_STRINGS ||
	       gretl_array_get_length(a) < dset->v - 1) {
	err = E_TYPES;
    }

    for (i=1; i<dset->v && !err; i++) {
	const char *s = gretl_array_get_data(a, i-1);

	series_set_label(dset, i, s);
    }

    return err;
}

int read_or_write_var_labels (gretlopt opt, DATASET *dset, PRN *prn)
{
    const char *lname = NULL;
    int err;

    err = incompatible_options(opt, OPT_D | OPT_T | OPT_F |
			       OPT_A | OPT_R);
    if (err) {
	return err;
    }

    if (opt & (OPT_T | OPT_F | OPT_A | OPT_R)) {
	lname = get_optval_string(LABELS, opt);
	if (lname == NULL) {
	    return E_BADOPT;
	} else if (opt & (OPT_T | OPT_F)) {
	    gretl_maybe_switch_dir(lname);
	}
    }

    if (opt & OPT_D) {
	/* delete */
	int i;

	for (i=1; i<dset->v; i++) {
	    series_set_label(dset, i, "");
	}
    } else if (opt & (OPT_T | OPT_R)) {
	/* to-file, to-array */
	if (!dataset_has_var_labels(dset)) {
	    pprintf(prn, _("No labels are available for writing\n"));
	    err = E_DATA;
	} else {
	    if (opt & OPT_T) {
		err = save_var_labels_to_file(dset, lname);
	    } else {
		err = save_var_labels_to_array(dset, lname);
	    }
	    if (!err && gretl_messages_on()) {
		pprintf(prn, _("Labels written OK\n"));
	    }
	}
    } else if (opt & (OPT_F | OPT_A)) {
	/* from-file, from-array */
	if (opt & OPT_F) {
	    err = add_var_labels_from_file(dset, lname);
	} else {
	    err = add_var_labels_from_array(dset, lname);
	}
	if (!err && gretl_messages_on()) {
	    pprintf(prn, _("Labels loaded OK\n"));
	}
    }

    return err;
}

int read_or_write_dset_description (gretlopt opt, DATASET *dset, PRN *prn)
{
    const char *fname = NULL;
    int err;

    if (dset == NULL) {
	return E_NODATA;
    }

    err = incompatible_options(opt, OPT_F | OPT_T);
    if (err) {
	return err;
    }

    fname = get_optval_string(INFO, opt);
    if (fname == NULL) {
	err = E_BADOPT;
    }

    gretl_maybe_switch_dir(fname);

    if (opt & OPT_T) {
	/* to-file */
	if (dset->descrip == NULL) {
	    pputs(prn, _("No description is available\n"));
	    err = E_DATA;
	} else {
	    err = save_dset_description_to_file(dset, fname);
	    if (!err && gretl_messages_on()) {
		pputs(prn, _("dataset description written OK\n"));
	    }
	}
    } else if (opt & OPT_F) {
	/* from-file */
	int have_desc = dset->descrip != NULL;

	err = add_dset_description_from_file(dset, fname);
	if (!err && gretl_messages_on()) {
	    if (have_desc) {
		pputs(prn, _("description replaced OK\n"));
	    } else {
		pputs(prn, _("description loaded OK\n"));
	    }
	}
    }

    return err;
}

static int save_obs_markers_to_file (DATASET *dset, const char *fname)
{
    FILE *fp = gretl_fopen(fname, "w");
    int err = 0;

    if (fp == NULL) {
	err = E_FOPEN;
    } else {
	int i;

	for (i=0; i<dset->n; i++) {
	    fprintf(fp, "%s\n", dset->S[i]);
	}
	fclose(fp);
    }

    return err;
}

int read_or_write_obs_markers (gretlopt opt, DATASET *dset, PRN *prn)
{
    const char *fname = NULL;
    int err;

    err = incompatible_options(opt, OPT_D | OPT_T | OPT_F |
			       OPT_A | OPT_R);
    if (err) {
	return err;
    }

    if (opt & (OPT_T | OPT_F)) {
	fname = get_optval_string(MARKERS, opt);
	if (fname == NULL) {
	    return E_BADOPT;
	} else {
	    fname = gretl_maybe_switch_dir(fname);
	}
    }

    if (opt & (OPT_A | OPT_T)) {
	/* writing to file or array */
	if (dset->S == NULL) {
	    gretl_errmsg_set(_("No markers are available for writing"));
	    return E_DATA;
	}
    }

    if (opt & OPT_D) {
	/* delete */
	dataset_destroy_obs_markers(dset);
    } else if (opt & OPT_T) {
	/* to-file */
	err = save_obs_markers_to_file(dset, fname);
	if (!err && gretl_messages_on()) {
	    pprintf(prn, _("Markers written OK\n"));
	}
    } else if (opt & OPT_F) {
	/* from-file */
	err = add_obs_markers_from_file(dset, fname);
	if (!err && gretl_messages_on()) {
	    pprintf(prn, _("Markers loaded OK\n"));
	}
    } else if (opt & OPT_A) {
	/* to-array */
	const char *aname = get_optval_string(MARKERS, OPT_A);

	err = save_obs_markers_to_array(dset, aname);
    } else if (opt & OPT_R) {
	/* from-array */
	const char *aname = get_optval_string(MARKERS, OPT_R);

	err = add_obs_markers_from_array(dset, aname);
    } else if (opt & OPT_S) {
	/* from-series */
	const char *vname = get_optval_string(MARKERS, OPT_S);

	err = add_obs_markers_from_series(dset, vname);
    }

    return err;
}

static void
octave_varname (char *name, const char *s, int nnum, int v)
{
    char nstr[12];
    int len, tr;

    if (nnum == 0) {
	strcpy(name, s);
    } else {
	sprintf(nstr, "%d", nnum);
	len = strlen(nstr);
	tr = VNAMELEN - len;

	if (tr > 0) {
	    strncat(name, s, tr);
	    strcat(name, nstr);
	} else {
	    sprintf(name, "v%d", v);
	}
    }
}

static int get_max_line_length (FILE *fp, PRN *prn)
{
    int c, c1, cc = 0;
    int maxlen = 0;

    while ((c = fgetc(fp)) != EOF) {
	if (c == 0x0d) {
	    /* CR */
	    c1 = fgetc(fp);
	    if (c1 == EOF) {
		break;
	    } else if (c1 == 0x0a) {
		/* CR + LF -> LF */
		c = c1;
	    } else {
		/* Mac-style: CR not followed by LF */
		c = 0x0a;
		ungetc(c1, fp);
	    }
	}
	if (c == 0x0a) {
	    if (cc > maxlen) {
		maxlen = cc;
	    }
	    cc = 0;
	    continue;
	}
	if (!isspace((unsigned char) c) && !isprint((unsigned char) c) &&
	    !(c == CTRLZ)) {
	    pprintf(prn, _("Binary data (%d) encountered: this is not a valid "
			   "text file\n"), c);
	    return -1;
	}
	cc++;
    }

    if (maxlen == 0) {
	pprintf(prn, _("Data file is empty\n"));
    }

    if (maxlen > 0) {
	/* allow for newline and null terminator */
	maxlen += 3;
    }

    return maxlen;
}

static int import_octave (const char *fname, DATASET *dset,
			  gretlopt opt, PRN *prn)
{
    DATASET *octset = NULL;
    FILE *fp = NULL;
    char *line = NULL;
    char tmp[8], fmt[16], name[32];
    int nrows = 0, ncols = 0, nblocks = 0;
    int brows = 0, bcols = 0, oldbcols = 0;
    int maxlen, got_type = 0, got_name = 0;
    int i, t, err = 0;

    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	return E_FOPEN;
    }

    pprintf(prn, "%s %s...\n", _("parsing"), fname);

    maxlen = get_max_line_length(fp, prn);
    if (maxlen <= 0) {
	err = E_DATA;
	goto oct_bailout;
    }

    line = malloc(maxlen);
    if (line == NULL) {
	err = E_ALLOC;
	goto oct_bailout;
    }

    pprintf(prn, _("   longest line: %d characters\n"), maxlen - 1);

    rewind(fp);

    while (fgets(line, maxlen, fp) && !err) {
	if (*line == '#') {
	    if (!got_name) {
		if (sscanf(line, "# name: %31s", name) == 1) {
		    got_name = 1;
		    nblocks++;
		    continue;
		}
	    }
	    if (!got_type) {
		if (sscanf(line, "# type: %7s", tmp) == 1) {
		    if (!got_name || strcmp(tmp, "matrix")) {
			err = 1;
		    } else {
			got_type = 1;
		    }
		    continue;
		}
	    }
	    if (brows == 0) {
		if (sscanf(line, "# rows: %d", &brows) == 1) {
		    if (!got_name || !got_type || brows <= 0) {
			err = 1;
		    } else if (nrows > 0 && brows != nrows) {
			err = 1;
		    } else {
			nrows = brows;
		    }
		    continue;
		}
	    }
	    if (bcols == 0) {
		if (sscanf(line, "# columns: %d", &bcols) == 1) {
		    if (!got_name || !got_type || bcols <= 0) {
			err = 1;
		    } else {
			ncols += bcols;
			pprintf(prn, _("   Found matrix '%s' with "
				       "%d rows, %d columns\n"), name, brows, bcols);
		    }
		    continue;
		}
	    }
	} else if (string_is_blank(line)) {
	    continue;
	} else {
	    got_name = 0;
	    got_type = 0;
	    brows = 0;
	    bcols = 0;
	}
    }

    if (err || nrows == 0 || ncols == 0) {
	pputs(prn, _("Invalid data file\n"));
	err = E_DATA;
	goto oct_bailout;
    }

    /* initialize datainfo and Z */

    octset = datainfo_new();
    if (octset == NULL) {
	pputs(prn, _("Out of memory!\n"));
	err = E_ALLOC;
	goto oct_bailout;
    }

    octset->n = nrows;
    octset->v = ncols + 1;

    if (start_new_Z(octset, 0)) {
	pputs(prn, _("Out of memory!\n"));
	err = E_ALLOC;
	goto oct_bailout;
    }

    rewind(fp);

    pprintf(prn, _("   number of variables: %d\n"), ncols);
    pprintf(prn, _("   number of observations: %d\n"), nrows);
    pprintf(prn, _("   number of data blocks: %d\n"), nblocks);

    i = 1;
    t = 0;

    sprintf(fmt, "# name: %%%ds", VNAMELEN - 1);

    while (fgets(line, maxlen, fp) && !err) {
	char *s = line;
	int j;

	if (*s == '#') {
	    if (sscanf(line, fmt, name) == 1) {
		;
	    } else if (sscanf(line, "# rows: %d", &brows) == 1) {
		t = 0;
	    } else if (sscanf(line, "# columns: %d", &bcols) == 1) {
		i += oldbcols;
		oldbcols = bcols;
	    }
	}

	if (*s == '#' || string_is_blank(s)) {
	    continue;
	}

	if (t >= octset->n) {
	    err = 1;
	}

	for (j=0; j<bcols && !err; j++) {
	    double x;
	    int v = i + j;

	    if (t == 0) {
		int nnum = (bcols > 1)? j + 1 : 0;

		octave_varname(octset->varname[i+j], name, nnum, v);
	    }

	    while (isspace(*s)) s++;
	    if (sscanf(s, "%lf", &x) != 1) {
		fprintf(stderr, "error: '%s', didn't get double\n", s);
		err = 1;
	    } else {
		octset->Z[v][t] = x;
		while (!isspace(*s)) s++;
	    }
	}
	t++;
    }

    if (err) {
	pputs(prn, _("Invalid data file\n"));
	err = E_DATA;
    } else {
	int merge = dset->Z != NULL;
	gretlopt merge_opt = 0;

	if (merge && (opt & OPT_T)) {
	    merge_opt = OPT_T;
	}
	err = merge_or_replace_data(dset, &octset, merge_opt, prn);
    }

 oct_bailout:

    if (fp != NULL) {
	fclose(fp);
    }

    if (line != NULL) {
	free(line);
    }

    if (octset != NULL) {
	clear_datainfo(octset, CLEAR_FULL);
    }

    return err;
}

/**
 * import_other:
 * @fname: name of file.
 * @ftype: type of data file.
 * @dset: pointer to dataset struct.
 * @opt: option flag; see gretl_get_data().
 * @prn: gretl printing struct.
 *
 * Open a data file of a type that requires a special plugin.
 *
 * Returns: 0 on successful completion, non-zero otherwise.
 */

int import_other (const char *fname, GretlFileType ftype,
		  DATASET *dset, gretlopt opt, PRN *prn)
{
    FILE *fp;
    int (*importer) (const char *, DATASET *,
		     gretlopt, PRN *);
    int err = 0;

    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	pprintf(prn, _("Couldn't open %s\n"), fname);
	err = E_FOPEN;
	goto bailout;
    }

    fclose(fp);

    if (ftype == GRETL_OCTAVE) {
	/* plugin not needed */
	return import_octave(fname, dset, opt, prn);
    }

    if (ftype == GRETL_WF1) {
	importer = get_plugin_function("wf1_get_data");
    } else if (ftype == GRETL_DTA) {
	importer = get_plugin_function("dta_get_data");
    } else if (ftype == GRETL_SAV) {
	importer = get_plugin_function("sav_get_data");
    } else if (ftype == GRETL_SAS) {
	importer = get_plugin_function("xport_get_data");
    } else if (ftype == GRETL_JMULTI) {
	importer = get_plugin_function("jmulti_get_data");
    } else if (ftype == GRETL_MAP) {
	importer = get_plugin_function("map_get_data");
    } else {
	pprintf(prn, _("Unrecognized data type"));
	pputc(prn, '\n');
	return E_DATA;
    }

    if (importer == NULL) {
        err = 1;
    } else {
	err = (*importer)(fname, dset, opt, prn);
    }

 bailout:

    return err;
}

/**
 * import_spreadsheet:
 * @fname: name of file.
 * @ftype: type of data file.
 * @list: list of parameters for spreadsheet import, or NULL.
 * @sheetname: name of specific worksheet, or NULL.
 * @dset: dataset struct.
 * @opt: option flag; see gretl_get_data().
 * @prn: gretl printing struct.
 *
 * Open a data file of a type that requires a special plugin.
 * Acceptable values for @ftype are %GRETL_GNUMERIC,
 * %GRETL_XLS, %GRETL_XLSX and %GRETL_ODS.
 *
 * Returns: 0 on successful completion, non-zero otherwise.
 */

int import_spreadsheet (const char *fname, GretlFileType ftype,
			int *list, char *sheetname,
			DATASET *dset, gretlopt opt, PRN *prn)
{
    FILE *fp;
    int (*importer) (const char*, int *, char *,
		     DATASET *, gretlopt, PRN *,
		     void *);
    int err = 0;

    import_na_init();

    fp = gretl_fopen(fname, "r");

    if (fp == NULL) {
	pprintf(prn, _("Couldn't open %s\n"), fname);
	err = E_FOPEN;
	goto bailout;
    }

    fclose(fp);

    if (ftype == GRETL_GNUMERIC) {
	importer = get_plugin_function("gnumeric_get_data");
    } else if (ftype == GRETL_XLS) {
	importer = get_plugin_function("xls_get_data");
    } else if (ftype == GRETL_XLSX) {
	importer = get_plugin_function("xlsx_get_data");
    } else if (ftype == GRETL_ODS) {
	importer = get_plugin_function("ods_get_data");
    } else {
	pprintf(prn, _("Unrecognized data type"));
	pputc(prn, '\n');
	return E_DATA;
    }

    if (importer == NULL) {
        err = 1;
    } else {
	gchar *thisdir = g_get_current_dir();

	err = (*importer)(fname, list, sheetname, dset, opt, prn, NULL);

	if (thisdir != NULL) {
	    /* come back out of dotdir? */
	    gretl_chdir(thisdir);
	    g_free(thisdir);
	}
    }

 bailout:

    return err;
}

static int is_jmulti_datafile (const char *fname)
{
    FILE *fp;
    int ret = 0;

    fp = gretl_fopen(fname, "r");

    if (fp != NULL) {
	char test[128] = {0};
	int gotobs = 0;
	int gotcomm = 0;
	int incomm = 0;

	/* look for characteristic C-style comment and
	   <obs stuff> field, outside of comment */

	while (fgets(test, sizeof test, fp)) {
	    if (!incomm && strstr(test, "/*")) {
		gotcomm = 1;
		incomm = 1;
	    }
	    if (incomm && strstr(test, "*/")) {
		incomm = 0;
	    }
	    if (!incomm && *test == '<' && strchr(test, '>')) {
		gotobs = 1;
	    }
	    if (gotcomm && gotobs) {
		ret = 1;
		break;
	    }
	}
	fclose(fp);
    }

    return ret;
}

/**
 * gretl_is_pkzip_file:
 * @fname: name of file to examine.
 *
 * Returns: 1 if @fname is readable and is a PKZIP file,
 * else 0.
 */

int gretl_is_pkzip_file (const char *fname)
{
    FILE *fp = gretl_fopen(fname, "rb");
    int ret = 0;

    if (fp != NULL) {
	char test[3] = {0};

	if (fread(test, 1, 2, fp) == 2 && !strcmp(test, "PK")) {
	    ret = 1;
	}
	fclose(fp);
    }

    return ret;
}

/**
 * detect_filetype:
 * @fname: the name of the file to test.
 * @opt: OPT_P may be included to permit path-searching if @fname
 * is not an absolute path; in that case the @fname argument
 * may be modified, otherwise it will be left unchanged.
 *
 * Attempts to determine the type of a file to be opened in gretl:
 * data file (of various formats), or command script. If OPT_P
 * is given, the @fname argument must be an array of length
 * at least %MAXLEN: a path may be prepended and in some cases
 * an extension may be appended.
 *
 * Returns: integer code indicating the type of file.
 */

GretlFileType detect_filetype (char *fname, gretlopt opt)
{
    const char *ext = get_filename_extension(fname);
    GretlFileType ftype = GRETL_UNRECOGNIZED;

    if (ext != NULL) {
	/* First try judging the type by extension */
	if (!strcmp(ext, ".inp")) {
	    ftype = GRETL_SCRIPT;
	} else if (!strcmp(ext, ".gretl")) {
	    if (gretl_is_pkzip_file(fname)) {
		ftype = GRETL_SESSION;
	    } else {
		ftype = GRETL_SCRIPT;
	    }
	} else {
	    ftype = data_file_type_from_extension(ext, fname);
	    if (ftype == GRETL_UNRECOGNIZED) {
		/* check for database types */
		if (!strcmp(ext, ".bin")) {
		    ftype = GRETL_NATIVE_DB;
		} else if (!strcmp(ext, ".rat")) {
		    ftype = GRETL_RATS_DB;
		} else if (!strcmp(ext, ".bn7")) {
		    ftype = GRETL_PCGIVE_DB;
		}
	    }
	}
	if (ftype != GRETL_UNRECOGNIZED) {
	    /* We got a type from the extension, but can we find
	       the file "as is"? If so, we're done.
	    */
	    if (gretl_test_fopen(fname, "r") == 0) {
		return ftype;
	    }
	}
    }

    if ((opt & OPT_P) && gretl_addpath(fname, 0) != NULL) {
	ext = get_filename_extension(fname);
	if (ext != NULL) {
	    /* check again for known data file types */
	    ftype = data_file_type_from_extension(ext, fname);
	}
    }

    if (ftype == GRETL_UNRECOGNIZED) {
	/* last gasp */
	if (gretl_is_xml_file(fname)) {
	    ftype = GRETL_XML_DATA;
	} else if (has_suffix(fname, ".dat") && is_jmulti_datafile(fname)) {
	    ftype = GRETL_JMULTI;
	} else {
	    /* default to assuming plain text data */
	    ftype = GRETL_CSV;
	}
    }

    return ftype;
}

/**
 * check_atof:
 * @numstr: string to check.
 *
 * Returns: 0 if @numstr is blank, or is a valid string representation
 * of a floating point number, else 1.
 */

int check_atof (const char *numstr)
{
    char *test;

    /* accept blank entries */
    if (*numstr == '\0') return 0;

    errno = 0;

    strtod(numstr, &test);

    if (*test == '\0' && errno != ERANGE) return 0;

    if (!strcmp(numstr, test)) {
	gretl_errmsg_sprintf(_("'%s' -- no numeric conversion performed!"), numstr);
	return 1;
    }

    if (*test != '\0') {
	if (isprint(*test)) {
	    gretl_errmsg_sprintf(_("Extraneous character '%c' in data"), *test);
	} else {
	    gretl_errmsg_sprintf(_("Extraneous character (0x%x) in data"), *test);
	}
	return 1;
    }

    if (errno == ERANGE) {
	gretl_errmsg_sprintf(_("'%s' -- number out of range!"), numstr);
    }

    return 1;
}

/**
 * check_atoi:
 * @numstr: string to check.
 *
 * Returns: 0 if @numstr is blank, or is a valid string representation
 * of an int, else 1.
 */

int check_atoi (const char *numstr)
{
    long int val;
    char *test;

    /* accept blank entries */
    if (*numstr == '\0') return 0;

    errno = 0;

    val = strtol(numstr, &test, 10);

    if (*test == '\0' && errno != ERANGE) return 0;

    if (!strcmp(numstr, test)) {
	gretl_errmsg_sprintf(_("'%s' -- no numeric conversion performed!"), numstr);
	return 1;
    }

    if (*test != '\0') {
	if (isprint(*test)) {
	    gretl_errmsg_sprintf(_("Extraneous character '%c' in data"), *test);
	} else {
	    gretl_errmsg_sprintf(_("Extraneous character (0x%x) in data"), *test);
	}
	return 1;
    }

    if (errno == ERANGE || val <= INT_MIN || val >= INT_MAX) {
	gretl_errmsg_sprintf(_("'%s' -- number out of range!"), numstr);
    }

    return 1;
}

static int transpose_varname_used (const char *vname,
				   DATASET *dinfo,
				   int imax)
{
    int i;

    for (i=0; i<imax; i++) {
	if (!strcmp(vname, dinfo->varname[i])) {
	    return 1;
	}
    }

    return 0;
}

/**
 * transpose_data:
 * @dset: pointer to dataset information struct.
 *
 * Attempts to transpose the current dataset, so that each
 * variable becomes interpreted as an observation and each
 * observation as a variable.
 *
 * Returns: 0 on success, non-zero error code on error.
 */

int transpose_data (DATASET *dset)
{
    DATASET *tset;
    int k = dset->n + 1;
    int T = dset->v - 1;
    int i, t;

    tset = create_new_dataset(k, T, 0);
    if (tset == NULL) {
	return E_ALLOC;
    }

    for (i=1; i<dset->v; i++) {
	for (t=0; t<dset->n; t++) {
	    tset->Z[t+1][i-1] = dset->Z[i][t];
	}
    }

    for (t=0; t<dset->n; t++) {
	int k = t + 1;
	char *targ = tset->varname[k];

	if (dset->S != NULL && dset->S[t][0] != '\0') {
	    int err;

	    *targ = '\0';
	    strncat(targ, dset->S[t], VNAMELEN - 1);
	    gretl_charsub(targ, ' ', '_');
	    err = check_varname(targ);
	    if (err) {
		sprintf(targ, "v%d", k);
		gretl_error_clear();
	    } else if (transpose_varname_used(targ, tset, k)) {
		sprintf(targ, "v%d", k);
	    }
	} else {
	    sprintf(targ, "v%d", k);
	}
    }

    free_Z(dset);
    dset->Z = tset->Z;

    clear_datainfo(dset, CLEAR_FULL);

    dset->v = k;
    dset->n = T;
    dset->t1 = 0;
    dset->t2 = dset->n - 1;

    dset->varname = tset->varname;
    dset->varinfo = tset->varinfo;

    dataset_obs_info_default(dset);

    free(tset);

    return 0;
}

void dataset_set_regular_markers (DATASET *dset)
{
    dset->markers = REGULAR_MARKERS;
}

struct filetype_info {
    GretlFileType type;
    const char *src;
};

/**
 * dataset_add_import_info:
 * @dset: pointer to dataset information struct.
 * @fname: the name of a file from which data have been imported.
 * @type: code representing the type of the file identified by
 * @fname.
 *
 * On successful import of data from some "foreign" format,
 * add a note to the "descrip" member of the new dataset
 * saying where it came from and when.
 */

void dataset_add_import_info (DATASET *dset, const char *fname,
			      GretlFileType type)
{
    struct filetype_info ftypes[] = {
	{ GRETL_CSV,      "CSV" },
	{ GRETL_GNUMERIC, "Gnumeric" },
	{ GRETL_XLS,      "Excel" },
	{ GRETL_XLSX,     "Excel" },
	{ GRETL_ODS,      "Open Document" },
	{ GRETL_WF1,      "Eviews" },
	{ GRETL_DTA,      "Stata" },
	{ GRETL_SAV,      "SPSS" },
	{ GRETL_SAS,      "SAS" },
	{ GRETL_JMULTI,   "JMulTi" }
    };
    int i, nt = G_N_ELEMENTS(ftypes);
    const char *src = NULL;
    gchar *note = NULL;
    char tstr[128];

    for (i=0; i<nt; i++) {
	if (type == ftypes[i].type) {
	    src = ftypes[i].src;
	    break;
	}
    }

    if (src == NULL) {
	return;
    }

    print_time(tstr);

    if (g_utf8_validate(fname, -1, NULL)) {
	const char *p = strrslash(fname);

	if (p != NULL) {
	    fname = p + 1;
	}
	note = g_strdup_printf(_("Data imported from %s file '%s', %s\n"),
			       src, fname, tstr);
    } else {
	note = g_strdup_printf(_("Data imported from %s, %s\n"),
			       src, tstr);
    }

    if (note != NULL) {
	if (dset->descrip == NULL) {
	    dset->descrip = gretl_strdup(note);
	} else {
	    int dlen = strlen(dset->descrip);
	    int nlen = strlen(note);
	    char *tmp = realloc(dset->descrip, dlen + nlen + 5);

	    if (tmp != NULL) {
		dset->descrip = tmp;
		strcat(dset->descrip, "\n\n");
		strncat(dset->descrip, note, nlen);
	    }
	}
	g_free(note);
    }
}

static int is_weekend (int t, int pd, int sat0, int sun0)
{
    int sat = 0, sun = 0;

    /* this is intended to identify weekend days for
       both 7- and 6-day data */

    if (sat0 >= 0) {
	sat = (t - sat0) % pd == 0;
    }

    if (!sat && sun0 >= 0) {
	sun = (t - sun0) % pd == 0;
    }

    return sat || sun;
}

/* Scan imported daily data for missing values, so as to
   be able to offer the user some options.

   Return values:

   0 : no missing values
   1 : all weekend data are missing (or possibly just all
       Saturdays, or just all Sundays), but no weekday
       data missing
   2 : as with 1, but also some weekdays missing
   3 : scattering of weekend and/or weekday data missing
*/

int analyse_daily_import (const DATASET *dset, PRN *prn)
{
    int all_weekends_blank = 0;
    int blank_weekends = 0;
    int blank_weekdays = 0;
    int n_weekdays = 0;
    int n_weekend_days = 0;
    int sat0 = -1, sun0 = -1;
    int i, t, pd = dset->pd;
    int all_missing, weekend;
    int ret = 0;

    if (pd > 5) {
	char datestr[OBSLEN];
	int wkday;

	/* start by finding first Sat and/or Sun */
	for (t=0; t<dset->n; t++) {
	    ntolabel(datestr, t, dset);
	    wkday = weekday_from_date(datestr);
	    if (wkday == 6 && sat0 < 0) {
		sat0 = t;
	    } else if (wkday == 7 && sun0 < 0) {
		sun0 = t;
	    }
	    if (sat0 >= 0 && sun0 >= 0) {
		break;
	    } else if (dset->pd == 6 && (sat0 >= 0 || sun0 >= 0)) {
		break;
	    }
	}
	all_weekends_blank = 1; /* may be revised below */
    } else {
	/* there are no weekend days in 5-day data */
	weekend = 0;
    }

    for (t=0; t<dset->n; t++) {
	if (pd > 5) {
	    weekend = is_weekend(t, pd, sat0, sun0);
	}
	all_missing = 1;
	for (i=1; i<dset->v; i++) {
	    if (!na(dset->Z[i][t])) {
		all_missing = 0;
		break;
	    }
	}
	if (weekend) {
	    n_weekend_days++;
	    if (!all_missing) {
		/* not all weekend data are missing */
		all_weekends_blank = 0;
	    }
	} else {
	    n_weekdays++;
	}
	if (all_missing) {
	    if (weekend) {
		blank_weekends++;
	    } else {
		blank_weekdays++;
	    }
	}
    }

    if (all_weekends_blank) {
	double misspc = 100.0 * blank_weekdays / (double) n_weekdays;

	if (pd == 7) {
	    pputs(prn, _("This dataset is on 7-day calendar, but weekends are blank."));
	} else {
	    pprintf(prn, _("This dataset is on 6-day calendar, but %s are blank."),
		    sat0 >= 0 ? _("Sundays") : _("Saturdays"));
	}
	ret = 1;
	if (misspc > 0.0) {
	    pputc(prn, '\n');
	    pprintf(prn, "%s ", _("In addition,"));
	    if (misspc >= 0.01) {
		pprintf(prn, _("%.2f percent of weekday observations are missing."),
			misspc);
	    } else {
		pprintf(prn, _("%g percent of weekday observations are missing."),
			misspc);
	    }
	    if (misspc < 10.0) {
		ret = 2;
	    }
	}
    } else if (blank_weekdays || blank_weekends) {
	int ndays = n_weekdays + n_weekend_days;
	int nmiss = blank_weekdays + blank_weekends;
	double misspc = 100.0 * nmiss / (double) ndays;

	if (misspc >= 0.01) {
	    pprintf(prn, _("%.2f percent of daily observations are missing."),
		    misspc);
	} else {
	    pprintf(prn, _("%g percent of daily observations are missing."),
		    misspc);
	}
	if (misspc < 10) {
	    ret = 3;
	}
    }

    return ret;
}
