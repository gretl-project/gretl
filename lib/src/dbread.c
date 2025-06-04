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

/* dbread.c for gretl */

#include "libgretl.h"
#include "swap_bytes.h"
#include "libset.h"
#include "uservar.h"
#include "matrix_extra.h"
#include "usermat.h"
#include "dbread.h"
#include "gretl_midas.h"
#include "gretl_typemap.h"
#ifdef USE_CURL
# include "gretl_www.h"
#endif

#include <glib.h>
#include <unistd.h>
#include <errno.h>

#if G_BYTE_ORDER == G_BIG_ENDIAN
# include <netinet/in.h>
#endif

/**
 * SECTION:dbread
 * @short_description: reading from databases
 * @title: DB read
 * @include: gretl/libgretl.h, gretl/dbread.h
 *
 * Functions that read data from native gretl databases as
 * well as RATS 4.0 and PcGive databases. As you will see,
 * this area is mostly undocumented at present, but since it
 * may ultimately be useful for third-party coders we will
 * try to remedy this!
 */

#define DB_DEBUG 0

#define RECNUM gint32
#define NAMELENGTH 16
#define RATSCOMMENTLENGTH 80
#define RATSCOMMENTS 2
#define RATS_PARSE_ERROR -999

typedef struct {
    gint32 daynumber;              /* Number of days from 1-1-90
				      to year, month, day */
    short panel;                   /* 1 for panel set, 2 for intraday
				      date set , 0 o.w. */
#define LINEAR    0                /* Single time direction */
#define RATSPANEL 1                /* panel:period */
#define INTRADAY  2                /* date:intraday period */
    gint32 panelrecord;            /* Size of panel or
				      number of periods per day */
    short dclass;                  /* See definitions below */
#define UNDATEDCLASS   0           /* No time series properties */
#define IRREGULARCLASS 1           /* Time series (irregular) */
#define PERYEARCLASS   2           /* x periods / year */
#define PERWEEKCLASS   3           /* x periods / week */
#define DAILYCLASS     4           /* x days / period */
    gint32 info;                   /* Number of periods per year or
				      per week */
    short digits;                  /* Digits for representing panel
				      or intraday period */
    short year,month,day;          /* Starting year, month, day */
} DATEINFO;

typedef struct {
    RECNUM back_point;             /* Pointer to previous series */
    RECNUM forward_point;          /* Pointer to next series */
    short back_class;              /* Reserved.  Should be 0 */
    short forward_class;           /* Reserved.  Should be 0 */
    RECNUM first_data;             /* First data record */
    char series_name[NAMELENGTH];  /* Series name */
    DATEINFO date_info;            /* Dating scheme for this series */
    gint32 datapoints;             /* Number of data points */
    short data_type;               /* real, char, complex.
                                      Reserved.  Should be 0 */
    short digits;                  /* . + digit count for representation
				      (0 = unspecified) */
    short misc1;                   /* For future expansion should be 0 */
    short misc2;
    short comment_lines;           /* Number of comment lines (0,1,2) */
    char series_class[NAMELENGTH]; /* Series class. Not used, blank */
    char comments[RATSCOMMENTS][RATSCOMMENTLENGTH];
    char pad[10];
} RATSDirect;

typedef struct {
    RECNUM back_point;             /* Previous record (0 for first) */
    RECNUM forward_point;          /* Next record (0 for last) */
    double data[31];               /* Data */
} RATSData;

static char saved_db_name[MAXLEN];
static int saved_db_type;

#if G_BYTE_ORDER == G_BIG_ENDIAN
float retrieve_float (netfloat nf)
{
    short exp = ntohs(nf.exp);
    long frac = ntohl(nf.frac);
    double receive = frac / 10e6;

    return ldexp(receive, exp);
}
#endif

static int lib_add_db_data (double **dbZ, SERIESINFO *sinfo,
			    DATASET *dset, CompactMethod cmethod,
			    int dbv, PRN *prn);

static int do_compact_spread (DATASET *dset, int newpd);

static FILE *open_binfile (const char *dbbase, int code, int offset, int *err)
{
    char dbbin[MAXLEN];
    FILE *fp = NULL;

    strcpy(dbbin, dbbase);
    if (code == GRETL_NATIVE_DB) {
	if (strstr(dbbin, ".bin") == NULL) {
	    strcat(dbbin, ".bin");
	}
    } else {
	if (strstr(dbbin, ".bn7") == NULL) {
	    strcat(dbbin, ".bn7");
	}
    }

    fp = gretl_fopen(dbbin, "rb");

    if (fp == NULL) {
	*err = E_FOPEN;
    } else if (fseek(fp, (long) offset, SEEK_SET)) {
	*err = DB_PARSE_ERROR;
	fclose(fp);
	fp = NULL;
    }

    return fp;
}

/**
 * get_native_db_data:
 * @dbbase:
 * @sinfo:
 * @Z: data array.
 *
 * Returns: 0 on success, non-zero code on failure.
 */

int get_native_db_data (const char *dbbase, SERIESINFO *sinfo,
			double **Z)
{
    char numstr[32];
    FILE *fp;
    dbnumber x;
    int v = sinfo->v;
    int t, t2, err = 0;

    fp = open_binfile(dbbase, GRETL_NATIVE_DB, sinfo->offset, &err);
    if (err) {
	return err;
    }

    t2 = (sinfo->t2 > 0)? sinfo->t2 : sinfo->nobs - 1;

    for (t=sinfo->t1; t<=t2 && !err; t++) {
	if (fread(&x, sizeof x, 1, fp) != 1) {
	    err = DB_PARSE_ERROR;
	} else {
	    sprintf(numstr, "%.7g", (double) x); /* N.B. converting a float */
	    Z[v][t] = atof(numstr);
	    if (Z[v][t] == DBNA) {
		Z[v][t] = NADBL;
	    }
	}
    }

    fclose(fp);

    return err;
}

#ifdef USE_CURL

/**
 * get_remote_db_data:
 * @dbbase:
 * @sinfo:
 * @Z: data array.
 *
 * Returns: 0 on success, non-zero code on failure.
 */

int get_remote_db_data (const char *dbbase, SERIESINFO *sinfo,
			double **Z)
{
    char *getbuf = NULL;
    int t, t2, err;
    int v = sinfo->v;
    dbnumber x;
    size_t offset;
#if G_BYTE_ORDER == G_BIG_ENDIAN
    netfloat nf;
#endif

    err = retrieve_remote_db_data(dbbase, sinfo->varname, &getbuf);
    if (err) {
	free(getbuf);
	return E_FOPEN;
    }

    t2 = (sinfo->t2 > 0)? sinfo->t2 : sinfo->nobs - 1;

    offset = 0L;
    for (t=sinfo->t1; t<=t2; t++) {
#if G_BYTE_ORDER == G_BIG_ENDIAN
	/* go via network byte order */
	memcpy(&(nf.frac), getbuf + offset, sizeof nf.frac);
	offset += sizeof nf.frac;
	memcpy(&(nf.exp), getbuf + offset, sizeof nf.exp);
	offset += sizeof nf.exp;
	x = retrieve_float(nf);
#else
	/* just read floats */
	memcpy(&x, getbuf + offset, sizeof x);
	offset += sizeof x;
#endif
	Z[v][t] = (x == DBNA)? NADBL : x;
    }

    free(getbuf);

    return 0;
}

#endif /* USE_CURL */

/**
 * get_pcgive_db_data:
 * @dbbase:
 * @sinfo:
 * @Z: data array.
 *
 *
 * Returns: 0 on success, non-zero code on failure.
 */

int get_pcgive_db_data (const char *dbbase, SERIESINFO *sinfo,
			double **Z)
{
    FILE *fp;
    double x;
    int v = sinfo->v;
    int t, t2, err = 0;

    fp = open_binfile(dbbase, GRETL_PCGIVE_DB, sinfo->offset, &err);
    if (err) {
	return err;
    }

    t2 = (sinfo->t2 > 0)? sinfo->t2 : sinfo->nobs - 1;

    for (t=sinfo->t1; t<=t2; t++) {
	if (fread(&x, sizeof x, 1, fp) != 1) {
	    err = E_DATA;
	    break;
	}
#if G_BYTE_ORDER == G_BIG_ENDIAN
	reverse_double(x);
#endif
	if (x == -9999.99 || isnan(x)) {
	    Z[v][t] = NADBL;
	    err = DB_MISSING_DATA;
	} else {
	    Z[v][t] = x;
	}
    }

    fclose(fp);

    return err;
}

static void get_native_series_comment (SERIESINFO *sinfo, const char *s)
{
    s += strcspn(s, " "); /* skip varname */
    s += strspn(s, " ");  /* skip space */

    series_info_set_description(sinfo, s);
    tailstrip(sinfo->descrip);
}

static int get_native_series_pd (SERIESINFO *sinfo, char pdc)
{
    sinfo->pd = 1;
    sinfo->undated = 0;

    if (pdc == 'M') sinfo->pd = 12;
    else if (pdc == 'Q') sinfo->pd = 4;
    else if (pdc == 'B') sinfo->pd = 5;
    else if (pdc == 'S') sinfo->pd = 6;
    else if (pdc == 'D') sinfo->pd = 7;
    else if (pdc == 'U') sinfo->undated = 1;
    else return 1;

    return 0;
}

static int get_native_series_obs (SERIESINFO *sinfo,
				  const char *stobs,
				  const char *endobs)
{
    char dc = 0;

    if (strchr(stobs, '-')) {
	dc = '-';
    } else if (strchr(stobs, '/')) {
	dc = '/';
    }

    if (dc != 0) {
	/* calendar data */
	const char *q = stobs;
	const char *p = strchr(stobs, dc);

	if (p - q == 4) {
	    strcpy(sinfo->stobs, q);
	}
	q = endobs;
	p = strchr(endobs, dc);
	if (p && p - q == 4) {
	    strcpy(sinfo->endobs, q);
	}
    } else {
	*sinfo->stobs = '\0';
	*sinfo->endobs = '\0';
	strncat(sinfo->stobs, stobs, OBSLEN-1);
	strncat(sinfo->endobs, endobs, OBSLEN-1);
    }

    return 0;
}

static int
open_native_db_files (const char *dname, FILE **f1, char *name1,
		      FILE **f2, char *name2)
{
    char dbbase[FILENAME_MAX];
    char fname[FILENAME_MAX];
    FILE *fidx = NULL, *fbin = NULL;
    int err = 0;

    if (dname != NULL) {
	strcpy(dbbase, dname);
    } else {
	strcpy(dbbase, saved_db_name);
    }

    if (has_suffix(dbbase, ".bin")) {
	dbbase[strlen(dbbase) - 4] = '\0';
    }

    if (f1 != NULL) {
	strcpy(fname, dbbase);
	strcat(fname, ".idx");

	if (name1 != NULL) {
	    err = gretl_write_access(fname);
	    if (!err) {
		strcpy(name1, fname);
	    }
	}

	if (!err) {
	    fidx = gretl_fopen(fname, "r");
	    if (fidx == NULL) {
		gretl_errmsg_set(_("Couldn't open database index file"));
		err = E_FOPEN;
	    }
	}
    }

    if (f2 != NULL && !err) {
	strcpy(fname, dbbase);
	strcat(fname, ".bin");

	if (name2 != NULL) {
	    err = gretl_write_access(fname);
	    if (!err) {
		strcpy(name2, fname);
	    }
	}

	if (!err) {
	    fbin = gretl_fopen(fname, "rb");
	    if (fbin == NULL) {
		gretl_errmsg_set(_("Couldn't open database binary file"));
		err = E_FOPEN;
	    }
	}
    }

    if (err) {
	if (fidx != NULL) {
	    fclose(fidx);
	}
    } else {
	if (f1 != NULL) {
	    *f1 = fidx;
	}
	if (f2 != NULL) {
	    *f2 = fbin;
	}
    }

    return err;
}

static char *native_db_index_name (void)
{
    char *fname;

    if (has_suffix(saved_db_name, ".bin")) {
	fname = g_strdup(saved_db_name);
	strcpy(fname + strlen(fname) - 3, "idx");
    } else {
	fname = g_strdup_printf("%s.idx", saved_db_name);
    }

    return fname;
}

static int db_match_glob (FILE *fp,
			  char *line, int linelen,
			  GPatternSpec *pspec,
			  char **S, int *err)
{
    char vname[VNAMELEN], l2[72];
    int n = 0;

    while (fgets(line, linelen, fp) && !*err) {
	if (*line == '#') {
	    continue;
	}
	if (gretl_scan_varname(line, vname) != 1) {
	    break;
	}
	if (g_pattern_match_string(pspec, vname)) {
	    if (S != NULL) {
		S[n] = gretl_strdup(vname);
	    }
	    n++;
	}
	if (fgets(l2, sizeof l2, fp) == NULL) {
	    *err = DB_PARSE_ERROR;
	}
    }

    return n;
}

static char **native_db_match_series (const char *glob, int *nmatch,
				      const char *idxname, int *err)
{
    GPatternSpec *pspec;
    char **S = NULL;
    FILE *fp = NULL;
    char line[256];

    fp = gretl_fopen(idxname, "rb");
    if (fp == NULL) {
	*err = E_FOPEN;
	*nmatch = 0;
	return NULL;
    }

    pspec = g_pattern_spec_new(glob);

    *nmatch = db_match_glob(fp, line, sizeof line, pspec, NULL, err);

    if (!*err && *nmatch > 0) {
	S = strings_array_new(*nmatch);
	if (S == NULL) {
	    *nmatch = 0;
	    *err = E_ALLOC;
	} else {
	    rewind(fp);
	    db_match_glob(fp, line, sizeof line, pspec, S, err);
	}
    }

    g_pattern_spec_free(pspec);
    fclose(fp);

    return S;
}

static int get_native_series_info (const char *series,
				   SERIESINFO *sinfo,
				   const char *idxname)
{
    FILE *fp = NULL;
    char sername[VNAMELEN];
    /* 2019-01-08: enlarge @s1 from 256 to 1024 */
    char s1[1024], s2[72];
    char stobs[OBSLEN], endobs[OBSLEN];
    char pdc;
    int offset = 0;
    int gotit = 0, err = 0;
    int n;

    fp = gretl_fopen(idxname, "rb");
    if (fp == NULL) {
	return E_FOPEN;
    }

    while (fgets(s1, sizeof s1, fp) && !gotit) {
	if (*s1 == '#') {
	    continue;
	}
	if (gretl_scan_varname(s1, sername) != 1) {
	    break;
	}
	if (!strcmp(series, sername)) {
	    gotit = 1;
	    strcpy(sinfo->varname, sername);
	}
	if (fgets(s2, sizeof s2, fp) == NULL) {
	    err = DB_PARSE_ERROR;
	    break;
	}
	if (gotit) {
	    get_native_series_comment(sinfo, s1);
	    if (sscanf(s2, "%c %10s %*s %10s %*s %*s %d",
		       &pdc, stobs, endobs, &sinfo->nobs) != 4) {
		gretl_errmsg_set(_("Failed to parse series information"));
		err = DB_PARSE_ERROR;
	    } else {
		get_native_series_pd(sinfo, pdc);
		get_native_series_obs(sinfo, stobs, endobs);
		sinfo->offset = offset;
		sinfo->t2 = sinfo->nobs - 1;
	    }
	} else {
	    if (sscanf(s2, "%*c %*s %*s %*s %*s %*s %d", &n) != 1) {
		gretl_errmsg_set(_("Failed to parse series information"));
		err = DB_PARSE_ERROR;
	    } else {
		offset += n * sizeof(dbnumber);
	    }
	}
    }

    fclose(fp);

    if (!err && !gotit) {
	gretl_errmsg_sprintf(_("Series not found, '%s'"), series);
	err = DB_NO_SUCH_SERIES;
    }

    return err;
}

#ifdef USE_CURL

static int remote_db_index_to_file (const char *fname)
{
    char *buf = NULL;
    int err;

    err = retrieve_remote_db_index(saved_db_name, &buf);

    if (!err) {
	FILE *fp = gretl_fopen(fname, "wb");

	if (fp == NULL) {
	    err = E_FOPEN;
	} else {
	    fputs(buf, fp);
	    fclose(fp);
#if 1
	    fprintf(stderr, "remote db index saved\n");
#endif
	}
	free(buf);
    }

    return err;
}

#endif /* USE_CURL */

static int in7_get_obs (int y0, int p0, int y1, int p1,
			SERIESINFO *sinfo)
{
    int pd = sinfo->pd;
    int n = (y1 - y0 + 1) * pd - (p0 - 1) - (pd - p1);
    int err = 0;

    if (n <= 0) {
	err = 1;
    } else {
	sinfo->nobs = n;
	sinfo->t2 = n - 1;
    }

    return err;
}

static int
pcgive_set_stobs_endobs (int y0, int p0, int y1, int p1,
			 SERIESINFO *sinfo)
{
    int err = 0;

    if (sinfo->pd == 1) {
	sprintf(sinfo->stobs, "%d", y0);
	sprintf(sinfo->endobs, "%d", y1);
	if (y0 == 1) {
	    sinfo->undated = 1;
	}
    } else if (sinfo->pd == 4) {
	sprintf(sinfo->stobs, "%d:%d", y0, p0);
	sprintf(sinfo->endobs, "%d:%d", y1, p1);
    } else if (sinfo->pd == 12 || sinfo->pd == 52) {
	sprintf(sinfo->stobs, "%d:%02d", y0, p0);
	sprintf(sinfo->endobs, "%d:%02d", y1, p1);
    } else {
	err = E_DATA; /* FIXME? */
    }

    return err;
}

static int
get_pcgive_series_info (const char *series, SERIESINFO *sinfo)
{
    FILE *fp;
    char dbidx[MAXLEN];
    char line[1024];
    char fmt[24];
    char *p;
    int y0, p0, y1, p1;
    int nf, gotit = 0;
    int err = 0;

    strcpy(dbidx, saved_db_name);
    p = strstr(dbidx, ".bn7");
    if (p != NULL) {
	strcpy(p, ".in7");
    } else {
	strcat(dbidx, ".in7");
    }

#if DB_DEBUG
    fprintf(stderr, "get_pcgive_series_info: dbidx = '%s'\n", dbidx);
#endif

    fp = gretl_fopen(dbidx, "r");
    if (fp == NULL) {
	gretl_errmsg_set(_("Couldn't open database index file"));
	return E_FOPEN;
    }

    sprintf(fmt, "%%%ds %%d %%d %%d %%d %%d %%d", VNAMELEN - 1);

    while (fgets(line, sizeof line, fp) && !gotit) {
	if (*line == '>') {
	    *sinfo->varname = 0;
	    nf = sscanf(line + 1, fmt, sinfo->varname, &y0, &p0,
			&y1, &p1, &sinfo->pd, &sinfo->offset);
	    fprintf(stderr, "in7: varname='%s'\n", sinfo->varname);
	    if (!strcmp(sinfo->varname, series)) {
		gotit = 1;
	    } else {
		continue;
	    }
	    if (nf == 7 && y0 > 0 && p0 > 0 && y1 > 0 && p1 > 0 &&
		sinfo->pd >= 1 && sinfo->offset > 0) {
		while (fgets(line, sizeof line, fp)) {
		    if (*line == ';') {
			gretl_strstrip(line);
			series_info_set_description(sinfo, line + 1);
		    } else {
			break;
		    }
		}
		/* transcribe info */
		err = in7_get_obs(y0, p0, y1, p1, sinfo);
		if (!err) {
		    err = pcgive_set_stobs_endobs(y0, p0, y1, p1, sinfo);
		}
	    } else {
		err = E_DATA;
	    }
	}
    }

    fclose(fp);

    if (!err && !gotit) {
	gretl_errmsg_sprintf(_("Series not found, '%s'"), series);
	err = DB_NO_SUCH_SERIES;
    }

    return err;
}

/* Figure the ending observation date of a series */

static int get_endobs (char *datestr, int startyr, int startfrac,
		       int pd, int n)
{
    int endyr, endfrac;

    endyr = startyr + n / pd;
    endfrac = startfrac - 1 + n % pd;

    if (endfrac >= pd) {
	endyr++;
	endfrac -= pd;
    }

    if (endfrac == 0) {
	endyr--;
	endfrac = pd;
    }

    if (pd == 1) {
	sprintf(datestr, "%d", endyr);
    } else if (pd == 4) {
	sprintf(datestr, "%d.%d", endyr, endfrac);
    } else if (pd == 12 || pd == 52) {
	sprintf(datestr, "%d.%02d", endyr, endfrac);
    }

    return 0;
}

static int dinfo_sanity_check (const DATEINFO *dinfo)
{
    int err = 0;

    if (dinfo->info < 0 || dinfo->info > 365) {
	err = 1;
    } else if (dinfo->day < 0 || dinfo->day > 365) {
	err = 1;
    } else if (dinfo->year < 0 || dinfo->year > 3000) {
	err = 1;
    } else if (dinfo->info == 52) {
	/* note: "month" = week */
	if (dinfo->month < 0 || dinfo->month > 52) {
	    err = 1;
	}
    } else {
	/* annual, quarterly, monthly */
	if (dinfo->month < 0 || dinfo->month > 12) {
	    err = 1;
	}
    }

    if (err) {
	gretl_errmsg_set(_("This is not a valid RATS 4.0 database"));
	fprintf(stderr, "rats database: failed dinfo_sanity_check:\n"
		" info=%d, year=%d, month=%d, day=%d\n",
		(int) dinfo->info, (int) dinfo->year, (int) dinfo->month,
		(int) dinfo->day);
    }

    return err;
}

static int dinfo_to_sinfo (const DATEINFO *dinfo, SERIESINFO *sinfo,
			   const char *varname, const char *comment,
			   int n, int offset)
{
    int startfrac = 0;
    char pdstr[8] = {0};
    int err = 0;

    if (dinfo_sanity_check(dinfo)) {
	return 1;
    }

    sprintf(sinfo->stobs, "%d", dinfo->year);

    if (dinfo->info == 4) {
	sprintf(pdstr, ".%d", dinfo->month);
	if (dinfo->month == 1) {
	    startfrac = 1;
	} else if (dinfo->month > 1 && dinfo->month <= 4) {
	    startfrac = 2;
	} else if (dinfo->month > 4 && dinfo->month <= 7) {
	    startfrac = 3;
	} else {
	    startfrac = 4;
	}
    } else if (dinfo->info == 12 || dinfo->info == 52) {
	sprintf(pdstr, ".%02d", dinfo->month);
	startfrac = dinfo->month;
    } else if (dinfo->info == 1) {
	startfrac = 0;
    } else {
	fprintf(stderr, "frequency (%d) does not seem to make sense\n",
		(int) dinfo->info);
	gretl_errmsg_sprintf(_(("frequency (%d) does not seem to make sense")),
			     (int) dinfo->info);
	err = 1;
    }

    if (*pdstr) {
	strcat(sinfo->stobs, pdstr);
    }

    get_endobs(sinfo->endobs, dinfo->year, startfrac, dinfo->info, n);

    sinfo->pd = dinfo->info;
    sinfo->nobs = n;
    sinfo->t2 = n - 1;
    sinfo->offset = offset;

    strncat(sinfo->varname, varname, VNAMELEN - 1);
    series_info_set_description(sinfo, comment);

#if DB_DEBUG
    fprintf(stderr, "dinfo_to_sinfo: '%s': set sinfo->offset = %d\n", varname,
	    (int) offset);
#endif

    return err;
}

static int in7_to_sinfo (const char *varname, const char *comment,
			 int y0, int p0, int y1, int p1, int pd,
			 int offset, SERIESINFO *sinfo)
{
    int err = 0;

    if (pd == 4) {
	sprintf(sinfo->stobs, "%d.%d", y0, p0);
	sprintf(sinfo->endobs, "%d.%d", y1, p1);
    } else if (pd == 12 || pd == 52) {
	sprintf(sinfo->stobs, "%d.%02d", y0, p0);
	sprintf(sinfo->endobs, "%d.%02d", y1, p1);
    } else if (pd == 1) {
	sprintf(sinfo->stobs, "%d", y0);
	sprintf(sinfo->endobs, "%d", y1);
    } else {
	fprintf(stderr, "frequency %d is not supported\n", pd);
	gretl_errmsg_sprintf(_("frequency %d is not supported"), pd);
	err = 1;
    }

    if (!err) {
	sinfo->pd = pd;
	err = in7_get_obs(y0, p0, y1, p1, sinfo);
    }

    if (!err) {
	strcpy(sinfo->varname, varname);
	if (comment != NULL && *comment != 0) {
	    series_info_set_description(sinfo, comment);
	}
	sinfo->pd = pd;
	sinfo->offset = offset;
    }

    return err;
}

static RECNUM read_rats_directory (FILE *fp, const char *series_name,
				   SERIESINFO *sinfo)
{
    RATSDirect rdir;
    DATEINFO dinfo;
    RECNUM ret;
    int nread;
    int i, err = 0;

    memset(rdir.series_name, 0, NAMELENGTH);

    if (fread(&rdir.back_point, sizeof(RECNUM), 1, fp) != 1) {
	err = 1;
    } else if (fread(&rdir.forward_point, sizeof(RECNUM), 1, fp) != 1) {
	err = 1;
    }
    if (!err) {
	fseek(fp, 4L, SEEK_CUR); /* skip two shorts */
	if (fread(&rdir.first_data, sizeof(RECNUM), 1, fp) != 1) {
	    err = 1;
	} else if (fread(rdir.series_name, NAMELENGTH, 1, fp) != 1) {
	    err = 1;
	}
    }

    if (!err) {
	rdir.series_name[NAMELENGTH-1] = '\0';
	gretl_strstrip(rdir.series_name);
#if DB_DEBUG
	fprintf(stderr, "read_rats_directory: name='%s'\n", rdir.series_name);
#endif
	if (!isprint(rdir.series_name[0])) {
	    err = 1;
	}
    }

    if (err) {
	return RATS_PARSE_ERROR;
    }

    if (series_name != NULL && strcmp(series_name, rdir.series_name)) {
	/* specific series not found yet: keep going */
	return rdir.forward_point;
    }

    /* Now the dateinfo: we can't read this in one go either :-( */

    /* skip long, short, long, short */
    fseek(fp, 12, SEEK_CUR);
    nread = 0;
    nread += fread(&dinfo.info, sizeof(gint32), 1, fp);
    nread += fread(&dinfo.digits, sizeof(short), 1, fp);
    nread += fread(&dinfo.year, sizeof(short), 1, fp);
    nread += fread(&dinfo.month, sizeof(short), 1, fp);
    nread += fread(&dinfo.day, sizeof(short), 1, fp);
    nread += fread(&rdir.datapoints, sizeof(gint32), 1, fp);

    if (nread != 6) {
	return RATS_PARSE_ERROR;
    }

    fseek(fp, sizeof(short) * 4L, SEEK_CUR);  /* skip 4 shorts */

#if DB_DEBUG
    fprintf(stderr, "info=%d, digits=%d, year=%d, mon=%d, day=%d\n",
	    (int) dinfo.info, (int) dinfo.digits, (int) dinfo.year,
	    (int) dinfo.month, (int) dinfo.day);
    fprintf(stderr, "datapoints = %d\n", (int) rdir.datapoints);
#endif

    if (fread(&rdir.comment_lines, sizeof(short), 1, fp) != 1) {
	err = 1;
    } else {
	fseek(fp, 1L, SEEK_CUR); /* skip one char */
	for (i=0; i<2 && !err; i++) {
	    if (i < rdir.comment_lines) {
		memset(rdir.comments[i], 0, 80);
		err = (fread(rdir.comments[i], 80, 1, fp) != 1);
		if (!err) {
		    rdir.comments[i][79] = '\0';
		    gretl_strstrip(rdir.comments[i]);
		}
	    } else {
		rdir.comments[i][0] = 0;
		fseek(fp, 80, SEEK_CUR);
	    }
	}
    }

#if DB_DEBUG
    if (!err) {
	fprintf(stderr, "comment_lines = %d\n", (int) rdir.comment_lines);
	fprintf(stderr, "comment[0] = '%s'\n", rdir.comments[0]);
	fprintf(stderr, "comment[1] = '%s'\n", rdir.comments[1]);
    }
#endif

    if (!err) {
	err = dinfo_to_sinfo(&dinfo, sinfo, rdir.series_name, rdir.comments[0],
			     rdir.datapoints, rdir.first_data);
    }

    ret = (err)? RATS_PARSE_ERROR : rdir.forward_point;

#if DB_DEBUG
    fprintf(stderr, "read_rats_directory: err = %d, forward_point=%d, first_data=%d\n",
	    err, (int) rdir.forward_point, (int) rdir.first_data);
    fprintf(stderr, "returning %d\n", (int) ret);
#endif

    return ret;
}

static void series_info_init (SERIESINFO *sinfo)
{
    sinfo->t1 = sinfo->t2 = 0;
    sinfo->nobs = 0;
    sinfo->v = 1;
    sinfo->pd = 1;
    sinfo->offset = -1;
    sinfo->err = 0;
    sinfo->undated = 0;

    sinfo->varname[0] = '\0';
    sinfo->stobs[0] = '\0';
    sinfo->endobs[0] = '\0';

    sinfo->descrip = NULL;
    sinfo->data = NULL;
}

void series_info_set_description (SERIESINFO *sinfo,
				  const char *s)
{
    if (sinfo->descrip != NULL) {
	free(sinfo->descrip);
	sinfo->descrip = NULL;
    }
    if (s != NULL && *s != '\0') {
	sinfo->descrip = gretl_strdup(s);
    }
}

static void series_info_clear (SERIESINFO *sinfo)
{
    free(sinfo->descrip);
}

#define DB_INIT_ROWS 32

/**
 * dbwrapper_destroy:
 * @dw: database series wrapper.
 *
 * Frees all resources associated with @dw as well as the pointer
 * itself.
 */

void dbwrapper_destroy (dbwrapper *dw)
{
    if (dw != NULL) {
	free(dw->fname);
	free(dw->sinfo);
	free(dw);
    }
}

/**
 * dbwrapper_new:
 * @n: initial number of series.
 * @fname: database filename.
 * @dbtype: database type code.
 *
 * Returns: an allocated database series wrapper.
 */

dbwrapper *dbwrapper_new (int n, const char *fname, int dbtype)
{
    dbwrapper *dw;
    int i;

    if (n == 0) {
	n = DB_INIT_ROWS;
    }

    dw = malloc(sizeof *dw);
    if (dw == NULL) {
	return NULL;
    }

    dw->fname = gretl_strdup(fname);
    dw->dbtype = dbtype;

    dw->sinfo = malloc(n * sizeof *dw->sinfo);
    if (dw->sinfo == NULL) {
	free(dw);
	return NULL;
    }

    for (i=0; i<n; i++) {
	series_info_init(&dw->sinfo[i]);
    }

    dw->nv = 0;
    dw->nalloc = n;

    return dw;
}

static int dbwrapper_expand (dbwrapper *dw)
{
    SERIESINFO *sinfo;
    int i, newsz;

    newsz = (dw->nv / DB_INIT_ROWS) + 1;
    newsz *= DB_INIT_ROWS;

    sinfo = realloc(dw->sinfo, newsz * sizeof *sinfo);
    if (sinfo == NULL) {
	free(dw->sinfo);
	dw->sinfo = NULL;
	return 1;
    }

    dw->sinfo = sinfo;

    for (i=dw->nalloc; i<newsz; i++) {
	series_info_init(&dw->sinfo[i]);
    }

    dw->nalloc = newsz;

    return 0;
}

static int read_in7_series_info (FILE *fp, dbwrapper *dw)
{
    char line[1024];
    char sname[VNAMELEN];
    char desc[MAXLABEL];
    char fmt[24];
    int y0, p0, y1, p1;
    int pd, offset, pos;
    int i, nf;
    int err = 0;

    sprintf(fmt, "%%%ds %%d %%d %%d %%d %%d %%d", VNAMELEN - 1);

    i = 0;
    while (fgets(line, sizeof line, fp) && !err) {
	if (*line == '>') {
	    nf = sscanf(line + 1, fmt, sname, &y0, &p0, &y1,
			&p1, &pd, &offset);
	    if (nf == 7 && y0 > 0 && p0 > 0 && y1 > 0 && p1 > 0 &&
		pd >= 1 && offset > 0) {
		*desc = 0;
		pos = ftell(fp);
		while (fgets(line, sizeof line, fp)) {
		    if (*line == ';') {
			/* following series description */
			int rem = MAXLABEL - strlen(desc) - 1;

			if (rem > 0) {
			    gretl_strstrip(line);
			    strncat(desc, line + 1, rem);
			}
			pos = ftell(fp);
		    } else {
			/* not a description: throw the line back */
			fseek(fp, pos, SEEK_SET);
			break;
		    }
		}
		/* record info */
		err = in7_to_sinfo(sname, desc, y0, p0, y1, p1,
				   pd, offset, &dw->sinfo[i++]);
		if (!err) {
		    dw->nv += 1;
		}
	    }
	}
    }

    return err;
}

static int count_in7_series (FILE *fp, int *err)
{
    char line[1024];
    char sname[VNAMELEN];
    char fmt[24];
    int y0, p0, y1, p1;
    int pd, offset;
    int nf, i = 0, nseries = 0;

    sprintf(fmt, "%%%ds %%d %%d %%d %%d %%d %%d", VNAMELEN - 1);

    while (fgets(line, sizeof line, fp)) {
	if (i == 0 && strncmp(line, "pcgive 700", 10)) {
	    *err = 1;
	    gretl_errmsg_set(_("This is not a PcGive 700 data file"));
	    return 0;
	}
	if (*line == '>') {
	    nf = sscanf(line + 1, fmt, sname, &y0, &p0, &y1,
			&p1, &pd, &offset);
	    if (nf < 7 || y0 < 0 || p0 < 0 || y1 < 0 || p1 < 0 ||
		pd < 1 || offset < 0) {
		fprintf(stderr, "Error reading series info\n");
	    } else {
		nseries++;
	    }
	}
	i++;
    }

    return nseries;
}

/**
 * read_pcgive_db:
 * @fname: name of database file.
 * @fp: pre-opened stream (caller to close it)
 *
 * Read the series info from a PcGive database, .in7 file
 *
 * Returns: pointer to a #dbwrapper containing the series info,
 * or NULL in case of failure.
 */

dbwrapper *read_pcgive_db (const char *fname, FILE *fp)
{
    dbwrapper *dw;
    int ns, err = 0;

    gretl_error_clear();

    ns = count_in7_series(fp, &err);
    if (ns == 0) {
	if (!err) {
	    gretl_errmsg_set(_("No valid series found"));
	}
	return NULL;
    }

#if DB_DEBUG
    fprintf(stderr, "in7: found %d series\n", ns);
#endif

    /* allocate table for series rows */
    dw = dbwrapper_new(ns, fname, GRETL_PCGIVE_DB);
    if (dw == NULL) {
	gretl_errmsg_set(_("Out of memory!"));
	return NULL;
    }

    rewind(fp);

    /* Go find the series info */
    err = read_in7_series_info(fp, dw);

    if (err) {
	dbwrapper_destroy(dw);
	dw = NULL;
    }

    return dw;
}

/**
 * read_rats_db:
 * @fname: database filename.
 * @fp: pre-opened stream (caller to close it)
 *
 * Read the series info from a RATS 4.0 database: read the base
 * block at offset 0 in the data file, and recurse through the
 * directory entries.
 *
 * Returns: pointer to a #dbwrapper containing the series info,
 * or NULL in case of failure.
 */

dbwrapper *read_rats_db (const char *fname, FILE *fp)
{
    dbwrapper *dw;
    long forward = 0;
    int i, err = 0;

    gretl_error_clear();

    /* get into position */
    fseek(fp, 30L, SEEK_SET); /* skip unneeded fields */
    if (fread(&forward, sizeof forward, 1, fp) == 1) {
	fseek(fp, 4L, SEEK_CUR);
    }

    /* basic check */
    if (forward <= 0) {
	gretl_errmsg_set(_("This is not a valid RATS 4.0 database"));
	fprintf(stderr, "rats database: got forward = %ld\n", forward);
	return NULL;
    }

    /* allocate table for series rows */
    dw = dbwrapper_new(0, fname, GRETL_RATS_DB);
    if (dw == NULL) {
	gretl_errmsg_set(_("Out of memory!"));
	return NULL;
    }

    /* Go find the series */
    i = 0;
    while (forward && !err) {
	dw->nv += 1;
#if DB_DEBUG
	fprintf(stderr, "read_rats_db: forward = %d, nv = %d\n",
		(int) forward, dw->nv);
#endif
	if (dw->nv > 0 && dw->nv % DB_INIT_ROWS == 0) {
	    err = dbwrapper_expand(dw);
	    if (err) {
		gretl_errmsg_set(_("Out of memory!"));
	    }
	}
	if (!err) {
	    err = fseek(fp, (forward - 1) * 256L, SEEK_SET);
	    if (!err) {
		forward = read_rats_directory(fp, NULL, &dw->sinfo[i++]);
		if (forward == RATS_PARSE_ERROR) {
		    err = 1;
		}
	    }
	}
#if DB_DEBUG
	fprintf(stderr, "bottom of loop, err = %d\n", err);
#endif
    }

#if DB_DEBUG
    fprintf(stderr, "read_rats_db: err = %d, dw = %p\n",
	    err, (void *) dw);
#endif

    if (err) {
	dbwrapper_destroy(dw);
	return NULL;
    }

    return dw;
}

/* retrieve the actual data values from the data blocks */

static int get_rats_series (int offset, SERIESINFO *sinfo, FILE *fp,
			    double **Z)
{
    RATSData rdata;
    double x;
    int v = sinfo->v;
    int i, t, T;
    int miss = 0;
    int err = 0;

    fprintf(stderr, "get_rats_series: starting from offset %d\n", offset);

    if (sinfo->t2 > 0) {
	T = sinfo->t2 + 1;
    } else {
	T = sinfo->nobs;
    }

    rdata.forward_point = offset;
    t = sinfo->t1;

    while (rdata.forward_point) {
	fseek(fp, (rdata.forward_point - 1) * 256L, SEEK_SET);
	/* the RATSData struct is actually 256 bytes.  Yay! */
	if (fread(&rdata, sizeof rdata, 1, fp) != 1) {
	    err = E_DATA;
	    break;
	}
	for (i=0; i<31 && t<T; i++) {
	    x = rdata.data[i];
#if G_BYTE_ORDER == G_BIG_ENDIAN
	    reverse_double(x);
#endif
	    if (isnan(x)) {
		x = NADBL;
		miss = 1;
	    }
	    Z[v][t++] = x;
	}
    }

    if (miss && !err) {
	err = DB_MISSING_DATA;
    }

    return err;
}

/**
 * get_rats_db_data:
 * @fname: name of RATS 4.0 database to read from
 * @sinfo: holds info on the given series (input)
 * @Z: data matrix
 *
 * Read the actual data values for a series from a RATS database.
 *
 * Returns: 0 on successful completion, E_FOPEN if
 * the data could not be read, and DB_MISSING_DATA if the
 * data were found but there were some missing values.
 */

int get_rats_db_data (const char *fname, SERIESINFO *sinfo,
		      double **Z)
{
    FILE *fp;
    int err = 0;

    fp = gretl_fopen(fname, "rb");
    if (fp == NULL) {
	err = E_FOPEN;
    } else {
	err = get_rats_series(sinfo->offset, sinfo, fp, Z);
	fclose(fp);
    }

    return err;
}

static int get_rats_series_info (const char *series_name, SERIESINFO *sinfo)
{
    FILE *fp;
    long forward = 0;
    int err = 0;

    gretl_error_clear();

    fp = gretl_fopen(saved_db_name, "rb");
    if (fp == NULL) {
	return E_FOPEN;
    }

#if DB_DEBUG
    fprintf(stderr, "Opened %s\n", saved_db_name);
#endif

    /* get into position */
    fseek(fp, 30L, SEEK_SET);
    if (fread(&forward, sizeof forward, 1, fp) == 1) {
	fseek(fp, 4L, SEEK_CUR);
    }

    /* basic check */
    if (forward <= 0) {
	gretl_errmsg_set(_("This is not a valid RATS 4.0 database"));
	fprintf(stderr, "rats database: got forward = %ld\n", forward);
	return DB_PARSE_ERROR;
    }

    sinfo->offset = 0;

    /* Go find the series */
    while (forward) {
	fseek(fp, (forward - 1) * 256L, SEEK_SET);
	forward = read_rats_directory(fp, series_name, sinfo);
	if (forward == RATS_PARSE_ERROR) {
	    sinfo->offset = -1;
	}
	if (sinfo->offset != 0) {
	    break;
	}
    }

    fclose(fp);

    if (sinfo->offset < 0) {
	err = DB_NO_SUCH_SERIES;
    }

#if DB_DEBUG
    fprintf(stderr, "get_rats_series_info: offset = %d\n", sinfo->offset);
    fprintf(stderr, " pd = %d, nobs = %d\n", sinfo->pd, sinfo->nobs);
#endif

    return err;
}

/* For importation of database series */

static double *get_compacted_xt (const double *src, int n,
				 CompactMethod method,
				 int compfac, int skip)
{
    int p, t;
    double *x;

    x = malloc(n * sizeof *x);
    if (x == NULL) {
	return NULL;
    }

    for (t=0; t<n; t++) {
	p = (t + 1) * compfac;
	x[t] = 0.0;
	if (method == COMPACT_AVG || method == COMPACT_SUM) {
	    int i, st;

	    for (i=1; i<=compfac; i++) {
		st = p - i + skip;
		if (na(src[st])) {
		    x[t] = NADBL;
		    break;
		} else {
		    x[t] += src[st];
		}
	    }
	    if (method == COMPACT_AVG) {
		x[t] /= (double) compfac;
	    }
	} else if (method == COMPACT_EOP) {
	    x[t] = src[p - 1 + skip];
	} else if (method == COMPACT_SOP) {
	    x[t] = src[p - compfac + skip];
	}
    }

    return x;
}

/* Compact a single series from a database, for importation
   into a working dataset of lower frequency.  At present
   this is permitted only for the cases:

     quarterly -> annual
     monthly   -> quarterly
     monthly   -> annual
*/

static double *compact_db_series (const double *src,
				  int pd, int *pnobs,
				  char *stobs,
				  int target_pd,
				  CompactMethod method)
{
    int p0, y0, endskip, goodobs;
    int skip = 0, compfac = pd / target_pd;
    double *x;

    if (target_pd == 1) {
	/* figure the annual dates */
	y0 = atoi(stobs);
	p0 = atoi(stobs + 5);
	if (p0 != 1) {
	    ++y0;
	    skip = compfac - (p0 + 1);
	}
	sprintf(stobs, "%d", y0);
    } else if (target_pd == 4) {
	/* figure the quarterly dates */
	float q;
	int q0;

	y0 = atoi(stobs);
	p0 = atoi(stobs + 5);
	q = 1.0 + p0 / 3.;
	q0 = q + .5;
	skip = ((q0 - 1) * 3) + 1 - p0;
	if (q0 == 5) {
	    y0++;
	    q0 = 1;
	}
	sprintf(stobs, "%d.%d", y0, q0);
    } else {
	return NULL;
    }

    endskip = (*pnobs - skip) % compfac;
    goodobs = (*pnobs - skip - endskip) / compfac;
    *pnobs = goodobs;

#if DB_DEBUG
    fprintf(stderr, "startskip = %d\n", skip);
    fprintf(stderr, "endskip = %d\n", endskip);
    fprintf(stderr, "goodobs = %d\n", goodobs);
    fprintf(stderr, "compfac = %d\n", compfac);
    fprintf(stderr, "starting date = %s\n", stobs);
#endif

    x = get_compacted_xt(src, goodobs, method, compfac, skip);

    return x;
}

#define EXPAND_DEBUG 0

/* Expand a single series from a database, for importation
   into a working dataset of higher frequency.  At present
   this is permitted only for the cases:

   1) annual    -> quarterly
   2) annual    -> monthly
   3) quarterly -> monthly
*/

static double *expand_db_series (const double *src,
				 int pd, int *pnobs,
				 char *stobs,
				 DATASET *dset)
{
    char new_stobs[OBSLEN] = {0};
    int target_pd = dset->pd;
    int oldn = *pnobs;
    int mult, newn;
    double *x = NULL;
    int j, t;
    int err = 0;

    mult = target_pd / pd;
    newn = mult * oldn;

    x = malloc(newn * sizeof *x);
    if (x == NULL) {
	err = E_ALLOC;
    } else {
	int s = 0;

	for (t=0; t<oldn; t++) {
	    for (j=0; j<mult; j++) {
		x[s++] = src[t];
	    }
	}
    }

#if EXPAND_DEBUG
    fprintf(stderr, "expand_db_series 1: mult=%d, newn=%d, stobs='%s'\n",
	    mult, newn, stobs);
#endif

    if (err) {
	return NULL;
    }

    if (pd == 1) {
	strcpy(new_stobs, stobs);
	if (target_pd == 4) {
	    strcat(new_stobs, ":1");
	} else {
	    strcat(new_stobs, ":01");
	}
    } else {
	int yr, qtr, mo;

	if (strchr(stobs, '.')) {
	    sscanf(stobs, "%d.%d", &yr, &qtr);
	} else {
	    sscanf(stobs, "%d:%d", &yr, &qtr);
	}
	mo = (qtr - 1) * 3 + 1;
	sprintf(new_stobs, "%d:%02d", yr, mo);
    }

    /* revise incoming values */
    strcpy(stobs, new_stobs);
    *pnobs = newn;

#if EXPAND_DEBUG
    fprintf(stderr, "expand_db_series 2: pd=%d, stobs='%s'\n",
	    pd, stobs);
#endif

    return x;
}

int set_db_name (const char *fname, int filetype, PRN *prn)
{
    FILE *fp;
    int err = 0;

    *saved_db_name = '\0';
    if (fname != NULL) {
	strncat(saved_db_name, fname, MAXLEN - 1);
    }

    if (filetype == GRETL_DBNOMICS || filetype == 0) {
	saved_db_type = filetype;
	return 0;
    }

    if (filetype == GRETL_NATIVE_DB_WWW) {
#ifdef USE_CURL
	int n = strlen(saved_db_name);

	if (n > 4) {
	    n -= 4;
	    if (!strcmp(saved_db_name + n, ".bin")) {
		saved_db_name[n] = '\0';
	    }
	}
	err = check_remote_db(saved_db_name);
	if (!err) {
	    saved_db_type = filetype;
	    pprintf(prn, "%s\n", saved_db_name);
	}
#else
	pprintf(prn, _("Internet access not supported"));
	pputc(prn, '\n');
	err = E_DATA;
#endif
	return err;
    }

    fp = gretl_fopen(saved_db_name, "rb");

    if (fp == NULL && !g_path_is_absolute(saved_db_name) &&
	filetype == GRETL_NATIVE_DB) {
	/* try looking a bit more */
	const char *path = gretl_binbase();

	if (path != NULL && *path != '\0') {
	    gretl_build_path(saved_db_name, path, fname, NULL);
	    fp = gretl_fopen(saved_db_name, "rb");
	}

#ifdef __APPLE__
	if (fp == NULL) {
	    gchar *tmp = g_build_filename(gretl_app_support_dir(), "db",
					  fname, NULL);

	    fp = gretl_fopen(tmp, "rb");
	    if (fp != NULL) {
		strcpy(saved_db_name, tmp);
	    }
	    g_free(tmp);
	}
#endif
    }

    if (fp == NULL) {
	*saved_db_name = '\0';
	pprintf(prn, _("Couldn't open %s\n"), fname);
	err = E_FOPEN;
    } else {
	fclose(fp);
	saved_db_type = filetype;
	pprintf(prn, "%s\n", saved_db_name);
    }

    return err;
}

const char *get_db_name (void)
{
    return saved_db_name;
}

int dbnomics_is_open (void)
{
    return !strcmp(saved_db_name, "dbnomics");
}

/* Handling of DSN setup for ODBC: grab the dsn, username
   and password strings.
*/

static char *get_dsn_field (const char *tag, const char *src)
{
    const char *p;
    char needle[12];
    char *ret = NULL;

    sprintf(needle, "%s=", tag);
    p = strstr(src, needle);

    if (p != NULL) {
	p += strlen(needle);
	if (*p == '"' || *p == '\'') {
	    ret = gretl_quoted_string_strdup(p, NULL);
	} else {
	    ret = gretl_strndup(p, strcspn(p, " "));
	}
    }

    return ret;
}

static ODBC_info gretl_odinfo;

static void ODBC_info_clear_read (void)
{
    int i;

    free(gretl_odinfo.query);
    gretl_odinfo.query = NULL;

    doubles_array_free(gretl_odinfo.X, gretl_odinfo.nvars);
    gretl_odinfo.X = NULL;

    strings_array_free(gretl_odinfo.S, gretl_odinfo.nrows);
    gretl_odinfo.S = NULL;

    gretl_string_table_destroy(gretl_odinfo.gst);
    gretl_odinfo.gst = NULL;

    for (i=0; i<ODBC_OBSCOLS; i++) {
	gretl_odinfo.coltypes[i] = 0;
    }

    if (gretl_odinfo.fmts != NULL) {
	strings_array_free(gretl_odinfo.fmts, gretl_odinfo.obscols);
	gretl_odinfo.fmts = NULL;
    }

    gretl_odinfo.nrows = 0;
    gretl_odinfo.obscols = 0;
    gretl_odinfo.nvars = 0;
}

static void gretl_odbc_cleanup (void)
{
    free(gretl_odinfo.dsn);
    gretl_odinfo.dsn = NULL;

    free(gretl_odinfo.username);
    gretl_odinfo.username = NULL;

    free(gretl_odinfo.password);
    gretl_odinfo.password = NULL;

    ODBC_info_clear_read();
}

int set_odbc_dsn (const char *line, PRN *prn)
{
    int (*check_dsn) (ODBC_info *);
    char *dbname = NULL;
    char *uname = NULL;
    char *pword = NULL;
    int got_plugin = 0;
    int err = 0;

    gretl_odbc_cleanup();

    dbname = get_dsn_field("dsn", line);
    if (dbname == NULL) {
	pputs(prn, _("You must specify a DSN using 'dsn=...'\n"));
	return E_DATA;
    }

    uname = get_dsn_field("user", line);
    pword = get_dsn_field("password", line);

    gretl_odinfo.dsn = dbname;
    gretl_odinfo.username = uname;
    gretl_odinfo.password = pword;

    gretl_error_clear();

    check_dsn = get_plugin_function("gretl_odbc_check_dsn");

    if (check_dsn == NULL) {
        err = 1;
    } else {
	got_plugin = 1;
        err = (*check_dsn) (&gretl_odinfo);
    }

    if (err) {
	if (!got_plugin) {
	    pprintf(prn, _("Couldn't open the gretl ODBC plugin\n"));
	} else {
	    pprintf(prn, _("Failed to connect to ODBC data source '%s'\n"),
		    gretl_odinfo.dsn);
	}
	gretl_odbc_cleanup();
    } else if (gretl_messages_on()) {
	pprintf(prn, _("Connected to ODBC data source '%s'\n"),
		gretl_odinfo.dsn);
    }

    return err;
}

int db_set_sample (const char *start, const char *stop, DATASET *dset)
{
    int t1 = 0, t2 = 0;

    if (strcmp(start, ";")) {
	t1 = dateton(start, dset);
	if (t1 < 0) {
	    return 1;
	}
    }

    t2 = dateton(stop, dset);
    if (t2 < 0) {
	return 1;
    }

    if (t1 > t2) {
	gretl_errmsg_set(_("Invalid null sample"));
	return 1;
    }

    dset->t1 = t1;
    dset->t2 = t2;
    dset->n = t2 - t1 + 1;
    strcpy(dset->endobs, stop);

#if DB_DEBUG
    fprintf(stderr, "db_set_sample: t1=%d, t2=%d, stobs='%s', endobs='%s' "
	    "sd0 = %g, n = %d\n",
	    dset->t1, dset->t2,
	    dset->stobs, dset->endobs,
	    dset->sd0, dset->n);
#endif

    return 0;
}

static const char *
get_word_and_advance (const char *s, char *word, size_t maxlen)
{
    size_t i = 0;

    while (isspace(*s)) s++;

    *word = '\0';

    while (*s && !isspace(*s)) {
	if (i < maxlen) word[i++] = *s;
	s++;
    }

    word[i] = '\0';

    return (*word != '\0')? s : NULL;
}

static const char *
get_compact_method_and_advance (const char *s, CompactMethod *method)
{
    const char *p;

    *method = COMPACT_UNSET;

    if ((p = strstr(s, "(compact")) != NULL) {
	char comp[8];
	int i = 0;

	p += 8;
	while (*p && *p != ')' && i < 7) {
	    if (!isspace(*p) && *p != '=') {
		comp[i++] = *p;
	    }
	    p++;
	}
	comp[i] = '\0';

	if (!strcmp(comp, "average")) {
	    *method = COMPACT_AVG;
	} else if (!strcmp(comp, "sum")) {
	    *method = COMPACT_SUM;
	} else if (!strcmp(comp, "first")) {
	    *method = COMPACT_SOP;
	} else if (!strcmp(comp, "last")) {
	    *method = COMPACT_EOP;
	} else if (!strcmp(comp, "spread")) {
	    *method = COMPACT_SPREAD;
	}

	p = strchr(p, ')');
	if (p != NULL) p++;
    } else if ((p = strstr(s, "data ")) != NULL) {
	p += 5;
    } else {
	p = s;
    }

    return p;
}

static CompactMethod compact_method_from_option (int *err)
{
    const char *s = get_optval_string(DATA, OPT_C);
    CompactMethod method = COMPACT_UNSET;

    if (s == NULL || *s == '\0') {
	*err = E_PARSE;
    } else if (!strcmp(s, "average")) {
	method = COMPACT_AVG;
    } else if (!strcmp(s, "sum")) {
	method = COMPACT_SUM;
    } else if (!strcmp(s, "first")) {
	method = COMPACT_SOP;
    } else if (!strcmp(s, "last")) {
	method = COMPACT_EOP;
    } else if (!strcmp(s, "spread")) {
	method = COMPACT_SPREAD;
    } else {
	gretl_errmsg_sprintf(_("field '%s' in command is invalid"), s);
	*err = E_PARSE;
    }

    return method;
}

/* 2-D array of doubles, allocated space in second
   position (as in a DATASET) */

static double **new_dbZ (int n)
{
    double **Z;
    int t;

    Z = malloc(2 * sizeof *Z);
    if (Z == NULL) return NULL;

    Z[0] = NULL;
    Z[1] = malloc(n * sizeof **Z);

    if (Z[1] == NULL) {
	free(Z);
	return NULL;
    }

    for (t=0; t<n; t++) {
	Z[1][t] = NADBL;
    }

    return Z;
}

static void free_dbZ (double **dbZ)
{
    if (dbZ != NULL) {
	free(dbZ[1]);
	free(dbZ);
    }
}

static int parse_odbc_format_chunk (char **ps, int i)
{
    const char *numchars = "0123456789";
    char *chunk = NULL;
    char *p = *ps;
    int n, err = 0;

    /* advance to '%' */
    while (*p && *p != '%') p++;
    if (*p == '\0') {
	return E_PARSE;
    }

    p++; /* move past '%' */

    /* zero padding? */
    if (*p == '0') {
	p++;
    }

    /* optional width? */
    n = strspn(p, numchars);
    if (n == 1) {
	p++;
    } else if (n > 0) {
	return E_PARSE;
    }

    /* optional dot plus precision? */
    if (*p == '.') {
	p++;
	n = strspn(p, numchars);
	if (n == 1) {
	    p++;
	} else {
	    return E_PARSE;
	}
    }

    /* now we should have a conversion character */
    if (*p == 'd') {
	gretl_odinfo.coltypes[i] = GRETL_TYPE_INT;
    } else if (*p == 's') {
	gretl_odinfo.coltypes[i] = GRETL_TYPE_STRING;
    } else if (*p == 'f' || *p == 'g') {
	gretl_odinfo.coltypes[i] = GRETL_TYPE_DOUBLE;
    } else if (*p == 'D') {
	*p = 's';
	gretl_odinfo.coltypes[i] = GRETL_TYPE_DATE;
    } else {
	return E_PARSE;
    }

    /* append any trailing fixed chars */
    p++;
    while (*p && *p != '%') p++;
    n = p - *ps;

    chunk = gretl_strndup(*ps, n);
    if (chunk == NULL) {
	err = E_ALLOC;
    } else {
	err = strings_array_add(&gretl_odinfo.fmts,
				&gretl_odinfo.obscols,
				chunk);
	free(chunk);
    }

    *ps = p;

#if 1
    fprintf(stderr, "set obs coltype[%d] = %d (%s), fmt='%s'\n", i,
	    gretl_odinfo.coltypes[i],
	    gretl_type_get_name(gretl_odinfo.coltypes[i]),
	    gretl_odinfo.fmts[i]);
#endif

    return err;
}

static int parse_odbc_format (char *fmt)
{
    char *s = fmt;
    int i, err = 0;

    for (i=0; i<ODBC_OBSCOLS && !err && *s; i++) {
	err = parse_odbc_format_chunk(&s, i);
    }

    if (!err && *s != '\0') {
	err = E_PARSE;
    }

    free(fmt);

    return err;
}

static char *odbc_get_query (const char *s, int *err)
{
    char *query = NULL;
    const char *p;

    if (*s == '"') {
	query = gretl_quoted_string_strdup(s, NULL);
    } else {
	p = get_string_by_name(s);
	if (p != NULL) {
	    query = gretl_strdup(p);
	} else {
	    query = gretl_strdup(s);
	}
    }

    if (query == NULL) {
	*err = E_ALLOC;
    } else if (*query == '\0') {
	gretl_errmsg_set(_("Expected an SQL query string"));
	*err = E_PARSE;
    }

    return query;
}

/* Grab the series name(s) out of an ODBC "data" command.  If the SQL
   query is marked by "query=" (which was not required in the original
   gretl ODBC setup) we're able to get multiple series names,
   otherwise we're restricted to one.
*/

static char **odbc_get_varnames (const char **line, int *err)
{
    char **vnames = NULL;
    char vname[VNAMELEN];
    const char *s = *line;
    int len, loop_ok = 0, nv = 0;

    if (strstr(s, "query=")) {
	/* we know where the SQL query starts */
	loop_ok = 1;
    }

    while (!*err) {
	*vname = '\0';
	*err = extract_varname(vname, s, &len);

	if (!*err && len == 0) {
	    gretl_errmsg_set(_("Expected a valid variable name"));
	    *err = E_PARSE;
	}

	if (!*err) {
	    *err = check_varname(vname);
	}

	if (!*err) {
	    *err = strings_array_add(&vnames, &nv, vname);
	}

	if (!*err) {
	    s += len;
	    s += strspn(s, " ");
	}

	if (!loop_ok || *s == '\0' || !strncmp(s, "obs-", 4) ||
	    !strncmp(s, "query=", 6)) {
	    /* got to the end of the varnames section */
	    break;
	}
    }

    if (*err) {
	strings_array_free(vnames, nv);
	vnames = NULL;
    } else {
	gretl_odinfo.nvars = nv;
    }

    *line = s;

    return vnames;
}

static double s_tab_get (int i, int t, series_table *stl, series_table *str)
{
    const double *x = gretl_odinfo.X[i];
    const char *sr;
    double ret = NADBL;

    /* get the string value for the imported obs */
    sr = series_table_get_string(str, x[t]);
    /* look up its index "on the left" */
    ret = series_table_get_value(stl, sr);
    if (na(ret)) {
	/* not found: so try adding it to the LHS table */
	series_table_add_string(stl, sr);
	ret = series_table_get_value(stl, sr);
    }

    return ret;
}

static int m2q (int m)
{
    if (m == 1) return 1;
    else if (m == 4) return 2;
    else if (m == 7) return 3;
    else if (m == 10) return 4;
    else return -1;
}

static int try_iso_8601 (const char *s, DATASET *dset)
{
    int t = -1;

    if (dataset_is_time_series(dset)) {
	char obsstr[OBSLEN] = {0};
	int y, m, d;

	if (sscanf(s, "%d-%d-%d", &y, &m, &d) == 3) {
	    if (dset->pd == 4 && d == 1) {
		sprintf(obsstr, "%04d:%d", y, m2q(m));
	    } else if (dset->pd == 12 && d == 1) {
		sprintf(obsstr, "%04d:%02d", y, m);
	    } else if (dset->pd == 1 && m == 1 && d == 1) {
		sprintf(obsstr, "%04d", y);
	    }
	    t = dateton(obsstr, dset);
	}
    }

    return t;
}

static int odbc_transcribe_data (char **vnames, DATASET *dset,
				 int vmin, int newvars,
				 gretlopt opt, PRN *prn)
{
    char label[MAXLABEL];
    int *gstlist = NULL;
    int *gstlnew = NULL;
    int nv = gretl_odinfo.nvars;
    int n = gretl_odinfo.nrows;
    int nrepl = nv - newvars;
    int simple_fill = (opt & OPT_F);
    int i, s, t, v;
    int spos = 1;
    int err = 0;

    if (gretl_odinfo.gst != NULL) {
	gstlist = string_table_copy_list(gretl_odinfo.gst);
	gstlnew = gretl_list_new(gstlist[0]);
	gstlnew[0] = 0;
    }

    for (i=0; i<nv && !err; i++) {
	series_table *str = NULL;
	series_table *stl = NULL;
	int vnew = 1; /* is this a new series? */
	int obs_used = 0;

	if (nrepl > 0) {
	    /* we're replacing some series */
	    v = current_series_index(dset, vnames[i]);
	} else {
	    /* all the series are new */
	    v = -1;
	}

	if (v < 0) {
	    /* a new series */
	    v = vmin++;
	    strcpy(dset->varname[v], vnames[i]);
	    sprintf(label, "ODBC series %d", i + 1);
	    series_set_label(dset, v, label);
	} else {
	    /* an existing series */
	    vnew = 0;
	    stl = series_get_string_table(dset, v);
	}

	if (in_gretl_list(gstlist, i+1)) {
	    /* the imported data are string-valued */
	    if (vnew) {
		gstlnew[spos++] = v;
		gstlnew[0] += 1;
	    } else if (stl == NULL) {
		gretl_errmsg_sprintf(_("%s: can't mix numeric and string data"),
				     dset->varname[v]);
		err = E_TYPES;
	    } else {
		str = gretl_string_table_detach_col(gretl_odinfo.gst, i+1);
	    }
	    if (!err && gretl_messages_on()) {
		pprintf(prn, _("%s: string-valued\n"), dset->varname[v]);
	    }
	} else if (stl != NULL) {
	    /* string-valued in dataset, numeric data from ODBC */
	    gretl_errmsg_sprintf(_("%s: can't mix numeric and string data"),
				 dset->varname[v]);
	    err = E_TYPES;
	}

	if (err) {
	    break;
	}

	if (gretl_odinfo.S != NULL) {
	    /* got obs identifiers via ODBC */
	    if (vnew) {
		for (t=0; t<dset->n; t++) {
		    dset->Z[v][t] = NADBL;
		}
	    }
	    for (s=0; s<n; s++) {
		t = dateton(gretl_odinfo.S[s], dset);
		if (t < 0) {
		    t = try_iso_8601(gretl_odinfo.S[s], dset);
		}
		if (t >= dset->t1 && t <= dset->t2) {
		    if (str != NULL) {
			dset->Z[v][t] = s_tab_get(i, s, stl, str);
		    } else {
			dset->Z[v][t] = gretl_odinfo.X[i][s];
		    }
		    obs_used++;
		} else {
		    fprintf(stderr, "Rejecting obs '%s'\n", gretl_odinfo.S[s]);
		}
	    }
	} else {
	    /* no obs identifiers via ODBC */
	    int ns = dset->t2 - dset->t1 + 1;

	    if (n == ns || simple_fill) {
		s = 0;
	    } else if (n == dset->n) {
		s = dset->t1;
	    } else {
		gretl_errmsg_sprintf(_("%s: don't know how to align the data!"),
				     dset->varname[v]);
		err = E_DATA;
	    }
	    for (t=0; t<dset->n && !err; t++) {
		if (t >= dset->t1 && t <= dset->t2 && s < n) {
		    if (str != NULL) {
			dset->Z[v][t] = s_tab_get(i, s++, stl, str);
		    } else {
			dset->Z[v][t] = gretl_odinfo.X[i][s++];
		    }
		    obs_used++;
		} else if (vnew) {
		    dset->Z[v][t] = NADBL;
		}
	    }
	}

	if (str != NULL) {
	    series_table_destroy(str);
	}

	if (!err && vnew && obs_used == 0) {
	    gretl_warnmsg_sprintf(_("ODBC import: '%s': no valid observations in sample range"),
				  vnames[i]);
	}
    }

    if (err) {
	dataset_drop_last_variables(dset, newvars);
	if (gretl_odinfo.gst != NULL) {
	    gretl_string_table_destroy(gretl_odinfo.gst);
	    gretl_odinfo.gst = NULL;
	}
    } else if (gretl_odinfo.gst != NULL) {
	if (gstlnew[0] == 0) {
	    /* no series tables to transfer */
	    gretl_string_table_destroy(gretl_odinfo.gst);
	    gretl_odinfo.gst = NULL;
	} else {
	    string_table_replace_list(gretl_odinfo.gst, gstlnew);
	    gstlnew = NULL; /* donated to table */
	    gretl_string_table_save(gretl_odinfo.gst, dset);
	}
    }

    free(gstlist);
    free(gstlnew);

    return err;
}

static int odbc_count_new_vars (char **vnames, int nv,
				const DATASET *dset)
{
    int newv = nv;

    if (dset->v > 0) {
	int i;

	for (i=0; i<nv; i++) {
	    if (current_series_index(dset, vnames[i]) > 0) {
		newv--;
	    }
	}
    }

    return newv;
}

/* data series [obs-format=format-string] [query=]query-string */

static int odbc_get_series (const char *line, DATASET *dset,
			    gretlopt opt, PRN *prn)
{
    int (*get_data) (ODBC_info *, gretlopt, PRN *);
    char **vnames = NULL;
    char *format = NULL;
    int err = 0;

    if (gretl_odinfo.dsn == NULL) {
	gretl_errmsg_set(_("No database has been opened"));
	return 1;
    } else if (dset->n == 0) {
	return E_NODATA;
    }

    /* get "series" field */
    vnames = odbc_get_varnames(&line, &err);
    if (err) {
	return err;
    }

    /* optional "obs-format" field */
    if (!strncmp(line, "obs-format=", 11)) {
	line += 11;
	format = gretl_quoted_string_strdup(line, (const char **) &line);
	if (format == NULL) {
	    err = E_PARSE;
	} else {
	    err = parse_odbc_format(format);
	}
    }

    /* now the query to pass to the database */
    if (!err) {
	line += strspn(line, " ");
	if (!strncmp(line, "query=", 6)) {
	    line += 6;
	}
	gretl_odinfo.query = odbc_get_query(line, &err);
    }

    if (!err) {
	if (opt & OPT_V) {
	    pprintf(prn, "SQL query: '%s'\n", gretl_odinfo.query);
	}
	gretl_error_clear();

	get_data = get_plugin_function("gretl_odbc_get_data");

	if (get_data == NULL) {
	    err = 1;
	} else {
	    err = (*get_data) (&gretl_odinfo, opt, prn);
	}
    }

    if (!err) {
	int n = gretl_odinfo.nrows;
	int nv = gretl_odinfo.nvars;
	int newvars, vmin = 1;

	if (gretl_messages_on()) {
	    pprintf(prn, _("Retrieved %d observations on %d series via ODBC\n"),
		    n, nv);
	}

	if (dset->v == 0) {
	    /* the data array is still empty */
	    newvars = nv;
	    dset->v = 1 + nv;
	    err = start_new_Z(dset, 0);
	} else {
	    newvars = odbc_count_new_vars(vnames, nv, dset);
	    vmin = dset->v;
	    if (newvars > 0) {
		err = dataset_add_series(dset, newvars);
	    }
	}

	if (!err) {
	    err = odbc_transcribe_data(vnames, dset, vmin, newvars, opt, prn);
	}
    }

    strings_array_free(vnames, gretl_odinfo.nvars);
    ODBC_info_clear_read();

    return err;
}

/* dbnomics function in separate file */

#include "dbnread.c"

/* called from loop in db_get_series() */

static int get_one_db_series (const char *sername,
			      const char *altname,
			      DATASET *dset,
			      CompactMethod cmethod,
			      const char *idxname,
			      PRN *prn)
{
    CompactMethod this_method = cmethod;
    const char *impname;
    SERIESINFO sinfo;
    double **dbZ;
    int v, err = 0;

    series_info_init(&sinfo);

    /* are we using a specified name for importation? */
    impname = (*altname == '\0')? sername : altname;

    /* see if the series is already in the dataset */
    v = series_index(dset, impname);
    if (v < dset->v && cmethod == COMPACT_UNSET) {
	this_method = series_get_compact_method(dset, v);
    }

#if DB_DEBUG
    fprintf(stderr, "get_one_db_series: dset->v=%d, v=%d, name='%s'\n",
	    dset->v, v, impname);
    fprintf(stderr, "this_var_method = %d\n", this_method);
#endif

    /* find the series information in the database */
    if (saved_db_type == GRETL_RATS_DB) {
	err = get_rats_series_info(sername, &sinfo);
    } else if (saved_db_type == GRETL_PCGIVE_DB) {
	err = get_pcgive_series_info(sername, &sinfo);
    } else {
	err = get_native_series_info(sername, &sinfo, idxname);
    }

    if (err) {
	fprintf(stderr, "get_one_db_series: failed to get series info\n");
	return err;
    }

    /* temporary data array */
    dbZ = new_dbZ(sinfo.nobs);
    if (dbZ == NULL) {
	gretl_errmsg_set(_("Out of memory!"));
	return E_ALLOC;
    }

#if DB_DEBUG
    fprintf(stderr, "get_one_db_series: offset=%d, nobs=%d\n",
	    sinfo.offset, sinfo.nobs);
#endif

    if (saved_db_type == GRETL_RATS_DB) {
	err = get_rats_db_data(saved_db_name, &sinfo, dbZ);
    } else if (saved_db_type == GRETL_PCGIVE_DB) {
	err = get_pcgive_db_data(saved_db_name, &sinfo, dbZ);
#ifdef USE_CURL
    } else if (saved_db_type == GRETL_NATIVE_DB_WWW) {
	err = get_remote_db_data(saved_db_name, &sinfo, dbZ);
#endif
    } else {
	err = get_native_db_data(saved_db_name, &sinfo, dbZ);
    }

#if DB_DEBUG
    fprintf(stderr, "get_one_db_series: get_db_data gave %d\n", err);
#endif

    if (err == DB_MISSING_DATA) {
	fprintf(stderr, "There were missing data\n");
	err = 0;
    }

#if DB_DEBUG
    fprintf(stderr, "sinfo.varname='%s', this_method=%d\n",
	    sinfo.varname, this_method);
#endif

    if (!err) {
	if (*altname != '\0') {
	    /* switch the recorded name now */
	    strcpy(sinfo.varname, altname);
	}
	if (this_method == COMPACT_SPREAD) {
	    err = lib_spread_db_data(dbZ, &sinfo, dset, prn);
	} else {
	    err = lib_add_db_data(dbZ, &sinfo, dset, this_method,
				  v, prn);
	}
    }

    series_info_clear(&sinfo);
    free_dbZ(dbZ);

    return err;
}

static int is_glob (const char *s)
{
    return strchr(s, '*') || strchr(s, '?');
}

static int process_import_name_option (char *vname)
{
    const char *s = get_optval_string(DATA, OPT_N);
    int err = 0;

    if (s == NULL) {
	err = E_DATA;
    } else {
	err = check_varname(s);
    }

    if (!err) {
	strcpy(vname, s);
    }

    return err;
}

/* Remedial code to handle dicontinued FRED series whose
   names end with "96" and which are no longer included in
   gretl's fedstl database, but are continued in series
   whose names (mostly) end in "1".
*/

static void maybe_update_name (char *argname,
                               char *altname)
{
    const char *dups[] = {
        "gdpc", "fpic", "fgcec", "gcec", "slcec9", "gpdic",
        "prfic", "expgsc", "impgsc96", "cbic", "finslc", NULL
    };
    char test[16];
    int i, done = 0;

    for (i=0; dups[i] != NULL; i++) {
        sprintf(test, "%s96", dups[i]);
        if (!strcmp(test, argname)) {
            if (*altname == '\0') {
                strcpy(altname, argname);
            }
            sprintf(argname, "%s1", dups[i]);
            done = 1;
            break;
        }
    }

    if (!done && !strcmp(argname, "netexc96")) {
        /* special case: netexc96 -> netexc */
        if (*altname == '\0') {
            strcpy(altname, argname);
        }
        strcpy(argname, "netexc");
    }
}

/* main function for getting one or more series out of a
   database (including ODBC) via command-line/script
*/

int db_get_series (const char *line, DATASET *dset,
		   gretlopt opt, PRN *prn)
{
    char altname[VNAMELEN] = {0};
    char **vnames = NULL;
    char *idxname = NULL;
    CompactMethod cmethod;
    int i, nnames = 0;
    int from_scratch = 0;
    int FRED = 0;
    int err = 0;

    if (opt & OPT_O) {
	return odbc_get_series(line, dset, opt, prn);
    }

    if (opt & OPT_N) {
	/* --name=whatever to rename a db import */
	err = process_import_name_option(altname);
	if (err) {
	    return err;
	}
    }

#if DB_DEBUG
    fprintf(stderr, "db_get_series: line='%s', dset=%p\n",
	    line, (void *) dset);
    fprintf(stderr, "db_name = '%s'\n", saved_db_name);
#endif

    if (*saved_db_name == '\0') {
	gretl_errmsg_set(_("No database has been opened"));
	return 1;
    }

    from_scratch = (dset->n == 0);
    FRED = strstr(saved_db_name, "fedstl") != NULL;

    if (opt & OPT_C) {
	/* new-style: compaction method supplied as option */
	cmethod = compact_method_from_option(&err);
    } else {
	/* legacy */
	line = get_compact_method_and_advance(line, &cmethod);
    }

    if (!err) {
	if (string_is_blank(line)) {
	    err = E_DATA;
	} else {
	    /* get the variable names on the line */
	    vnames = gretl_string_split(line, &nnames, NULL);
	    if (vnames == NULL) {
		err = E_ALLOC;
	    }
	}
    }

    if (!err && nnames > 1 && *altname != '\0') {
	/* altname only works for a single series? */
	err = E_BADOPT;
    }

    if (!err) {
	if (saved_db_type == GRETL_NATIVE_DB) {
	    idxname = native_db_index_name();
	} else if (saved_db_type == GRETL_NATIVE_DB_WWW) {
#ifdef USE_CURL
	    idxname = g_strdup_printf("%sdbtmp.idx", gretl_dotdir());
	    err = remote_db_index_to_file(idxname);
#endif
	}
    }

    /* now process the imports individually */

    for (i=0; i<nnames && !err; i++) {
	if (is_glob(vnames[i])) {
	    /* globbing works only for native databases */
	    if (*altname != '\0') {
		/* can't do it */
		err = E_BADOPT;
	    } else if (saved_db_type == GRETL_NATIVE_DB ||
		       saved_db_type == GRETL_NATIVE_DB_WWW) {
		char **tmp;
		int j, nmatch;

		tmp = native_db_match_series(vnames[i], &nmatch,
					     idxname, &err);
		for (j=0; j<nmatch && !err; j++) {
		    err = get_one_db_series(tmp[j], altname, dset,
					    cmethod, idxname, prn);
		}
		strings_array_free(tmp, nmatch);
	    } else {
		err = E_INVARG;
	    }
	} else if (saved_db_type == GRETL_DBNOMICS) {
	    err = get_one_dbnomics_series(vnames[i], altname, dset,
					  cmethod, prn);
	} else {
            if (FRED && strstr(vnames[i], "96")) {
                maybe_update_name(vnames[i], altname);
            }
	    err = get_one_db_series(vnames[i], altname, dset,
				    cmethod, idxname, prn);
            *altname = '\0'; /* zero after use */
	}
    }

    strings_array_free(vnames, nnames);

    if (!err && !(opt & OPT_Q) && gretl_messages_on()) {
	pprintf(prn, _("Series imported OK"));
	pputc(prn, '\n');
	if (from_scratch) {
	    print_smpl(dset, 0, OPT_NONE, prn);
	}
    }

    if (idxname != NULL) {
	if (saved_db_type == GRETL_NATIVE_DB_WWW) {
	    /* this file is a temporary download */
	    gretl_remove(idxname);
	}
	free(idxname);
    }

    return err;
}

static FILE *tempfile_open (char *fname, int *err)
{
    FILE *fp;

    strcat(fname, ".XXXXXX");
    fp = gretl_mktemp(fname, "w+");
    if (fp == NULL) {
	*err = E_FOPEN;
    }

    return fp;
}

static void maybe_fclose (FILE *fp)
{
    if (fp != NULL) {
	fclose(fp);
    }
}

#define DBUFLEN 1024

static int db_delete_series (const char *line, const int *list,
			     const char *fname, PRN *prn)
{
    dbnumber buf[DBUFLEN];
    char src1[FILENAME_MAX];
    char src2[FILENAME_MAX];
    char tmp1[FILENAME_MAX];
    char tmp2[FILENAME_MAX];
    char series[VNAMELEN];
    char *p, s[512];
    char **snames = NULL;
    FILE *fidx = NULL, *fbin = NULL;
    FILE *f1 = NULL, *f2 = NULL;
    int i, j, k, print, n, ns;
    int ndel = 0;
    int err = 0;

    if (fname == NULL) {
	if (*saved_db_name == '\0') {
	    gretl_errmsg_set(_("No database has been opened"));
	    err = 1;
	} else if (saved_db_type != GRETL_NATIVE_DB) {
	    gretl_errmsg_set(_("This only works for gretl databases"));
	    err = 1;
	} else {
	    err = open_native_db_files(saved_db_name, &fidx, src1, &fbin, src2);
	}
    } else {
	err = open_native_db_files(fname, &fidx, src1, &fbin, src2);
    }

    if (err) {
	return err;
    }

    strcpy(tmp1, gretl_dotdir());
    strcat(tmp1, "tmpidx");
    f1 = tempfile_open(tmp1, &err);
    if (err) {
	goto bailout;
    }

    strcpy(tmp2, gretl_dotdir());
    strcat(tmp2, "tmpbin");
    f2 = tempfile_open(tmp2, &err);
    if (err) {
	goto bailout;
    }

    if (line != NULL) {
	/* extract the variable names given on the line */
	ns = 0;
	while ((line = get_word_and_advance(line, series, VNAMELEN-1))
	       && !err) {
	    err = strings_array_add(&snames, &ns, series);
	}
	if (!err && ns == 0) {
	    fprintf(stderr, "Found no series names\n");
	    err = E_PARSE;
	}
    }

    print = k = 1;
    i = j = 0;

    while (fgets(s, sizeof s, fidx) && !err) {
	if (i == 0) {
	    /* always reprint the header */
	    fputs(s, f1);
	    i++;
	    continue;
	}

	if (i % 2 != 0) {
	    /* odd lines contain varnames */
	    print = 1;
	    if (snames != NULL) {
		sscanf(s, "%s", series);
		for (j=0; j<ns; j++) {
		    if (!strcmp(series, snames[j])) {
			print = 0;
			ndel++;
			break;
		    }
		}
	    } else {
		if (k <= list[0] && list[k] == j) {
		    k++;
		    print = 0;
		    ndel++;
		}
		j++;
	    }
	    if (print) {
		fputs(s, f1);
	    }
	} else {
	    /* even lines have obs information */
	    p = strstr(s, "n = ");
	    if (p != NULL) {
		sscanf(p + 4, "%d", &n);
	    } else {
		err = E_DATA;
		fprintf(stderr, "couldn't find obs for series\n");
	    }

	    if (!print) {
		fseek(fbin, n * sizeof(dbnumber), SEEK_CUR);
	    } else {
		int get, got, rem = n;

		fputs(s, f1);

		while (rem > 0 && !err) {
		    get = (rem > DBUFLEN)? DBUFLEN : rem;
		    got = fread(buf, sizeof(dbnumber), get, fbin);
		    if (got != get) {
			fprintf(stderr, "error reading binary data\n");
			err = E_DATA;
		    } else {
			fwrite(buf, sizeof(dbnumber), got, f2);
			rem -= got;
		    }
		}
	    }
	}
	i++;
    }

    if (snames != NULL) {
	strings_array_free(snames, ns);
    }

 bailout:

    maybe_fclose(fidx);
    maybe_fclose(fbin);
    maybe_fclose(f1);
    maybe_fclose(f2);

    if (!err && ndel > 0) {
	err = gretl_rename(tmp1, src1);
	if (!err) {
	    err = gretl_rename(tmp2, src2);
	}
    } else {
	gretl_remove(tmp1);
	gretl_remove(tmp2);
    }

    if (!err && prn != NULL) {
	pprintf(prn, _("Deleted %d series from %s\n"), ndel, src2);
    }

    return err;
}

int db_delete_series_by_name (const char *line, PRN *prn)
{
    return db_delete_series(line, NULL, NULL, prn);
}

int db_delete_series_by_number (const int *list, const char *fname)
{
    return db_delete_series(NULL, list, fname, NULL);
}

static void obs_to_ymd (const char *obs, int pd, int *y, int *m, int *d)
{
    *y = atoi(obs);
    *d = 1;

    if (pd == 12) {
	*m = atoi(obs + 5);
    } else if (pd == 4) {
	int q = atoi(obs + 5);

	*m = q * 3 - 2;
    } else {
	*m = 1;
    }
}

int db_range_check (int db_pd,
		    const char *db_stobs,
		    const char *db_endobs,
		    const char *varname,
		    DATASET *dset)
{
    double sd0_orig, sdn_orig, sd0, sdn;
    int err = 0;

    sd0 = get_date_x(db_pd, db_stobs);
    sdn = get_date_x(db_pd, db_endobs);

    if (db_pd >= 5 && db_pd <= 7 && !dated_daily_data(dset)) {
	/* convert 'orig' info to daily dates */
	int y, m, d;

	obs_to_ymd(dset->stobs, dset->pd, &y, &m, &d);
	sd0_orig = epoch_day_from_ymd(y, m, d);
	obs_to_ymd(dset->endobs, dset->pd, &y, &m, &d);
	sdn_orig = epoch_day_from_ymd(y, m, d);
    } else {
	sd0_orig = dset->sd0;
	sdn_orig = get_date_x(dset->pd, dset->endobs);
    }

    if (sd0 > sdn_orig || sdn < sd0_orig) {
	gretl_errmsg_sprintf(_("%s: observation range does not overlap\n"
			       "with the working data set"),
			     varname);
	err = 1;
    }

    return err;
}

int check_db_import_conversion (int pd, DATASET *dset)
{
    int target = dset->pd;
    int err = 0;

    if (pd == target) {
	; /* no conversion needed */
    } else if (pd == 1 && target == 4) {
	; /* annual to quarterly expansion */
    } else if (pd == 1 && target == 12) {
	; /* annual to monthly expansion */
    } else if (pd == 4 && target == 12) {
	; /* quarterly to monthly expansion */
    } else if (pd == 12 && target == 1) {
	; /* monthly to annual compaction */
    } else if (pd == 4 && target == 1) {
	; /* quarterly to annual compaction */
    } else if (pd == 12 && target == 4) {
	; /* monthly to quarterly compaction */
    } else {
	fprintf(stderr, "db import fail: pd = %d, target %d\n", pd, target);
	err = E_DATA;
    }

    return err;
}

static int check_db_import_full (int pd,
				 const char *stobs,
				 const char *endobs,
				 const char *varname,
				 DATASET *dset)
{
    int err = check_db_import_conversion(pd, dset);

    if (err) {
	gretl_errmsg_sprintf(_("%s: can't handle conversion"),
			     varname);
    } else {
	err = db_range_check(pd, stobs, endobs, varname, dset);
    }

#if DB_DEBUG
    if (err) {
	fprintf(stderr, "check_db_import_full: err = %d\n", err);
	fprintf(stderr, "(dset->n = %d)\n", dset->n);
    }
#endif

    return err;
}

/* We'll do "spread" compaction for monthly to quarterly or annual,
   quarterly to annual, or daily to monthly or quarterly. Other
   cases are rejected.
*/

static int compact_spread_pd_check (int high, int low)
{
    if ((low == 12 || low == 4) &&
	(high == 5 || high == 6 || high == 7)) {
	/* daily to monthly or quarterly */
	return 0;
    }

    if (!(high == 12 && low == 1) &&
	!(high == 12 && low == 4) &&
	!(high == 4 && low == 1)) {
	gretl_errmsg_set(_("Unsupported conversion"));
	return E_DATA;
    }

    return 0;
}

static void
init_datainfo_from_sinfo (DATASET *dset, SERIESINFO *sinfo)
{
    dset->pd = sinfo->pd;

    strcpy(dset->stobs, sinfo->stobs);
    strcpy(dset->endobs, sinfo->endobs);
    colonize_obs(dset->stobs);
    colonize_obs(dset->endobs);

    dset->sd0 = get_date_x(dset->pd, dset->stobs);
    dset->n = sinfo->nobs;
    dset->v = 2;

    dset->t1 = 0;
    dset->t2 = dset->n - 1;
}

/* construct a little dataset as a temporary wrapper for an
   import using compact=spread
*/

static DATASET *make_import_tmpset (const DATASET *dset,
				    SERIESINFO *sinfo,
				    double **dbZ,
				    int *err)
{
    DATASET *tmpset = NULL;

    *err = compact_spread_pd_check(sinfo->pd, dset->pd);
    if (*err) {
	return NULL;
    }

    tmpset = datainfo_new();
    if (tmpset == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    tmpset->v = 2;
    tmpset->n = sinfo->nobs;

    tmpset->Z = malloc(2 * sizeof *tmpset->Z);
    if (tmpset->Z == NULL) {
	*err = E_ALLOC;
	free(tmpset);
	return NULL;
    }

    *err = dataset_allocate_varnames(tmpset);
    if (*err) {
	free(tmpset->Z[1]);
	free(tmpset->Z);
	free(tmpset);
	return NULL;
    }

    tmpset->Z[0] = NULL;
    tmpset->Z[1] = dbZ[1];
    dbZ[1] = NULL; /* note: stolen! */

    tmpset->t1 = sinfo->t1;
    tmpset->t2 = sinfo->t2;
    tmpset->pd = sinfo->pd;
    strcpy(tmpset->stobs, sinfo->stobs);
    strcpy(tmpset->endobs, sinfo->endobs);
    tmpset->structure = TIME_SERIES;
    tmpset->sd0 = get_date_x(tmpset->pd, tmpset->stobs);

    strcpy(tmpset->varname[1], sinfo->varname);

#if 0
    PRN *prn = gretl_print_new(GRETL_PRINT_STDERR, NULL);
    fprintf(stderr, "import_tmpset: t1=%d, t2=%d, nobs=%d, pd=%d, offset=%d\n",
	    sinfo->t1, sinfo->t2, sinfo->nobs, sinfo->pd, sinfo->offset);
    printdata(NULL, NULL, tmpset, OPT_O, prn);
    gretl_print_destroy(prn);
#endif

    return tmpset;
}

static int
real_transcribe_db_data (const char *stobs, int nobs,
			 const DATASET *dset, int dbv,
			 const double *xvec)
{
    int t, pad1, pad2;
    int start, stop;
    double x;

    pad1 = dateton(stobs, dset);
    pad2 = dset->n - nobs - pad1;

    if (pad1 > 0) {
	fprintf(stderr, "Padding at start, %d obs\n", pad1);
	for (t=0; t<pad1; t++) {
	    dset->Z[dbv][t] = NADBL;
	}
	start = pad1;
    } else {
	start = 0;
    }
    if (pad2 > 0) {
	int n = dset->n;

	fprintf(stderr, "Padding at end, %d obs\n", pad2);
	for (t=n-1; t>=n-1-pad2; t--) {
	    dset->Z[dbv][t] = NADBL;
	}
	stop = n - pad2;
    } else {
	stop = dset->n;
    }

    fprintf(stderr, "Filling in values from %d to %d\n", start, stop - 1);
    for (t=start; t<stop; t++) {
	x = xvec[t - pad1];
	dset->Z[dbv][t] = (x == DBNA)? NADBL : x;
    }

    return 0;
}

int transcribe_db_data (DATASET *dset, int targv,
			const double *src, int pd,
			int nobs, char *stobs,
			CompactMethod cmethod)
{
    double *xvec = (double *) src;
    int free_xvec = 0;

    if (pd != dset->pd) {
	if (pd < dset->pd) {
	    /* the series needs to be expanded */
	    xvec = expand_db_series(src, pd, &nobs, stobs, dset);
	} else {
	    /* the series needs to be compacted */
	    xvec = compact_db_series(src, pd, &nobs, stobs, dset->pd,
				     cmethod);
	}
	if (xvec == NULL) {
	    return E_ALLOC;
	}
	free_xvec = 1;
    }

    real_transcribe_db_data(stobs, nobs, dset, targv, xvec);

    if (free_xvec) {
	free(xvec);
    }

    return 0;
}

/* Processes a single db series in "spread" mode, meaning
   that multiple series are added to the target dataset,
   @dset. This variant is associated with gretl databases.
*/

int lib_spread_db_data (double **dbZ, SERIESINFO *sinfo,
			DATASET *dset, PRN *prn)
{
    int err = 0;

    if (dset == NULL || dset->v == 0) {
	gretl_errmsg_set(_("\"compact=spread\": requires a dataset in place"));
	err = E_DATA;
    } else {
	DATASET *tmpset = make_import_tmpset(dset, sinfo, dbZ, &err);

	if (!err) {
	    err = do_compact_spread(tmpset, dset->pd);
	}
	if (!err) {
	    err = merge_or_replace_data(dset, &tmpset, OPT_X | OPT_U, prn);
	}
    }

    return err;
}

/* Processes a single db series in "spread" mode, meaning
   that multiple series are added to the target dataset,
   @dset. This variant is associated with dbnomics import.
*/

int lib_spread_dbnomics_data (DATASET *dset, DATASET *dbset,
			      PRN *prn)
{
    int err = 0;

    if (dset == NULL || dset->v == 0) {
	gretl_errmsg_set(_("\"compact=spread\": requires a dataset in place"));
	err = E_DATA;
    } else {
	err = do_compact_spread(dbset, dset->pd);
	if (!err) {
	    /* we add OPT_K ("keep") to prevent destruction of @dbset:
	       we're bypassing get_merge_opts(), so we'd better know
	       what we're doing!
	    */
	    gretlopt merge_opt = (OPT_X | OPT_U | OPT_K);

	    err = merge_or_replace_data(dset, &dbset, merge_opt, prn);
	}
    }

    return err;
}

/* Processes a single db series, adding it to @dset if
   possible (perhaps after compaction or expansion).
*/

static int lib_add_db_data (double **dbZ, SERIESINFO *sinfo,
			    DATASET *dset, CompactMethod cmethod,
			    int dbv, PRN *prn)
{
    int new = (dbv == dset->v);
    int err = 0;

    if (sinfo == NULL && dbZ == NULL) {
	fprintf(stderr, "lib_add_db_data: broken call!\n");
	return E_DATA;
    }

    if (cmethod == COMPACT_UNSET) {
	/* impose default if need be */
	cmethod = COMPACT_AVG;
    }

    if (dset->n == 0) {
	/* if the existing dataset is empty, initialize it
	   using info from the database series
	*/
	init_datainfo_from_sinfo(dset, sinfo);
	dset->v = 0; /* trigger for creating data array below */
	if (dset->pd != 1 || strcmp(dset->stobs, "1")) {
	    dset->structure = TIME_SERIES;
	}
    } else {
	err = check_db_import_full(sinfo->pd, sinfo->stobs, sinfo->endobs,
				   sinfo->varname, dset);
	if (err) {
	    return err;
	}
    }

    if (dset->v == 0) {
	/* the data array is still empty */
	dset->v = 2;
	dbv = 1;
	if (start_new_Z(dset, 0)) {
	    return E_ALLOC;
	}
    } else if (new && dataset_add_series(dset, 1)) {
	return E_ALLOC;
    }

#if DB_DEBUG
    fprintf(stderr, "dset->Z=%p\n", (void *) dset->Z);
    fprintf(stderr, "dset->n = %d, dset->v = %d, dbv = %d\n",
	    dset->n, dset->v, dbv);
#endif

    err = transcribe_db_data(dset, dbv, dbZ[1], sinfo->pd, sinfo->nobs,
			     sinfo->stobs, cmethod);

    if (!err) {
	/* common stuff for adding a var */
	strcpy(dset->varname[dbv], sinfo->varname);
	series_set_label(dset, dbv, sinfo->descrip);
	series_set_compact_method(dset, dbv, cmethod);
	if (sinfo->pd < dset->pd) {
	    series_set_orig_pd(dset, dbv, sinfo->pd);
	}
    } else if (new) {
	/* we added a series that has not been filled */
	dataset_drop_last_variables(dset, 1);
    }

    return err;
}

/* Compact an individual series, in the context of converting an
   entire working dataset to a lower frequency: used in all cases
   except conversion from daily to monthly.
*/

static double *compact_series (const DATASET *dset, int i, int oldn,
			       int startskip, int min_startskip,
			       int compfac, CompactMethod method)
{
    const double *src = dset->Z[i];
    double *x;
    int lead = startskip - min_startskip;
    int to_weekly = (compfac >= 5 && compfac <= 7);
    int t, idx;

#if DB_DEBUG
    fprintf(stderr, "compact_series: startskip=%d, min_startskip=%d, compfac=%d "
	    "lead=%d\n", startskip, min_startskip, compfac, lead);
#endif

    x = malloc(dset->n * sizeof *x);
    if (x == NULL) {
	return NULL;
    }

    for (t=0; t<dset->n; t++) {
	x[t] = NADBL;
    }

    idx = startskip;

    for (t=lead; t<dset->n && idx<oldn; t++) {
	if (method == COMPACT_SOP) {
	    if (to_weekly && na(src[idx]) && idx < oldn - 1) {
		/* allow one day's slack */
		x[t] = src[idx + 1];
	    } else {
		x[t] = src[idx];
	    }
	} else if (method == COMPACT_EOP) {
	    if (to_weekly && na(src[idx]) && idx > startskip) {
		/* one day's slack */
		x[t] = src[idx - 1];
	    } else {
		x[t] = src[idx];
	    }
	} else if (method == COMPACT_SUM || method == COMPACT_AVG) {
	    int j, st, den = compfac;
	    int n_ok = 0;

	    if (idx + compfac - 1 > oldn - 1) {
		break;
	    }

	    x[t] = 0.0;

	    for (j=0; j<compfac; j++) {
		st = idx + j;
		if (na(src[st])) {
		    if (to_weekly) {
			den--;
		    } else {
			x[t] = NADBL;
			break;
		    }
		} else {
		    /* got a valid observation */
		    n_ok++;
		    x[t] += src[st];
		}
	    }

	    if (n_ok == 0) {
		x[t] = NADBL;
	    }

	    if (method == COMPACT_AVG && !na(x[t])) {
		if (den > 0) {
		    x[t] /= den;
		} else {
		    x[t] = NADBL;
		}
	    }
	}
	idx += compfac;
    }

    return x;
}

/* Determine year and period (either month or quarter,
   depending on the value of @pd) for observation @t in
   daily dataset @dset.
*/

static int daily_yp (const DATASET *dset, int t,
		     int pd, int *y, int *p)
{
    char obs[12];
    int mon, day;

    ntolabel(obs, t, dset);

    if (sscanf(obs, YMD_READ_FMT, y, &mon, &day) != 3) {
	return E_DATA;
    }

    if (pd == 12) {
	*p = mon;
    } else {
	/* convert month to quarter */
	*p = 1 + (mon - 1) / 3;
    }

    return 0;
}

#define DAYDBG 0

/* For a single row, @cset_t, of a compacted dataset,
   write daily values into the set of monthly or
   quarterly series that will represent them. The
   daily data are drawn from @dset and transcribed to
   @cset.
*/

static void fill_cset_t (const DATASET *dset,
			 int *startday,
			 DATASET *cset,
			 int cset_t,
			 int compfac,
			 int qmonth)
{
    char obs[OBSLEN];
    double cvec[30];
    int idx[30];
    const double *z;
    int y, p, pstart = 0;
    int effn, ndays = 0;
    int i, j, k, s, t, t0;
    double zsum = 0.0;

    t0 = *startday;
    y = p = 0;

    /* how many daily obs do we have in this month? */
    for (t=t0; t<dset->n; t++) {
	daily_yp(dset, t, 12, &y, &p);
	if (t == t0) {
	    pstart = p;
	} else if (p != pstart) {
	    break;
	}
	ndays++;
    }

#if 0
    fprintf(stderr, "fill_cset_t: ndays = %d, compfac = %d\n",
	    ndays, compfac);
#endif

    /* construct array of month-day indices */
    for (j=0; j<compfac && j<ndays; j++) {
	ntolabel(obs, t0 + j, dset);
	idx[j] = date_to_daily_index(obs, dset->pd);
    }

    /* the outer loop is over the daily series in the
       source dataset */

    k = 1 + qmonth * compfac;

    for (i=1; i<dset->v; i++) {
	z = dset->Z[i] + t0;
	for (j=0; j<compfac; j++) {
	    cvec[j] = NADBL;
	}
	effn = 0;
	zsum = 0.0;
	for (j=0; j<compfac && j<ndays; j++) {
	    s = idx[j];
	    cvec[s] = z[j];
	    if (!na(cvec[s])) {
		zsum += cvec[s];
		effn++;
	    }
	}
	if (effn < compfac) {
	    /* we have some padding to do */
	    double zbar = zsum / effn;

	    for (j=0; j<compfac; j++) {
		if (na(cvec[j])) {
		    cvec[j] = zbar;
		}
	    }
	}
	/* transcribe into target dataset */
	for (j=0; j<compfac; j++) {
	    cset->Z[k+j][cset_t] = cvec[j];
	}

	k += compfac;
    }

    *startday += ndays;
}

#define SPREAD_DEBUG 0

/* compact daily data to monthly or quarterly using the
   "spread" method */

static DATASET *compact_daily_spread (const DATASET *dset,
				      int newpd,
				      int *nv,
				      int *err)
{
    const char *periods[] = {
	"month",
	"quarter"
    };
    const char *period;
    DATASET *cset = NULL;
    char label[MAXLABEL];
    int oldpd = dset->pd;
    int compfac;
    int v, i, j, k, t, T;
    int startyr, startper = 0;
    int endyr, endper = 0;
    int startday;

#if SPREAD_DEBUG
    fprintf(stderr, "*** compact_daily_spread (newpd=%d) ***\n", newpd);
#endif

    daily_yp(dset, 0, newpd, &startyr, &startper);
    daily_yp(dset, dset->n - 1, newpd, &endyr, &endper);
    compfac = midas_days_per_period(dset->pd, newpd);

    if (newpd == 12) {
	period = periods[0];
    } else if (newpd == 4) {
	period = periods[1];
    } else {
	*err = E_DATA;
	return NULL;
    }

    T = newpd * (endyr - startyr) + (endper - startper + 1);

    if (T <= 1) {
	*err = E_DATA;
	return NULL;
    }

    /* the number of series, after compaction */
    v = 1 + (dset->v - 1) * compfac;

#if SPREAD_DEBUG
    fprintf(stderr, "oldpd %d, newpd %d, nvars=%d, T=%d, start=%d:%d, end=%d:%d\n",
	    dset->pd, newpd, v, T, startyr, startper, endyr, endper);
#endif

    cset = create_new_dataset(v, T, 0);
    if (cset == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    if (newpd == 12) {
	sprintf(cset->stobs, "%d:%02d", startyr, startper);
	sprintf(cset->endobs, "%d:%02d", endyr, endper);
    } else {
	sprintf(cset->stobs, "%d:%d", startyr, startper);
	sprintf(cset->endobs, "%d:%d", endyr, endper);
    }

    cset->pd = newpd;
    cset->structure = TIME_SERIES;
    cset->sd0 = get_date_x(cset->pd, cset->stobs);

    /* ensure no uninitialized data */
    for (i=1; i<v; i++) {
	for (t=0; t<T; t++) {
	    cset->Z[i][t] = NADBL;
	}
    }

    /* do the actual data transcription first */
    startday = 0;
    for (t=0; t<T; t++) {
	if (newpd == 4) {
	    fill_cset_t(dset, &startday, cset, t, compfac/3, 0);
	    fill_cset_t(dset, &startday, cset, t, compfac/3, 1);
	    fill_cset_t(dset, &startday, cset, t, compfac/3, 2);
	} else {
	    fill_cset_t(dset, &startday, cset, t, compfac, 0);
	}
    }

    /* then name the series and reorganize */

    k = 1;
    for (i=1; i<dset->v; i++) {
	double *xtmp;
	char sfx[16];
	int p;

	/* switch data order */
	for (j=0; j<compfac/2; j++) {
	    p = k + compfac - j - 1;
	    xtmp = cset->Z[k+j];
	    cset->Z[k+j] = cset->Z[p];
	    cset->Z[p] = xtmp;
	}

	/* names and labels */
	for (j=0; j<compfac; j++) {
	    strcpy(cset->varname[k+j], dset->varname[i]);
	    gretl_trunc(cset->varname[k+j], VNAMELEN - 5);
	    sprintf(sfx, "_d%02d", compfac - j);
	    strcat(cset->varname[k+j], sfx);
	    sprintf(label, "%s in day %d of %s", dset->varname[i],
		    compfac - j, period);
	    series_record_label(cset, k+j, label);
	    series_set_midas_period(cset, k+j, compfac - j);
	    series_set_midas_freq(cset, k+j, oldpd);
	    if (j == 0) {
		series_set_midas_anchor(cset, k+j);
	    }
	}

	/* advance column write position for next source series */
	k += compfac;
    }

#if SPREAD_DEBUG > 1
    PRN *prn = gretl_print_new(GRETL_PRINT_STDERR, NULL);
    printdata(NULL, NULL, cset, OPT_O, prn);
    gretl_print_destroy(prn);
#endif

    return cset;
}

/* compact an entire dataset, transcribing from each higher-frequency
   series to a set of lower-frequency series, each of which holds the
   observations from a given sub-period
*/

static DATASET *compact_data_spread (const DATASET *dset, int newpd,
				     int startmaj, int startmin,
				     int endmaj, int endmin,
				     int *nv, int *err)
{
    const char *subper[] = {
	"month",
	"quarter"
    };
    const char *period[] = {
	"year",
	"quarter"
    };
    const char *p0, *p1;
    DATASET *cset = NULL;
    char sfx[16];
    char label[MAXLABEL];
    int oldpd = dset->pd;
    int compfac = oldpd / newpd;
    int v, i, j, k, t, T;
    int q0 = 0, qT = 0;

    /* calculate @T, the number of observations that the compacted
       dataset should comprise
    */
    if (newpd == 1) {
	T = endmaj - startmaj + 1;
    } else if (newpd == 4) {
	T = oldpd * (endmaj - startmaj + 1) / compfac;
	q0 = 1 + (startmin - 1) / 3;
	qT = 1 + (endmin - 1) / 3;
	T += qT - q0 - 3;
    } else {
	*err = E_DATA;
	return NULL;
    }

    if (T <= 1) {
	*err = E_DATA;
	return NULL;
    }

    /* calculate @v, the number of series after compaction */
    v = 1 + (dset->v - 1) * compfac;

#if SPREAD_DEBUG
    fprintf(stderr, "oldpd %d, newpd %d, v=%d, T=%d, start=%d:%d, end=%d:%d\n",
	    oldpd, newpd, v, T, startmaj, startmin, endmaj, endmin);
#endif

    cset = create_new_dataset(v, T, 0);
    if (cset == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    if (newpd == 1) {
	sprintf(cset->stobs, "%d", startmaj);
	sprintf(cset->endobs, "%d", endmaj);
	p1 = period[0];
    } else {
	/* newpd must be 4 */
	sprintf(cset->stobs, "%d:%d", startmaj, q0);
	sprintf(cset->endobs, "%d:%d", endmaj, qT);
	p1 = period[1];
    }

    p0 = (oldpd == 12)? subper[0] : subper[1];

    cset->pd = newpd;
    cset->structure = TIME_SERIES;
    cset->sd0 = get_date_x(cset->pd, cset->stobs);

#if SPREAD_DEBUG
    fprintf(stderr, "stobs '%s', endobs '%s', sd0=%g, q0=%d\n",
	    cset->stobs, cset->endobs, cset->sd0, q0);
#endif

    k = 1; /* the first new series */

    for (i=1; i<dset->v; i++) {
	/* loop across original data series */
	double *xtmp;
	int offset;
	int p, s = 0;

	/* how many initial observations should be set to NA? */
	if (newpd == 1) {
	    offset = startmin - 1;
	} else {
	    offset = startmin - (1 + (q0 - 1) * compfac);
	}

	for (t=0; t<T; t++) {
	    /* loop across new time periods */
	    for (j=0; j<compfac; j++) {
		/* loop across new series <- sub-periods */
		while (s < offset) {
		    cset->Z[k+j][t] = NADBL;
		    offset--;
		    j++;
		}
		if (s < dset->n) {
		    cset->Z[k+j][t] = dset->Z[i][s];
		} else {
		    cset->Z[k+j][t] = NADBL;
		}
		s++;
	    }
	}

	/* reverse the new columns: most recent first */
	for (j=0; j<compfac/2; j++) {
	    p = k + compfac - j - 1;
	    xtmp = cset->Z[k+j];
	    cset->Z[k+j] = cset->Z[p];
	    cset->Z[p] = xtmp;
	}

	/* names and labels */
	for (j=0; j<compfac; j++) {
	    strcpy(cset->varname[k+j], dset->varname[i]);
	    if (oldpd == 12 && newpd == 4) {
		gretl_trunc(cset->varname[k+j], VNAMELEN - 4);
		sprintf(sfx, "_m%d", compfac - j);
	    } else if (oldpd == 12) {
		/* going to annual */
		gretl_trunc(cset->varname[k+j], VNAMELEN - 5);
		sprintf(sfx, "_m%02d", compfac - j);
	    } else {
		gretl_trunc(cset->varname[k+j], VNAMELEN - 4);
		sprintf(sfx, "_q%d", compfac - j);
	    }
	    strcat(cset->varname[k+j], sfx);
	    sprintf(label, "%s in %s %d of %s", dset->varname[i],
		    p0, compfac - j, p1);
	    series_record_label(cset, k+j, label);
	    series_set_midas_period(cset, k+j, compfac - j);
	    series_set_midas_freq(cset, k+j, oldpd);
	    if (j == 0) {
		series_set_midas_anchor(cset, k+j);
	    }
	}

	/* advance column write position for next source series */
	k += compfac;
    }

#if SPREAD_DEBUG > 1
    PRN *prn = gretl_print_new(GRETL_PRINT_STDERR, NULL);
    printdata(NULL, NULL, cset, OPT_O, prn);
    gretl_print_destroy(prn);
#endif

    return cset;
}

/* specific apparatus for converting daily time series to monthly */

static double *extend_series (const double *z, int n)
{
    double *x = malloc(n * sizeof *x);

    if (x != NULL) {
	int t;

	x[0] = NADBL;
	for (t=1; t<n; t++) {
	    x[t] = z[t-1];
	}
    }

    return x;
}

#define DMDEBUG 0

static double *
daily_series_to_monthly (DATASET *dset, int i,
			 int nm, int yr, int mon, int offset,
			 int any_eop, CompactMethod method)
{
    double *x;
    const double *src = dset->Z[i];
    const double *z;
    double *tmp = NULL;
    int t, sop_t, eop_t;

    x = malloc(nm * sizeof *x);
    if (x == NULL) {
	return NULL;
    }

    if (offset < 0) {
	tmp = extend_series(src, dset->n + 1);
	if (tmp == NULL) {
	    free(x);
	    return NULL;
	}
	/* this permits use of a negative offset */
	z = tmp + 1;
    } else {
	z = src;
    }

    /* Note: we can't necessarily assume that the first obs
       is the first day of a month. The @offset value gives the
       number of daily observations (allowing for the number of
       observed days in the week) in the first month of the daily
       data, prior to the data actually starting.
    */

    /* The "one day's slack" business below, with start-of-period and
       end-of-period compaction, is designed to allow for the
       possibility that the first (or last) day of the month may not
       have been a trading day.
    */

    /* first obs for start-of-period */
    sop_t = offset;

    /* first obs for end-of-period */
    if (sop_t > 0) {
	eop_t = offset - 1;
    } else {
	/* the first obs starts a month */
	eop_t = get_days_in_month(mon, yr, dset->pd, 0) - 1;
    }

#if DMDEBUG
    fprintf(stderr, "starting: offset=%d, any_eop=%d, sop_t=%d, eop_t=%d\n",
	    offset, any_eop, sop_t, eop_t);
#endif

    for (t=0; t<nm; t++) {
	/* loop across the months in the compacted data */
	int mdays = get_days_in_month(mon, yr, dset->pd, 0);

	if (t > 0) {
	    eop_t += mdays;
	}

#if DMDEBUG
	fprintf(stderr, "t=%d: mon=%d, mdays=%d, sop_t=%d, eop_t=%d\n",
		t, mon, mdays, sop_t, eop_t);
#endif

	if (t == 0 && offset > 0 && any_eop && method != COMPACT_EOP) {
	    /* we started with an incomplete month: so any
	       method other than EOP yields an NA */
	    x[t] = NADBL;
	} else if (method == COMPACT_SOP) {
	    /* allow one days's slack */
	    if (na(z[sop_t]) && sop_t < dset->n - 1) {
		x[t] = z[sop_t + 1];
	    } else {
		x[t] = z[sop_t];
	    }
	} else if (method == COMPACT_EOP) {
	    if (eop_t >= dset->n) {
		x[t] = NADBL;
	    } else {
		/* allow one days's slack */
		if (na(z[eop_t]) && eop_t > 0) {
		    x[t] = z[eop_t - 1];
		} else {
		    x[t] = z[eop_t];
		}
	    }
	} else if (method == COMPACT_SUM ||
		   method == COMPACT_AVG) {
	    int j, dayt, den = mdays;
	    int n_ok = 0;

	    x[t] = 0.0;

	    for (j=0; j<mdays; j++) {
		dayt = sop_t + j;
		if (dayt >= dset->n) {
		    x[t] = NADBL;
		    break;
		} else if (na(z[dayt])) {
		    if (method == COMPACT_AVG) {
			den--;
		    }
		} else {
		    /* got a valid observation */
		    x[t] += z[dayt];
		    n_ok++;
		}
	    }

	    if (n_ok == 0) {
		x[t] = NADBL;
	    }

	    if (method == COMPACT_AVG && !na(x[t])) {
		if (den > 0) {
		    x[t] /= (double) den;
		} else {
		    x[t] = NADBL;
		}
	    }
	}

	sop_t += mdays;

	if (mon == 12) {
	    mon = 1;
	    yr++;
	} else {
	    mon++;
	}
    }

    if (tmp != NULL) {
	free(tmp);
    }

    return x;
}

static void
get_startskip_etc (int compfac, int startmin, int endmin,
		   DATASET *dset, CompactMethod method,
		   int *startskip, int *newn)
{
    int oldn = dset->n;
    int ss = 0, n = 0;

    if (method == COMPACT_EOP) {
	int unused;

	ss = (compfac - (startmin % compfac)) % compfac;
	n = oldn / compfac;
	unused = oldn - 1 - ss - (n-1) * compfac;
	if (unused >= compfac) {
	    n++;
	}
    } else if (method == COMPACT_SOP) {
	int unused;

	ss = (compfac - (startmin % compfac) + 1) % compfac;
	n = oldn / compfac;
	unused = oldn - 1 - ss - (n-1) * compfac;
	if (unused >= compfac) {
	    n++;
	}
    } else if (dated_daily_data(dset)) {
        int es = endmin % compfac;

        /* in this case @startmin is actually startskip */
        ss = startmin;
        n = (oldn - ss - es) / compfac;
    } else {
	int es = endmin % compfac;

	ss = (compfac - (startmin % compfac) + 1) % compfac;
	n = (oldn - ss - es) / compfac;
    }

#if DB_DEBUG
    fprintf(stderr, "get_startskip_etc: startskip %d, startmin %d, newn %d\n",
            ss, startmin, n);
#endif

    *startskip = ss;
    *newn = n;
}

/* Note: min_startskip is the minimum, across the series to
   be compacted in case we're doing multiple series, of the
   initial number of observations that cannot be used in the
   compacted series. The per-series startskip will depend on
   the compaction method (which may differ across series).

   For example, consider quarterly to annual compaction where
   the quarterly data start in a third quarter. When doing
   COMPACT_AVG or COMPACT_SUM the initial quarters III and IV
   would be unusable (startskip 2), but with COMPACT_EOP only
   the initial quarter III would be unusable (startskip 1).

   The minimum startskip in effect tells us where the result
   dataset will start in relation to the source dataset.
*/

typedef struct compact_params_ {
    int min_startskip;  /* see Note above */
    int newn;           /* number of observations in compacted dataset */
    int any_eop;        /* boolean: 1 iff any series wants COMPACT_EOP */
    int any_sop;        /* boolean: 1 iff any series wants COMPACT_SOP */
    int all_same;       /* boolean: 1 iff all series want the same method */
} compact_params;

static void compact_params_init (compact_params *cp, int oldpd)
{
    cp->min_startskip = oldpd;
    cp->newn = 0;
    cp->any_sop = 0;
    cp->any_eop = 0;
    cp->all_same = 1;
}

/* used only when compacting daily data to monthly */

static void
get_daily_compact_params (CompactMethod default_method,
			  compact_params *cp,
			  const DATASET *dset)
{
    int i, n_not_eop = 0, n_not_sop = 0;

    cp->all_same = 1;
    cp->any_eop = (default_method == COMPACT_EOP)? 1 : 0;
    cp->any_sop = (default_method == COMPACT_SOP)? 1 : 0;

    for (i=1; i<dset->v; i++) {
	CompactMethod method = series_get_compact_method(dset, i);

	if (method != default_method && method != COMPACT_UNSET) {
	    cp->all_same = 0;
	    if (method == COMPACT_EOP) {
		cp->any_eop = 1;
	    } else {
		n_not_eop++;
	    }
	    if (method == COMPACT_SOP) {
		cp->any_sop = 1;
	    } else {
		n_not_sop++;
	    }
	}
    }

    if (n_not_eop == dset->v - 1) {
	cp->any_eop = 0;
    }
    if (n_not_sop == dset->v - 1) {
	cp->any_sop = 0;
    }
}

/* used in cases other than daily to monthly */

static void
get_global_compact_params (int compfac, int startmin, int endmin,
			   CompactMethod default_method,
			   compact_params *cp, DATASET *dset)
{
    CompactMethod method = default_method;
    int nvars = dset->v - 1;
    int i, startskip, n = 0;
    int n_not_eop = 0;

    get_startskip_etc(compfac, startmin, endmin, dset,
		      method, &startskip, &n);
    if (method == COMPACT_EOP) {
	cp->any_eop = 1;
    }

    for (i=1; i<=nvars; i++) {
	method = series_get_compact_method(dset, i);
	if (method != default_method && method != COMPACT_UNSET) {
	    get_startskip_etc(compfac, startmin, endmin, dset,
			      method, &startskip, &n);
	    cp->all_same = 0;
	    if (method == COMPACT_EOP) {
		cp->any_eop = 1;
	    } else {
		n_not_eop++;
	    }
	}
	if (startskip < cp->min_startskip) {
	    cp->min_startskip = startskip;
	}
	if (n > cp->newn) {
	    cp->newn = n;
	}
    }

    if (n_not_eop == nvars) {
	cp->any_eop = 0;
    }
}

static int get_obs_maj_min (const char *obs, int *maj, int *min)
{
    int np = sscanf(obs, "%d:%d", maj, min);

    if (np < 2) {
	np = sscanf(obs, "%d.%d", maj, min);
    }

    return (np == 2);
}

/* for daily to monthly compaction, figure the number of
   observations to be skipped at the start of each series
*/

static int get_daily_offset (const DATASET *dset,
			     int y, int m, int d,
			     int skip, int any_eop)
{
    int ret = 0;

    if (skip) {
	/* moving to start of next month: offset = no. of
	   observations in the first month */
	ret = days_in_month_after(y, m, d, dset->pd) + 1;
    } else if (any_eop && !day_starts_month(d, m, y, dset->pd, NULL)) {
	ret = days_in_month_after(y, m, d, dset->pd) + 1;
    } else {
	/* offset = no. of obs missing at start of first month */
	ret = days_in_month_before(y, m, d, dset->pd);
#if DMDEBUG
	fprintf(stderr, "days_in_month_before %d-%02d-%02d = %d "
		"for pd=%d\n", y, m, d, ret, dset->pd);
#endif
    }

    return ret;
}

/* for daily to monthly compaction, figure the number of
   monthly observations that can be constructed
*/

static int get_n_ok_months (const DATASET *dset,
			    CompactMethod default_method,
			    int *startyr, int *startmon,
			    int *endyr, int *endmon,
			    int *offset, int *any_eop)
{
    compact_params cp;
    int y1, m1, d1, y2, m2, d2;
    int skip = 0, pad = 0, nm = -1;

    if (sscanf(dset->stobs, YMD_READ_FMT, &y1, &m1, &d1) != 3) {
	return -1;
    }
    if (sscanf(dset->endobs, YMD_READ_FMT, &y2, &m2, &d2) != 3) {
	return -1;
    }

    if (y1 < 100) {
	y1 = FOUR_DIGIT_YEAR(y1);
    }
    if (y2 < 100) {
	y2 = FOUR_DIGIT_YEAR(y2);
    }

    nm = 12 * (y2 - y1) + m2 - m1 + 1;

    compact_params_init(&cp, 0);
    get_daily_compact_params(default_method, &cp, dset);

    *startyr = y1;
    *startmon = m1;
    *endyr = y2;
    *endmon = m2;

#if DMDEBUG
    fprintf(stderr, "get_n_ok_months: any_sop=%d, any_eop=%d, "
	    "all_same=%d\n", cp->any_sop,  cp->any_eop,  cp->all_same);
    fprintf(stderr, "y1=%d m1=%d d1=%d; y2=%d m2=%d d2=%d\n",
	    y1, m1, d1, y2, m2, d2);
#endif

    if (!day_starts_month(d1, m1, y1, dset->pd, &pad) && ! cp.any_eop) {
	if (*startmon == 12) {
	    *startmon = 1;
	    *startyr += 1;
	} else {
	    *startmon += 1;
	}
	skip = 1;
	nm--;
    }

    if (!day_ends_month(d2, m2, y2, dset->pd) && ! cp.any_sop) {
	if (*endmon == 1) {
	    *endmon = 12;
	    *endyr -= 1;
	} else {
	    *endmon -= 1;
	}
	nm--;
    }

#if DMDEBUG
    fprintf(stderr, "after adjustment: range %d:%02d to %d:%02d, "
	    "pad=%d, skip=%d\n", *startyr, *startmon, *endyr, *endmon,
	    pad, skip);
#endif

    if (pad) {
	*offset = -1;
    } else {
	*offset = get_daily_offset(dset, y1, m1, d1, skip, cp.any_eop);
    }

    *any_eop = cp.any_eop;

    return nm;
}

static int
weeks_to_months_exec (double **mZ, const DATASET *dset,
		      CompactMethod method)
{
    char obsstr[OBSLEN];
    int *mn = NULL;
    int yr, mon, day;
    int monbak = 0;
    int i, s, t = 0;
    int err = 0;

    mn = malloc(dset->v * sizeof *mn);
    if (mn == NULL) {
	return E_ALLOC;
    }

    for (i=1; i<dset->v; i++) {
	/* initialize all series, first obs */
	mZ[i][0] = NADBL;
	mn[i] = 0;
    }

    for (s=0; s<dset->n; s++) {
	/* loop across the weekly obs in this month */
	ntolabel(obsstr, s, dset);
	sscanf(obsstr, YMD_READ_FMT, &yr, &mon, &day);
	if (monbak > 0 && mon != monbak) {
	    /* new month: finalize the previous one */
	    for (i=1; i<dset->v; i++) {
		if (method == COMPACT_EOP) {
		    if (s > 0) {
			mZ[i][t] = dset->Z[i][s-1];
		    }
		} else if (method == COMPACT_AVG) {
		    if (mn[i] > 0) {
			mZ[i][t] /= (double) mn[i];
		    }
		}
	    }
	    /* and start another? */
	    if (s < dset->n - 1) {
		t++;
		for (i=1; i<dset->v; i++) {
		    /* initialize all series, current obs */
		    if (method == COMPACT_SOP) {
			mZ[i][t] = dset->Z[i][s];
		    } else {
			mZ[i][t] = NADBL;
		    }
		    mn[i] = 0;
		}
	    }
	}

	/* cumulate non-missing weekly observations? */
	for (i=1; i<dset->v; i++) {
	    if (method == COMPACT_SOP) {
		; /* handled above */
	    } else if (method == COMPACT_EOP) {
		mZ[i][t] = dset->Z[i][s];
	    } else if (!na(dset->Z[i][s])) {
		if (na(mZ[i][t])) {
		    mZ[i][t] = dset->Z[i][s];
		} else {
		    mZ[i][t] += dset->Z[i][s];
		}
		mn[i] += 1;
	    }
	    if (mon == monbak && s == dset->n - 1) {
		/* reached the end: ship out last obs */
		if (method == COMPACT_EOP) {
		    mZ[i][t] = NADBL;
		} else if (method == COMPACT_AVG && mn[i] > 0) {
		    mZ[i][t] /= (double) mn[i];
		}
	    }
	}
	monbak = mon;
    }

    free(mn);

    return err;
}

static int
weeks_to_months_check (const DATASET *dset, int *startyr, int *endyr,
		       int *startmon, int *endmon)
{
    char obsstr[OBSLEN];
    int yr, mon, day;
    int mcount = 0;
    int monbak = 0;
    int t, err = 0;

    for (t=0; t<dset->n; t++) {
	ntolabel(obsstr, t, dset);
	if (sscanf(obsstr, YMD_READ_FMT, &yr, &mon, &day) != 3) {
	    err = 1;
	    break;
	}
	if (monbak == 0) {
	    /* first obs */
	    fprintf(stderr, "starting month = '%d'\n", mon);
	    *startyr = yr;
	    *startmon = mon;
	    mcount++;
	} else if (mon != monbak) {
	    /* got a new month */
	    mcount++;
	}
	monbak = mon;
    }

    if (err) {
	mcount = 0;
    } else {
	/* flush the last observation */
	*endyr = yr;
	*endmon = mon;
    }

    return mcount;
}

static int weekly_dataset_to_monthly (DATASET *dset,
				      CompactMethod method)
{
    DATASET mset;
    int startyr = 1, endyr = 1;
    int startmon = 1, endmon = 1;
    int err = 0;

    mset.n = weeks_to_months_check(dset, &startyr, &endyr, &startmon, &endmon);
    fprintf(stderr, "Weekly data: found %d months\n", mset.n);
    if (mset.n <= 0) {
	return E_DATA;
    }

    mset.v = dset->v;
    err = allocate_Z(&mset, 0);
    if (err) {
	return err;
    }

    /* compact series */
    if (!err && dset->v > 1) {
	err = weeks_to_months_exec(mset.Z, dset, method);
    }

    if (err) {
	free_Z(&mset);
    } else {
	free_Z(dset);
	dset->Z = mset.Z;
	dset->n = mset.n;
	dset->pd = 12;
	sprintf(dset->stobs, "%04d:%02d", startyr, startmon);
	sprintf(dset->endobs, "%04d:%02d", endyr, endmon);
	dset->sd0 = get_date_x(dset->pd, dset->stobs);
	dset->t1 = 0;
	dset->t2 = dset->n - 1;
    }

    return err;
}

static int shorten_the_constant (double **Z, int n)
{
    double *tmp = realloc(Z[0], n * sizeof *tmp);

    if (tmp == NULL) {
	return E_ALLOC;
    } else {
	Z[0] = tmp;
	return 0;
    }
}

/* compaction of daily to weekly data in the special case of
   COMPACT_WDAY -- that is, via use of a "representative day",
   for example the weekly value is Wednesday's value. Note
   that on input @repday is 0-based on Sunday.
*/

static int daily_dataset_to_weekly_special (DATASET *dset,
					    int repday)
{
    int y1, m1, d1;
    char obs[OBSLEN];
    double *x = NULL;
    double *tmp;
    int n = 0, n_ok = 0;
    int wday, ok;
    int i, t, err = 0;

    int repmin = G_DATE_MONDAY;
    int repmax = dset->pd == 5 ? G_DATE_FRIDAY :
	dset->pd == 6 ? G_DATE_SATURDAY : G_DATE_SUNDAY;

    if (repday < repmin || repday > repmax) {
	gretl_errmsg_set("Invalid repday value");
	return E_INVARG;
    }

    for (t=0; t<dset->n; t++) {
	ntolabel(obs, t, dset);
	wday = weekday_from_date(obs);
	if (wday == repday) {
	    ok = 0;
	    for (i=1; i<dset->v; i++) {
		if (!na(dset->Z[i][t])) {
		    ok = 1;
		    break;
		}
	    }
	    if (ok) {
		n_ok++;
	    }
	    if (n == 0) {
		sscanf(obs, YMD_READ_FMT, &y1, &m1, &d1);
	    }
	    n++;
	}
    }

    if (n_ok == 0) {
	gretl_errmsg_set(_("Compacted dataset would be empty"));
	return 1;
    }

    fprintf(stderr, "n=%d, n_ok=%d, y1=%d, m1=%d, d1=%d\n",
	    n, n_ok, y1, m1, d1);

    x = malloc(n * sizeof *x);
    if (x == NULL) {
	return E_ALLOC;
    }

    err = shorten_the_constant(dset->Z, n);

    for (i=1; i<dset->v && !err; i++) {
	int s = 0;

	for (t=0; t<dset->n; t++) {
	    ntolabel(obs, t, dset);
	    wday = weekday_from_date(obs);
	    if (wday == repday) {
		x[s++] = dset->Z[i][t];
	    }
	}
	tmp = realloc(dset->Z[i], n * sizeof *tmp);
	if (tmp == NULL) {
	    err = E_ALLOC;
	} else {
	    dset->Z[i] = tmp;
	    for (t=0; t<n; t++) {
		dset->Z[i][t] = x[t];
	    }
	}
    }

    free(x);

    if (!err) {
	dset->n = n;
	dset->pd = 52;
	sprintf(dset->stobs, YMD_WRITE_Y4_FMT, y1, m1, d1);
	dset->sd0 = get_date_x(dset->pd, dset->stobs);
	dset->t1 = 0;
	dset->t2 = dset->n - 1;
	ntolabel(dset->endobs, dset->t2, dset);
	dataset_destroy_obs_markers(dset);
    }

    return err;
}

static int daily_dataset_to_monthly (DATASET *dset,
				     CompactMethod default_method)
{
    int nm, startyr, startmon, endyr, endmon;
    int offset, any_eop;
    CompactMethod method;
    double *x;
    int i, err = 0;

    nm = get_n_ok_months(dset, default_method, &startyr, &startmon,
			 &endyr, &endmon, &offset, &any_eop);

    if (nm <= 0) {
	gretl_errmsg_set(_("Compacted dataset would be empty"));
	return E_DATA;
    }

    err = shorten_the_constant(dset->Z, nm);

    for (i=1; i<dset->v && !err; i++) {
	method = series_get_compact_method(dset, i);
	if (method == COMPACT_UNSET) {
	    method = default_method;
	}
	x = daily_series_to_monthly(dset, i, nm,
				    startyr, startmon,
				    offset, any_eop, method);
	if (x == NULL) {
	    err = E_ALLOC;
	} else {
	    free(dset->Z[i]);
	    dset->Z[i] = x;
	}
    }

    if (!err) {
	dset->n = nm;
	dset->pd = 12;
	sprintf(dset->stobs, "%04d:%02d", startyr, startmon);
	sprintf(dset->endobs, "%04d:%02d", endyr, endmon);
	dset->sd0 = get_date_x(dset->pd, dset->stobs);
	dset->t1 = 0;
	dset->t2 = dset->n - 1;
	dataset_destroy_obs_markers(dset);
    }

    return err;
}

static int get_daily_skip (const DATASET *dset, int t)
{
    int dd = calendar_obs_number(dset->S[t], dset, 0) -
	calendar_obs_number(dset->S[t-1], dset, 0);

    if (dd == 0) {
	fprintf(stderr, "get_daily_skip: S[%d] = '%s', S[%d] = '%s'\n",
		t, dset->S[t], t-1, dset->S[t-1]);
    }

    return dd - 1;
}

static int insert_missing_hidden_obs (DATASET *dset, int nmiss)
{
    int oldn = dset->n;
    double *tmp, **Z;
    int i, s, t, skip;
    int err = 0;

    err = dataset_add_observations(dset, nmiss, OPT_NONE);
    if (err) {
	return err;
    }

#if DB_DEBUG
    fprintf(stderr, "daily data: expanded n from %d to %d\n",
	    oldn, dset->n);
#endif

    Z = dset->Z;
    tmp = Z[0];

    for (i=1; i<dset->v && !err; i++) {
	for (s=0; s<oldn; s++) {
	    tmp[s] = Z[i][s];
	}

	Z[i][0] = tmp[0];
	t = 1;
	for (s=1; s<oldn; s++) {
	    skip = get_daily_skip(dset, s);
	    if (skip < 0) {
		err = E_DATA;
		break;
	    }
	    while (skip--) {
		Z[i][t++] = NADBL;
	    }
	    Z[i][t++] = tmp[s];
	}
    }

    for (t=0; t<dset->n; t++) {
	Z[0][t] = 1.0;
	if (dset->S != NULL) {
	    calendar_date_string(dset->S[t], t, dset);
	}
    }

    if (!err) {
	dset->t2 = dset->n - 1;
	ntolabel(dset->endobs, dset->n - 1, dset);
    }

#if DB_DEBUG > 1
    fprintf(stderr, "insert_missing_hidden_obs, done, err = %d\n", err);
    for (t=0; t<dset->n; t++) {
	fprintf(stderr, "Z[1][%d] = %14g\n", t, Z[1][t]);
    }
#endif

    return err;
}

int maybe_expand_daily_data (DATASET *dset)
{
    int nmiss = n_hidden_missing_obs(dset, 0, dset->n - 1);
    int err = 0;

    fprintf(stderr, "n_hidden_missing_obs: nmiss = %d\n", nmiss);

    if (nmiss < 0) {
	err = 1;
    } else if (nmiss > 0) {
	err = insert_missing_hidden_obs(dset, nmiss);
    }

    return err;
}

static int do_compact_spread (DATASET *dset, int newpd)
{
    DATASET *cset = NULL;
    int nv = 0;
    int err;

    err = compact_spread_pd_check(dset->pd, newpd);
    if (err) {
	return err;
    }

    if (dated_daily_data(dset)) {
	err = maybe_expand_daily_data(dset);
	if (err) {
	    gretl_errmsg_set(_("Error expanding daily data with missing observations"));
	} else {
	    cset = compact_daily_spread(dset, newpd, &nv, &err);
	}
    } else {
	int startmaj, startmin;
	int endmaj, endmin;

	/* get starting obs major and minor components */
	if (!get_obs_maj_min(dset->stobs, &startmaj, &startmin)) {
	    return E_DATA;
	}

	/* get ending obs major and minor components */
	if (!get_obs_maj_min(dset->endobs, &endmaj, &endmin)) {
	    return E_DATA;
	}

	cset = compact_data_spread(dset, newpd, startmaj, startmin,
				   endmaj, endmin, &nv, &err);
    }

    if (!err) {
	free_Z(dset);
	clear_datainfo(dset, CLEAR_FULL);
	*dset = *cset;
	free(cset);
    }

    return err;
}

static int dated_daily_startmin (DATASET *dset, int wkstart)
{
    char obs[OBSLEN];
    int t, smin = 0;

    for (t=0; t<dset->pd; t++) {
        ntolabel(obs, t, dset);
        if (weekday_from_date(obs) == wkstart) {
            smin = t;
            break;
        }
    }

    return smin;
}

#define is_daily(p) (p==5 || p==6 || p==7)

/**
 * compact_dataset:
 * @dset: dataset struct.
 * @newpd: target data frequency.
 * @default_method: code for the default compaction method.
 * @wkstart: the day on which a week is deemed to start:
 * 1 for Monday or 7 for Sunday; relevant only for daily
 * to weekly compaction.
 * @repday: "representative day" for conversion from daily
 * to weekly data (with method %COMPACT_WDAY only). Value
 * 1 to 7 with 1 = Monday, or 0 to ignore it.
 *
 * Compact the data set from higher to lower frequency.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int compact_dataset (DATASET *dset, int newpd,
		     CompactMethod default_method,
		     int wkstart, int repday)
{
    compact_params cp;
    int oldn = dset->n;
    int oldpd = dset->pd;
    int compfac = 0;
    int start_major;
    int start_minor;
    int end_major = 0;
    int end_minor = 0;
    int from_dated_daily;
    char stobs[OBSLEN];
    int i, err = 0;

    gretl_error_clear();
    compact_params_init(&cp, oldpd);

    if (wkstart == 0) {
	; /* OK, "ignore me" */
    } else if (wkstart != G_DATE_MONDAY && wkstart != G_DATE_SUNDAY) {
	gretl_errmsg_set("Invalid weekstart value");
	return E_INVARG;
    }

    /* catch a couple of special cases */
    if (default_method == COMPACT_SPREAD) {
	return do_compact_spread(dset, newpd);
    } else if (oldpd == 52) {
	return weekly_dataset_to_monthly(dset, default_method);
    }

    from_dated_daily = dated_daily_data(dset);

    if (from_dated_daily) {
	/* allow for the possibility that the daily dataset
	   contains "hidden" or suppressed missing observations
	   (holidays are just skipped, not marked as NA)
	*/
	err = maybe_expand_daily_data(dset);
	if (err) {
	    gretl_errmsg_set(_("Error expanding daily data with missing observations"));
	    return err;
	} else {
	    oldn = dset->n;
	}
    }

    if (newpd == 52 && is_daily(oldpd) && default_method == COMPACT_WDAY) {
	/* daily to weekly, using "representative day" */
	return daily_dataset_to_weekly_special(dset, repday);
    } else if (newpd == 12 && is_daily(oldpd)) {
	/* daily to monthly: special */
	return daily_dataset_to_monthly(dset, default_method);
    } else if (is_daily(oldpd)) {
	/* daily to weekly, not COMPACT_WDAY */
	compfac = oldpd;
	if (from_dated_daily) {
            /* here we're borrowing @start_minor for a daily-
               specific purpose
            */
	    start_minor = dated_daily_startmin(dset, wkstart);
	} else {
	    start_minor = 1;
	}
    } else if (oldpd == 24 && is_daily(newpd)) {
	/* hourly to daily */
	compfac = 24;
	if (!get_obs_maj_min(dset->stobs, &start_major, &start_minor)) {
	    return 1;
	}
    } else {
	compfac = oldpd / newpd;
	/* get starting obs major and minor components */
	if (!get_obs_maj_min(dset->stobs, &start_major, &start_minor)) {
	    return 1;
	}
	/* get ending obs major and minor components */
	if (!get_obs_maj_min(dset->endobs, &end_major, &end_minor)) {
	    return 1;
	}
    }

    get_global_compact_params(compfac, start_minor, end_minor,
                              default_method, &cp, dset);

    if (cp.newn == 0 && default_method != COMPACT_SPREAD) {
	gretl_errmsg_set(_("Compacted dataset would be empty"));
	return 1;
    }

    if (newpd == 1) {
	if (cp.min_startskip > 0 && !cp.any_eop) {
	    start_major++;
	}
	sprintf(stobs, "%d", start_major);
    } else if (newpd == 52) {
	if (dated_daily_data(dset)) {
	    ntolabel(stobs, cp.min_startskip, dset);
	} else {
	    strcpy(stobs, "1");
	}
    } else {
	int m0 = start_minor + cp.min_startskip;
	int minor = m0 / compfac + (m0 % compfac > 0);

	if (minor > newpd) {
	    start_major++;
	    minor -= newpd;
	}
	format_obs(stobs, start_major, minor, newpd);
    }

    /* revise datainfo members */
    strcpy(dset->stobs, stobs);
    dset->pd = newpd;
    dset->n = cp.newn;
    dset->sd0 = get_date_x(dset->pd, dset->stobs);
    dset->t1 = 0;
    dset->t2 = dset->n - 1;
    ntolabel(dset->endobs, dset->t2, dset);

    if (from_dated_daily) {
	/* delete any daily date strings */
	if (dset->S != NULL) {
	    dataset_destroy_obs_markers(dset);
	}
	/* revise endobs */
	ntolabel(dset->endobs, dset->t2, dset);
    }

    err = shorten_the_constant(dset->Z, dset->n);

    /* compact the individual data series */
    for (i=1; i<dset->v && !err; i++) {
	CompactMethod this_method = default_method;
	int startskip = cp.min_startskip;
	double *x;

	if (!cp.all_same) {
	    CompactMethod m_i = series_get_compact_method(dset, i);

	    if (m_i != COMPACT_UNSET) {
		this_method = m_i;
	    }

	    startskip = compfac - (start_minor % compfac) + 1;
	    startskip = startskip % compfac;
	    if (this_method == COMPACT_EOP) {
		if (startskip > 0) {
		    startskip--;
		} else {
		    startskip = compfac - 1;
		}
	    }
	}
	x = compact_series(dset, i, oldn, startskip, cp.min_startskip,
			   compfac, this_method);
	if (x == NULL) {
	    err = E_ALLOC;
	} else {
	    free(dset->Z[i]);
	    dset->Z[i] = x;
	}
    }

    return err;
}

/**
 * expand_dataset:
 * @dset: dataset struct.
 * @newpd: target data frequency.
 *
 * Expand the data set from lower to higher frequency: an "expert"
 * option.  This is supported only for expansion from annual
 * to quarterly or monthly, or from quarterly to monthly.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int expand_dataset (DATASET *dset, int newpd)
{
    char stobs[OBSLEN];
    int oldn = dset->n;
    int oldpd = dset->pd;
    int t1 = dset->t1;
    int t2 = dset->t2;
    int mult, newn, nadd;
    double *x = NULL;
    size_t sz;
    int i, j, t, s;
    int err = 0;

    if (oldpd != 1 && oldpd != 4) {
	return E_PDWRONG;
    } else if (oldpd == 1 && newpd != 4 && newpd != 12) {
	return E_DATA;
    } else if (oldpd == 4 && newpd != 12) {
	return E_DATA;
    }

    x = malloc(oldn * sizeof *x);
    if (x == NULL) {
	return E_ALLOC;
    }

    mult = newpd / oldpd;  /* frequency increase factor */
    newn = mult * dset->n; /* revised number of observations */
    nadd = newn - oldn;    /* number of obs to add */

    err = dataset_add_observations(dset, nadd, OPT_D);
    if (err) {
	goto bailout;
    }

    sz = oldn * sizeof *x;
    for (i=1; i<dset->v; i++) {
	memcpy(x, dset->Z[i], sz);
	s = 0;
	for (t=0; t<oldn; t++) {
	    for (j=0; j<mult; j++) {
		dset->Z[i][s++] = x[t];
	    }
	}
	series_set_orig_pd(dset, i, oldpd);
    }

    if (dset->pd == 1) {
	/* starting with annual data */
	strcpy(stobs, dset->stobs);
	if (newpd == 4) {
	    strcat(stobs, ":1");
	} else {
	    strcat(stobs, ":01");
	}
    } else {
	/* starting with quarterly data */
	int yr, qtr, mo;

	sscanf(dset->stobs, "%d:%d", &yr, &qtr);
	mo = (qtr - 1) * 3 + 1;
	sprintf(stobs, "%d:%02d", yr, mo);
    }

    /* revise the sample range, if set */
    if (dset->t1 > 0) {
	dset->t1 *= mult;
    }
    if (dset->t2 < oldn - 1) {
	dset->t2 = dset->t1 + (t2 - t1 + 1) * mult - 1;
    }

    strcpy(dset->stobs, stobs);
    dset->pd = newpd;
    dset->sd0 = get_date_x(dset->pd, dset->stobs);
    ntolabel(dset->endobs, dset->n - 1, dset);

 bailout:

    free(x);

    return err;
}
