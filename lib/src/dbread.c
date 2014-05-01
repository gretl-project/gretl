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

static int cli_add_db_data (double **dbZ, SERIESINFO *sinfo, 
			    DATASET *dset, CompactMethod method, 
			    int dbv, int *newdata, PRN *prn);

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

static int db_row_wanted (const gretl_matrix *mask, 
			  int masklen, int t)
{
    if (t >= masklen) {
	return 0;
    } else {
	return gretl_vector_get(mask, t) != 0;
    }
}

static int get_native_db_data_masked (const char *dbbase, 
				      SERIESINFO *sinfo,
				      double **dbZ,
				      const gretl_matrix *mask)
{
    char numstr[32];
    FILE *fp;
    dbnumber x;
    int s, t, masklen;
    int err = 0;

    fp = open_binfile(dbbase, GRETL_NATIVE_DB, sinfo->offset, &err);
    if (err) {
	return err;
    }

    masklen = gretl_vector_get_length(mask);

    s = 0;
    for (t=0; t<=sinfo->nobs && !err; t++) {
	if (fread(&x, sizeof x, 1, fp) != 1) {
	    err = DB_PARSE_ERROR;
	} else if (db_row_wanted(mask, masklen, t)) {
	    sprintf(numstr, "%.7g", (double) x);
	    dbZ[1][s] = atof(numstr);
	    if (dbZ[1][s] == DBNA) {
		dbZ[1][s] = NADBL;
	    }
	    s++;
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
 * Returns:  0 on success, non-zero code on failure.
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

#if G_BYTE_ORDER == G_BIG_ENDIAN
    err = retrieve_remote_db_data(dbbase, sinfo->varname, &getbuf,
				  GRAB_NBO_DATA);
#else
    err = retrieve_remote_db_data(dbbase, sinfo->varname, &getbuf,
				  GRAB_DATA);
#endif

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

    *sinfo->descrip = 0;
    strncat(sinfo->descrip, s, MAXLABEL - 1);
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
    if (strchr(stobs, '/')) { 
	/* calendar data */
	const char *q = stobs;
	const char *p = strchr(stobs, '/');

	if (p - q == 4) {
	    strcpy(sinfo->stobs, q + 2);
	}
	q = endobs;
	p = strchr(endobs, '/');
	if (p && p - q == 4) {
	    strcpy(sinfo->endobs, q + 2);
	}
    } else {
	*sinfo->stobs = 0;
	*sinfo->endobs = 0;
	strncat(sinfo->stobs, stobs, 8);
	strncat(sinfo->endobs, endobs, 8);
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

static int 
get_native_series_info (const char *series, SERIESINFO *sinfo)
{
    FILE *fp = NULL;
    char sername[VNAMELEN];
    char s1[256], s2[72];
    char stobs[16], endobs[16];
    char pdc;
    int offset = 0;
    int gotit = 0, err = 0;
    int n;

    err = open_native_db_files(saved_db_name, &fp, NULL, NULL, NULL);
    if (err) {
	return err;
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

static int 
get_remote_series_info (const char *series, SERIESINFO *sinfo)
{
    char *buf = NULL;
    char sername[VNAMELEN];
    char s1[256], s2[72];
    char stobs[16], endobs[16];
    char pdc;
    int gotit = 0;
    int err = 0;

    err = retrieve_remote_db_index(saved_db_name, &buf);
    if (err) {
	return err;
    }

    sinfo->offset = 0;

    bufgets_init(buf);

    while (bufgets(s1, sizeof s1, buf) && !gotit) {
	if (*s1 == '#') {
	    continue;
	}

	if (gretl_scan_varname(s1, sername) != 1) {
	    break;
	}

	if (strcmp(series, sername)) {
	    continue;
	}

	gotit = 1;

	if (bufgets(s2, sizeof s2, buf) == NULL) {
	    err = DB_PARSE_ERROR;
	    break;
	} 

	strcpy(sinfo->varname, sername);
	get_native_series_comment(sinfo, s1);

	if (sscanf(s2, "%c %10s %*s %10s %*s %*s %d", 
		   &pdc, stobs, endobs, &sinfo->nobs) != 4) {
	    gretl_errmsg_set(_("Failed to parse series information"));
	    err = DB_PARSE_ERROR;
	} else {
	    get_native_series_pd(sinfo, pdc);
	    get_native_series_obs(sinfo, stobs, endobs);
	    sinfo->t2 = sinfo->nobs - 1;
	}
    }

    bufgets_finalize(buf);

    free(buf);

    if (!gotit) {
	gretl_errmsg_sprintf(_("Series not found, '%s'"), series);
	err = DB_NO_SUCH_SERIES;
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
	    *sinfo->descrip = 0;
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
			int rem = MAXLABEL - strlen(sinfo->descrip) - 1;

			if (rem > 0) {
			    gretl_strstrip(line);
			    strncat(sinfo->descrip, line + 1, rem);
			}
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
    char pdstr[4] = {0}; 
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
	fprintf(stderr, I_("frequency (%d) does not make seem to make sense"),
		(int) dinfo->info);
	fputc('\n', stderr);
	gretl_errmsg_sprintf(("frequency (%d) does not make seem to make sense"), 
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
    strncat(sinfo->descrip, comment, MAXLABEL - 1);

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
	fprintf(stderr, I_("frequency %d is not supported"), pd);
	fputc('\n', stderr);
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
	    strcpy(sinfo->descrip, comment);
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

    *sinfo->varname = '\0';
    *sinfo->descrip = '\0';
    *sinfo->stobs = '\0';
    *sinfo->endobs = '\0';
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
	    gretl_errmsg_set("This is not a PcGive 700 data file");
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
				 CompactMethod method, int compfac,
				 int skip)
{
    int p, t;
    double *x;

    x = malloc(n * sizeof *x);
    if (x == NULL) return NULL;

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

double *compact_db_series (const double *src, SERIESINFO *sinfo,
			   int target_pd, CompactMethod method)
{
    int p0, y0, endskip, goodobs;
    int skip = 0, compfac = sinfo->pd / target_pd;
    double *x;

    if (target_pd == 1) {
	/* figure the annual dates */
	y0 = atoi(sinfo->stobs);
	p0 = atoi(sinfo->stobs + 5);
	if (p0 != 1) {
	    ++y0;
	    skip = compfac - (p0 + 1);
	}
	sprintf(sinfo->stobs, "%d", y0);
    } else if (target_pd == 4) {
	/* figure the quarterly dates */
	float q;
	int q0;

	y0 = atoi(sinfo->stobs);
	p0 = atoi(sinfo->stobs + 5);
	q = 1.0 + p0 / 3.;
	q0 = q + .5;
	skip = ((q0 - 1) * 3) + 1 - p0;
	if (q0 == 5) {
	    y0++;
	    q0 = 1;
	}
	sprintf(sinfo->stobs, "%d.%d", y0, q0);
    } else {
	return NULL;
    }

    endskip = (sinfo->nobs - skip) % compfac;
    goodobs = (sinfo->nobs - skip - endskip) / compfac;
    sinfo->nobs = goodobs;

#if DB_DEBUG
    fprintf(stderr, "startskip = %d\n", skip);
    fprintf(stderr, "endskip = %d\n", endskip);
    fprintf(stderr, "goodobs = %d\n", goodobs);
    fprintf(stderr, "compfac = %d\n", compfac);
    fprintf(stderr, "starting date = %s\n", sinfo->stobs);
#endif

    x = get_compacted_xt(src, goodobs, method, compfac, skip);

    sinfo->pd = target_pd;

    return x;
}

static double *interpolate_db_series (const double *src,
				      int oldn, int mult,
				      int *err)
{
    gretl_matrix *yx;
    gretl_matrix *y;
    double *ret = NULL;
    int t;

    y = gretl_column_vector_alloc(oldn);
    if (y == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    for (t=0; t<oldn; t++) {
	y->val[t] = src[t];
    }

    yx = matrix_chowlin(y, NULL, mult, err);
    gretl_matrix_free(y);

    if (!*err) {
	ret = yx->val;
	yx->val = NULL;
    }

    gretl_matrix_free(yx);

    return ret;
}

/* Expand a single series from a database, for importation 
   into a working dataset of higher frequency.  At present
   this is permitted only for the cases:

   1) annual    -> quarterly
   2) annual    -> monthly
   3) quarterly -> monthly

   Interpolation is supported for cases 1 and 3 only.
*/

double *expand_db_series (const double *src, SERIESINFO *sinfo,
			  int target_pd, int interpol)
{
    char stobs[12] = {0};
    int oldn = sinfo->nobs;
    int mult, newn;
    double *x = NULL;
    int j, t;
    int err = 0;

    mult = target_pd / sinfo->pd;
    newn = mult * sinfo->nobs;

    if (!((target_pd == 4 && sinfo->pd == 1) ||
	  (target_pd == 12 && sinfo->pd == 4))) {
	interpol = 0;
    }

    if (interpol) {
	x = interpolate_db_series(src, oldn, mult, &err);
    } else {
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
    }

    if (err) {
	return NULL;
    }

    if (sinfo->pd == 1) {
	strcpy(stobs, sinfo->stobs);
	if (target_pd == 4) {
	    strcat(stobs, ":1");
	} else {
	    strcat(stobs, ":01");
	}
    } else {
	int yr, qtr, mo;

	sscanf(sinfo->stobs, "%d.%d", &yr, &qtr);
	mo = (qtr - 1) * 3 + 1;
	sprintf(stobs, "%d:%02d", yr, mo);
    }

    strcpy(sinfo->stobs, stobs);
    sinfo->pd = target_pd;
    sinfo->nobs = newn;

    return x;
}

int set_db_name (const char *fname, int filetype, PRN *prn)
{
    FILE *fp;
    int err = 0;

    *saved_db_name = '\0';
    strncat(saved_db_name, fname, MAXLEN - 1);

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
	    build_path(saved_db_name, path, fname, NULL);
	    fp = gretl_fopen(saved_db_name, "rb");
	}

#ifdef OS_OSX
	if (fp == NULL) {
	    gchar *tmp = g_build_path(SLASHSTR, gretl_app_support_dir(), "db",
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

static void ODBC_info_clear_all (void)
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
    void *handle = NULL;
    int (*check_dsn) (ODBC_info *);
    char *dbname = NULL;
    char *uname = NULL;
    char *pword = NULL;
    int got_plugin = 0;
    int err = 0;

    /* skip command word */
    line += strcspn(line, " ");
    line += strspn(line, " ");

    ODBC_info_clear_all();

    dbname = get_dsn_field("dsn", line);
    if (dbname == NULL) {
	pputs(prn, "You must specify a DSN using 'dsn=...'\n");
	return E_DATA;
    }

    uname = get_dsn_field("user", line);
    pword = get_dsn_field("password", line);

    gretl_odinfo.dsn = dbname;
    gretl_odinfo.username = uname;
    gretl_odinfo.password = pword;

    gretl_error_clear();

    check_dsn = get_plugin_function("gretl_odbc_check_dsn", &handle);

    if (check_dsn == NULL) {
        err = 1;
    } else {
	got_plugin = 1;
        err = (* check_dsn) (&gretl_odinfo);
        close_plugin(handle);
    }

    if (err) {
	if (!got_plugin) {
	    pprintf(prn, "Couldn't open the gretl ODBC plugin\n");
	} else {	    
	    pprintf(prn, "Failed to connect to ODBC data source '%s'\n", 
		    gretl_odinfo.dsn);
	} 
	ODBC_info_clear_all();
    } else if (gretl_messages_on()) {
	pprintf(prn, "Connected to ODBC data source '%s'\n", 
		gretl_odinfo.dsn);
    }

    return err;
}

int db_set_sample (const char *s, DATASET *dset)
{
    char start[OBSLEN], stop[OBSLEN];
    int t1 = 0, t2 = 0;

    if (sscanf(s, "%10s %10s", start, stop) != 2) {
	gretl_errmsg_set(_("error reading smpl line"));
	return 1;
    }

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

static char *
get_word_and_advance (char *s, char *word, size_t maxlen)
{
    size_t i = 0;

    while (isspace(*s)) s++;

    *word = 0;

    while (*s && !isspace(*s)) {
	if (i < maxlen) word[i++] = *s;
	s++;
    }

    word[i] = 0;

    if (*word) return s;
    else return NULL;
}

static char *
get_compact_method_and_advance (char *s, CompactMethod *method)
{
    char *p;

    *method = COMPACT_NONE;

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
    fprintf(stderr, "set coltype[%d] = %d, fmt='%s'\n", i, 
	    gretl_odinfo.coltypes[i], gretl_odinfo.fmts[i]);
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

static char *odbc_get_query (char *s, int *err)
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

static char **odbc_get_varnames (char **line, int *err)
{
    char **vnames = NULL;
    char vname[VNAMELEN];
    char *s = *line;
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

static int odbc_transcribe_data (char **vnames, DATASET *dset, 
				 int vmin, int newvars)
{
    char label[MAXLABEL];
    int nv = gretl_odinfo.nvars;
    int n = gretl_odinfo.nrows;
    int nrepl = nv - newvars;
    int i, s, t, v;

    for (i=0; i<nv; i++) {
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
	    v = vmin++;
	    strcpy(dset->varname[v], vnames[i]);
	    sprintf(label, "ODBC series %d", i + 1);
	    series_set_label(dset, v, label);
	} else {
	    vnew = 0;
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
		if (t >= 0 && t < dset->n) {
		    dset->Z[v][t] = gretl_odinfo.X[i][s];
		    obs_used++;
		} else {
		    fprintf(stderr, "Rejecting obs '%s'\n", gretl_odinfo.S[s]);
		}
	    }
	} else {
	    /* no obs identifiers via ODBC */
	    s = 0;
	    for (t=0; t<dset->n; t++) {
		if (t >= dset->t1 && t <= dset->t2 && s < n) {
		    dset->Z[v][t] = gretl_odinfo.X[i][s++];
		    obs_used++;
		} else if (vnew) {
		    dset->Z[v][t] = NADBL;
		}
	    }
	}

	if (vnew && obs_used == 0) {
	    gretl_warnmsg_sprintf("ODBC import: '%s': no valid observations",
				  vnames[i]);
	}
    }

    return 0;
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

static int odbc_get_series (char *line, DATASET *dset, 
			    PRN *prn)
{
    void *handle = NULL;
    int (*get_data) (ODBC_info *);
    char **vnames = NULL;
    char *format = NULL;
    int err = 0;

    if (gretl_odinfo.dsn == NULL) {
	gretl_errmsg_set(_("No database has been opened"));
	return 1;
    } 

    if (dset->n == 0) {
	gretl_errmsg_set(_("No series length has been defined"));
	return 1;
    }

    /* skip "data" plus following space */
    line += strcspn(line, " ");
    line += strspn(line, " ");

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
	fprintf(stderr, "SQL query: '%s'\n", gretl_odinfo.query);
	gretl_error_clear();

	get_data = get_plugin_function("gretl_odbc_get_data", &handle);

	if (get_data == NULL) {
	    err = 1;
	} else {
	    err = (*get_data) (&gretl_odinfo);
	    close_plugin(handle);
	}
    }

    if (!err) {
	int n = gretl_odinfo.nrows;
	int nv = gretl_odinfo.nvars;
	int newvars, vmin = 1;

	if (gretl_messages_on()) {
	    pprintf(prn, "Retrieved %d observations on %d series via ODBC\n", 
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
	    odbc_transcribe_data(vnames, dset, vmin, newvars);
	}
    }

    strings_array_free(vnames, gretl_odinfo.nvars);
    ODBC_info_clear_read();

    return err;
}

static int db_get_row_mask (const gretl_matrix **pmat)
{
    const char *rows;
    int err = 0;

    rows = get_optval_string(DATA, OPT_M);
    if (rows == NULL || *rows == '\0') {
	return E_PARSE;
    }   

    *pmat = get_matrix_by_name(rows);
    if (*pmat == NULL) {
	gretl_errmsg_sprintf(_("'%s': no such matrix"), rows);
	err = E_DATA;
    } else if (gretl_vector_get_length(*pmat) == 0) {
	err = E_NONCONF;
    }

    return err;
}

static int db_n_from_row_mask (const gretl_matrix *mask,
			       int n_max)
{
    int mlen = gretl_vector_get_length(mask);
    int i, n = 0;

    for (i=0; i<mlen && i<n_max; i++) {
	if (gretl_vector_get(mask, i) != 0) {
	    n++;
	}
    }

    return n;
}

/* when a row mask is applied in reading from a gretl
   database, we rejig the series information such that
   it's flat (undated) with the observations index
   running from 1 to the number of rows selected
*/

static int update_sinfo_masked (SERIESINFO *sinfo, int nobs)
{
    int err = 0;

    sinfo->nobs = nobs;
    sinfo->t1 = 0;
    sinfo->t2 = nobs - 1;

    if (sinfo->pd != 1) {
	err = E_PDWRONG;
    } else if (strcmp(sinfo->stobs, "1")) {
	err = E_DATA;
    } else {
	sprintf(sinfo->endobs, "%d", nobs);
    }

    return err;
}

/* main function for getting a series out of a database, using the
   command-line client or in script or console mode
*/

int db_get_series (char *line, DATASET *dset, 
		   gretlopt opt, PRN *prn)
{
    const gretl_matrix *rowmask = NULL;
    char series[VNAMELEN];
    CompactMethod method;
    SERIESINFO sinfo;
    double **dbZ;
    int err = 0;

    if (opt & OPT_O) {
	return odbc_get_series(line, dset, prn);
    }

    if (opt & OPT_M) {
	err = db_get_row_mask(&rowmask);
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

    line = get_compact_method_and_advance(line, &method);

    /* now loop over variable names given on the line */

    while ((line = get_word_and_advance(line, series, VNAMELEN-1)) && !err) {
	int v, this_var_method = method; 
	int nobs, newdata = 0;

	series_info_init(&sinfo);

	/* see if the series is already in the dataset */
	v = series_index(dset, series);
	if (v < dset->v && method == COMPACT_NONE) {
	    this_var_method = series_get_compact_method(dset, v);
	}

#if DB_DEBUG
	fprintf(stderr, "db_get_series: dset->v = %d, v = %d, series = '%s'\n", 
		dset->v, v, series);
	fprintf(stderr, "this_var_method = %d\n", this_var_method);
#endif

	/* find the series information in the database */
	if (saved_db_type == GRETL_RATS_DB) {
	    err = get_rats_series_info(series, &sinfo);
	} else if (saved_db_type == GRETL_PCGIVE_DB) {
	    err = get_pcgive_series_info(series, &sinfo);
#ifdef USE_CURL
	} else if (saved_db_type == GRETL_NATIVE_DB_WWW) {
	    err = get_remote_series_info(series, &sinfo);
#endif
	} else {
	    err = get_native_series_info(series, &sinfo);
	} 

	if (err) {
	    fprintf(stderr, "db_get_series: failed to get series info\n");
	    return 1;
	}

	if (rowmask != NULL) {
	    nobs = db_n_from_row_mask(rowmask, sinfo.nobs);
	} else {
	    nobs = sinfo.nobs;
	}

	/* temporary dataset */
	dbZ = new_dbZ(nobs);
	if (dbZ == NULL) {
	    gretl_errmsg_set(_("Out of memory!"));
	    return 1;
	}

#if DB_DEBUG
	fprintf(stderr, "db_get_series: offset=%d, nobs=%d\n", 
		sinfo.offset, nobs);
#endif

	if (saved_db_type == GRETL_RATS_DB) {
	    err = get_rats_db_data(saved_db_name, &sinfo, dbZ);
	} else if (saved_db_type == GRETL_PCGIVE_DB) {
	    err = get_pcgive_db_data(saved_db_name, &sinfo, dbZ);
	} else if (saved_db_type == GRETL_NATIVE_DB_WWW) {
	    err = get_remote_db_data(saved_db_name, &sinfo, dbZ);
	} else if (rowmask != NULL) {
	    err = get_native_db_data_masked(saved_db_name, &sinfo, dbZ, 
					    rowmask);
	} else {
	    err = get_native_db_data(saved_db_name, &sinfo, dbZ);
	} 

#if DB_DEBUG
	fprintf(stderr, "db_get_series: get_XXX_db_data gave %d\n", err);
#endif

	if (err == DB_MISSING_DATA) {
	    fprintf(stderr, "There were missing data\n");
	    err = 0;
	}

	if (!err && rowmask != NULL) {
	    err = update_sinfo_masked(&sinfo, nobs);
	}

	if (!err) {
	    err = cli_add_db_data(dbZ, &sinfo, dset, 
				  this_var_method, v, 
				  &newdata, prn);
	}

	/* free up temp stuff */
	free_dbZ(dbZ);

	if (!err && !(opt & OPT_Q)) {
	    pprintf(prn, _("Series imported OK"));
	    pputc(prn, '\n');
	    if (newdata) {
		print_smpl(dset, 0, prn);
	    }
	}
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

static int db_delete_series (char *line, const int *list,
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
	    gretl_errmsg_set("This only works for gretl databases");
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
	pprintf(prn, "Deleted %d series from %s\n", ndel, src2);
    }

    return err;
}

int db_delete_series_by_name (char *line, PRN *prn)
{
    return db_delete_series(line, NULL, NULL, prn);
}

int db_delete_series_by_number (const int *list, const char *fname)
{
    return db_delete_series(NULL, list, fname, NULL);
}

void get_db_padding (SERIESINFO *sinfo, DATASET *dset, 
		     int *pad1, int *pad2)
{
    *pad1 = dateton(sinfo->stobs, dset); 
    fprintf(stderr, "pad1 = dateton(%s) = %d\n", sinfo->stobs, *pad1);

    *pad2 = dset->n - sinfo->nobs - *pad1;
    fprintf(stderr, "pad2 = dset->n - sinfo->nobs - pad1 = %d - %d - %d = %d\n", 
	    dset->n, sinfo->nobs, *pad1, *pad2);
} 

int db_range_check (SERIESINFO *sinfo, DATASET *dset)
{
    double sdn_orig = get_date_x(dset->pd, dset->endobs);
    double sd0 = get_date_x(sinfo->pd, sinfo->stobs);
    double sdn = get_date_x(sinfo->pd, sinfo->endobs);
    int err = 0;

    if (sd0 > sdn_orig || sdn < dset->sd0) {
	gretl_errmsg_sprintf(_("%s: observation range does not overlap\n"
			       "with the working data set"),
			     sinfo->varname);
	err = 1;
    }

    return err;
}

int check_db_import (SERIESINFO *sinfo, DATASET *dset)
{
    int err = 0;

    if (sinfo->pd < dset->pd) {
	if (sinfo->pd != 1 && sinfo->pd != 4 && 
	    dset->pd != 4 && dset->pd != 12) {
	    gretl_errmsg_sprintf(_("%s: can't handle conversion"),
				 sinfo->varname);
	    err = 1;
	} 
    }

    if (!err) {
	err = db_range_check(sinfo, dset);
    }

#if DB_DEBUG
    if (err) {
	fprintf(stderr, "check_db_import: err = %d\n", err);
	fprintf(stderr, "(dset->n = %d)\n", dset->n);
    }
#endif
    
    return err;
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

static int cli_add_db_data (double **dbZ, SERIESINFO *sinfo, 
			    DATASET *dset, CompactMethod method, 
			    int dbv, int *newdata, PRN *prn)
{
    double *xvec = NULL;
    int pad1 = 0, pad2 = 0;
    int t, start, stop;
    int free_xvec = 0;
    int new = (dbv == dset->v);

    if (dset->n == 0) {
	/* if the length of the dataset is not defined yet,
	   initialize using sinfo */
	init_datainfo_from_sinfo(dset, sinfo);
	dset->v = 0; /* trigger for creating data array below */
	*newdata = 1;
	if (dset->pd != 1 || strcmp(dset->stobs, "1")) { 
	    dset->structure = TIME_SERIES;
	}
    } else if (check_db_import(sinfo, dset)) {
	return 1;
    }

    if (dset->v == 0) {
	/* the data matrix is still empty */
	dset->v = 2;
	dbv = 1;
	if (start_new_Z(dset, 0)) {
	    gretl_errmsg_set(_("Out of memory!"));
	    return 1;
	}
    } else if (new && dataset_add_series(dset, 1)) {
	gretl_errmsg_set(_("Out of memory!"));
	return 1;
    }

#if DB_DEBUG
    fprintf(stderr, "dset->Z=%p\n", (void *) dset->Z);
    fprintf(stderr, "dset->n = %d, dset->v = %d, dbv = %d\n", 
	    dset->n, dset->v, dbv);
#endif

    /* is the frequency of the new var higher? */
    if (sinfo->pd > dset->pd) {
	if (dset->pd != 1 && dset->pd != 4 &&
	    sinfo->pd != 12) {
	    gretl_errmsg_set(_("Sorry, can't handle this conversion yet!"));
	    if (new) {
		dataset_drop_last_variables(dset, 1);
	    }
	    return 1;
	}
	if (method == COMPACT_NONE) {
	    method = COMPACT_AVG;
	    pprintf(prn, _("%s: using default compaction method: averaging"), 
		    sinfo->varname);
	    pputc(prn, '\n');
	}
	xvec = compact_db_series(dbZ[1], sinfo, dset->pd, method);
	if (xvec == NULL) {
	    gretl_errmsg_set(_("Out of memory!"));
	    if (new) {
		dataset_drop_last_variables(dset, 1);
	    }
	    return 1;
	}
	free_xvec = 1;
    } else {  
	/* series does not need compacting */
	xvec = dbZ[1];
    }

    /* common stuff for adding a var */
    strcpy(dset->varname[dbv], sinfo->varname);
    series_set_label(dset, dbv, sinfo->descrip);
    series_set_compact_method(dset, dbv, method);
    get_db_padding(sinfo, dset, &pad1, &pad2);

    if (pad1 > 0) {
#if DB_DEBUG
	fprintf(stderr, "Padding at start, %d obs\n", pad1);
#endif
	for (t=0; t<pad1; t++) {
	    dset->Z[dbv][t] = NADBL;
	}
	start = pad1;
    } else {
	start = 0;
    }

    if (pad2 > 0) {
#if DB_DEBUG
	fprintf(stderr, "Padding at end, %d obs\n", pad2);
#endif
	for (t=dset->n - 1; t>=dset->n - 1 - pad2; t--) {
	    dset->Z[dbv][t] = NADBL;
	}
	stop = dset->n - pad2;
    } else {
	stop = dset->n;
    }

    /* fill in actual data values */
#if DB_DEBUG
    fprintf(stderr, "Filling in values from %d to %d\n", start, stop - 1);
#endif
    for (t=start; t<stop; t++) {
	dset->Z[dbv][t] = xvec[t-pad1];
    }

    if (free_xvec) {
	free(xvec);
    }

    return 0;
}

/* compact an individual series, in the context of converting an
   entire working dataset to a lower frequency: used in all cases
   except conversion from daily to monthly
*/

static double *compact_series (const DATASET *dset, int i, int oldn,
			       int startskip, int min_startskip, int compfac,
			       CompactMethod method)
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
       end-of-period compaction is designed to allow for the
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
	eop_t = get_days_in_month(mon, yr, dset->pd) - 1;
    }

#if DMDEBUG
    fprintf(stderr, "starting: offset=%d, any_eop=%d, sop_t=%d, eop_t=%d\n", 
	    offset, any_eop, sop_t, eop_t);
#endif

    for (t=0; t<nm; t++) {
	/* loop across the months in the compacted data */
	int mdays = get_days_in_month(mon, yr, dset->pd);

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
		   int oldn, CompactMethod method, 
		   int *startskip, int *newn) 
{
    int ss, es, n;

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
    } else {
	ss = (compfac - (startmin % compfac) + 1) % compfac;
	es = endmin % compfac;
	n = (oldn - ss - es) / compfac;
    }

    *startskip = ss;
    *newn = n;
}

/* specific to compaction of daily time series */

static void 
get_daily_compact_params (CompactMethod default_method, 
			  int *any_eop, int *any_sop,
			  int *all_same,
			  const DATASET *dset)
{
    int i, n_not_eop = 0, n_not_sop = 0;

    *all_same = 1;
    *any_eop = (default_method == COMPACT_EOP)? 1 : 0;
    *any_sop = (default_method == COMPACT_SOP)? 1 : 0;

    for (i=1; i<dset->v; i++) {
	CompactMethod method = series_get_compact_method(dset, i);

	if (method != default_method && method != COMPACT_NONE) {
	    *all_same = 0;
	    if (method == COMPACT_EOP) {
		*any_eop = 1;
	    } else {
		n_not_eop++;
	    }
	    if (method == COMPACT_SOP) {
		*any_sop = 1;
	    } else {
		n_not_sop++;
	    }
	}
    }

    if (n_not_eop == dset->v - 1) {
	*any_eop = 0;
    }

    if (n_not_sop == dset->v - 1) {
	*any_sop = 0;
    }
}

/* specific to non-daily time series (monthly or quarterly) */

static void 
get_global_compact_params (int compfac, int startmin, int endmin,
			   CompactMethod default_method, 
			   int *min_startskip, int *max_n,
			   int *any_eop, int *all_same,
			   DATASET *dset)
{
    CompactMethod method;
    int i, startskip, n;
    int n_not_eop = 0;

    for (i=0; i<dset->v; i++) {
	if (i == 0) {
	    get_startskip_etc(compfac, startmin, endmin, dset->n, 
			      default_method, &startskip, &n);
	    if (default_method == COMPACT_EOP) {
		*any_eop = 1;
	    }
	} else {
	    method = series_get_compact_method(dset, i);
	    if (method != default_method && method != COMPACT_NONE) {
		get_startskip_etc(compfac, startmin, endmin, dset->n, 
				  method, &startskip, &n);
		*all_same = 0;
		if (method == COMPACT_EOP) {
		    *any_eop = 1;
		} else {
		    n_not_eop++;
		}
	    } 
	}
	if (startskip < *min_startskip) {
	    *min_startskip = startskip;
	}
	if (n > *max_n) {
	    *max_n = n;
	}
    }

    if (n_not_eop == dset->v - 1) {
	*any_eop = 0;
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

/* for daily data, figure the number of observations to
   be skipped at the start of each series 
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

/* for daily data, figure the number of valid monthly 
   observations that can be constructed by compaction
*/

static int get_n_ok_months (const DATASET *dset, 
			    CompactMethod default_method,
			    int *startyr, int *startmon,
			    int *endyr, int *endmon,
			    int *offset, int *p_any_eop)
{
    int y1, m1, d1, y2, m2, d2;
    int any_eop, any_sop, all_same;
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

    get_daily_compact_params(default_method, &any_eop, &any_sop,
			     &all_same, dset);

    *startyr = y1;
    *startmon = m1;
    *endyr = y2;
    *endmon = m2;

#if DMDEBUG
    fprintf(stderr, "get_n_ok_months: any_sop=%d, any_eop=%d, "
	    "all_same=%d\n", any_sop, any_eop, all_same);
    fprintf(stderr, "y1=%d m1=%d d1=%d; y2=%d m2=%d d2=%d\n", 
	    y1, m1, d1, y2, m2, d2);
#endif

    if (!day_starts_month(d1, m1, y1, dset->pd, &pad) && !any_eop) {
	if (*startmon == 12) {
	    *startmon = 1;
	    *startyr += 1;
	} else {
	    *startmon += 1;
	}
	skip = 1;
	nm--;
    }

    if (!day_ends_month(d2, m2, y2, dset->pd) && !any_sop) {
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
	*offset = get_daily_offset(dset, y1, m1, d1, skip, any_eop);
    }

    *p_any_eop = any_eop;

    return nm;
}

#define WEEKLY_DEBUG 0

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
	ntodate(obsstr, s, dset);
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
    int wcount = 0, mcount = 0;
    int monbak = 0;
    int t, err = 0;

    for (t=0; t<dset->n; t++) {
	ntodate(obsstr, t, dset);
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
	    wcount = 1;	
	} else if (mon != monbak) {
	    /* got a new month: report on previous one */
#if WEEKLY_DEBUG
	    fprintf(stderr, "month %d ('%d'), weekly obs = %d\n", 
		    mcount, monbak, wcount);
#endif
	    mcount++;
	    wcount = 1;
	} else {
	    /* continuation of current month */
	    wcount++;
	}
	monbak = mon;
    }

    if (err) {
	mcount = 0;
    } else {
	/* flush the last observation */
#if WEEKLY_DEBUG
	fprintf(stderr, "month %d ('%d'), weekly obs = %d\n", 
		mcount, monbak, wcount);
#endif
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

/* conversion to weekly using a "representative day", e.g. use
   each Wednesday value: @repday is 0-based on Sunday.
*/

static int daily_dataset_to_weekly (DATASET *dset, int repday)
{
    int y1, m1, d1;
    char obs[OBSLEN];
    double *x = NULL;
    double *tmp;
    int n = 0, n_ok = 0;
    int wday, ok;
    int i, t, err = 0;

    fprintf(stderr, "daily_dataset_to_weekly: repday = %d\n", repday);

    for (t=0; t<dset->n; t++) {
	ntodate(obs, t, dset);
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
	    ntodate(obs, t, dset);
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
	ntodate(dset->endobs, dset->t2, dset);

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
	if (method == COMPACT_NONE) {
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
    int dd = calendar_obs_number(dset->S[t], dset) -
	calendar_obs_number(dset->S[t-1], dset);

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
	ntodate(dset->endobs, dset->n - 1, dset);
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
    int nmiss = n_hidden_missing_obs(dset);
    int err = 0;

    fprintf(stderr, "n_hidden_missing_obs: nmiss = %d\n", nmiss);

    if (nmiss < 0) {
	err = 1;
    } else if (nmiss > 0) {
	err = insert_missing_hidden_obs(dset, nmiss);
    }

    return err;
}

/**
 * compact_data_set:
 * @dset: dataset struct.
 * @newpd: target data frequency.
 * @default_method: code for the default compaction method.
 * @monstart: FIXME add explanation.
 * @repday: "representative day" for conversion from daily
 * to weekly data (with method %COMPACT_WDAY only).
 * 
 * Compact the data set from higher to lower frequency.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int compact_data_set (DATASET *dset, int newpd,
		      CompactMethod default_method, 
		      int monstart, int repday)
{
    int newn, oldn = dset->n, oldpd = dset->pd;
    int compfac;
    int startmaj, startmin;
    int endmaj, endmin;
    int any_eop, all_same;
    int min_startskip = 0;
    char stobs[OBSLEN];
    int i, err = 0;

    gretl_error_clear();

    if (oldpd == 52) {
	return weekly_dataset_to_monthly(dset, default_method);
    }

    if (dated_daily_data(dset)) {
	/* allow for the possibility that the daily dataset
	   contains "hidden" or suppressed missing observations
	   (holidays are just skipped, not marked as NA)
	*/
	err = maybe_expand_daily_data(dset);
	if (err) {
	    gretl_errmsg_set("Error expanding daily data with missing observations");
	    return err;
	} else {
	    oldn = dset->n;
	}
    }

    if (newpd == 52 && oldpd >= 5 && oldpd <= 7 && 
	default_method == COMPACT_WDAY) {
	/* daily to weekly, using "representative day" */
	return daily_dataset_to_weekly(dset, repday);
    } else if (newpd == 12 && oldpd >= 5 && oldpd <= 7) {
	/* daily to monthly: special */
	return daily_dataset_to_monthly(dset, default_method);
    } else if (oldpd >= 5 && oldpd <= 7) {
	/* daily to weekly */
	compfac = oldpd;
	if (dated_daily_data(dset)) {
	    startmin = weekday_from_date(dset->stobs);
	    if (oldpd == 7) {
		if (monstart) {
		    if (startmin == 0) startmin = 7;
		} else {
		    startmin++;
		}
	    }
	} else {
	    startmin = 1;
	}
    } else if (oldpd == 24 && newpd >= 5 && newpd <= 7) {
	/* hourly to daily */
	compfac = 24;
	if (!get_obs_maj_min(dset->stobs, &startmaj, &startmin)) {
	    return 1;
	}
    } else {
	compfac = oldpd / newpd;
	/* get starting obs major and minor components */
	if (!get_obs_maj_min(dset->stobs, &startmaj, &startmin)) {
	    return 1;
	}
	/* get ending obs major and minor components */
	if (!get_obs_maj_min(dset->endobs, &endmaj, &endmin)) {
	    return 1;
	} 
    }

    min_startskip = oldpd;
    newn = 0;
    any_eop = 0;
    all_same = 1;
    get_global_compact_params(compfac, startmin, endmin, default_method,
			      &min_startskip, &newn, &any_eop, &all_same, 
			      dset);
    if (newn == 0) {
	gretl_errmsg_set(_("Compacted dataset would be empty"));
	return 1;
    }    

    if (newpd == 1) {
	if (min_startskip > 0 && !any_eop) { 
	    startmaj++;
	}
	sprintf(stobs, "%d", startmaj);
    } else if (newpd == 52) {
	if (oldpd >= 5 && oldpd <= 7 && dset->S != NULL) {
	    strcpy(stobs, dset->S[min_startskip]);
	} else {
	    strcpy(stobs, "1");
	}
    } else {
	int m0 = startmin + min_startskip;
	int minor = m0 / compfac + (m0 % compfac > 0);

	if (minor > newpd) {
	    startmaj++;
	    minor -= newpd;
	}
	format_obs(stobs, startmaj, minor, newpd);
    }

    /* revise datainfo members */
    strcpy(dset->stobs, stobs);
    dset->pd = newpd;
    dset->n = newn;
    dset->sd0 = get_date_x(dset->pd, dset->stobs);
    dset->t1 = 0;
    dset->t2 = dset->n - 1;
    ntodate(dset->endobs, dset->t2, dset);
    
    if (oldpd >= 5 && oldpd <= 7 && dset->markers) {
	/* remove any daily date strings; revise endobs */
	dataset_destroy_obs_markers(dset);
	ntodate(dset->endobs, dset->t2, dset);
    }

    err = shorten_the_constant(dset->Z, dset->n);

    /* compact the individual data series */
    for (i=1; i<dset->v && !err; i++) {
	CompactMethod this_method = default_method;
	int startskip = min_startskip;
	double *x;

	if (!all_same) {
	    CompactMethod m_i = series_get_compact_method(dset, i);

	    if (m_i != COMPACT_NONE) {
		this_method = m_i;
	    }

	    startskip = compfac - (startmin % compfac) + 1;
	    startskip = startskip % compfac;
	    if (this_method == COMPACT_EOP) {
		if (startskip > 0) {
		    startskip--;
		} else {
		    startskip = compfac - 1;
		}
	    }
	}

	x = compact_series(dset, i, oldn, startskip, min_startskip, 
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

static gretl_matrix *interpol_expand_dataset (const DATASET *dset, 
					      int newpd, int *err)
{
    gretl_matrix *Y0, *Y1 = NULL;
    int *list;

    list = gretl_consecutive_list_new(1, dset->v - 1);
    if (list == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    Y0 = gretl_matrix_data_subset(list, dset, dset->t1, dset->t2,
				  M_MISSING_ERROR, err);

    if (!*err) {
	int f = newpd / dset->pd;

	Y1 = matrix_chowlin(Y0, NULL, f, err);
	gretl_matrix_free(Y0);
    }

    free(list);

    return Y1;
}

/**
 * expand_data_set:
 * @dset: dataset struct.
 * @newpd: target data frequency.
 * @interpol: use interpolation (0/1).
 * 
 * Expand the data set from lower to higher frequency: an "expert"
 * option.  This is supported at present only for expansion from
 * annual to quarterly or monthly, or from quarterly to monthly.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int expand_data_set (DATASET *dset, int newpd, int interpol)
{
    char stobs[12];
    int oldn = dset->n;
    int oldpd = dset->pd;
    int t1 = dset->t1;
    int t2 = dset->t2;
    int mult, newn, nadd;
    gretl_matrix *X = NULL;
    double *x = NULL;
    int i, j, s, t;
    int err = 0;

    if (oldpd != 1 && oldpd != 4) {
	return E_PDWRONG;
    } else if (oldpd == 1 && newpd != 4 && newpd != 12) {
	return E_DATA;
    } else if (oldpd == 4 && newpd != 12) {
	return E_DATA;
    } else if (oldpd == 1 && newpd == 12 && interpol) {
	return E_DATA;
    }

    if (interpol) {
	X = interpol_expand_dataset(dset, newpd, &err);
    } else {
	x = malloc(oldn * sizeof *x);
	if (x == NULL) {
	    err = E_ALLOC;
	}
    }

    if (err) {
	return err;
    }

    mult = newpd / oldpd;
    newn = mult * dset->n;
    nadd = newn - oldn;

    err = dataset_add_observations(dset, nadd, OPT_NONE);
    if (err) {
	goto bailout;
    }

    if (interpol) {
	for (i=1; i<dset->v; i++) {
	    for (t=0; t<newn; t++) {
		dset->Z[i][t] = gretl_matrix_get(X, t, i-1);
	    }
	}
    } else {
	for (i=1; i<dset->v; i++) {
	    for (t=0; t<oldn; t++) {
		x[t] = dset->Z[i][t];
	    }
	    s = 0;
	    for (t=0; t<oldn; t++) {
		for (j=0; j<mult; j++) {
		    dset->Z[i][s++] = x[t];
		}
	    }
	}
    }

    if (dset->pd == 1) {
	strcpy(stobs, dset->stobs);
	if (newpd == 4) {
	    strcat(stobs, ":1");
	} else {
	    strcat(stobs, ":01");
	}
    } else {
	int yr, qtr, mo;

	sscanf(dset->stobs, "%d:%d", &yr, &qtr);
	mo = (qtr - 1) * 3 + 1;
	sprintf(stobs, "%d:%02d", yr, mo);
    }

    if (dset->t1 > 0) {
	dset->t1 *= mult;
    }

    if (dset->t2 < oldn - 1) {
	dset->t2 = dset->t1 + (t2 - t1 + 1) * mult - 1;
    }    

    strcpy(dset->stobs, stobs);
    dset->pd = newpd;
    dset->sd0 = get_date_x(dset->pd, dset->stobs);
    ntodate(dset->endobs, dset->n - 1, dset);

    if (dset->markers) {
	dataset_destroy_obs_markers(dset);
    }

 bailout:
    
    free(x);
    gretl_matrix_free(X);

    return err;
}
