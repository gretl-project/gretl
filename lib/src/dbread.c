/*
 *   Copyright (c) by Allin Cottrell
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

/* dbread.c for gretl */

#include <glib.h>

#include "libgretl.h"

/* #define DB_DEBUG */

#define RECNUM long
#define NAMELENGTH 16
#define RATSCOMMENTLENGTH 80
#define RATSCOMMENTS 2

typedef struct {
    long daynumber;                /* Number of days from 1-1-90
				      to year, month, day */
    short panel;                   /* 1 for panel set, 2 for intraday
				      date set , 0 o.w. */
#define LINEAR   0                 /* Single time direction */
#define PANEL    1                 /* panel:period */    
#define INTRADAY 2                 /* date:intraday period */
    long panelrecord;              /* Size of panel or 
				      number of periods per day */
    short dclass;                  /* See definitions below */
#define UNDATEDCLASS   0           /* No time series properties */
#define IRREGULARCLASS 1           /* Time series (irregular) */
#define PERYEARCLASS   2           /* x periods / year */
#define PERWEEKCLASS   3           /* x periods / week */
#define DAILYCLASS     4           /* x days / period */
    long info;                     /* Number of periods per year or
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
    long datapoints;               /* Number of data points */
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

static char db_name[MAXLEN];
static int db_type;

static int cli_add_db_data (double **dbZ, SERIESINFO *sinfo, 
			    double ***pZ, DATAINFO *pdinfo,
			    int compact_method, int dbv);

/* ........................................................... */

int get_native_db_data (const char *dbbase, SERIESINFO *sinfo, 
			double **Z)
{
    char dbbin[MAXLEN], numstr[16];
    FILE *fp;
    int t, n = sinfo->nobs;
    dbnumber val;

    strcpy(dbbin, dbbase);
    if (strstr(dbbin, ".bin") == NULL) {
	strcat(dbbin, ".bin");
    }

    fp = fopen(dbbin, "rb");
    if (fp == NULL) return 1;
    
    fseek(fp, (long) sinfo->offset, SEEK_SET);
    for (t=0; t<n; t++) {
	fread(&val, sizeof val, 1, fp);
	sprintf(numstr, "%g", val);
	Z[1][t] = atof(numstr);
    }

    fclose(fp);

    return 0;
}

static void get_native_series_comment (SERIESINFO *sinfo, const char *s)
{
    size_t n = strlen(sinfo->varname);
    const char *p = s + n + 1;
    int i;

    while (*p) {
	if (isspace(*p)) p++;
	else break;
    }

    *sinfo->descrip = 0;
    strncat(sinfo->descrip, p, MAXLABEL - 1);

    n = strlen(sinfo->descrip) - 1;
    
    for (i=n; i>0; i--) {
	if (isspace(sinfo->descrip[i])) sinfo->descrip[i] = 0;
	else if (sinfo->descrip[i] == '\r' || sinfo->descrip[i] == '\n') {
	    sinfo->descrip[i] = 0;
	}
	else break;
    }
}

static int get_native_series_pd (SERIESINFO *sinfo, char pdc)
{
    sinfo->pd = 1;
    sinfo->undated = 0;

    if (pdc == 'M') sinfo->pd = 12;
    else if (pdc == 'Q') sinfo->pd = 4;
    else if (pdc == 'B') sinfo->pd = 5;
    else if (pdc == 'D') sinfo->pd = 7;
    else if (pdc == 'U') sinfo->undated = 1;
    else return 1;

    return 0;
}

static int get_native_series_obs (SERIESINFO *sinfo, 
				  const char *stobs,
				  const char *endobs)
{
    if (strchr(stobs, '/')) { /* daily data */
	const char *q = stobs;
	const char *p = strchr(stobs, '/');

	if (p - q == 4) strcpy(sinfo->stobs, q + 2);
	q = endobs;
	p = strchr(endobs, '/');
	if (p && p - q == 4) strcpy(sinfo->endobs, q + 2);
    } else {
	*sinfo->stobs = 0;
	*sinfo->endobs = 0;
	strncat(sinfo->stobs, stobs, 8);
	strncat(sinfo->endobs, endobs, 8);
    }

    return 0;
}

static int 
get_native_series_info (const char *series, SERIESINFO *sinfo)
{
    FILE *fp;
    char dbidx[MAXLEN];
    char sername[9];
    char line1[256], line2[72];
    char stobs[11], endobs[11];
    char *p, pdc;
    int offset = 0;
    int gotit = 0, err = DB_OK;
    int n;

    strcpy(dbidx, db_name);
    p = strstr(dbidx, ".bin");
    if (p != NULL) {
	strcpy(p, ".idx");
    } else {
	strcat(dbidx, ".idx");
    }

    fp = fopen(dbidx, "r");
    if (fp == NULL) {
	strcpy(gretl_errmsg, _("Couldn't open database index file"));
	return DB_NOT_FOUND;
    }

    while (!gotit) {

	if (fgets(line1, 255, fp) == NULL) break;

	if (*line1 == '#') continue;
	line1[255] = 0;

	if (sscanf(line1, "%8s", sername) != 1) break;
	sername[8] = 0;

	if (!strcmp(series, sername)) {
	    gotit = 1;
	    strcpy(sinfo->varname, sername);
	}

	fgets(line2, 71, fp);
	line2[71] = 0;

	if (gotit) {
	    get_native_series_comment(sinfo, line1);
	    if (sscanf(line2, "%c %10s %*s %10s %*s %*s %d", 
		       &pdc, stobs, endobs, &(sinfo->nobs)) != 4) {
		strcpy(gretl_errmsg,
		       _("Failed to parse series information"));
		err = DB_PARSE_ERROR;
	    } else {
		get_native_series_pd(sinfo, pdc);
		get_native_series_obs(sinfo, stobs, endobs);
		sinfo->offset = offset;
	    }
	} else {
	    if (sscanf(line2, "%*c %*s %*s %*s %*s %*s %d", &n) != 1) {
		strcpy(gretl_errmsg,
		       _("Failed to parse series information"));
		err = DB_PARSE_ERROR;
	    } else {
		offset += n * sizeof(dbnumber);
	    }
	}
    }

    fclose(fp);

    if (!gotit) {
	sprintf(gretl_errmsg, _("Series not found, '%s'"), series);
	err = DB_NO_SUCH_SERIES;
    }

    return err;
}

/* ........................................................... */

static int get_endobs (char *datestr, int startyr, int startfrac, 
		       int pd, int n)
/* Figure the ending observation date of a series */
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
    } else if (pd == 12) {
	sprintf(datestr, "%d.%02d", endyr, endfrac);
    }

    return 0;
}

static int dinfo_sanity_check (const DATEINFO *dinfo)
{
    if (dinfo->info < 0 || dinfo->info > 365 ||
	dinfo->year < 0 || dinfo->year > 3000 ||
	dinfo->month < 0 || dinfo->month > 12 ||
	dinfo->day < 0 || dinfo->day > 31) {
	strcpy(gretl_errmsg, _("This is not a valid RATS 4.0 database"));
	return 1;
    }

    return 0;
}

static int dinfo_to_sinfo (const DATEINFO *dinfo, SERIESINFO *sinfo,
			   const char *varname, const char *comment,
			   int n, int offset)
{
    int startfrac = 0;
    char pdstr[3]; 
    int err = 0;

    if (dinfo_sanity_check(dinfo)) return 1;

    sprintf(sinfo->stobs, "%d", dinfo->year); 
    *pdstr = 0;

    if (dinfo->info == 4) {
	sprintf(pdstr, ".%d", dinfo->month);
	if (dinfo->month == 1) startfrac = 1;
	else if (dinfo->month > 1 && dinfo->month <= 4) startfrac = 2;
	else if (dinfo->month > 4 && dinfo->month <= 7) startfrac = 3;
	else startfrac = 4;
    }
    else if (dinfo->info == 12) {
	sprintf(pdstr, ".%02d", dinfo->month);
	startfrac = dinfo->month;
    }
    else if (dinfo->info == 1) {
	startfrac = 0;
    } else {
	fprintf(stderr, _("frequency (%d) does not make seem to make sense"),
		(int) dinfo->info);
	err = 1;
    }   

    if (*pdstr) strcat(sinfo->stobs, pdstr);
    get_endobs(sinfo->endobs, dinfo->year, startfrac, dinfo->info, n);

    sinfo->pd = dinfo->info;
    sinfo->nobs = n;
    *sinfo->varname = 0;
    strncat(sinfo->varname, varname, 8);
    *sinfo->descrip = 0;
    strncat(sinfo->descrip, comment, MAXLABEL - 1);

    sinfo->offset = offset;

    return err;
}

static int dinfo_to_tbl_row (const DATEINFO *dinfo, db_table_row *row,
			     const char *varname, const char *comment,
			     int n)
{
    char pd = 0, pdstr[3], endobs[9];
    int startfrac = 0;
    int err = 0;

    if (dinfo_sanity_check(dinfo)) return 1;

    *pdstr = 0;

    if (dinfo->info == 4) {
	pd = 'Q';
	sprintf(pdstr, ".%d", dinfo->month);
	if (dinfo->month == 1) startfrac = 1;
	else if (dinfo->month > 1 && dinfo->month <= 4) startfrac = 2;
	else if (dinfo->month > 4 && dinfo->month <= 7) startfrac = 3;
	else startfrac = 4;
    }
    else if (dinfo->info == 12) {
	pd = 'M';
	sprintf(pdstr, ".%02d", dinfo->month);
	startfrac = dinfo->month;
    }
    else if (dinfo->info == 1) {
	pd = 'A';
	startfrac = 0;
    } else {
	fprintf(stderr, _("frequency (%d) does not make seem to make sense"),
		(int) dinfo->info);
	err = 1;
    }

    if (!err) {
	get_endobs(endobs, dinfo->year, startfrac, dinfo->info, n);
	row->varname = g_strdup(varname);
	row->comment = g_strdup(comment);
	row->obsinfo = g_strdup_printf("%c  %d%s - %s  n = %d", pd, 
				       (int) dinfo->year, pdstr, endobs, n);
    }

    return err;
}

static RECNUM read_rats_directory (FILE *fp, db_table_row *row,
				   const char *series_name,
				   SERIESINFO *sinfo) 
{
    RATSDirect rdir;
    DATEINFO dinfo;
    int err = 0;

    fread(&rdir.back_point, sizeof(RECNUM), 1, fp);
    fread(&rdir.forward_point, sizeof(RECNUM), 1, fp);
    fseek(fp, 4L, SEEK_CUR); /* skip two shorts */
    fread(&rdir.first_data, sizeof(RECNUM), 1, fp);
    fread(rdir.series_name, 16, 1, fp);  
    rdir.series_name[8] = '\0';
    chopstr(rdir.series_name);

    if (series_name != NULL && strcmp(series_name, rdir.series_name)) {
	/* specific series not found yet: keep going */
	return rdir.forward_point;
    }

    /* Now the dateinfo: we can't read this in one go either :-( */
    fseek(fp, 12, SEEK_CUR); /* skip long, short, long, short */
    fread(&dinfo.info, sizeof(long), 1, fp);
    fread(&dinfo.digits, sizeof(short), 1, fp);
    fread(&dinfo.year, sizeof(short), 1, fp);
    fread(&dinfo.month, sizeof(short), 1, fp);
    fread(&dinfo.day, sizeof(short), 1, fp);

    fread(&rdir.datapoints, sizeof(long), 1, fp);
    fseek(fp, sizeof(short) * 4L, SEEK_CUR);  /* skip 4 shorts */

    fread(&rdir.comment_lines, sizeof(short), 1, fp);
    fseek(fp, 1L, SEEK_CUR); /* skip one char */

    fread(rdir.comments[0], 80, 1, fp);
    rdir.comments[0][79] = '\0';
    chopstr(rdir.comments[0]);

    fread(rdir.comments[1], 80, 1, fp);
    rdir.comments[1][79] = '\0';
    chopstr(rdir.comments[1]);

    if (sinfo != NULL) {
	err = dinfo_to_sinfo(&dinfo, sinfo, rdir.series_name, rdir.comments[0],
			     rdir.datapoints, rdir.first_data);
    } else if (row != NULL) {
	err = dinfo_to_tbl_row(&dinfo, row, rdir.series_name, rdir.comments[0],
			       rdir.datapoints);
    } else {
	err = 1;
    }

    if (err) {
	return -999;
    } else {
	return rdir.forward_point;
    }
}

#define DB_INIT_ROWS 32

static db_table *db_table_new (void)
{
    db_table *tbl;

    tbl = malloc(sizeof *tbl);
    if (tbl == NULL) return NULL;

    tbl->rows = malloc(DB_INIT_ROWS * sizeof *tbl->rows);

    if (tbl->rows == NULL) {
	free(tbl);
	return NULL;
    }

    tbl->nrows = 0;

    return tbl;
}

static int db_table_expand (db_table *tbl)
{
    db_table_row *rows;
    int newsz;

    newsz = (tbl->nrows / DB_INIT_ROWS) + 1;
    newsz *= DB_INIT_ROWS;

    rows = realloc(tbl->rows, newsz * sizeof *rows);
    if (rows == NULL) {
	free(tbl->rows);
	tbl->rows = NULL;
	return 1;
    }

    tbl->rows = rows;

    return 0;
}

static void db_table_free (db_table *tbl)
{
    int i;

    for (i=0; i<tbl->nrows; i++) {
	free(tbl->rows[i].varname);
	free(tbl->rows[i].comment);
	free(tbl->rows[i].obsinfo);
    }

    free(tbl->rows);
    free(tbl);
}

/**
 * read_rats_db:
 * @fp: pre-opened stream (caller to close it)
 *
 * Read the series info from a RATS 4.0 database.
 * 
 * Returns: pointer to a #db_table containing the series info,
 * or NULL in case of failure.
 *
 */

db_table *read_rats_db (FILE *fp) 
/* read the base block at offset 0 in the data file,
   and recurse through the directory entries */
{
    db_table *tbl;
    long forward;
    int i, err = 0;

    *gretl_errmsg = 0;
    
    /* get into position */
    fseek(fp, 30L, SEEK_SET); /* skip unneeded fields */
    fread(&forward, sizeof forward, 1, fp);
    fseek(fp, 4L, SEEK_CUR);

    /* basic check */
    if (forward <= 0) {
	strcpy(gretl_errmsg, _("This is not a valid RATS 4.0 database"));
	return NULL;
    }

    /* allocate table for series rows */
    tbl = db_table_new();
    if (tbl == NULL) {
	strcpy(gretl_errmsg, _("Out of memory!"));
	return NULL;
    }
    
    /* Go find the series */
    i = 0;
    while (forward && !err) {
	tbl->nrows += 1;
	if (tbl->nrows > 0 && tbl->nrows % DB_INIT_ROWS == 0) {
	    err = db_table_expand(tbl);
	    if (err) {
		strcpy(gretl_errmsg, _("Out of memory!"));
	    }
	}
	if (!err) {
	    fseek(fp, (forward - 1) * 256L, SEEK_SET);
	    forward = read_rats_directory(fp, &tbl->rows[i++], NULL, NULL);
	    if (forward == -999) err = 1;
	}
    }

    if (err) {
	db_table_free(tbl);
	return NULL;
    }

    return tbl;
}

/* ........................................................... */

static int find_rats_dir_by_number (FILE *fp, int first_dir, 
				    int series_number)
{
    long forward;
    int count = 1;

    forward = first_dir;

    while (forward && count < series_number) {
	fseek(fp, (forward - 1) * 256L, SEEK_SET);
	fseek(fp, 4L, SEEK_CUR);
	fread(&forward, 4L, 1, fp);
	count++;
    }

    return (int) forward;
}

/* ........................................................... */

static int get_rats_series_offset_by_number (FILE *fp, 
					     int series_number)
{
    long num_series, first_dir;
    int offset;

    fseek(fp, 6L, SEEK_SET);
    fread(&num_series, sizeof num_series, 1, fp);
    if (series_number > num_series) return -1;

    fseek(fp, sizeof(long) * 5L, SEEK_CUR);  
    fread(&first_dir, sizeof first_dir, 1, fp);

    offset = find_rats_dir_by_number(fp, first_dir, series_number); 

    return offset;
}

/* ........................................................... */

static int get_rats_series (int offset, SERIESINFO *sinfo, FILE *fp, 
			    double **Z)
/* retrieve the actual data values from the data blocks */
{
    RATSData rdata;
    char numstr[16];
    int miss = 0, i, t = 0;
    double val;
    
    rdata.forward_point = offset;

    while (rdata.forward_point) {
	fseek(fp, (rdata.forward_point - 1) * 256L, SEEK_SET);
	/* the RATSData struct is actually 256 bytes.  Yay! */
	fread(&rdata, sizeof(RATSData), 1, fp);
	for (i=0; i<31 && t<sinfo->nobs; i++) {
	    sprintf(numstr, "%g", rdata.data[i]);
	    val = atof(numstr);
	    if (isnan(val)) {
		val = NADBL;
		miss = 1;
	    }
	    Z[1][t] = val;
	    t++;
	}
    }

    return miss;
}

/**
 * get_rats_data_by_series_number:
 * @fname: name of RATS 4.0 database to read from
 * @series_number: number of the series within the database
 * @sinfo: holds info on the given series (input)
 * @pZ: pointer to data matrix
 *
 * Read the actual data values for a series from a RATS database.
 * 
 * Returns: DB_OK on successful completion, DB_NOT_FOUND if
 * the data could not be read, and DB_MISSING_DATA if the
 * data were found but there were some missing values.
 *
 */

int get_rats_data_by_series_number (const char *fname, 
				    int series_number,
				    SERIESINFO *sinfo, 
				    double **Z)
     /* series are numbered from 1 for this function */
{
    FILE *fp;
    int offset;
    long first_data;
    int ret = DB_OK;

    fp = fopen(fname, "rb");
    if (fp == NULL) return DB_NOT_FOUND;
    
    offset = get_rats_series_offset_by_number(fp, series_number);
    if (offset < 0) return DB_NOT_FOUND;

    fprintf(stderr, "get_rats_data_by_series_number: offset=%d\n", offset);
    
    fseek(fp, (offset - 1) * 256 + 12, SEEK_SET); 
    fread(&first_data, sizeof(RECNUM), 1, fp);
    if (get_rats_series(first_data, sinfo, fp, Z)) {
	ret = DB_MISSING_DATA;
    }

    fclose(fp);

    return ret;
}

static long 
get_rats_series_offset_by_name (FILE *fp,
				const char *series_name,
				SERIESINFO *sinfo)
{
    long forward;

    *gretl_errmsg = 0;
    
    /* get into position */
    fseek(fp, 30L, SEEK_SET); 
    fread(&forward, sizeof forward, 1, fp);
    fseek(fp, 4L, SEEK_CUR);

    /* basic check */
    if (forward <= 0) {
	strcpy(gretl_errmsg, _("This is not a valid RATS 4.0 database"));
	return -1;
    }

    sinfo->offset = 0;

    /* Go find the series */
    while (forward) {
	fseek(fp, (forward - 1) * 256L, SEEK_SET);
	forward = read_rats_directory(fp, NULL, series_name, sinfo);
	if (forward == -999) sinfo->offset = -1;
	if (sinfo->offset != 0) break;
    }

    return sinfo->offset;
}

static int 
get_rats_series_info_by_name (const char *series_name,
			      SERIESINFO *sinfo)
{
    FILE *fp;
    int offset;
    int ret = DB_OK;

    fp = fopen(db_name, "rb");
    if (fp == NULL) return DB_NOT_FOUND;
    
    offset = get_rats_series_offset_by_name(fp, series_name, sinfo);
    fclose(fp);

    if (offset <= 0) return DB_NOT_FOUND;

#ifdef DB_DEBUG
    fprintf(stderr, "get_rats_series_info_by_name: offset=%d\n", offset);
    fprintf(stderr, " pd = %d, nobs = %d\n", sinfo->pd, sinfo->nobs);
#endif

    return ret;
}

/* ........................................................... */

static int 
get_rats_data_by_offset (const char *fname, 
			 SERIESINFO *sinfo,
			 double **Z)
{
    FILE *fp;
    int ret = DB_OK;

    fp = fopen(fname, "rb");
    if (fp == NULL) return DB_NOT_FOUND;

    if (get_rats_series(sinfo->offset, sinfo, fp, Z)) {
	ret = DB_MISSING_DATA;
    }

    fclose(fp);

    return ret;
}

/* ........................................................... */

#ifdef LIMIT_DIGITS
static int get_places (double x)
{
    char numstr[16];
    int i, n, p = 0;

    sprintf(numstr, "%.3f", x);
    n = strlen(numstr);
    for (i=n-3; i<n; i++) {
	if (numstr[i] != '0') p++;
	else break;
    }

    return p;
}
#endif

int mon_to_quart (double **pq, double *mvec, SERIESINFO *sinfo,
		  int method)
{
    int t, p, m0, q0, y0, skip = 0, endskip, goodobs;
    float q;
    double val = 0.;
#ifdef LIMIT_DIGITS
    int pmax = 0;
    char numstr[16];

    /* record the precision of the original data */
    for (t=0; t<sinfo->nobs; t++) {
	p = get_places(mvec[t]);
	if (p > pmax) pmax = p;
    }
#endif

    /* figure the quarterly dates */
    y0 = atoi(sinfo->stobs);
    m0 = atoi(sinfo->stobs + 5);
    q = 1.0 + m0/3.;
    q0 = q + .5;
    skip = ((q0 - 1) * 3) + 1 - m0;
    if (q0 == 5) {
	y0++;
	q0 = 1;
    }

    endskip = (sinfo->nobs - skip) % 3;
    goodobs = (sinfo->nobs - skip - endskip) / 3;
    sinfo->nobs = goodobs;
    sprintf(sinfo->stobs, "%d.%d", y0, q0);

#ifdef DB_DEBUG
    fprintf(stderr, "startskip = %d\n", skip);
    fprintf(stderr, "endskip = %d\n", endskip);
    fprintf(stderr, "goodobs = %d\n", goodobs);
    fprintf(stderr, "starting date = %s\n", sinfo->stobs);
#endif

    *pq = malloc(goodobs * sizeof **pq);
    if (*pq == NULL) return 1;

    for (t=0; t<goodobs; t++) {
	p = (t + 1) * 3;
	if (method == COMPACT_AVG) 
	    val = (mvec[p-3+skip] + mvec[p-2+skip] + mvec[p-1+skip]) / 3.0;
	else if (method == COMPACT_SUM)
	    val = mvec[p-3+skip] + mvec[p-2+skip] + mvec[p-1+skip];
	else if (method == COMPACT_EOP)
	    val = mvec[p-1+skip];
	else if (method == COMPACT_SOP)
	    val = mvec[p-3+skip];
	/* do we want to limit the precision of the compacted
	   data to that of the original data? */
#ifdef LIMIT_DIGITS
	sprintf(numstr, "%.*f", pmax, val);
	(*pq)[t] = atof(numstr);
#else
	(*pq)[t] = val;
#endif
    }

    sinfo->pd = 4;

    return 0;
}

/* ........................................................... */

int to_annual (double **pq, double *mvec, SERIESINFO *sinfo,
	       int method)
{
    int i, t, p, p0, y0, skip = 0, endskip, goodobs;
    int pd = sinfo->pd;
    double val;
#ifdef LIMIT_DIGITS
    int pmax = 0;
    char numstr[16];

    /* record the precision of the original data */
    for (t=0; t<sinfo->nobs; t++) {
	p = get_places(mvec[t]);
	if (p > pmax) pmax = p;
    }
#endif

    /* figure the annual dates */
    y0 = atoi(sinfo->stobs);
    p0 = atoi(sinfo->stobs + 5);
    if (p0 != 1) {
	++y0;
	skip = pd - (p0 + 1);
    }

    fprintf(stderr, "startskip = %d\n", skip);
    endskip = (sinfo->nobs - skip) % pd;
    fprintf(stderr, "endskip = %d\n", endskip);
    goodobs = (sinfo->nobs - skip - endskip) / pd;
    fprintf(stderr, "goodobs = %d\n", goodobs);
    sinfo->nobs = goodobs;
    sprintf(sinfo->stobs, "%d", y0);
    fprintf(stderr, "starting date = %s\n", sinfo->stobs);

    *pq = malloc(goodobs * sizeof **pq);
    if (*pq == NULL) return 1;

    for (t=0; t<goodobs; t++) {
	p = (t + 1) * pd;
	val = 0.;

	if (method == COMPACT_AVG || method == COMPACT_SUM) { 
	    for (i=1; i<=pd; i++) val += mvec[p-i+skip];
	    if (method == COMPACT_AVG) {
		val /= (double) pd;
	    }
	}
	else if (method == COMPACT_EOP) {
	    val = mvec[p-1+skip];
	}
	else if (method == COMPACT_SOP) {
	    val = mvec[p-pd+skip];
	}

#ifdef LIMIT_DIGITS
	sprintf(numstr, "%.*f", pmax, val);
	(*pq)[t] = atof(numstr);
#else
	(*pq)[t] = val;
#endif
    }

    sinfo->pd = 1;

    return 0;
}

int set_db_name (const char *fname, int filetype, const PATHS *ppaths, 
		 PRN *prn)
{
    FILE *fp;
    int err = 0;

    *db_name = 0;
    strncat(db_name, fname, MAXLEN - 1);

    fp = fopen(db_name, "rb");

    if (fp == NULL) {
	/* try looking a bit more */
	if (filetype == GRETL_NATIVE_DB && 
	    strstr(db_name, ppaths->binbase) == NULL) {
	    strcpy(db_name, ppaths->binbase);
	    strcat(db_name, fname);
	}
	else if (filetype == GRETL_RATS_DB && 
		 strstr(db_name, ppaths->ratsbase) == NULL) {
	    strcpy(db_name, ppaths->ratsbase);
	    strcat(db_name, fname);
	}
	fp = fopen(db_name, "rb");
    }
	    
    if (fp == NULL) {
	*db_name = 0;
	pprintf(prn, _("Couldn't open %s\n"), fname);
	err = 1;
    } else {
	fclose(fp);
	db_type = filetype;
	pprintf(prn, "%s\n", db_name);
    }

    return err;
}

int db_set_sample (const char *line, DATAINFO *pdinfo)
{
    char cmd[5], start[9], stop[9];
    int t1 = 0, t2 = 0;

    if (sscanf(line, "%4s %8s %8s", cmd, start, stop) != 3) {
	sprintf(gretl_errmsg, _("error reading smpl line"));
	return 1;
    }

    if (strcmp(start, ";")) {
	t1 = dateton(start, pdinfo);
	if (t1 < 0 || strlen(gretl_errmsg)) {
	    return 1;
	}
    }

    t2 = dateton(stop, pdinfo);
    if (strlen(gretl_errmsg)) return 1;

    if (t1 > t2) {
	sprintf(gretl_errmsg, _("Invalid null sample"));
	return 1;
    }

    pdinfo->t1 = t1;
    pdinfo->t2 = t2;
    pdinfo->n = t2 - t1 + 1;
    strcpy(pdinfo->endobs, stop);

#ifdef DB_DEBUG
    fprintf(stderr, "db_set_sample: t1=%d, t2=%d, stobs='%s', endobs='%s' "
	    "sd0 = %g, n = %d\n", 
	    pdinfo->t1, pdinfo->t2, 
	    pdinfo->stobs, pdinfo->endobs,
	    pdinfo->sd0, pdinfo->n);
#endif

    return 0;
}

static const char *
get_word_and_advance (const char *s, char *word, size_t maxlen)
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

static const char *
get_compact_method_and_advance (const char *s, int *method)
{
    const char *p;

    *method = COMPACT_NONE;

    if ((p = strstr(s, "(compact"))) {
	char comp[8];
	int i;

	p += 8;
	i = 0;
	while (*p && *p != ')' && i < 7) {
	    if (!isspace(*p) && *p != '=') {
		comp[i++] = *p;
	    }
	    p++;
	}
	comp[i] = 0;

	if (!strcmp(comp, "average")) *method = COMPACT_AVG;
	else if (!strcmp(comp, "sum")) *method = COMPACT_SUM;
	else if (!strcmp(comp, "first")) *method = COMPACT_SOP;
	else if (!strcmp(comp, "last")) *method = COMPACT_EOP;

	p = strchr(p, ')');
	if (p != NULL) p++;
    } else {
	/* no compaction method given */
	if ((p = strstr(s, "data "))) p += 5;
    }

    return p;
}

static double **new_dbZ (int n)
{
    double **Z;

    Z = malloc(2 * sizeof *Z);
    if (Z == NULL) return NULL;

    Z[0] = NULL;
    Z[1] = malloc(n * sizeof **Z);

    if (Z[1] == NULL) {
	free(Z);
	return NULL;
    }

    return Z;
}

/* main function for getting a series out of a database, using the
   command-line client or in script or console mode
*/

int db_get_series (const char *line, double ***pZ, DATAINFO *pdinfo, 
		   PRN *prn)
{
    char series[16];
    int comp_method;
    SERIESINFO sinfo;
    double **dbZ;
    int err = 0;

    if (*db_name == 0) {
	strcpy(gretl_errmsg, _("No database has been opened"));
	return 1;
    }   

    line = get_compact_method_and_advance(line, &comp_method);

    /* now loop over variable names given on the line */

    while ((line = get_word_and_advance(line, series, 8))) {
	int v, this_var_method = comp_method; 

	/* see if the series is already in the dataset */
	v = varindex(pdinfo, series);
	fprintf(stderr, "db_get_series: pdinfo->v = %d, v = %d\n", pdinfo->v, v);
	if (v < pdinfo->v && comp_method == COMPACT_NONE) {
	    this_var_method = COMPACT_METHOD(pdinfo, v);
	}

	/* find the series information in the database */
	if (db_type == GRETL_RATS_DB) {
	    err = get_rats_series_info_by_name (series, &sinfo);
	} else {	
	    err = get_native_series_info (series, &sinfo);
	} 

	if (err) {
	    return 1;
	}

	/* temporary dataset */
	dbZ = new_dbZ(sinfo.nobs);
	if (dbZ == NULL) {
	    strcpy(gretl_errmsg, _("Out of memory!"));
	    return 1;
	}

	if (db_type == GRETL_RATS_DB) {
	    err = get_rats_data_by_offset (db_name, &sinfo, dbZ);
	} else {
	    get_native_db_data (db_name, &sinfo, dbZ);
	}

	if (!err) {
	    err = cli_add_db_data(dbZ, &sinfo, pZ, pdinfo, this_var_method, v);
	}

	/* free up temp stuff */
	free(dbZ[1]);
	free(dbZ);

	if (!err) {
	    pprintf(prn, _("Series imported OK"));
	    pputs(prn, "\n");
	}
    }
    
    return err;
}

void get_db_padding (SERIESINFO *sinfo, DATAINFO *pdinfo, 
			 int *pad1, int *pad2)
{
    *pad1 = dateton(sinfo->stobs, pdinfo); 
    *pad2 = pdinfo->n - sinfo->nobs - *pad1;
} 


int check_db_import (SERIESINFO *sinfo, DATAINFO *pdinfo)
{
    double sd0, sdn_new, sdn_old;

    if (sinfo->pd < pdinfo->pd) {
	strcpy(gretl_errmsg, _("You can't add a lower frequency series to a\nhigher "
	       "frequency working data set."));
	return 1;
    }

    sd0 = get_date_x(sinfo->pd, sinfo->stobs);
    sdn_new = get_date_x(sinfo->pd, sinfo->endobs);
    sdn_old = get_date_x(pdinfo->pd, pdinfo->endobs);
    if (sd0 > sdn_old || sdn_new < pdinfo->sd0) {
	strcpy(gretl_errmsg, _("Observation range does not overlap\nwith the working "
	       "data set"));
	return 1;
    }

    return 0;
}

static int cli_add_db_data (double **dbZ, SERIESINFO *sinfo, 
			    double ***pZ, DATAINFO *pdinfo,
			    int compact_method, int dbv)
{
    int err = 0;
    double *xvec;
    int n, t, start, stop, pad1 = 0, pad2 = 0;
    int new = (dbv == pdinfo->v);
    

    if (check_db_import(sinfo, pdinfo)) {
	return 1;
    }

    /* the data matrix may still be empty */
    if (pdinfo->v == 0) {
	pdinfo->v = 1;
	dbv = 1;
	if (start_new_Z(pZ, pdinfo, 0)) {
	    strcpy(gretl_errmsg, _("Out of memory adding series"));
	    return 1;
	}
    }

#ifdef DB_DEBUG
    fprintf(stderr, "Z=%p\n", (void *) *pZ);
    fprintf(stderr, "pdinfo->n = %d, pdinfo->v = %d, pdinfo->varname = %p\n",
	    pdinfo->n, pdinfo->v, (void *) pdinfo->varname);
#endif

    if (new && dataset_add_vars(1, pZ, pdinfo)) {
	strcpy(gretl_errmsg, _("Out of memory adding series"));
	return 1;
    }

    n = pdinfo->n;

    /* is the frequency of the new var higher? */
    if (sinfo->pd > pdinfo->pd) {
	if (pdinfo->pd != 1 && pdinfo->pd != 4 &&
	    sinfo->pd != 12) {
	    strcpy(gretl_errmsg, _("Sorry, can't handle this conversion yet!"));
	    if (new) dataset_drop_vars(1, pZ, pdinfo);
	    return 1;
	}
	if (compact_method == COMPACT_NONE) {
	    strcpy(gretl_errmsg, _("You must specify a compaction method"));
	    if (new) dataset_drop_vars(1, pZ, pdinfo);
	    return 1;
	}
	if (sinfo->pd == 12 && pdinfo->pd == 4) {
	    mon_to_quart(&xvec, dbZ[1], sinfo, compact_method);
	}
	else if (pdinfo->pd == 1) {
	    to_annual(&xvec, dbZ[1], sinfo, compact_method);
	}
    } else {  
	/* series does not need compacting */
	xvec = malloc(sinfo->nobs * sizeof *xvec);
	if (xvec == NULL) err = 1;
	else {
	    for (t=0; t<sinfo->nobs; t++) {
		xvec[t] = dbZ[1][t];
	    }
	} 
    }

    if (err) {
	strcpy(gretl_errmsg, _("Out of memory adding series"));
	if (new) dataset_drop_vars(1, pZ, pdinfo);
	return 1;
    }

    /* common stuff for adding a var */
    strcpy(pdinfo->varname[dbv], sinfo->varname);
    strcpy(VARLABEL(pdinfo, dbv), sinfo->descrip);
    COMPACT_METHOD(pdinfo, dbv) = compact_method;
    get_db_padding(sinfo, pdinfo, &pad1, &pad2);

    if (pad1 > 0) {
#ifdef DB_DEBUG
	fprintf(stderr, "Padding at start, %d obs\n", pad1);
#endif
	for (t=0; t<pad1; t++) {
	    (*pZ)[dbv][t] = NADBL;
	}
	start = pad1;
    } else start = 0;

    if (pad2 > 0) {
#ifdef DB_DEBUG
	fprintf(stderr, "Padding at end, %d obs\n", pad2);
#endif
	for (t=n-1; t>=n-1-pad2; t--) {
	    (*pZ)[dbv][t] = NADBL;
	}
	stop = n - pad2;
    } else stop = n;

    /* fill in actual data values */
#ifdef DB_DEBUG
    fprintf(stderr, "Filling in values from %d to %d\n", start, stop - 1);
#endif
    for (t=start; t<stop; t++) {
	(*pZ)[dbv][t] = xvec[t-pad1];
    }

    free(xvec);

    return 0;
}

