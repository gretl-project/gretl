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
#include "dbread.h"

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

/* ........................................................... */

static int get_rats_series (int offset, SERIESINFO *sinfo, FILE *fp, 
			    double ***pZ)
/* print the actual data values from the data blocks */
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
	    sprintf(numstr, "%.3f ", rdata.data[i]);
	    val = atof(numstr);
	    if (isnan(val)) {
		val = NADBL;
		miss = 1;
	    }
	    (*pZ)[1][t] = val;
	    t++;
	}
    }

    return miss;
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
    if (pd == 1)
	sprintf(datestr, "%d", endyr);
    else if (pd == 4)
	sprintf(datestr, "%d.%d", endyr, endfrac);
    else if (pd == 12)
	sprintf(datestr, "%d.%02d", endyr, endfrac);
    return 0;
}

/* ........................................................... */

static RECNUM read_RATS_directory (FILE *fp, db_table_row *row) 
/* read the RATS directory struct.  Note that we can't do this
   in one gulp, since the info is packed to 256 bytes in the RATS
   file, which is more compact than the C struct we're reading
   the info into, due to padding in the latter. */
{
    RATSDirect rdir;
    DATEINFO dinfo;
    char pd = 0, pdstr[3], endobs[9], datestuff[48];    
    int startfrac = 0;

    fread(&rdir.back_point, sizeof(RECNUM), 1, fp);
    fread(&rdir.forward_point, sizeof(RECNUM), 1, fp);
    fseek(fp, 4L, SEEK_CUR); /* skip two shorts */
    fread(&rdir.first_data, sizeof(RECNUM), 1, fp);
    fread(rdir.series_name, 16, 1, fp);  
    rdir.series_name[8] = '\0';
    chopstr(rdir.series_name);

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

    if ((int) dinfo.info == 4) {
	pd = 'Q';
	sprintf(pdstr, ".%d", dinfo.month);
	if (dinfo.month == 1) startfrac = 1;
	else if (dinfo.month > 1 && dinfo.month <= 4) startfrac = 2;
	else if (dinfo.month > 4 && dinfo.month <= 7) startfrac = 3;
	else startfrac = 4;
    }
    else if ((int) dinfo.info == 12) {
	pd = 'M';
	sprintf(pdstr, ".%02d", dinfo.month);
	startfrac = dinfo.month;
    }
    else if ((int) dinfo.info == 1) {
	pd = 'A';
	strcpy(pdstr, "");
	startfrac = 0;
    }

    get_endobs(endobs, dinfo.year, startfrac, dinfo.info, 
	       rdir.datapoints);

    /* add info to trawl */
    row->varname = g_strdup(rdir.series_name);
    row->comment = g_strdup(rdir.comments[0]);
    sprintf(datestuff, "%c  %d%s - %s  n = %d", pd, (int) dinfo.year, 
	   pdstr, endobs, (int) rdir.datapoints);
    row->obsinfo = g_strdup(datestuff);

    return rdir.forward_point;
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

    rows = realloc(tbl->rows, 3 * sizeof *rows);
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

/* ........................................................... */

db_table *read_RATS_db (FILE *fp) 
/* read the base block at offset 0 in the data file,
   and recurse through the directory entries */
{
    long forward;
    db_table *tbl;
    int i, err = 0;
    
    tbl = db_table_new();
    if (tbl == NULL) return NULL;

    fseek(fp, 30L, SEEK_SET); /* skip unneeded fields */
    fread(&forward, sizeof forward, 1, fp);
    fseek(fp, 4L, SEEK_CUR);

    /* Go find the series */
    i = 0;
    while (forward && !err) {
	tbl->nrows += 1;
	if (tbl->nrows > 0 && tbl->nrows % DB_INIT_ROWS == 0) {
	    err = db_table_expand(tbl);
	}
	if (!err) {
	    fseek(fp, (forward - 1) * 256L, SEEK_SET);
	    forward = read_RATS_directory(fp, &tbl->rows[i++]);
	}
    }

    if (err) {
	db_table_free(tbl);
	return NULL;
    }

    return tbl;
}

/* ........................................................... */

static int find_RATSDirect (FILE *fp, const int first_dir, 
			    const int series_number)
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

static int get_rats_series_offset (FILE *fp, const int series_number)
{
    long num_series, first_dir;
    int offset;

    fseek(fp, 6L, SEEK_SET);
    fread(&num_series, sizeof num_series, 1, fp);
    if (series_number > num_series) return -1;
    fseek(fp, sizeof(long) * 5L, SEEK_CUR);  
    fread(&first_dir, sizeof first_dir, 1, fp);
    offset = find_RATSDirect(fp, first_dir, series_number); 

    return offset;
}

/* ........................................................... */

int get_rats_data (const char *fname, const int series_number,
		   SERIESINFO *sinfo, double ***pZ)
/* series are numbered from 1 for this function.
   We need to know the specific filename. */
{
    FILE *fp;
    int offset;
    long first_data;
    int ret = DB_OK;

    fp = fopen(fname, "rb");
    if (fp == NULL) return DB_NOT_FOUND;
    
    offset = get_rats_series_offset(fp, series_number);
    if (offset < 0) return DB_NOT_FOUND;
    
    fseek(fp, (offset - 1) * 256 + 12, SEEK_SET); 
    fread(&first_data, sizeof(RECNUM), 1, fp);
    if (get_rats_series(first_data, sinfo, fp, pZ)) {
	ret = DB_MISSING_DATA;
    }

    fclose(fp);

    return ret;
}

/* ........................................................... */

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

int mon_to_quart (double **pq, double *mvec, SERIESINFO *sinfo,
		  gint method)
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

    fprintf(stderr, "startskip = %d\n", skip);
    endskip = (sinfo->nobs - skip) % 3;
    fprintf(stderr, "endskip = %d\n", endskip);
    goodobs = (sinfo->nobs - skip - endskip) / 3;
    fprintf(stderr, "goodobs = %d\n", goodobs);
    sinfo->nobs = goodobs;
    sprintf(sinfo->stobs, "%d.%d", y0, q0);
    fprintf(stderr, "starting date = %s\n", sinfo->stobs);

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
	       gint method)
{
    int i, t, p, pmax = 0, p0, y0, skip = 0, endskip, goodobs;
    int pd = sinfo->pd;
    double val;
    char numstr[16];

    /* record the precision of the original data */
    for (t=0; t<sinfo->nobs; t++) {
	p = get_places(mvec[t]);
	if (p > pmax) pmax = p;
    }

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
	else if (method == COMPACT_EOP) 
	    val = mvec[p-1+skip];
	else if (method == COMPACT_SOP)
	    val = mvec[p-pd+skip];
	sprintf(numstr, "%.*f", pmax, val);
	(*pq)[t] = atof(numstr);
    }

    sinfo->pd = 1;

    return 0;
}
