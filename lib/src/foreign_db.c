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

/**
 * SECTION:foreign_db
 * @short_description: reading from RATS and PcGive databases
 * @title: Foreign DB
 * @include: gretl/libgretl.h, gretl/dbread.h
 *
 * Functions that read data from RATS 4.0 and PcGive databases.
 */

#include "libgretl.h"
#include "dbread.h"
#include "swap_bytes.h"
#include "foreign_db.h"

#define RATSCOMMENTLENGTH 80
#define RATSCOMMENTS 2
#define RATS_PARSE_ERROR -999

#define RECNUM gint32
#define NAMELENGTH 16

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

int get_rats_series_info (const char *series_name, SERIESINFO *sinfo)
{
    const char *dbn = get_db_name();
    FILE *fp;
    long forward = 0;
    int err = 0;

    gretl_error_clear();

    fp = gretl_fopen(dbn, "rb");
    if (fp == NULL) {
	return E_FOPEN;
    }

#if DB_DEBUG
    fprintf(stderr, "Opened %s\n", dbn);
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

/* end of RATS database material, below: PcGive databases */

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

int get_pcgive_series_info (const char *series, SERIESINFO *sinfo)
{
    FILE *fp;
    char dbidx[MAXLEN];
    char line[1024];
    char fmt[24];
    char *p;
    int y0, p0, y1, p1;
    int nf, gotit = 0;
    int err = 0;

    strcpy(dbidx, get_db_name());
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

static FILE *open_bn7 (const char *dbbase, int offset, int *err)
{
    char dbbin[MAXLEN];
    FILE *fp = NULL;

    strcpy(dbbin, dbbase);
    if (strstr(dbbin, ".bn7") == NULL) {
        strcat(dbbin, ".bn7");
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

    fp = open_bn7(dbbase, sinfo->offset, &err);
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


