/*
 *  Copyright (c) by Ramu Ramanathan and Allin Cottrell
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
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

#include "libgretl.h"
#include "internal.h"

#include <zlib.h>
#include <ctype.h>
#include <time.h>
#include <errno.h>
#include <sys/stat.h>
#include <unistd.h>

#include <libxml/xmlmemory.h>
#include <libxml/parser.h>

#define QUOTE                  '\''
#define SCALAR_DIGITS           12

#define IS_DATE_SEP(c) (c == '.' || c == ':' || c == ',')

static int prepZ (double ***pZ, const DATAINFO *pdinfo);
static int writelbl (const char *lblfile, const int *list, 
		     const DATAINFO *pdinfo);
static int writehdr (const char *hdrfile, const int *list, 
		     const DATAINFO *pdinfo, int opt);
static double obs_float (const DATAINFO *pdinfo, int end);
static int write_xmldata (const char *fname, const int *list, 
			  double **Z, const DATAINFO *pdinfo, 
			  int opt, PATHS *ppaths);
static int xmlfile (const char *fname);

static char STARTCOMMENT[3] = "(*";
static char ENDCOMMENT[3] = "*)";

/* from strutils.c */
extern void check_for_console (PRN *prn);
extern void console_off (void);

#define PROGRESS_BAR "progress_bar"

enum {
    NO_MARKERS = 0,
    REGULAR_MARKERS,
    DAILY_DATE_STRINGS
};

/**
 * free_Z:
 * @Z: data matrix.
 * @pdinfo: data information struct.
 *
 * Do a deep free on the data matrix.
 * 
 */

void free_Z (double **Z, DATAINFO *pdinfo)
{
    int i;

    if (Z == NULL) return;

    for (i=0; i<pdinfo->v; i++) free(Z[i]);
    free(Z);
}

/**
 * clear_datainfo:
 * @pdinfo: data information struct.
 * @code: either CLEAR_FULL or CLEAR_SUBSAMPLE.
 *
 * Free the allocated content of a data information struct.
 * 
 */

void clear_datainfo (DATAINFO *pdinfo, int code)
{
    int i;

    if (pdinfo->S != NULL) {
	for (i=0; i<pdinfo->n; i++) 
	   free(pdinfo->S[i]); 
	free(pdinfo->S);
	pdinfo->S = NULL;
	pdinfo->markers = NO_MARKERS;
    } 

    /* if this is not a sub-sample datainfo, free varnames, labels, etc. */
    if (code == CLEAR_FULL) {
	if (pdinfo->varname != NULL) {
	    for (i=0; i<pdinfo->v; i++) 
		free(pdinfo->varname[i]); 
	    free(pdinfo->varname);
	    pdinfo->varname = NULL;
	}
	if (pdinfo->varinfo != NULL) {
	    for (i=0; i<pdinfo->v; i++) 
		free(pdinfo->varinfo[i]); 
	    free(pdinfo->varinfo);
	    pdinfo->varinfo = NULL;
	}
	if (pdinfo->descrip) {
	    free(pdinfo->descrip);
	    pdinfo->descrip = NULL;
	}
	if (pdinfo->vector) {
	    free(pdinfo->vector);
	    pdinfo->vector = NULL;
	}
    }
}

/* ......................................................... */

static void dataset_dates_defaults (DATAINFO *pdinfo)
{
    strcpy(pdinfo->stobs, "1");
    sprintf(pdinfo->endobs, "%d", pdinfo->n);
    pdinfo->sd0 = 1.0;
    pdinfo->pd = 1;
    pdinfo->time_series = 0;
    pdinfo->extra = 0; 
    pdinfo->decpoint = '.';
}

/* ......................................................... */

double get_date_x (int pd, const char *obs)
{
    double x = 1.0;

    if ((pd == 5 || pd == 7) && strlen(obs) > 4) { 
	/* daily data */
	long ed = get_epoch_day(obs);

	if (ed >= 0) x = ed;
    } else 
	x = obs_str_to_double(obs); 

    return x;
}

/* ......................................................... */

static int skipcomments (FILE *fp, const char *str)
/* Skips past comments in .hdr file.  Returns 
   0 if comments found, otherwise 1.
*/
{
    char commentword[MAXLEN];  /* should be big enough to accommodate
			          strings among the comments? */
    commentword[0] = '\0';
    if (strncmp(str, STARTCOMMENT, 2) == 0) {
        while (strcmp(commentword, ENDCOMMENT)) {
            fscanf(fp, "%s", commentword);
        }
        return 0;
    } else return 1;
}

/* ................................................. */

static int comment_lines (FILE *fp, char **pbuf)
{
    char s[MAXLEN], *mybuf = NULL;
    int count = 0, bigger = 1, bufsize;

    if (fgets(s, MAXLEN-1, fp) == NULL) return 0;

    if (!strncmp(s, STARTCOMMENT, 2)) {
	*pbuf = malloc(20 * MAXLEN);
	if (*pbuf == NULL) return -1;
	**pbuf = '\0';
	do {
	    if (fgets(s, MAXLEN-1, fp) == NULL) break;
	    if (!strncmp(s, ENDCOMMENT, 2)) break;
	    count++;
	    if (count > 20*bigger) {
		bigger++;
		bufsize = 20 * MAXLEN * bigger;
		mybuf = realloc(*pbuf, bufsize);
		if (mybuf == NULL) return -1;
		else *pbuf = mybuf;
	    }
	    strcat(*pbuf, s);
	} while (s != NULL);
    }

    return count;
}

/* ................................................. */

static int dataset_allocate_markers (DATAINFO *pdinfo)
{
    int i, k;

    pdinfo->S = malloc(pdinfo->n * sizeof *pdinfo->S);

    if (pdinfo->S == NULL) return 1; 

    for (i=0; i<pdinfo->n; i++) {
	pdinfo->S[i] = malloc(OBSLEN);
	if (pdinfo->S[i] == NULL) {
	    for (k=0; k<i; k++) free(pdinfo->S[k]);
	    free(pdinfo->S);
	    pdinfo->S = NULL;
	    return 1; 
	}
    }

    return 0;
}

/* ................................................. */

void gretl_varinfo_init (VARINFO *vinfo)
{
    *vinfo->label = '\0';
    *vinfo->display_name = '\0';
    vinfo->compact_method = COMPACT_NONE;
}

/* ................................................. */

static char *compact_method_to_string (int method)
{
    if (method == COMPACT_SUM) return "COMPACT_SUM";
    else if (method == COMPACT_AVG) return "COMPACT_AVG";
    else if (method == COMPACT_SOP) return "COMPACT_SOP";
    else if (method == COMPACT_EOP) return "COMPACT_EOP";
    else return "COMPACT_NONE";
}

static int compact_string_to_int (const char *str)
{
    if (!strcmp(str, "COMPACT_SUM")) return COMPACT_SUM;
    else if (!strcmp(str, "COMPACT_AVG")) return COMPACT_AVG;
    else if (!strcmp(str, "COMPACT_SOP")) return COMPACT_SOP;
    else if (!strcmp(str, "COMPACT_EOP")) return COMPACT_EOP;
    else return COMPACT_NONE;
}

/* ................................................. */

static int dataset_allocate_varnames (DATAINFO *pdinfo)
{
    int i, v = pdinfo->v;
    
    pdinfo->varname = malloc(v * sizeof *pdinfo->varname);
    pdinfo->varinfo = malloc(v * sizeof *pdinfo->varinfo);
    pdinfo->vector = malloc(v * sizeof *pdinfo->vector);

    if (pdinfo->varname == NULL || 
	pdinfo->varinfo == NULL ||
	pdinfo->vector == NULL) return 1;

    for (i=0; i<v; i++) {
	pdinfo->varname[i] = malloc(VNAMELEN);
	if (pdinfo->varname[i] == NULL) return 1;

	pdinfo->varname[i][0] = '\0';
	pdinfo->varinfo[i] = malloc(sizeof **pdinfo->varinfo);
	if (pdinfo->varinfo[i] == NULL) return 1;

	gretl_varinfo_init(pdinfo->varinfo[i]);
	pdinfo->vector[i] = 1;
    }

    strcpy(pdinfo->varname[0], "const");
    strcpy(VARLABEL(pdinfo, 0), _("auto-generated constant"));
    return 0;
}

#ifdef notdef
static char locale_friendly_delim (void)
{
    if (',' == get_local_decpoint())
	return ' ';
    else
	return ',';
}
#endif

/**
 * datainfo_new:
 *
 * Create a new data information struct pointer from scratch,
 * properly initialized.
 * 
 * Returns: pointer to data information struct, or NULL on error.
 *
 */

DATAINFO *datainfo_new (void)
{
    DATAINFO *dinfo;

    dinfo = malloc(sizeof *dinfo);
    if (dinfo == NULL) return NULL;

    dinfo->v = 0;
    dinfo->n = 0;
    dinfo->pd = 1;
    dinfo->bin = 0;
    dinfo->extra = 0;
    dinfo->sd0 = 1.0;
    dinfo->t1 = 0;
    dinfo->t2 = 0;
    *dinfo->stobs = '\0';
    *dinfo->endobs = '\0';
    dinfo->varname = NULL;
    dinfo->varinfo = NULL;    
    dinfo->markers = NO_MARKERS;  
    dinfo->delim = ',';
    dinfo->decpoint = '.';
    dinfo->S = NULL;
    dinfo->descrip = NULL;
    dinfo->vector = NULL;
    dinfo->data = NULL;
    dinfo->time_series = 0;

    return dinfo;
}

/**
 * create_new_dataset:
 * @pZ: pointer to data matrix.
 * @nvar: number of variables.
 * @nobs: number of observations per variable 
 * @markers: 1 if there are case markers for the observations, 0
 * otherwise.
 *
 * Create a new data information struct corresponding to a given
 * data matrix.
 * 
 * Returns: pointer to data information struct, or NULL on error.
 *
 */

DATAINFO *create_new_dataset (double ***pZ,  
			      int nvar,
			      int nobs,
			      int markers
			      )
{
    DATAINFO *pdinfo;

    pdinfo = malloc(sizeof *pdinfo);
    if (pdinfo == NULL) return NULL;

    pdinfo->v = nvar;
    pdinfo->n = nobs;
    *pZ = NULL;

    if (start_new_Z(pZ, pdinfo, 0)) {
	free(pdinfo);
	return NULL;
    }

    pdinfo->markers = (unsigned char) markers;
    if (pdinfo->markers) {
	if (dataset_allocate_markers(pdinfo)) {
	    free_datainfo(pdinfo);
	    return NULL;
	}
    } 

    dataset_dates_defaults(pdinfo);
    pdinfo->descrip = NULL;

    return pdinfo;
}

/**
 * start_new_Z:
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 * @resample: 1 if we're sub-sampling from a full data set, 0 otherwise.
 *
 * Initialize data matrix (add the constant) and the data information
 * struct.
 * 
 * Returns: 0 on successful completion, 1 on error.
 *
 */

int start_new_Z (double ***pZ, DATAINFO *pdinfo, int resample)
{
    if (prepZ(pZ, pdinfo)) return 1;

    pdinfo->t1 = 0; 
    pdinfo->t2 = pdinfo->n - 1;

    if (resample) {
	pdinfo->varname = NULL;
	pdinfo->varinfo = NULL;
    } 
    else if (dataset_allocate_varnames(pdinfo)) {
	return 1;
    }

    pdinfo->S = NULL;
    pdinfo->markers = NO_MARKERS;
    pdinfo->delim = ',';
    pdinfo->descrip = NULL;
    pdinfo->data = NULL;
    
    return 0;
}

/* ................................................. */

static int prepZ (double ***pZ, const DATAINFO *pdinfo)
    /* allocate data array and put in constant */
{
    int i, t;

    if (*pZ != NULL) free(*pZ);
    *pZ = malloc(pdinfo->v * sizeof **pZ);
    if (*pZ == NULL) return 1;

    for (i=0; i<pdinfo->v; i++) {
	(*pZ)[i] = malloc(pdinfo->n * sizeof ***pZ);
	if (*pZ == NULL) return 1;
    }

    for (t=0; t<pdinfo->n; t++) (*pZ)[0][t] = 1.0; 
    return 0;
}

/* ................................................. */

static void eatspace (FILE *fp)
{
    char c;

    while (1) {
	c = fgetc(fp);
	if (!isspace((unsigned char) c)) {
	    ungetc(c, fp);
	    return;
	}
    }
}

/* ................................................. */

static int readdata (FILE *fp, const DATAINFO *pdinfo, double **Z)
{
    int i, t, n = pdinfo->n;
    char c, marker[OBSLEN];
    int err = 0;
    float x;

    gretl_errmsg[0] = '\0';

    if (pdinfo->bin == 1) { /* single-precision binary data */
	for (i=1; i<pdinfo->v; i++) {
	    for (t=0; t<n; t++) {
		if (!fread(&x, sizeof x, 1, fp)) {
		    sprintf(gretl_errmsg, _("WARNING: binary data read error at "
			    "var %d"), i);
		    return 1;
		}
		Z[i][t] = (double) x;
	    }
	}
    }
    else if (pdinfo->bin == 2) { /* double-precision binary data */
	for (i=1; i<pdinfo->v; i++) {
	    if (!fread(Z[i], sizeof(double), n, fp)) {
		sprintf(gretl_errmsg, 
			_("WARNING: binary data read error at var %d"), i);
		return 1;
	    }
	}
    } else { /* ascii data file */
#ifdef ENABLE_NLS
	setlocale(LC_NUMERIC, "C");
#endif
	for (t=0; t<n && !err; t++) {
	    eatspace(fp);
	    c = fgetc(fp);  /* test for a #-opened comment line */
	    if (c == '#') {
		while (c != '\n') c = fgetc(fp);
	    } else {
		ungetc(c, fp);
	    }
	    if (pdinfo->markers) {
		fscanf(fp, "%8s", marker);
		strcpy(pdinfo->S[t], marker);
	    }
	    for (i=1; i<pdinfo->v; i++) {
		if ((fscanf(fp, "%lf", &Z[i][t])) != 1) {
		    sprintf(gretl_errmsg, 
			    _("WARNING: ascii data read error at var %d, "
			    "obs %d"), i, t + 1);
		    err = 1;
		    break;
		}
	    }
	}
#ifdef ENABLE_NLS
	setlocale(LC_NUMERIC, "");
#endif
    }

    return err;
}

/* ................................................. */

static int gz_readdata (gzFile fz, const DATAINFO *pdinfo, double **Z)
{
    int i, t, n = pdinfo->n;
    int err = 0;
    
    gretl_errmsg[0] = '\0';

    if (pdinfo->bin == 1) { /* single-precision binary data */
	float xx;

	for (i=1; i<pdinfo->v; i++) {
	    for (t=0; t<n; t++) {
		if (!gzread(fz, &xx, sizeof xx)) {
		    sprintf(gretl_errmsg, _("WARNING: binary data read error at "
			    "var %d"), i);
		    return 1;
		}
		Z[i][t] = (double) xx;
	    }
	}
    }
    else if (pdinfo->bin == 2) { /* double-precision binary data */
	for (i=1; i<pdinfo->v; i++) {
	    if (!gzread(fz, &Z[i][0], n * sizeof(double))) {
		sprintf(gretl_errmsg, 
			_("WARNING: binary data read error at var %d"), i);
		return 1;
	    }
	}
    } else { /* ascii data file */
	char *line, numstr[24];
	int llen = pdinfo->v * 32;
	size_t offset;

	line = malloc(llen);
	if (line == NULL) return E_ALLOC;

#ifdef ENABLE_NLS
	setlocale(LC_NUMERIC, "C");
#endif

	for (t=0; t<n; t++) {
	    offset = 0L;
	    if (!gzgets(fz, line, llen - 1)) {
		sprintf(gretl_errmsg, _("WARNING: ascii data read error at "
			"obs %d"), t + 1);
		err = 1;
		break;
	    }
	    chopstr(line);
	    compress_spaces(line);
	    if (line[0] == '#') {
		t--;
		continue;
	    }
	    if (pdinfo->markers) {
		if (sscanf(line, "%8s", pdinfo->S[t]) != 1) {
		   sprintf(gretl_errmsg, 
			   _("WARNING: failed to read case marker for "
			   "obs %d"), t + 1);
		   err = 1;
		   break;
		}
		pdinfo->S[t][8] = 0;
		offset += strlen(pdinfo->S[t]) + 1;
	    }
	    for (i=1; i<pdinfo->v; i++) {
		if (sscanf(line + offset, "%23s", numstr) != 1) {
		    sprintf(gretl_errmsg, 
			    _("WARNING: ascii data read error at var %d, "
			    "obs %d"), i, t + 1);
		    err = 1;
		    break;
		}
		numstr[23] = 0;
		Z[i][t] = atof(numstr);
		if (i < pdinfo->v - 1) offset += strlen(numstr) + 1;
	    }
	    if (err) break;
	}
	free(line);

#ifdef ENABLE_NLS
	setlocale(LC_NUMERIC, "");
#endif

    }
    return err;
}

/**
 * check_varname:
 * @varname: putative name for variable.
 * 
 * Check a variable name for legality.
 * 
 * Returns: 0 if name is OK, 1 if not.
 *
 */

int check_varname (const char *varname)
{
    int i, n = strlen(varname);

    *gretl_errmsg = '\0';

    if (_reserved(varname)) return 1;
    
    if (!(isalpha((unsigned char) *varname))) {
        sprintf(gretl_errmsg, _("First char of varname ('%c') is bad\n"
               "(first must be alphabetical)"), *varname);
        return 1;
    }

    for (i=1; i<n; i++) {
        if (!(isalpha((unsigned char) varname[i]))  
            && !(isdigit((unsigned char) varname[i]))
            && varname[i] != '_') {
	    if (isprint((unsigned char) varname[i]))
		sprintf(gretl_errmsg, _("Varname contains illegal character '%c'\n"
			"Use only letters, digits and underscore"), 
			varname[i]);
	    else
		sprintf(gretl_errmsg, _("Varname contains illegal character 0x%x\n"
			"Use only letters, digits and underscore"), 
			(unsigned) varname[i]);
            return 1;
        }
    }

    return 0;
}   

/* ................................................ */

static int readhdr (const char *hdrfile, DATAINFO *pdinfo)
{
    FILE *fp;
    int n, i = 0, panel = 0, descrip = 0;
    char str[MAXLEN], byobs[6], option[8];

    *gretl_errmsg = '\0';

    fp = fopen(hdrfile, "r");
    if (fp == NULL) {
	sprintf(gretl_errmsg, _("Couldn't open file %s"),  hdrfile);
	return E_FOPEN;
    }
    fscanf(fp, "%s", str);
    i += skipcomments(fp, str); 
    while (1) { /* find number of variables */
        if (fscanf(fp, "%s", str) != 1) {
	    fclose(fp);
	    sprintf(gretl_errmsg, _("Opened header file %s\n"
		    "Couldn't find list of variables (must "
		    "be terminated with a semicolon)"), hdrfile);
	    return 1;
	}
	n = strlen(str);
	if (str[n-1] == ';') {
	    if (n > 1) i++;
	    break;
	} else i++;
    }
    pdinfo->v = i + 1;
    fclose(fp);

    pdinfo->S = NULL;
    if (dataset_allocate_varnames(pdinfo)) return E_ALLOC;

    i = 1;
    fp = fopen(hdrfile, "r");
    str[0] = 0;
    fscanf(fp, "%s", str);
    if (skipcomments(fp, str)) {
        safecpy(pdinfo->varname[i], str, 8);
	if (check_varname(pdinfo->varname[i++])) 
	    goto varname_error;
    } else {
	descrip = 1; /* comments were found */
    }

    while (1) {
        fscanf(fp, "%s", str);
	n = strlen(str);
	if (str[n-1] != ';') {
            safecpy(pdinfo->varname[i], str, 8);
	    if (check_varname(pdinfo->varname[i++])) 
		goto varname_error;
        } else {
	    if (n > 1) {
		safecpy(pdinfo->varname[i], str, n-1);
		pdinfo->varname[i][n] = '\0';
		if (check_varname(pdinfo->varname[i]))
		    goto varname_error; 
	    }
	    break;
	}
    }

    fscanf(fp, "%d", &pdinfo->pd);
    fscanf(fp, "%s", pdinfo->stobs);
    fscanf(fp, "%s", pdinfo->endobs);

    colonize_obs(pdinfo->stobs);
    colonize_obs(pdinfo->endobs);

    pdinfo->sd0 = get_date_x(pdinfo->pd, pdinfo->stobs);

    if (pdinfo->sd0 >= 2.0) 
        pdinfo->time_series = TIME_SERIES; /* actual time series? */
    else if (pdinfo->sd0 > 1.0) {
	pdinfo->time_series = STACKED_TIME_SERIES; /* panel data? */
    }
    else pdinfo->time_series = 0;

    pdinfo->n = -1;
    pdinfo->n = dateton(pdinfo->endobs, pdinfo) + 1;
    pdinfo->extra = 0;      

    pdinfo->bin = 0;
    pdinfo->markers = NO_MARKERS;

    if (fscanf(fp, "%5s %7s", byobs, option) == 2) {
	if (strcmp(option, "SINGLE") == 0)
	    pdinfo->bin = 1;
	else if (strcmp(option, "BINARY") == 0)
	    pdinfo->bin = 2;
	else if (strcmp(option, "MARKERS") == 0) 
	    pdinfo->markers = 1;
	else if (strcmp(option, "PANEL2") == 0) {
	    panel = 1;
	    pdinfo->time_series = STACKED_TIME_SERIES;
	} else if (strcmp(option, "PANEL3") == 0) {
	    panel = 1;
	    pdinfo->time_series = STACKED_CROSS_SECTION;
	}
    }

    if (!panel && fscanf(fp, "%6s", option) == 1) {
	if (strcmp(option, "PANEL2") == 0)
	    pdinfo->time_series = STACKED_TIME_SERIES;
	else if (strcmp(option, "PANEL3") == 0)
	    pdinfo->time_series = STACKED_CROSS_SECTION;
    }

    if (fp != NULL) fclose(fp);

    /* last pass, to pick up comments */
    pdinfo->descrip = NULL;
    if (descrip) {
	char *dbuf = NULL;
	int lines;

	fp = fopen(hdrfile, "r");
	if (fp == NULL) return 0;
	if ((lines = comment_lines(fp, &dbuf)) > 0) {
	    delchar('\r', dbuf);
	    pdinfo->descrip = malloc(strlen(dbuf) + 1);
	    if (pdinfo->descrip != NULL) {
		strcpy(pdinfo->descrip, dbuf);
	    }
	    free(dbuf);
	}
	else if (lines < 0) {
	    fprintf(stderr, I_("Failed to store data comments\n"));
	}
	fclose(fp);
    } 
	
    return 0;

    varname_error:
    fclose(fp);
    clear_datainfo(pdinfo, CLEAR_FULL);
    return E_DATA;
}

/* ................................................ */

static int check_date (const char *date)
{
    int i, n = strlen(date);

    *gretl_errmsg = 0;

    for (i=0; i<n; i++) {
	if (!isdigit((unsigned char) date[i]) && !IS_DATE_SEP(date[i])) {
	    if (isprint((unsigned char) date[i]))
		sprintf(gretl_errmsg, 
			_("Bad character '%c' in date string"), date[i]);
	    else 
		sprintf(gretl_errmsg, 
			_("Bad character %d in date string"), date[i]);
	    return 1;
	}
    }
    return 0;
}

/**
 * dateton:
 * @date: string representation of date for processing.
 * @pdinfo: pointer to data information struct.
 * 
 * Given a "current" date string, a periodicity, and a starting
 * date string, returns the observation number corresponding to
 * the current date string, counting from zero.
 * 
 * Returns: integer observation number.
 *
 */

int dateton (const char *date, const DATAINFO *pdinfo)
{
    int dotpos1 = 0, dotpos2 = 0, maj = 0, min = 0, n, i;
    char majstr[5], minstr[3];
    char startmajstr[5], startminstr[3];
    int startmaj, startmin;

    if (dated_daily_data(pdinfo)) {
	if (pdinfo->markers && pdinfo->S != NULL) {
	    for (i=0; i<pdinfo->n; i++) {
		if (!strcmp(date, pdinfo->S[i])) {
		    return i;
		}
	    }
	    return -1;
	} else {
	    return daily_obs_number(date, pdinfo);
	}
    }

    if (check_date(date)) {
	return -1;
    }

    n = strlen(date);
    for (i=1; i<n; i++) {
        if (IS_DATE_SEP(date[i])) {
	    dotpos1 = i;
	    break;
	}
    }

    if (dotpos1) {
        safecpy(majstr, date, dotpos1);
        maj = atoi(majstr);
        strcpy(minstr, date + dotpos1 + 1);
        min = atoi(minstr);
    }

    n = strlen(pdinfo->stobs);
    for (i=1; i<n; i++) {
        if (IS_DATE_SEP(pdinfo->stobs[i])) {
	    dotpos2 = i;
	    break;
	}
    }

    if ((dotpos1 && !dotpos2) || (dotpos2 && !dotpos1)) {
	sprintf(gretl_errmsg, _("Date strings inconsistent"));
	return -1;  
    }

    if (!dotpos1 && !dotpos2) {
	n = atoi(date) - atoi(pdinfo->stobs);
	if (n < 0 || (pdinfo->n != -1 && n > pdinfo->n)) {
	    sprintf(gretl_errmsg, _("Observation number out of bounds"));
	    return -1; 
	}
        else return n;
    }

    safecpy(startmajstr, pdinfo->stobs, dotpos2);
    startmaj = atoi(startmajstr);
    strcpy(startminstr, pdinfo->stobs + dotpos2 + 1);
    startmin = atoi(startminstr);
    n = pdinfo->pd * (maj - startmaj);
    n += min - startmin;
   
    return n;
}

/**
 * ntodate:
 * @datestr: string to which date is to be printed.
 * @nt: an observation number (zero-based).
 * @pdinfo: data information struct.
 * 
 * print to @datestr the calendar representation of observation
 * number nt.
 * 
 * Returns: the observation string.
 */

char *ntodate (char *datestr, int t, const DATAINFO *pdinfo)
{
    double x;
    static int decpoint;

    decpoint = get_local_decpoint();

    if (dated_daily_data(pdinfo)) {
	if (pdinfo->markers && pdinfo->S != NULL) {
	    strcpy(datestr, pdinfo->S[t]);
	} else {
	    daily_date_string(datestr, t, pdinfo);
	}
	return datestr;
    }

    x = date(t, pdinfo->pd, pdinfo->sd0);

    if (pdinfo->pd == 1) {
        int n = (int) x;

        sprintf(datestr, "%d", n);
    } else {
	if (pdinfo->pd < 10) sprintf(datestr, "%.1f", x);
	else sprintf(datestr, "%.2f", x);
	charsub(datestr, decpoint, ':');
    }
    
    return datestr;
}

/* .......................................................... */

static int blank_check (FILE *fp)
{
    int i, deflt = 1;
    char s[MAXLEN];

    for (i=0; i<3 && deflt && fgets(s, MAXLEN-1, fp) ; i++) {
	if (i == 0 && strncmp(s, "(*", 2)) deflt = 0;
	else if (i == 1 && strncmp(s, _("space for comments"), 18)) deflt = 0;
	else if (i == 2 && strncmp(s, "*)", 2)) deflt = 0;
    }

    fclose(fp);

    return deflt;
}

/**
 * get_info:
 * @hdrfile: name of data header file
 * @prn: gretl printing struct.
 * 
 * print to @prn the informative comments contained in the given
 * data file (if any).
 * 
 * Returns: 0 on successful completion, non-zero on error or if there
 * are no informative comments.
 * 
 */

int get_info (const char *hdrfile, PRN *prn)
{      
    char s[MAXLEN];
    int i = 0;
    FILE *hdr;

    if ((hdr = fopen(hdrfile, "r")) == NULL) {
	pprintf(prn, _("Couldn't open %s\n"), hdrfile); 
	return 1;
    }

    /* see if it's just the default "space for comments" */
    if (blank_check(hdr)) { /* yes */
	pprintf(prn, _("No info in %s\n"), hdrfile);
	return 2;
    } 

    /* no, so restart the read */
    if ((hdr = fopen(hdrfile, "r")) == NULL) {
	pprintf(prn, _("Couldn't open %s\n"), hdrfile); 
	return 1;
    }    

    pprintf(prn, _("Data info in file %s:\n\n"), hdrfile);

    if (fgets(s, MAXLEN-1, hdr) != NULL && strncmp(s, STARTCOMMENT, 2) == 0) {
	do {
	    if (fgets(s, MAXLEN-1, hdr) != NULL && strncmp(s, "*)", 2)) {
#ifndef WIN32
		delchar('\r', s);
#endif
		pputs(prn, s);
		i++;
	    }
	} while (s != NULL && strncmp(s, ENDCOMMENT, 2));
    }

    if (i == 0) pputs(prn, _(" (none)\n"));
    pputc(prn, '\n');

    if (hdr != NULL) fclose(hdr);
    return 0;
}

/* ................................................ */

static int writehdr (const char *hdrfile, const int *list, 
		     const DATAINFO *pdinfo, int opt)
{
    FILE *fp;
    int bin = 0, i;
    char startdate[OBSLEN], enddate[OBSLEN];

    if (opt == GRETL_DATA_FLOAT) bin = 1;
    else if (opt == GRETL_DATA_DOUBLE) bin = 2;

    ntodate(startdate, pdinfo->t1, pdinfo);
    ntodate(enddate, pdinfo->t2, pdinfo);

    fp = fopen(hdrfile, "w");
    if (fp == NULL) return 1;

    /* write description of data set, if any */
    if (pdinfo->descrip != NULL) {
	size_t len = strlen(pdinfo->descrip);

	if (len > 2) {
	    fprintf(fp, "(*\n%s%s*)\n", pdinfo->descrip,
		    (pdinfo->descrip[len-1] == '\n')? "" : "\n");
	}
    }

    /* then list of variables */
    for (i=1; i<=list[0]; i++) {
	if (list[i] == 0) continue;
	fprintf(fp, "%s ", pdinfo->varname[list[i]]);
	if (i && i <list[0] && (i+1) % 8 == 0) fprintf(fp, "\n");
    }    
    fputs(";\n", fp);

    /* then obs line */
    fprintf(fp, "%d %s %s\n", pdinfo->pd, startdate, enddate);
    
    /* and flags as required */
    if (bin == 1) fputs("BYVAR\nSINGLE\n", fp);
    else if (bin == 2) fputs("BYVAR\nBINARY\n", fp);
    else { 
	fputs("BYOBS\n", fp);
	if (pdinfo->markers) fputs("MARKERS\n", fp);
    }
    if (pdinfo->time_series == STACKED_TIME_SERIES) 
	fprintf(fp, "PANEL2\n");
    else if (pdinfo->time_series == STACKED_CROSS_SECTION) 
	fprintf(fp, "PANEL3\n");
    
    if (fp != NULL) fclose(fp);
    return 0;
}

/**
 * get_precision:
 * @x: data vector.
 * @n: length of x.
 * @placemax: maximum number of decimal places to try.
 *
 * Find the number of decimal places required to represent a given
 * data series uniformly.
 * 
 * Returns: the required number of decimal places.
 *
 */

int get_precision (double *x, int n, int placemax)
{
    int i, p, pmax = 0;
    char *s, numstr[48];

    for (i=0; i<n; i++) {
	if (na(x[i])) continue;
	/* escape clause: numbers are too big or too small for
	   this treatment */
	if (x[i] < 1.0e-6 || x[i] > 1.0e+8) {
	    return PMAX_NOT_AVAILABLE;
	}
	p = placemax;
	sprintf(numstr, "%.*f", p, x[i]);
	s = numstr + strlen(numstr) - 1;
	while (*s-- == '0') p--;
	if (p > pmax) pmax = p;
    }

    return pmax;
}

/**
 * write_data:
 * @fname: name of file to write.
 * @list: list of variables to write.
 * @Z: data matrix.
 * @pdinfo: data information struct.
 * @opt: code for format in which to write the data (see #data_options).
 * @ppaths: pointer to paths information (should be NULL when not
 * called from gui).
 * 
 * Write out a data file containing the values of the given set
 * of variables.
 * 
 * Returns: 0 on successful completion, non-zero on error.
 * 
 */

int write_data (const char *fname, const int *list, 
		double **Z, const DATAINFO *pdinfo, 
		int opt, PATHS *ppaths)
{
    int i = 0, t, l0, n = pdinfo->n;
    char datfile[MAXLEN], hdrfile[MAXLEN], lblfile[MAXLEN];
    FILE *fp = NULL;
    int *pmax = NULL, tsamp = pdinfo->t2 - pdinfo->t1 + 1;
    double xx;

    *gretl_errmsg = 0;

    if (list == NULL) return 1;

    l0 = list[0];
    if (l0 == 0) return 1;

    if (opt == 0 || opt == GRETL_DATA_GZIPPED) {
	return write_xmldata(fname, list, Z, pdinfo, opt, ppaths);
    }

    if (opt == GRETL_DATA_CSV && pdinfo->delim == ',' && 
	',' == pdinfo->decpoint) {
	sprintf(gretl_errmsg, _("You can't use the same character for "
				"the column delimiter and the decimal point"));
	return 1;
    }

    strcpy(datfile, fname);

    if (opt == GRETL_DATA_R && pdinfo->time_series == TIME_SERIES) {
	opt = GRETL_DATA_R_ALT;
    }

    /* write header and label files if not exporting to other formats */
    if (opt != GRETL_DATA_R && opt != GRETL_DATA_R_ALT && 
	opt != GRETL_DATA_CSV && opt != GRETL_DATA_OCTAVE) {
	if (!has_gz_suffix(datfile)) {
	    switch_ext(hdrfile, datfile, "hdr");
	    switch_ext(lblfile, datfile, "lbl");
	} else {
	    gz_switch_ext(hdrfile, datfile, "hdr");
	    gz_switch_ext(lblfile, datfile, "lbl");
	}
	if (writehdr(hdrfile, list, pdinfo, opt)) {
	    fprintf(stderr, I_("Write of header file failed"));
	    return 1;
	}
	if (writelbl(lblfile, list, pdinfo)) {
	    fprintf(stderr, I_("Write of labels file failed"));
	    return 1;
	}
    }

    /* open files, other than for gzipped output */
    if (opt == GRETL_DATA_FLOAT || opt == GRETL_DATA_DOUBLE) {
	fp = fopen(datfile, "wb");
    } else {
	fp = fopen(datfile, "w");
    }
    if (fp == NULL) return 1;

    if (opt == GRETL_DATA_FLOAT) { /* single-precision binary */
	float x;

	for (i=1; i<=l0; i++) {
	    for (t=0; t<n; t++) {
		x = (float) (pdinfo->vector[list[i]])? 
			     Z[list[i]][t] : Z[list[i]][0];
		fwrite(&x, sizeof(float), 1, fp);
	    }
	}
    }
    else if (opt == GRETL_DATA_DOUBLE) { /* double-precision binary */
	for (i=1; i<=l0; i++) {
	    if (pdinfo->vector[list[i]])
		fwrite(&Z[list[i]][0], sizeof(double), n, fp);
	    else {
		for (t=0; t<n; t++) 
		    fwrite(&Z[list[i]][0], sizeof(double), 1, fp);
	    }
	}
    }

    if (opt == GRETL_DATA_CSV || opt == GRETL_DATA_OCTAVE || 
	GRETL_DATA_R || opt == GRETL_DATA_TRAD) { 
	/* an ASCII variant of some sort */
	pmax = malloc(l0 * sizeof *pmax);
	if (pmax == NULL) return 1;
	for (i=1; i<=l0; i++) {
	    if (pdinfo->vector[list[i]]) {
		pmax[i-1] = get_precision(&Z[list[i]][pdinfo->t1], tsamp, 10);
	    } else {
		pmax[i-1] = SCALAR_DIGITS;
	    }
	}	
    }

#ifdef ENABLE_NLS
    if (opt == GRETL_DATA_CSV && pdinfo->decpoint == ',') ;
    else setlocale(LC_NUMERIC, "C");
#endif

    if (opt == GRETL_DATA_TRAD) { /* plain ASCII */
	for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	    if (pdinfo->markers && pdinfo->S != NULL) {
		fprintf(fp, "%s ", pdinfo->S[t]);
	    }
	    for (i=1; i<=l0; i++) {
		if (na(Z[list[i]][t])) {
		    fprintf(fp, "-999 ");
		} else if (pmax[i-1] == PMAX_NOT_AVAILABLE) {
		    fprintf(fp, "%.12g ", 
			    (pdinfo->vector[list[i]])? 
			    Z[list[i]][t] : Z[list[i]][0]);
		} else {
		    fprintf(fp, "%.*f ", pmax[i-1], 
			    (pdinfo->vector[list[i]])? 
			    Z[list[i]][t] : Z[list[i]][0]);
		}
	    }
	    fputs("\n", fp);
	}
	fputs("\n", fp);
    }
    else if (opt == GRETL_DATA_CSV || opt == GRETL_DATA_R) { 
	/* export CSV or GNU R (dataframe) */
	char delim;
	
	if (opt == GRETL_DATA_CSV) delim = pdinfo->delim;
	else delim = ' ';

	/* variable names */
	if (opt == GRETL_DATA_CSV) fprintf(fp, "obs%c", delim);
	for (i=1; i<l0; i++) {
	    fprintf(fp, "%s%c", pdinfo->varname[list[i]], delim);
	}
	fprintf(fp, "%s\n", pdinfo->varname[list[l0]]);
	
	for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	    if (pdinfo->S != NULL) {
		fprintf(fp, "%s%c", pdinfo->S[t], delim);
	    } else {
		char tmp[OBSLEN];

		ntodate(tmp, t, pdinfo);
		fprintf(fp, "\"%s\"%c", tmp, delim);
	    }
	    for (i=1; i<=l0; i++) { 
		xx = (pdinfo->vector[list[i]])? Z[list[i]][t] : Z[list[i]][0];
		if (na(xx)) {
		    fprintf(fp, "NA");
		} else if (pmax[i-1] == PMAX_NOT_AVAILABLE) {
		    fprintf(fp, "%.12g", xx);;
		} else {
		    fprintf(fp, "%.*f", pmax[i-1], xx);
		}
		if (i < l0) fputc(delim, fp);
		else fputc('\n', fp);
	    }
	}
	fputc('\n', fp);
    }
    else if (opt == GRETL_DATA_R_ALT && pdinfo->time_series == TIME_SERIES) {
	/* new (October, 2003) attempt at improved R time-series structure */
	char *p, datestr[OBSLEN];
	int subper = 1;

	fprintf(fp, "\"%s\" <- ts (t (matrix (data = c(\n", "gretldata");
	for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	    for (i=1; i<=l0; i++) {
		xx = (pdinfo->vector[list[i]])? Z[list[i]][t] : Z[list[i]][0];
		if (na(xx)) fputs("NA", fp);
		else fprintf(fp, "%g", xx);
		if (i == l0) {
		    if (t == pdinfo->t2) fputs("),\n", fp);
		    else fputs(" ,\n", fp);
		} else {
		    fputs(" , ", fp);
		}
	    }
	}
	ntodate(datestr, pdinfo->t1, pdinfo);
	p = strchr(datestr, ':');
	if (p != NULL) subper = atoi(p + 1);
	fprintf(fp, "nrow = %d, ncol = %d)), start = c(%d,%d), frequency = %d)\n",
		l0, pdinfo->t2 - pdinfo->t1 + 1, 
		atoi(datestr), subper, pdinfo->pd);
	fprintf(fp, "colnames(%s) <- c(", "gretldata");
	for (i=1; i<=l0; i++) {
	    fprintf(fp, "\"%s\"", pdinfo->varname[list[i]]);
	    if (i < l0) fputs(", ", fp);
	    else fputs(")\n", fp);
	}
    }
    else if (opt == GRETL_DATA_R_ALT) { 
	/* export GNU R (structure) */
	for (i=1; i<=l0; i++) {
	    fprintf(fp, "\"%s\" <-\n", pdinfo->varname[list[i]]);
	    fprintf(fp, "structure(c(");
	    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
		xx = (pdinfo->vector[list[i]])? Z[list[i]][t] : Z[list[i]][0];
		if (na(xx)) fprintf(fp, "NA");
		else fprintf(fp, "%g", xx);
		if (t < pdinfo->t2) fprintf(fp, ", "); 
		if (t > pdinfo->t1 && (t - pdinfo->t1) % 8 == 0 && t < pdinfo->t2) {
		    fputc('\n', fp);
		}
	    }
	    fputc(')', fp);
	    if (pdinfo->time_series == TIME_SERIES) 
		fprintf(fp, ",\n.Tsp = c(%f, %f, %d), class = \"ts\"",
			obs_float(pdinfo, 0), obs_float(pdinfo, 1), 
			pdinfo->pd);
	    fprintf(fp, ")\n");
	}
    }
    else if (opt == GRETL_DATA_OCTAVE) { /* export GNU Octave */
	/* write out info on dependent variable */
	fprintf(fp, "# name: %s\n# type: matrix\n# rows: %d\n# columns: 1\n", 
		pdinfo->varname[list[1]], n);
	/* write out column of values of dep. var. */
	for (t=pdinfo->t1; t<=pdinfo->t2; t++) { 
	    if (pmax[0] == PMAX_NOT_AVAILABLE) {
		fprintf(fp, "%.12g\n", Z[list[1]][t]);
	    } else {
		fprintf(fp, "%.*f\n", pmax[0], 
			(pdinfo->vector[list[i]])? 
			Z[list[1]][t] : Z[list[1]][0]);
	    }
	}
	/* write out info for indep vars matrix */
	fprintf(fp, "# name: X\n# type: matrix\n# rows: %d\n# columns: %d\n", 
		n, list[0] - 1);
	/* write out indep. var. matrix */
	for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	    for (i=2; i<=list[0]; i++) {
		if (pmax[i-1] == PMAX_NOT_AVAILABLE) {
		    fprintf(fp, "%.12g ", Z[list[i]][t]);
		} else {
		    fprintf(fp, "%.*f ", pmax[i-1], 
			    (pdinfo->vector[list[i]])? 
			    Z[list[i]][t] : Z[list[i]][0]);
		}
	    }
	    fputc('\n', fp);
	}
	fputc('\n', fp);
    }

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif

    if (pmax) free(pmax);
    if (fp != NULL) fclose(fp);

    return 0;
}

static void type_string (char *str, const DATAINFO *pdinfo)
{
    if (dataset_is_time_series(pdinfo)) {
	strcpy(str, _("time series"));
    } else if (dataset_is_panel(pdinfo)) {
        strcpy(str, _("panel"));
    } else {
        strcpy(str, _("undated"));
    }
}

static void pd_string (char *str, const DATAINFO *pdinfo)
{
    switch (pdinfo->pd) {
    case 1:
	strcpy(str, _("annual")); break;
    case 4:
	strcpy(str, _("quarterly")); break;
    case 12:
	strcpy(str, _("monthly")); break;
    case 24:
	strcpy(str, _("hourly")); break;
    case 52:
	strcpy(str, _("weekly")); break;
    case 5:
	strcpy(str, _("daily")); break;
    case 7:
	strcpy(str, _("daily")); break;
    default:
	strcpy(str, _("unknown")); break;
    }
}

/**
 * data_report:
 * @pdinfo: data information struct.
 * @ppaths: path information struct.
 * @prn: gretl printing struct.
 * 
 * Write out a summary of the content of the current data set.
 * 
 * Returns: 0 on successful completion, non-zero on error.
 * 
 */

int data_report (const DATAINFO *pdinfo, PATHS *ppaths, PRN *prn)
{
    char startdate[OBSLEN], enddate[OBSLEN], tmp[MAXLEN];
    time_t prntime = time(NULL);
    int i;

    ntodate(startdate, 0, pdinfo);
    ntodate(enddate, pdinfo->n - 1, pdinfo);

    sprintf(tmp, _("Data file %s\nas of"), 
	    strlen(ppaths->datfile)? ppaths->datfile : _("(unsaved)"));

    pprintf(prn, "%s %s\n\n", tmp, print_time(&prntime));

    if (pdinfo->descrip != NULL && strlen(pdinfo->descrip)) {
	pprintf(prn, "%s:\n\n", _("Description"));
	pprintf(prn, "%s\n\n", pdinfo->descrip);
    }

    type_string(tmp, pdinfo);
    pprintf(prn, "%s: %s\n", _("Type of data"), tmp);
    
    if (dataset_is_time_series(pdinfo)) {
	pd_string(tmp, pdinfo);
	pprintf(prn, "%s: %s\n", _("Frequency"), tmp);
    }	

    pprintf(prn, "%s: %s - %s (n = %d)\n\n", _("Range"),
	    startdate, enddate, pdinfo->n);

    pprintf(prn, "%s:\n\n", _("Listing of variables"));

    for (i=1; i<pdinfo->v; i++) {
	pprintf(prn, "%*s  %s\n", VNAMELEN, pdinfo->varname[i], 
		VARLABEL(pdinfo, i));
    }

    return 0;
}

/* ................................................. */

static double obs_float (const DATAINFO *pdinfo, int end)
{
    double xx, xx2 = 0.;
    int i, x1, x2 = 0;

    if (end) {
	xx = obs_str_to_double(pdinfo->endobs);
	if ((i = haschar(':', pdinfo->endobs)) > 0)
	   x2 = atoi(pdinfo->endobs + i + 1) - 1;
    } else {
	xx = obs_str_to_double(pdinfo->stobs);
	if ((i = haschar(':', pdinfo->stobs)) > 0)
	   x2 = atoi(pdinfo->stobs + i + 1) - 1;
    }
    x1 = (int) xx;
    if (x2 > 0)
	xx2 = (double) x2 / pdinfo->pd;
    
    return (double) x1 + xx2;
}

/* ................................................. */

static int readlbl (const char *lblfile, DATAINFO *pdinfo)
     /* read data "labels" from file */
{
    FILE * fp;
    char line[MAXLEN], *label, varname[VNAMELEN];
    int v;
    
    *gretl_errmsg = '\0';

    fp = fopen(lblfile, "r");
    if (fp == NULL) return 0;
    while (1) {
        if (fgets(line, MAXLEN-1, fp) == NULL) {
            fclose(fp);
            return 0;
        }
        if (sscanf(line, "%s", varname) != 1) {
            fclose(fp);
	    sprintf(gretl_errmsg, _("Bad data label in %s"), lblfile); 
            return 0;
        }
        label = line + strlen(varname);
        if (top_n_tail(label) == E_ALLOC) {
            fclose(fp);
            return E_ALLOC;
        }
	v = varindex(pdinfo, varname);
	if (v < pdinfo->v) {
	    strcpy(VARLABEL(pdinfo, v), label);
	} else {
	    fprintf(stderr, I_("extraneous label for var '%s'\n"), varname);
	}
    }

    if (fp != NULL) fclose(fp);

    return 0;
}

/* ................................................ */

static int writelbl (const char *lblfile, const int *list, 
		     const DATAINFO *pdinfo)
{
    FILE *fp;
    int i, lblcount = 0;

    for (i=1; i<=list[0]; i++) {
	if (list[i] == 0) continue;
	if (strlen(VARLABEL(pdinfo, list[i])) > 2) {
	    lblcount++;
	    break;
	}
    }
    if (lblcount == 0) return 0;

    fp = fopen(lblfile, "w");
    if (fp == NULL) return 1;

    /* spit out varnames and labels (if filled out) */
    for (i=1; i<=list[0]; i++) {
	if (list[i] == 0) continue;
	if (strlen(VARLABEL(pdinfo, list[i])) > 2) {
	    fprintf(fp, "%s %s\n", pdinfo->varname[list[i]],
		    VARLABEL(pdinfo, list[i]));
	}
    }
    
    if (fp != NULL) fclose(fp);

    return 0;
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

    if (fname == NULL || *fname == '\0') return 0;

    fp = fopen(fname, "rb");
    if (fp == NULL) return 0;

    if (fgetc(fp) == 037 && fgetc(fp) == 0213) 
	gz = 1;

    fclose(fp);
    return gz;
}

/**
 * has_gz_suffix:
 * @fname: filename to examine.
 * 
 * Determine if the given filename ends with ".gz".
 * 
 * Returns: 1 in case of a ".gz" suffix, otherwise 0.
 * 
 */

int has_gz_suffix (const char *fname)
{
    size_t n = strlen(fname);
	
    if (n < 4 || strncmp(fname + n - 3, ".gz", 3))
	return 0;
    else
	return 1;
}

/**
 * gz_switch_ext:
 * @targ: target or "output" filename (must be pre-allocated).
 * @src: "source or "input" filename.
 * @ext: suffix to add to filename.
 * 
 * Copy @src filename to @targ, without the existing suffix (if any),
 * and adding the supplied extension or suffix.
 * 
 */

void gz_switch_ext (char *targ, char *src, char *ext)
{
    size_t i = dotpos(src), j = slashpos(src), k;

    strcpy(targ, src);
    targ[i] = '\0';
    k = dotpos(targ);
    if (j > 0 && k < strlen(targ) && k > j) i = k;
    targ[i] = '.';
    targ[i + 1] = '\0';
    strcat(targ, ext);
}

/* ................................................ */

static void try_gdt (char *fname)
{
    char *suff;

    if (fname != NULL) {
	suff = strrchr(fname, '.');
	if (suff != NULL && !strcmp(suff, ".dat")) {
	    strcpy(suff, ".gdt");
	} else {
	    strcat(fname, ".gdt");
	}
    }
}

/**
 * get_data:
 * @pZ: pointer to data set.
 * @pdinfo: data information struct.
 * @datfile: name of file to try.
 * @ppaths: path information struct.
 * @data_status: indicator for whether a data file is currently open
 * in gretl's work space (not-0) or not (0).
 * @prn: where messages should be written.
 * 
 * Read data from file into gretl's work space, allocating space as
 * required.
 * 
 * Returns: 0 on successful completion, non-zero otherwise.
 *
 */

int get_data (double ***pZ, DATAINFO *pdinfo, char *datfile, PATHS *ppaths, 
	      int data_status, PRN *prn) 
{

    FILE *dat = NULL;
    gzFile fz = NULL;
    int err, gzsuff = 0, add_gdt = 0;
    char hdrfile[MAXLEN], lblfile[MAXLEN];

    *gretl_errmsg = '\0';

    /* get filenames organized */
    *hdrfile = '\0';
    gzsuff = has_gz_suffix(datfile);

    if (addpath(datfile, ppaths, 0) == NULL) { /* not found yet */
	char tryfile[MAXLEN];
	int found = 0;

	/* try using the .gdt suffix? */
	*tryfile = '\0';
	strncat(tryfile, datfile, MAXLEN-1);
	try_gdt(tryfile); 
	found = (addpath(tryfile, ppaths, 0) != NULL);
	if (found) add_gdt = 1;

	/* or maybe the file is gzipped but lacks a .gz extension? */
	if (!found && !gzsuff) { 
	    sprintf(tryfile, "%s.gz", datfile);
	    if (addpath(tryfile, ppaths, 0) != NULL) {
		gzsuff = 1;
		found = 1;
	    }
	}
	if (!found) return E_FOPEN;
	else strcpy(datfile, tryfile);
    }

    /* catch XML files that have strayed in here? */
    if (add_gdt && xmlfile(datfile)) {
	return get_xmldata(pZ, pdinfo, datfile, ppaths, 
			   data_status, prn, 0);
    }
	
    if (!gzsuff) {
	switch_ext(hdrfile, datfile, "hdr");
	switch_ext(lblfile, datfile, "lbl");
    } else {
	gz_switch_ext(hdrfile, datfile, "hdr");
	gz_switch_ext(lblfile, datfile, "lbl");
    }

    /* clear any existing data info */
    if (data_status) clear_datainfo(pdinfo, CLEAR_FULL);

    /* read data header file */
    err = readhdr(hdrfile, pdinfo);
    if (err) {
	return err;
    } else { 
	pprintf(prn, I_("\nReading header file %s\n"), hdrfile);
    }

    /* deal with case where first col. of data file contains
       "marker" strings */
    pdinfo->S = NULL;
    if (pdinfo->markers && dataset_allocate_markers(pdinfo)) {
	return E_ALLOC; 
    }
    
    /* allocate dataset */
    if (prepZ(pZ, pdinfo)) return E_ALLOC;

    /* Invoke data (Z) reading function */
    if (gzsuff) {
	fz = gzopen(datfile, "rb");
	if (fz == NULL) return E_FOPEN;
    } else {
	if (pdinfo->bin) {
	    dat = fopen(datfile, "rb");
	} else {
	    dat = fopen(datfile, "r");
	}
	if (dat == NULL) return E_FOPEN;
    }

    /* print out basic info from the files read */
    pprintf(prn, I_("periodicity: %d, maxobs: %d, "
	   "observations range: %s-%s\n"), pdinfo->pd, pdinfo->n,
	   pdinfo->stobs, pdinfo->endobs);

    pputs(prn, I_("\nReading "));
    pputs(prn, (pdinfo->time_series == TIME_SERIES) ? 
	    I_("time-series") : _("cross-sectional"));
    pputs(prn, I_(" datafile"));
    if (strlen(datfile) > 40) pputs(prn, "\n");
    pprintf(prn, " %s\n\n", datfile);

    if (gzsuff) {
	err = gz_readdata(fz, pdinfo, *pZ); 
	gzclose(fz);
    } else {
	err = readdata(dat, pdinfo, *pZ); 
	fclose(dat);
    }

    if (err) return err;

    /* Set sample range to entire length of dataset by default */
    pdinfo->t1 = 0; 
    pdinfo->t2 = pdinfo->n - 1;

    strcpy(ppaths->datfile, datfile);

    err = readlbl(lblfile, pdinfo);
    if (err) return err;    

    return 0;
}

/**
 * open_nulldata:
 * @pZ: pointer to data set.
 * @pdinfo: data information struct.
 * @data_status: indicator for whether a data file is currently open
 * in gretl's work space (1) or not (0).
 * @length: desired length of data series.
 * @prn: gretl printing struct.
 * 
 * Create an empty "dummy" data set, suitable for Monte Carlo simulations.
 * 
 * Returns: 0 on successful completion, non-zero otherwise.
 *
 */

int open_nulldata (double ***pZ, DATAINFO *pdinfo, 
		   int data_status, int length,
		   PRN *prn) 
{
    int t;

    /* clear any existing data info */
    if (data_status) clear_datainfo(pdinfo, CLEAR_FULL);

    /* dummy up the data info */
    pdinfo->n = length;
    pdinfo->v = 2;
    dataset_dates_defaults(pdinfo);

    if (dataset_allocate_varnames(pdinfo)) return E_ALLOC;

    /* no observation markers */
    pdinfo->markers = NO_MARKERS;
    pdinfo->S = NULL; 

    /* no descriptive comments */
    pdinfo->descrip = NULL;

    /* allocate dataset */
    if (prepZ(pZ, pdinfo)) return E_ALLOC;

    /* add an index var */
    strcpy(pdinfo->varname[1], "index");
    strcpy(VARLABEL(pdinfo, 1), _("index variable"));
    for (t=0; t<pdinfo->n; t++) (*pZ)[1][t] = (double) (t + 1);

    /* print out basic info */
    pprintf(prn, I_("periodicity: %d, maxobs: %d, "
	   "observations range: %s-%s\n"), pdinfo->pd, pdinfo->n,
	   pdinfo->stobs, pdinfo->endobs);

    /* Set sample range to entire length of data-set by default */
    pdinfo->t1 = 0; 
    pdinfo->t2 = pdinfo->n - 1;

    return 0;
}

/* .......................................................... */

#ifdef PAD_DAILY_DATA

static int daily_data_resize (double ***pZ, DATAINFO *pdinfo, 
			      char *missvec)
{
    int i, j, t;
    double **tmp;

    tmp = malloc(pdinfo->v * sizeof *tmp);
    if (tmp == NULL) return 1;

    for (i=0; i<pdinfo->v; i++) {
	tmp[i] = malloc(pdinfo->n * sizeof **tmp);
	if (tmp[i] == NULL) {
	    for (j=0; j<i; j++) {
		free(tmp[j]);
	    }
	    free(tmp);
	    return 1;
	}
    }

    j = 0;
    for (t=0; t<pdinfo->n; t++) {
	int missing = missvec[t];

	tmp[0][t] = 1.0;
	for (i=1; i<pdinfo->v; i++) {
	    if (missing) tmp[i][t] = NADBL;
	    else tmp[i][t] = (*pZ)[i][j];
	}
	if (!missing) j++;
    }

    free_Z(*pZ, pdinfo);
    *pZ = tmp;

    ntodate(pdinfo->stobs, 0, pdinfo);
    ntodate(pdinfo->endobs, pdinfo->n - 1, pdinfo);
    pdinfo->t2 = pdinfo->n - 1;

    return 0;
}

#endif /* PAD_DAILY_DATA */

static int check_daily_dates (DATAINFO *pdinfo, char **pmiss)
{
    int fulln = 0, n, t;
    int pd = pdinfo->pd;
    double sd0 = pdinfo->sd0;
    long ed1, ed2;
    int nmiss = 0, err = 0;
    char *missvec = NULL;
    
    pdinfo->pd = 5; /* let's try it */
    pdinfo->time_series = TIME_SERIES;

    ed1 = get_epoch_day(pdinfo->S[0]);
    if (ed1 < 0) err = 1;

    if (!err) {
	ed2 = get_epoch_day(pdinfo->S[pdinfo->n - 1]);
	if (ed2 <= ed1) err = 1;
	else pdinfo->sd0 = ed1;
    }

    if (!err) {
	int n1 = daily_obs_number(pdinfo->S[0], pdinfo);
	int n2 = daily_obs_number(pdinfo->S[pdinfo->n - 1], pdinfo);

	fulln = n2 - n1 + 1;
	if (pdinfo->n > fulln) {
	    pdinfo->pd = 7;
	    /* FIXME need to do more here */
	} else {
	    nmiss = fulln - pdinfo->n;
	    fprintf(stderr, "Observations: %d; days in sample: %d\n", 
		    pdinfo->n, fulln);
	    fprintf(stderr, "Missing daily observations: %d\n", nmiss);
	}
    }

#ifdef PAD_DAILY_DATA
    if (!err && nmiss > 0) {
	missvec = malloc(fulln * sizeof *missvec);
	if (missvec == NULL) err = 1;
	else {
	    for (t=0; t<fulln; t++) {
		missvec[t] = 1;
	    }
	}
    }
#endif

    for (t=0; t<pdinfo->n && !err; t++) {
	n = daily_obs_number(pdinfo->S[t], pdinfo);
	if (n < t) {
	    fprintf(stderr, "Error: n = %d < t = %d\n", n, t);
	    err = 1;
	} 
	else if (n > fulln - 1) {
	    fprintf(stderr, "Error: n = %d >= fulln = %d\n", n, fulln);
	    err = 1;
	}
#ifdef PAD_DAILY_DATA 
	else {
	    missvec[n] = 0;
	}
#endif
    }

    if (err) {
	pdinfo->pd = pd;
	pdinfo->sd0 = sd0;
	pdinfo->time_series = 0;
	if (missvec) free(missvec);
    } else if (pmiss != NULL && missvec != NULL) {
	*pmiss = missvec;
	for (t=0; t<pdinfo->n; t++) {
	    free(pdinfo->S[t]);
	}
	free(pdinfo->S);
	pdinfo->S = NULL;
	pdinfo->markers = NO_MARKERS;
	pdinfo->n = fulln;
    } else {
	strcpy(pdinfo->stobs, pdinfo->S[0]);
	strcpy(pdinfo->endobs, pdinfo->S[pdinfo->n - 1]);
	pdinfo->t2 = pdinfo->n - 1;
	if (nmiss > 0) pdinfo->markers = DAILY_DATE_STRINGS;
    }

    if (!err) return pdinfo->pd;
    else return -1;
}

/* ......................................................... */

#if 0
static int labels_all_numeric (DATAINFO *pdinfo)
{
    int t, ret = 1;

    for (t=0; t<pdinfo->n; t++) {
	if (!_isnumber(pdinfo->S[t])) {
	    ret = 0;
	    break;
	}
    }

    return  ret;
}
#endif

static int complete_year_labels (DATAINFO *pdinfo)
{
    int t, yr, yrbak = atoi(pdinfo->S[0]);
    int ret = 1;

    for (t=1; t<pdinfo->n; t++) {
	yr = atoi(pdinfo->S[t]);
	if (yr != yrbak + 1) {
	    ret = 0;
	    break;
	}
	yrbak = yr;
    }

    return  ret;
}

/* ......................................................... */

static int test_label (DATAINFO *pdinfo, char **missvec, PRN *prn)
     /* attempt to parse csv row labels as dates.  Return -1 if
	this doesn't work out, 0 if the labels seem to be just
	integer observation numbers, else return the inferred data
	frequency */
{
    int n1, n2, try;
    char year[5], subper[3], endobs[OBSLEN];
    char lbl1[OBSLEN], lbl2[OBSLEN];

    *lbl1 = *lbl2 = 0;
    strncat(lbl1, pdinfo->S[0], 8);
    strncat(lbl2, pdinfo->S[pdinfo->n - 1], 8);
    n1 = strlen(lbl1);
    n2 = strlen(lbl2);

    pprintf(prn, M_("   first row label \"%s\", last label \"%s\"\n"), 
	   lbl1, lbl2);

    /* are the labels (probably) just 1, 2, 3 etc.? */
    sprintf(endobs, "%d", pdinfo->n);
    if (strcmp(lbl1, "1") == 0 && strcmp(lbl2, endobs) == 0)
	return 0;

    if (n1 != n2) {
	pputs(prn, M_("   label strings can't be consistent dates\n"));
	return -1;
    }

    pputs(prn, M_("trying to parse row labels as dates...\n"));

    /* daily data? */
    if (n1 == 8) {
	int yr1, mon1, day1;
	int yr2, mon2, day2;

	if (sscanf(lbl1, "%d/%d/%d", &yr1, &mon1, &day1) == 3 &&
	    sscanf(lbl2, "%d/%d/%d", &yr2, &mon2, &day2) == 3) {
	    if (yr2 >= yr1 && 
		mon1 > 0 && mon1 < 13 &&
		mon2 > 0 && mon2 < 13 && 
		day1 > 0 && day1 < 32 &&
		day2 > 0 && day2 < 32) {
		/* looks promising for calendar dates */
		return check_daily_dates(pdinfo, missvec);
	    }
	} else {
	    return -1;
	}
    }

    else if (n1 >= 4) {
	*year = 0;
	if (isdigit((unsigned char) lbl1[0]) &&
	    isdigit((unsigned char) lbl1[1]) &&
	    isdigit((unsigned char) lbl1[2]) && 
	    isdigit((unsigned char) lbl1[3])) {
	    strncat(year, lbl1, 4);
	    try = atoi(year);
	    if (try > 0 && try < 3000) {
		pprintf(prn, M_("   %s: probably a year... "), year);
	    } else {
		pprintf(prn, M_("   %s: out of bounds for a year?\n"), year);
	    }
	    if (n1 == 5) {
		pputs(prn, M_("   but I can't make sense of the extra bit\n"));
		return -1;
	    }
	    if (n1 == 4) {
		pputs(prn, M_("and just a year\n"));
		if (complete_year_labels(pdinfo)) {
		    strcpy(pdinfo->stobs, year);
		    pdinfo->sd0 = atof(pdinfo->stobs);
		    strcpy(pdinfo->endobs, lbl2);
		    pdinfo->pd = 1;
		    return 1;
		} else {
		    pputs(prn, M_("   but the dates are not complete and consistent\n"));
		    return -1;
		}
	    }
	    if (lbl1[4] == '.' || lbl1[4] == ':' || lbl1[4] == 'Q') {
		strcpy(subper, lbl1 + 5);
		if (n1 == 6) {
		    pprintf(prn, M_("quarter %s?\n"), subper);
		    sprintf(pdinfo->stobs, "%s:%s", year, subper);
		    pdinfo->sd0 = obs_str_to_double(pdinfo->stobs);
		    *pdinfo->endobs = 0;
		    strncat(pdinfo->endobs, lbl2, 4);
		    strcat(pdinfo->endobs, ":");
		    strcat(pdinfo->endobs, lbl2 + 5);
		    pdinfo->pd = 4;
		    return 4;
		}
		if (n1 == 7) {
		    pprintf(prn, M_("month %s?\n"), subper);
		    sprintf(pdinfo->stobs, "%s:%s", year, subper);
		    pdinfo->sd0 = obs_str_to_double(pdinfo->stobs);
		    *pdinfo->endobs = 0;
		    strncat(pdinfo->endobs, lbl2, 4);
		    strcat(pdinfo->endobs, ":");
		    strcat(pdinfo->endobs, lbl2 + 5);
		    pdinfo->pd = 12;
		    return 12;
		}
	    }
	} else pputs(prn, M_("   definitely not a four-digit year\n"));
    }

    return -1;
}

static int var_overlap (const DATAINFO *pdinfo, const DATAINFO *addinfo)
{
    int i, j;
    int match = 0;

    for (i=1; i<pdinfo->v; i++) {
	for (j=1; j<addinfo->v; j++) {
	    if (!strcmp(pdinfo->varname[i], addinfo->varname[j])) {
		match++;
		break;
	    }
	}
    }

    return match;
}

static int new_data_offset_ok (const DATAINFO *pdinfo,
			       const DATAINFO *addinfo,
			       int *offset)
{
    int sd0, sd1;
    int ed0, ed1;
    int ok = 0;

    sd0 = dateton(pdinfo->stobs, pdinfo);
    sd1 = dateton(addinfo->stobs, pdinfo);
    ed0 = dateton(pdinfo->endobs, pdinfo);
    ed1 = dateton(addinfo->endobs, pdinfo);

    /* case: starting obs the same */
    if (sd1 == sd0) {
	*offset = 0;
	ok = 1;
    }

    /* case: new data start later */
    else if (sd1 > sd0) {
	if (sd1 < ed0) {
	    /* but there's some overlap */
	    *offset = sd1 - sd0;
	    ok = 1;
	}
    }

    /* case: new data start earlier */
    else if (sd1 < sd0) {
	if (ed1 > sd0) {
	    /* but there's some overlap */
	    *offset = sd1 - sd0;
	    ok = 1;
	}
    }

    return ok;
}

/**
 * merge_data:
 * @pZ: pointer to data set.
 * @pdinfo: data information struct.
 * @addZ: new data set to be merged in.
 * @addinfo: data information associated with @addZ.
 * @prn: print struct to accept messages.
 * 
 * Attempt to merge the content of a newly opened data file into
 * gretl's current working data set.
 * 
 * Returns: 0 on successful completion, non-zero otherwise.
 *
 */

int merge_data (double ***pZ, DATAINFO *pdinfo,
		double **addZ, DATAINFO *addinfo,
		PRN *prn)
{
    int err = 0, addrows = 0, addcols = 0;
    int offset = 0;
    int match_obs = 0, match_vars = 0;

    /* first check for conformability */

    if (pdinfo->pd != addinfo->pd) {
	pprintf(prn, _("Data frequency does not match\n"));
	err = 1;
    }

    if (!err) {
	match_obs = new_data_offset_ok(pdinfo, addinfo, &offset);
	match_vars = var_overlap(pdinfo, addinfo);
    }

    if (!err && !match_obs && !match_vars) {
	pputs(prn, _("New data not conformable for appending\n"));
	err = 1;
    }

    else if (!err && match_obs && !match_vars) {
	/* adding variable(s) */
	addcols = 1;
    }

    else if (!err && match_vars) {
	/* adding observations */
	addrows = 1;
    }

    if (addrows && !addcols) {
	if (pdinfo->time_series && 
	    dateton(addinfo->stobs, pdinfo) != pdinfo->n) {
	    pputs(prn, _("Starting point of new data does not fit\n"));
	    err = 1;
	}
	else if (pdinfo->markers != addinfo->markers) {
	    pputs(prn, _("Inconsistency in observation markers\n"));
	    err = 1;
	}
	if (err) addrows = 0;
    }

    /* if checks are passed, try merging the data */

   if (!err && addcols) { 
       int orig_vars = pdinfo->v;
       int i, t, nvars = pdinfo->v + addinfo->v - 1;

       if (dataset_add_vars(addinfo->v - 1, pZ, pdinfo)) {
	   pputs(prn, _("Out of memory adding data\n"));
	   err = 1;
       }

       for (i=orig_vars; i<nvars && !err; i++) {
	   strcpy(pdinfo->varname[i], addinfo->varname[i - orig_vars + 1]);
	   if (offset >= 0) {
	       for (t=0; t<pdinfo->n; t++) {
		   if (t < offset) {
		       (*pZ)[i][t] = NADBL;
		   } else {
		       (*pZ)[i][t] = addZ[i - orig_vars + 1][t - offset];
		   }
	       }
	   } else {
	       /* negative offset: new data start earlier */
	       for (t=0; t<pdinfo->n; t++) {
		   if (t < addinfo->n + offset) {
		       (*pZ)[i][t] = addZ[i - orig_vars + 1][t - offset];
		   } else {
		       (*pZ)[i][t] = NADBL;
		   }
	       }
	   }
       }
   }

   else if (!err && addrows) { 
       int i, t, tnew = pdinfo->n + addinfo->n;
       double *xx;

       if (pdinfo->markers) {
	   char **S = realloc(pdinfo->S, tnew * sizeof *S);
	   
	   if (S == NULL) {
	       err = 1;
	   } else {
	       for (t=pdinfo->n; t<tnew && !err; t++) {
		   S[t] = malloc(OBSLEN);
		   if (S[t] == NULL) err = 1;
		   else strcpy(S[t], addinfo->S[t - pdinfo->n]);
	       }
	       pdinfo->S = S;
	   }
       }

       for (i=0; i<pdinfo->v && !err; i++) {
	   xx = realloc((*pZ)[i], tnew * sizeof *xx);
	   if (xx == NULL) {
	       err = 1;
	       break;
	   } else {
	       for (t=pdinfo->n; t<tnew; t++) {
		   xx[t] = addZ[i][t - pdinfo->n];
	       }
	       (*pZ)[i] = xx;
	   }
       }

       if (err) { 
	   pputs(prn, _("Out of memory adding data\n"));
       } else {
	   pdinfo->n = tnew;
	   ntodate(pdinfo->endobs, tnew - 1, pdinfo);
	   pdinfo->t2 = pdinfo->n - 1;
       }
   }

   if (!err && (addcols || addrows)) {
       pputs(prn, _("Data appended OK\n"));
   }

   free_Z(addZ, addinfo);
   clear_datainfo(addinfo, CLEAR_FULL);

   return err;
}

/* ......................................................... */

/* The function below checks for the maximum line length in the given
   file.  It also checks for extraneous binary data (the file is 
   supposed to be plain text), and checks whether the 'delim'
   character is present in the file, on a non-comment line (where
   a comment line is one that starts with '#').
*/

static int get_max_line_length (FILE *fp, char delim, int *gotdelim, 
				int *gottab, PRN *prn)
{
    int c, cc = 0;
    int comment = 0, maxlen = 0;

    while (1) {
	c = fgetc(fp);
	if (c == '\n') {
	    if (cc > maxlen) maxlen = cc;
	    cc = 0;
	    continue;
	}
	if (c == EOF) break;
	if (!isspace((unsigned char) c) && !isprint((unsigned char) c) &&
	    !(c == CTRLZ)) {
	    pprintf(prn, M_("Binary data (%d) encountered: this is not a valid "
			   "text file\n"), c);
	    return -1;
	}
	if (cc == 0) {
	    if (c == '#') comment = 1;
	    else comment = 0;
	}
	if (!comment && *gottab == 0 && c == '\t') *gottab = 1;
	if (!comment && *gotdelim == 0 && c == delim) *gotdelim = 1;
	cc++;
    }

    if (maxlen == 0) {
	pprintf(prn, M_("Data file is empty\n"));
    }	

    return maxlen;
}

static int count_fields (const char *line, char delim)
{
    int cbak, nf = 0;
    const char *p = line;

    if (*p == delim && *p == ' ') p++;

    while (*p) {
	if (*p == delim) nf++;
	cbak = *p;
	p++;
	if (*p == '\0' && cbak == delim) {
	    nf--;
	}
    }

    return nf + 1;
}

static void compress_csv_line (char *line, char delim)
{
    int n = strlen(line);
    char *p = line + n - 1;

    if (*p == '\n') {
	*p = '\0';
	p--;
    }
    if (*p == '\r') *p = '\0';

    if (delim != ' ') {
	delchar(' ', line);
    } else {
	compress_spaces(line);
    }
    delchar('"', line);
}

static void check_first_field (const char *line, char delim, 
			       int *blank_1, int *obs_1, PRN *prn)
{
    *blank_1 = 0;
    *obs_1 = 0;
    
    if (delim != ' ' && *line == delim) {
	*blank_1 = 1;
    } else {
	char field1[16];
	int i = 0;

	if (delim == ' ' && *line == ' ') line++;

	while (*line && i < 15) {
	    if (*line == delim) break;
	    field1[i++] = *line++;
	}
	field1[i] = '\0';
	pprintf(prn, M_("   first field: '%s'\n"), field1);
	lower(field1);
	if (!strcmp(field1, "obs") || !strcmp(field1, "date") ||
	    !strcmp(field1, "year")) {
	    pputs(prn, M_("   seems to be observation label\n"));
	    *obs_1 = 1;
	}
    }
}

#define ISNA(s) (strcmp(s, "NA") == 0 || \
                 strcmp(s, "N.A.") == 0 || \
                 strcmp(s, "n.a.") == 0 || \
                 strcmp(s, "na") == 0 || \
                 strcmp(s, ".") == 0 || \
                 strncmp(s, "-999", 4) == 0)

static int csv_missval (const char *str, int k, int t, PRN *prn)
{
    int miss = 0;

    if (strlen(str) == 0) {
	pprintf(prn, M_("   the cell for variable %d, obs %d "
		       "is empty: treating as missing value\n"), 
		k, t);
	miss = 1;
    }

    if (ISNA(str)) {
	pprintf(prn, M_("   warning: missing value for variable "
		       "%d, obs %d\n"), k, t);
	miss = 1;
    }

    return miss;
}

static int dataset_add_obs (double ***pZ, DATAINFO *pdinfo)
{
    int i;

    for (i=0; i<pdinfo->v; i++) {
	double *tmp = realloc((*pZ)[i], (pdinfo->n + 1) * sizeof ***pZ);

	if (tmp != NULL) {
	    (*pZ)[i] = tmp;
	} else {
	    return 1;
	}
    }

    pdinfo->n += 1;

    return 0;
}

/**
 * import_csv:
 * @pZ: pointer to data set.
 * @pdinfo: data information struct.
 * @fname: name of CSV file.
 * @prn: gretl printing struct.
 * 
 * Open a Comma-Separated Values data file and read the data into
 * the current work space.
 * 
 * Returns: 0 on successful completion, non-zero otherwise.
 *
 */

int import_csv (double ***pZ, DATAINFO **ppdinfo, 
		const char *fname, PRN *prn)
{
    int ncols, chkcols;
    int gotdata = 0, gotdelim = 0, gottab = 0, markertest = -1;
    int i, k, t, blank_1 = 0, obs_1 = 0, maxlen;
    char csvstr[32];
    FILE *fp = NULL;
    DATAINFO *csvinfo = NULL;
    double **csvZ = NULL;
    char *line = NULL, *p = NULL;
    const char *msg = M_("\nPlease note:\n"
	"- The first row of the CSV file should contain the "
	"names of the variables.\n"
	"- The first column may optionally contain date "
	"strings or other 'markers':\n  in that case its row 1 entry "
	"should be blank, or should say 'obs' or 'date'.\n"
	"- The remainder of the file must be a rectangular "
	"array of data.\n");
    char delim = (*ppdinfo)->delim;
    int numcount, auto_name_vars = 0;
    char *missvec = NULL;

    check_for_console(prn);

    fp = fopen(fname, "r");
    if (fp == NULL) {
	pprintf(prn, M_("Couldn't open %s\n"), fname);
	goto csv_bailout;
    }

    csvinfo = datainfo_new();
    if (csvinfo == NULL) {
	fclose(fp);
	pputs(prn, M_("Out of memory\n"));
	goto csv_bailout;
    }
    csvinfo->delim = delim;

    pprintf(prn, "%s %s...\n", M_("parsing"), fname);

    /* get line length, also check for binary data */
    maxlen = get_max_line_length(fp, delim, &gotdelim, &gottab, prn);    
    if (maxlen <= 0) {
	goto csv_bailout;
    }
    if (!gotdelim) {
	if (gottab) delim = csvinfo->delim = '\t';
	else delim = csvinfo->delim = ' ';
    }

    pprintf(prn, M_("using delimiter '%c'\n"), delim);
    pprintf(prn, M_("   longest line: %d characters\n"), maxlen + 1);

    /* create buffer to hold lines */
    line = malloc(maxlen + 1);
    if (line == NULL) {
	pputs(prn, M_("Out of memory\n"));
	goto csv_bailout;
    }  
    
    rewind(fp);
    
    /* read lines, check for consistency in number of fields */
    chkcols = ncols = gotdata = 0;
    while (fgets(line, maxlen + 1, fp)) {
	/* skip comment lines */
	if (*line == '#') continue;
	/* skip blank lines */
	if (string_is_blank(line)) continue;
	csvinfo->n += 1;
	compress_csv_line(line, delim);
	if (!gotdata) {
	    /* scrutinize first "real" line */
	    check_first_field(line, delim, &blank_1, &obs_1, prn);
	    gotdata = 1;
	} 
	chkcols = count_fields(line, delim);
	if (ncols == 0) {
	    ncols = chkcols;
	    pprintf(prn, M_("   number of columns = %d\n"), ncols);	    
	} else {
	    if (chkcols != ncols) {
		pprintf(prn, M_("   ...but row %d has %d fields: aborting\n"),
			csvinfo->n, chkcols);
		pputs(prn, msg);
		goto csv_bailout;
	    }
	}
    }

    /* need to decrement csvinfo->n to allow for var headings */
    csvinfo->n -= 1;

    csvinfo->v = (blank_1 || obs_1)? ncols: ncols + 1;
    pprintf(prn, M_("   number of variables: %d\n"), csvinfo->v - 1);
    pprintf(prn, M_("   number of non-blank lines: %d\n"), 
	    csvinfo->n + 1);

    /* end initial checking */
    fclose(fp);
    fp = NULL;

    if (csvinfo->n == 0) {
	pputs(prn, M_("Invalid data file\n"));
	goto csv_bailout;
    }

    /* initialize datainfo and Z */
    if (start_new_Z(&csvZ, csvinfo, 0)) {
	pputs(prn, M_("Out of memory\n"));
	goto csv_bailout;
    }

    if (blank_1 || obs_1) {
	csvinfo->markers = REGULAR_MARKERS;
	if (dataset_allocate_markers(csvinfo)) {
	    pputs(prn, M_("Out of memory\n"));
	    goto csv_bailout;
	}
    }

    /* second pass */
    fp = fopen(fname, "r");
    if (fp == NULL) {
	goto csv_bailout;
    }

    /* parse the line containing variable names */
    pputs(prn, M_("scanning for variable names...\n"));

    while (fgets(line, maxlen + 1, fp)) {
	if (*line == '#') continue;
	if (string_is_blank(line)) continue;
	else break;
    }
    compress_csv_line(line, delim);    

    p = line;
    if (delim == ' ' && *p == ' ') p++;
    pprintf(prn, M_("   line: %s\n"), p);
    
    numcount = 0;
    for (k=0; k<ncols; k++) {
	int nv = 0;

	i = 0;
	while (*p && *p != delim) {
	    if (i < 31) csvstr[i++] = *p;
	    p++;
	}
	if (*p == delim) p++;
	csvstr[i] = 0;
	if (k == 0 && (blank_1 || obs_1)) {
	    ;
	} else {
	    nv = (blank_1 || obs_1)? k : k + 1;

	    if (*csvstr == '\0') {
		pprintf(prn, M_("   variable name %d is missing: aborting\n"), nv);
		pputs(prn, msg);
		goto csv_bailout;
	    } else {
		csvinfo->varname[nv][0] = 0;
		strncat(csvinfo->varname[nv], csvstr, 8);
		if (isdigit(*csvstr)) {
		    numcount++;
		} else {
		    if (check_varname(csvinfo->varname[nv])) {
			pprintf(prn, "%s\n", gretl_errmsg);
			*gretl_errmsg = '\0';
			goto csv_bailout;
		    }
		}
	    }
	}
	if (nv == csvinfo->v - 1) break;
    }

    if (numcount == csvinfo->v - 1) {
	pputs(prn, M_("it seems there are no variable names\n"));
	if (dataset_add_obs(&csvZ, csvinfo)) {
	    pputs(prn, _("Out of memory\n"));
	    goto csv_bailout;
	}
	auto_name_vars = 1;
	rewind(fp);
    } else if (numcount > 0) {
	for (i=0; i<csvinfo->v; i++) {
	    if (check_varname(csvinfo->varname[i])) {
		pprintf(prn, "%s\n", gretl_errmsg);
		*gretl_errmsg = '\0';
		break;
	    }
	}	    
	goto csv_bailout;
    }
    
#ifdef ENABLE_NLS
    if ((*ppdinfo)->decpoint != ',') setlocale(LC_NUMERIC, "C");
#endif

    pputs(prn, M_("scanning for row labels and data...\n"));

    t = 0;
    while (fgets(line, maxlen + 1, fp)) {
	int nv;

	if (*line == '#') continue;
	if (string_is_blank(line)) continue;
	compress_csv_line(line, delim);
	p = line;
	if (delim == ' ' && *p == ' ') p++;

	for (k=0; k<ncols; k++) {
	    i = 0;
	    while (*p && *p != delim) {
		if (i < 31) csvstr[i++] = *p;
		else {
		    pprintf(prn, M_("warning: truncating data at row %d, column %d\n"),
			    t+1, k+1);
		}
		p++;
	    }
	    if (*p == delim) p++;
	    csvstr[i] = 0;
	    if (k == 0 && (blank_1 || obs_1)) {
		csvinfo->S[t][0] = 0;
		strncat(csvinfo->S[t], csvstr, 8);
	    } else {
		nv = (blank_1 || obs_1)? k : k + 1;
		if (csv_missval(csvstr, k, t+1, prn)) {
		    csvZ[nv][t] = NADBL;
		} else {
		    if (check_atof(csvstr)) {
			pprintf(prn, "%s\n", gretl_errmsg);
			*gretl_errmsg = '\0';
			goto csv_bailout;
		    } else {
			csvZ[nv][t] = atof(csvstr);
		    }
		}
	    }
	}
	t++;
	if (t == csvinfo->n) break;
    }

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif

    csvinfo->t1 = 0;
    csvinfo->t2 = csvinfo->n - 1;
    if (blank_1 || obs_1) {
	markertest = test_label(csvinfo, &missvec, prn);
    }
    if (markertest > 0) {
	pputs(prn, M_("taking date information from row labels\n\n"));
    } else {
	pputs(prn, M_("treating these as undated data\n\n"));
	dataset_dates_defaults(csvinfo);
    }

#ifdef PAD_DAILY_DATA
    if (missvec != NULL) {
	daily_data_resize(&csvZ, csvinfo, missvec);
	free(missvec);
    }
#endif

    if (csvinfo->pd != 1 || strcmp(csvinfo->stobs, "1")) 
        csvinfo->time_series = TIME_SERIES;

    /* If there were observation labels and they were not interpretable
       as dates, and they weren't simply "1, 2, 3, ...", then they 
       should probably be preserved; otherwise discard them. */

    if (csvinfo->S != NULL && markertest >= 0 && 
	csvinfo->markers != DAILY_DATE_STRINGS) {
	csvinfo->markers = NO_MARKERS;
	for (i=0; i<csvinfo->n; i++) {
	    free(csvinfo->S[i]);
	}
	free(csvinfo->S);
	csvinfo->S = NULL;
    }

    if (auto_name_vars) {
	for (i=1; i<csvinfo->v; i++) {
	    sprintf(csvinfo->varname[i], "v%d", i);
	}
    }

    if (*pZ == NULL) {
	*pZ = csvZ;
	if (*ppdinfo != NULL) free(*ppdinfo);
	*ppdinfo = csvinfo;
    } else {
	if (merge_data(pZ, *ppdinfo, csvZ, csvinfo, prn))
	    goto csv_bailout;
    }

    fclose(fp); 
    free(line);

    console_off();

    return 0;

 csv_bailout:
    if (fp != NULL) fclose(fp);
    if (line != NULL) free(line);
    if (csvinfo != NULL) 
	clear_datainfo(csvinfo, CLEAR_FULL);

    console_off();

    return 1;
}

/**
 * add_case_markers:
 * @pdinfo: data information struct.
 * @fname: name of file containing case markers.
 * 
 * Read case markers (strings of 8 characters or less that identify
 * the observations) from a file, and associate tham with the 
 * current data set.
 * 
 * Returns: 0 on successful completion, non-zero otherwise.
 *
 */

int add_case_markers (DATAINFO *pdinfo, const char *fname)
{
    FILE *fp;
    char **S, marker[OBSLEN];
    int i, t;

    fp = fopen(fname, "r");
    if (fp == NULL) return E_FOPEN;
    
    S = malloc(pdinfo->n * sizeof *S);
    if (S == NULL) return E_ALLOC; 
    for (i=0; i<pdinfo->n; i++) {
	S[i] = malloc(OBSLEN);
	if (S[i] == NULL) return E_ALLOC; 
    }

    for (t=0; t<pdinfo->n; t++) {
	eatspace(fp);
	if (fscanf(fp, "%8s", marker) != 1) {
	    for (i=0; i<pdinfo->n; i++) free (S[i]);
	    free(S);
	    fclose(fp);
	    return 1;
	}
	marker[8] = '\0';
	strcpy(S[t], marker);
    }
    fclose(fp);

    if (pdinfo->S != NULL) {
	for (i=0; i<pdinfo->n; i++) free (pdinfo->S[i]);
	free(pdinfo->S);
    }

    pdinfo->S = S;
    pdinfo->markers = REGULAR_MARKERS;

    return 0;
}

/* ................................................. */

static char *unspace (char *s)
{
    int i;
    size_t n = strlen(s);

    for (i=n-1; i>=0; i--) { 
	if (s[i] == ' ') s[i] = '\0';
	else break;
    }
    return s;
}

/* #define BOX_DEBUG 1 */

/**
 * import_box:
 * @pZ: pointer to data set.
 * @pdinfo: data information struct.
 * @fname: name of CSV file.
 * @prn: gretl printing struct.
 * 
 * Open a BOX1 data file (as produced by the US Census Bureau's
 * Data Extraction Service) and read the data into
 * the current work space.
 * 
 * Returns: 0 on successful completion, non-zero otherwise.
 *
 */

int import_box (double ***pZ, DATAINFO **ppdinfo, 
		const char *fname, PRN *prn)
{
    int c, cc, i, t, v, realv, gotdata;
    int maxline, dumpvars;
    char tmp[48];
    unsigned *varsize = NULL, *varstart = NULL;
    char *line = NULL;
    double x;
    FILE *fp;
    DATAINFO *boxinfo;
    double **boxZ = NULL;
    int err = 0;

    check_for_console(prn);

    fp = fopen(fname, "r");
    if (fp == NULL) {
	pprintf(prn, M_("Couldn't open %s\n"), fname);
	err = E_FOPEN;
	goto box_bailout;
    }

    boxinfo = datainfo_new();
    if (boxinfo == NULL) {
	pputs(prn, M_("Out of memory\n"));
	err = E_ALLOC;
	goto box_bailout;
    }

    pprintf(prn, "%s %s...\n", M_("parsing"), fname);

    /* first pass: find max line length, number of vars and number
       of observations, plus basic sanity check */
    cc = maxline = 0;
    boxinfo->v = 1;
    do {
	c = getc(fp); 
	if (c != EOF && c != 10 && !isprint((unsigned char) c)) {
	    pprintf(prn, M_("Binary data (%d) encountered: this is not a valid "
		   "BOX1 file\n"), c);
	    fclose(fp);
	    err = 1;
	    goto box_bailout;
	}
	if (c == '\n') {
	    if (cc > maxline) maxline = cc;
	    cc = 0;
	    if ((c = getc(fp)) != EOF) {
		tmp[0] = c; cc++;
	    } else break;
	    if ((c = getc(fp)) != EOF) {
		tmp[1] = c; cc++;
	    } else break;
	    tmp[2] = '\0';
	    if (!strcmp(tmp, "03")) boxinfo->v += 1;
	    else if (!strcmp(tmp, "99")) boxinfo->n += 1;
	} else
	    cc++;
    } while (c != EOF);

    fclose(fp);

    pprintf(prn, M_("   found %d variables\n"), boxinfo->v - 1);
    pprintf(prn, M_("   found %d observations\n"), boxinfo->n);
    pprintf(prn, M_("   longest line = %d characters\n"), maxline); 
    maxline += 2;

    /* allocate space for data etc */
    pputs(prn, M_("allocating memory for data... "));

    if (start_new_Z(&boxZ, boxinfo, 0)) {
	err = E_ALLOC;
	goto box_bailout;
    }

    varstart = malloc((boxinfo->v - 1) * sizeof *varstart);
    varsize = malloc((boxinfo->v - 1) * sizeof *varsize);
    line = malloc(maxline);

    if (varstart == NULL || varsize == NULL || line == NULL) {
	free(varstart);
	free(varsize);
	free(line);
	err = E_ALLOC;
	goto box_bailout;
    }

    pputs(prn, M_("done\n"));

    fp = fopen(fname, "r");
    if (fp == NULL) {
	err = E_FOPEN;
	goto box_bailout;
    }
    pputs(prn, M_("reading variable information...\n"));

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "C");
#endif

    /* second pass: get detailed info on variables */
    v = 0; realv = 1; t = 0;
    dumpvars = 0; gotdata = 0;
    while (fgets(line, maxline, fp)) {
	strncpy(tmp, line, 2);
	tmp[2] = '\0';
	switch (atoi(tmp)) {
	case 0: /* comment */
	    break;
	case 1: /* BOX info (ignored for now) */
	    break;
	case 2: /* raw data records types (ignored for now) */
	    break;
	case 3: /* variable info */
	    strncpy(boxinfo->varname[realv], line+11, 8);
	    boxinfo->varname[realv][8] = '\0';
	    unspace(boxinfo->varname[realv]);
	    lower(boxinfo->varname[realv]);
	    pprintf(prn, M_(" variable %d: '%s'\n"), v+1, boxinfo->varname[realv]);
#ifdef notdef  
	    /* This is wrong!  How do you identify character data? */
	    if (line[51] != '2') {
		pputs(prn, M_("   Non-numeric data: will be skipped\n"));
		varstart[v] = 0;
		varsize[v] = 0;
		v++;
		break;
	    }
#endif
	    strncpy(tmp, line+52, 6);
	    tmp[6] = '\0';
	    varstart[v] = atoi(tmp) - 1;
	    pprintf(prn, M_("   starting col. %d, "), varstart[v]);
	    strncpy(tmp, line+58, 4);
	    tmp[4] = '\0';
	    varsize[v] = atoi(tmp);
	    pprintf(prn, M_("field width %d, "), varsize[v]);
	    strncpy(tmp, line+62, 2);
	    tmp[2] = '\0';
	    pprintf(prn, M_("decimal places %d\n"), atoi(tmp));
	    tmp[0] = '\0';
	    strncpy(tmp, line+64, 20);
	    tmp[20] = '\0';
	    unspace(tmp);
	    if (strlen(tmp))
		pprintf(prn, M_("   Warning: coded variable (format '%s' "
			"in BOX file)\n"), tmp);
	    *VARLABEL(boxinfo, realv) = 0;
	    strncat(VARLABEL(boxinfo, realv), line+87, 99);
	    unspace(VARLABEL(boxinfo, realv));
	    pprintf(prn, M_("   definition: '%s'\n"), VARLABEL(boxinfo, realv));
	    realv++;
	    v++;
	    break;
	case 4: /* category info (ignored for now) */
	    break;
	case 99: /* data line */
	    realv = 1;
 	    for (i=0; i<v; i++) {
		if (varstart[i] == 0 && varsize[i] == 0) {
		    if (!gotdata) dumpvars++;
		    continue;
		}
		strncpy(tmp, line + varstart[i], varsize[i]);
		tmp[varsize[i]] = '\0';
		top_n_tail(tmp);

		if (check_atof(tmp)) {
		    pprintf(prn, "%s\n", gretl_errmsg);
		    x = NADBL;
		} else {
		    x = atof(tmp);
		}
#ifdef BOX_DEBUG
		fprintf(stderr, "read %d chars from pos %d: '%s' -> %g\n",
			varsize[i], varstart[i], tmp, x); 
#endif
		boxZ[realv][t] = x;
#ifdef BOX_DEBUG
		fprintf(stderr, "setting Z[%d][%d] = %g\n", realv, t, x);
#endif
		realv++;
	    }
	    t++;
	    gotdata = 1;
	    break;
	default:
	    break;
	}
    }

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif

    pputs(prn, M_("done reading data\n"));
    fclose(fp);

    free(varstart);
    free(varsize);
    free(line);

    dataset_dates_defaults(boxinfo);

    if (dumpvars) {
	dataset_drop_vars(dumpvars, &boxZ, boxinfo);
	pprintf(prn, M_("Warning: discarded %d non-numeric variable(s)\n"), 
		dumpvars);
    }

    if (*pZ == NULL) {
	*pZ = boxZ;
	if (*ppdinfo != NULL) free(*ppdinfo);
	*ppdinfo = boxinfo;
    }

 box_bailout:

    console_off();

    return err;
}

static int xmlfile (const char *fname)
{
    gzFile fz;
    char test[6];
    int ret = 0;

    fz = gzopen(fname, "rb");
    if (fz != Z_NULL) {
	if (gzread(fz, test, 5)) {
	    test[5] = '\0';
	    if (!strcmp(test, "<?xml")) ret = 1;
	} 
	gzclose(fz);
    } 
    return ret;
} 

static int file_has_suffix (const char *fname, const char *sfx)
{
    const char *p = strrchr(fname, '.');

    if (p != NULL && !strcmp(p + 1, sfx)) return 1;
    else return 0;
}

/**
 * detect_filetype:
 * @fname: name of file to examine.
 * @ppaths: path information struct.
 * @prn: gretl printing struct.
 * 
 * Attempt to determine the type of a file to be opened in gretl:
 * data file (native, CSV or BOX), or command script.
 * 
 * Returns: integer code indicating the type of file (see #gretl_filetypes).
 *
 */

int detect_filetype (char *fname, PATHS *ppaths, PRN *prn)
{
    int i, c, ftype = GRETL_NATIVE_DATA;
    char teststr[5];
    FILE *fp;

    /* might be a script file? (watch out for DOS-mangled names) */
    if (file_has_suffix(fname, "inp") ||
	file_has_suffix(fname, "gre") ||
	file_has_suffix(fname, "GRE"))
	return GRETL_SCRIPT;
    if (file_has_suffix(fname, "gretl"))
	return GRETL_SCRIPT; 
    if (file_has_suffix(fname, "gnumeric"))
	return GRETL_GNUMERIC;
    if (file_has_suffix(fname, "xls"))
	return GRETL_EXCEL;
    if (file_has_suffix(fname, "des"))
	return GRETL_DES_DATA;
    if (file_has_suffix(fname, "bin"))
	return GRETL_NATIVE_DB;
    if (file_has_suffix(fname, "rat"))
	return GRETL_RATS_DB;

    addpath(fname, ppaths, 0); 

    if (xmlfile(fname))
	return GRETL_XML_DATA;   

    fp = fopen(fname, "r");
    if (fp == NULL) { 
	return GRETL_NATIVE_DATA; /* may be native file in different location */
    }

    /* look at extension */
    if (file_has_suffix(fname, "csv")) 
	ftype = GRETL_CSV_DATA;
    else if (file_has_suffix(fname, "txt"))
	ftype = GRETL_CSV_DATA;
    else if (file_has_suffix(fname, "box"))
	ftype = GRETL_BOX_DATA;

    /* then take a peek at content */
    for (i=0; i<80; i++) {
	c = getc(fp);
	if (c == EOF || c == '\n') break;
	if (!isprint(c) && c != '\r' && c != '\t') {
	    ftype = GRETL_NATIVE_DATA; /* native binary data? */
	    break;
	}
	if (i < 4) teststr[i] = c;
    }

    fclose(fp);
    teststr[4] = 0;

    switch (ftype) {
    case GRETL_NATIVE_DATA: 
	return GRETL_NATIVE_DATA;
    case GRETL_CSV_DATA: 
	return GRETL_CSV_DATA;
    case GRETL_BOX_DATA: 
	if (strcmp(teststr, "00**") == 0) return GRETL_BOX_DATA;
	else {
	    pputs(prn, M_("box file seems to be malformed\n"));
	    return GRETL_UNRECOGNIZED;
	}
    }
    return GRETL_NATIVE_DATA; /* FIXME? */
}

#define UTF const xmlChar *

/* #define XML_DEBUG */

static char *simple_fname (char *dest, const char *src)
{
    char *p;
    const char *s;

    s = strrchr(src, SLASH);

    /* take last part of src filename */
    if (s != NULL) {
        strcpy(dest, s + 1);
    } else {
        strcpy(dest, src);
    }

    /* trash any extension */
    p = strrchr(dest, '.');
    if (p != NULL && strlen(dest) > 3) {
	*p = '\0';
    }

    return dest;
}

/**
 * write_xmldata:
 * @fname: name of file to write.
 * @list: list of variables to write.
 * @Z: data matrix.
 * @pdinfo: data information struct.
 * @opt: if non-zero, write gzipped data, else plain.
 * @ppaths: pointer to paths information (or NULL).
 * 
 * Write out in xml a data file containing the values of the given set
 * of variables.
 * 
 * Returns: 0 on successful completion, non-zero on error.
 * 
 */

static int write_xmldata (const char *fname, const int *list, 
			  double **Z, const DATAINFO *pdinfo, 
			  int opt, PATHS *ppaths)
{
    int err, i, t;
    FILE *fp = NULL;
    gzFile *fz = Z_NULL;
    int *pmax = NULL, tsamp = pdinfo->t2 - pdinfo->t1 + 1;
    char startdate[OBSLEN], enddate[OBSLEN], datname[MAXLEN], type[32];
    char *xmlbuf = NULL;
    long sz = 0L;
    void *handle = NULL;
    int (*show_progress) (long, long, int) = NULL;
#ifdef USE_GTK2
    const char *enc = "UTF-8";
#else
    const char *enc = "ISO-8859-1";
#endif

    err = 0;

    if (opt) {
	fz = gzopen(fname, "wb");
	if (fz == Z_NULL) err = 1;
    } else {
	fp = fopen(fname, "wb");
	if (fp == NULL) err = 1;
    }
    if (err) {
	sprintf(gretl_errmsg, _("Couldn't open %s for writing"), fname);
	return 1;
    }

    pmax = malloc(list[0] * sizeof *pmax);
    if (pmax == NULL) {
	sprintf(gretl_errmsg, _("Out of memory"));
	return 1;
    } 

    sz = (tsamp * pdinfo->v * sizeof(double));
    if (sz > 100000) {
	fprintf(stderr, I_("Writing %ld Kbytes of data\n"), sz / 1024);
	if (ppaths == NULL) sz = 0L;
    } else {
	sz = 0L;
    }

    if (sz) {
	show_progress = get_plugin_function("show_progress", &handle);
	if (show_progress == NULL) sz = 0L;
    }

    if (sz) (*show_progress)(0, sz, SP_SAVE_INIT); 

    for (i=1; i<=list[0]; i++) {
	if (pdinfo->vector[list[i]]) {
	    pmax[i-1] = get_precision(&Z[list[i]][pdinfo->t1], tsamp, 10);
	} else {
	    pmax[i-1] = SCALAR_DIGITS;
	}
    }

    ntodate(startdate, pdinfo->t1, pdinfo);
    ntodate(enddate, pdinfo->t2, pdinfo);

    simple_fname(datname, fname);
    xmlbuf = gretl_xml_encode(datname);
    if (xmlbuf == NULL) return 1;

    if (opt) {
	gzprintf(fz, "<?xml version=\"1.0\" encoding=\"%s\"?>\n"
		 "<!DOCTYPE gretldata SYSTEM \"gretldata.dtd\">\n\n"
		 "<gretldata name=\"%s\" frequency=\"%d\" "
		 "startobs=\"%s\" endobs=\"%s\" ", 
		 enc, datname, pdinfo->pd, startdate, enddate);
    } else {
	fprintf(fp, "<?xml version=\"1.0\" encoding=\"%s\"?>\n"
		"<!DOCTYPE gretldata SYSTEM \"gretldata.dtd\">\n\n"
		"<gretldata name=\"%s\" frequency=\"%d\" "
		"startobs=\"%s\" endobs=\"%s\" ", 
		enc, datname, pdinfo->pd, startdate, enddate);
    }

    free(xmlbuf);

    switch (pdinfo->time_series) {
    case 0:
	strcpy(type, "cross-section"); break;
    case TIME_SERIES:
	strcpy(type, "time-series"); break;
    case STACKED_TIME_SERIES:
	strcpy(type, "stacked-time-series"); break;
    case STACKED_CROSS_SECTION:
	strcpy(type, "stacked-cross-section"); break;
    default:
	strcpy(type, "cross-section"); break;
    }

    if (opt) {
	gzprintf(fz, "type=\"%s\">\n", type);
    } else {
	fprintf(fp, "type=\"%s\">\n", type);
    }

    /* first deal with description, if any */
    if (pdinfo->descrip != NULL) {
	xmlbuf = gretl_xml_encode(pdinfo->descrip);
	if (xmlbuf == NULL) return 1;
	else {
	    if (opt) {
		gzputs(fz, "<description>\n");
		gzputs(fz, xmlbuf);
		gzputs(fz, "</description>\n");
	    }
	    else {
		fprintf(fp, "<description>\n%s</description>\n", xmlbuf);
	    }
	    free(xmlbuf);
#ifdef XML_DEBUG
	    fprintf(stderr, "xmlbuf encoded buffer freed\n");
#endif
	}
    }

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "C");
#endif

    /* then listing of variable names and labels */
    if (opt) {
	gzprintf(fz, "<variables count=\"%d\">\n", list[0]);
    } else {
	fprintf(fp, "<variables count=\"%d\">\n", list[0]);
    }

    for (i=1; i<=list[0]; i++) {
	xmlbuf = gretl_xml_encode(pdinfo->varname[list[i]]);
	if (xmlbuf == NULL) return 1;
	else {
	    if (opt) gzprintf(fz, "<variable name=\"%s\"", xmlbuf);
	    else fprintf(fp, "<variable name=\"%s\"", xmlbuf);
	    free(xmlbuf);
	}
	if (!pdinfo->vector[list[i]]) {
	    if (opt) { 
		if (pmax[i-1] == PMAX_NOT_AVAILABLE) {
		    gzprintf(fz, "\n role=\"scalar\" value=\"%.12g\"",
			     Z[list[i]][0]);
		} else {
		    gzprintf(fz, "\n role=\"scalar\" value=\"%.*g\"",
			     pmax[i-1], Z[list[i]][0]);
		}
	    } else {
		if (pmax[i-1] == PMAX_NOT_AVAILABLE) {
		    fprintf(fp, "\n role=\"scalar\" value=\"%.12g\"",
			    Z[list[i]][0]);
		} else {
		    fprintf(fp, "\n role=\"scalar\" value=\"%.*g\"",
			    pmax[i-1], Z[list[i]][0]);
		}
	    }
	}
	if (*VARLABEL(pdinfo, list[i])) {
	    xmlbuf = gretl_xml_encode(VARLABEL(pdinfo, list[i]));
	    if (xmlbuf == NULL) return 1;
	    else {
		if (opt) gzprintf(fz, "\n label=\"%s\"", xmlbuf);
		else fprintf(fp, "\n label=\"%s\"", xmlbuf);
		free(xmlbuf);
	    }
	} 
	if (*DISPLAYNAME(pdinfo, list[i])) {
	    xmlbuf = gretl_xml_encode(DISPLAYNAME(pdinfo, list[i]));
	    if (xmlbuf == NULL) return 1;
	    else {
		if (opt) gzprintf(fz, "\n displayname=\"%s\"", xmlbuf);
		else fprintf(fp, "\n displayname=\"%s\"", xmlbuf);
		free(xmlbuf);
	    }
	} 
	if (COMPACT_METHOD(pdinfo, list[i]) != COMPACT_NONE) {
	    xmlbuf = compact_method_to_string(COMPACT_METHOD(pdinfo, list[i]));
	    if (opt) gzprintf(fz, "\n compact-method=\"%s\"", xmlbuf);
	    else fprintf(fp, "\n compact-method=\"%s\"", xmlbuf);
	} 
	if (opt) gzputs(fz, "\n/>\n");
	else fputs("\n/>\n", fp);

    }
    if (opt) gzputs(fz, "</variables>\n");
    else fputs("</variables>\n", fp);

    /* then listing of observations */
    if (opt) {
	gzprintf(fz, "<observations count=\"%d\" labels=\"%s\">\n",
		tsamp, (pdinfo->markers && pdinfo->S != NULL)? "true" : "false");
    } else {
	fprintf(fp, "<observations count=\"%d\" labels=\"%s\">\n",
		tsamp, (pdinfo->markers && pdinfo->S != NULL)? "true" : "false");
    }

    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	if (opt) {
	    gzputs(fz, "<obs");
	} else {
	    fputs("<obs", fp);
	}
	if (pdinfo->markers && pdinfo->S != NULL) {
	    if (opt) gzprintf(fz, " label=\"%s\">", pdinfo->S[t]);
	    else fprintf(fp, " label=\"%s\">", pdinfo->S[t]);
	} else {
	    if (opt) gzputs(fz, ">");
	    else fputc('>', fp);
	}
	for (i=1; i<=list[0]; i++) {
	    if (!pdinfo->vector[list[i]]) continue;
	    if (na(Z[list[i]][t])) {
		if (opt) gzputs(fz, "NA ");
		else fputs("NA ", fp);
	    } else {
		if (opt) {
		    if (pmax[i-1] == PMAX_NOT_AVAILABLE) {
			gzprintf(fz, "%.12g ", Z[list[i]][t]);
		    } else {
			gzprintf(fz, "%.*f ", pmax[i-1], Z[list[i]][t]);
		    }
		}
		else {
		    if (pmax[i-1] == PMAX_NOT_AVAILABLE) {
			fprintf(fp, "%.12g ", Z[list[i]][t]);
		    } else {
			fprintf(fp, "%.*f ", pmax[i-1], Z[list[i]][t]);
		    }
		}
	    }
	}
	if (opt) gzputs(fz, "</obs>\n");
	else fputs("</obs>\n", fp);
	if (sz && t && ((t - pdinfo->t1) % 50 == 0)) { 
	    (*show_progress) (50, tsamp, SP_NONE);
	}
    }

    if (opt) {
	gzprintf(fz, "</observations>\n</gretldata>\n");
    } else {
	fprintf(fp, "</observations>\n</gretldata>\n");
    }

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif

    if (sz) {
	(*show_progress)(0, pdinfo->t2 - pdinfo->t1 + 1, SP_FINISH);
	close_plugin(handle);
    }  

    if (pmax) free(pmax);
    if (fp != NULL) fclose(fp);
    if (fz != Z_NULL) gzclose(fz);

    return 0;
}

static int process_varlist (xmlNodePtr node, DATAINFO *pdinfo, double ***pZ)
{
    xmlNodePtr cur;
    int i;
    char *tmp = xmlGetProp(node, (UTF) "count");

    if (tmp) {
	int v, err = 0;

	if (sscanf(tmp, "%d", &v) == 1) {
	    pdinfo->v = v + 1;
	} else {
	    sprintf(gretl_errmsg, _("Failed to parse count of variables"));
	    err = 1;
	}
	if (!err && dataset_allocate_varnames(pdinfo)) {
	    sprintf(gretl_errmsg, _("Out of memory reading data file"));
	    err = 1;
	}
	if (!err) {
	    *pZ = malloc(pdinfo->v * sizeof **pZ);
	    if (*pZ == NULL) {
		sprintf(gretl_errmsg, _("Out of memory reading data file"));
		err = 1;
	    }
	}		
	free(tmp);
	if (err) return 1;
    }
    else {
	sprintf(gretl_errmsg, _("Got no variables"));
	return 1;
    }

    /* now get individual variable info: names and labels */
    cur = node->xmlChildrenNode;
    while (cur && xmlIsBlankNode(cur)) {
	cur = cur->next;
    }
    if (cur == 0) {
	sprintf(gretl_errmsg, _("Got no variables"));
	return 1;
    }

    i = 1;
    while (cur != NULL) {
        if (!xmlStrcmp(cur->name, (UTF) "variable")) {
	    tmp = xmlGetProp(cur, (UTF) "name");
	    if (tmp) {
		pdinfo->varname[i][0] = 0;
		strncat(pdinfo->varname[i], tmp, 8);
		free(tmp);
	    } else {
		sprintf(gretl_errmsg, _("Variable %d has no name"), i);
		return 1;
	    }
	    tmp = xmlGetProp(cur, (UTF) "label");
	    if (tmp) {
		*VARLABEL(pdinfo, i) = 0;
		strncat(VARLABEL(pdinfo, i), tmp, MAXLABEL-1);
		free(tmp);
	    }
	    tmp = xmlGetProp(cur, (UTF) "displayname");
	    if (tmp) {
		*DISPLAYNAME(pdinfo, i) = 0;
		strncat(DISPLAYNAME(pdinfo, i), tmp, MAXDISP-1);
		free(tmp);
	    }
	    tmp = xmlGetProp(cur, (UTF) "compact-method");
	    if (tmp) {
		COMPACT_METHOD(pdinfo, i) = compact_string_to_int(tmp);
		free(tmp);
	    }
	    tmp = xmlGetProp(cur, (UTF) "role");
	    if (tmp) {
		if (!strcmp(tmp, "scalar")) {
		    char *val = xmlGetProp(cur, (UTF) "value");
		    
		    if (val) {
			double xx = atof(val);

			free(val);
			(*pZ)[i] = malloc(sizeof ***pZ);
			(*pZ)[i][0] = xx;
			pdinfo->vector[i] = 0;
		    }
		}
		free(tmp);
	    }
	    i++;
	}	    
	cur = cur->next;
    }
   
    if (i != pdinfo->v) {
	sprintf(gretl_errmsg, _("Number of variables does not match declaration"));
	return 1;
    }
    else return 0;
}

static int process_values (double **Z, DATAINFO *pdinfo, int t, char *s)
{
    int i;
    double x;

    for (i=1; i<pdinfo->v; i++) {
	if (!pdinfo->vector[i]) continue;
	s = strpbrk(s, "01234567890+-NA");
	if (!strncmp(s, "NA", 2)) x = NADBL;
	else if (*s && (sscanf(s, "%lf", &x) != 1)) {
	    sprintf(gretl_errmsg, _("Failed to parse data values at obs %d"), t+1);
	    return 1;
	}
	Z[i][t] = x;
	s = strpbrk(s, " \t\n\r");
    }

    return 0;
}

static int process_observations (xmlDocPtr doc, xmlNodePtr node, 
				 double ***pZ, DATAINFO *pdinfo,
				 long progress)
{
    xmlNodePtr cur;
    char *tmp = xmlGetProp(node, (UTF) "count");
    int i, t;
    void *handle;
    int (*show_progress) (long, long, int) = NULL;

    if (progress > 0) {
	show_progress = get_plugin_function("show_progress", &handle);
	if (show_progress == NULL) progress = 0L;
    }

    if (tmp) {
	int n;

	if (sscanf(tmp, "%d", &n) == 1) 
	    pdinfo->n = n;
	else {
	    sprintf(gretl_errmsg, _("Failed to parse number of observations"));
	    return 1;
	}
	free(tmp);
    } else 
	return 1;

    tmp = xmlGetProp(node, (UTF) "labels");
    if (tmp) {
	if (!strcmp(tmp, "true")) {
	    pdinfo->markers = REGULAR_MARKERS;
	    if (dataset_allocate_markers(pdinfo)) {
		sprintf(gretl_errmsg, "Out of memory");
		return 1;
	    }
	} else if (strcmp(tmp, "false")) {
	    sprintf(gretl_errmsg, _("labels attribute for observations must be "
		    "'true' or 'false'"));
	    return 1;
	}
	free(tmp);
    } else
	return 1;

    if (pdinfo->endobs[0] == '\0') 
	sprintf(pdinfo->endobs, "%d", pdinfo->n);

    pdinfo->t2 = pdinfo->n - 1;

    for (i=0; i<pdinfo->v; i++) {
	if (!pdinfo->vector[i]) continue;
	(*pZ)[i] = malloc(pdinfo->n * sizeof ***pZ);
	if ((*pZ)[i] == NULL) return 1;
    }

    for (t=0; t<pdinfo->n; t++) (*pZ)[0][t] = 1.0;

    /* now get individual obs info: labels and values */
    cur = node->xmlChildrenNode;
    while (cur && xmlIsBlankNode(cur)) {
	cur = cur->next;
    }
    if (cur == 0) {
	sprintf(gretl_errmsg, _("Got no observations\n"));
	return 1;
    }

    if (progress) (*show_progress)(0L, progress, SP_LOAD_INIT);

    t = 0;
    while (cur != NULL) {
        if (!xmlStrcmp(cur->name, (UTF) "obs")) {
	    if (pdinfo->markers) {
		tmp = xmlGetProp(cur, (UTF) "label");
		if (tmp) {
		    strncpy(pdinfo->S[t], tmp, 8);
		    pdinfo->S[t][8] = '\0';
		    free(tmp);
		} else {
		    sprintf(gretl_errmsg, _("Case marker missing at obs %d"), t+1);
		    return 1;
		}
	    }
	    tmp = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
	    if (tmp) {
		if (process_values(*pZ, pdinfo, t, tmp))
		    return 1;
		free(tmp);
		t++;
	    } else {
		sprintf(gretl_errmsg, _("Values missing at observation %d"), t+1);
		return 1;
	    }
	}	    
	cur = cur->next;
	if (progress && t > 0 && t % 50 == 0) {
	    (*show_progress) (50L, (long) pdinfo->n, SP_NONE);
	}
    }

    if (progress) {
	(*show_progress)(0L, (long) pdinfo->n, SP_FINISH);
	close_plugin(handle);
    }

    if (t != pdinfo->n) {
	sprintf(gretl_errmsg, _("Number of observations does not match declaration"));
	return 1;
    }

    else return 0;
}

static long get_filesize (const char *fname)
{
    struct stat buf;

    if (stat(fname, &buf) == 0) {
        return buf.st_size;
    } else {
        return -1;
    }
}

static void colon_to_point (char *str)
{
    char *p = str;

    while (*p) {
	if (*p == ':') {
	    *p = '.';
	    break;
	}
	p++;
    }
}

/**
 * get_xmldata:
 * @pZ: pointer to data set.
 * @pdinfo: data information struct.
 * @fname: name of file to try.
 * @ppaths: path information struct.
 * @data_status: indicator for whether a data file is currently open
 * in gretl's work space (not-0) or not (0).
 * @prn: where messages should be written.
 * @gui: should = 1 if the function is launched from the GUI, else 0.
 * 
 * Read data from file into gretl's work space, allocating space as
 * required.
 * 
 * Returns: 0 on successful completion, non-zero otherwise.
 *
 */

int get_xmldata (double ***pZ, DATAINFO *pdinfo, char *fname,
		 PATHS *ppaths, int data_status, PRN *prn, int gui) 
{
    xmlDocPtr doc;
    xmlNodePtr cur;
    char *tmp;
    int gotvars = 0, gotobs = 0, err = 0;
    long fsz, progress = 0L;

    *gretl_errmsg = '\0';

    check_for_console(prn);

    /* COMPAT: Do not generate nodes for formatting spaces */
    LIBXML_TEST_VERSION
	xmlKeepBlanksDefault(0);

    fsz = get_filesize(fname);
    if (fsz > 100000) {
	fprintf(stderr, "%s %ld bytes %s...\n", 
		(is_gzipped(fname))? I_("Uncompressing") : I_("Reading"),
		fsz, I_("of data"));
	if (gui) progress = fsz;
    }

    doc = xmlParseFile(fname);
    if (doc == NULL) {
	sprintf(gretl_errmsg, _("xmlParseFile failed on %s"), fname);
	err = 1;
	goto xmldata_bailout;
    }

    /* clear any existing data info */
    if (data_status) clear_datainfo(pdinfo, CLEAR_FULL);

    cur = xmlDocGetRootElement(doc);
    if (cur == NULL) {
        sprintf(gretl_errmsg, _("%s: empty document"), fname);
	xmlFreeDoc(doc);
	err = 1;
	goto xmldata_bailout;
    }

    if (xmlStrcmp(cur->name, (UTF) "gretldata")) {
        sprintf(gretl_errmsg, _("File of the wrong type, root node not gretldata"));
	xmlFreeDoc(doc);
	err = 1;
	goto xmldata_bailout;
    }

    /* set some datainfo parameters */
    tmp = xmlGetProp(cur, (UTF) "type");
    if (tmp == NULL) {
	sprintf(gretl_errmsg, 
		_("Required attribute 'type' is missing from data file"));
	err = 1;
	goto xmldata_bailout;
    } else {
	if (!strcmp(tmp, "cross-section")) 
	    pdinfo->time_series = 0;
	else if (!strcmp(tmp, "time-series"))
	    pdinfo->time_series = TIME_SERIES;
	else if (!strcmp(tmp, "stacked-time-series"))
	    pdinfo->time_series = STACKED_TIME_SERIES;
	else if (!strcmp(tmp, "stacked-cross-section"))
	    pdinfo->time_series = STACKED_CROSS_SECTION;
	else {
	    sprintf(gretl_errmsg, _("Unrecognized type attribute for data file"));
	    free(tmp);
	    err = 1;
	    goto xmldata_bailout;
	}
	free(tmp);
    }

    pdinfo->pd = 1;
    tmp = xmlGetProp(cur, (UTF) "frequency");
    if (tmp) {
	int pd = 0;

	if (sscanf(tmp, "%d", &pd) == 1) {
	    pdinfo->pd = pd;
	} else {
	    strcpy(gretl_errmsg, _("Failed to parse data frequency"));
	    free(tmp);
	    err = 1;
	    goto xmldata_bailout;
	}
	free(tmp);
    }

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "C");
#endif

    strcpy(pdinfo->stobs, "1");
    tmp = xmlGetProp(cur, (UTF) "startobs");
    if (tmp) {
	char obstr[16];

	strcpy(obstr, tmp);
	colon_to_point(obstr);
	
	if (dataset_is_daily(pdinfo)) {
	    long ed = get_epoch_day(tmp);

	    if (ed < 0) err = 1;
	    else pdinfo->sd0 = ed;
	} else {
	    double x;

	    if (sscanf(obstr, "%lf", &x) != 1) err = 1;
	    else pdinfo->sd0 = x;
	}
	if (err) {
	    strcpy(gretl_errmsg, _("Failed to parse startobs"));
	    free(tmp);
	    err = 1;
	    goto xmldata_bailout;
	}
	strncpy(pdinfo->stobs, tmp, 8);
	pdinfo->stobs[8] = '\0';
	colonize_obs(pdinfo->stobs);
	free(tmp);
    }

    *pdinfo->endobs = '\0';
    tmp = xmlGetProp(cur, (UTF) "endobs");
    if (tmp) {
	if (dataset_is_daily(pdinfo)) {
	    long ed = get_epoch_day(tmp);

	    if (ed < 0) err = 1;
	} else {
	    double x;

	    if (sscanf(tmp, "%lf", &x) != 1) err = 1;
	} 
	if (err) {
	    strcpy(gretl_errmsg, _("Failed to parse endobs"));
	    free(tmp);
	    err = 1;
	    goto xmldata_bailout;
	}
	strncpy(pdinfo->endobs, tmp, 8);
	pdinfo->endobs[8] = '\0';
	colonize_obs(pdinfo->endobs);
	free(tmp);
    }

    /* Now walk the tree */
    cur = cur->xmlChildrenNode;
    while (cur != NULL && !err) {
        if (!xmlStrcmp(cur->name, (UTF) "description")) {
	    pdinfo->descrip = 
		xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
        } else if (!xmlStrcmp(cur->name, (UTF) "variables")) {
	    if (process_varlist(cur, pdinfo, pZ)) {
		err = 1;
	    } else {
		gotvars = 1;
	    }
	}
        else if (!xmlStrcmp(cur->name, (UTF) "observations")) {
	    if (!gotvars) {
		sprintf(gretl_errmsg, _("Variables information is missing"));
		err = 1;
	    }
	    if (process_observations(doc, cur, pZ, pdinfo, progress)) {
		err = 1;
	    } else {
		gotobs = 1;
	    }
	}
	if (!err) cur = cur->next;
    }

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif

    xmlFreeDoc(doc);
    xmlCleanupParser();

    if (err) goto xmldata_bailout;

    if (!gotvars) {
	sprintf(gretl_errmsg, _("Variables information is missing"));
	err = 1;
	goto xmldata_bailout;
    }
    if (!gotobs) {
	sprintf(gretl_errmsg, _("No observations were found"));
	err = 1;
	goto xmldata_bailout;
    }

    if (fname != ppaths->datfile) {
	strcpy(ppaths->datfile, fname);
    }

    pprintf(prn, M_("\nRead datafile %s\n"), fname);
    pprintf(prn, M_("periodicity: %d, maxobs: %d, "
		    "observations range: %s-%s\n\n"), pdinfo->pd, pdinfo->n,
	    pdinfo->stobs, pdinfo->endobs);

 xmldata_bailout:

    console_off();

    return err;
}

/**
 * get_xml_description:
 * @fname: name of file to try.
 * 
 * Read data description for gretl xml data file.
 * 
 * Returns: buffer containing description, or NULL on failure.
 *
 */

char *get_xml_description (const char *fname)
{
    xmlDocPtr doc;
    xmlNodePtr cur;
    char *buf = NULL;

    *gretl_errmsg = '\0';

    LIBXML_TEST_VERSION
	xmlKeepBlanksDefault(0);

    doc = xmlParseFile(fname);
    if (doc == NULL) {
	sprintf(gretl_errmsg, _("xmlParseFile failed on %s"), fname);
	return NULL;
    }

    cur = xmlDocGetRootElement(doc);
    if (cur == NULL) {
        sprintf(gretl_errmsg, _("%s: empty document"), fname);
	xmlFreeDoc(doc);
	return NULL;
    }

    if (xmlStrcmp(cur->name, (UTF) "gretldata")) {
        sprintf(gretl_errmsg, _("File of the wrong type, root node not gretldata"));
	xmlFreeDoc(doc);
	return NULL;
    }

    cur = cur->xmlChildrenNode;
    while (cur != NULL) {
        if (!xmlStrcmp(cur->name, (UTF) "description")) {
	    buf = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
	    break;
        }
	cur = cur->next;
    }

    xmlFreeDoc(doc);
    xmlCleanupParser();

    return buf;
}

int check_atof (const char *numstr)
{
    char *test;
    extern int errno;

    errno = 0;

    /* accept blank entries */
    if (*numstr == '\0') return 0;

    (void) strtod(numstr, &test);

    if (!strcmp(numstr, test)) {
	sprintf(gretl_errmsg, M_("'%s' -- no numeric conversion performed!"), numstr);
	return 1;
    }

    if (*test != '\0') {
	if (isprint(*test)) {
	    sprintf(gretl_errmsg, M_("Extraneous character '%c' in data"), *test);
	} else {
	    sprintf(gretl_errmsg, M_("Extraneous character (0x%x) in data"), *test);
	}
	return 1;
    }

    if (errno == ERANGE) {
	sprintf(gretl_errmsg, M_("'%s' -- number out of range!"), numstr);
	return 1;
    }

    return 0;
}
