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
#include <errno.h>

#define QUOTE '\''

static int prepZ (double **pZ, const DATAINFO *pdinfo);
static int writelbl (const char *lblfile, const int *list, 
		     const DATAINFO *pdinfo);
static int writehdr (const char *hdrfile, const int *list, 
		     const DATAINFO *pdinfo, const int opt);
static double obs_float (const DATAINFO *pdinfo, const int end);

static char STARTCOMMENT[3] = "(*";
static char ENDCOMMENT[3] = "*)";

/**
 * clear_datainfo:
 * @pdinfo: data information struct.
 * @subsample: 1 if @pdinfo relates to a sub-sample of a full data set,
 * 0 otherwise.
 *
 * Free the allocated content of a data information struct.
 * 
 */

void clear_datainfo (DATAINFO *pdinfo, int subsample)
{
    int i;

    if (pdinfo->S != NULL) {
	for (i=0; i<pdinfo->n; i++) 
	   free(pdinfo->S[i]); 
	free(pdinfo->S);
	pdinfo->S = NULL;
    } 
    /* if this is not a sub-sample datainfo, free varnames and labels */
    if (!subsample) {
	if (pdinfo->varname != NULL) {
	    for (i=0; i<pdinfo->v; i++) 
		free(pdinfo->varname[i]); 
	    free(pdinfo->varname);
	    pdinfo->varname = NULL;
	}
	if (pdinfo->label != NULL) {
	    for (i=0; i<pdinfo->v; i++) 
		free(pdinfo->label[i]); 
	    free(pdinfo->label);
	    pdinfo->label = NULL;
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
    char s[MAXLEN], *mybuf;
    int count = 0, bigger = 1, bufsize;

    if (fgets(s, MAXLEN-1, fp) == NULL) return 0;
    if (strncmp(s, STARTCOMMENT, 2) == 0) {
	*pbuf = malloc(20 * MAXLEN);
	if (*pbuf == NULL) return -1;
	**pbuf = '\0';
	strcat(*pbuf, s);
	do {
	    if (fgets(s, MAXLEN-1, fp) == NULL) break;
	    count++;
	    if (count > 20*bigger) {
		bigger++;
		bufsize = 20 * MAXLEN * bigger;
		mybuf = realloc(*pbuf, bufsize);
		if (mybuf == NULL) return -1;
		else *pbuf = mybuf;
	    }
	    strcat(*pbuf, s);
	} while (s != NULL && strncmp(s, ENDCOMMENT, 2));
    }
    return count;
}

/* ................................................. */

static int dataset_allocate_markers (DATAINFO *pdinfo)
{
    int i;

    pdinfo->S = malloc(pdinfo->n * sizeof(char *));
    if (pdinfo->S == NULL) return 1; 
    for (i=0; i<pdinfo->n; i++) {
	pdinfo->S[i] = malloc(9);
	if (pdinfo->S[i] == NULL) {
	    free(pdinfo->S);
	    return 1; 
	}
    }
    return 0;
}

/* ................................................. */

static int dataset_allocate_varnames (DATAINFO *pdinfo)
{
    int i, v = pdinfo->v;
    
    pdinfo->varname = malloc(v * sizeof(char *));
    pdinfo->label = malloc(v * sizeof(char *));
    if (pdinfo->varname == NULL || pdinfo->label == NULL) return 1;
    for (i=0; i<v; i++) {
	pdinfo->varname[i] = malloc(9);
	if (pdinfo->varname[i] == NULL) return 1;
	pdinfo->label[i] = malloc(MAXLABEL);
	if (pdinfo->label[i] == NULL) return 1;
	pdinfo->label[i][0] = '\0';
    }
    strcpy(pdinfo->varname[0], "const");
    strcpy(pdinfo->label[0], "auto-generated constant");
    return 0;
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

DATAINFO *create_new_dataset (double **pZ,  
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
    pdinfo->markers = (short) markers;
    if (pdinfo->markers) {
	if (dataset_allocate_markers(pdinfo)) {
	    free_datainfo(pdinfo);
	    return NULL;
	}
    } 
    dataset_dates_defaults(pdinfo);
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

int start_new_Z (double **pZ, DATAINFO *pdinfo, int resample)
{
    if (prepZ(pZ, pdinfo)) return 1;

    pdinfo->t1 = 0; 
    pdinfo->t2 = pdinfo->n - 1;

    if (resample) {
	pdinfo->varname = NULL;
	pdinfo->label = NULL;
    } else {
	if (dataset_allocate_varnames(pdinfo))
	    return 1;
    }

    pdinfo->S = NULL;
    pdinfo->markers = 0;
    
    return 0;
}

/* ................................................. */

static int prepZ (double **pZ, const DATAINFO *pdinfo)
    /* allocate data array and put in constant */
{
    int t;

    if (*pZ != NULL) free(*pZ);
    *pZ = malloc(pdinfo->v * pdinfo->n * sizeof **pZ);
    if (*pZ == NULL) return 1;

    for (t=0; t<pdinfo->n; t++) (*pZ)[t] = 1.0; 
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

static int readdata (FILE *fp, const DATAINFO *pdinfo, double *Z)
{
    int i, t, n = pdinfo->n;
    char c, marker[9];
    float x;

    gretl_errmsg[0] = '\0';

    if (pdinfo->bin == 1) { /* single-precision binary data */
	for (i=1; i<pdinfo->v; i++) {
	    for (t=0; t<n; t++) {
		if (!fread(&x, sizeof(float), 1, fp)) {
		    sprintf(gretl_errmsg, "WARNING: binary data read error at "
			    "var %d\n", i);
		    return 1;
		}
		Z(i, t) = (double) x;
	    }
	}
    }
    else if (pdinfo->bin == 2) { /* double-precision binary data */
	for (i=1; i<pdinfo->v; i++) {
	    if (!fread(&Z(i, 0), sizeof(double), n, fp)) {
		sprintf(gretl_errmsg, 
			"WARNING: binary data read error at var %d\n", i);
		return 1;
	    }
	}
    } else { /* ascii data file */
	for (t=0; t<n; t++) {
	    eatspace(fp);
	    c = fgetc(fp);  /* test for a #-opened comment line */
	    if (c == '#') while (c != '\n') c = fgetc(fp);
	    else ungetc(c, fp);
	    if (pdinfo->markers) {
		fscanf(fp, "%s", marker);
		strcpy(pdinfo->S[t], marker);
	    }
	    for (i=1; i<pdinfo->v; i++) {
		if ((fscanf(fp, "%lf", &Z(i, t))) != 1) {
		    sprintf(gretl_errmsg, 
			    "WARNING: ascii data read error at var %d, "
			    "obs %d\n", i, t + 1);
		    return 1;
		}
	    }
	}
    }
    return 0;
}

/* ................................................. */

static int gz_readdata (gzFile fgz, const DATAINFO *pdinfo, double *Z)
{
    int i, t, n = pdinfo->n;
    
    gretl_errmsg[0] = '\0';

    if (pdinfo->bin == 1) { /* single-precision binary data */
	float xx;

	for (i=1; i<pdinfo->v; i++) {
	    for (t=0; t<n; t++) {
		if (!gzread(fgz, &xx, sizeof xx)) {
		    sprintf(gretl_errmsg, "WARNING: binary data read error at "
			    "var %d", i);
		    return 1;
		}
		Z(i, t) = (double) xx;
	    }
	}
    }
    else if (pdinfo->bin == 2) { /* double-precision binary data */
	for (i=1; i<pdinfo->v; i++) {
	    if (!gzread(fgz, &Z(i, 0), n * sizeof(double))) {
		sprintf(gretl_errmsg, 
			"WARNING: binary data read error at var %d", i);
		return 1;
	    }
	}
    } else { /* ascii data file */
	char *line, numstr[24];
	int llen = pdinfo->v * 32;
	size_t offset;

	line = malloc(llen);
	if (line == NULL) return E_ALLOC;
	for (t=0; t<n; t++) {
	    offset = 0L;
	    if (!gzgets(fgz, line, llen - 1)) {
		sprintf(gretl_errmsg, "WARNING: ascii data read error at "
			"obs %d", t + 1);
		free(line);
		return 1;
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
			   "WARNING: failed to read case marker for "
			   "obs %d", t + 1);
		   free(line);
		   return 1;
		}
		pdinfo->S[t][8] = 0;
		offset += strlen(pdinfo->S[t]) + 1;
	    }
	    for (i=1; i<pdinfo->v; i++) {
		if (sscanf(line + offset, "%23s", numstr) != 1) {
		    sprintf(gretl_errmsg, 
			    "WARNING: ascii data read error at var %d, "
			    "obs %d", i, t + 1);
		    return 1;
		}
		numstr[23] = 0;
		Z(i, t) = atof(numstr);
		if (i < pdinfo->v - 1)
		    offset += strlen(numstr) + 1;
	    }
	}
	free(line);
    }
    return 0;
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
    static char varerr[] = "Reading data header file\n";

    gretl_errmsg[0] = '\0';
    
    if (!(isalpha((unsigned char) varname[0]))) {
        sprintf(gretl_errmsg, "%sfirst char of varname ('%c') is bad\n"
               "(first must be alphabetical)", varerr, varname[0]);
        return 1;
    }
    for (i=1; i<n; i++) {
        if (!(isalpha((unsigned char) varname[i]))  
            && !(isdigit((unsigned char) varname[i]))
            && varname[i] != '_') {
	    if (isprint((unsigned char) varname[i]))
		sprintf(gretl_errmsg, "%svarname contains illegal character '%c'\n"
			"Use only letters, digits and underscore", 
			varerr, varname[i]);
	    else
		sprintf(gretl_errmsg, "%svarname contains illegal character 0x%x\n"
			"Use only letters, digits and underscore", 
			varerr, (unsigned) varname[i]);
            return 1;
        }
    }
    return 0;
}       

/* ................................................ */

static int readhdr (const char *hdrfile, DATAINFO *pdinfo)
{
    FILE *fp;
    int n, i = 0, panel = 0;
    char str[MAXLEN], byobs[6], option[8];

    gretl_errmsg[0] = '\0';

    fp = fopen(hdrfile, "r");
    if (fp == NULL) {
	sprintf(gretl_errmsg, "\nCouldn't open data header file %s\n",  hdrfile);
	return E_FOPEN;
    }
    fscanf(fp, "%s", str);
    i += skipcomments(fp, str); 
    while (1) { /* find number of variables */
        if (fscanf(fp, "%s", str) != 1) {
	    fclose(fp);
	    sprintf(gretl_errmsg, "\nOpened header file %s\n"
		    "Couldn't find list of variables (must "
		    "be terminated with a semicolon)", hdrfile);
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
    pdinfo->sd0 = atof(pdinfo->stobs);
    pdinfo->n = dateton(pdinfo->endobs, pdinfo->pd, 
			pdinfo->stobs) + 1;
    pdinfo->extra = 0;    

    if (pdinfo->sd0 >= 2.0) 
        pdinfo->time_series = TIME_SERIES; /* actual time series? */
    else if (pdinfo->sd0 > 1.0)
	pdinfo->time_series = STACKED_TIME_SERIES; /* panel data? */
    else pdinfo->time_series = 0;

    pdinfo->bin = 0;
    pdinfo->markers = 0;
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
    return 0;

    varname_error:
    fclose(fp);
    clear_datainfo(pdinfo, 0);
    return E_DATA;
}

/* .......................................................... */

static int check_date (const char *date)
{
    int i, n = strlen(date);

    gretl_errmsg[0] = 0;

    for (i=0; i<n; i++) {
	if (!isdigit((unsigned char) date[i]) && date[i] != '.' 
	    && date[i] != ':') {
	    if (isprint((unsigned char) date[i]))
		sprintf(gretl_errmsg, 
			"Bad character '%c' in date string", date[i]);
	    else 
		sprintf(gretl_errmsg, 
			"Bad character %d in date string", date[i]);
	    return 1;
	}
    }
    return 0;
}

/**
 * dateton:
 * @date: string representation of date for processing.
 * @pd: periodicity or frequency of data.
 * @startdate: string representing starting date.
 * 
 * Given a "current" date string, a periodicity, and a starting
 * date string, returns the observation number corresponding to
 * the current date string, counting from zero.
 * 
 * Returns: integer observation number.
 *
 */

int dateton (const char *date, const int pd, const char *startdate)
{
    int dotpos1 = 0, dotpos2 = 0, maj = 0, min = 0, n, i;
    char majstr[5], minstr[3];
    char startmajstr[5], startminstr[3];
    int startmaj, startmin;

    if (check_date(date)) return -1;

    n = strlen(date);
    for (i=1; i<n; i++) {
        if (date[i] == '.' || date[i] == ':') {
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
    n = strlen(startdate);
    for (i=1; i<n; i++) {
        if (startdate[i] == '.' || startdate[i] == ':') {
	    dotpos2 = i;
	    break;
	}
    }
    if ((dotpos1 && !dotpos2) || (dotpos2 && !dotpos1)) {
	sprintf(gretl_errmsg, "date strings inconsistent");
	return -1;  
    }
    if (!dotpos1 && !dotpos2) {
        return (atoi(date) - atoi(startdate));
    }
    safecpy(startmajstr, startdate, dotpos2);
    startmaj = atoi(startmajstr);
    strcpy(startminstr, startdate + dotpos2 + 1);
    startmin = atoi(startminstr);
    n = pd * (maj - startmaj);
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
 */

void ntodate (char *datestr, const int nt, const DATAINFO *pdinfo)
/* print to datestr the calendar representation of nt */
{
    double xn;
    int n;

    xn = date(nt, pdinfo->pd, pdinfo->sd0);
    if (pdinfo->pd == 1) {
        n = (int) xn;
        sprintf(datestr, "%d", n);
    }
    else if (pdinfo->pd < 10) sprintf(datestr, "%.1f", xn);
    else sprintf(datestr, "%.2f", xn);
}

/* .......................................................... */

static int blank_check (FILE *fp)
{
    int i, deflt = 1;
    char s[MAXLEN];

    for (i=0; i<3 && deflt && fgets(s, MAXLEN-1, fp) ; i++) {
	if (i == 0 && strncmp(s, "(*", 2)) deflt = 0;
	else if (i == 1 && strncmp(s, "space for comments", 18)) deflt = 0;
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
	pprintf(prn, "Couldn't open %s\n", hdrfile); 
	return 1;
    }

    /* see if it's just the default "space for comments" */
    if (blank_check(hdr)) { /* yes */
	pprintf(prn, "No info in %s\n", hdrfile);
	return 2;
    } 

    /* no, so restart the read */
    if ((hdr = fopen(hdrfile, "r")) == NULL) {
	pprintf(prn, "Couldn't open %s\n", hdrfile); 
	return 1;
    }    

    pprintf(prn, "Data info in file %s:\n\n", hdrfile);
    if (fgets(s, MAXLEN-1, hdr) != NULL && strncmp(s, STARTCOMMENT, 2) == 0) {
	do {
	    if (fgets(s, MAXLEN-1, hdr) != NULL && strncmp(s, "*)", 2)) {
#ifndef OS_WIN32
		delchar('\r', s);
#endif
		pprintf(prn, "%s", s);
		i++;
	    }
	} while (s != NULL && strncmp(s, ENDCOMMENT, 2));
    }
    if (i == 0) pprintf(prn, " (none)\n");
    pprintf(prn, "\n");

    if (hdr != NULL) fclose(hdr);
    return 0;
}

/* ................................................ */

static int writehdr (const char *hdrfile, const int *list, 
		     const DATAINFO *pdinfo, const int opt)
{
    FILE *fp;
    int bin = 0, i;
    int old_comments = 0;
    char *comment_buf = NULL;
    char startdate[8], enddate[8];

    if (opt == GRETL_DATA_FLOAT) bin = 1;
    else if (opt == GRETL_DATA_DOUBLE) bin = 2;

    ntodate(startdate, pdinfo->t1, pdinfo);
    ntodate(enddate, pdinfo->t2, pdinfo);

    /* If there's already a file of this name, and it contains
       comments on the data set, then preserve them */
    fp = fopen(hdrfile, "r");
    if (fp != NULL) {
	old_comments = comment_lines(fp, &comment_buf);
	if (old_comments < 0) return 1;
	if (fp != NULL) fclose (fp);
    }
	
    fp = fopen(hdrfile, "w");
    if (fp == NULL) return 1;

    /* write original comments, if any */
    if (old_comments) 
	fprintf(fp, "%s", comment_buf);
    else
	fprintf(fp, "(*\nspace for comments\n*)\n");

    /* then list of variables */
    for (i=1; i<=list[0]; i++) {
	if (list[i] == 0) continue;
	fprintf(fp, "%s ", pdinfo->varname[list[i]]);
	if (i && i <list[0] && (i+1) % 8 == 0) fprintf(fp, "\n");
    }    
    fputs(";\n", fp);

    /* then obs line */
    fprintf(fp, "%d ", pdinfo->pd);
    fprintf(fp, "%s ", startdate);
    fprintf(fp, "%s\n", enddate);
    
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
    if (comment_buf != NULL) free(comment_buf);
    return 0;
}

/**
 * _get_precision:
 * @x: data vector.
 * @n: length of x.
 *
 * Find the number of decimal places required to represent a given
 * data series uniformly.
 * 
 * Returns: the required number of decimal places.
 *
 */

int _get_precision (double *x, int n)
{
    int i, j, p, dot, len, pmax = 0;
    char numstr[48];

    for (i=0; i<n; i++) {
	sprintf(numstr, "%f", x[i]);
	numstr[31] = '\0';
	len = strlen(numstr);
	dot = dotpos(numstr);
	p = len - dot - 1;
	if (p > 0) { /* figures after decimal point */
	    for (j=len-1; j>1; j--) {
		if (numstr[j] == '0') p--;
		else break;
	    }
	    if (p > pmax) pmax = p;
	}
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
 * 
 * Write out a data file containing the values of the given set
 * of variables.
 * 
 * Returns: 0 on successful completion, non-zero on error.
 * 
 */

int write_data (const char *fname, const int *list, 
		double *Z, const DATAINFO *pdinfo, int opt)
{
    int i = 0, t, l0 = list[0], n = pdinfo->n;
    char datfile[MAXLEN], hdrfile[MAXLEN], lblfile[MAXLEN];
    FILE *fp = NULL;
    int *pmax = NULL, tsamp = pdinfo->t2 - pdinfo->t1 + 1;
    double xx;

    strcpy(datfile, fname);

    if (opt == GRETL_DATA_GZIPPED && !has_gz_suffix(datfile))
	strcat(datfile, ".gz");

    if (opt == GRETL_DATA_R && pdinfo->time_series == TIME_SERIES) 
	opt = GRETL_DATA_R_ALT;

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
	    fprintf(stderr, "Write of header file failed");
	    return 1;
	}
	if (writelbl(lblfile, list, pdinfo)) {
	    fprintf(stderr, "Write of labels file failed");
	    return 1;
	}
    }

    /* open files, other than for gzipped output */
    if (opt == GRETL_DATA_FLOAT || opt == GRETL_DATA_DOUBLE) 
	fp = fopen(datfile, "wb");
    else if (opt != GRETL_DATA_GZIPPED) 
	fp = fopen(datfile, "w");
    if (opt != GRETL_DATA_GZIPPED && fp == NULL) return 1;

    if (opt == GRETL_DATA_FLOAT) { /* single-precision binary */
	float x;

	for (i=1; i<=l0; i++) {
	    for (t=0; t<n; t++) {
		x = (float) Z(list[i], t);
		fwrite(&x, sizeof(float), 1, fp);
	    }
	}
    }
    else if (opt == GRETL_DATA_DOUBLE) { /* double-precision binary */
	for (i=1; i<=l0; i++) 
	    fwrite(&Z(list[i], 0), sizeof(double), n, fp);
    }

    if (opt == GRETL_DATA_GZIPPED || opt == GRETL_DATA_CSV || 
	opt == GRETL_DATA_OCTAVE || GRETL_DATA_R || opt == 0) { 
	/* an ASCII variant of some sort */
	pmax = malloc(l0 * sizeof *pmax);
	if (pmax == NULL) return 1;
	for (i=1; i<=l0; i++) 
	    pmax[i-1] = _get_precision(&Z(list[i], pdinfo->t1), tsamp);
    }

    if (opt == GRETL_DATA_GZIPPED) { /* compressed ASCII */
	gzFile fgz;

	fgz = gzopen(datfile, "wb");
	if (fgz == Z_NULL) return 1;
	for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	    if (pdinfo->markers && pdinfo->S != NULL) 
		gzprintf(fgz, "%s ", pdinfo->S[t]);
	    for (i=1; i<=l0; i++) {
		if (na(Z(list[i], t)))
		    gzprintf(fgz, "-999 ");
		else 
		    gzprintf(fgz, "%.*f ", pmax[i-1], Z(list[i], t));
	    }
	    gzprintf(fgz, "\n");
	}
	gzprintf(fgz, "\n");
	gzclose(fgz);
    }
    else if (opt == 0) { /* plain ASCII */
	for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	    if (pdinfo->markers && pdinfo->S != NULL) 
		fprintf(fp, "%s ", pdinfo->S[t]);
	    for (i=1; i<=l0; i++) {
		if (na(Z(list[i], t)))
		    fprintf(fp, "-999 ");
		else 
		    fprintf(fp, "%.*f ", pmax[i-1], Z(list[i], t));
	    }
	    fputs("\n", fp);
	}
	fputs("\n", fp);
    }
    else if (opt == GRETL_DATA_CSV || opt == GRETL_DATA_R) { 
	/* export CSV or GNU R (dataframe) */
	int idate;
	double xdate;
	char comma[2];
	
	if (opt == GRETL_DATA_CSV) strcpy(comma, ",");
	else strcpy(comma, " ");

	/* variable names */
	if (opt == GRETL_DATA_CSV) fprintf(fp, "obs,");
	for (i=1; i<l0; i++) 
	    fprintf(fp, "%s%s", pdinfo->varname[list[i]], comma);
	fprintf(fp, "%s\n", pdinfo->varname[list[l0]]);
	
	for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	    if (pdinfo->S != NULL) 
		fprintf(fp, "%s%s", pdinfo->S[t], comma);
	    else {
		xdate = date(t, pdinfo->pd, pdinfo->sd0);
		idate = (int) xdate;
		if (pdinfo->pd == 1) 
		    fprintf(fp, "\"%d\"%s", idate, comma);
		else if (pdinfo->pd < 10) 
		    fprintf(fp, "\"%.1f\"%s", xdate, comma);
		else fprintf(fp, "\"%.2f\"%s", xdate, comma);
	    }
	    for (i=1; i<=l0; i++) { 
		xx = Z(list[i], t);
		if (na(xx))
		    fprintf(fp, "NA");
		else
		    fprintf(fp, "%.*f", pmax[i-1], Z(list[i], t));
		if (i < l0) fprintf(fp, "%s", comma);
		else fprintf(fp, "\n");
	    }
	}
	fputs("\n", fp);
    }
    else if (opt == GRETL_DATA_R_ALT) { /* export GNU R (structure) */
	for (i=1; i<=l0; i++) {
	    fprintf(fp, "\"%s\" <-\n", pdinfo->varname[list[i]]);
	    fprintf(fp, "structure(c(");
	    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
		xx = Z(list[i], t);
		if (na(xx)) fprintf(fp, "NA");
		else fprintf(fp, "%g", Z(list[i], t));
		if (t < pdinfo->t2) fprintf(fp, ", "); 
		if (t > pdinfo->t1 && (t - pdinfo->t1) % 8 == 0 &&
		    t < pdinfo->t2)
		    fprintf(fp, "\n");
	    }
	    fprintf(fp, ")");
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
	for (t=pdinfo->t1; t<=pdinfo->t2; t++) 
	    fprintf(fp, "%.*f\n", pmax[0], Z(list[1], t));
	/* write out info for indep vars matrix */
	fprintf(fp, "# name: X\n# type: matrix\n# rows: %d\n# columns: %d\n", 
		n, list[0] - 1);
	/* write out indep. var. matrix */
	for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	    for (i=2; i<=list[0]; i++) {
		fprintf(fp, "%.*f ", pmax[i-1], Z(list[i], t));
	    }
	    fputs("\n", fp);
	}
	fputs("\n", fp);
    }

    if (pmax) free(pmax);
    if (fp != NULL) fclose(fp);
    return 0;
}

/* ................................................. */

static double obs_float (const DATAINFO *pdinfo, const int end)
{
    double xx, xx2 = 0.;
    int i, x1, x2 = 0;

    if (end) {
	xx = atof(pdinfo->endobs);
	if ((i = haschar('.', pdinfo->endobs)) > 0)
	   x2 = atoi(pdinfo->endobs + i + 1) - 1;
    } else {
	xx = atof(pdinfo->stobs);
	if ((i = haschar('.', pdinfo->stobs)) > 0)
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
    char line[MAXLEN], *label, varname[9];
    int len, v;
    
    gretl_errmsg[0] = '\0';

    fp = fopen(lblfile, "r");
    if (fp == NULL) return 0;
    while (1) {
        if (fgets(line, MAXLEN-1, fp) == NULL) {
            fclose(fp);
            return 0;
        }
        if (sscanf(line, "%s", varname) != 1) {
            fclose(fp);
	    sprintf(gretl_errmsg, "Bad data label in %s", lblfile); 
            return 0;
        }
        len = strlen(varname);
        label = line + len;
        if (top_n_tail(label) == E_ALLOC) {
            fclose(fp);
            return E_ALLOC;
        }
	v = varindex(pdinfo, varname);
	if (v < pdinfo->v) strcpy(pdinfo->label[v], label);
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
	if (strlen(pdinfo->label[list[i]]) > 2) {
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
	if (strlen(pdinfo->label[list[i]]) > 2) 
	    fprintf(fp, "%s %s\n", pdinfo->varname[list[i]],
		    pdinfo->label[list[i]]);
    }    
    if (fp != NULL) fclose(fp);
    return 0;
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
    int i = dotpos(src), j = slashpos(src), k;

    strcpy(targ, src);
    targ[i] = '\0';
    k = dotpos(targ);
    if (j > 0 && k < strlen(targ) && k > j) i = k;
    targ[i] = '.';
    targ[i + 1] = '\0';
    strcat(targ, ext);
}

/* ................................................ */

static int try_gdt (char *fname)
{
    char *suff;
    int ret = 0;

    if (fname != NULL) {
	suff = strrchr(fname, '.');
	if (suff != NULL && !strcmp(suff, ".dat")) {
	    strcpy(suff, ".gdt");
	    ret = 1;
	}
    }
    return ret;
}

/**
 * get_data:
 * @pZ: pointer to data set.
 * @pdinfo: data information struct.
 * @datfile: name of file to try.
 * @ppaths: path information struct.
 * @data_status: indicator for whether a data file is currently open
 * in gretl's work space (not-0) or not (0).
 * @fp: file pointer from which data may be read.
 * 
 * Read data from file into gretl's work space, allocating space as
 * required.
 * 
 * Returns: 0 on successful completion, non-zero otherwise.
 *
 */

int get_data (double **pZ, DATAINFO *pdinfo, char *datfile, PATHS *ppaths, 
	      const int data_status, FILE *fp) 
{

    FILE *dat = NULL;
    gzFile fgz = NULL;
    int err, gzsuff = 0;
    char hdrfile[MAXLEN], lblfile[MAXLEN];

    gretl_errmsg[0] = '\0';

    /* get filenames organized */
    hdrfile[0] = '\0';
    gzsuff = has_gz_suffix(datfile);

    if (addpath(datfile, ppaths, 0) == NULL) { /* not found yet */
	char tryfile[MAXLEN];
	int found = 0;

	if (!gzsuff) { /* maybe the file is gzipped but we 
			  didn't get the .gz extension? */
	    sprintf(tryfile, "%s.gz", datfile);
	    if (addpath(tryfile, ppaths, 0) != NULL) {
		gzsuff = 1;
		found = 1;
	    }
	}
	
	if (!found) { /* no -- then try using the .gdt suffix? */
	    strcpy(tryfile, datfile);
	    if (try_gdt(tryfile)) 
		found = (addpath(tryfile, ppaths, 0) != NULL);
	} 
	    
	if (!found) return E_FOPEN;
	else strcpy(datfile, tryfile);
    }
	
    if (!gzsuff) {
	switch_ext(hdrfile, datfile, "hdr");
	switch_ext(lblfile, datfile, "lbl");
    } else {
	gz_switch_ext(hdrfile, datfile, "hdr");
	gz_switch_ext(lblfile, datfile, "lbl");
    }

    /* clear any existing data info */
    if (data_status) clear_datainfo(pdinfo, 0);

    /* read data header file */
    err = readhdr(hdrfile, pdinfo);
    if (err) return err;
    else 
	fprintf(fp, "\nReading header file %s\n", hdrfile);

    /* deal with case where first col. of data file contains
       "marker" strings */
    if (pdinfo->markers) {
	if (dataset_allocate_markers(pdinfo)) return E_ALLOC; 
    } else pdinfo->S = NULL;
    
    /* allocate dataset */
    if (prepZ(pZ, pdinfo)) return E_ALLOC;

    /* Invoke data (Z) reading function */
    if (gzsuff) {
	fgz = gzopen(datfile, "rb");
	if (fgz == NULL) return E_FOPEN;
    } else {
	if (pdinfo->bin)
	    dat = fopen(datfile, "rb");
	else
	    dat = fopen(datfile, "r");
	if (dat == NULL) return E_FOPEN;
    }

    /* print out basic info from the files read */
    fprintf(fp, "periodicity: %d, maxobs: %d, "
	   "observations range: %s-%s\n", pdinfo->pd, pdinfo->n,
	   pdinfo->stobs, pdinfo->endobs);

    fprintf(fp, "\nReading ");
    fprintf(fp, (pdinfo->time_series == TIME_SERIES) ? 
	    "time-series" : "cross-sectional");
    fprintf(fp, " datafile");
    if (strlen(datfile) > 40) putc('\n', fp);
    fprintf(fp, " %s\n\n", datfile);

    if (gzsuff) {
	err = gz_readdata(fgz, pdinfo, *pZ); 
	gzclose(fgz);
    } else {
	err = readdata(dat, pdinfo, *pZ); 
	fclose(dat);
    }
    if (err) return err;

    /* Set sample range to entire length of dataset by default */
    pdinfo->t1 = 0; 
    pdinfo->t2 = pdinfo->n - 1;

    strcpy(ppaths->datfile, datfile);
    strcpy(ppaths->hdrfile, hdrfile);
    strcpy(ppaths->lblfile, lblfile);

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

int open_nulldata (double **pZ, DATAINFO *pdinfo, 
		   const int data_status, const int length,
		   PRN *prn) 
{
    /* clear any existing data info */
    if (data_status) clear_datainfo(pdinfo, 0);

    /* dummy up the data info */
    pdinfo->n = length;
    pdinfo->v = 1;
    dataset_dates_defaults(pdinfo);

    if (dataset_allocate_varnames(pdinfo)) return E_ALLOC;

    /* no observation markers */
    pdinfo->markers = 0;
    pdinfo->S = NULL; 

    /* allocate dataset */
    if (prepZ(pZ, pdinfo)) return E_ALLOC;

    /* print out basic info */
    pprintf(prn, "periodicity: %d, maxobs: %d, "
	   "observations range: %s-%s\n", pdinfo->pd, pdinfo->n,
	   pdinfo->stobs, pdinfo->endobs);

    /* Set sample range to entire length of data-set by default */
    pdinfo->t1 = 0; 
    pdinfo->t2 = pdinfo->n - 1;

    return 0;
}

/* ......................................................... */

static int test_label (DATAINFO *pdinfo, PRN *prn)
     /* attempt to parse csv row labels as dates.  Return -1 if
	this doesn't work out, 0 if the labels seem to be just
	integer observation numbers, else return the inferred data
	frequency */
{
    int n1, n2, try;
    char year[5], subper[3], endobs[9];
    char lbl1[9], lbl2[9];

    strcpy(lbl1, pdinfo->S[0]);
    strcpy(lbl2, pdinfo->S[pdinfo->n - 1]);
    n1 = strlen(lbl1);
    n2 = strlen(lbl2);

    pprintf(prn, "   first row label \"%s\", last label \"%s\"\n", 
	   lbl1, lbl2);

    /* are the labels (probably) just 1, 2, 3 etc.? */
    sprintf(endobs, "%d", pdinfo->n);
    if (strcmp(lbl1, "1") == 0 && strcmp(lbl2, endobs) == 0)
	return 0;

    if (n1 > 7) {
	pprintf(prn, "   label strings too long for dates?\n");
	pdinfo->pd = 1;
	pdinfo->sd0 = 1.0;
	return -1;
    }
    if (n1 != n2) {
	pprintf(prn, "   label strings can't be consistent dates\n");
	return -1;
    }

    /* does it look like it starts with a year? */
    pprintf(prn, "trying to parse row labels as dates...\n");
    if (n1 >= 4) {
	if (isdigit((unsigned char) lbl1[0]) 
	    && isdigit((unsigned char) lbl1[1]) &&
	    isdigit((unsigned char) lbl1[2]) && 
	    isdigit((unsigned char) lbl1[3])) {
	    safecpy(year, lbl1, 4);
	    try = atoi(year);
	    if (try > 0 && try < 3000) {
		pprintf(prn, "   %s: probably a year... ", year);
	    } else {
		pprintf(prn, "   %s: out of bounds for a year?\n", year);
	    }
	    if (n1 == 5) {
		pprintf(prn, "   but I can't make sense of the extra bit\n");
		return -1;
	    }
	    if (n1 == 4) {
		pprintf(prn, "and just a year\n");
		strcpy(pdinfo->stobs, year);
		pdinfo->sd0 = atof(pdinfo->stobs);
		/* need more checking!! FIXME */
		strcpy(pdinfo->endobs, lbl2);
		pdinfo->pd = 1;
		return 1;
	    }
	    if (lbl1[4] == '.' || lbl1[4] == ':' || lbl1[4] == 'Q') {
		strcpy(subper, lbl1+5);
		if (n1 == 6) {
		    pprintf(prn, "quarter %s?\n", subper);
		    sprintf(pdinfo->stobs, "%s.%s", year, subper);
		    pdinfo->sd0 = atof(pdinfo->stobs);
		    strncpy(pdinfo->endobs, lbl2, 4);
		    pdinfo->endobs[4] = '.';
		    strcat(pdinfo->endobs, lbl2+5);
		    pdinfo->pd = 4;
		    return 4;
		}
		if (n1 == 7) {
		    pprintf(prn, "month %s?\n", subper);
		    sprintf(pdinfo->stobs, "%s.%s", year, subper);
		    pdinfo->sd0 = atof(pdinfo->stobs);
		    strncpy(pdinfo->endobs, lbl2, 4);
		    pdinfo->endobs[4] = '.';
		    strcat(pdinfo->endobs, lbl2+5);
		    pdinfo->pd = 12;
		    return 12;
		}
	    }
	} else pprintf(prn, "   definitely not a four-digit year\n");
    }

    return -1;
}

/* ......................................................... */

static int check_csv_merge (DATAINFO *pdinfo, DATAINFO *pcinfo, 
			    PRN *prn)
{
    pprintf(prn, "Checking for conformability with present data set...\n");
    if (pdinfo->n != pcinfo->n) {
	pprintf(prn, "   Number of observations does not match.\n");
	return 1;
    }
    if (pdinfo->pd != pcinfo->pd) {
	pprintf(prn, "   Frequency does not match.\n");
	return 1;
    }
    if (strcmp(pdinfo->stobs, pcinfo->stobs)) {
	pprintf(prn, "   Starting observation does not match.\n");
	return 1;
    } 
    if (strcmp(pdinfo->endobs, pcinfo->endobs)) {
	pprintf(prn, "   Ending observation does not match.\n");
	return 1;
    }
    pprintf(prn, "   OK.\n");
    return 0;
}

/* ......................................................... */

static int do_csv_merge (DATAINFO *pdinfo, DATAINFO *pcinfo,
			 double **pZ, double **csvZ, PRN *prn)
{
    int i, t, newvars = pcinfo->v, oldvars = pdinfo->v;
    int n = pdinfo->n;

    pprintf(prn, "Attempting data merge...\n");
    if (_grow_Z(newvars - 1, pZ, pdinfo)) {
	pprintf(prn, "   Out of memory.\n");
	return E_ALLOC;
    }
    for (i=1; i<newvars; i++) {
	for (t=0; t<pdinfo->n; t++) 
	    (*pZ)[n*(oldvars+i-1) + t] = (*csvZ)[n*i + t];
	strcpy(pdinfo->varname[oldvars+i-1], pcinfo->varname[i]);
    }  
    clear_datainfo(pcinfo, 0);
    free(*csvZ);
    pprintf(prn, "   OK, I think.\n");
    return 0;
}

/* ......................................................... */

static void trim_csv_line (char *line)
{
    size_t n;

    if (line == NULL) return;
    n = strlen(line);
    if (n == 0) return;

    if (line[n-1] == '\n') line[n-1] = '\0';
    n = strlen(line);
    if (n == 0) return;
    if (line[n-1] == '\r') line[n-1] = '\0';
    n = strlen(line);
    if (n == 0) return;
    if (line[n-1] == ',') line[n-1] = '\0';
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

int import_csv (double **pZ, DATAINFO *pdinfo, 
		const char *fname, PRN *prn)
{
    int n, nv, missval = 0, ncols = 0, chkcols = 0;
    char c, cbak;
    int bad_commas = 0, skipvar = 0, markertest = -1;
    int i, j, k, t, blank_1 = 0, obs_1 = 0, len = 0, maxlen = 0, ok = 0;
    char *line, varname[9], numstr[32], field_1[32];
    FILE *fp;
    DATAINFO csvinfo;
    double *csvZ = NULL;
    const char *msg = "\nPlease note:\n"
	"- The first row of the CSV file should contain the "
	"names of the variables.\n"
	"- The first column may optionally contain date "
	"strings or other 'markers':\n  in that case its row 1 entry "
	"should be blank, or should say 'obs' or 'date'.\n"
	"- The remainder of the file must be a rectangular "
	"array of data.\n";

    fp = fopen(fname, "r");
    if (fp == NULL) {
	pprintf(prn, "couldn't open %s\n", fname);
	return 1;
    }

    csvinfo.markers = 0;

    pprintf(prn, "parsing %s...\n", fname);

    /* count chars and fields in first line */
    if (fread(&cbak, 1, 1, fp) == 0 || cbak == '\n') {
	pprintf(prn, "   empty first line!\n");
	fclose(fp);
	return 1;
    }
    if (cbak == ',') {
	blank_1 = 1;
	pprintf(prn, "   first field is blank (dates?)\n");
	ncols++;
    }
    maxlen++;
    while (fread(&c, 1, 1, fp)) {
	if ((c == '\n' || c == '\r') && cbak == ',') {
	    bad_commas = 1;
	    pprintf(prn, "   file has trailing commas (lame)\n");
	}
	if (c == '\n') break;
	cbak = c;
	maxlen++;
	if (c == ',') ncols++;
    }
    if (!bad_commas) ncols++;

    pprintf(prn, "   number of columns = %d\n", ncols);

    /* now count remaining non-blank rows, checking for fields */
    csvinfo.n = 0;
    chkcols = (bad_commas)? -1: 0;
    while (fread(&c, 1, 1, fp)) {
	if (!(isspace((unsigned char) c))) ok = 1;
	if (c == ',') chkcols += 1;
	if (c != '\n') len++;
	else {
	    if (len > maxlen) maxlen = len;
	    len = 0;
	    if (ok) {
		chkcols += 1; 
		csvinfo.n += 1;
		if (chkcols != ncols) {
		    pprintf(prn, "   ...but row %d has %d fields: aborting\n",
			    csvinfo.n, chkcols);
		    fclose(fp);
		    pprintf(prn, msg);
		    return 1;
		}
	    }
	    ok = 0; 
	    chkcols = (bad_commas)? -1: 0;
	}
    }
    pprintf(prn, "   longest line: %d characters\n", maxlen + 1);

    if (!blank_1) {
	rewind(fp);
	c = 0; i = 0;
	while (c != ',') {
	    fread(&c, 1, 1, fp);
	    field_1[i++] = c;
	}
	field_1[i-1] = '\0';
	delchar('"', field_1);
	pprintf(prn, "   first field: '%s'\n", field_1);
	lower(field_1);
	if (strcmp(field_1, "obs") == 0 || strcmp(field_1, "date") == 0) {
	    pprintf(prn, "   seems to be observation label\n");
	    obs_1 = 1;
	    skipvar = 1;
	}
    }

    csvinfo.v = (blank_1 || obs_1)? ncols: ncols + 1;
    pprintf(prn, "   number of variables: %d\n", csvinfo.v - 1);
    pprintf(prn, "   number of non-blank lines: %d\n", 
	    csvinfo.n + 1);

    fclose(fp);
    /* end initial checking */

    /* initialize datainfo and Z */
    if (start_new_Z(&csvZ, &csvinfo, 0)) return E_ALLOC;

    if (blank_1 || obs_1) {
	csvinfo.markers = 1;
	if (dataset_allocate_markers(&csvinfo)) return E_ALLOC;
    }

    /* second pass */
    fp = fopen(fname, "r");
    if (fp == NULL) return 1;
    line = malloc(maxlen + 1);
    if (line == NULL) {
	fclose(fp);
	clear_datainfo(&csvinfo, 0);
	return E_ALLOC;
    }

    /* parse the variable names, truncating to 8 chars */
    pprintf(prn, "scanning for variable names...\n");
    fgets(line, maxlen + 1, fp);
    trim_csv_line(line);
    delchar(' ', line);
    delchar('"', line);
    pprintf(prn, "   line: %s\n", line);
    n = strlen(line);
    k = 0;
    for (i=0; i<n; i++) {
	while (line[i] == ',' || line[i] == '"' || line[i] == QUOTE) i++;
	j = 0; 
	while (line[i] != '"' && line[i] != QUOTE && line[i] != ',' && j < 8) 
	    varname[j++] = line[i++];
	varname[j] = '\0';
	k++;
	if (strlen(varname) == 0 || varname[0] == '\n') {
	    pprintf(prn, "   variable name %d is missing: aborting\n", k);
	    pprintf(prn, msg);
	    fclose(fp);
	    free(line);
	    clear_datainfo(&csvinfo, 0);
	    return 1;
	}
	if (k == 1 && skipvar) {
	    k--;
	    skipvar = 0;
	} else {
	    pprintf(prn, "   variable %d: '%s'\n", k, varname); 
	    strcpy(csvinfo.varname[k], varname);
	} 
	if (k == csvinfo.v - 1) break;
	if ((k) && j == 8 && line[i] != ',') 
	    while (line[i+1] != ',') i++;
    }

    pprintf(prn, "scanning for row labels and data...\n");
    for (t=0; t<csvinfo.n; t++) {
	ok = 0;
	fgets(line, maxlen + 1, fp);
	trim_csv_line(line);
	delchar(' ', line);
	delchar('"', line);
	n = strlen(line);
	for (i=0; i<n; i++) {
	    if (!(isspace((unsigned char) line[i]))) { 
		ok = 1; 
		break; 
	    }
	}
	if (!ok) { 
	    t--; 
	    continue; 
	}
	k = 0;
	/* printf("t=%d: ", t); */
	for (i=0; i<n; i++) {   /* parse line */
	    if (i != 0 || line[i] != ',') {
		if (k == 0 && line[i] == ',') i++;
		if (line[i] == '"' || line[i] == QUOTE) i++;
		j = 0; 
		while (line[i] != ',' && line[i] != '\n' 
		       && line[i] != '"' && line[i] != QUOTE
		       && line[i] != '\0') {
		    numstr[j++] = line[i++];
		}
		numstr[j] = '\0';
	    } else {
		numstr[0] = '\0';
	    }
	    if (strlen(numstr) == 0 || numstr[0] == '\n' ||
		strcmp(numstr, "NA") == 0) {
		if (numstr[0] == 'N') 
		    pprintf(prn, "   warning: missing value for variable "
			    "%d, obs %d\n", k, t+1);
		else 
		    pprintf(prn, "   the cell for variable %d, obs %d "
			    "is empty: treating as missing value\n", 
			    k, t+1);
		missval = 1;
	    } 
	    k++;
	    if (blank_1 || obs_1) nv = k - 1;
	    else nv = k;
	    if ((blank_1 || obs_1) && k == 1) {
		_esl_trunc(numstr, 8);
		strcpy(csvinfo.S[t], numstr);
	    } else {
		if (missval) 
		    csvZ[csvinfo.n * nv + t] = NADBL;
		else
		    csvZ[csvinfo.n * nv + t] = atof(numstr);
		missval = 0;
	    }
	    if (k == ncols) break;
	}
    }
    fclose(fp); 
    free(line);

    csvinfo.t1 = 0;
    csvinfo.t2 = csvinfo.n - 1;
    if (blank_1 || obs_1) markertest = test_label(&csvinfo, prn);
    if ((blank_1 || obs_1) && (markertest > 0))
	pprintf(prn, "taking date information from row labels\n\n");
    else {
	pprintf(prn, "treating these as undated data\n\n");
	dataset_dates_defaults(&csvinfo);
    }	
    if (csvinfo.pd != 1 || strcmp(csvinfo.stobs, "1")) 
        csvinfo.time_series = TIME_SERIES;

    /* If there were observation labels and they were not interpretable
       as dates, and they weren't simply "1, 2, 3, ...", then they 
       should probably be preserved. */

    if (csvinfo.S != NULL && markertest >= 0) {
	csvinfo.markers = 0;
	for (i=0; i<csvinfo.n; i++) free(csvinfo.S[i]);
	free(csvinfo.S);
	csvinfo.S = NULL;
    }

    if (*pZ == NULL) {
	*pZ = csvZ;
	*pdinfo = csvinfo;
    } else {
	if (check_csv_merge(pdinfo, &csvinfo, prn)) return 1;
	if (do_csv_merge(pdinfo, &csvinfo, pZ, &csvZ, prn)) return 1;
    }

    return 0;
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
    char **S, marker[9];
    int i, t;

    fp = fopen(fname, "r");
    if (fp == NULL) return E_FOPEN;
    
    S = malloc(pdinfo->n * sizeof *S);
    if (S == NULL) return E_ALLOC; 
    for (i=0; i<pdinfo->n; i++) {
	S[i] = malloc(9);
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
	for (i=0; i<pdinfo->n; i++) 
	    free (pdinfo->S[i]);
	free(pdinfo->S);
    }
    pdinfo->S = S;
    pdinfo->markers = 1;
    return 0;
}

/* ................................................. */

static char *unspace (char *s)
{
    size_t i, n = strlen(s);

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

int import_box (double **pZ, DATAINFO *pdinfo, 
		const char *fname, PRN *prn)
{
    int c, cc, i, t, v, realv, gotdata;
    int maxline, dumpvars;
    char tmp[48];
    unsigned *varsize, *varstart;
    char *test, *line;
    double x;
    FILE *fp;
    DATAINFO boxinfo;
    double *boxZ = NULL;
    extern int errno;

    fp = fopen(fname, "r");
    if (fp == NULL) {
	pprintf(prn, "couldn't open %s\n", fname);
	return 1;
    }

    pprintf(prn, "parsing %s...\n", fname);

    /* first pass: find max line length, number of vars and number
       of observations, plus basic sanity check */
    cc = maxline = 0;
    boxinfo.n = 0;
    boxinfo.v = 1;
    do {
	c = getc(fp); 
	if (c != EOF && c != 10 && !isprint((unsigned char) c)) {
	    pprintf(prn, "Binary data (%d) encountered: this is not a valid "
		   "BOX1 file\n", c);
	    fclose(fp);
	    return 1;
	}
	if (c == '\n') {
	    if (cc > maxline) maxline = cc;
	    cc = 0;
	    if ((c=getc(fp)) != EOF) {
		tmp[0] = c; cc++;
	    } else break;
	    if ((c=getc(fp)) != EOF) {
		tmp[1] = c; cc++;
	    } else break;
	    tmp[2] = '\0';
	    if (!strcmp(tmp, "03")) boxinfo.v += 1;
	    else if (!strcmp(tmp, "99")) boxinfo.n += 1;
	} else
	    cc++;
    } while (c != EOF);
    fclose(fp);

    pprintf(prn, "   found %d variables\n", boxinfo.v - 1);
    pprintf(prn, "   found %d observations\n", boxinfo.n);
    pprintf(prn, "   longest line = %d characters\n", maxline); 
    maxline += 2;

    /* allocate space for data etc */
    pprintf(prn, "allocating memory for data... ");
    if (start_new_Z(&boxZ, &boxinfo, 0)) return E_ALLOC;
    boxinfo.markers = 0;
    varstart = malloc((boxinfo.v - 1) * sizeof *varstart);
    if (varstart == NULL) return E_ALLOC;
    varsize = malloc((boxinfo.v - 1) * sizeof *varsize);
    if (varsize == NULL) return E_ALLOC;
    line = malloc(maxline);
    if (line == NULL) return E_ALLOC;
    pprintf(prn, "done\n");

    fp = fopen(fname, "r");
    if (fp == NULL) return 1;
    pprintf(prn, "reading variable information...\n");

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
	    strncpy(boxinfo.varname[realv], line+11, 8);
	    boxinfo.varname[realv][8] = '\0';
	    unspace(boxinfo.varname[realv]);
	    lower(boxinfo.varname[realv]);
	    pprintf(prn, " variable %d: '%s'\n", v+1, boxinfo.varname[realv]);
#ifdef notdef  
	    /* This is wrong!  How do you identify character data? */
	    if (line[51] != '2') {
		pprintf(prn, "   Non-numeric data: will be skipped\n");
		varstart[v] = 0;
		varsize[v] = 0;
		v++;
		break;
	    }
#endif
	    strncpy(tmp, line+52, 6);
	    tmp[6] = '\0';
	    varstart[v] = atoi(tmp) - 1;
	    pprintf(prn, "   starting col. %d, ", varstart[v]);
	    strncpy(tmp, line+58, 4);
	    tmp[4] = '\0';
	    varsize[v] = atoi(tmp);
	    pprintf(prn, "field width %d, ", varsize[v]);
	    strncpy(tmp, line+62, 2);
	    tmp[2] = '\0';
	    pprintf(prn, "decimal places %d\n", atoi(tmp));
	    tmp[0] = '\0';
	    strncpy(tmp, line+64, 20);
	    tmp[20] = '\0';
	    unspace(tmp);
	    if (strlen(tmp))
		pprintf(prn, "   Warning: coded variable (format '%s' "
			"in BOX file)\n", tmp);
	    strncpy(boxinfo.label[realv], line+87, 99);
	    boxinfo.label[realv][99] = '\0';
	    unspace(boxinfo.label[realv]);
	    pprintf(prn, "   definition: '%s'\n", boxinfo.label[realv]);
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
		x = strtod(tmp, &test);
#ifdef BOX_DEBUG
		fprintf(stderr, "read %d chars from pos %d: '%s' -> %g\n",
			varsize[i], varstart[i], tmp, x); 
#endif
		if (!strcmp(tmp, test)) {
		    pprintf(prn, "'%s' -- no numeric conversion performed!\n", 
			    tmp);
		    x = -999.0;
		}
		if (test[0] != '\0') {
		    if (isprint(test[0]))
			pprintf(prn, "Extraneous character '%c' in data\n", 
				test[0]);
		    else
			pprintf(prn, "Extraneous character (0x%x) in data\n", 
				test[0]);
		    x = -999.0;
		}
		if (errno == ERANGE) {
		    pprintf(prn, "'%s' -- number out of range!\n", tmp);
		    x = -999.0;
		}
		boxZ[boxinfo.n * realv + t] = x;
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

    pprintf(prn, "done reading data\n");
    fclose(fp);

    free(varstart);
    free(varsize);
    free(line);

    dataset_dates_defaults(&boxinfo);

    if (dumpvars) {
	dataset_drop_vars(dumpvars, &boxZ, &boxinfo);
	pprintf(prn, "Warning: discarded %d non-numeric variable(s)\n", 
		dumpvars);
    }

    if (*pZ == NULL) {
	*pZ = boxZ;
	*pdinfo = boxinfo;
    }

    return 0;
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
    size_t n = strlen(fname);
    int i, c, comma, ftype = GRETL_NATIVE_DATA;
    char teststr[5];
    FILE *fp;

    /* might be a script file? */
    if (n > 4 && strcmp(fname + n - 4, ".inp") == 0) 
	return GRETL_SCRIPT;

    addpath(fname, ppaths, 0);    

    fp = fopen(fname, "r");
    if (fp == NULL)  
	return GRETL_NATIVE_DATA; /* may be native file in different location */

    /* first look at extension */
    n = strlen(fname);
    if (n >= 5) {
	if (strcmp(fname + n - 4, ".csv") == 0) 
	    ftype = GRETL_CSV_DATA;
	else if (strcmp(fname + n - 4, ".box") == 0) 
	    ftype = GRETL_BOX_DATA;
    }
    /* then take a peek at content */
    comma = 0;
    for (i=0; i<80; i++) {
	c = getc(fp);
	if (c == EOF || c == '\n') break;
	if (c == ',') comma = 1;
	if (!isprint(c) && c != '\r') {
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
	if (comma) return GRETL_CSV_DATA;
	else {
	    pprintf(prn, "csv file seems to be malformed\n");
	    return GRETL_UNRECOGNIZED;
	}
    case GRETL_BOX_DATA: 
	if (strcmp(teststr, "00**") == 0) return GRETL_BOX_DATA;
	else {
	    pprintf(prn, "box file seems to be malformed\n");
	    return GRETL_UNRECOGNIZED;
	}
    }
    return GRETL_NATIVE_DATA; /* FIXME? */
}

