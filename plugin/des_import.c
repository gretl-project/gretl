/*
 *  Copyright (c) by Allin Cottrell 2002
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

#include <gtk/gtk.h>

#ifdef G_OS_WIN32
# include "../winconfig.h"
#else
# include "../config.h"
#endif

#include "libgretl.h"

/* #define VERBOSE 1 */

#define MAXIGNORE 4
#define DESLINE 1024

typedef struct {
    int ts;
    int pd;
    double sd0;
    int ignore[MAXIGNORE+1];
} DESINFO;

/* Reader for the datasets that accompany Wooldridge's "Introductory
   Econometrics: A Modern Approach" (South-Western)
*/

static int blankline (const char *line)
{
    while (*line) 
	if (!isspace((unsigned char) *line++)) return 0;
    return 1;
}

static int count_varnames (const char *line, int *v)
{
    char varname[9];

    while (1) {
	line += strspn(line, " \t");
	if (*line == 0 || sscanf(line, "%8s", varname) != 1) break;
	*v += 1;
	line += strlen(varname);
    }

    return 0;
}

static char *go_to_field (char *s, int field)
{
    int i;

    while (isspace((unsigned char) *s)) s++;
    
    for (i=0; i<field-1; i++) {
	while (!isspace((unsigned char) *s)) s++;
	while (isspace((unsigned char) *s)) s++;
    }

    return s;
}

static int des_ignore (int varnum, DESINFO *des)
{
    int i;

    for (i=1; i<=des->ignore[0]; i++) {
	if (varnum == des->ignore[i]) return 1;
    }
    return 0;
}

static void desinfo_init (DESINFO *des)
{
    des->ts = 0;
    des->pd = 0;
    des->sd0 = 0.0;
    des->ignore[0] = 0;
}

static void parse_ignore (DESINFO *des, const char *ignore)
{
    int ig, ic = 0;
    char numstr[4];
    const char *p = ignore;
    size_t sz, len = strlen(ignore);

    if (strstr(p, "ignore=") == NULL) return;
    p += strlen("ignore=");
    while (ic < MAXIGNORE && p - ignore < len) {
	sz = strcspn(p, ",");
	if (sz == 0) break;
	*numstr = 0;
	strncat(numstr, p, sz);
	ig = atoi(numstr);
	ic++;
	des->ignore[ic] = ig;
	p += sz + 1;
    }
    
    des->ignore[0] = ic;
#ifdef VERBOSE
    printlist(des->ignore, "ignore list");
#endif
}

static void get_jwfile (const char *fname, char *jwfile)
{
    char *p;

    strcpy(jwfile, fname);
    p = strrchr(jwfile, '/');
    if (p == NULL) p = strrchr(jwfile, '\\');
    if (p != NULL) strcpy(p + 1, "jw_structure");
}

static void get_dname (const char *fname, char *dname)
{
    char *p;

    *dname = 0;
    p = strrchr(fname, '/');
    if (p != NULL) strcpy(dname, p + 1);
    else {
	p = strrchr(fname, '\\');
	if (p != NULL) strcpy(dname, p + 1);
    }
    p = strrchr(dname, '.');
    if (p) *p = 0;
}

static int read_jw_structure (const char *fname, DESINFO *des, char *errbuf)
{
    FILE *fp;
    char line[80];
    char gotname[12], ts[32], ignore[32];
    char jwfile[MAXLEN], dname[32];
    int pd;
    double sd0;

    get_jwfile(fname, jwfile);

    fp = fopen(jwfile, "r");
    if (fp == NULL) {
	sprintf(errbuf, _("Couldn't open %s"), jwfile);
	return 1;
    }

    desinfo_init(des);

    get_dname(fname, dname);
    if (*dname == 0) return 0;

    while (fgets(line, 79, fp)) {
	if (*line == '#') continue;
	if (strstr(line, "difficult")) break;
	if (strstr(line, dname)) {
	    int f;
	    char *p;

	    if ((p = strchr(line, '\n'))) *p = 0;
	    if ((p = strchr(line, '\r'))) *p = 0;
#ifdef VERBOSE
	    fprintf(stderr, "%s\ngot line: '%s'\n", jwfile, line);
#endif
	    f = sscanf(line, "%s %s %d %lf %s", gotname, ts, &pd, &sd0, ignore);
	    if (f != 5) break;
	    if (!strcmp(ts, "TIME_SERIES")) des->ts = 1;
	    des->pd = pd;
	    des->sd0 = sd0; 
	    parse_ignore(des, ignore);
	}
    }

    fclose(fp);

    return 0;
}

static int read_des (const char *fname, DESINFO *des, double ***pZ,
		     DATAINFO **ppinfo, char *errbuf)
{
    FILE *fp;
    char *p, line[DESLINE];
    int n = 0, v = 0;
    int blocknum = 0, err = 0, v2 = 0;
    int blankbak = 0, donealloc = 0;
    DATAINFO *dinfo = NULL;

    fp = fopen(fname, "r");
    if (fp == NULL) return 1;

    /* order of events: (1) filename (upper case);
       (2) listing of vars; (3) Obs line; (4) Numbered listing
       of varnames plus short descriptions.  Blocks separated
       by vertical whitespace (usually just one blank line)
    */

    while (!err && fgets(line, DESLINE-1, fp)) {
	if (blankline(line)) {
	    /* allow for extra blank lines,
	       found in a few of the files */
	    if (!blankbak) blocknum++;
	    blankbak = 1;
	    continue;
	} else
	    blankbak = 0;
	if ((p = strchr(line, '\n'))) *p = 0;
	if ((p = strchr(line, '\r'))) *p = 0;
        if (blocknum == 0) continue;
	else if (blocknum == 1) {
	    count_varnames(line, &v);
	}
	else if (blocknum == 2) {
	    if (sscanf(line, "%*s %d", &n) != 1) {
		sprintf(errbuf, "Failed to read number of observations");
		err = 1;
	    }
	}
	else if (blocknum == 3) {
	    int i;
	    char varname[9];

	    v -= des->ignore[0];
	    v++; /* allow for constant */

	    if (!donealloc) {
		dinfo = create_new_dataset(pZ, v, n, 0);
		if (dinfo != NULL) {
		    *ppinfo = dinfo;
		    donealloc = 1;
		} else {
		    sprintf(errbuf, "Out of memory");
		    err = 1;
		}
	    }

	    if (!err && sscanf(line, "%d. %9s", &i, varname) != 2) {
		sprintf(errbuf, "Failed to read variables info block");
		err = 1;
	    } else {
		char *p = go_to_field(line, 3);

#ifdef VERBOSE
		printf("got variable #%d: name='%s'", i, varname);
		if (*p) printf("comment='%s'\n", p);
		else printf("\n");
#endif
		if (!des_ignore(i, des)) {
		    v2++;
		    strcpy(dinfo->varname[v2], varname);
		    strcpy(dinfo->label[v2], p);
		} else 
		    printf("ignoring variable '%s'\n", varname);
	    }
	}
	else break;
    }

    if (!err)
	fprintf(stderr, "Read %s:\ngot %d real variables, %d observations\n", 
	       fname, dinfo->v, dinfo->n);

    fclose(fp);
    return err;
}

static int read_raw (const char *fname, DESINFO *des,
		     DATAINFO *pdinfo, double **Z,
		     char *errbuf)
{
    FILE *fp;
    int i, t, des_i;
    char value[16];
    int err = 0;
    int literal = 0, litcol = 0;

    fp = fopen(fname, "r");
    if (fp == NULL) return 1;

    for (t=0; t<pdinfo->n; t++) {
	des_i = 1;
	for (i=1; i<pdinfo->v; i++) {
	    while (1) {
		if (des_ignore(des_i, des)) {
		    fscanf(fp, "%s", value);
#ifdef VVERBOSE
		    fprintf(stderr, "ignoring data(%d,%d) '%s'\n", 
			    des_i, t, value);
#endif	
		    des_i++;
		} else break;
	    }
	    if (fscanf(fp, "%15s", value) != 1) {
		sprintf(errbuf, "scan fail on data[%d][%d]", i, t);
		err = 1;
	    } else {
#ifdef VVERBOSE
		fprintf(stderr, "got data[%d][%d] = '%s'\n", i, t, value);
#endif
		if (*value == '"') {
		    literal = 1;
		    if (i > litcol) litcol = i;
		}
		else {
		    if (strcmp(value, ".") == 0) {
#ifdef VVERBOSE
			fprintf(stderr, "got a missing value\n"); 
#endif
			Z[i][t] = NADBL;
		    } else {
			Z[i][t] = atof(value);
		    }
		}
	    }
	    if (err) break;
	    des_i++;
	}
	if (err) break;
    }

    if (literal)
	sprintf(errbuf, "Data file contains non-numeric entries");

    if (litcol)
	sprintf(errbuf, "%s: found literal data in col %d",
		fname, litcol + 1);

    fclose(fp);
    return err;
}


int des_get_data (const char *fname, double ***pZ, DATAINFO *pdinfo,
		  char *errbuf)
{
    int err = 0;
    double **newZ = NULL;
    DATAINFO *newinfo = NULL;
    DESINFO des;
    char rawfile[MAXLEN];

    *errbuf = '\0';

    /* get any time-series info, cols to be ignored */
    err = read_jw_structure (fname, &des, errbuf);
    if (err) return err;

    /* read des file for number of obs, varnames and descriptions */
    err = read_des (fname, &des, &newZ, &newinfo, errbuf);

    /* modify datainfo based on time-series info gathered */
    if (des.pd) {
	newinfo->pd = des.pd;
	newinfo->sd0 = des.sd0;
	ntodate(newinfo->stobs, newinfo->t1, newinfo);
	ntodate(newinfo->endobs, newinfo->t2, newinfo);
	newinfo->time_series = des.ts;
    }

    switch_ext(rawfile, fname, "raw");
    err = read_raw(rawfile, &des, newinfo, newZ, errbuf);
    if (err) {
	sprintf(errbuf, "Problem reading %s", rawfile);
	return 1;
    }    

    if (*pZ == NULL) { /* FIXME */
	*pZ = newZ;
	*pdinfo = *newinfo;
    } 

    return err;

}



