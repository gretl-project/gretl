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

#define MAXIGNORE 4
#define DESLINE 1024
#define MAXLEN 512

typedef struct {
    int ts;
    int pd;
    double sd0;
    int ignore[MAXIGNORE+1];
} DESINFO;

/* Reader for the datasets that accompany Wooldridge's "Introductory
   Econometrics: A Modern Approach" (South-Western)
*/

static int names_differ (const char *s1, const char *s2)
{
    char *p = strrchr(s1, '/');

    if (p != NULL && *(p + 1)) s1 = p + 1;

    while (*s1 && *s2) 
	if (*s1++ != tolower(*s2++)) return 1;
    if (*s1 || *s2) return 1;
    return 0;
}

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

static int read_des (const char *fname, double ***pZ,
		     DATAINFO **ppinfo)
{
    FILE *fp;
    char *p, line[DESLINE];
    int n = 0, v = 0;
    int blocknum = 0, err = 0, v2 = 0;
    int blankbak = 0, donealloc = 0;
    DATAINFO *dinfo;

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
	if (blocknum == 0) {
	    if (names_differ(fname, line)) {
		fprintf(stderr, "filename check failed: '%s' versus '%s'\n", 
			fname, line);
		err = 1;
	    }
	}
	else if (blocknum == 1) {
	    count_varnames(line, &v);
	}
	else if (blocknum == 2) {
	    if (sscanf(line, "%*s %d", &n) != 1) {
		fprintf(stderr, "Failed to read number of observations\n");
		err = 1;
	    }
	}
	else if (blocknum == 3) {
	    int i;
	    char varname[9];

	    if (!donealloc) {
		dinfo = create_new_dataset(pZ, v, n, 0);
		if (dinfo != NULL) {
		    *ppinfo = dinfo;
		    donealloc = 1;
		} else {
		    fprintf(stderr, "Out of memory\n");
		    err = 1;
		}
	    }

	    if (!err && sscanf(line, "%d. %9s", &i, varname) != 2) {
		fprintf(stderr, "Failed to read variables info block\n");
		err = 1;
	    } else {
		char *p = go_to_field(line, 3);

		printf("variable #%d: name='%s'", i, varname);
		if (*p) printf("comment='%s'\n", p);
		else printf("\n");
		v2++;
	    }
	}
	else break;
    }

    if (v2 != v) {
	fprintf(stderr, "Numbers of vars in block 1 and block 3 differ\n");
	err = 1;
    }

    if (!err)
	fprintf(stderr, "%s: got %d vars, %d observations\n", 
		fname, dinfo->v, dinfo->n);

    fclose(fp);
    return err;
}

static int read_raw (const char *fname, DATAINFO *pdinfo)
{
    FILE *fp;
    int i, t;
    char value[16];
    int err = 0;
    int literal = 0, litcol = 0;

    fp = fopen(fname, "r");
    if (fp == NULL) return 1;

    for (t=0; t<pdinfo->n; t++) {
	for (i=0; i<pdinfo->v; i++) {
	    if (fscanf(fp, "%15s", value) != 1) {
		fprintf(stderr, "scan fail on data[%d][%d]\n", i, t);
		err = 1;
	    } else {
		printf("got data[%d][%d] = '%s'\n", i, t, value);
		if (*value == '"') {
		    literal = 1;
		    if (i > litcol) litcol = i;
		}
		else if (strcmp(value, ".") == 0)
		    printf("got a missing value\n"); 
	    }
	    if (err) break;
	}
	if (err) break;
    }

    if (literal)
	printf("Data file contains non-numeric entries\n");

    if (litcol)
	fprintf(stderr, "%s: found literal data in col %d\n",
		fname, litcol + 1);

    fclose(fp);
    return err;
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
    size_t sz;

    if (strstr(ignore, "ignore=") == NULL) return;
    ignore += strlen("ignore=");
    while (ic < MAXIGNORE) {
	sz = strcspn(ignore, ",");
	if (sz == 0) break;
	*numstr = 0;
	strncat(numstr, ignore, sz);
	ig = atoi(numstr);
	ic++;
	des->ignore[ic] = ig;
	ignore += sz + 1;
    }
    
    des->ignore[0] = ic;
}


static int get_jw_structure (const char *fname, DESINFO *des, char *errbuf)
{
    FILE *fp;
    char line[80];
    char gotname[12], ts[32], ignore[32];
    char *p, *dname;
    int pd;
    double sd0;

    fp = fopen("jw_structure", "r");
    if (fp == NULL) {
	sprintf(errbuf, _("bad"));
	return 1;
    }

    desinfo_init(des);

    p = strrchr(fname, '/');
    if (p != NULL) dname = p + 1;
    else {
	p = strrchr(fname, '\\');
	if (p != NULL) dname = p + 1;
	else dname = fname;
    }

    while (fgets(line, 79, fp)) {
	if (*line == '#') continue;
	if (strstr(line, "difficult")) break;
	if (strstr(line, dname)) {
	    int f;

	    f = sscanf(line, "%s %s %d %lf %s", gotname, ts, &pd, &sd0, ignore);
	    if (f != 5) break;
	    if (!strcmp(ts, "TIME_SERIES")) des->ts = 1;
	    des->pd = pd;
	    des->sd0 = sd0; 
	    parse_ignore(des, ignore);
	}
    }

    return 0;
}

int des_get_data (const char *fname, double ***pZ, DATAINFO *pdinfo,
		  char *errbuf)
{
    int err = 0;
    double **newZ;
    DATAINFO *newinfo;
    DESINFO des;

    *errbuf = '\0';

    newinfo = datainfo_new();
    if (newinfo == NULL) {
	sprintf(errtext, _("Out of memory\n"));
	return 1;
    }

    err = get_jw_structure (fname, &des, errbuf);

    if (des.pd) {
	newinfo->pd = des.pd;
	newinfo->sd0 = des.sd0;
	/* newinfo->stobs?? */
	newinfo->time_series = des.ts;
    } 

    start_new_Z(&newZ, newinfo, 0);

    if (!des.ts) {
	strcpy(newinfo->stobs, "1");
	sprintf(newinfo->endobs, "%d", newinfo->n);
	newinfo->sd0 = 1.0;
	newinfo->pd = 1;
	newinfo->time_series = 0;
    } else {
	ntodate(newinfo->endobs, newinfo->n - 1, newinfo);
    }
    newinfo->extra = 0; 


    if (*pZ == NULL) { /* FIXME */
	*pZ = newZ;
	*pdinfo = *newinfo;
    } 

    return err;

}

int main (int argc, char *argv[])
{
    char desfile[MAXLEN], rawfile[MAXLEN];
    int err = 0;
    double **Z = NULL;
    DATAINFO *pdinfo = NULL;

    pdinfo = datainfo_new();
    if (pdinfo == NULL) {
	fprintf(stderr, "Out of memory\n");
	exit(EXIT_FAILURE);
    }

    strcpy(desfile, argv[1]);
    strcat(desfile, ".des");
    strcpy(rawfile, argv[1]);
    strcat(rawfile, ".raw");

    err = read_des(desfile, &Z, &pdinfo, markers);
    if (err) {
	sprintf(errbuf, "Problem reading %s\n", desfile);
	return 1;
    }

    err = read_raw(rawfile, pdinfo);
    if (err) {
	sprintf(errbuf, "Problem reading %s\n", rawfile);
	return 1;
    }

    return 0;
}


