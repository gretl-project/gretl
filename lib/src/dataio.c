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
#include "gretl_string_table.h"
#include "dbwrite.h"
#include "libset.h"
#include "gretl_xml.h"
#include "gretl_panel.h"
#include "csvdata.h"
#include "usermat.h"

#include <ctype.h>
#include <time.h>
#include <errno.h>

#include <glib.h>

typedef enum {
    GRETL_FMT_FLOAT = 1, /* single-precision binary data */
    GRETL_FMT_DOUBLE,    /* double-precision binary data */
    GRETL_FMT_OCTAVE,    /* data in Gnu Octave format */
    GRETL_FMT_CSV,       /* data in Comma Separated Values format */
    GRETL_FMT_R,         /* data in Gnu R format */
    GRETL_FMT_GZIPPED,   /* gzipped data */
    GRETL_FMT_TRAD,      /* traditional (ESL-style) data */
    GRETL_FMT_DAT,       /* data in PcGive format */
    GRETL_FMT_DB,        /* gretl native database format */
    GRETL_FMT_JM         /* JMulti ascii data */
} GretlDataFormat;

#define IS_DATE_SEP(c) (c == '.' || c == ':' || c == ',')

static int writelbl (const char *lblfile, const int *list, 
		     const DATAINFO *pdinfo);
static int writehdr (const char *hdrfile, const int *list, 
		     const DATAINFO *pdinfo, int opt);

static char STARTCOMMENT[3] = "(*";
static char ENDCOMMENT[3] = "*)";

#define PROGRESS_BAR "progress_bar"

double get_date_x (int pd, const char *obs)
{
    double x = 1.0;

    if ((pd == 5 || pd == 6 || pd == 7 || pd == 52) && strlen(obs) > 4) { 
	/* calendar data */
	long ed = get_epoch_day(obs);

	if (ed >= 0) {
	    x = ed;
	}
    } else {
	x = obs_str_to_double(obs); 
    }

    return x;
}

/* Skip past comments in .hdr file.  Return 0 if comments found,
   otherwise 1.
*/

static int skipcomments (FILE *fp, const char *str)
{
    char word[MAXLEN];  /* should be big enough to accommodate
			   strings among the comments? */

    *word = '\0';

    if (strncmp(str, STARTCOMMENT, 2) == 0) {
        while (strcmp(word, ENDCOMMENT)) {
            fscanf(fp, "%s", word);
        }
        return 0;
    } 

    return 1;
}

static int comment_lines (FILE *fp, char **pbuf)
{
    char s[MAXLEN], *mybuf = NULL;
    int count = 0, bigger = 1;

    if (fgets(s, sizeof s, fp) == NULL) {
	return 0;
    }

    if (!strncmp(s, STARTCOMMENT, 2)) {
	*pbuf = malloc(20 * MAXLEN);

	if (*pbuf == NULL) {
	    return -1;
	}

	**pbuf = '\0';

	while (fgets(s, sizeof s, fp)) {
	    if (!strncmp(s, ENDCOMMENT, 2)) {
		break;
	    }
	    if (++count > 20 * bigger) {
		size_t bufsize = 20 * MAXLEN * ++bigger;

		mybuf = realloc(*pbuf, bufsize);
		if (mybuf == NULL) {
		    return -1;
		} else {
		    *pbuf = mybuf;
		}
	    }
	    strcat(*pbuf, s);
	} 
    }

    return count;
}

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

static int readdata (FILE *fp, const DATAINFO *pdinfo, double **Z,
		     int binary, int old_byvar)
{
    int i, t, n = pdinfo->n;
    char c, marker[OBSLEN];
    int err = 0;

    gretl_error_clear();

    if (binary == 1) { 
	/* single-precision binary data */
	float x;

	for (i=1; i<pdinfo->v; i++) {
	    for (t=0; t<n; t++) {
		if (!fread(&x, sizeof x, 1, fp)) {
		    sprintf(gretl_errmsg, _("WARNING: binary data read error at "
			    "var %d"), i);
		    return 1;
		}
		if (x == -999.0) {
		    Z[i][t] = NADBL;
		} else {
		    Z[i][t] = (double) x;
		}
	    }
	}
    } else if (binary == 2) { 
	/* double-precision binary data */
	double x;

	for (i=1; i<pdinfo->v; i++) {
	    for (t=0; t<n; t++) {
		if (!fread(&x, sizeof x, 1, fp)) {
		    sprintf(gretl_errmsg, 
			    _("WARNING: binary data read error at var %d"), i);
		    return 1;
		}
		if (x == -999.0) {
		    Z[i][t] = NADBL;
		} else {
		    Z[i][t] = x;
		}
	    }
	}
    } else if (old_byvar) {
	/* ascii data by variable */
	for (i=1; i<pdinfo->v; i++) {
	   for (t=0; t<n && !err; t++) {
		if ((fscanf(fp, "%lf", &Z[i][t])) != 1) {
		    sprintf(gretl_errmsg, 
			    _("WARNING: ascii data read error at var %d, "
			    "obs %d"), i, t + 1);
		    err = 1;
		    break;
		}
		if (Z[i][t] == -999.0) {
		    Z[i][t] = NADBL;
		} 
	   }
	}	       
    } else { 
	/* ascii data by observation */
	char sformat[8];

	sprintf(sformat, "%%%ds", OBSLEN - 1);

	gretl_push_c_numeric_locale();

	for (t=0; t<n && !err; t++) {
	    eatspace(fp);
	    c = fgetc(fp);  /* test for a #-opened comment line */
	    if (c == '#') {
		while (c != '\n') {
		    c = fgetc(fp);
		}
	    } else {
		ungetc(c, fp);
	    }
	    if (pdinfo->markers) {
		*marker = '\0';
		fscanf(fp, sformat, marker);
		if (*marker == '"' || *marker == '\'') {
		    strcpy(pdinfo->S[t], marker + 1);
		} else {
		    strcpy(pdinfo->S[t], marker);
		}
	    }
	    for (i=1; i<pdinfo->v; i++) {
		if ((fscanf(fp, "%lf", &Z[i][t])) != 1) {
		    sprintf(gretl_errmsg, 
			    _("WARNING: ascii data read error at var %d, "
			    "obs %d"), i, t + 1);
		    err = 1;
		    break;
		}
		if (Z[i][t] == -999.0) {
		    Z[i][t] = NADBL;
		} 
	    }
	}

	gretl_pop_c_numeric_locale();
    }

    return err;
}

static int gz_readdata (gzFile fz, const DATAINFO *pdinfo, double **Z,
			int binary)
{
    int i, t, n = pdinfo->n;
    int err = 0;
    
    gretl_error_clear();

    if (binary == 1) { 
	/* single-precision binary data */
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
    } else if (binary == 2) { 
	/* double-precision binary data */
	for (i=1; i<pdinfo->v; i++) {
	    if (!gzread(fz, &Z[i][0], n * sizeof(double))) {
		sprintf(gretl_errmsg, 
			_("WARNING: binary data read error at var %d"), i);
		return 1;
	    }
	}
    } else { 
	/* ascii data */
	char *line, numstr[24], sformat[8];
	int llen = pdinfo->v * 32;
	size_t offset;

	line = malloc(llen);
	if (line == NULL) {
	    return E_ALLOC;
	}

	sprintf(sformat, "%%%ds", OBSLEN - 1);

	gretl_push_c_numeric_locale();

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
		if (sscanf(line, sformat, pdinfo->S[t]) != 1) {
		   sprintf(gretl_errmsg, 
			   _("WARNING: failed to read case marker for "
			   "obs %d"), t + 1);
		   err = 1;
		   break;
		}
		pdinfo->S[t][OBSLEN-1] = 0;
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
		if (i < pdinfo->v - 1) {
		    offset += strlen(numstr) + 1;
		}
	    }

	    if (err) break;
	}

	free(line);

	gretl_pop_c_numeric_locale();
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
    int testchar = 'a';
    int ret = 0;

    gretl_error_clear();

    if (gretl_reserved_word(varname)) {
	ret = VARNAME_RESERVED;
    } else if (!(isalpha((unsigned char) *varname))) {
	testchar = *varname;
        ret = VARNAME_FIRSTCHAR;
    } else {
	const char *p = varname;

	while (*p && testchar == 'a') {
	    if (!(isalpha((unsigned char) *p))  
		&& !(isdigit((unsigned char) *p))
		&& *p != '_') {
		testchar = *p;
		ret = VARNAME_BADCHAR;
	    }
	    p++;
	}
    }

    if (testchar != 'a') {
	if (isprint((unsigned char) testchar)) {
	    if (ret == VARNAME_FIRSTCHAR) {
		sprintf(gretl_errmsg, _("First char of varname ('%c') is bad\n"
					"(first must be alphabetical)"), 
			(unsigned char) testchar);
	    } else {
		sprintf(gretl_errmsg, _("Varname contains illegal character '%c'\n"
					"Use only letters, digits and underscore"), 
			(unsigned char) testchar);
	    }
	} else {
	    if (ret == VARNAME_FIRSTCHAR) {
		sprintf(gretl_errmsg, _("First char of varname (0x%x) is bad\n"
					"(first must be alphabetical)"), 
			(unsigned) testchar);
	    } else {
		sprintf(gretl_errmsg, _("Varname contains illegal character 0x%x\n"
					"Use only letters, digits and underscore"), 
			(unsigned) testchar);
	    }
	}
    }

    return ret;
}   

static int readhdr (const char *hdrfile, DATAINFO *pdinfo, 
		    int *binary, int *old_byvar)
{
    FILE *fp;
    int n, i = 0, panel = 0, descrip = 0;
    char str[MAXLEN], byobs[6], option[8];

    gretl_error_clear();

    fp = gretl_fopen(hdrfile, "r");
    if (fp == NULL) {
	sprintf(gretl_errmsg, _("Couldn't open file %s"),  hdrfile);
	return E_FOPEN;
    }

    fscanf(fp, "%s", str);
    i += skipcomments(fp, str); 

    /* find number of variables */

    while (1) {
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

    if (dataset_allocate_varnames(pdinfo)) {
	return E_ALLOC;
    }

    i = 1;
    fp = gretl_fopen(hdrfile, "r");

    str[0] = 0;
    fscanf(fp, "%s", str);
    if (skipcomments(fp, str)) {
        safecpy(pdinfo->varname[i], str, VNAMELEN - 1);
	if (check_varname(pdinfo->varname[i++])) {
	    goto varname_error;
	}
    } else {
	descrip = 1; /* comments were found */
    }

    while (1) {
        fscanf(fp, "%s", str);
	n = strlen(str);
	if (str[n-1] != ';') {
            safecpy(pdinfo->varname[i], str, VNAMELEN - 1);
	    if (check_varname(pdinfo->varname[i++])) {
		goto varname_error;
	    }
        } else {
	    if (n > 1) {
		safecpy(pdinfo->varname[i], str, n-1);
		pdinfo->varname[i][n] = '\0';
		if (check_varname(pdinfo->varname[i])) {
		    goto varname_error; 
		}
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

    if (pdinfo->sd0 >= 2.0) {
        pdinfo->structure = TIME_SERIES; /* actual time series? */
    } else if (pdinfo->sd0 > 1.0) {
	pdinfo->structure = STACKED_TIME_SERIES; /* panel data? */
    } else {
	pdinfo->structure = CROSS_SECTION;
    }

    pdinfo->n = -1;
    pdinfo->n = dateton(pdinfo->endobs, pdinfo) + 1;

    *binary = 0;
    pdinfo->markers = NO_MARKERS;

    n = fscanf(fp, "%5s %7s", byobs, option);

    if (n == 1 && strcmp(byobs, "BYVAR") == 0) {
	*old_byvar = 1;
    } else if (n == 2) {
	if (strcmp(option, "SINGLE") == 0) {
	    *binary = 1;
	} else if (strcmp(option, "BINARY") == 0) {
	    *binary = 2;
	} else if (strcmp(option, "MARKERS") == 0) {
	    pdinfo->markers = 1;
	} else if (strcmp(option, "PANEL2") == 0) {
	    panel = 1;
	    pdinfo->structure = STACKED_TIME_SERIES;
	} else if (strcmp(option, "PANEL3") == 0) {
	    panel = 1;
	    pdinfo->structure = STACKED_CROSS_SECTION;
	}
    } 

    if (!panel && fscanf(fp, "%6s", option) == 1) {
	if (strcmp(option, "PANEL2") == 0) {
	    pdinfo->structure = STACKED_TIME_SERIES;
	} else if (strcmp(option, "PANEL3") == 0) {
	    pdinfo->structure = STACKED_CROSS_SECTION;
	}
    }

    if (fp != NULL) {
	fclose(fp);
    }

    /* last pass, to pick up data description */
    pdinfo->descrip = NULL;
    if (descrip) {
	char *dbuf = NULL;
	int lines;

	fp = gretl_fopen(hdrfile, "r");
	if (fp == NULL) return 0;
	if ((lines = comment_lines(fp, &dbuf)) > 0) {
	    delchar('\r', dbuf);
	    pdinfo->descrip = malloc(strlen(dbuf) + 1);
	    if (pdinfo->descrip != NULL) {
		strcpy(pdinfo->descrip, dbuf);
	    }
	    free(dbuf);
	} else if (lines < 0) {
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

static int bad_date_string (const char *s)
{
    int err = 0;

    gretl_error_clear();

    while (*s && !err) {
	if (!isdigit((unsigned char) *s) && !IS_DATE_SEP(*s)) {
	    if (isprint((unsigned char) *s)) {
		sprintf(gretl_errmsg, 
			_("Bad character '%c' in date string"), *s);
	    } else {
		sprintf(gretl_errmsg, 
			_("Bad character %d in date string"), *s);
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

#define DATES_DEBUG 0

static int match_obs_marker (const char *s, const DATAINFO *pdinfo)
{
    char test[OBSLEN];
    int t;

#if DATES_DEBUG
    fprintf(stderr, "dateton: checking marker strings\n");
#endif

    maybe_unquote_label(test, s);

    for (t=0; t<pdinfo->n; t++) {
	if (!strcmp(test, pdinfo->S[t])) {
	    /* handled */
	    return t;
	}
    }

    if (isalpha(*s)) {
	/* try harder */
	int k = strlen(test);

	for (t=0; t<pdinfo->n; t++) {
	    if (!strncmp(test, pdinfo->S[t], k)) {
		return t;
	    }
	}
    }

    return -1;
}

static int 
real_dateton (const char *date, const DATAINFO *pdinfo, int nolimit)
{
    int handled = 0;
    int t, n = -1;

    /* first check if this is calendar data and if so,
       treat accordingly */

    if (calendar_data(pdinfo)) {
#if DATES_DEBUG
	fprintf(stderr, "dateton: treating as calendar data\n");
#endif
	if (pdinfo->markers && pdinfo->S != NULL) {
	    /* "hard-wired" calendar dates as strings */
	    for (t=0; t<pdinfo->n; t++) {
		if (!strcmp(date, pdinfo->S[t])) {
		    /* handled */
		    return t;
		}
	    }
	    /* try allowing for 2- versus 4-digit years? */
	    if (strlen(pdinfo->S[0]) == 10 &&
		(!strncmp(pdinfo->S[0], "19", 2) || 
		 !strncmp(pdinfo->S[0], "20", 2))) {
		for (t=0; t<pdinfo->n; t++) {
		    if (!strcmp(date, pdinfo->S[t] + 2)) {
			/* handled */
			return t;
		    }
		}		
	    }
	    /* out of options: abort */
	    return -1;
	} else {
	    /* automatic calendar dates */
	    n = calendar_obs_number(date, pdinfo);
	    handled = 1;
	} 
    } else if (dataset_is_daily(pdinfo) ||
	       dataset_is_weekly(pdinfo)) {
#if DATES_DEBUG
	fprintf(stderr, "dateton: trying undated time series\n");
#endif
	t = positive_int_from_string(date);
	if (t > 0) {
	    n = t - 1;
	    handled = 1;
	}
    } else if (dataset_is_decennial(pdinfo)) {
	t = positive_int_from_string(date);
	if (t > 0) {
	    n = (t - pdinfo->sd0) / 10;
	    handled = 1;
	}	
    } else if (pdinfo->markers && pdinfo->S != NULL) {
	t = match_obs_marker(date, pdinfo);
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
	int dotpos1, dotpos2;

#if DATES_DEBUG
	fprintf(stderr, "dateton: treating as regular numeric obs\n");
#endif
	if (bad_date_string(date)) {
	    return -1;
	}

	dotpos1 = get_dot_pos(date);
	dotpos2 = get_dot_pos(pdinfo->stobs);

	if ((dotpos1 && !dotpos2) || (dotpos2 && !dotpos1)) {
	    sprintf(gretl_errmsg, _("Date strings inconsistent"));
	} else if (!dotpos1 && !dotpos2) {
	    n = atoi(date) - atoi(pdinfo->stobs);
	} else {
	    char majstr[5] = {0};
	    char minstr[3] = {0};
	    char majstr0[5] = {0};
	    char minstr0[3] = {0};

	    int maj, min;
	    int maj0, min0;

	    strncat(majstr, date, dotpos1);
	    maj = atoi(majstr);
	    strncat(minstr, date + dotpos1 + 1, 2);
	    min = atoi(minstr);	    

	    strncat(majstr0, pdinfo->stobs, dotpos2);
	    maj0 = atoi(majstr0);
	    strncat(minstr0, pdinfo->stobs + dotpos2 + 1, 2);
	    min0 = atoi(minstr0);
    
	    n = pdinfo->pd * (maj - maj0) + (min - min0);
	}
    }

    if (!nolimit && pdinfo->n > 0 && n >= pdinfo->n) {
	fprintf(stderr, "n = %d, pdinfo->n = %d: out of bounds\n", n, pdinfo->n);
	gretl_errmsg_set(_("Observation number out of bounds"));
	n = -1; 
    }

    return n;
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
    return real_dateton(date, pdinfo, 0);
}

/* special for appending data: allow the date to be outside of
   the range of the current dataset */

static int merge_dateton (const char *date, const DATAINFO *pdinfo)
{
    return real_dateton(date, pdinfo, 1);
}

static char *
out_of_range_panel_obs (char *s, int t, const DATAINFO *pdinfo)
{
    int i = t / pdinfo->pd + 1;
    int j = (t + 1) % pdinfo->pd;

    if (j == 0) {
	j = pdinfo->pd;
    }

    sprintf(s, "%d:%0*d", i, pdinfo->paninfo->olen, j);
    return s;
}

static char *
real_ntodate (char *datestr, int t, const DATAINFO *pdinfo, int full)
{
    double x;

#if 0
    fprintf(stderr, "real_ntodate: t=%d, pd=%d, sd0=%g\n",
	    t, pdinfo->pd, pdinfo->sd0);
#endif

    if (calendar_data(pdinfo)) {
	/* handles both daily and dated weekly data */
	if (pdinfo->markers && pdinfo->S != NULL) {
	    strcpy(datestr, pdinfo->S[t]);
	} else {
	    calendar_date_string(datestr, t, pdinfo);
	}
	if (!full && strlen(datestr) > 8) {
	    char tmp[12];

	    strcpy(tmp, datestr);
	    strcpy(datestr, tmp + 2);
	}
	return datestr;
    } else if (dataset_is_daily(pdinfo) || 
	       dataset_is_weekly(pdinfo)) {
	/* undated time series */
	x = date(t, 1, pdinfo->sd0);
	sprintf(datestr, "%d", (int) x);
	return datestr;
    } else if (dataset_is_decennial(pdinfo)) {
	x = pdinfo->sd0 + 10 * t;
	sprintf(datestr, "%d", (int) x);
	return datestr;
    } else if (pdinfo->paninfo != NULL) {
	/* indexed panel data */
	if (t < 0 || t >= pdinfo->n) {
	    return out_of_range_panel_obs(datestr, t, pdinfo);
	}
	sprintf(datestr, "%d:%0*d", pdinfo->paninfo->unit[t] + 1,
		pdinfo->paninfo->olen,
		pdinfo->paninfo->period[t] + 1);
	return datestr;
    }

    x = date(t, pdinfo->pd, pdinfo->sd0);

    if (pdinfo->pd == 1) {
        sprintf(datestr, "%d", (int) x);
    } else {
	int pdp = pdinfo->pd, len = 1;
	char fmt[8];

	while ((pdp = pdp / 10)) len++;
	sprintf(fmt, "%%.%df", len);
	sprintf(datestr, fmt, x);
	colonize_obs(datestr);
    }
    
    return datestr;
}

/**
 * ntodate:
 * @datestr: string to which date is to be printed.
 * @t: an observation number (zero-based).
 * @pdinfo: data information struct.
 * 
 * print to @datestr the calendar representation of observation
 * number @t.
 * 
 * Returns: the observation string.
 */

char *ntodate (char *datestr, int t, const DATAINFO *pdinfo)
{
    return real_ntodate(datestr, t, pdinfo, 0);
}

char *ntodate_full (char *datestr, int t, const DATAINFO *pdinfo)
{
    return real_ntodate(datestr, t, pdinfo, 1);
}

#define xround(x) (((x-floor(x))>.5)? ceil(x) : floor(x))

/* for "seasonal" time series data (broad sense): given
   the 0-based observation number, t, determine the
   sub-period at that obs. The "sub-period" might
   be the quarter, month, hour or whatever.  The value
   returned is zero-based (e.g. first quarter = 0).
*/

int get_subperiod (int t, const DATAINFO *pdinfo, int *err)
{
    int ret = 0;

    if (!dataset_is_seasonal(pdinfo)) {
	if (err != NULL) {
	    *err = E_PDWRONG;
	}
	return 0;
    }

    if (dataset_is_weekly(pdinfo)) {
	/* bodge -- what else to do? */
	ret = t % pdinfo->pd;
    } else if (calendar_data(pdinfo)) {
	/* dated daily data */
	char datestr[12];

	calendar_date_string(datestr, t, pdinfo);
	ret = get_day_of_week(datestr); 
    } else if (dataset_is_daily(pdinfo)) {
	/* bodge, again */
	ret = t % pdinfo->pd;
    } else {
	/* quarterly, monthly, hourly... */
	double x = date(t, pdinfo->pd, pdinfo->sd0);
	int i, d = ceil(log10(pdinfo->pd));

	x -= floor(x);
	for (i=0; i<d; i++) {
	    x *= 10;
	}
	ret = xround(x) - 1;
    }
    
    return ret;    
}

/* .......................................................... */

static int blank_check (FILE *fp)
{
    int i, deflt = 1;
    char s[MAXLEN];

    for (i=0; i<3 && deflt && fgets(s, MAXLEN-1, fp); i++) {
	if (i == 0 && strncmp(s, "(*", 2)) {
	    deflt = 0;
	} else if (i == 1 && strncmp(s, _("space for comments"), 18)) {
	    deflt = 0;
	} else if (i == 2 && strncmp(s, "*)", 2)) {
	    deflt = 0;
	}
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

    if ((hdr = gretl_fopen(hdrfile, "r")) == NULL) {
	pprintf(prn, _("Couldn't open %s\n"), hdrfile); 
	return 1;
    }

    /* see if it's just the default "space for comments" */
    if (blank_check(hdr)) { /* yes */
	pprintf(prn, _("No info in %s\n"), hdrfile);
	return 2;
    } 

    /* no, so restart the read */
    if ((hdr = gretl_fopen(hdrfile, "r")) == NULL) {
	pprintf(prn, _("Couldn't open %s\n"), hdrfile); 
	return 1;
    }    

    pprintf(prn, _("Data info in file %s:\n\n"), hdrfile);

    if (fgets(s, MAXLEN-1, hdr) != NULL && !strncmp(s, STARTCOMMENT, 2)) {
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

    if (i == 0) {
	pputs(prn, _(" (none)\n"));
    }

    pputc(prn, '\n');

    if (hdr != NULL) {
	fclose(hdr);
    }

    return 0;
}

static int writehdr (const char *hdrfile, const int *list, 
		     const DATAINFO *pdinfo, int opt)
{
    FILE *fp;
    char startdate[OBSLEN], enddate[OBSLEN];
    int i, binary = 0;

    if (opt == GRETL_FMT_FLOAT) {
	binary = 1;
    } else if (opt == GRETL_FMT_DOUBLE) {
	binary = 2;
    }

    ntodate_full(startdate, pdinfo->t1, pdinfo);
    ntodate_full(enddate, pdinfo->t2, pdinfo);

    fp = gretl_fopen(hdrfile, "w");
    if (fp == NULL) {
	return 1;
    }

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
	if (list[i] == 0) {
	    continue;
	}
	fprintf(fp, "%s ", pdinfo->varname[list[i]]);
	if (i && i <list[0] && (i+1) % 8 == 0) {
	    fputc('\n', fp);
	}
    }  
  
    fputs(";\n", fp);

    /* then obs line */
    fprintf(fp, "%d %s %s\n", pdinfo->pd, startdate, enddate);
    
    /* and flags as required */
    if (binary == 1) {
	fputs("BYVAR\nSINGLE\n", fp);
    } else if (binary == 2) {
	fputs("BYVAR\nBINARY\n", fp);
    } else { 
	fputs("BYOBS\n", fp);
	if (pdinfo->markers) {
	    fputs("MARKERS\n", fp);
	}
    }

    if (pdinfo->structure == STACKED_TIME_SERIES) {
	fputs("PANEL2\n", fp);
    } else if (pdinfo->structure == STACKED_CROSS_SECTION) {
	fputs("PANEL3\n", fp);
    }
    
    fclose(fp);

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
 */

int get_precision (const double *x, int n, int placemax)
{
    int t, p, pmax = 0;
    char *s, numstr[48];
    int n_ok = 0;
    double z;

    for (t=0; t<n; t++) {
	if (na(x[t])) {
	    continue;
	}

	n_ok++;

	z = fabs(x[t]);

	/* escape clause: numbers are too big or too small for
	   this treatment */
	if (z > 0 && (z < 1.0e-6 || z > 1.0e+8)) {
	    return PMAX_NOT_AVAILABLE;
	}

	p = placemax;
	sprintf(numstr, "%.*f", p, z);
	s = numstr + strlen(numstr) - 1;
	while (*s-- == '0') {
	    p--;
	}
	if (p > pmax) {
	    pmax = p;
	}
    }

    if (n_ok == 0) {
	pmax = PMAX_NOT_AVAILABLE;
    }

    return pmax;
}

gretlopt data_save_opt_from_suffix (const char *fname)
{
    gretlopt opt = OPT_NONE;

    if (has_suffix(fname, ".R")) {
	opt = OPT_R;
    } else if (has_suffix(fname, ".m")) {
	opt = OPT_M;
    } else if (has_suffix(fname, ".csv") ||
	       has_suffix(fname, ".txt") ||
	       has_suffix(fname, ".asc")) {
	opt = OPT_C;
    } 

    return opt;
}

static GretlDataFormat 
format_from_opt_or_name (gretlopt opt, const char *fname,
			 char *delim)
{
    GretlDataFormat fmt = 0;
    
    if (opt & OPT_T) {
	fmt = GRETL_FMT_TRAD;
    } else if (opt & OPT_M) {
	fmt = GRETL_FMT_OCTAVE;
    } else if (opt & OPT_R) {
	fmt = GRETL_FMT_R;
    } else if (opt & OPT_C) {
	fmt = GRETL_FMT_CSV;
    } else if (opt & OPT_Z) {
	fmt = GRETL_FMT_GZIPPED;
    } else if (opt & OPT_D) {
	fmt = GRETL_FMT_DB;
    } else if (opt & OPT_G) {
	fmt = GRETL_FMT_DAT;
    } else if (opt & OPT_J) {
	fmt = GRETL_FMT_JM;
    }

    if (fmt == 0) {
	if (has_suffix(fname, ".R")) {
	    fmt = GRETL_FMT_R;
	} else if (has_suffix(fname, ".csv")) {
	    fmt = GRETL_FMT_CSV;
	} else if (has_suffix(fname, ".m")) {
	    fmt = GRETL_FMT_OCTAVE;
	} else if (has_suffix(fname, ".txt") ||
		   has_suffix(fname, ".asc")) {
	    fmt = GRETL_FMT_CSV;
	    *delim = ' ';
	} 
    }

    return fmt;
}

static void date_maj_min (int t, const DATAINFO *pdinfo, int *maj, int *min)
{
    char obs[OBSLEN];
    char *s;

    ntodate(obs, t, pdinfo);

    *maj = atoi(obs);
    s = strchr(obs, ':');
    if (s != NULL && strlen(s) > 1) {
	*min = atoi(s + 1);
    } else {
	*min = 1;
    }
}

#define annual_data(p) (p->structure == TIME_SERIES && p->pd == 1)

/**
 * write_data:
 * @fname: name of file to write.
 * @list: list of variables to write (or %NULL to write all series).
 * @Z: data matrix.
 * @pdinfo: data information struct.
 * @opt: option flag indicating format in which to write the data.
 * @ppaths: pointer to paths information (should be NULL when not
 * called from gui).
 * 
 * Write out a data file containing the values of the given set
 * of variables.
 * 
 * Returns: 0 on successful completion, non-zero on error.
 */

int write_data (const char *fname, int *list, 
		const double **Z, const DATAINFO *pdinfo, 
		gretlopt opt, PATHS *ppaths)
{
    int i, t, v, l0;
    GretlDataFormat fmt;
    char datfile[MAXLEN], hdrfile[MAXLEN], lblfile[MAXLEN];
    int tsamp = sample_size(pdinfo);
    int n = pdinfo->n;
    char delim = 0;
    FILE *fp = NULL;
    int *pmax = NULL;
    int freelist = 0;
    double xx;
    int err = 0;

    gretl_error_clear();

    if (list != NULL && list[0] == 0) {
	return E_ARGS;
    }

    if (list == NULL) {
	list = full_var_list(pdinfo, &l0);
	if (l0 == 0) {
	    return E_ARGS;
	} else if (list == NULL) {
	    return E_ALLOC;
	} else {
	    freelist = 1;
	}
    }

    l0 = list[0];

    fmt = format_from_opt_or_name(opt, fname, &delim);

    fname = gretl_maybe_switch_dir(fname);

    if (fmt == 0 || fmt == GRETL_FMT_GZIPPED) {
	err = gretl_write_gdt(fname, list, Z, pdinfo, 
			      (fmt == GRETL_FMT_GZIPPED)? OPT_Z : OPT_NONE,
			      ppaths);
	goto write_exit;
    }

    if (fmt == GRETL_FMT_DB) {
	err = write_db_data(fname, list, opt, Z, pdinfo);
	goto write_exit;
    }

    if (fmt == GRETL_FMT_CSV && get_csv_delim(pdinfo) == ',' && 
	',' == pdinfo->decpoint) {
	sprintf(gretl_errmsg, _("You can't use the same character for "
				"the column delimiter and the decimal point"));
	err = E_DATA;
	goto write_exit;
    }

    strcpy(datfile, fname);

    /* write header and label files if not exporting to other formats */
    if (fmt != GRETL_FMT_R && fmt != GRETL_FMT_CSV && 
	fmt != GRETL_FMT_OCTAVE && fmt != GRETL_FMT_DAT && 
	fmt != GRETL_FMT_JM) {
	if (!has_suffix(datfile, ".gz")) {
	    switch_ext(hdrfile, datfile, "hdr");
	    switch_ext(lblfile, datfile, "lbl");
	} else {
	    gz_switch_ext(hdrfile, datfile, "hdr");
	    gz_switch_ext(lblfile, datfile, "lbl");
	}
	if (writehdr(hdrfile, list, pdinfo, fmt)) {
	    fprintf(stderr, I_("Write of header file failed"));
	    err = E_FOPEN;
	    goto write_exit;
	}
	if (writelbl(lblfile, list, pdinfo)) {
	    fprintf(stderr, I_("Write of labels file failed"));
	    err = E_FOPEN;
	    goto write_exit;
	}
    }

    /* open file for output */
    fp = gretl_fopen(datfile, "w");
    if (fp == NULL) {
	err = E_FOPEN;
	goto write_exit;
    }

    if (fmt == GRETL_FMT_CSV || fmt == GRETL_FMT_OCTAVE || 
	GRETL_FMT_R || fmt == GRETL_FMT_TRAD || 
	fmt == GRETL_FMT_DAT || fmt == GRETL_FMT_JM) { 
	/* an ASCII variant of some sort */
	pmax = malloc(l0 * sizeof *pmax);
	if (pmax == NULL) {
	    fclose(fp);
	    err = E_ALLOC;
	    goto write_exit;
	}
	for (i=1; i<=l0; i++) {
	    v = list[i];
	    pmax[i-1] = get_precision(&Z[v][pdinfo->t1], tsamp, 10);
	}	
    }

    if (fmt != GRETL_FMT_CSV || pdinfo->decpoint != ',') {
	gretl_push_c_numeric_locale();
    }

    if (fmt == GRETL_FMT_TRAD) { 
	/* plain ASCII */
	for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	    if (pdinfo->markers && pdinfo->S != NULL) {
		fprintf(fp, "%s ", pdinfo->S[t]);
	    }
	    for (i=1; i<=l0; i++) {
		v = list[i];
		if (na(Z[v][t])) {
		    fprintf(fp, "-999 ");
		} else if (pmax[i-1] == PMAX_NOT_AVAILABLE) {
		    fprintf(fp, "%.12g ", Z[v][t]);
		} else {
		    fprintf(fp, "%.*f ", pmax[i-1], Z[v][t]);
		}
	    }
	    fputc('\n', fp);
	}
    } else if (fmt == GRETL_FMT_CSV || fmt == GRETL_FMT_R) { 
	/* export CSV or GNU R (dataframe) */
	int print_obs = 0;

	if (fmt == GRETL_FMT_CSV) {
	    if ((pdinfo->structure == TIME_SERIES || pdinfo->S != NULL)
		&& !(opt & OPT_X)) {
		print_obs = 1;
	    }
	    if (!delim) {
		delim = get_csv_delim(pdinfo);
	    }
	} else {
	    print_obs = (pdinfo->S != NULL);
	    delim = ' ';
	}

	if (fmt == GRETL_FMT_R && dataset_is_time_series(pdinfo)) {
	    char datestr[OBSLEN];

	    ntodate_full(datestr, pdinfo->t1, pdinfo);
	    fprintf(fp, "# time-series data: start = %s, frequency = %d\n",
		    datestr, pdinfo->pd);
	}

	/* variable names */
	if (fmt == GRETL_FMT_CSV && print_obs && 
	    (pdinfo->S != NULL || pdinfo->structure != CROSS_SECTION)) {
	    fprintf(fp, "obs%c", delim);
	}
	for (i=1; i<l0; i++) {
	    fprintf(fp, "%s%c", pdinfo->varname[list[i]], delim);
	}
	fprintf(fp, "%s\n", pdinfo->varname[list[l0]]);
	
	for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	    if (print_obs) {
		if (pdinfo->S != NULL) {
		    fprintf(fp, "\"%s\"%c", pdinfo->S[t], delim);
		} else {
		    char tmp[OBSLEN];

		    ntodate_full(tmp, t, pdinfo);
		    if (quarterly_or_monthly(pdinfo)) {
			modify_date_for_csv(tmp, pdinfo->pd);
		    }
		    fprintf(fp, "%s%c", tmp, delim);
		}
	    }
	    for (i=1; i<=l0; i++) { 
		v = list[i];
		xx = Z[v][t];
		if (na(xx)) {
		    fprintf(fp, "NA");
		} else if (pmax[i-1] == PMAX_NOT_AVAILABLE) {
		    fprintf(fp, "%.12g", xx);
		} else {
		    fprintf(fp, "%.*f", pmax[i-1], xx);
		}
		if (i < l0) {
		    fputc(delim, fp);
		} else {
		    fputc('\n', fp);
		}
	    }
	}
    } else if (fmt == GRETL_FMT_OCTAVE) { 
	/* GNU Octave: write out data as a matrix */
	fprintf(fp, "# name: X\n# type: matrix\n# rows: %d\n# columns: %d\n", 
		n, list[0]);
	for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	    for (i=1; i<=list[0]; i++) {
		v = list[i];
		if (pmax[i-1] == PMAX_NOT_AVAILABLE) {
		    fprintf(fp, "%.12g ", Z[v][t]);
		} else {
		    fprintf(fp, "%.*f ", pmax[i-1], Z[v][t]); 
		}
	    }
	    fputc('\n', fp);
	}
    } else if (fmt == GRETL_FMT_DAT) { 
	/* PcGive: data file with load info */
	int pd = pdinfo->pd;

	for (i=1; i<=list[0]; i++) {
	    fprintf(fp, ">%s ", pdinfo->varname[list[i]]);
	    if (pdinfo->structure == TIME_SERIES &&
		(pd == 1 || pd == 4 || pd == 12)) {
		int maj, min;

		date_maj_min(pdinfo->t1, pdinfo, &maj, &min);
		fprintf(fp, "%d %d ", maj, min);
		date_maj_min(pdinfo->t2, pdinfo, &maj, &min);
		fprintf(fp, "%d %d %d", maj, min, pd);
	    } else {
		fprintf(fp, "%d 1 %d 1 1", pdinfo->t1, pdinfo->t2);
	    }
			   
	    fputc('\n', fp);

	    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
		v = list[i];
		xx = Z[v][t];
		if (na(xx)) {
		    fprintf(fp, "-9999.99");
		} else if (pmax[i-1] == PMAX_NOT_AVAILABLE) {
		    fprintf(fp, "%.12g", xx);;
		} else {
		    fprintf(fp, "%.*f", pmax[i-1], xx);
		}
		fputc('\n', fp);
	    }
	    fputc('\n', fp);
	}
    } else if (fmt == GRETL_FMT_JM) { 
	/* JMulti: ascii with comments and date info */
	int maj, min;

	fputs("/*\n", fp);
	for (i=1; i<=list[0]; i++) {
	    fprintf(fp, " %s: %s\n", pdinfo->varname[list[i]], VARLABEL(pdinfo, i));
	}
	fputs("*/\n", fp);
	date_maj_min(pdinfo->t1, pdinfo, &maj, &min);
	if (pdinfo->pd == 4 || pdinfo->pd == 12) {
	    fprintf(fp, "<%d %c%d>\n", maj, (pdinfo->pd == 4)? 'Q' : 'M', min);
	} else if (pdinfo->pd == 1) {
	    fprintf(fp, "<%d>\n", maj);
	} else {
	    fputs("<1>\n", fp);
	}
	for (i=1; i<=list[0]; i++) {
	    fprintf(fp, " %s", pdinfo->varname[list[i]]);
	}
	fputc('\n', fp);
	for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	    for (i=1; i<=list[0]; i++) {
		v = list[i];
		if (na(Z[v][t])) {
		    fputs("NaN ", fp);
		} else if (pmax[i-1] == PMAX_NOT_AVAILABLE) {
		    fprintf(fp, "%.12g ", Z[v][t]);
		} else {
		    fprintf(fp, "%.*f ", pmax[i-1], Z[v][t]);
		}
	    }
	    fputc('\n', fp);
	}
    }

    if (fmt != GRETL_FMT_CSV || pdinfo->decpoint != ',') {
	gretl_pop_c_numeric_locale();
    }

    if (pmax != NULL) {
	free(pmax);
    }

    if (fp != NULL) {
	fclose(fp);
    }

 write_exit:

    if (freelist) {
	free(list);
    }

    return err;
}

static void dataset_type_string (char *str, const DATAINFO *pdinfo)
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
    if (custom_time_series(pdinfo)) {
	strcpy(str, _("special"));
    } else {
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
	case 6:
	case 7:
	    strcpy(str, _("daily")); break;
	case 10:
	    strcpy(str, _("decennial")); break;
	default:
	    strcpy(str, _("unknown")); break;
	}
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
    char tstr[48];
    int i;

    ntodate_full(startdate, 0, pdinfo);
    ntodate_full(enddate, pdinfo->n - 1, pdinfo);

    sprintf(tmp, _("Data file %s\nas of"), 
	    strlen(ppaths->datfile)? ppaths->datfile : _("(unsaved)"));

    print_time(tstr);
    pprintf(prn, "%s %s\n\n", tmp, tstr);

    if (pdinfo->descrip != NULL && *pdinfo->descrip != '\0') {
	pprintf(prn, "%s:\n\n", _("Description"));
	pputs(prn, pdinfo->descrip);
	pputs(prn, "\n\n");
    }

    dataset_type_string(tmp, pdinfo);
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

/* read data "labels" from file */

static int readlbl (const char *lblfile, DATAINFO *pdinfo)
{
    FILE * fp;
    char line[MAXLEN], varname[VNAMELEN];
    char *p;
    int v;
    
    gretl_error_clear();

    fp = gretl_fopen(lblfile, "r");
    if (fp == NULL) {
	return 0;
    }

    while (fgets(line, MAXLEN, fp)) {
	tailstrip(line);
        if (sscanf(line, "%s", varname) != 1) {
	    sprintf(gretl_errmsg, _("Bad data label in %s"), lblfile); 
            break;
        }
	v = series_index(pdinfo, varname);
	if (v < pdinfo->v) {
	    p = line + strlen(varname);
	    p += strspn(p, " \t");
	    VARLABEL(pdinfo, v)[0] = '\0';
	    strncat(VARLABEL(pdinfo, v), p, MAXLABEL - 1);
	} else {
	    fprintf(stderr, I_("extraneous label for var '%s'\n"), varname);
	}
    }

    fclose(fp);

    return 0;
}

static int writelbl (const char *lblfile, const int *list, 
		     const DATAINFO *pdinfo)
{
    FILE *fp;
    int i, lblcount = 0;

    for (i=1; i<=list[0]; i++) {
	if (list[i] == 0) {
	    continue;
	}
	if (strlen(VARLABEL(pdinfo, list[i])) > 2) {
	    lblcount++;
	    break;
	}
    }

    if (lblcount == 0) return 0;

    fp = gretl_fopen(lblfile, "w");
    if (fp == NULL) return 1;

    /* spit out varnames and labels (if filled out) */
    for (i=1; i<=list[0]; i++) {
	if (list[i] == 0) {
	    continue;
	}
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
    if (j > 0 && k < strlen(targ) && k > j) {
	i = k;
    }

    targ[i] = '.';
    targ[i + 1] = '\0';
    strcat(targ, ext);
}

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
 * gretl_get_data:
 * @datfile: name of file to try.
 * @ppaths: path information struct.
 * @pZ: pointer to data set.
 * @pdinfo: pointer to data information struct.
 * @opt: for use with "append".
 * @prn: where messages should be written.
 * 
 * Read data from file into gretl's work space, allocating space as
 * required.
 * 
 * Returns: 0 on successful completion, non-zero otherwise.
 */

int gretl_get_data (char *datfile, PATHS *ppaths,
		    double ***pZ, DATAINFO *pdinfo, 
		    gretlopt opt, PRN *prn) 
{
    DATAINFO *tmpdinfo = NULL;
    double **tmpZ = NULL;
    FILE *dat = NULL;
    gzFile fz = NULL;
    char hdrfile[MAXLEN], lblfile[MAXLEN];
    int newdata = (*pZ == NULL);
    int gdtsuff, gzsuff = 0;
    int binary = 0, old_byvar = 0;
    int err = 0;

    gretl_error_clear();

    *hdrfile = '\0';

    gdtsuff = has_suffix(datfile, ".gdt");
    if (!gdtsuff) {
	gzsuff = has_suffix(datfile, ".gz");
    }

    if (addpath(datfile, ppaths, 0) == NULL) { /* not found yet */
	char tryfile[MAXLEN];
	int found = 0;

	if (!gdtsuff) {
	    /* try using the .gdt suffix? */
	    *tryfile = '\0';
	    strncat(tryfile, datfile, MAXLEN-1);
	    try_gdt(tryfile); 
	    found = (addpath(tryfile, ppaths, 0) != NULL);
	    if (found) {
		gdtsuff = 1;
	    }
	}

	/* or maybe the file is gzipped but lacks a .gz extension? */
	if (!found && !gzsuff) { 
	    sprintf(tryfile, "%s.gz", datfile);
	    if (addpath(tryfile, ppaths, 0) != NULL) {
		gzsuff = 1;
		found = 1;
	    }
	}

	if (!found) {
	    sprintf(gretl_errmsg, _("Couldn't open file %s"),  datfile);
	    return E_FOPEN;
	} else {
	    strcpy(datfile, tryfile);
	}
    }

    /* catch XML files that have strayed in here? */
    if (gdtsuff && gretl_is_xml_file(datfile)) {
	return gretl_read_gdt(datfile, ppaths, pZ, pdinfo, OPT_NONE, prn);
    }

    tmpdinfo = datainfo_new();
    if (tmpdinfo == NULL) {
	return E_ALLOC;
    }
	
    if (!gzsuff) {
	switch_ext(hdrfile, datfile, "hdr");
	switch_ext(lblfile, datfile, "lbl");
    } else {
	gz_switch_ext(hdrfile, datfile, "hdr");
	gz_switch_ext(lblfile, datfile, "lbl");
    }

    /* read data header file */
    err = readhdr(hdrfile, tmpdinfo, &binary, &old_byvar);
    if (err == E_FOPEN) {
	/* no header file, so maybe it's just an ascii datafile */
	return import_csv(datfile, pZ, pdinfo, OPT_NONE, prn);
    } else if (err) {
	return err;
    } else { 
	pprintf(prn, I_("\nReading header file %s\n"), hdrfile);
    }

    /* deal with case where first col. of data file contains
       "marker" strings */
    tmpdinfo->S = NULL;
    if (tmpdinfo->markers && dataset_allocate_obs_markers(tmpdinfo)) {
	return E_ALLOC; 
    }
    
    /* allocate dataset */
    if (allocate_Z(&tmpZ, tmpdinfo)) {
	err = E_ALLOC;
	goto bailout;
    }

    /* Invoke data (Z) reading function */
    if (gzsuff) {
	fz = gretl_gzopen(datfile, "rb");
	if (fz == NULL) {
	    err = E_FOPEN;
	    goto bailout;
	}
    } else {
	if (binary) {
	    dat = gretl_fopen(datfile, "rb");
	} else {
	    dat = gretl_fopen(datfile, "r");
	}
	if (dat == NULL) {
	    err = E_FOPEN;
	    goto bailout;
	}
    }

    if (gzsuff) {
	err = gz_readdata(fz, tmpdinfo, tmpZ, binary); 
	gzclose(fz);
    } else {
	err = readdata(dat, tmpdinfo, tmpZ, binary, old_byvar); 
	fclose(dat);
    }

    if (err) goto bailout;

    if (tmpdinfo->structure == STACKED_CROSS_SECTION) {
	err = switch_panel_orientation(tmpZ, tmpdinfo);
    }

    if (!err && tmpdinfo->structure == STACKED_TIME_SERIES) {
	err = dataset_add_default_panel_indices(tmpdinfo);
    }

    if (err) goto bailout;

    /* print out basic info from the files read */
    pprintf(prn, I_("periodicity: %d, maxobs: %d,\n"
	   "observations range: %s-%s\n"), tmpdinfo->pd, tmpdinfo->n,
	   tmpdinfo->stobs, tmpdinfo->endobs);

    pputs(prn, I_("\nReading "));
    pputs(prn, (tmpdinfo->structure == TIME_SERIES) ? 
	    I_("time-series") : _("cross-sectional"));
    pputs(prn, I_(" datafile"));
    if (strlen(datfile) > 40) {
	pputc(prn, '\n');
    }
    pprintf(prn, " %s\n\n", datfile);

    /* Set sample range to entire length of dataset by default */
    tmpdinfo->t1 = 0; 
    tmpdinfo->t2 = tmpdinfo->n - 1;

    err = readlbl(lblfile, tmpdinfo);
    if (err) goto bailout;

    err = merge_or_replace_data(pZ, pdinfo, &tmpZ, &tmpdinfo, opt, prn);

    if (!err && newdata && ppaths != NULL && datfile != ppaths->datfile) {
	strcpy(ppaths->datfile, datfile);
    }

 bailout:

    if (err && tmpdinfo != NULL) {
	destroy_dataset(tmpZ, tmpdinfo);
    }

    return err;
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
    if (data_status) {
	clear_datainfo(pdinfo, CLEAR_FULL);
    }

    /* dummy up the data info */
    pdinfo->n = length;
    pdinfo->v = 2;
    dataset_obs_info_default(pdinfo);

    if (dataset_allocate_varnames(pdinfo)) {
	return E_ALLOC;
    }

    /* allocate dataset */
    if (allocate_Z(pZ, pdinfo)) {
	return E_ALLOC;
    }

    /* add an index var */
    strcpy(pdinfo->varname[1], "index");
    strcpy(VARLABEL(pdinfo, 1), _("index variable"));
    for (t=0; t<pdinfo->n; t++) {
	(*pZ)[1][t] = (double) (t + 1);
    }

    /* print out basic info */
    pprintf(prn, M_("periodicity: %d, maxobs: %d,\n"
	   "observations range: %s-%s\n"), pdinfo->pd, pdinfo->n,
	   pdinfo->stobs, pdinfo->endobs);

    /* Set sample range to entire length of data-set by default */
    pdinfo->t1 = 0; 
    pdinfo->t2 = pdinfo->n - 1;

    return 0;
}

static int extend_markers (DATAINFO *pdinfo, int old_n, int new_n)
{
    char **S = realloc(pdinfo->S, new_n * sizeof *S);
    int t, err = 0;
	   
    if (S == NULL) {
	err = 1;
    } else {
	pdinfo->S = S;
	for (t=old_n; t<new_n && !err; t++) {
	    S[t] = malloc(OBSLEN);
	    if (S[t] == NULL) {
		err = 1;
	    } 
	}
    }

    return err;
}

static void merge_error (char *msg, PRN *prn)
{
    pputs(prn, msg);
    strcpy(gretl_errmsg, msg);
}

static int count_new_vars (const DATAINFO *pdinfo, const DATAINFO *addinfo,
			   PRN *prn)
{
    const char *newname;
    /* default to all new, and subtract */
    int addvars = addinfo->v - 1;
    int i, j;

    for (i=1; i<addinfo->v && addvars >= 0; i++) {
	newname = addinfo->varname[i];
	if (get_matrix_by_name(newname)) {
	    merge_error("can't replace matrix with series\n", prn);
	    addvars = -1;
	} else if (get_string_by_name(newname)) {
	    merge_error("can't replace string with series\n", prn);
	    addvars = -1;
	} else {
	    for (j=1; j<pdinfo->v; j++) {
		/* FIXME collision with scalar, matrix names */
		if (!strcmp(newname, pdinfo->varname[j])) {
		    addvars--;
		}
	    }
	}
    }

    return addvars;
}

static int compare_ranges (const DATAINFO *targ,
			   const DATAINFO *src,
			   int *offset)
{
    int ed0, sd1, ed1;
    int addobs = -1;

    ed0 = dateton(targ->endobs, targ);
    sd1 = merge_dateton(src->stobs, targ);
    ed1 = merge_dateton(src->endobs, targ);

#if 0
    fprintf(stderr, "compare_ranges:\n"
	    " targ->n = %d, src->n = %d\n"
	    " targ->stobs = '%s', src->stobs = '%s'\n" 
	    " sd1 = %d, ed1 = %d\n",
	    targ->n, src->n, targ->stobs, src->stobs,
	    sd1, ed1);
#endif

    if (sd1 < 0) {
	/* case: new data start earlier than old */
	if (ed1 < 0) {
	    fprintf(stderr, "no overlap in ranges, can't merge\n");
	} else if (ed1 > ed0) {
	    fprintf(stderr, "new data start earlier, end later, can't handle\n");
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

    if (addobs < 0) {
	fputs("compare_ranges: returning error\n", stderr);
    }

    return addobs;
}

static int panel_expand_ok (DATAINFO *pdinfo, DATAINFO *addinfo,
			    gretlopt opt)
{
    int n = pdinfo->paninfo->nunits;
    int T = pdinfo->paninfo->Tmax;
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
				 double ***pZ, 
				 DATAINFO *pdinfo, 
				 double **addZ, 
				 DATAINFO *addinfo,
				 gretlopt opt,
				 PRN *prn)
{
    int n = pdinfo->paninfo->nunits;
    int T = pdinfo->paninfo->Tmax;
    int k = pdinfo->v;
    int tsdata;
    int i, j, s, p, t;
    int err = 0;

    if (addvars > 0 && dataset_add_series(addvars, pZ, pdinfo)) {
	merge_error(_("Out of memory!\n"), prn);
	err = E_ALLOC;
    }

    tsdata = ((opt & OPT_T) || addinfo->n != n);

    for (i=1; i<addinfo->v && !err; i++) {
	int v = series_index(pdinfo, addinfo->varname[i]);

	if (v >= k) {
	    /* a new variable */
	    v = k++;
	    strcpy(pdinfo->varname[v], addinfo->varname[i]);
	    copy_varinfo(pdinfo->varinfo[v], addinfo->varinfo[i]);
	} 

	s = 0;
	for (j=0; j<n; j++) {
	    /* loop across units */
	    for (t=0; t<T; t++) {
		/* loop across periods */
		p = (tsdata)? t : j;
		(*pZ)[v][s++] = addZ[i][p]; 
	    }
	}
    }

    return err;
}

static int 
just_append_rows (const DATAINFO *targ, const DATAINFO *src,
		  int *offset)
{
    if (targ->structure == CROSS_SECTION &&
	src->structure == CROSS_SECTION &&
	targ->markers == 0 && src->markers == 0 &&
	targ->sd0 == 1 && src->sd0 == 1) {
	*offset = targ->n;
	return src->n;
    } else {
	return 0;
    }
}

static int simple_range_match (const DATAINFO *targ, const DATAINFO *src,
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

#if 0
static int markers_are_ints (const DATAINFO *pdinfo)
{
    char *test;
    int i;

    errno = 0;

    for (i=0; i<pdinfo->n; i++) {
	strtol(pdinfo->S[i], &test, 10);
	if (*test || errno) {
	    errno = 0;
	    return 0;
	}
    }

    return 1;
}
#endif

#define simple_structure(p) (p->structure == TIME_SERIES ||		\
			     p->structure == SPECIAL_TIME_SERIES ||	\
			     (p->structure == CROSS_SECTION &&		\
			      p->S == NULL))

/**
 * merge_data:
 * @pZ: pointer to data set.
 * @pdinfo: data information struct.
 * @addZ: new data set to be merged in.
 * @addinfo: data information associated with @addZ.
 * @opt: may include %OPT_T to force a time-series interpretation
 * when appending to a panel dataset.
 * @prn: print struct to accept messages.
 * 
 * Attempt to merge the content of a newly opened data file into
 * gretl's current working data set.  
 * 
 * Returns: 0 on successful completion, non-zero otherwise.
 */

static int merge_data (double ***pZ, DATAINFO *pdinfo,
		       double **addZ, DATAINFO *addinfo,
		       gretlopt opt, PRN *prn)
{
    int addsimple = 0;
    int addpanel = 0;
    int addvars = 0;
    int addobs = 0;
    int offset = 0;
    int err = 0;

    /* first see how many new vars we have */
    addvars = count_new_vars(pdinfo, addinfo, prn);
    if (addvars < 0) {
	return 1;
    }

    /* below: had additional condition: simple_structure(pdinfo)
       relaxed this on 2009-05-15 */

    if (simple_range_match(pdinfo, addinfo, &offset)) {
	/* we'll allow undated data to be merged with the existing
	   dateset, sideways, provided the number of observations
	   matches OK */
	addsimple = 1;
    } else if (dataset_is_panel(pdinfo) && 
	       panel_expand_ok(pdinfo, addinfo, opt)) {
	/* allow appending to panel when the number of obs matches
	   either the cross-section size or the time-series length */
	addpanel = 1;
    } else if (pdinfo->pd != addinfo->pd) {
	merge_error(_("Data frequency does not match\n"), prn);
	err = 1;
    }

    if (!err) {
	if (!addsimple && !addpanel) {
	    addobs = compare_ranges(pdinfo, addinfo, &offset);
	}
	if (addobs <= 0 && addvars == 0) {
	    addobs = just_append_rows(pdinfo, addinfo, &offset);
	}
    }

    if (!err && (addobs < 0 || addvars < 0)) {
	merge_error(_("New data not conformable for appending\n"), prn);
	err = 1;
    }

    if (!err && !addpanel && pdinfo->markers != addinfo->markers) {
	if (addinfo->n != pdinfo->n) {
	    merge_error(_("Inconsistency in observation markers\n"), prn);
	    err = 1;
	} else if (addinfo->markers && !pdinfo->markers) {
	    dataset_destroy_obs_markers(addinfo);
	}
    }

#if 0
    fprintf(stderr, "merge_data: addvars = %d, addobs = %d\n",
	    addvars, addobs);
#endif

    /* if checks are passed, try merging the data */

    if (!err && addobs > 0) { 
	int i, t, new_n = pdinfo->n + addobs;

	if (pdinfo->markers) {
	    err = extend_markers(pdinfo, pdinfo->n, new_n);
	    if (!err) {
		for (t=pdinfo->n; t<new_n; t++) {
		    strcpy(pdinfo->S[t], addinfo->S[t - offset]);
		}
	    }
	}

	for (i=0; i<pdinfo->v && !err; i++) {
	    double *x;

	    x = realloc((*pZ)[i], new_n * sizeof *x);
	    if (x == NULL) {
		err = 1;
		break;
	    }

	    for (t=pdinfo->n; t<new_n; t++) {
		if (i == 0) {
		    x[t] = 1.0;
		} else {
		    x[t] = NADBL;
		}
	    }
	    (*pZ)[i] = x;
	}

	if (err) { 
	    merge_error(_("Out of memory!\n"), prn);
	} else {
	    pdinfo->n = new_n;
	    ntodate_full(pdinfo->endobs, new_n - 1, pdinfo);
	    pdinfo->t2 = pdinfo->n - 1;
	}
    }

    if (!err && addpanel) {
	err = panel_append_special(addvars, pZ, pdinfo, addZ, addinfo, 
				   opt, prn);
    } else if (!err) { 
	int k = pdinfo->v;
	int i, t;

	if (addvars > 0 && dataset_add_series(addvars, pZ, pdinfo)) {
	    merge_error(_("Out of memory!\n"), prn);
	    err = E_ALLOC;
	}

	for (i=1; i<addinfo->v && !err; i++) {
	    int v = series_index(pdinfo, addinfo->varname[i]);
	    int newvar = 0;

	    if (v >= k) {
		/* a new variable */
		v = k++;
		newvar = 1;
		strcpy(pdinfo->varname[v], addinfo->varname[i]);
		copy_varinfo(pdinfo->varinfo[v], addinfo->varinfo[i]);
	    } 

	    for (t=0; t<pdinfo->n; t++) {
		if (t >= offset && t - offset < addinfo->n) {
		    (*pZ)[v][t] = addZ[i][t - offset];
		} else if (newvar) {
		    (*pZ)[v][t] = NADBL;
		}
	    }
	}
    }

    if (!err && (addvars || addobs) && gretl_messages_on()) {
	pputs(prn, _("Data appended OK\n"));
    }

    return err;
}

/**
 * merge_or_replace_data:
 * @pZ0: pointer to original data set.
 * @pdinfo0: original dataset information struct.
 * @pZ1: new data set.
 * @ppdinfo1: pointer to dataset information associated with @pZ1.
 * @opt: may include %OPT_T when appending to a panel dataset,
 * to force a time-series interpretation of the added data.
 * @prn: print struct to accept messages.
 *
 * Given a newly-created dataset, pointed to by @pZ1 and
 * @ppdinfo1, either attempt to merge it with @pZ0, if the
 * original dataset is non-NULL, or replace the content of
 * the original pointers with the new dataset.
 * In case merging is not successful, the new dataset is
 * destroyed.
 * 
 * Returns: 0 on successful completion, non-zero otherwise.
 */

int merge_or_replace_data (double ***pZ0, DATAINFO *pdinfo0,
			   double ***pZ1, DATAINFO **ppdinfo1,
			   gretlopt opt, PRN *prn)
{
    int err = 0;

    if (*pZ0 != NULL) {
	err = merge_data(pZ0, pdinfo0, *pZ1, *ppdinfo1, opt, prn);
	destroy_dataset(*pZ1, *ppdinfo1);
    } else {
	*pdinfo0 = **ppdinfo1;
	free(*ppdinfo1);
	*pZ0 = *pZ1;
    }

    *pZ1 = NULL;
    *ppdinfo1 = NULL;

    return err;
}

static int check_marker (char *src, int i)
{
    int err = 0;

    if (!g_utf8_validate(src, -1, NULL)) {
	gchar *trstr = NULL;
	gsize bytes;

	trstr = g_locale_to_utf8(src, -1, NULL, &bytes, NULL);

	if (trstr == NULL) {
	    sprintf(gretl_errmsg, "Invalid characters in marker, line %d", i);
	    err = E_DATA;
	} else {
	    *src = '\0';
	    strncat(src, trstr, OBSLEN - 1);
	    g_free(trstr);
	}
    }

    return err;
}

/**
 * add_obs_markers_from_file:
 * @pdinfo: data information struct.
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

int add_obs_markers_from_file (DATAINFO *pdinfo, const char *fname)
{
    char **S = NULL;
    FILE *fp;
    char line[128], marker[32];
    int t, err = 0;

    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	return E_FOPEN;
    }

    S = strings_array_new_with_length(pdinfo->n, OBSLEN);
    if (S == NULL) {
	fclose(fp);
	return E_ALLOC;
    }
    
    for (t=0; t<pdinfo->n && !err; t++) {
	if (fgets(line, sizeof line, fp) == NULL) {
	    sprintf(gretl_errmsg, "Expected %d markers; found %d\n", 
		    pdinfo->n, t);
	    err = E_DATA;
	} else if (sscanf(line, "%31[^\n\r]", marker) != 1) {
	    sprintf(gretl_errmsg, "Couldn't read marker on line %d", t+1);
	    err = E_DATA;
	} else {
	    strncat(S[t], marker, OBSLEN - 1);
	    err = check_marker(S[t], t+1);
	}
    }

    if (err) {
	free_strings_array(S, pdinfo->n);
    } else {
	if (pdinfo->S != NULL) {
	    free_strings_array(pdinfo->S, pdinfo->n);
	} 
	pdinfo->markers = REGULAR_MARKERS;
	pdinfo->S = S;
    }

    return err;
}

static void 
octave_varname (char *name, const char *s, int nnum, int v)
{
    char nstr[8];
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
    int c, c1, cbak = 0, cc = 0;
    int comment = 0, maxlen = 0;

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
	cbak = c;
	if (!isspace((unsigned char) c) && !isprint((unsigned char) c) &&
	    !(c == CTRLZ)) {
	    pprintf(prn, M_("Binary data (%d) encountered: this is not a valid "
			   "text file\n"), c);
	    return -1;
	}
	if (cc == 0) {
	    comment = (c == '#');
	}
	cc++;
    }

    if (maxlen == 0) {
	pprintf(prn, M_("Data file is empty\n"));
    } 

    if (maxlen > 0) {
	/* allow for newline and null terminator */
	maxlen += 3;
    }

    return maxlen;
}

static int 
import_octave (const char *fname, double ***pZ, DATAINFO *pdinfo, 
	       gretlopt opt, PRN *prn)
{
    DATAINFO *octinfo = NULL;
    double **octZ = NULL;
    FILE *fp = NULL;
    char *line = NULL;
    char tmp[8], name[32];
    int nrows = 0, ncols = 0, nblocks = 0;
    int brows = 0, bcols = 0, oldbcols = 0;
    int maxlen, got_type = 0, got_name = 0;
    int i, t, err = 0;

    pprintf(prn, "%s %s...\n", M_("parsing"), fname);

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

    pprintf(prn, M_("   longest line: %d characters\n"), maxlen - 1);

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
			pprintf(prn, M_("   Found matrix '%s' with "
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
	pputs(prn, M_("Invalid data file\n"));
	err = E_DATA;
	goto oct_bailout;
    } 

    /* initialize datainfo and Z */

    octinfo = datainfo_new();
    if (octinfo == NULL) {
	pputs(prn, M_("Out of memory!\n"));
	err = E_ALLOC;
	goto oct_bailout;
    }

    octinfo->n = nrows;
    octinfo->v = ncols + 1;

    if (start_new_Z(&octZ, octinfo, 0)) {
	pputs(prn, M_("Out of memory!\n"));
	err = E_ALLOC;
	goto oct_bailout;
    }  

    rewind(fp);

    pprintf(prn, M_("   number of variables: %d\n"), ncols);
    pprintf(prn, M_("   number of observations: %d\n"), nrows);
    pprintf(prn, M_("   number of data blocks: %d\n"), nblocks); 

    i = 1;
    t = 0;

    while (fgets(line, maxlen, fp) && !err) {
	char *s = line;
	int j;

	if (*s == '#') {
	    if (sscanf(line, "# name: %15s", name) == 1) {
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

	if (t >= octinfo->n) {
	    err = 1;
	}

	for (j=0; j<bcols && !err; j++) {
	    double x;
	    int v = i + j;

	    if (t == 0) {
		int nnum = (bcols > 1)? j + 1 : 0;

		octave_varname(octinfo->varname[i+j], name, nnum, v);
	    }

	    while (isspace(*s)) s++;
	    if (sscanf(s, "%lf", &x) != 1) {
		fprintf(stderr, "error: '%s', didn't get double\n", s);
		err = 1;
	    } else {
		octZ[v][t] = x;
		while (!isspace(*s)) s++;
	    }	
	}
	t++;
    }

    if (err) {
	pputs(prn, M_("Invalid data file\n"));
	err = E_DATA;
	goto oct_bailout;
    } 

    err = merge_or_replace_data(pZ, pdinfo, &octZ, &octinfo, opt, prn);

 oct_bailout:

    if (fp != NULL) {
	fclose(fp);
    }

    if (line != NULL) {
	free(line);
    }

    if (octinfo != NULL) {
	clear_datainfo(octinfo, CLEAR_FULL);
    }

    console_off();

    return err;
}

/**
 * import_other:
 * @fname: name of file.
 * @ftype: type of data file.
 * @pZ: pointer to data set.
 * @pdinfo: pointer to data information struct.
 * @opt: for use with "append".
 * @prn: gretl printing struct.
 * 
 * Open a data file of a type that requires a special plugin.
 * 
 * Returns: 0 on successful completion, non-zero otherwise.
 */

int import_other (const char *fname, int ftype,
		  double ***pZ, DATAINFO *pdinfo, 
		  gretlopt opt, PRN *prn)
{
    void *handle;
    FILE *fp;
    int (*importer) (const char *, 
		     double ***, DATAINFO *, 
		     gretlopt, PRN *);
    int err = 0;

    check_for_console(prn);

    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	pprintf(prn, M_("Couldn't open %s\n"), fname);
	err = E_FOPEN;
	goto bailout;
    }

    fclose(fp);

    if (ftype == GRETL_OCTAVE) {
	/* plugin not needed */
	return import_octave(fname, pZ, pdinfo, opt, prn);
    }

    if (ftype == GRETL_WF1) {
	importer = get_plugin_function("wf1_get_data", &handle);
    } else if (ftype == GRETL_DTA) {
	importer = get_plugin_function("dta_get_data", &handle);
    } else if (ftype == GRETL_SAV) {
	importer = get_plugin_function("sav_get_data", &handle);
    } else if (ftype == GRETL_JMULTI) {
	importer = get_plugin_function("jmulti_get_data", &handle);
    } else {
	pprintf(prn, M_("Unrecognized data type"));
	pputc(prn, '\n');
	return E_DATA;
    }

    if (importer == NULL) {
        err = 1;
    } else {
	err = (*importer)(fname, pZ, pdinfo, opt, prn);
	close_plugin(handle);
    }

 bailout:

    console_off();

    return err;
}

/**
 * import_spreadsheet:
 * @fname: name of file.
 * @ftype: type of data file.
 * @list: list of parameters for spreadsheet import, or %NULL.
 * @sheetname: name of worksheet, or %NULL.
 * @pZ: pointer to data set.
 * @pdinfo: pointer to data information struct.
 * @opt: fir use with "append".
 * @prn: gretl printing struct.
 * 
 * Open a data file of a type that requires a special plugin.
 * 
 * Returns: 0 on successful completion, non-zero otherwise.
 */

int import_spreadsheet (const char *fname, int ftype, 
			int *list, char *sheetname,
			double ***pZ, DATAINFO *pdinfo, 
			gretlopt opt, PRN *prn)
{
    void *handle;
    FILE *fp;
    int (*importer) (const char*, int *, char *,
		     double ***, DATAINFO *, 
		     gretlopt, PRN *);
    int err = 0;

    check_for_console(prn);

    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	pprintf(prn, M_("Couldn't open %s\n"), fname);
	err = E_FOPEN;
	goto bailout;
    }

    fclose(fp);

    if (ftype == GRETL_GNUMERIC) {
	importer = get_plugin_function("gnumeric_get_data", &handle);
    } else if (ftype == GRETL_XLS) {
	importer = get_plugin_function("xls_get_data", &handle);
    } else if (ftype == GRETL_ODS) {
	importer = get_plugin_function("ods_get_data", &handle);
    } else {
	pprintf(prn, M_("Unrecognized data type"));
	pputc(prn, '\n');
	return E_DATA;
    }

    if (importer == NULL) {
        err = 1;
    } else {
	err = (*importer)(fname, list, sheetname, pZ, pdinfo, opt, prn);
	close_plugin(handle);
    }

 bailout:

    console_off();

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

int gretl_is_pkzip_file (const char *fname)
{
    FILE *fp;
    char test[3] = {0};
    int ret = 0;

    fp = gretl_fopen(fname, "rb");
    if (fp != NULL) {
	if (fread(test, 1, 2, fp) == 2) {
	    if (!strcmp(test, "PK")) ret = 1;
	} 
	fclose(fp);
    } 

    return ret;
}

/**
 * detect_filetype:
 * @fname: name of file to examine.
 * @ppaths: path information struct.
 * 
 * Attempt to determine the type of a file to be opened in gretl:
 * data file (of various formats), or command script.
 * 
 * Returns: integer code indicating the type of file.
 */

GretlFileType detect_filetype (char *fname, PATHS *ppaths)
{
    int i, c, ftype = GRETL_NATIVE_DATA;
    char teststr[5];
    FILE *fp;

    /* might be a script file? (watch out for DOS-mangled names) */
    if (has_suffix(fname, ".inp")) { 
	return GRETL_SCRIPT;
    }

    if (has_suffix(fname, ".gretl")) {
	if (gretl_is_pkzip_file(fname)) {
	    return GRETL_SESSION;
	} else {
	    return GRETL_SCRIPT;
	}
    }

    if (has_suffix(fname, ".gnumeric"))
	return GRETL_GNUMERIC;
    if (has_suffix(fname, ".xls"))
	return GRETL_XLS;
    if (has_suffix(fname, ".ods"))
	return GRETL_ODS;
    if (has_suffix(fname, ".wf1"))
	return GRETL_WF1;
    if (has_suffix(fname, ".dta"))
	return GRETL_DTA;
    if (has_suffix(fname, ".sav"))
	return GRETL_SAV;
    if (has_suffix(fname, ".bin"))
	return GRETL_NATIVE_DB;
    if (has_suffix(fname, ".rat"))
	return GRETL_RATS_DB;
    if (has_suffix(fname, ".csv"))
	return GRETL_CSV;
    if (has_suffix(fname, ".txt"))
	return GRETL_CSV;
    if (has_suffix(fname, ".m"))
	return GRETL_OCTAVE;
    if (has_suffix(fname, ".bn7"))
	return GRETL_PCGIVE_DB;

    if (ppaths == NULL) {
	return GRETL_NATIVE_DATA; 
    }

    addpath(fname, ppaths, 0); 

    if (gretl_is_xml_file(fname)) {
	return GRETL_XML_DATA;  
    } 

    if (has_suffix(fname, ".dat") && is_jmulti_datafile(fname)) {
	return GRETL_JMULTI; 
    }

    fp = gretl_fopen(fname, "r");
    if (fp == NULL) { 
	/* may be native file in different location */
	return GRETL_NATIVE_DATA; 
    }

    /* take a peek at content */
    for (i=0; i<80; i++) {
	c = getc(fp);
	if (c == EOF || c == '\n') {
	    break;
	}
	if (!isprint(c) && c != '\r' && c != '\t') {
	    ftype = GRETL_NATIVE_DATA; /* native binary data? */
	    break;
	}
	if (i < 4) {
	    teststr[i] = c;
	}
    }

    fclose(fp);
    teststr[4] = 0;

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

    if (errno == ERANGE || val <= INT_MIN || val >= INT_MAX) {
	sprintf(gretl_errmsg, M_("'%s' -- number out of range!"), numstr);
    }

    return 1;
}

/**
 * transpose_data:
 * @pZ: pointer to data array.
 * @pdinfo: pointer to dataset information struct.
 *
 * Attempts to transpose the current dataset, so that each
 * variable becomes interpreted as an observation and each
 * observation as a variable.  This will not work if the
 * dataset contains scalar variables.
 *
 * Returns: 0 on success, non-zero error code on error.
 */

int transpose_data (double ***pZ, DATAINFO *pdinfo)
{
    double **tZ = NULL;
    DATAINFO *tinfo;
    int k = pdinfo->n + 1;
    int T = pdinfo->v - 1;
    int i, t;

    tinfo = create_new_dataset(&tZ, k, T, 0);
    if (tinfo == NULL) {
	return E_ALLOC;
    }

    for (i=1; i<pdinfo->v; i++) {
	for (t=0; t<pdinfo->n; t++) {
	    tZ[t+1][i-1] = (*pZ)[i][t];
	}
    }

    for (t=0; t<pdinfo->n; t++) {
	if (pdinfo->S != NULL && pdinfo->S[t][0] != '\0') {
	    tinfo->varname[t+1][0] = '\0';
	    strncat(tinfo->varname[t+1], pdinfo->S[t], VNAMELEN - 1);
	} else {
	    sprintf(tinfo->varname[t+1], "v%d", t+1);
	}
    }

    free_Z(*pZ, pdinfo);
    *pZ = tZ;

    clear_datainfo(pdinfo, CLEAR_FULL);

    pdinfo->v = k;
    pdinfo->n = T;
    pdinfo->t1 = 0;
    pdinfo->t2 = pdinfo->n - 1;

    pdinfo->varname = tinfo->varname;
    pdinfo->varinfo = tinfo->varinfo;

    dataset_obs_info_default(pdinfo);

    free(tinfo);

    return 0;
}

void dataset_set_regular_markers (DATAINFO *pdinfo)
{
    pdinfo->markers = REGULAR_MARKERS;
}

struct filetype_info {
    GretlFileType type;
    const char *src;
};

#define NOTELEN 256

/* on successful import of data from some "foreign" format,
   add a note to the "descrip" member of the new dataset
   saying ehere it came from and when
*/

void dataset_add_import_info (DATAINFO *pdinfo, const char *fname,
			      GretlFileType type)
{
    struct filetype_info ftypes[] = {
	{ GRETL_CSV,      "CSV" },
	{ GRETL_GNUMERIC, "Gnumeric" },
	{ GRETL_XLS,      "Excel" },
	{ GRETL_WF1,      "Eviews" },
	{ GRETL_DTA,      "Stata" },
	{ GRETL_SAV,      "SPSS" },
	{ GRETL_JMULTI,   "JMulTi" }
    };
    int i, nt = sizeof ftypes / sizeof ftypes[0];
    const char *p, *src = NULL;
    gchar *basename = NULL;
    char note[NOTELEN], tstr[48];

    for (i=0; i<nt; i++) {
	if (type == ftypes[i].type) {
	    src = ftypes[i].src;
	    break;
	}
    }

    p = strrchr(fname, SLASH);
    if (p != NULL) {
	basename = g_strdup(p + 1);
    } 

    print_time(tstr); 

    if (src != NULL) {
	snprintf(note, NOTELEN-1, "Data imported from %s file '%s', %s\n",
		 src, (basename == NULL)? fname : basename, tstr);
    } else {
	snprintf(note, NOTELEN-1, "Data imported from '%s', %s\n",
		 (basename == NULL)? fname : basename, tstr);
    }

    if (pdinfo->descrip == NULL) {
	pdinfo->descrip = gretl_strdup(note);
    } else {
	int dlen = strlen(pdinfo->descrip);
	int nlen = strlen(note);
	char *tmp = realloc(pdinfo->descrip, dlen + nlen + 3);

	if (tmp != NULL) {
	    pdinfo->descrip = tmp;
	    strcat(pdinfo->descrip, "\n\n");
	    strncat(pdinfo->descrip, note, strlen(note));
	}
    }

    g_free(basename);
}

