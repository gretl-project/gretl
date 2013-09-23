/*
   Reader for JMulTi data files.
   Allin Cottrell, Aigust 2006
*/

#include <stdio.h>
#include <stdlib.h>

#include "libgretl.h"
#include "version.h"

static int try_fix_varname (char *name)
{
    char test[VNAMELEN];
    int err = 0;

    *test = 0;
    strncat(test, name, VNAMELEN - 2);
    strcat(test, "1");
    
    err = check_varname(test);
    if (!err) {
	fprintf(stderr, "Warning: illegal name '%s' changed to '%s'\n",
		name, test);
	strcpy(name, test);
    } else {
	/* get the right error message in place */
	check_varname(name);
    }

    return err;
}

static int get_max_line_length (FILE *fp)
{
    int c, cc = 0;
    int maxlen = 0;

    while ((c = fgetc(fp)) != EOF) {
	if (c == '\n') {
	    if (cc > maxlen) {
		maxlen = cc;
	    }
	    cc = 0;
	}
	cc++;
    }

    if (maxlen > 0) {
	/* allow for newline and null terminator */
	maxlen += 2;
    }

    return maxlen;
}

static int push_desc (char *line, char **desc)
{
    char *p = strstr(line, "/*");
    char *q = strstr(line, "*/");
    char *s = NULL;
    int len, err = 0;

    if (p != NULL && q != NULL && q - p < 0) {
	return E_DATA;
    }

    if (p == NULL) {
	p = line;
    } else {
	p += 2;
    }

    if (q == NULL) {
	q = p + strlen(p);
    }

    len = q - p;

    if (*desc == NULL) {
	s = malloc(len + 2);
	if (s != NULL) {
	    *s = 0;
	}
    } else {
	s = realloc(*desc, strlen(*desc) + len + 2);
    }

    if (s == NULL) {
	free(*desc);
	*desc = NULL;
	err = E_ALLOC;
    } else {
	strncat(s, p, len);
	strcat(s, "\n");
	*desc = s;
    }
	
    return err;
}

static char *get_jmulti_comment (FILE *fp, char *line, int len,
				 int *err)
{
    char *desc = NULL;
    int outcomm = 0;
    int incomm = 0;

    while (fgets(line, len, fp) && !*err) {
	if (!incomm && strstr(line, "/*") != NULL) {
	    incomm = 1;
	}
	if (incomm) {
	    tailstrip(line);
	    *err = push_desc(line, &desc);
	    if (strstr(line, "*/") != NULL) {
		outcomm = 1;
		break;
	    }
	}
    } 

    if (incomm && !outcomm) {
	if (desc != NULL) {
	    free(desc);
	    desc = NULL;
	}
	*err = E_DATA;
    } else if (!incomm) {
	rewind(fp);
    }

    return desc;
}

static int parse_obs_bits (const char *s1, const char *s2,
			   int *pd, char *stobs)
{
    int y0, p0;
    int err = 0;

    if ((*s2 == 'Q' || *s2 == 'M') && strlen(s2) < 2) {
	return E_DATA;
    }

    if (strchr(s1, '/')) {
	/* dated daily */
	int y, m, d;
	int nf;

	if (strlen(s2) < 3) {
	    err = E_DATA;
	} else {
	    nf = sscanf(s1, "%d/%d/%d", &m, &d, &y); /* ?? */
	    if (nf != 3) {
		err = E_DATA;
	    } else {
		nf = sscanf(s2, "(%d)", pd);
		if (nf != 1) {
		   err = E_DATA;
		} 
	    } 
	}
	if (!err) {
	    sprintf(stobs, YMD_WRITE_FMT, y, m, d);
	    fprintf(stderr, "Daily data, starting '%s', pd = %d\n", stobs, *pd);
	}
    } else {
	y0 = atoi(s1);
	if (*s2 == 0) {
	    /* annual? */
	    sprintf(stobs, "%d", y0);
	    *pd = 1;
	} else if (*s2 == 'Q') {
	    /* quarterly */
	    p0 = atoi(s2 + 1);
	    sprintf(stobs, "%d:%d", y0, p0);
	    *pd = 4;
	} else if (*s2 == 'M') {
	    /* monthly */
	    p0 = atoi(s2 + 1);
	    sprintf(stobs, "%d:%02d", y0, p0);
	    *pd = 12;
	} else {
	    err = E_DATA;
	}
    }

    return err;
}

static int get_jmulti_obs_info (FILE *fp, char *line, int len,
				int *pd, char *stobs)
{
    char s1[16], s2[8];
    int nf, err = 0;

    while (fgets(line, len, fp) && !err) {
	if (*line == '<' && strchr(line, '>') != NULL) {
	    tailstrip(line);
	    fprintf(stderr, "jmulti obs line: %s\n", line);
	    if (!strncmp(line, "<1>", 3)) {
		/* undated data */
		*pd = 1;
		strcpy(stobs, "1");
	    } else {
		*s1 = *s2 = 0;
		nf = sscanf(line, "<%15s %7[^> ]", s1, s2);
		if (nf < 1) {
		    err = E_DATA;
		} else {
		    err = parse_obs_bits(s1, s2, pd, stobs);
		}
	    }
	    break;
	}
    } 

    return err;
}

static int get_jmulti_var_info (FILE *fp, char *line, int len,
				int *nv, int *nobs)
{
    long offset = ftell(fp);
    int err = 0;

    if (fgets(line, len, fp) == NULL) {
	return E_DATA;
    }

    tailstrip(line);
    fprintf(stderr, "jmulti var line: %s\n", line);

    *nv = count_fields(line, NULL);

    while (fgets(line, len, fp)) {
	if (!string_is_blank(line)) {
	    *nobs += 1;
	}
    }

    fprintf(stderr, "jmulti: %d variables, %d observations\n",
	    *nv, *nobs);

    fseek(fp, offset, SEEK_SET);

    return err;
}

static int read_jmulti_data (FILE *fp, char *line, int len,
			     DATASET *dset)
{
    char *s, numstr[32];
    int i, t, err = 0;

    if (fgets(line, len, fp) == NULL) {
	return E_DATA;
    }

    s = line;
    for (i=1; i<dset->v; i++) {
	s += strspn(s, " ");
	gretl_scan_varname(s, dset->varname[i]);
	if (check_varname(dset->varname[i]) && 
	    try_fix_varname(dset->varname[i])) {
	    err = 1;
	}
	s += strcspn(s, " \t\n");
    }

    gretl_push_c_numeric_locale();

    t = 0;
    while (fgets(line, len, fp) && !err) {
	if (string_is_blank(line)) {
	    continue;
	}
	s = line;
	for (i=1; i<dset->v && !err; i++) {
	    s += strspn(s, " \t");
	    if (sscanf(s, "%31s", numstr) != 1) {
		err = 1;
	    } else {
		if (!strcmp(numstr, "NaN")) {
		    dset->Z[i][t] = NADBL;
		} else {
		    dset->Z[i][t] = atof(numstr);
		}
		s += strcspn(s, " \t");
	    }
	}
	t++;
    } 

    gretl_pop_c_numeric_locale();

    return err;
}

int jmulti_get_data (const char *fname, DATASET *dset,
		     gretlopt opt, PRN *prn)
{
    int maxlen = 0;
    int nvar = 0, nobs = 0;
    FILE *fp;
    char *descrip = NULL;
    char *line = NULL;
    DATASET *newset = NULL;
    char stobs[OBSLEN];
    int pd = 0;
    int err = 0;

    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	return E_FOPEN;
    }

    maxlen = get_max_line_length(fp);

    if (maxlen == 0) {
	fclose(fp);
	return E_DATA;
    }

    line = malloc(maxlen);
    if (line == NULL) {
	fclose(fp);
	return E_ALLOC;
    }

    rewind(fp);

    descrip = get_jmulti_comment(fp, line, maxlen, &err);
    if (err) {
	free(line);
	fclose(fp);
	return err;
    }

    err = get_jmulti_obs_info(fp, line, maxlen, &pd, stobs);
    if (err) {
	goto bailout;
    }

    err = get_jmulti_var_info(fp, line, maxlen, &nvar, &nobs);
    if (err) {
	goto bailout;
    }

    newset = datainfo_new();
    if (newset == NULL) {
	pputs(prn, _("Out of memory\n"));
	err = E_ALLOC;
	goto bailout;
    }

    newset->v = nvar + 1;
    newset->n = nobs;
    newset->pd = pd;
    strcpy(newset->stobs, stobs);

    if (strcmp(newset->stobs, "1")) {
	newset->structure = TIME_SERIES;
	newset->sd0 = get_date_x(newset->pd, newset->stobs);
    }

    err = start_new_Z(newset, 0);
    if (err) {
	pputs(prn, _("Out of memory\n"));
	free_datainfo(newset);
	goto bailout;
    }	

    err = read_jmulti_data(fp, line, maxlen, newset);
    if (err) {
	destroy_dataset(newset);
    } else {
	int merge = (dset->Z != NULL);

	if (fix_varname_duplicates(newset)) {
	    pputs(prn, _("warning: some variable names were duplicated\n"));
	}

	newset->descrip = descrip;
	descrip = NULL;

	err = merge_or_replace_data(dset, &newset, opt, prn);

	if (!err && !merge) {
	    dataset_add_import_info(dset, fname, GRETL_JMULTI);
	}
    }

 bailout:

    fclose(fp);
    free(line);
    free(descrip);

    return err;
}  





