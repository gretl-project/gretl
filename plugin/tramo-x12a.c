/* TRAMO/SEATS, X-12-ARIMA plugin for gretl */

#include "libgretl.h"

#ifdef OS_WIN32
# include <windows.h>
#endif

typedef enum {
    D11,     /* seasonally adjusted series */
    D12,     /* trend/cycle */
    D13      /* irregular component */
} x12a_series;

const char *x12a_series_strings[] = {
    "d11", "d12", "d13"
};

const char *x12a_descrip_formats[] = {
    N_("seasonally adjusted %s"),
    N_("trend/cycle for %s"),
    N_("irregular component of %s")
};

int write_tramo_data (char *fname, int varnum, 
		      double ***pZ, DATAINFO *pdinfo, 
		      const char *tramodir)
{
    int i, t;
    int tsamp = pdinfo->t2 - pdinfo->t1 + 1;
    double x;
    char tmp[9], varname[9], cmd[MAXLEN];
    int startyr, startper;
    FILE *fp = NULL;

    *gretl_errmsg = 0;

    if (!pdinfo->vector[varnum]) {
	sprintf(gretl_errmsg, "%s %s", pdinfo->varname[varnum], 
		_("is a scalar"));
	return 1;
    }

    sprintf(varname, pdinfo->varname[varnum]);
    lower(varname);
    sprintf(fname, "%s%c%s", tramodir, SLASH, varname);

    fp = fopen(fname, "w");
    if (fp == NULL) return 1;

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "C");
#endif   

    x = date(pdinfo->t1, pdinfo->pd, pdinfo->sd0);
    startyr = (int) x;
    sprintf(tmp, "%g", x);
    startper = atoi(strchr(tmp, '.') + 1);

    fprintf(fp, "%s\n", pdinfo->varname[varnum]);
    fprintf(fp, "%d %d %d %d\n", tsamp, startyr, startper, pdinfo->pd);

    for (t=pdinfo->t1; t<=pdinfo->t2; ) {
	for (i=0; i<pdinfo->pd; i++) {
	    if (t == pdinfo->t1) {
		i += startper - 1;
	    }
	    if (na((*pZ)[varnum][t])) {
		fprintf(fp, "-99999 ");
	    } else {
		fprintf(fp, "%g ", (*pZ)[varnum][t]);
	    }
	    t++;
	}
	fputs("\n", fp);
    }

    /* FIXME: make these values configurable */
    fprintf(fp, "$INPUT lam=-1,iatip=1,aio=2,va=3.3,noadmiss=1,seats=2,$\n");

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif

    if (fp != NULL) {
	fclose(fp);
    }

    /* testing */
#ifdef notdef
    sprintf(cmd, "cd %s && ./tramo -i %s >/dev/null && mv seats.itr serie && ./seats -OF %s", 
	    tramodir, varname, varname);
#endif
#ifdef OS_WIN32 /* FIXME */
    sprintf(cmd, "tramo -i %s", 
	    tramodir, varname);
    WinExec(cmd, SW_SHOWMINIMIZED);
#else
    sprintf(cmd, "cd %s && ./tramo -i %s >/dev/null", 
	    tramodir, varname);
    system(cmd);
#endif

    strcpy(tmp, pdinfo->varname[varnum]);
    lower(tmp);
    sprintf(fname, "%s%coutput%c%s.out", tramodir, SLASH, SLASH, varname); /* .OUT for seats */

    return 0;
}

static void truncate (char *str, int n)
{
    int len = strlen(str);

    if (len > n) str[n] = 0;
}

static int graph_x12a_series (double **Z, DATAINFO *pdinfo, int varno)
{
    FILE *fp;
    int i, t;

    fp = fopen("x12a_gp.dat", "w");
    if (fp == NULL) return 1;

    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	fprintf(fp, "%g ", Z[varno][t]);
	for (i=1; i<4; i++) {
	    fprintf(fp, "%g ", Z[pdinfo->v - i][t]);
	}
	fprintf(fp, "\n");
    }
    fclose(fp);

    fp = fopen("x12a.gp", "w");
    if (fp == NULL) return 1;

    /* fixme tics */
    fprintf(fp, "set multiplot\n"
	    "set size 1.0,0.32\n"
	    "set origin 0.0,0.0\n"
	    "plot 'x12a_gp.dat' using 2 w l t 'irregular'\n"
	    "set origin 0.0,0.33\n"
	    "plot 'x12a_gp.dat' using 1 w l t '%s', \\\n"
	    " 'x12a_gp.dat' using 3 w l t 'trend/cycle'\n"
	    "set origin 0.0,0.66\n"
	    "plot 'x12a_gp.dat' using 1 w l t '%s', \\\n"
	    " 'x12a_gp.dat' using 4 w l t 'adjusted'\n"
	    "unset multiplot\n", 
	    pdinfo->varname[varno], pdinfo->varname[varno]);

    fclose(fp);
    
    system("gnuplot -persist x12a.gp");

    return 0;
}

static int add_x12a_series (const char *fname, int code,
			    double ***pZ, DATAINFO *pdinfo,
			    int varno)
{
    FILE *fp;
    char *p, line[128], varname[16], sfname[MAXLEN], date[8];
    double x;
    int d, yr, per, err = 0;
    int t, v = pdinfo->v;

    strcpy(sfname, fname);
    p = strrchr(sfname, '.');
    if (p != NULL) strcpy(p + 1, x12a_series_strings[code]);

    fp = fopen(sfname, "r");
    if (fp == NULL) {
	sprintf(gretl_errmsg, "%s %s", _("Couldn't open"), sfname);
	return 1;
    }

    /* formulate name of new variable to add */
    strcpy(varname, pdinfo->varname[varno]);
    truncate(varname, 4);
    strcat(varname, "_");
    strcat(varname, x12a_series_strings[code]);

    /* expand the dataset */
    if (dataset_add_vars(1, pZ, pdinfo)) {
	sprintf(gretl_errmsg, _("Out of memory adding data"));
	fclose(fp);
	return 1;
    }

    /* copy varname and label into place */
    strcpy(pdinfo->varname[v], varname);
    sprintf(pdinfo->label[v], _(x12a_descrip_formats[code]), pdinfo->varname[varno]);

    for (t=0; t<pdinfo->n; t++) (*pZ)[v][t] = NADBL;

    /* grab the data from the x12arima file */
    while (fgets(line, 127, fp)) {
	if (*line == 'd' || *line == '-') continue;
	if (sscanf(line, "%d %lf", &d, &x) != 2) {
	    err = 1; 
	    break;
	}
	yr = d / 100;
	per = d % 100;
	sprintf(date, "%d.%d", yr, per);
	t = dateton(date, pdinfo);
	if (t < 0 || t >= pdinfo->n) {
	    err = 1;
	    break;
	}
	(*pZ)[v][t] = x;
    }

    fclose(fp);

    return err;
}

int write_x12a_data (char *fname, int varnum, 
		     double ***pZ, DATAINFO *pdinfo, 
		     const char *x12adir)
{
    int i, t, err = 0;
    char tmp[8], varname[9], cmd[MAXLEN];
    int startyr, startper;
    int *savelist = NULL; /* FIXME configuration */
    double x;
    FILE *fp = NULL;

    *gretl_errmsg = 0;

    if (!pdinfo->vector[varnum]) {
	sprintf(gretl_errmsg, "%s %s", pdinfo->varname[varnum], 
		_("is a scalar"));
	return 1;
    }

    /* make a default x12a.mdl file if it doesn't already exist */
    sprintf(fname, "%s%cx12a.mdl", x12adir, SLASH);
    fp = fopen(fname, "r");
    if (fp == NULL) {
	fp = fopen(fname, "w");
	if (fp == NULL) return 1;
	fprintf(fp, "(0 1 1)(0 1 1) X\n"
		"(0 1 2)(0 1 1) X\n"
		"(2 1 0)(0 1 1) X\n"
		"(0 2 2)(0 1 1) X\n"
		"(2 1 2)(0 1 1)\n");
	fclose(fp);
    } else {
	fclose(fp);
    }

    sprintf(varname, pdinfo->varname[varnum]);
    sprintf(fname, "%s%c%s.spc", x12adir, SLASH, varname);

    fp = fopen(fname, "w");
    if (fp == NULL) return 1;

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "C");
#endif 

    x = date(pdinfo->t1, pdinfo->pd, pdinfo->sd0);
    startyr = (int) x;
    sprintf(tmp, "%g", x);
    startper = atoi(strchr(tmp, '.') + 1);

    fprintf(fp, "series{\n period=%d\n title=\"%s\"\n", pdinfo->pd, varname);
    fprintf(fp, " start=%d.%d\n", startyr, startper);
    fprintf(fp, " data=(\n");

    i = 0;
    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	if (na((*pZ)[varnum][t])) {
	    fprintf(fp, "-99999 "); /* FIXME? */
	} else {
	    fprintf(fp, "%g ", (*pZ)[varnum][t]);
	}
	if ((i + 1) % 7 == 0) fputs("\n", fp);
	i++;
    }
    fputs(" )\n}\n", fp);

    /* FIXME: make these values configurable */
    fputs("automdl{}\nx11{", fp);
    
    if (savelist != NULL) {
	if (savelist[0] == 1) {
	    fprintf(fp, " save=%s ", x12a_series_strings[savelist[1]]); 
	} else {
	    fputs(" save=( ", fp);
	    for (i=1; i<=savelist[0]; i++) {
		fprintf(fp, "%s ", x12a_series_strings[savelist[i]]);
	    }
	    fputs(") ", fp);
	}
    }

    fputs("}\n", fp);

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif

    if (fp != NULL) {
	fclose(fp);
    }

    /* FIXME win32 */
    sprintf(cmd, "cd %s && x12a %s -r -p -q >/dev/null", x12adir, varname);
    system(cmd);

    sprintf(fname, "%s%c%s.out", x12adir, SLASH, varname); 

    if (savelist != NULL) {
	for (i=1; i<=savelist[0]; i++) {
	    err = add_x12a_series (fname, savelist[i], pZ, pdinfo, varnum);
	}
	/* testing */
	graph_x12a_series(*pZ, pdinfo, varnum);
    }

    return err;
}

