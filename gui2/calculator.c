/*
 *   Copyright (C) Allin Cottrell
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

/* calculator.c for gretl */

#define NTESTS 6
#define NLOOKUPS 5
#define NTESTENTRY 7
#define NLOOKUPENTRY 4

#include "gretl.h"
#include "dlgutils.h"
#include <ctype.h>

#ifdef G_OS_WIN32
# include <windows.h>
#endif

#if GTK_MAJOR_VERSION < 2
# define OLD_GTK
#endif

typedef struct _GretlChild GretlChild;

struct _GretlChild {
    GtkWidget *win;
    GtkWidget *vbox;
    GtkWidget *action_area;
    gpointer data;
};

typedef struct {
    GtkWidget *entry[NTESTENTRY];
    GtkWidget *book;
    GtkWidget *check;
    GtkWidget *graph;
} test_t;

typedef struct {
    GtkWidget *entry[NLOOKUPENTRY];
    GtkWidget *book;
} lookup_t;

enum {
    NORMAL_DIST,
    T_DIST,
    CHISQ_DIST,
    F_DIST,
    DW_DIST
};

enum {
    NORMAL_PVAL,
    T_PVAL,
    CHISQ_PVAL,
    F_PVAL,
    GAMMA_PVAL
};

/* ........................................................... */

static int printnum (char *dest, const char *s, int d) 
{
    char numstr[16];

    if (s == NULL || *s == '\0') {
	errbox(_("Incomplete entry for p-value"));
	return 1;
    }

    if (isalpha(*s)) {
	int v; 
	double xx;

	if (data_status && (v = varindex(datainfo, s)) < datainfo->v) {
	    xx = Z[v][0];
	    if (na(xx)) {
		sprintf(errtext, _("Data missing for variable '%s'"), s);
		errbox(errtext);
		return 1;
	    }
	    if (d) sprintf(numstr, "%f", xx);
	    else sprintf(numstr, "%d", (int) xx);
	} else {
	    sprintf(errtext, 
		    _("Unrecognized variable '%s' in p-value command"), s);
	    errbox(errtext);
	    return 1;
	}
    } else if (d != 0) {
	sprintf(numstr, "%f ", atof(s));
    } else {
	sprintf(numstr, "%d ", atoi(s));
    }

    strcat(dest, numstr);

    return 0;
}

/* ........................................................... */

static double getval (const char *s, PRN *prn, int pos)
     /* if pos != 0, value must be positive or it is invalid */ 
{
    double x = NADBL;

    if (s == NULL || *s == '\0') {
	errbox(_("Incomplete entry"));
    } else {
	if (check_atof(s)) {
	    errbox(get_gretl_errmsg());
	} else {
	    x = atof(s);
	}

	if (pos && !na(x) && x <= 0.0) {
	    errbox(_("Invalid entry"));
	    x = NADBL;
	} 
    }

    if (x == NADBL && prn != NULL) {
	gretl_print_destroy(prn);
    }

    return x;
}

/* ........................................................... */

static int getint (const char *s, PRN *prn) 
{
    int n;

    if (s == NULL || *s == '\0') {
	errbox(_("Incomplete entry for hypothesis test"));
	gretl_print_destroy(prn);
	return -1;
    }

    n = atoi(s);

    if (n <= 0) {
	errbox(_("Invalid entry for hypothesis test"));
	return -1;
    } else {
	return n;
    }
}

/* ........................................................... */

static void get_critical (GtkWidget *w, gpointer data)
{
    lookup_t **look = (lookup_t **) data;
    void *handle = NULL;
    void *funp = NULL;
    void (*norm_table)(PRN *, int) = NULL;
    void (*dw)(int, PRN *) = NULL;
    void (*tcrit)(int, PRN *, int) = NULL;
    void (*chicrit)(int, PRN *, int) = NULL;
    int i, n = -1, df = -1, err = 0;
    int winheight = 300;
    PRN *prn;

    i = gtk_notebook_get_current_page(GTK_NOTEBOOK(look[0]->book));

    if (bufopen(&prn)) {
	return;
    }	

    switch (i) {

    case NORMAL_DIST:
	funp = norm_table = gui_get_plugin_function("norm_lookup", &handle);
	break;

    case T_DIST:
	df = atoi(gtk_entry_get_text(GTK_ENTRY(look[i]->entry[0])));
	funp = tcrit = gui_get_plugin_function("t_lookup", &handle);
	break;

    case CHISQ_DIST:
	df = atoi(gtk_entry_get_text(GTK_ENTRY(look[i]->entry[0])));
	funp = chicrit = gui_get_plugin_function("chisq_lookup", &handle);
	break;

    case F_DIST:
	df = atoi(gtk_entry_get_text(GTK_ENTRY(look[i]->entry[0])));
	n = atoi(gtk_entry_get_text(GTK_ENTRY(look[i]->entry[1])));
	break;

    case DW_DIST:
	n = atoi(gtk_entry_get_text(GTK_ENTRY(look[i]->entry[0])));
	funp = dw = gui_get_plugin_function("dw_lookup", &handle);
	break;

    default:
	break;
    }

    /* sanity */
    if ((0 < i && i < 4 && df <= 0) || (i == 3 && n <= 0)) {
	errbox(_("Invalid degrees of freedom"));
	err = 1;
    }
    else if (i == 4 && n <= 0) {
	errbox(_("Invalid sample size"));
	err = 1;
    }
    else if (i != 3 && funp == NULL)  {
	err = 1;
    }

    if (!err) {
	switch (i) {

	case NORMAL_DIST:
	    (*norm_table)(prn, 1);
	    winheight = 340;
	    break;

	case T_DIST:
	    (*tcrit)(df, prn, 1);
	    winheight = 340;
	    break;

	case CHISQ_DIST:
	    (*chicrit)(df, prn, 1);
	    break;
	
	case F_DIST:
	    pprintf(prn, _("Approximate critical values of F(%d, %d)\n\n"),
		    df, n);
	    pprintf(prn, _(" 10%% in right tail %.2f\n"), f_crit_a(.10, df, n));
	    pprintf(prn, "  5%%               %.2f\n", f_crit_a(.05, df, n));	
	    pprintf(prn, "  1%%               %.2f\n", f_crit_a(.01, df, n));
	    break;

	case DW_DIST:
	    (*dw)(n, prn);
	    break;

	default:
	    break;
	}
    }

    if (handle != NULL) {
	close_plugin(handle);
    }

    if (err) {
	gretl_print_destroy(prn);
    } else {
	view_buffer(prn, 77, winheight, _("gretl: statistical table"), 
		    STAT_TABLE, NULL);
    }
}

/* ........................................................... */

static void get_pvalue (GtkWidget *w, gpointer data)
{
    lookup_t **pval = (lookup_t **) data;
    gint i, j;
    int df;
    double val, xx, zz;
    const gchar *tmp;
    gchar cmd[128];
    PRN *prn;

    i = gtk_notebook_get_current_page(GTK_NOTEBOOK(pval[0]->book));

    sprintf(cmd, "pvalue %d ", i + 1);
    
    switch (i) {

    case NORMAL_PVAL:
	tmp = gtk_entry_get_text(GTK_ENTRY(pval[i]->entry[0]));
	xx = getval(tmp, NULL, 0); /* value */
	if (na(xx)) return;
	tmp = gtk_entry_get_text(GTK_ENTRY(pval[i]->entry[1]));
	zz = getval(tmp, NULL, 0); /* mean */
	if (na(zz)) return;
	xx -= zz;
	tmp = gtk_entry_get_text(GTK_ENTRY(pval[i]->entry[2]));
	val = getval(tmp, NULL, 0); /* std. deviation */
	if (na(val)) return;
	if (val <= 0) {
	    errbox(_("Invalid standard deviation"));
	    return;
	}
	xx /= val; 
	sprintf(cmd, "pvalue 1 %g", xx);
	break;

    case T_PVAL: 
	tmp = gtk_entry_get_text(GTK_ENTRY(pval[i]->entry[0]));
	df = atoi(tmp);   /* df */
	if (df <= 0) {
	    errbox(_("Invalid degrees of freedom"));
	    return;
	}
	tmp = gtk_entry_get_text(GTK_ENTRY(pval[i]->entry[1]));
	xx = getval(tmp, NULL, 0); /* value */
	if (na(xx)) return;
	tmp = gtk_entry_get_text(GTK_ENTRY(pval[i]->entry[2]));
	zz = getval(tmp, NULL, 0); /* mean */
	if (na(zz)) return;
	xx -= zz;
	tmp = gtk_entry_get_text(GTK_ENTRY(pval[i]->entry[3]));
	val = getval(tmp, NULL, 0); /* std. deviation */
	if (na(val)) return;
	if (val <= 0) {
	    errbox(_("Invalid standard deviation"));
	    return;
	}
	xx /= val; 
	sprintf(cmd, "pvalue 2 %d %g", df, xx);
	break;

    case CHISQ_PVAL:
	tmp = gtk_entry_get_text(GTK_ENTRY(pval[i]->entry[0]));
	df = atoi(tmp);   /* df */
	if (df <= 0) {
	    errbox(_("Invalid degrees of freedom"));
	    return;
	}	
	if (printnum(cmd, tmp, 0)) return;
	tmp = gtk_entry_get_text(GTK_ENTRY(pval[i]->entry[1]));
	if (printnum(cmd, tmp, 1)) return;
	break;

    case F_PVAL:
	for (j=0; j<2; j++) {
	    tmp = gtk_entry_get_text(GTK_ENTRY(pval[i]->entry[j]));
	    df = atoi(tmp);
	    if (df <= 0) {
		errbox(_("Invalid degrees of freedom"));
		return;
	    }	    
	    if (printnum(cmd, tmp, 0)) return;
	}
	tmp = gtk_entry_get_text(GTK_ENTRY(pval[i]->entry[2]));
	if (printnum(cmd, tmp, 1)) return;
	break;

    case GAMMA_PVAL: 
	for (j=0; j<3; j++) {
	    tmp = gtk_entry_get_text(GTK_ENTRY(pval[i]->entry[j]));
	    if (printnum(cmd, tmp, 1)) return;
	}
	break;

    default:
	break;
    }

    if (bufopen(&prn)) return;

    batch_pvalue(cmd, (const double **) Z, datainfo, prn);
    view_buffer(prn, 78, 200, _("gretl: p-value"), PVALUE, NULL);

    return;
}

/* ........................................................... */

static void print_pv (PRN *prn, double p1, double p2)
{
    pprintf(prn, _("Two-tailed p-value = %.4g\n(one-tailed = %.4g)\n"),
	    p1, p2);
}

/* ...................................................................	*/

static double chi_crit (const int df)
     /* v. rough 99.9 percent critical value of chi-square */
{
    double x = 10.0;
    
    while (chisq(x, df) > .001) x += .5;
    return x;
}

/* ...................................................................	*/

static double f_crit (const int df1, const int df2)
     /* v. rough 99.9 percent critical value of F */
{
    double x = 2.0;

    while (fdist(x, df1, df2) > .001) x += .5;
    return x;
}

/* ........................................................... */

static void htest_graph (int dist, double x, int df1, int df2)
{
    double xx, prange, spike = 0.0;
    FILE *fp = NULL;

    if (gnuplot_init(PLOT_SAMPLING_DIST, &fp)) return;

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "C");
#endif

    fprintf(fp, "set key right top\n");
    if (df1) fprintf(fp, "df1=%.1f\n", (double) df1);
    if (df2) fprintf(fp, "df2=%.1f\n", (double) df2);

    if (dist == NORMAL_DIST || dist == T_DIST) {
	xx = fabs(x);
	prange = ((xx > 3.5)? xx + .5 : 3.5);
	spike = .25;
	fprintf(fp, "set xrange [%.3f:%.3f]\n", -prange, prange);
	fprintf(fp, "set yrange [0:.50]\n");
	fprintf(fp, "set xlabel \"%s\"\n", I_("Standard errors"));
    }

    if (dist == T_DIST || dist == F_DIST) {
	fprintf(fp, "Binv(p,q)=exp(lgamma(p+q)-lgamma(p)-lgamma(q))\n");
    }

    if (dist == CHISQ_DIST) {
	prange = chi_crit(df1);
	if (x > prange) prange = 1.1 * x;
	spike = 1.0/prange;
	fprintf(fp, "set xrange [0:%.3f]\n", prange);
	if (df1 > 69) {
	    fprintf(fp, "log2=log(2.0)\n");
	    fprintf(fp, "chi(x)=exp((0.5*df1-1.0)*log(x)-0.5*x-"
		    "lgamma(0.5*df1)-df1*0.5*log2)\n");
	} else 
	    fprintf(fp, "chi(x)=x**(0.5*df1-1.0)*exp(-0.5*x)/gamma(0.5*df1)"
		    "/2**(0.5*df1)\n");
    }

    if (dist == F_DIST) { 
	prange = f_crit(df1, df2);
	if (x > prange) prange = 1.1 * x;
	spike = 1.0/prange;
	fprintf(fp, "set xrange [0:%.3f]\n", prange);
	fprintf(fp, "f(x)=Binv(0.5*df1,0.5*df2)*(df1/df2)**(0.5*df1)"
		"*x**(0.5*df1-1.0)/(1.0+df1/df2*x)**(0.5*(df1+df2))\n");
    }

    fprintf(fp, "plot \\\n");

    if (dist == NORMAL_DIST) {
	fprintf(fp, "(1/(sqrt(2*pi))*exp(-(x)**2/2)) "
		"title '%s' w lines , \\\n",
		I_("Gaussian sampling distribution"));
    }
    else if (dist == T_DIST) {
	char tmp[64];

	sprintf(tmp, I_("t(%d) sampling distribution"), df1);
	fprintf(fp, "Binv(0.5*df1,0.5)/sqrt(df1)*(1.0+(x*x)/df1)"
		"**(-0.5*(df1+1.0)) "
		"title '%s' w lines , \\\n", tmp);
    }
    else if (dist == CHISQ_DIST) {
	char tmp[64];
	
	sprintf(tmp, I_("Chi-square(%d) sampling distribution"), df1);
	fprintf(fp, "chi(x) title '%s' w lines , \\\n", tmp);
    }
    else if (dist == F_DIST) {
	char tmp[64];

	sprintf(tmp, I_("F(%d, %d) sampling distribution"), df1, df2);
	fprintf(fp, "f(x) title '%s' w lines , \\\n", tmp);
    }

    fprintf(fp, "'-' using 1:2 title '%s' w impulses\n",
	    I_("test statistic"));
    fprintf(fp, "%f %f\n", x, spike);
    fprintf(fp, "e\n");

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif

    fclose(fp);

    if (gnuplot_make_graph()) {
	errbox(_("gnuplot command failed"));
    } else {
	register_graph();
    }
}

/* ........................................................... */

static void h_test (GtkWidget *w, gpointer data)
{
    test_t **test = (test_t **) data;
    int i, j, n1, n2, grf = 0;
    double x[5], sderr, ts, pv;
    const gchar *tmp;
    PRN *prn;

    i = gtk_notebook_get_current_page(GTK_NOTEBOOK(test[0]->book));

    if (bufopen(&prn)) return;

    if (GTK_TOGGLE_BUTTON(test[i]->graph)->active) grf = 1;

    x[4] = 0.0;

    switch (i) {

    case 0: /* mean */
	for (j=0; j<2; j++) {
	    /* get sample mean and std dev */
	    tmp = gtk_entry_get_text(GTK_ENTRY(test[i]->entry[j]));
	    x[j] = getval(tmp, prn, (j==1)? 1 : 0);
	    if (na(x[j])) return;
	}
	tmp = gtk_entry_get_text(GTK_ENTRY(test[i]->entry[2]));
	if ((n1 = getint(tmp, prn)) == -1) return;
	tmp = gtk_entry_get_text(GTK_ENTRY(test[i]->entry[3]));
	x[2] = getval(tmp, prn, 0);
	if (na(x[2])) return;
	sderr = x[1]/sqrt((double) n1);
	ts = (x[0] - x[2])/sderr;
	pprintf(prn, _("Null hypothesis: population mean = %g\n"), x[2]);
	pprintf(prn, _("Sample size: n = %d\n"), n1);
	pprintf(prn, _("Sample mean = %g, std. deviation = %g\n"), 
		x[0], x[1]);
	if (GTK_TOGGLE_BUTTON(test[i]->check)->active) {
	    pprintf(prn, _("Test statistic: z = (%g - %g)/%g = %g\n"), 
		    x[0], x[2], sderr, ts);
	    pv = normal_pvalue_2(ts);
	    print_pv(prn, pv, pv / 2.0);
	    if (grf) htest_graph(0, ts, 0, 0);
	} else {
	    pprintf(prn, _("Test statistic: t(%d) = (%g - %g)/%g = %g\n"), n1-1,
		    x[0], x[2], sderr, ts);
	    pv = t_pvalue_2(ts, n1 - 1);
	    print_pv(prn, pv, 0.5 * pv);
	    if (grf) htest_graph(1, ts, n1-1, 0);
	}
	break;

    case 1: /* variance */
	tmp = gtk_entry_get_text(GTK_ENTRY(test[i]->entry[0]));
	x[0] = getval(tmp, prn, 1);
	if (na(x[0])) return;
	tmp = gtk_entry_get_text(GTK_ENTRY(test[i]->entry[1]));
	if ((n1 = getint(tmp, prn)) == -1) return;
	tmp = gtk_entry_get_text(GTK_ENTRY(test[i]->entry[2]));
	x[1] = getval(tmp, prn, 1);
	if (na(x[1])) return;
	ts = (n1 - 1) * x[0] / x[1];
	pprintf(prn, _("Null hypothesis: population variance = %g\n"), x[1]);
	pprintf(prn, _("Sample size: n = %d\n"), n1);
	pprintf(prn, _("Sample variance = %g\n"), x[0]);
	pprintf(prn, _("Test statistic: chi-square(%d) = %d * %g/%g = %g\n"), 
		n1-1, n1-1, x[0], x[1], ts);
	if (x[0] > x[1])
	    pv = chisq(ts, n1-1);
	else
	    pv = 1.0 - chisq(ts, n1-1);
	print_pv(prn, 2.0 * pv, pv);
	if (grf) htest_graph(2, ts, n1-1, 0);
	break;

    case 2: /* proportion */
	tmp = gtk_entry_get_text(GTK_ENTRY(test[i]->entry[0]));
	x[0] = getval(tmp, prn, 1);
	if (na(x[0])) return;
	tmp = gtk_entry_get_text(GTK_ENTRY(test[i]->entry[1]));
	if ((n1 = getint(tmp, prn)) == -1) return;
	tmp = gtk_entry_get_text(GTK_ENTRY(test[i]->entry[2]));
	x[1] = getval(tmp, prn, 1);
	if (na(x[1])) return;
	if (n1 * x[1] < 5.0 || n1 * (1.0 - x[1]) < 5.0) {
	    infobox(_("The assumption of a normal sampling distribution\n"
		      "is not justified here.  Abandoning the test."));
	    gretl_print_destroy(prn);
	    return;
	}
	sderr = sqrt(x[1] * (1.0 - x[1]) / n1);
	ts = (x[0] - x[1]) / sderr;
	pprintf(prn, _("Null hypothesis: population proportion = %g\n"), x[1]);
	pprintf(prn, _("Sample size: n = %d\n"), n1);
	pprintf(prn, _("Sample proportion = %g\n"), x[0]);
	pprintf(prn, _("Test statistic: z = (%g - %g)/%g = %g\n"), 
		x[0], x[1], sderr, ts);
	pv = normal_pvalue_2(ts);
	print_pv(prn, pv, pv / 2.0);
	if (grf) htest_graph(0, ts, 0, 0);
	break;

    case 3: /* two means */
	for (j=0; j<2; j++) {
	    /* mean and std dev, sample 1 */
	    tmp = gtk_entry_get_text(GTK_ENTRY(test[i]->entry[j]));
	    x[j] = getval(tmp, prn, (j==1)? 1 : 0);
	    if (na(x[j])) return;
	}
	tmp = gtk_entry_get_text(GTK_ENTRY(test[i]->entry[2]));
	if ((n1 = getint(tmp, prn)) == -1) return;
	for (j=2; j<4; j++) {
	    /* mean and std dev, sample 2 */
	    tmp = gtk_entry_get_text(GTK_ENTRY(test[i]->entry[j+1]));
	    x[j] = getval(tmp, prn, (j==3)? 1 : 0);
	    if (na(x[j])) return;
	}
	tmp = gtk_entry_get_text(GTK_ENTRY(test[i]->entry[5]));
	if ((n2 = getint(tmp, prn)) == -1) return;
	tmp = gtk_entry_get_text(GTK_ENTRY(test[i]->entry[6]));
	x[4] = getval(tmp, prn, 0);
	if (na(x[4])) return;
	pprintf(prn, _("Null hypothesis: Difference of means = %g\n"), x[4]);
	pprintf(prn, _("Sample 1:\n n = %d, mean = %g, s.d. = %g\n"),
		n1, x[0], x[1]);
	pprintf(prn, _("Sample 2:\n n = %d, mean = %g, s.d. = %g\n"),
		n2, x[2], x[3]);
	/* are we assuming a common variance? */
	j = 0;
	if (GTK_TOGGLE_BUTTON(test[i]->check)->active) j = 1;
	if (n1 < 30 || n2 < 30) j = 2;
	if (j > 0) {
	    ts = ((n1-1)*x[1]*x[1] + (n2-1)*x[3]*x[3])/(n1 + n2 - 2);
	    sderr = sqrt(ts/n1 + ts/n2);
	} else {
	    sderr = sqrt(x[1]*x[1]/n1 + x[3]*x[3]/n2);
	}
	ts = (x[0] - x[2] - x[4]) / sderr;
	if (j) {
	    if (j == 2)
		pprintf(prn, _("Small samples: assuming normality and common "
			       "variance\n"));
	    pprintf(prn, _("Test statistic: t(%d) = (%g - %g)/%g = %g\n"),
		    n1+n2-2, x[0], x[2], sderr, ts);
	    if (ts > 0) {
		pv = t_pvalue_2(ts, n1+n2-2);
	    } else {
		pv = t_pvalue_2(-ts, n1+n2-2);
	    }
	    print_pv(prn, pv, 0.5 * pv);
	    if (grf) htest_graph(1, ts, n1+n2-2, 0);
	} else {
	    if (x[4] > 0.0) {
		pprintf(prn, _("Test statistic: z = (%g - %g - %g)/%g = %g\n"),
			x[0], x[2], x[4], sderr, ts);
	    } else if (x[4] < 0.0) {
		pprintf(prn, _("Test statistic: z = [(%g - %g) - (%g)]/%g = %g\n"),
			x[0], x[2], x[4], sderr, ts);
	    } else {
		pprintf(prn, _("Test statistic: z = (%g - %g)/%g = %g\n"),
			x[0], x[2], sderr, ts);
	    }
	    pv = normal_pvalue_2(ts);
	    print_pv(prn, pv, pv / 2.0);
	    if (grf) htest_graph(0, ts, 0, 0);
	}
	break;

    case 4: /* two variances */
	tmp = gtk_entry_get_text(GTK_ENTRY(test[i]->entry[0]));
	x[0] = getval(tmp, prn, 1);
	if (na(x[0])) return;
	tmp = gtk_entry_get_text(GTK_ENTRY(test[i]->entry[1]));
	if ((n1 = getint(tmp, prn)) == -1) return;
	tmp = gtk_entry_get_text(GTK_ENTRY(test[i]->entry[2]));
	x[1] = getval(tmp, prn, 1);
	if (na(x[1])) return;
	tmp = gtk_entry_get_text(GTK_ENTRY(test[i]->entry[3]));
	if ((n2 = getint(tmp, prn)) == -1) return;
	pprintf(prn, _("Null hypothesis: The population variances are "
		       "equal\n"));
	pprintf(prn, _("Sample 1:\n n = %d, variance = %g\n"), n1, x[0]);
	pprintf(prn, _("Sample 2:\n n = %d, variance = %g\n"), n2, x[1]);
	if (x[0] > x[1]) {
	    ts = x[0]/x[1];
	    pprintf(prn, _("Test statistic: F(%d, %d) = %g\n"), 
		    n1-1, n2-1, ts);
	    pv = fdist(ts, n1-1, n2-1);
	} else {
	    ts = x[1]/x[0];
	    pprintf(prn, _("Test statistic: F(%d, %d) = %g\n"), 
		    n2-1, n1-1, ts);
	    pv = fdist(ts, n2-1, n1-1);
	}
	print_pv(prn, 2.0 * pv, pv);
	if (grf) htest_graph(3, ts, n1-1, n2-1);
	break;

    case 5: /* two proportions */
	tmp = gtk_entry_get_text(GTK_ENTRY(test[i]->entry[0]));
	x[0] = getval(tmp, prn, 1);
	if (na(x[0])) return;

	tmp = gtk_entry_get_text(GTK_ENTRY(test[i]->entry[1]));
	if ((n1 = getint(tmp, prn)) == -1) return;
	tmp = gtk_entry_get_text(GTK_ENTRY(test[i]->entry[2]));
	x[1] = getval(tmp, prn, 1);
	if (na(x[1])) return;

	tmp = gtk_entry_get_text(GTK_ENTRY(test[i]->entry[3]));
	if ((n2 = getint(tmp, prn)) == -1) return;

	pprintf(prn, _("Null hypothesis: the population proportions are "
		       "equal\n"));
	pprintf(prn, _("Sample 1:\n n = %d, proportion = %g\n"), n1, x[0]);
	pprintf(prn, _("Sample 2:\n n = %d, proportion = %g\n"), n2, x[1]);
	x[2] = (n1*x[0] + n2*x[1]) / (n1 + n2);
	sderr = sqrt((x[2] * (1.0-x[2])) * (1.0/n1 + 1.0/n2));
	ts = (x[0] - x[1]) / sderr;
	pprintf(prn, _("Test statistic: z = (%g - %g) / %g = %g\n"),
		x[0], x[1], sderr, ts);
	pv = normal_pvalue_2(ts);
	print_pv(prn, pv, pv / 2.0);
	if (grf) htest_graph(0, ts, 0, 0);
	break;

    default:
	break;
    }

    view_buffer(prn, 78, 300, _("gretl: hypothesis test"), H_TEST,
                NULL);
}

/* ........................................................... */

static void add_lookup_entry (GtkWidget *tbl, gint *tbl_len, 
			      const gchar *label, lookup_t **look, 
			      int code, int pval)
{
    GtkWidget *tempwid;

    *tbl_len += 1;

    /* label */
    gtk_table_resize(GTK_TABLE (tbl), *tbl_len, 2);
    tempwid = gtk_label_new(_(label));
    gtk_misc_set_alignment (GTK_MISC (tempwid), 1, 0.5);
    gtk_table_attach_defaults(GTK_TABLE (tbl), 
			      tempwid, 0, 1, *tbl_len - 1, *tbl_len);
    gtk_widget_show(tempwid);

    /* numeric entry */
    tempwid = gtk_entry_new();
    gtk_table_attach_defaults(GTK_TABLE (tbl), 
			      tempwid, 1, 2, *tbl_len - 1, *tbl_len);
    gtk_widget_show (tempwid);
    look[code]->entry[*tbl_len - 2] = tempwid;

    g_signal_connect(G_OBJECT(tempwid), "activate", 
		     (pval)? G_CALLBACK(get_pvalue) : 
		     G_CALLBACK(get_critical),
		     look);
}

/* .................................................................. */

static void make_lookup_tab (GtkWidget *notebook, int code, lookup_t **look) 
{
    GtkWidget *tempwid, *box, *tbl;
    gint tbl_len;
    const gchar *titles[] = {
	N_("normal"), 
	N_(" t "), 
	N_("chi-square"), 
	N_(" F "), 
	N_(" DW ")
    };
   
    box = gtk_vbox_new (FALSE, 0);
    gtk_container_set_border_width (GTK_CONTAINER (box), 10);
    gtk_widget_show (box);

    tempwid = gtk_label_new (_(titles[code]));
    gtk_widget_show (tempwid);
    gtk_notebook_append_page (GTK_NOTEBOOK (notebook), box, tempwid);   

    tbl_len = 1;
    tbl = gtk_table_new (tbl_len, 2, FALSE);
    gtk_table_set_row_spacings (GTK_TABLE (tbl), 5);
    gtk_table_set_col_spacings (GTK_TABLE (tbl), 5);
    gtk_box_pack_start (GTK_BOX (box), tbl, FALSE, FALSE, 0);
    gtk_widget_show (tbl);
   
    switch (code) {

    case NORMAL_DIST:
	gtk_table_resize(GTK_TABLE (tbl), ++tbl_len, 1);
	tempwid = gtk_label_new(_("Critical values for\n"
				  "standard normal distribution"));
	gtk_table_attach_defaults(GTK_TABLE (tbl), 
				  tempwid, 0, 2, tbl_len - 1, tbl_len);
	gtk_widget_show(tempwid);
	break;

    case T_DIST:
	add_lookup_entry(tbl, &tbl_len, N_("df"), look, code, 0);
	break;

    case CHISQ_DIST:
	add_lookup_entry(tbl, &tbl_len, N_("df"), look, code, 0);
	break;

    case F_DIST:
	add_lookup_entry(tbl, &tbl_len, N_("dfn"), look, code, 0);
	add_lookup_entry(tbl, &tbl_len, N_("dfd"), look, code, 0);
	break;	

    case DW_DIST:
	add_lookup_entry(tbl, &tbl_len, N_("n"), look, code, 0);
	break;

    default:
	break;
    } 
}

/* .................................................................. */

static void make_dist_tab (GtkWidget *notebook, int code, lookup_t **pval) 
{
    GtkWidget *tempwid, *box, *tbl;
    gint tbl_len;
    const gchar *titles[] = {
	N_("normal"), 
	N_(" t "), 
	N_("chi-square"), 
	N_(" F "), 
	N_("gamma")
    };
   
    box = gtk_vbox_new (FALSE, 0);
    gtk_container_set_border_width (GTK_CONTAINER (box), 10);
    gtk_widget_show (box);

    tempwid = gtk_label_new (_(titles[code]));
    gtk_widget_show (tempwid);
    gtk_notebook_append_page (GTK_NOTEBOOK (notebook), box, tempwid);   

    tbl_len = 1;
    tbl = gtk_table_new (tbl_len, 2, FALSE);
    gtk_table_set_row_spacings (GTK_TABLE (tbl), 5);
    gtk_table_set_col_spacings (GTK_TABLE (tbl), 5);
    gtk_box_pack_start (GTK_BOX (box), tbl, FALSE, FALSE, 0);
    gtk_widget_show (tbl);
   
    switch (code) {

    case NORMAL_PVAL: 
	add_lookup_entry(tbl, &tbl_len, N_("value"), pval, code, 1);
	add_lookup_entry(tbl, &tbl_len, N_("mean"), pval, code, 1);
	gtk_entry_set_text(GTK_ENTRY(pval[0]->entry[1]), "0");
	add_lookup_entry(tbl, &tbl_len, N_("std. deviation"), pval, code, 1);
	gtk_entry_set_text(GTK_ENTRY(pval[0]->entry[2]), "1");
	break;

    case T_PVAL:
	add_lookup_entry(tbl, &tbl_len, N_("df"), pval, code, 1);
	add_lookup_entry(tbl, &tbl_len, N_("value"), pval, code, 1);
	add_lookup_entry(tbl, &tbl_len, N_("mean"), pval, code, 1);
	gtk_entry_set_text(GTK_ENTRY(pval[1]->entry[2]), "0");
	add_lookup_entry(tbl, &tbl_len, N_("std. deviation"), pval, code, 1);
	gtk_entry_set_text(GTK_ENTRY(pval[1]->entry[3]), "1");
	break;

    case CHISQ_PVAL:
	add_lookup_entry(tbl, &tbl_len, N_("df"), pval, code, 1);
	add_lookup_entry(tbl, &tbl_len, N_("x-value"), pval, code, 1);
	break;

    case F_PVAL:
	add_lookup_entry(tbl, &tbl_len, N_("dfn"), pval, code, 1);
	add_lookup_entry(tbl, &tbl_len, N_("dfd"), pval, code, 1);
	add_lookup_entry(tbl, &tbl_len, N_("F-value"), pval, code, 1);
	break;

    case GAMMA_PVAL:
	add_lookup_entry(tbl, &tbl_len, N_("mean"), pval, code, 1);
	add_lookup_entry(tbl, &tbl_len, N_("variance"), pval, code, 1);
	add_lookup_entry(tbl, &tbl_len, N_("x-value"), pval, code, 1);
	break;

    default:
	break;
    } 
}

/* ........................................................... */
 
static void trash_look (GtkWidget *w, gpointer data)
{
    lookup_t **look = (lookup_t **) data;
    int i;

    for (i=0; i<NLOOKUPS; i++) free(look[i]);
    free(look);
}

/* ........................................................... */
 
static void trash_test (GtkWidget *w, gpointer data)
{
    test_t **test = (test_t **) data;
    int i;

    for (i=0; i<NTESTS; i++) free(test[i]);
    free(test);
}

/* ........................................................... */

static void add_test_entry (GtkWidget *tbl, gint *tbl_len, 
			    const gchar *label, test_t **test, int code)
{
    GtkWidget *tempwid;

    *tbl_len += 1;
    gtk_table_resize (GTK_TABLE (tbl), *tbl_len, 2);
    tempwid = gtk_label_new (label);
    gtk_misc_set_alignment (GTK_MISC (tempwid), 1, 0.5);
    gtk_table_attach_defaults (GTK_TABLE (tbl), 
			       tempwid, 0, 1, *tbl_len - 1, *tbl_len);
    gtk_widget_show (tempwid);
    tempwid = gtk_entry_new ();
    gtk_table_attach_defaults (GTK_TABLE (tbl), 
			       tempwid, 1, 2, *tbl_len - 1, *tbl_len);
    gtk_widget_show (tempwid);
    test[code]->entry[*tbl_len - 2] = tempwid;

    g_signal_connect(G_OBJECT(tempwid), "activate", 
		     G_CALLBACK(h_test), test);
}

/* ........................................................... */

static void add_test_label (GtkWidget *tbl, gint *tbl_len, 
			    const gchar *label)
{
    GtkWidget *tempwid;

    *tbl_len += 1;
    gtk_table_resize (GTK_TABLE (tbl), *tbl_len, 2);
    tempwid = gtk_label_new (label);
    gtk_misc_set_alignment (GTK_MISC (tempwid), 0, 0.5);
    gtk_table_attach_defaults (GTK_TABLE (tbl), 
			       tempwid, 0, 2, *tbl_len - 1, *tbl_len);
    gtk_widget_show (tempwid);
}

/* ........................................................... */

static void add_test_check (GtkWidget *tbl, gint *tbl_len, 
			    const gchar *label, test_t **test, int code)
{
    GtkWidget *tempwid;

    *tbl_len += 1;
    gtk_table_resize (GTK_TABLE (tbl), *tbl_len, 2);
    tempwid = gtk_check_button_new_with_label (label);
    gtk_table_attach_defaults (GTK_TABLE (tbl), 
			       tempwid, 0, 2, *tbl_len - 1, *tbl_len);
    gtk_widget_show (tempwid);
    test[code]->check = tempwid;
}

/* .................................................................. */

static void make_test_tab (GtkWidget *notebook, int code, test_t **test) 
{
    GtkWidget *tempwid, *box, *tbl;
    gint tbl_len;
    const gchar *titles[] = {
	N_("mean"), 
	N_("variance"), 
	N_("proportion"),
	N_("2 means"), 
	N_("2 variances"), 
	N_("2 proportions")
    };
   
    box = gtk_vbox_new (FALSE, 0);
    gtk_container_set_border_width (GTK_CONTAINER (box), 10);
    gtk_widget_show (box);

    tempwid = gtk_label_new (_(titles[code]));
    gtk_widget_show (tempwid);
    gtk_notebook_append_page (GTK_NOTEBOOK (notebook), box, tempwid);   

    tbl_len = 1;
    tbl = gtk_table_new (tbl_len, 2, FALSE);
    gtk_table_set_row_spacings (GTK_TABLE (tbl), 5);
    gtk_table_set_col_spacings (GTK_TABLE (tbl), 5);
    gtk_box_pack_start (GTK_BOX (box), tbl, FALSE, FALSE, 0);
    gtk_widget_show (tbl);
   
    switch (code) {

    case 0: /* mean */
	add_test_entry(tbl, &tbl_len, _("sample mean"), test, code);
	add_test_entry(tbl, &tbl_len, _("std. deviation"), test, code);
	add_test_entry(tbl, &tbl_len, _("sample size"), test, code);
	add_test_entry(tbl, &tbl_len, _("H0: mean ="), test, code);
	add_test_check(tbl, &tbl_len, _("Assume standard deviation is "
		       "population value"), test, code);
	break;

    case 1: /* variance */
	add_test_entry(tbl, &tbl_len, _("sample variance"), test, code);
	add_test_entry(tbl, &tbl_len, _("sample size"), test, code);
	add_test_entry(tbl, &tbl_len, _("H0: variance ="), test, code);
	break;

    case 2: /* proportion */
	add_test_entry(tbl, &tbl_len, _("sample proportion"), test, code);
	add_test_entry(tbl, &tbl_len, _("sample size"), test, code);
	add_test_entry(tbl, &tbl_len, _("H0: proportion ="), test, code);
	break;

    case 3: /* two means */
	add_test_entry(tbl, &tbl_len, _("mean of sample 1"), test, code);
	add_test_entry(tbl, &tbl_len, _("std. deviation, sample 1"), test, code);
	add_test_entry(tbl, &tbl_len, _("size of sample 1"), test, code);
	add_test_entry(tbl, &tbl_len, _("mean of sample 2"), test, code);
	add_test_entry(tbl, &tbl_len, _("std. deviation, sample 2"), test, code);
	add_test_entry(tbl, &tbl_len, _("size of sample 2"), test, code);
	add_test_entry(tbl, &tbl_len, _("H0: Difference of means ="), test, code);
	gtk_entry_set_text(GTK_ENTRY(test[3]->entry[6]), "0");
	add_test_check(tbl, &tbl_len, _("Assume common population standard "
		       "deviation"), test, code);
	break;

    case 4: /* two variances */
	add_test_entry(tbl, &tbl_len, _("variance of sample 1"), test, code);
	add_test_entry(tbl, &tbl_len, _("size of sample 1"), test, code);
	add_test_entry(tbl, &tbl_len, _("variance of sample 2"), test, code);
	add_test_entry(tbl, &tbl_len, _("size of sample 2"), test, code);
	add_test_label(tbl, &tbl_len, _("H0: Ratio of variances = 1"));
	break;

    case 5: /* two proportions */
	add_test_entry(tbl, &tbl_len, _("proportion, sample 1"), test, code);
	add_test_entry(tbl, &tbl_len, _("size of sample 1"), test, code);
	add_test_entry(tbl, &tbl_len, _("proportion, sample 2"), test, code);
	add_test_entry(tbl, &tbl_len, _("size of sample 2"), test, code);
	add_test_label(tbl, &tbl_len, _("H0: Difference of proportions = 0"));
	break;

    default:
	break;
    } 

    /* add check box for showing graph of sampling dist. */
    tbl_len += 1;
    gtk_table_resize(GTK_TABLE(tbl), tbl_len, 2);
    tempwid = gtk_check_button_new_with_label(_("Show graph of sampling "
						"distribution"));
    gtk_table_attach_defaults(GTK_TABLE(tbl), tempwid, 0, 2, 
			      tbl_len - 1, tbl_len);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON (tempwid), TRUE);
    gtk_widget_show (tempwid);

    test[code]->graph = tempwid; 
}

static void gretl_child_destroy (GtkWidget *w, gpointer data)
{
    GretlChild *gchild = (GretlChild *) data;
    GtkWidget **wp = (GtkWidget **) gchild->data;

    *wp = NULL;
    free(gchild);
}

static GretlChild *gretl_child_new (const gchar *title)
{
    GretlChild *gchild;
    GtkWidget *base, *hsep;

    gchild = mymalloc(sizeof *gchild);
    if (gchild == NULL) return NULL;

    gchild->win = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    if (title != NULL)
	gtk_window_set_title(GTK_WINDOW(gchild->win), title);

    base = gtk_vbox_new(FALSE, 0);
    gtk_container_add(GTK_CONTAINER(gchild->win), base);
    gtk_widget_show(base);

    gchild->vbox = gtk_vbox_new(FALSE, 5);
    gtk_widget_show(gchild->vbox);
    gtk_container_add(GTK_CONTAINER(base), gchild->vbox);
    gtk_container_set_border_width(GTK_CONTAINER(gchild->vbox), 5);

    hsep = gtk_hseparator_new();
    gtk_widget_show(hsep);
    gtk_container_add(GTK_CONTAINER(base), hsep);
    
    gchild->action_area = gtk_hbox_new(TRUE, 5);
    gtk_widget_show(gchild->action_area);
    gtk_container_add(GTK_CONTAINER(base), gchild->action_area);
    gtk_container_set_border_width(GTK_CONTAINER(gchild->action_area), 5);

    gtk_window_set_position(GTK_WINDOW(gchild->win), GTK_WIN_POS_MOUSE);

    g_signal_connect(G_OBJECT(gchild->win), "destroy",
                     G_CALLBACK(gretl_child_destroy),
                     gchild);

    return gchild;
}

/* ........................................................... */

void stats_calculator (gpointer data, guint code, GtkWidget *widget) 
{
    GtkWidget *tempwid = NULL, *notebook;
    static GtkWidget *winptr[3];
    GtkWidget *thiswin;
    GretlChild *dialog;
    test_t **test = NULL;
    lookup_t **look = NULL;
    const char *window_titles[] = {
	N_("gretl: p-value finder"),
	N_("gretl: statistical tables"),
	N_("gretl: test calculator")
    };
    gpointer statp = NULL;
    gint i;

    g_return_if_fail(code == CALC_PVAL || 
		     code == CALC_DIST || 
		     code == CALC_TEST);

    thiswin = winptr[code];

    if (thiswin != NULL) {
        gdk_window_show(thiswin->window);
        gdk_window_raise(thiswin->window);
 	return;
    }

    if (code == CALC_TEST) {
	test = mymalloc(NTESTS * sizeof *test);
	if (test == NULL) return;
	statp = test;
    } else {
	look = mymalloc(NLOOKUPS * sizeof *look);
	if (look == NULL) return;
	statp = look;
    }	

    dialog = gretl_child_new(NULL);
    winptr[code] = dialog->win;
    dialog->data = &(winptr[code]);

    gtk_window_set_title(GTK_WINDOW(dialog->win), _(window_titles[code]));

    notebook = gtk_notebook_new ();
    gtk_box_pack_start(GTK_BOX(dialog->vbox), notebook, TRUE, TRUE, 0);
    gtk_widget_show (notebook);

    if (code == CALC_TEST) {
	for (i=0; i<NTESTS; i++) {
	    test[i] = mymalloc(sizeof **test);
	    if (test[i] == NULL) return;
	    test[i]->book = notebook;
	    make_test_tab(notebook, i, test);
	}
    } else {
	for (i=0; i<NLOOKUPS; i++) {
	    look[i] = mymalloc(sizeof **look);
	    if (look[i] == NULL) return;
	    look[i]->book = notebook;
	    if (code == 0) {
		make_dist_tab(notebook, i, look);
	    } else {
		make_lookup_tab(notebook, i, look);
	    }
	}
    }	

    tempwid = standard_button(GTK_STOCK_OK);
    GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX(dialog->action_area), 
			tempwid, TRUE, TRUE, 0);

    g_signal_connect (G_OBJECT (tempwid), "clicked", 
		      (code == CALC_PVAL)? G_CALLBACK(get_pvalue) :
		      (code == CALC_DIST)? G_CALLBACK(get_critical) :
		      G_CALLBACK(h_test),
		      statp);

    gtk_widget_show (tempwid);

    /* Close button */
    tempwid = standard_button(GTK_STOCK_CLOSE);
    GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX(dialog->action_area), 
			tempwid, TRUE, TRUE, 0);

    g_signal_connect(G_OBJECT(tempwid), "clicked", 
		     (code == CALC_TEST)? G_CALLBACK(trash_test) :
		     G_CALLBACK(trash_look),
		     statp);
    g_signal_connect(G_OBJECT(tempwid), "clicked", 
		     G_CALLBACK(delete_widget), 
		     dialog->win);

    gtk_widget_show(tempwid);

    gtk_widget_show(dialog->win);
}





