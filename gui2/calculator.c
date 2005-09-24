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
#define NPVAL 6
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

typedef struct GretlChild_ GretlChild;
typedef struct test_t_ test_t;
typedef struct lookup_t_ lookup_t;

struct GretlChild_ {
    GtkWidget *win;
    GtkWidget *vbox;
    GtkWidget *action_area;
    gpointer data;
};

struct test_t_ {
    int code;
    GtkWidget *entry[NTESTENTRY];
    GtkWidget *combo[2];
    GtkWidget *book;
    GtkWidget *check;
    GtkWidget *graph;
};

struct lookup_t_ {
    GtkWidget *entry[NLOOKUPENTRY];
    GtkWidget *book;
};

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
    GAMMA_PVAL,
    BINOMIAL_PVAL
};

enum {
    ONE_MEAN,
    ONE_VARIANCE,
    ONE_PROPN,
    TWO_MEANS,
    TWO_VARIANCES,
    TWO_PROPNS
};

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
	sprintf(numstr, "%g ", atof(s));
    } else {
	sprintf(numstr, "%d ", atoi(s));
    }

    strcat(dest, numstr);

    return 0;
}

/* if pos != 0, value must be positive or it is invalid */ 

static double getval (const char *s, PRN *prn, int pos)
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
    } else if (i == 4 && n <= 0) {
	errbox(_("Invalid sample size"));
	err = 1;
    } else if (i != 3 && funp == NULL)  {
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

static int do_binomial_pdf (const char *s, PRN *prn)
{
    int k, n;
    double p, x1, x2;
    int err = 0;

    if (sscanf(s, "pvalue 6 %lf %d %d", &p, &n, &k) != 3) {
	err = 1;
    } else if (p < 0.0 || p > 1.0) {
	err = 1;
    } else if (n <= 0) {
	err = 1;
    } else if (k < 0 || k > n) {
	err = 1;
    }

    if (err) {
	errbox(_("Invalid entry"));
	return err;
    }

    x1 = binomial_cdf(k, n, p);
    x2 = binomial_pvalue(k, n, p);

    if (na(x1) || na(x2)) {
	err = 1;
    } else {
	pprintf(prn, _("Binomial distribution: p = %g, %d trials"), p, n);
	pputc(prn, '\n');
	pputs(prn, _("x = number of 'successes'"));
	pputs(prn, "\n\n");
	pprintf(prn, " P(x <= %d) = %g\n", k, x1);
	pprintf(prn, " P(x > %d)  = %g\n", k, x2);
	pprintf(prn, " P(x = %d)  = %g", k, binomial_cdf(k+1, n, p) - x1);
    }

    return err;
}

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

    case BINOMIAL_PVAL: 
	for (j=0; j<3; j++) {
	    tmp = gtk_entry_get_text(GTK_ENTRY(pval[i]->entry[j]));
	    if (printnum(cmd, tmp, 1)) return;
	}
	break;

    default:
	break;
    }

    if (bufopen(&prn)) return;

    if (i == BINOMIAL_PVAL) {
	if (do_binomial_pdf(cmd, prn)) {
	    return;
	}
    } else {
	batch_pvalue(cmd, (const double **) Z, datainfo, prn);
    }

    view_buffer(prn, 78, 200, _("gretl: p-value"), PVALUE, NULL);
}

static void print_pv (PRN *prn, double p1, double p2)
{
    pprintf(prn, _("Two-tailed p-value = %.4g\n(one-tailed = %.4g)\n"),
	    p1, p2);
}

/* v. rough 99.9 percent critical value of chi-square */

static double chi_crit (const int df)
{
    double x = 10.0;
    
    while (chisq(x, df) > .001) x += .5;
    return x;
}

/* v. rough 99.9 percent critical value of F */

static double f_crit (const int df1, const int df2)
{
    double x = 2.0;

    while (fdist(x, df1, df2) > .001) x += .5;
    return x;
}

static void htest_graph (int dist, double x, int df1, int df2)
{
    double xx, prange, spike = 0.0;
    FILE *fp = NULL;

    if (gnuplot_init(PLOT_SAMPLING_DIST, &fp)) return;

    fprintf(fp, "set key right top\n");

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "C");
#endif

    if (dist == NORMAL_DIST || dist == T_DIST) {
	xx = fabs(x);
	prange = ((xx > 3.5)? xx + .5 : 3.5);
	spike = .25;
	fprintf(fp, "set xrange [%.3f:%.3f]\n", -prange, prange);
	fprintf(fp, "set yrange [0:.50]\n");
	fprintf(fp, "set xlabel '%s'\n", I_("Standard errors"));
    } else if (dist == CHISQ_DIST || dist == F_DIST) {
	prange = (dist == CHISQ_DIST)? chi_crit(df1) : f_crit(df1, df2);
	if (x > prange) 
	    prange = 1.1 * x;
	spike = 1.0 / prange;
	fprintf(fp, "set xrange [0:%.3f]\n", prange);
    } 

    /* required variables and formulae */
    if (dist == T_DIST) {
	fputs("# literal lines = 2\n", fp);
	fprintf(fp, "df1=%.1f\n", (double) df1);
	fprintf(fp, "Binv(p,q)=exp(lgamma(p+q)-lgamma(p)-lgamma(q))\n");
    } else if (dist == CHISQ_DIST) {
	fprintf(fp, "# literal lines = %d\n", (df1 > 69)? 3 : 2);
	fprintf(fp, "df1=%.1f\n", (double) df1);
	if (df1 > 69) {
	    fputs("log2=log(2.0)\n", fp);
	    fputs("chi(x)=exp((0.5*df1-1.0)*log(x)-0.5*x-"
		  "lgamma(0.5*df1)-df1*0.5*log2)\n", fp);
	} else {
	    fprintf(fp, "chi(x)=x**(0.5*df1-1.0)*exp(-0.5*x)/gamma(0.5*df1)"
		    "/2**(0.5*df1)\n");
	}
    } else if (dist == F_DIST) {
	fputs("# literal lines = 4\n", fp);
	fprintf(fp, "df1=%.1f\n", (double) df1);
	fprintf(fp, "df2=%.1f\n", (double) df2);
	fputs("Binv(p,q)=exp(lgamma(p+q)-lgamma(p)-lgamma(q))\n", fp);
	fputs("f(x)=Binv(0.5*df1,0.5*df2)*(df1/df2)**(0.5*df1)"
	      "*x**(0.5*df1-1.0)/(1.0+df1/df2*x)**(0.5*(df1+df2))\n", fp);
    }	

    fprintf(fp, "plot \\\n");

    if (dist == NORMAL_DIST) {
	fprintf(fp, "(1/(sqrt(2*pi))*exp(-(x)**2/2)) "
		"title '%s' w lines , \\\n",
		I_("Gaussian sampling distribution"));
    } else if (dist == T_DIST) {
	char tmp[64];

	sprintf(tmp, I_("t(%d) sampling distribution"), df1);
	fprintf(fp, "Binv(0.5*df1,0.5)/sqrt(df1)*(1.0+(x*x)/df1)"
		"**(-0.5*(df1+1.0)) "
		"title '%s' w lines , \\\n", tmp);
    } else if (dist == CHISQ_DIST) {
	char tmp[64];
	
	sprintf(tmp, I_("Chi-square(%d) sampling distribution"), df1);
	fprintf(fp, "chi(x) title '%s' w lines , \\\n", tmp);
    } else if (dist == F_DIST) {
	char tmp[64];

	sprintf(tmp, I_("F(%d, %d) sampling distribution"), df1, df2);
	fprintf(fp, "f(x) title '%s' w lines , \\\n", tmp);
    }

    fprintf(fp, "'-' using 1:2 title '%s' w impulses\n",
	    I_("test statistic"));
    fprintf(fp, "%g %g\n", x, spike);
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

static void h_test (GtkWidget *w, gpointer data)
{
    test_t *test = (test_t *) data;
    int j, n1, n2, grf;
    double x[5], sderr, ts, pv;
    const gchar *tmp;
    PRN *prn;

    if (bufopen(&prn)) return;

    grf = GTK_TOGGLE_BUTTON(test->graph)->active;

    x[4] = 0.0;

    switch (test->code) {

    case ONE_MEAN:
	for (j=0; j<2; j++) {
	    /* get sample mean and std dev */
	    tmp = gtk_entry_get_text(GTK_ENTRY(test->entry[j]));
	    x[j] = getval(tmp, prn, (j==1)? 1 : 0);
	    if (na(x[j])) return;
	}
	tmp = gtk_entry_get_text(GTK_ENTRY(test->entry[2]));
	if ((n1 = getint(tmp, prn)) == -1) return;
	tmp = gtk_entry_get_text(GTK_ENTRY(test->entry[3]));
	x[2] = getval(tmp, prn, 0);
	if (na(x[2])) return;
	sderr = x[1] / sqrt((double) n1);
	ts = (x[0] - x[2])/sderr;
	pprintf(prn, _("Null hypothesis: population mean = %g\n"), x[2]);
	pprintf(prn, _("Sample size: n = %d\n"), n1);
	pprintf(prn, _("Sample mean = %g, std. deviation = %g\n"), 
		x[0], x[1]);
	if (GTK_TOGGLE_BUTTON(test->check)->active) {
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

    case ONE_VARIANCE:
	tmp = gtk_entry_get_text(GTK_ENTRY(test->entry[0]));
	x[0] = getval(tmp, prn, 1);
	if (na(x[0])) return;
	tmp = gtk_entry_get_text(GTK_ENTRY(test->entry[1]));
	if ((n1 = getint(tmp, prn)) == -1) return;
	tmp = gtk_entry_get_text(GTK_ENTRY(test->entry[2]));
	x[1] = getval(tmp, prn, 1);
	if (na(x[1])) return;
	ts = (n1 - 1) * x[0] / x[1];
	pprintf(prn, _("Null hypothesis: population variance = %g\n"), x[1]);
	pprintf(prn, _("Sample size: n = %d\n"), n1);
	pprintf(prn, _("Sample variance = %g\n"), x[0]);
	pprintf(prn, _("Test statistic: chi-square(%d) = %d * %g/%g = %g\n"), 
		n1-1, n1-1, x[0], x[1], ts);
	if (x[0] > x[1]) {
	    pv = chisq(ts, n1 - 1);
	} else {
	    pv = 1.0 - chisq(ts, n1 -1);
	}
	print_pv(prn, 2.0 * pv, pv);
	if (grf) htest_graph(2, ts, n1 - 1, 0);
	break;

    case ONE_PROPN: /* proportion */
	tmp = gtk_entry_get_text(GTK_ENTRY(test->entry[0]));
	x[0] = getval(tmp, prn, 1);
	if (na(x[0])) return;
	tmp = gtk_entry_get_text(GTK_ENTRY(test->entry[1]));
	if ((n1 = getint(tmp, prn)) == -1) return;
	tmp = gtk_entry_get_text(GTK_ENTRY(test->entry[2]));
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

    case TWO_MEANS:
	for (j=0; j<2; j++) {
	    /* mean and std dev, sample 1 */
	    tmp = gtk_entry_get_text(GTK_ENTRY(test->entry[j]));
	    x[j] = getval(tmp, prn, (j==1)? 1 : 0);
	    if (na(x[j])) return;
	}
	tmp = gtk_entry_get_text(GTK_ENTRY(test->entry[2]));
	if ((n1 = getint(tmp, prn)) == -1) return;
	for (j=2; j<4; j++) {
	    /* mean and std dev, sample 2 */
	    tmp = gtk_entry_get_text(GTK_ENTRY(test->entry[j+1]));
	    x[j] = getval(tmp, prn, (j==3)? 1 : 0);
	    if (na(x[j])) return;
	}
	tmp = gtk_entry_get_text(GTK_ENTRY(test->entry[5]));
	if ((n2 = getint(tmp, prn)) == -1) return;
	tmp = gtk_entry_get_text(GTK_ENTRY(test->entry[6]));
	x[4] = getval(tmp, prn, 0);
	if (na(x[4])) return;
	pprintf(prn, _("Null hypothesis: Difference of means = %g\n"), x[4]);
	pprintf(prn, _("Sample 1:\n n = %d, mean = %g, s.d. = %g\n"),
		n1, x[0], x[1]);
	pprintf(prn, _("Sample 2:\n n = %d, mean = %g, s.d. = %g\n"),
		n2, x[2], x[3]);
	/* are we assuming a common variance? */
	j = 0;
	if (GTK_TOGGLE_BUTTON(test->check)->active) j = 1;
	if (n1 < 30 || n2 < 30) j = 2;
	if (j > 0) {
	    ts = ((n1-1)*x[1]*x[1] + (n2-1)*x[3]*x[3])/(n1 + n2 - 2);
	    sderr = sqrt(ts / n1 + ts /n2);
	} else {
	    sderr = sqrt(x[1] * x[1] / n1 + x[3] * x[3] / n2);
	}
	ts = (x[0] - x[2] - x[4]) / sderr;
	if (j) {
	    if (j == 2)
		pprintf(prn, _("Small samples: assuming normality and common "
			       "variance\n"));
	    pprintf(prn, _("Test statistic: t(%d) = (%g - %g)/%g = %g\n"),
		    n1 + n2 - 2, x[0], x[2], sderr, ts);
	    if (ts > 0) {
		pv = t_pvalue_2(ts, n1 + n2 - 2);
	    } else {
		pv = t_pvalue_2(-ts, n1 + n2 - 2);
	    }
	    print_pv(prn, pv, 0.5 * pv);
	    if (grf) htest_graph(1, ts, n1 + n2 - 2, 0);
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

    case TWO_VARIANCES:
	tmp = gtk_entry_get_text(GTK_ENTRY(test->entry[0]));
	x[0] = getval(tmp, prn, 1);
	if (na(x[0])) return;
	tmp = gtk_entry_get_text(GTK_ENTRY(test->entry[1]));
	if ((n1 = getint(tmp, prn)) == -1) return;
	tmp = gtk_entry_get_text(GTK_ENTRY(test->entry[2]));
	x[1] = getval(tmp, prn, 1);
	if (na(x[1])) return;
	tmp = gtk_entry_get_text(GTK_ENTRY(test->entry[3]));
	if ((n2 = getint(tmp, prn)) == -1) return;
	pprintf(prn, _("Null hypothesis: The population variances are "
		       "equal\n"));
	pprintf(prn, _("Sample 1:\n n = %d, variance = %g\n"), n1, x[0]);
	pprintf(prn, _("Sample 2:\n n = %d, variance = %g\n"), n2, x[1]);
	if (x[0] > x[1]) {
	    ts = x[0] / x[1];
	    pprintf(prn, _("Test statistic: F(%d, %d) = %g\n"), 
		    n1 - 1, n2 - 1, ts);
	    pv = fdist(ts, n1 - 1, n2 - 1);
	} else {
	    ts = x[1]/x[0];
	    pprintf(prn, _("Test statistic: F(%d, %d) = %g\n"), 
		    n2 - 1, n1 - 1, ts);
	    pv = fdist(ts, n2 - 1, n1 - 1);
	}
	print_pv(prn, 2.0 * pv, pv);
	if (grf) htest_graph(3, ts, n1 - 1, n2 - 1);
	break;

    case TWO_PROPNS: /* two proportions */
	tmp = gtk_entry_get_text(GTK_ENTRY(test->entry[0]));
	x[0] = getval(tmp, prn, 1);
	if (na(x[0])) return;

	tmp = gtk_entry_get_text(GTK_ENTRY(test->entry[1]));
	if ((n1 = getint(tmp, prn)) == -1) return;
	tmp = gtk_entry_get_text(GTK_ENTRY(test->entry[2]));
	x[1] = getval(tmp, prn, 1);
	if (na(x[1])) return;

	tmp = gtk_entry_get_text(GTK_ENTRY(test->entry[3]));
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

static void h_test_global (GtkWidget *w, gpointer data)
{
    test_t **test = (test_t **) data;
    int i;

    i = gtk_notebook_get_current_page(GTK_NOTEBOOK(test[0]->book));
    h_test(NULL, test[i]);
}

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

static void make_dist_tab (GtkWidget *notebook, int code, lookup_t **pval) 
{
    GtkWidget *tempwid, *box, *tbl;
    gint tbl_len;
    const gchar *titles[] = {
	N_("normal"), 
	N_(" t "), 
	N_("chi-square"), 
	N_(" F "), 
	N_("gamma"),
	N_("binomial")
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

    case BINOMIAL_PVAL:
	add_lookup_entry(tbl, &tbl_len, N_("Prob"), pval, code, 1);
	add_lookup_entry(tbl, &tbl_len, N_("trials"), pval, code, 1);
	add_lookup_entry(tbl, &tbl_len, N_("x-value"), pval, code, 1);
	break;

    default:
	break;
    } 
}

static void trash_pval (GtkWidget *w, gpointer data)
{
    lookup_t **look = (lookup_t **) data;
    int i;

    for (i=0; i<NPVAL; i++) {
	free(look[i]);
    }
    free(look);
}

static void trash_look (GtkWidget *w, gpointer data)
{
    lookup_t **look = (lookup_t **) data;
    int i;

    for (i=0; i<NLOOKUPS; i++) {
	free(look[i]);
    }
    free(look);
}

static void trash_test (GtkWidget *w, gpointer data)
{
    test_t **test = (test_t **) data;
    int i;

    for (i=0; i<NTESTS; i++) {
	free(test[i]);
    }
    free(test);
}

static int get_restriction_vxy (const char *s, int *vx, int *vy, 
				GretlOp *yop, double *yval)
{
    char test[16];
    char *p, *q = NULL;
    char *str = g_strdup(s);
    int err = 0;

    if (str == NULL) {
	return 1;
    }

    p = strchr(str, '(');
    *p = 0;
    p++;

    if (sscanf(str, "%8s", test) != 1) {
	err = 1;
    } else {
	*vx = varindex(datainfo, test);
	if (*vx >= datainfo->v) {
	    err = 1;
	}
    }

    if (!err) {
	int len = strcspn(p, "=<>!");

	q = p + len;
	if (*q == 0) {
	    err = 1;
	} else {
	    len = 1;
	    if (!strncmp(q, "!=", 2)) {
		*yop = OP_NEQ;
		len = 2;
	    } else if (!strncmp(q, ">=", 2)) {
		*yop = OP_GTE;
		len = 2;
	    } else if (!strncmp(q, "<=", 2)) {
		*yop = OP_LTE;
		len = 2;
	    } else if (*q == '=') {
		*yop = OP_EQ;
	    } else if (*q == '>') {
		*yop = OP_GT;
	    } else if (*q == '<') {
		*yop = OP_LT;
	    } else {
		err = 1;
	    }
	    if (!err) {
		*q = 0;
		q += len;
	    }
	}
    }

    if (!err) {
	if (sscanf(p, "%8s", test) != 1) {
	    err = 1;
	} else {
	    *vy = varindex(datainfo, test);
	    if (*vx >= datainfo->v) {
		err = 1;
	    }
	}
    }

    if (!err) {
	if (sscanf(q, "%lf", yval) != 1) {
	    err = 1;
	}
    }
    
    g_free(str);

    return err;
}

/* fill out the sample statistics boxes based on the user's
   choice or variable (or variable plus restriction) */

static void populate_stats (GtkWidget *w, gpointer p)
{
    test_t *test = g_object_get_data(G_OBJECT(p), "test");
    int pos = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(p), "pos"));
    gchar **pbuf = g_object_get_data(G_OBJECT(p), "pbuf");
    int t, n = datainfo->t2 - datainfo->t1 + 1;
    int vx, vy = -1;
    GretlOp yop;
    const gchar *buf;
    char numstr[16];
    double x1, x2, yval;

    g_return_if_fail(GTK_IS_COMBO(p));
    if (!GTK_WIDGET_SENSITIVE(p)) {
	return;
    }

    g_return_if_fail(GTK_IS_ENTRY(GTK_COMBO(p)->entry));
    buf = gtk_entry_get_text(GTK_ENTRY(GTK_COMBO(p)->entry));
    if (*buf == 0) {
	return;
    }

    if (pbuf != NULL) {
	if (*pbuf != NULL && !strcmp(buf, *pbuf)) {
	    /* no real change */
	    return;
	}
	if (*pbuf != NULL) {
	    free(*pbuf);
	    *pbuf = NULL;
	}
	*pbuf = g_strdup(buf);
    }

    if (strchr(buf, '(') != NULL) {
	/* e.g. "cholest (gender = 1)" */
	if (get_restriction_vxy(buf, &vx, &vy, &yop, &yval)) {
	    return;
	}
    } else {
	vx = varindex(datainfo, buf);
	if (vx >= datainfo->v) {
	    return;
	}
    }

    /* scalars are not valid input in this context */
    if (!datainfo->vector[vx] || (vy > 0 && !datainfo->vector[vy])) {
	errbox(_("Invalid entry"));
	return;
    }

    for (t=datainfo->t1; t<=datainfo->t2; t++) {
	if (na(Z[vx][t]) || (vy > 0 && eval_ytest(Z[vy][t], yop, yval))) {
	    n--;
	}
    }

    if (n == 0) {		
	sprintf(errtext, _("Data missing for variable '%s'"),
		datainfo->varname[vx]);
	errbox(errtext);
	return;
    }

    if (test->code == ONE_MEAN || test->code == TWO_MEANS) {
	if (vy < 0) {
	    x1 = gretl_mean(datainfo->t1, datainfo->t2, Z[vx]);
	    x2 = gretl_stddev(datainfo->t1, datainfo->t2, Z[vx]);
	} else {
	    x1 = gretl_restricted_mean(datainfo->t1, datainfo->t2, 
				      Z[vx], Z[vy], yop, yval);
	    x2 = gretl_restricted_stddev(datainfo->t1, datainfo->t2, 
					 Z[vx], Z[vy], yop, yval);
	}
	sprintf(numstr, "%.10g", x1);
	gtk_entry_set_text(GTK_ENTRY(test->entry[pos]), numstr);
	sprintf(numstr, "%.10g", x2);
	gtk_entry_set_text(GTK_ENTRY(test->entry[pos + 1]), numstr);
	sprintf(numstr, "%d", n);
	gtk_entry_set_text(GTK_ENTRY(test->entry[pos + 2]), numstr);
    } else if (test->code == ONE_VARIANCE || test->code == TWO_VARIANCES) {
	if (vy < 0) {
	    x1 = gretl_variance(datainfo->t1, datainfo->t2, Z[vx]);
	} else {
	    x1 = gretl_restricted_variance(datainfo->t1, datainfo->t2, 
					   Z[vx], Z[vy], yop, yval);
	}
	sprintf(numstr, "%.10g", x1);
	gtk_entry_set_text(GTK_ENTRY(test->entry[pos]), numstr);
	sprintf(numstr, "%d", n);
	gtk_entry_set_text(GTK_ENTRY(test->entry[pos + 1]), numstr);
    } 
}

static void add_vars_to_combo (GtkWidget *w, int pos)
{
    GList *vlist = NULL;
    int i, vmin = (pos > 0)? 2 : 1;

    for (i=vmin; i<datainfo->v; i++) {
	if (!is_hidden_variable(i, datainfo) && datainfo->vector[i]) {
	    vlist = g_list_append(vlist, datainfo->varname[i]);
	}
    }

    if (pos > 0) {
	/* add first variable at the end of the list */
	for (i=1; i<datainfo->v; i++) {
	    if (!is_hidden_variable(i, datainfo) && datainfo->vector[i]) {
		vlist = g_list_append(vlist, datainfo->varname[i]);
		break;
	    }
	}
    }	

    gtk_combo_set_popdown_strings(GTK_COMBO(w), vlist);
}

static void switch_combo_ok (GtkWidget *b, gpointer p)
{
    test_t *test = g_object_get_data(G_OBJECT(p), "test");
    int use_combo = GTK_TOGGLE_BUTTON(b)->active;
    int pos = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(p), "pos"));
    int maxent = 0;

    gtk_widget_set_sensitive(GTK_WIDGET(p), use_combo);

    if (test->code == ONE_MEAN || test->code == TWO_MEANS) {
	maxent = pos + 3;
    } else if (test->code == ONE_VARIANCE || test->code == TWO_VARIANCES) {
	maxent = pos + 2;
    }

    if (use_combo) {
	populate_stats(NULL, p);
    }
}

static gint catch_combo_key (GtkWidget *w, GdkEventKey *key, gpointer p)
{
    if (key->keyval == GDK_Return) { 
	populate_stats(NULL, p);
        return TRUE;
    } 

    return FALSE;
}

static void free_pbuf (GtkWidget *w, gpointer p)
{
    gchar **pbuf = g_object_get_data(G_OBJECT(w), "pbuf");

    if (pbuf != NULL) {
	if (*pbuf != NULL) {
	    free(*pbuf);
	}
	free(pbuf);
	g_object_set_data(G_OBJECT(w), "pbuf", NULL);
    }
}

static void select_child_callback (GtkList *l, GtkWidget *w, gpointer p)
{
    populate_stats(NULL, p);
}

static void add_test_combo (GtkWidget *tbl, gint *tbl_len, 
			    test_t *test, int pos)
{
    GtkWidget *button, *tmp;
    gchar **pbuf;

    *tbl_len += 1;
    gtk_table_resize(GTK_TABLE(tbl), *tbl_len, 2);
    button = gtk_check_button_new_with_label(_("Use variable from dataset"));
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      button, 0, 1, *tbl_len - 1, *tbl_len);
    gtk_widget_show(button);

    tmp = gtk_combo_new();
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      tmp, 1, 2, *tbl_len - 1, *tbl_len);
    gtk_widget_show(tmp);
    g_object_set_data(G_OBJECT(tmp), "test", test);
    g_object_set_data(G_OBJECT(tmp), "pos", GINT_TO_POINTER(pos));

    pbuf = malloc(sizeof *pbuf);
    if (pbuf != NULL) {
	*pbuf = NULL;
	g_object_set_data(G_OBJECT(tmp), "pbuf", pbuf);
	g_signal_connect(G_OBJECT(tmp), "destroy", G_CALLBACK(free_pbuf), NULL);
    }

    if (pos > 0) {
	test->combo[1] = tmp;
    } else {
	test->combo[0] = tmp;
    }

    add_vars_to_combo(tmp, pos);
    gtk_widget_set_sensitive(tmp, FALSE);
    gtk_combo_disable_activate(GTK_COMBO(tmp));

    g_signal_connect(G_OBJECT(GTK_ENTRY(GTK_COMBO(tmp)->entry)), "key_press_event",
		     G_CALLBACK(catch_combo_key), tmp);
    g_signal_connect(G_OBJECT(GTK_LIST(GTK_COMBO(tmp)->list)), "select-child",
		     G_CALLBACK(select_child_callback), tmp);

    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(switch_combo_ok), tmp);
}

static void add_test_entry (GtkWidget *tbl, gint *tbl_len, 
			    const gchar *label, test_t *test, int i)
{
    GtkWidget *tmp;

    *tbl_len += 1;
    gtk_table_resize(GTK_TABLE(tbl), *tbl_len, 2);
    tmp = gtk_label_new(label);
    gtk_misc_set_alignment(GTK_MISC(tmp), 1, 0.5);
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      tmp, 0, 1, *tbl_len - 1, *tbl_len);
    gtk_widget_show(tmp);
    tmp = gtk_entry_new();
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      tmp, 1, 2, *tbl_len - 1, *tbl_len);
    gtk_widget_show(tmp);
    test->entry[i] = tmp;

    g_signal_connect(G_OBJECT(tmp), "activate", 
		     G_CALLBACK(h_test), test);
}

static void add_test_label (GtkWidget *tbl, gint *tbl_len, 
			    const gchar *label)
{
    GtkWidget *tmp;

    *tbl_len += 1;
    gtk_table_resize(GTK_TABLE(tbl), *tbl_len, 2);
    tmp = gtk_label_new(label);
    gtk_misc_set_alignment(GTK_MISC(tmp), 0, 0.5);
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      tmp, 0, 2, *tbl_len - 1, *tbl_len);
    gtk_widget_show(tmp);
}

static void add_test_check (GtkWidget *tbl, gint *tbl_len, 
			    const gchar *label, test_t *test)
{
    GtkWidget *tempwid;

    *tbl_len += 1;
    gtk_table_resize(GTK_TABLE(tbl), *tbl_len, 2);
    tempwid = gtk_check_button_new_with_label(label);
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      tempwid, 0, 2, *tbl_len - 1, *tbl_len);
    gtk_widget_show(tempwid);
    test->check = tempwid;
}

static int n_non_hidden_vars (void)
{
    int i, nv = 0;

    if (datainfo != NULL) {
	for (i=1; i<datainfo->v; i++) {
	    if (!is_hidden_variable(i, datainfo)) {
		nv++;
	    }
	}
    }

    return nv;
}

static void make_test_tab (GtkWidget *notebook, int code, test_t *test) 
{
    GtkWidget *tempwid, *box, *tbl;
    int nv = n_non_hidden_vars();
    gint i, tbl_len;
    const gchar *titles[] = {
	N_("mean"), 
	N_("variance"), 
	N_("proportion"),
	N_("2 means"), 
	N_("2 variances"), 
	N_("2 proportions")
    };
   
    box = gtk_vbox_new(FALSE, 0);
    gtk_container_set_border_width(GTK_CONTAINER(box), 10);
    gtk_widget_show(box);

    tempwid = gtk_label_new(_(titles[code]));
    gtk_widget_show(tempwid);
    gtk_notebook_append_page(GTK_NOTEBOOK(notebook), box, tempwid);   

    tbl_len = 1;
    tbl = gtk_table_new(tbl_len, 2, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(tbl), 5);
    gtk_table_set_col_spacings(GTK_TABLE(tbl), 5);
    gtk_box_pack_start(GTK_BOX(box), tbl, FALSE, FALSE, 0);
    gtk_widget_show(tbl);

    test->combo[0] = test->combo[1] = NULL;

    for (i=0; i<NTESTENTRY; i++) {
	test->entry[i] = NULL;
    }

    if ((code == ONE_MEAN || code == ONE_VARIANCE) && nv > 0) {
	add_test_combo(tbl, &tbl_len, test, 0);
    } else if ((code == TWO_MEANS || code == TWO_VARIANCES) && nv > 1) {
	add_test_combo(tbl, &tbl_len, test, 0);
    }
   
    switch (code) {

    case ONE_MEAN: 
	add_test_entry(tbl, &tbl_len, _("sample mean"), test, 0);
	add_test_entry(tbl, &tbl_len, _("std. deviation"), test, 1);
	add_test_entry(tbl, &tbl_len, _("sample size"), test, 2);
	add_test_entry(tbl, &tbl_len, _("H0: mean ="), test, 3);
	add_test_check(tbl, &tbl_len, _("Assume standard deviation is "
		       "population value"), test);
	break;

    case ONE_VARIANCE: 
	add_test_entry(tbl, &tbl_len, _("sample variance"), test, 0);
	add_test_entry(tbl, &tbl_len, _("sample size"), test, 1);
	add_test_entry(tbl, &tbl_len, _("H0: variance ="), test, 2);
	break;

    case ONE_PROPN: /* proportion */
	add_test_entry(tbl, &tbl_len, _("sample proportion"), test, 0);
	add_test_entry(tbl, &tbl_len, _("sample size"), test, 1);
	add_test_entry(tbl, &tbl_len, _("H0: proportion ="), test, 2);
	break;

    case TWO_MEANS:
	add_test_entry(tbl, &tbl_len, _("mean of sample 1"), test, 0);
	add_test_entry(tbl, &tbl_len, _("std. deviation, sample 1"), test, 1);
	add_test_entry(tbl, &tbl_len, _("size of sample 1"), test, 2);
	if (nv > 1) {
	    add_test_combo(tbl, &tbl_len, test, 3);
	}
	add_test_entry(tbl, &tbl_len, _("mean of sample 2"), test, 3);
	add_test_entry(tbl, &tbl_len, _("std. deviation, sample 2"), test, 4);
	add_test_entry(tbl, &tbl_len, _("size of sample 2"), test, 5);
	add_test_entry(tbl, &tbl_len, _("H0: Difference of means ="), test, 6);
	gtk_entry_set_text(GTK_ENTRY(test->entry[6]), "0");
	add_test_check(tbl, &tbl_len, _("Assume common population standard "
		       "deviation"), test);
	break;

    case TWO_VARIANCES:
	add_test_entry(tbl, &tbl_len, _("variance of sample 1"), test, 0);
	add_test_entry(tbl, &tbl_len, _("size of sample 1"), test, 1);
	if (nv > 1) {
	    add_test_combo(tbl, &tbl_len, test, 2);
	}	
	add_test_entry(tbl, &tbl_len, _("variance of sample 2"), test, 2);
	add_test_entry(tbl, &tbl_len, _("size of sample 2"), test, 3);
	add_test_label(tbl, &tbl_len, _("H0: Ratio of variances = 1"));
	break;

    case TWO_PROPNS:
	add_test_entry(tbl, &tbl_len, _("proportion, sample 1"), test, 0);
	add_test_entry(tbl, &tbl_len, _("size of sample 1"), test, 1);
	add_test_entry(tbl, &tbl_len, _("proportion, sample 2"), test, 2);
	add_test_entry(tbl, &tbl_len, _("size of sample 2"), test, 3);
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
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tempwid), TRUE);
    gtk_widget_show(tempwid);

    test->graph = tempwid; 
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
    } else if (code == CALC_PVAL) {
	look = mymalloc(NPVAL * sizeof *look);
	if (look == NULL) return;
	statp = look;
    } else {
	look = mymalloc(NLOOKUPS * sizeof *look);
	if (look == NULL) return;
	statp = look;
    }	

    dialog = gretl_child_new(NULL);
    winptr[code] = dialog->win;
    dialog->data = &(winptr[code]);

    gtk_window_set_title(GTK_WINDOW(dialog->win), _(window_titles[code]));

    notebook = gtk_notebook_new();
    gtk_box_pack_start(GTK_BOX(dialog->vbox), notebook, TRUE, TRUE, 0);
    gtk_widget_show(notebook);

    if (code == CALC_TEST) {
	for (i=0; i<NTESTS; i++) {
	    test[i] = mymalloc(sizeof **test);
	    if (test[i] == NULL) return;
	    test[i]->book = notebook;
	    test[i]->code = i;
	    make_test_tab(notebook, i, test[i]);
	}
    } else if (code == CALC_PVAL) {
	for (i=0; i<NPVAL; i++) {
	    look[i] = mymalloc(sizeof **look);
	    if (look[i] == NULL) return;
	    look[i]->book = notebook;
	    make_dist_tab(notebook, i, look);
	}	
    } else {
	for (i=0; i<NLOOKUPS; i++) {
	    look[i] = mymalloc(sizeof **look);
	    if (look[i] == NULL) return;
	    look[i]->book = notebook;
	    make_lookup_tab(notebook, i, look);
	}
    }	

    /* OK button */
    tempwid = standard_button(GTK_STOCK_OK);
    GTK_WIDGET_SET_FLAGS(tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(dialog->action_area), 
		       tempwid, TRUE, TRUE, 0);

    g_signal_connect(G_OBJECT (tempwid), "clicked", 
		     (code == CALC_PVAL)? G_CALLBACK(get_pvalue) :
		     (code == CALC_DIST)? G_CALLBACK(get_critical) :
		     G_CALLBACK(h_test_global),
		     statp);
    gtk_widget_show(tempwid);

    /* Close button */
    tempwid = standard_button(GTK_STOCK_CLOSE);
    GTK_WIDGET_SET_FLAGS(tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(dialog->action_area), 
		       tempwid, TRUE, TRUE, 0);

    g_signal_connect(G_OBJECT(tempwid), "clicked", 
		     (code == CALC_PVAL)? G_CALLBACK(trash_pval) :
		     (code == CALC_TEST)? G_CALLBACK(trash_test) :
		     G_CALLBACK(trash_look),
		     statp);
    g_signal_connect(G_OBJECT(tempwid), "clicked", 
		     G_CALLBACK(delete_widget), 
		     dialog->win);
    gtk_widget_show(tempwid);

    /* Help button? */
    if (code == CALC_TEST) {
	tempwid = standard_button(GTK_STOCK_HELP);
	GTK_WIDGET_SET_FLAGS(tempwid, GTK_CAN_DEFAULT);
	gtk_box_pack_start(GTK_BOX(dialog->action_area), 
			   tempwid, TRUE, TRUE, 0);
	g_signal_connect(G_OBJECT(tempwid), "clicked", 
			 G_CALLBACK(context_help), 
			 GINT_TO_POINTER(HTEST));
	gtk_widget_show(tempwid);
    }	

    gtk_widget_show(dialog->win);
}





