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
#define NPTESTS 2
#define NPVAL 7
#define NLOOKUPS 5
#define NGRAPHS 4
#define NTESTENTRY 7
#define NLOOKUPENTRY 4

#include "gretl.h"
#include "calculator.h"
#include "dlgutils.h"
#include "gpt_control.h"

typedef struct CalcChild_ CalcChild;
typedef struct test_t_ test_t;
typedef struct lookup_t_ lookup_t;

struct CalcChild_ {
    int code;
    GtkWidget *dlg;
    GtkWidget *vbox;
    GtkWidget *bbox;
    GtkWidget *book;
    gpointer calcp;
    gpointer winp;
    png_plot *plot;
};

struct test_t_ {
    int code;
    GtkWidget *entry[NTESTENTRY];
    GtkWidget *combo[2];
    GtkWidget *check;
    GtkWidget *radio[3];
    GtkWidget *extra;
};

struct lookup_t_ {
    GtkWidget *entry[NLOOKUPENTRY];
};

enum {
    NORMAL_PVAL,
    T_PVAL,
    CHISQ_PVAL,
    F_PVAL,
    GAMMA_PVAL,
    BINOMIAL_PVAL,
    POISSON_PVAL
};

enum {
    ONE_MEAN,
    ONE_VARIANCE,
    ONE_PROPN,
    TWO_MEANS,
    TWO_VARIANCES,
    TWO_PROPNS
};

enum {
    NP_DIFF,
    NP_RUNS
};

/* if pos != 0, value must be positive or it is invalid */ 

static double getval (const char *s, PRN *prn, int pos)
{
    double x = NADBL;

    if (s == NULL || *s == '\0') {
	errbox(_("Incomplete entry"));
    } else {
	if (check_atof(s)) {
	    errbox(gretl_errmsg_get());
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

static void get_critical (GtkWidget *w, CalcChild *child)
{
    lookup_t **look = child->calcp;
    void *handle = NULL;
    void (*dw)(int, PRN *) = NULL;
    int d, n = -1, df = -1, err = 0;
    int winwidth = 60;
    int winheight = 200;
    double x = NADBL, a = 0.0;
    PRN *prn;

    d = gtk_notebook_get_current_page(GTK_NOTEBOOK(child->book));

    if (bufopen(&prn)) {
	return;
    }	

    if (d != DW_DIST) {
	a = atof(gtk_entry_get_text(GTK_ENTRY(look[d]->entry[0])));
    }

    switch (d) {
    case T_DIST:
    case CHISQ_DIST:
	df = atoi(gtk_entry_get_text(GTK_ENTRY(look[d]->entry[1])));
	break;
    case F_DIST:
	df = atoi(gtk_entry_get_text(GTK_ENTRY(look[d]->entry[1])));
	n = atoi(gtk_entry_get_text(GTK_ENTRY(look[d]->entry[2])));
	break;
    case DW_DIST:
	n = atoi(gtk_entry_get_text(GTK_ENTRY(look[d]->entry[0])));
	dw = gui_get_plugin_function("dw_lookup", &handle);
	winwidth = 77;
	winheight = 300;
	break;
    default:
	break;
    }

    /* sanity checks */
    if ((0 < d && d < 4 && df <= 0) || (d == F_DIST && n <= 0)) {
	errbox(_("Invalid degrees of freedom"));
	err = 1;
    } else if (d != DW_DIST && (a <= 0.0 || a > 1.0)) {
	errbox(_("Invalid right-tail probability"));
	err = 1;
    } else if (d == DW_DIST && n <= 0) {
	errbox(_("Invalid sample size"));
	err = 1;
    } else if (d == DW_DIST && dw == NULL)  {
	err = 1;
    }

    /* get the values */
    if (!err) {
	switch (d) {
	case NORMAL_DIST:
	    x = normal_critval(a);
	    break;
	case T_DIST:
	    x = t_critval(a, df);
	    break;
	case CHISQ_DIST:
	    x = chisq_critval(a, df);
	    break;
	case F_DIST:
	    x = f_critval(a, df, n);
	    break;
	case DW_DIST:
	    (*dw)(n, prn);
	    break;
	default:
	    break;
	}
    }

    if (d != DW_DIST) {
	if (na(x)) {
	    errbox(_("Failed to compute critical value"));
	    err = 1;
	} else {
	    switch (d) {
	    case NORMAL_DIST:
		pprintf(prn, "%s", _("Standard normal distribution"));
		break;
	    case T_DIST:
		pprintf(prn, "t(%d)", df);
		break;
	    case CHISQ_DIST:
		pprintf(prn, _("Chi-square(%d)"), df);
		break;
	    case F_DIST:
		pprintf(prn, "F(%d, %d)", df, n);
		break;
	    }

	    pputs(prn, "\n ");
	    pprintf(prn, _("right-tail probability = %g"), a);
	    pputs(prn, "\n ");
	    pprintf(prn, _("complementary probability = %g"), 1.0 - a);
	    if (a < 0.5 && (d == NORMAL_DIST || d == T_DIST)) {
		pputs(prn, "\n ");
		pprintf(prn, _("two-tailed probability = %g"), 2.0 * a);
	    }
	    pputs(prn, "\n\n ");
	    pprintf(prn, _("Critical value = %g"), x);
	    pputc(prn, '\n');
	}
    }	
	
    if (handle != NULL) {
	close_plugin(handle);
    }

    if (err) {
	gretl_print_destroy(prn);
    } else {
	view_buffer(prn, winwidth, winheight, _("gretl: critical values"), 
		    STAT_TABLE, NULL);
    }
}

static void get_pvalue (GtkWidget *w, CalcChild *child)
{
    lookup_t **pval = child->calcp;
    double pv, parm[3];
    char st = 0;
    int i, j, df;
    double val, xx, zz;
    const gchar *tmp;
    PRN *prn;

    i = gtk_notebook_get_current_page(GTK_NOTEBOOK(child->book));

    switch (i) {

    case NORMAL_PVAL:
	st = 'z';
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
	parm[0] = xx/val;
	break;

    case T_PVAL: 
	st = 't';
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
	parm[0] = df;
	parm[1] = xx / val;
	break;

    case CHISQ_PVAL:
	st = 'X';
	tmp = gtk_entry_get_text(GTK_ENTRY(pval[i]->entry[0]));
	df = atoi(tmp);   /* df */
	if (df <= 0) {
	    errbox(_("Invalid degrees of freedom"));
	    return;
	}	
	tmp = gtk_entry_get_text(GTK_ENTRY(pval[i]->entry[1]));
	xx = getval(tmp, NULL, 0); /* value */
	if (na(xx)) return;
	parm[0] = df;
	parm[1] = xx;
	break;

    case POISSON_PVAL: 
	st = 'P';
	tmp = gtk_entry_get_text(GTK_ENTRY(pval[i]->entry[0]));
	zz = getval(tmp, NULL, 0); /* mean */
	if (na(zz)) return;
	if (zz <= 0.0) {
	    errbox(_("Invalid mean"));
	    return;
	}	
	tmp = gtk_entry_get_text(GTK_ENTRY(pval[i]->entry[1]));
	df = atoi(tmp); /* value, actually */
	parm[0] = zz;
	parm[1] = df;
	break;

    case F_PVAL:
	st = 'F';
	for (j=0; j<2; j++) {
	    tmp = gtk_entry_get_text(GTK_ENTRY(pval[i]->entry[j]));
	    df = atoi(tmp);
	    if (df <= 0) {
		errbox(_("Invalid degrees of freedom"));
		return;
	    }
	    parm[j] = df;
	}
	tmp = gtk_entry_get_text(GTK_ENTRY(pval[i]->entry[2]));
	xx = getval(tmp, NULL, 0); /* value */
	if (na(xx)) return;
	parm[2] = xx;
	break;

    case GAMMA_PVAL: 
	st = 'G';
	for (j=0; j<3; j++) {
	    tmp = gtk_entry_get_text(GTK_ENTRY(pval[i]->entry[j]));
	    xx = getval(tmp, NULL, 0);
	    if (na(xx)) return;
	    parm[j] = xx;
	}
	break;

    case BINOMIAL_PVAL: 
	st = 'B';
	for (j=0; j<3; j++) {
	    tmp = gtk_entry_get_text(GTK_ENTRY(pval[i]->entry[j]));
	    xx = getval(tmp, NULL, 0);
	    if (na(xx)) return;
	    parm[j] = xx;
	}
	break;

    default:
	errbox(_("Failed to compute p-value"));
	return;
    }

    if (bufopen(&prn)) return;

    pv = gretl_get_pvalue(st, parm);

    if (na(pv)) {
	errbox(_("Failed to compute p-value"));
    } else {
	print_pvalue(st, parm, pv, prn);
	view_buffer(prn, 78, 200, _("gretl: p-value"), PVALUE, NULL);
    }
}

static void print_pv (PRN *prn, double p1, double p2)
{
    pprintf(prn, _("Two-tailed p-value = %.4g\n(one-tailed = %.4g)\n"),
	    p1, p2);
}

gchar *dist_graph_title (int dist, double x, int df1, int df2)
{
    gchar *s = NULL;

    if (na(x)) {
	if (dist == NORMAL_DIST) {
	    s = g_strdup(I_("Standard normal distribution"));
	} else if (dist == T_DIST) {
	    s = g_strdup_printf(I_("t(%d)"), df1);
	} else if (dist == CHISQ_DIST) {
	    s = g_strdup_printf(I_("Chi-square(%d)"), df1);
	} else if (dist == F_DIST) {
	    s = g_strdup_printf(I_("F(%d, %d)"), df1, df2);
	}	
    } else {
	if (dist == NORMAL_DIST) {
	    s = g_strdup(I_("Gaussian sampling distribution"));
	} else if (dist == T_DIST) {
	    s = g_strdup_printf(I_("t(%d) sampling distribution"), df1);
	} else if (dist == CHISQ_DIST) {
	    s = g_strdup_printf(I_("Chi-square(%d) sampling distribution"), df1);
	} else if (dist == F_DIST) {
	    s = g_strdup_printf(I_("F(%d, %d) sampling distribution"), df1, df2);
	}
    }

    return s;
}

gchar *dist_marker_line (int dist, int df1, int df2)
{
    gchar *s = NULL;

    if (dist == NORMAL_DIST) {
	s = g_strdup("# standard normal");
    } else if (dist == T_DIST) {
	s = g_strdup_printf("# t(%d)", df1);
    } else if (dist == CHISQ_DIST) {
	s = g_strdup_printf("# chi-square(%d)", df1);
    } else if (dist == F_DIST) {
	s = g_strdup_printf("# F(%d,%d)", df1, df2);
    }	

    return s;
}

static const char *formulae[] = {
    "Binv(p,q)=exp(lgamma(p+q)-lgamma(p)-lgamma(q))",
    "chi(x,m)=x**(0.5*m-1.0)*exp(-0.5*x)/gamma(0.5*m)/2**(0.5*m)",
    "log2=log(2.0)",
    "bigchi(x,m)=exp((0.5*m-1.0)*log(x)-0.5*x-lgamma(0.5*m)-df1*0.5*log2)",
    "f(x,m,n)=Binv(0.5*m,0.5*n)*(m/n)**(0.5*m)*"
    "x**(0.5*m-1.0)/(1.0+m/n*x)**(0.5*(m+n))"
};

const char *dist_formula (FormulaCode c)
{
    if (c <= F_F) {
	return formulae[c];
    }

    return "";
}

double dist_xmax (int d, int df1, int df2)
{
    double a = 0.005;

    if (d == F_DIST && df1 + df2 < 16) {
	a = 0.009;
    }

    return (d == CHISQ_DIST)? chisq_critval(a, df1) : 
		f_critval(a, df1, df2);
}

static void htest_graph (int d, double x, int df1, int df2)
{
    PlotType pt = (na(x))? PLOT_PROB_DIST : PLOT_H_TEST;
    double xx, prange, spike = 0.0;
    int bigchi = 0, nlit = 0;
    gchar *title = NULL;
    FILE *fp = NULL;

    if (gnuplot_init(pt, &fp)) {
	return;
    }

    if (d == CHISQ_DIST && df1 > 69) {
	bigchi = 1;
    }

    fputs("set key right top\n", fp);

    gretl_push_c_numeric_locale();

    if (na(x)) {
	/* no test statistic to be shown */
	if (d == NORMAL_DIST || d == T_DIST) {
	    prange = 5.0;
	    fprintf(fp, "set xrange [%.3f:%.3f]\n", -prange, prange);
	    fprintf(fp, "set yrange [0:.50]\n");
	} else if (d == CHISQ_DIST || d == F_DIST) {
	    prange = dist_xmax(d, df1, df2);
	    fprintf(fp, "set xrange [0:%.3f]\n", prange);
	}	    
    } else {
	/* set range based on test stat */
	if (d == NORMAL_DIST || d == T_DIST) {
	    xx = fabs(x);
	    prange = ((xx > 3.5)? xx + .5 : 3.5);
	    spike = .25;
	    fprintf(fp, "set xrange [%.3f:%.3f]\n", -prange, prange);
	    fprintf(fp, "set yrange [0:.50]\n");
	    fprintf(fp, "set xlabel '%s'\n", I_("Standard errors"));
	} else if (d == CHISQ_DIST || d == F_DIST) {
	    prange = dist_xmax(d, df1, df2);
	    if (x > prange) {
		prange = 1.1 * x;
	    }
	    spike = 1.0 / prange;
	    fprintf(fp, "set xrange [0:%.3f]\n", prange);
	} 
    }

    /* header */
    nlit = (d == NORMAL_DIST)? 1 :
	(d == T_DIST)? 3 : 
	(bigchi)? 4 :
	(d == CHISQ_DIST)? 3 : 5;
    fprintf(fp, "# literal lines = %d\n", nlit);
    title = dist_marker_line(d, df1, df2);
    fprintf(fp, "%s\n", title);

    /* required variables and formulae */
    if (d == T_DIST) {
	fprintf(fp, "df1=%.1f\n", (double) df1);
	fprintf(fp, "%s\n", formulae[F_BINV]);
    } else if (d == CHISQ_DIST) {
	fprintf(fp, "df1=%.1f\n", (double) df1);
	if (bigchi) {
	    fprintf(fp, "%s\n", formulae[F_LOG2]);
	    fprintf(fp, "%s\n", formulae[F_BIGCHI]);
	} else {
	    fprintf(fp, "%s\n", formulae[F_CHI]);
	}
    } else if (d == F_DIST) {
	fprintf(fp, "df1=%.1f\n", (double) df1);
	fprintf(fp, "df2=%.1f\n", (double) df2);
	fprintf(fp, "%s\n", formulae[F_BINV]);
	fprintf(fp, "%s\n", formulae[F_F]);
    }

    g_free(title);

    fprintf(fp, "plot \\\n");

    title = dist_graph_title(d, x, df1, df2);

    if (d == NORMAL_DIST) {
	fprintf(fp, "(1/(sqrt(2*pi))*exp(-(x)**2/2)) "
		"title '%s' w lines", title);
    } else if (d == T_DIST) {
	fprintf(fp, "Binv(0.5*df1,0.5)/sqrt(df1)*(1.0+(x*x)/df1)"
		"**(-0.5*(df1+1.0)) "
		"title '%s' w lines", title);
    } else if (d == CHISQ_DIST) {
	fprintf(fp, "%s(x,df1) title '%s' w lines", (bigchi)? "bigchi" : "chi",
		title);
    } else if (d == F_DIST) {
	fprintf(fp, "f(x,df1,df2) title '%s' w lines", title);
    }

    if (!na(x)) {
	fputs(" , \\\n", fp);
	fprintf(fp, "'-' using 1:2 title '%s' w impulses\n",
		I_("test statistic"));
	fprintf(fp, "%g %g\n", x, spike);
	fputs("e\n", fp);
    } else {
	fputc('\n', fp);
    }

    gretl_pop_c_numeric_locale();

    fclose(fp);

    if (gnuplot_make_graph()) {
	errbox(_("gnuplot command failed"));
    } else {
	register_graph();
    }
}

static void np_test (GtkWidget *w, test_t *test)
{
    gretlopt opt = OPT_NONE;
    const char *var1, *var2;
    int v1, v2 = 0;
    PRN *prn = NULL;
    int err = 0;

    var1 = gtk_entry_get_text(GTK_ENTRY(test->entry[0]));
    v1 = varindex(datainfo, var1);

    if (v1 == datainfo->v) {
	gui_errmsg(E_UNKVAR);
	return;
    }

    if (test->code == NP_DIFF) {
	var2 = gtk_entry_get_text(GTK_ENTRY(test->entry[1]));
	v2 = varindex(datainfo, var2);
	if (v2 == datainfo->v) {
	    gui_errmsg(E_UNKVAR);
	    return;
	}
    }

    if (bufopen(&prn)) {
	return;
    }

    if (test->code == NP_DIFF) {
	int list[3] = { 2, v1, v2 };

	if (test->extra != NULL &&
	    GTK_TOGGLE_BUTTON(test->extra)->active) {
	    opt |= OPT_V;
	}

	if (GTK_TOGGLE_BUTTON(test->radio[0])->active) {
	    opt |= OPT_G;
	} else if (GTK_TOGGLE_BUTTON(test->radio[1])->active) {
	    opt |= OPT_R;
	} else if (GTK_TOGGLE_BUTTON(test->radio[2])->active) {
	    opt |= OPT_I;
	}

	err = diff_test(list, (const double **) Z, datainfo, 
			opt, prn);
    } else if (test->code == NP_RUNS) {
	if (test->extra != NULL &&
	    GTK_TOGGLE_BUTTON(test->extra)->active) {
	    opt |= OPT_D;
	}
	err = runs_test(v1, (const double **) Z, datainfo, 
			opt, prn);
    }	

    if (err) {
	gui_errmsg(err);
    } else {
	view_buffer(prn, 78, 380, "gretl: nonparametric test", 
		    PRINT, NULL);
    }
}

/* FIXME : should we record a relevant command when a test is done
   using dataset variables? (Two means, or two variances)
*/

static void h_test (GtkWidget *w, test_t *test)
{
    int j, n1, n2, grf;
    double x[5], sderr, ts, pv, z;
    const gchar *tmp;
    PRN *prn;

    if (bufopen(&prn)) {
	return;
    }

    grf = GTK_TOGGLE_BUTTON(test->extra)->active;

    x[4] = 0.0;

    switch (test->code) {

    case ONE_MEAN:
	for (j=0; j<2; j++) {
	    /* get sample mean and std dev */
	    tmp = gtk_entry_get_text(GTK_ENTRY(test->entry[j]));
	    x[j] = getval(tmp, prn, (j==1)? 1 : 0);
	    if (na(x[j])) {
		return;
	    }
	}

	tmp = gtk_entry_get_text(GTK_ENTRY(test->entry[2]));
	if ((n1 = getint(tmp, prn)) == -1) {
	    return;
	}

	tmp = gtk_entry_get_text(GTK_ENTRY(test->entry[3]));
	x[2] = getval(tmp, prn, 0);
	if (na(x[2])) {
	    return;
	}

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
	    if (grf) {
		htest_graph(0, ts, 0, 0);
	    }
	} else {
	    pprintf(prn, _("Test statistic: t(%d) = (%g - %g)/%g = %g\n"), n1-1,
		    x[0], x[2], sderr, ts);
	    pv = t_pvalue_2(ts, n1 - 1);
	    print_pv(prn, pv, 0.5 * pv);
	    if (grf) {
		htest_graph(1, ts, n1-1, 0);
	    }
	}
	break;

    case ONE_VARIANCE:
	tmp = gtk_entry_get_text(GTK_ENTRY(test->entry[0]));
	x[0] = getval(tmp, prn, 1);
	if (na(x[0])) {
	    return;
	}

	tmp = gtk_entry_get_text(GTK_ENTRY(test->entry[1]));
	if ((n1 = getint(tmp, prn)) == -1) {
	    return;
	}

	tmp = gtk_entry_get_text(GTK_ENTRY(test->entry[2]));
	x[1] = getval(tmp, prn, 1);
	if (na(x[1])) {
	    return;
	}

	ts = (n1 - 1) * x[0] / x[1];

	pprintf(prn, _("Null hypothesis: population variance = %g\n"), x[1]);
	pprintf(prn, _("Sample size: n = %d\n"), n1);
	pprintf(prn, _("Sample variance = %g\n"), x[0]);
	pprintf(prn, _("Test statistic: chi-square(%d) = %d * %g/%g = %g\n"), 
		n1-1, n1-1, x[0], x[1], ts);

	if (x[0] > x[1]) {
	    pv = chisq_cdf_comp(ts, n1 - 1);
	} else {
	    pv = chisq_cdf(ts, n1 - 1);
	}
	print_pv(prn, 2.0 * pv, pv);
	if (grf) {
	    htest_graph(2, ts, n1 - 1, 0);
	}
	break;

    case ONE_PROPN: /* proportion */
	tmp = gtk_entry_get_text(GTK_ENTRY(test->entry[0]));
	x[0] = getval(tmp, prn, 1);
	if (na(x[0])) {
	    return;
	}

	tmp = gtk_entry_get_text(GTK_ENTRY(test->entry[1]));
	if ((n1 = getint(tmp, prn)) == -1) {
	    return;
	}

	tmp = gtk_entry_get_text(GTK_ENTRY(test->entry[2]));
	x[1] = getval(tmp, prn, 1);
	if (na(x[1])) {
	    return;
	}

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
	if (grf) {
	    htest_graph(0, ts, 0, 0);
	}
	break;

    case TWO_MEANS:
	for (j=0; j<2; j++) {
	    /* mean and std dev, sample 1 */
	    tmp = gtk_entry_get_text(GTK_ENTRY(test->entry[j]));
	    x[j] = getval(tmp, prn, (j==1)? 1 : 0);
	    if (na(x[j])) {
		return;
	    }
	}

	tmp = gtk_entry_get_text(GTK_ENTRY(test->entry[2]));
	if ((n1 = getint(tmp, prn)) == -1) {
	    return;
	}

	for (j=2; j<4; j++) {
	    /* mean and std dev, sample 2 */
	    tmp = gtk_entry_get_text(GTK_ENTRY(test->entry[j+1]));
	    x[j] = getval(tmp, prn, (j==3)? 1 : 0);
	    if (na(x[j])) {
		return;
	    }
	}

	tmp = gtk_entry_get_text(GTK_ENTRY(test->entry[5]));
	if ((n2 = getint(tmp, prn)) == -1) {
	    return;
	}

	tmp = gtk_entry_get_text(GTK_ENTRY(test->entry[6]));
	x[4] = getval(tmp, prn, 0);
	if (na(x[4])) {
	    return;
	}

	pprintf(prn, _("Null hypothesis: Difference of means = %g\n"), x[4]);
	pputc(prn, '\n');

	pprintf(prn, _("Sample 1:\n n = %d, mean = %g, s.d. = %g\n"),
		n1, x[0], x[1]);

	z = x[1] / sqrt((double) n1);

	pprintf(prn, _(" standard error of mean = %g\n"), z);
	z *= tcrit95(n1 - 1);
	pprintf(prn, _(" 95%% confidence interval for mean: %g to %g\n"), 
		x[0] - z, x[0] + z);
	pputc(prn, '\n');

	pprintf(prn, _("Sample 2:\n n = %d, mean = %g, s.d. = %g\n"),
		n2, x[2], x[3]);

	z = x[3] / sqrt((double) n2);

	pprintf(prn, _(" standard error of mean = %g\n"), z);
	z *= tcrit95(n2 - 1);
	pprintf(prn, _(" 95%% confidence interval for mean: %g to %g\n"), 
		x[2] - z, x[2] + z);
	pputc(prn, '\n');

	if (GTK_TOGGLE_BUTTON(test->check)->active) {
	    /* the user specified a common variance */
	    j = 1;
	} else if (n1 < 30 || n2 < 30) {
	    /* flag for warning: unequal variances and small samples */
	    j = 2;
	} else {
	    j = 0;
	}

	if (j == 1) {
	    ts = ((n1-1) * x[1] * x[1] + (n2-1) * x[3] * x[3]) / (n1 + n2 - 2);
	    sderr = sqrt(ts / n1 + ts / n2);
	} else {
	    double v1 = x[1] * x[1] / n1;
	    double v2 = x[3] * x[3] / n2;

	    sderr = sqrt(v1 + v2);
	}

	ts = (x[0] - x[2] - x[4]) / sderr;

	if (j == 1) {
	    pprintf(prn, _("Test statistic: t(%d) = (%g - %g)/%g = %g\n"),
		    n1 + n2 - 2, x[0], x[2], sderr, ts);
	    if (ts > 0) {
		pv = t_pvalue_2(ts, n1 + n2 - 2);
	    } else {
		pv = t_pvalue_2(-ts, n1 + n2 - 2);
	    }
	    print_pv(prn, pv, 0.5 * pv);
	    if (grf) {
		htest_graph(1, ts, n1 + n2 - 2, 0);
	    }
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
	    if (grf) {
		htest_graph(0, ts, 0, 0);
	    }
	}
	
	if (j == 2) {
	    pputc(prn, '\n');
	    pprintf(prn, _("Warning: with small samples, asymptotic "
			   "approximation may be poor\n"));
	}

	break;

    case TWO_VARIANCES:
	tmp = gtk_entry_get_text(GTK_ENTRY(test->entry[0]));
	x[0] = getval(tmp, prn, 1);
	if (na(x[0])) {
	    return;
	}

	tmp = gtk_entry_get_text(GTK_ENTRY(test->entry[1]));
	if ((n1 = getint(tmp, prn)) == -1) {
	    return;
	}

	tmp = gtk_entry_get_text(GTK_ENTRY(test->entry[2]));
	x[1] = getval(tmp, prn, 1);
	if (na(x[1])) {
	    return;
	}

	tmp = gtk_entry_get_text(GTK_ENTRY(test->entry[3]));
	if ((n2 = getint(tmp, prn)) == -1) {
	    return;
	}

	pprintf(prn, _("Null hypothesis: The population variances are "
		       "equal\n"));
	pprintf(prn, _("Sample 1:\n n = %d, variance = %g\n"), n1, x[0]);
	pprintf(prn, _("Sample 2:\n n = %d, variance = %g\n"), n2, x[1]);

	if (x[0] > x[1]) {
	    ts = x[0] / x[1];
	    pprintf(prn, _("Test statistic: F(%d, %d) = %g\n"), 
		    n1 - 1, n2 - 1, ts);
	    pv = f_cdf_comp(ts, n1 - 1, n2 - 1);
	} else {
	    ts = x[1] / x[0];
	    pprintf(prn, _("Test statistic: F(%d, %d) = %g\n"), 
		    n2 - 1, n1 - 1, ts);
	    pv = f_cdf_comp(ts, n2 - 1, n1 - 1);
	}

	print_pv(prn, 2.0 * pv, pv);
	if (grf) {
	    htest_graph(3, ts, n1 - 1, n2 - 1);
	}
	break;

    case TWO_PROPNS: /* two proportions */
	tmp = gtk_entry_get_text(GTK_ENTRY(test->entry[0]));
	x[0] = getval(tmp, prn, 1);
	if (na(x[0])) {
	    return;
	}

	tmp = gtk_entry_get_text(GTK_ENTRY(test->entry[1]));
	if ((n1 = getint(tmp, prn)) == -1) {
	    return;
	}

	tmp = gtk_entry_get_text(GTK_ENTRY(test->entry[2]));
	x[1] = getval(tmp, prn, 1);
	if (na(x[1])) {
	    return;
	}

	tmp = gtk_entry_get_text(GTK_ENTRY(test->entry[3]));
	if ((n2 = getint(tmp, prn)) == -1) {
	    return;
	}

	pprintf(prn, _("Null hypothesis: the population proportions are "
		       "equal\n"));
	pputc(prn, '\n');
	pprintf(prn, _("Sample 1:\n n = %d, proportion = %g\n"), n1, x[0]);
	pputc(prn, '\n');
	pprintf(prn, _("Sample 2:\n n = %d, proportion = %g\n"), n2, x[1]);
	pputc(prn, '\n');

	x[2] = (n1*x[0] + n2*x[1]) / (n1 + n2);
	sderr = sqrt((x[2] * (1.0-x[2])) * (1.0/n1 + 1.0/n2));
	ts = (x[0] - x[1]) / sderr;

	pprintf(prn, _("Test statistic: z = (%g - %g) / %g = %g\n"),
		x[0], x[1], sderr, ts);

	pv = normal_pvalue_2(ts);
	print_pv(prn, pv, pv / 2.0);
	if (grf) {
	    htest_graph(0, ts, 0, 0);
	}
	break;

    default:
	break;
    }

    view_buffer(prn, 78, 340, _("gretl: hypothesis test"), H_TEST,
                NULL);
}

static void h_test_global (GtkWidget *w, CalcChild *child)
{
    test_t **test = child->calcp;
    int i;

    i = gtk_notebook_get_current_page(GTK_NOTEBOOK(child->book));
    h_test(NULL, test[i]);
}

static void np_test_global (GtkWidget *w, CalcChild *child)
{
    test_t **test = child->calcp;
    int i;

    i = gtk_notebook_get_current_page(GTK_NOTEBOOK(child->book));
    np_test(NULL, test[i]);
}

static void dist_graph (GtkWidget *w, CalcChild *child)
{
    lookup_t **look = child->calcp;
    int m = 0, n = 0;
    int d, err = 0;

    d = gtk_notebook_get_current_page(GTK_NOTEBOOK(child->book));

    switch (d) {
    case NORMAL_DIST:
	break;
    case T_DIST:
    case CHISQ_DIST:
	m = atoi(gtk_entry_get_text(GTK_ENTRY(look[d]->entry[0])));
	if (m <= 0) {
	    err = 1;
	}
	break;
    case F_DIST:
	m = atoi(gtk_entry_get_text(GTK_ENTRY(look[d]->entry[0])));
	n = atoi(gtk_entry_get_text(GTK_ENTRY(look[d]->entry[1])));
	if (m <= 0 || n <= 0) {
	    err = 1;
	}
	break;
    }

    if (err) {
	errbox(_("Invalid degrees of freedom"));
	return;
    }

    if (child->plot != NULL) {
	revise_distribution_plotspec(child->plot, d, m, n);
    } else {
	htest_graph(d, NADBL, m, n);
    } 
}

static void add_lookup_entry (GtkWidget *tbl, gint *tbl_len, 
			      const gchar *label, CalcChild *child,
			      int dist)
{
    lookup_t **look = child->calcp;
    int c = child->code;
    GtkWidget *tmp;

    *tbl_len += 1;

    /* label */
    gtk_table_resize(GTK_TABLE(tbl), *tbl_len, 2);
    tmp = gtk_label_new(_(label));
    gtk_misc_set_alignment(GTK_MISC(tmp), 1, 0.5);
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      tmp, 0, 1, *tbl_len - 1, *tbl_len);
    gtk_widget_show(tmp);

    /* numeric entry */
    tmp = gtk_entry_new();
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      tmp, 1, 2, *tbl_len - 1, *tbl_len);
    gtk_widget_show(tmp);
    look[dist]->entry[*tbl_len - 2] = tmp;

    g_signal_connect(G_OBJECT(tmp), "activate", 
		     (c == CALC_PVAL)? G_CALLBACK(get_pvalue) : 
		     (c == CALC_GRAPH || c == CALC_GRAPH_ADD)? 
		     G_CALLBACK(dist_graph) :
		     G_CALLBACK(get_critical),
		     child);
}

static void make_lookup_tab (CalcChild *child, int d)
{
    GtkWidget *tmp, *box, *tbl;
    gint tbl_len;
    const gchar *titles[] = {
	N_("normal"), 
	N_(" t "), 
	N_("chi-square"), 
	N_(" F "), 
	N_(" DW ")
    };
   
    box = gtk_vbox_new(FALSE, 0);
    gtk_container_set_border_width(GTK_CONTAINER(box), 10);
    gtk_widget_show(box);

    tmp = gtk_label_new(_(titles[d]));
    gtk_widget_show(tmp);
    gtk_notebook_append_page(GTK_NOTEBOOK(child->book), box, tmp);   

    tbl_len = 1;
    tbl = gtk_table_new(tbl_len, 2, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE (tbl), 5);
    gtk_table_set_col_spacings(GTK_TABLE (tbl), 5);
    gtk_box_pack_start(GTK_BOX(box), tbl, FALSE, FALSE, 0);
    gtk_widget_show(tbl);

    if (child->code == CALC_DIST && d != DW_DIST) {
	add_lookup_entry(tbl, &tbl_len, N_("right-tail probability"), child, d);
    }
   
    switch (d) {
    case NORMAL_DIST:
	break;
    case T_DIST:
	add_lookup_entry(tbl, &tbl_len, N_("df"), child, d);
	break;
    case CHISQ_DIST:
	add_lookup_entry(tbl, &tbl_len, N_("df"), child, d);
	break;
    case F_DIST:
	add_lookup_entry(tbl, &tbl_len, N_("dfn"), child, d);
	add_lookup_entry(tbl, &tbl_len, N_("dfd"), child, d);
	break;	
    case DW_DIST:
	add_lookup_entry(tbl, &tbl_len, N_("n"), child, d);
	break;
    default:
	break;
    } 
}

static void make_dist_tab (CalcChild *child, int d) 
{
    lookup_t **pval = child->calcp;
    GtkWidget *tempwid, *box, *tbl;
    gint tbl_len;
    const gchar *titles[] = {
	N_("normal"), 
	N_(" t "), 
	N_("chi-square"), 
	N_(" F "), 
	N_("gamma"),
	N_("binomial"),
	N_("poisson")
    };
   
    box = gtk_vbox_new(FALSE, 0);
    gtk_container_set_border_width(GTK_CONTAINER(box), 10);
    gtk_widget_show(box);

    tempwid = gtk_label_new(_(titles[d]));
    gtk_widget_show(tempwid);
    gtk_notebook_append_page(GTK_NOTEBOOK(child->book), box, tempwid);   

    tbl_len = 1;
    tbl = gtk_table_new(tbl_len, 2, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(tbl), 5);
    gtk_table_set_col_spacings(GTK_TABLE(tbl), 5);
    gtk_box_pack_start(GTK_BOX(box), tbl, FALSE, FALSE, 0);
    gtk_widget_show(tbl);
   
    switch (d) {

    case NORMAL_PVAL: 
	add_lookup_entry(tbl, &tbl_len, N_("value"), child, d);
	add_lookup_entry(tbl, &tbl_len, N_("mean"), child, d);
	gtk_entry_set_text(GTK_ENTRY(pval[0]->entry[1]), "0");
	add_lookup_entry(tbl, &tbl_len, N_("std. deviation"), child, d);
	gtk_entry_set_text(GTK_ENTRY(pval[0]->entry[2]), "1");
	break;

    case T_PVAL:
	add_lookup_entry(tbl, &tbl_len, N_("df"), child, d);
	add_lookup_entry(tbl, &tbl_len, N_("value"), child, d);
	add_lookup_entry(tbl, &tbl_len, N_("mean"), child, d);
	gtk_entry_set_text(GTK_ENTRY(pval[1]->entry[2]), "0");
	add_lookup_entry(tbl, &tbl_len, N_("std. deviation"), child, d);
	gtk_entry_set_text(GTK_ENTRY(pval[1]->entry[3]), "1");
	break;

    case CHISQ_PVAL:
	add_lookup_entry(tbl, &tbl_len, N_("df"), child, d);
	add_lookup_entry(tbl, &tbl_len, N_("x-value"), child, d);
	break;

    case F_PVAL:
	add_lookup_entry(tbl, &tbl_len, N_("dfn"), child, d);
	add_lookup_entry(tbl, &tbl_len, N_("dfd"), child, d);
	add_lookup_entry(tbl, &tbl_len, N_("F-value"), child, d);
	break;

    case GAMMA_PVAL:
	add_lookup_entry(tbl, &tbl_len, N_("mean"), child, d);
	add_lookup_entry(tbl, &tbl_len, N_("variance"), child, d);
	add_lookup_entry(tbl, &tbl_len, N_("x-value"), child, d);
	break;

    case BINOMIAL_PVAL:
	add_lookup_entry(tbl, &tbl_len, N_("Prob"), child, d);
	add_lookup_entry(tbl, &tbl_len, N_("trials"), child, d);
	add_lookup_entry(tbl, &tbl_len, N_("x-value"), child, d);
	break;

    case POISSON_PVAL:
	add_lookup_entry(tbl, &tbl_len, N_("mean"), child, d);
	add_lookup_entry(tbl, &tbl_len, N_("x-value"), child, d);
	break;

    default:
	break;
    } 
}

static int get_restriction_vxy (const char *s, int *vx, int *vy, 
				GretlOp *yop, double *yval)
{
    char test[VNAMELEN];
    char *p, *q = NULL;
    char *str = g_strdup(s);
    int err = 0;

    if (str == NULL) {
	return 1;
    }

    p = strchr(str, '(');
    *p = 0;
    p++;

    if (sscanf(str, "%15s", test) != 1) {
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
	if (sscanf(p, "%15s", test) != 1) {
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
   choice of variable (or variable plus restriction) */

static void populate_stats (GtkWidget *w, gpointer p)
{
    test_t *test = g_object_get_data(G_OBJECT(p), "test");
    int pos = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(p), "pos"));
    gchar **pbuf = g_object_get_data(G_OBJECT(p), "pbuf");
    int t, n = datainfo->t2 - datainfo->t1 + 1;
    int vx = -1, vy = -1;
    GretlOp yop = 0;
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
    if (var_is_scalar(datainfo, vx) || (vy > 0 && var_is_scalar(datainfo, vy))) {
	errbox(_("Invalid entry"));
	return;
    }

    for (t=datainfo->t1; t<=datainfo->t2; t++) {
	if (na(Z[vx][t]) || (vy > 0 && !eval_ytest(Z[vy][t], yop, yval))) {
	    n--;
	}
    }

    if (n == 0) {		
	errbox(_("Data missing for variable '%s'"), datainfo->varname[vx]);
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
    } else if (test->code == ONE_PROPN || test->code == TWO_PROPNS) {
	if (vy < 0) {
	    x1 = gretl_mean(datainfo->t1, datainfo->t2, Z[vx]);
	} else {
	    x1 = gretl_restricted_mean(datainfo->t1, datainfo->t2, 
				      Z[vx], Z[vy], yop, yval);
	}
	sprintf(numstr, "%.10g", x1);
	gtk_entry_set_text(GTK_ENTRY(test->entry[pos]), numstr);
	sprintf(numstr, "%d", n);
	gtk_entry_set_text(GTK_ENTRY(test->entry[pos + 1]), numstr);
    }	
}

static int var_is_ok (int i, int code)
{
    int ret = 1;

    if (var_is_hidden(datainfo, i)) {
	ret = 0;
    } else if (var_is_scalar(datainfo, i)) {
	ret = 0;
    } else if ((code == ONE_PROPN || code == TWO_PROPNS) &&
	       !gretl_isdummy(datainfo->t1, datainfo->t2, Z[i])) {
	ret = 0;
    }

    return ret;
}

static void add_vars_to_combo (GtkWidget *w, int code, int pos)
{
    GList *vlist = NULL;
    int i, vmin = (pos > 0)? 2 : 1;

    for (i=vmin; i<datainfo->v; i++) {
	if (var_is_ok(i, code)) {
	    vlist = g_list_append(vlist, datainfo->varname[i]);
	}
    }

    if (pos > 0) {
	/* add first variable at the end of the list */
	for (i=1; i<datainfo->v; i++) {
	    if (var_is_ok(i, code)) {
		vlist = g_list_append(vlist, datainfo->varname[i]);
		break;
	    }
	}
    }	

    gtk_combo_set_popdown_strings(GTK_COMBO(w), vlist);
}

static void toggle_combo_ok (GtkWidget *toggle, gpointer p)
{
    if (GTK_TOGGLE_BUTTON(toggle)->active) {
	gtk_widget_set_sensitive(GTK_WIDGET(p), TRUE);
	populate_stats(NULL, p);
    } else {
	gtk_widget_set_sensitive(GTK_WIDGET(p), FALSE);
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

static void add_test_var_selector (GtkWidget *tbl, gint *tbl_len, 
				   test_t *test, int pos,
				   int labelit)
{
    GtkWidget *label, *tmp;
    gchar **pbuf;

    *tbl_len += 1;
    gtk_table_resize(GTK_TABLE(tbl), *tbl_len, 2);
    if (labelit) {
	gchar *tmp = g_strdup_printf(_("Variable %d"), pos + 1);
	label = gtk_label_new(tmp);
	g_free(tmp);
    } else {
	label = gtk_label_new(_("Variable"));
    }
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      label, 0, 1, *tbl_len - 1, *tbl_len);
    gtk_widget_show(label);

    tmp = gtk_combo_new();
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      tmp, 1, 2, *tbl_len - 1, *tbl_len);
    gtk_widget_show(tmp);
    g_object_set_data(G_OBJECT(tmp), "test", test);
    g_object_set_data(G_OBJECT(tmp), "pos", GINT_TO_POINTER(pos));
    test->entry[pos] = GTK_COMBO(tmp)->entry;

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

    add_vars_to_combo(tmp, test->code, pos);
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

    add_vars_to_combo(tmp, test->code, pos);
    gtk_widget_set_sensitive(tmp, FALSE);
    gtk_combo_disable_activate(GTK_COMBO(tmp));

    g_signal_connect(G_OBJECT(GTK_ENTRY(GTK_COMBO(tmp)->entry)), "key_press_event",
		     G_CALLBACK(catch_combo_key), tmp);
    g_signal_connect(G_OBJECT(GTK_LIST(GTK_COMBO(tmp)->list)), "select-child",
		     G_CALLBACK(select_child_callback), tmp);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(toggle_combo_ok), tmp);
}

static void add_test_entry (GtkWidget *tbl, gint *tbl_len, 
			    const gchar *label, test_t *test, 
			    int i)
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
			    const gchar *label, test_t *test,
			    gboolean val)
{
    GtkWidget *tmp;

    *tbl_len += 1;
    gtk_table_resize(GTK_TABLE(tbl), *tbl_len, 2);
    tmp = gtk_check_button_new_with_label(label);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), val);
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      tmp, 0, 2, *tbl_len - 1, *tbl_len);
    gtk_widget_show(tmp);
    test->check = tmp;
}

static int n_ok_series (void)
{
    int i, nv = 0;

    if (datainfo != NULL) {
	for (i=1; i<datainfo->v; i++) {
	    if (var_is_series(datainfo, i) && !var_is_hidden(datainfo, i)) {
		nv++;
	    }
	}
    }

    return nv;
}

static int n_ok_dummies (void)
{
    int i, nv = 0;

    if (datainfo != NULL) {
	for (i=1; i<datainfo->v; i++) {
	    if (var_is_series(datainfo, i) && 
		!var_is_hidden(datainfo, i) &&
		gretl_isdummy(datainfo->t1, datainfo->t2, Z[i])) {
		nv++;
	    }
	}
    }

    return nv;
}

static gint toggle_verbose_state (GtkWidget *w, GtkWidget *b)
{
    if (GTK_TOGGLE_BUTTON(w)->active) {
	gtk_widget_set_sensitive(b, FALSE);
    } else {
	gtk_widget_set_sensitive(b, TRUE);
    }

    return FALSE;
}

static void make_nptest_tab (CalcChild *child, int idx) 
{
    test_t **tests = child->calcp;
    test_t *test = tests[idx];
    GtkWidget *tmp, *box, *tbl;
    GSList *group;
    gint i, tbl_len;
    const gchar *titles[] = {
	N_("Difference test"),
	N_("Runs test")
    };

    box = gtk_vbox_new(FALSE, 0);
    gtk_container_set_border_width(GTK_CONTAINER(box), 10);
    gtk_widget_show(box);

    tmp = gtk_label_new(_(titles[idx]));
    gtk_widget_show(tmp);
    gtk_notebook_append_page(GTK_NOTEBOOK(child->book), box, tmp);   

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

    switch (idx) {

    case NP_DIFF: 
	add_test_var_selector(tbl, &tbl_len, test, 0, 1);
	add_test_var_selector(tbl, &tbl_len, test, 1, 1);

	/* option radios */
	tbl_len += 3;
	gtk_table_resize(GTK_TABLE(tbl), tbl_len, 2);

	test->radio[0] = gtk_radio_button_new_with_label(NULL, _("Sign test"));
	gtk_table_attach_defaults(GTK_TABLE(tbl), test->radio[0], 0, 2, 
				  tbl_len - 3, tbl_len - 2);
	gtk_widget_show(test->radio[0]);

	group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(test->radio[0]));
	test->radio[1] = gtk_radio_button_new_with_label(group, _("Wilcoxon rank sum test"));
	gtk_table_attach_defaults(GTK_TABLE(tbl), test->radio[1], 0, 2, 
				  tbl_len - 2, tbl_len - 1);
	gtk_widget_show(test->radio[1]);

	group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(test->radio[1]));
	test->radio[2] = gtk_radio_button_new_with_label(group, _("Wilcoxon signed rank test"));
	gtk_table_attach_defaults(GTK_TABLE(tbl), test->radio[2], 0, 2, 
				  tbl_len - 1, tbl_len);
	gtk_widget_show(test->radio[2]);
	break;

    case NP_RUNS: 
	add_test_var_selector(tbl, &tbl_len, test, 0, 0);
	break;

    default:
	break;
    } 

    /* check box for extra option */
    tbl_len += 1;
    gtk_table_resize(GTK_TABLE(tbl), tbl_len, 2);
    if (idx == NP_DIFF) {
	tmp = gtk_check_button_new_with_label(_("Show details"));
    } else {
	tmp = gtk_check_button_new_with_label(_("Use first difference"));
    }
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 2, 
			      tbl_len - 1, tbl_len);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), FALSE);

    if (idx == NP_DIFF) {
	gtk_widget_set_sensitive(tmp, FALSE);
	g_signal_connect(G_OBJECT(test->radio[0]), "toggled",
			 G_CALLBACK(toggle_verbose_state),
			 tmp);
    }

    gtk_widget_show(tmp);
    test->extra = tmp; 
}

static void make_test_tab (CalcChild *child, int idx) 
{
    test_t **tests = child->calcp;
    test_t *test = tests[idx];
    GtkWidget *tempwid, *box, *tbl;
    int nv = 0;
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

    tempwid = gtk_label_new(_(titles[idx]));
    gtk_widget_show(tempwid);
    gtk_notebook_append_page(GTK_NOTEBOOK(child->book), box, tempwid);   

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

    if (idx == ONE_MEAN || idx == ONE_VARIANCE) {
	nv = n_ok_series();
    } else if (idx == TWO_MEANS || idx == TWO_VARIANCES) {
	nv = (n_ok_series() > 1);
    } else if (idx == ONE_PROPN) {
	nv = n_ok_dummies();
    } else if (idx == TWO_PROPNS) {
	nv = (n_ok_dummies() > 1);
    }

    if (nv > 0) {
	add_test_combo(tbl, &tbl_len, test, 0);
    }
   
    switch (idx) {

    case ONE_MEAN: 
	add_test_entry(tbl, &tbl_len, _("sample mean"), test, 0);
	add_test_entry(tbl, &tbl_len, _("std. deviation"), test, 1);
	add_test_entry(tbl, &tbl_len, _("sample size"), test, 2);
	add_test_entry(tbl, &tbl_len, _("H0: mean ="), test, 3);
	add_test_check(tbl, &tbl_len, _("Assume standard deviation is "
		       "population value"), test, FALSE);
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
	if (nv > 0) {
	    add_test_combo(tbl, &tbl_len, test, 3);
	}
	add_test_entry(tbl, &tbl_len, _("mean of sample 2"), test, 3);
	add_test_entry(tbl, &tbl_len, _("std. deviation, sample 2"), test, 4);
	add_test_entry(tbl, &tbl_len, _("size of sample 2"), test, 5);
	add_test_entry(tbl, &tbl_len, _("H0: Difference of means ="), test, 6);
	gtk_entry_set_text(GTK_ENTRY(test->entry[6]), "0");
	add_test_check(tbl, &tbl_len, _("Assume common population standard "
		       "deviation"), test, TRUE);
	break;

    case TWO_VARIANCES:
	add_test_entry(tbl, &tbl_len, _("variance of sample 1"), test, 0);
	add_test_entry(tbl, &tbl_len, _("size of sample 1"), test, 1);
	if (nv > 0) {
	    add_test_combo(tbl, &tbl_len, test, 2);
	}	
	add_test_entry(tbl, &tbl_len, _("variance of sample 2"), test, 2);
	add_test_entry(tbl, &tbl_len, _("size of sample 2"), test, 3);
	add_test_label(tbl, &tbl_len, _("H0: Ratio of variances = 1"));
	break;

    case TWO_PROPNS:
	add_test_entry(tbl, &tbl_len, _("proportion, sample 1"), test, 0);
	add_test_entry(tbl, &tbl_len, _("size of sample 1"), test, 1);
	if (nv > 0) {
	    add_test_combo(tbl, &tbl_len, test, 2);
	}	
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

    test->extra = tempwid; 
}

static void gretl_child_destroy (GtkWidget *w, CalcChild *child)
{
    GtkWidget **wp = (GtkWidget **) child->winp;
    int c = child->code;
    int i, n;

    *wp = NULL;

    if (c == CALC_TEST) {
	test_t **test = child->calcp;

	for (i=0; i<NTESTS; i++) {
	    free(test[i]);
	}
	free(test);
    } else if (c == CALC_NPTEST) {
	test_t **test = child->calcp;

	for (i=0; i<NPTESTS; i++) {
	    free(test[i]);
	}
	free(test);
    } else {	
	lookup_t **look = child->calcp;

	n = (c == CALC_DIST)? NLOOKUPS :
	    (c == CALC_PVAL)? NPVAL : NGRAPHS;

	for (i=0; i<n; i++) {
	    free(look[i]);
	}
	free(look);	
    } 

    free(child);
}

static int child_allocate_calcp (CalcChild *child)
{
    int c = child->code;
    test_t **test = NULL;
    lookup_t **look = NULL;
    int i, n, err = 0;

    child->calcp = NULL;
    
    if (c == CALC_TEST) {
	test = mymalloc(NTESTS * sizeof *test);
	if (test != NULL) {
	    child->calcp = test;
	    for (i=0; i<NTESTS && !err; i++) {
		test[i] = mymalloc(sizeof **test);
		if (test[i] == NULL) {
		    err = E_ALLOC;
		} else {
		    test[i]->code = i;
		}
	    }
	} 
    } else if (c == CALC_NPTEST) {
	test = mymalloc(NPTESTS * sizeof *test);
	if (test != NULL) {
	    child->calcp = test;
	    for (i=0; i<NPTESTS && !err; i++) {
		test[i] = mymalloc(sizeof **test);
		if (test[i] == NULL) {
		    err = E_ALLOC;
		} else {
		    test[i]->code = i;
		}
	    }
	} 
    } else {
	n = (c == CALC_PVAL)? NPVAL : 
	    (c == CALC_DIST)? NLOOKUPS :
	    NGRAPHS;
	look = mymalloc(n * sizeof *look);
	if (look != NULL) {
	    child->calcp = look;
	    for (i=0; i<n && !err; i++) {
		look[i] = mymalloc(sizeof **look);
		if (look[i] == NULL) {
		    err = E_ALLOC;
		}
	    }
	} 
    }

    if (child->calcp == NULL) {
	err = E_ALLOC;
    }

    return err;
}

static CalcChild *gretl_child_new (int code, gpointer p)
{
    CalcChild *child;
    GtkWidget *base;

    child = mymalloc(sizeof *child);
    if (child == NULL) return NULL;
    
    child->code = code;
    child->plot = (code == CALC_GRAPH_ADD)? p : NULL;

    if (child_allocate_calcp(child)) {
	free(child);
	return NULL;
    }

    child->dlg = gtk_window_new(GTK_WINDOW_TOPLEVEL);

    g_object_set_data(G_OBJECT(child->dlg), "gchild", child);

    base = gtk_vbox_new(FALSE, 0);
    gtk_container_add(GTK_CONTAINER(child->dlg), base);
    gtk_widget_show(base);

    child->vbox = gtk_vbox_new(FALSE, 5);
    gtk_container_add(GTK_CONTAINER(base), child->vbox);
    gtk_container_set_border_width(GTK_CONTAINER(child->vbox), 5);
    gtk_widget_show(child->vbox);

    child->book = gtk_notebook_new();
    gtk_box_pack_start(GTK_BOX(child->vbox), child->book, TRUE, TRUE, 0);
    gtk_widget_show(child->book);

    child->bbox = gtk_hbutton_box_new();
    gtk_button_box_set_layout(GTK_BUTTON_BOX(child->bbox), 
			      GTK_BUTTONBOX_END);
    gtk_button_box_set_spacing(GTK_BUTTON_BOX(child->bbox),
			       10);
    gtk_widget_show(child->bbox);
    gtk_container_add(GTK_CONTAINER(base), child->bbox);
    gtk_container_set_border_width(GTK_CONTAINER(child->bbox), 5);

    gtk_window_set_position(GTK_WINDOW(child->dlg), GTK_WIN_POS_MOUSE);

    g_signal_connect(G_OBJECT(child->dlg), "destroy",
                     G_CALLBACK(gretl_child_destroy),
                     child);

    return child;
}

static void 
make_graph_window_transient (GtkWidget *win, png_plot *plot)
{
    GtkWidget *pshell = plot_get_shell(plot);

    gtk_window_set_transient_for(GTK_WINDOW(win), GTK_WINDOW(pshell));
    gtk_window_set_destroy_with_parent(GTK_WINDOW(win), TRUE);
}

static void switch_child_role (GtkWidget *win, png_plot *plot)
{
    CalcChild *child;

    child = g_object_get_data(G_OBJECT(win), "gchild");
    child->plot = plot;
    child->code = CALC_GRAPH_ADD;

    gtk_window_set_title(GTK_WINDOW(win), 
			 _("gretl: add distribution graph"));

    make_graph_window_transient(win, plot);

    if (gtk_notebook_get_current_page(GTK_NOTEBOOK(child->book)) ==
	NORMAL_DIST) {
	gtk_notebook_set_current_page(GTK_NOTEBOOK(child->book),
				      T_DIST);
    }
}

void stats_calculator (gpointer data, guint code, GtkWidget *widget) 
{
    GtkWidget *tmp = NULL;
    static GtkWidget *winptr[5];
    GtkWidget *oldwin;
    CalcChild *child;
    const char *calc_titles[] = {
	N_("gretl: p-value finder"),
	N_("gretl: statistical tables"),
	N_("gretl: test calculator"),
	N_("gretl: nonparametric tests"),
	N_("gretl: distribution graphs"),
	N_("gretl: add distribution graph"),
    };
    int i, nv = 0;

    g_return_if_fail(code == CALC_PVAL || 
		     code == CALC_DIST || 
		     code == CALC_TEST ||
		     code == CALC_NPTEST ||
		     code == CALC_GRAPH ||
		     code == CALC_GRAPH_ADD);

    oldwin = winptr[code];
    if (oldwin != NULL) {
	gtk_window_present(GTK_WINDOW(oldwin));
 	return;
    }

    if (code == CALC_GRAPH_ADD) {
	oldwin = winptr[CALC_GRAPH];
	if (oldwin != NULL) {
	    switch_child_role(oldwin, data);
	    gtk_window_present(GTK_WINDOW(oldwin));
	    return;
	}
    }

    if (code == CALC_NPTEST) {
	nv = n_ok_series();
	if (nv == 0) {
	    errbox(_("No suitable data are available"));
	    return;
	}
    }

    child = gretl_child_new(code, data);
    winptr[code] = child->dlg;
    child->winp = &(winptr[code]);

    gtk_window_set_title(GTK_WINDOW(child->dlg), _(calc_titles[code]));

    if (code == CALC_GRAPH_ADD) {
	make_graph_window_transient(child->dlg, data);
    }

    if (code == CALC_TEST) {
	for (i=0; i<NTESTS; i++) {
	    make_test_tab(child, i);
	}
    } else if (code == CALC_NPTEST) {
	for (i=0; i<NPTESTS; i++) {
	    make_nptest_tab(child, i);
	}
    } else if (code == CALC_PVAL) {
	for (i=0; i<NPVAL; i++) {
	    make_dist_tab(child, i);
	}	
    } else if (code == CALC_DIST) {
	for (i=0; i<NLOOKUPS; i++) {
	    make_lookup_tab(child, i);
	}
    } else {
	for (i=0; i<NGRAPHS; i++) {
	    make_lookup_tab(child, i);
	}
    }

    if (code == CALC_GRAPH || code == CALC_GRAPH_ADD) {
	gtk_notebook_set_current_page(GTK_NOTEBOOK(child->book), T_DIST);
    }

    if (code == CALC_NPTEST && nv < 2) {
	GtkWidget *p = gtk_notebook_get_nth_page(GTK_NOTEBOOK(child->book), NP_DIFF);

	gtk_widget_set_sensitive(p, FALSE);
	gtk_notebook_set_current_page(GTK_NOTEBOOK(child->book), NP_RUNS);
	
    }

    /* Close button */
    tmp = gtk_button_new_from_stock(GTK_STOCK_CLOSE);
    GTK_WIDGET_SET_FLAGS(tmp, GTK_CAN_DEFAULT);
    gtk_container_add(GTK_CONTAINER(child->bbox), tmp);
    g_signal_connect(G_OBJECT(tmp), "clicked", 
		     G_CALLBACK(delete_widget), 
		     child->dlg);
    gtk_widget_show(tmp);

    /* OK button */
    tmp = gtk_button_new_from_stock(GTK_STOCK_OK);
    GTK_WIDGET_SET_FLAGS(tmp, GTK_CAN_DEFAULT);
    gtk_container_add(GTK_CONTAINER(child->bbox), tmp);
    g_signal_connect(G_OBJECT (tmp), "clicked", 
		     (code == CALC_PVAL)? G_CALLBACK(get_pvalue) :
		     (code == CALC_DIST)? G_CALLBACK(get_critical) :
		     (code == CALC_NPTEST)? G_CALLBACK(np_test_global) :
		     (code == CALC_GRAPH || code == CALC_GRAPH_ADD)? 
		     G_CALLBACK(dist_graph) :
		     G_CALLBACK(h_test_global),
		     child);
    gtk_widget_show(tmp);

    /* Help button? */
    if (code == CALC_TEST || code == CALC_NPTEST) { 
	int hcode = (code == CALC_TEST)? HTEST : HTESTNP;

	tmp = gtk_button_new_from_stock(GTK_STOCK_HELP);
	GTK_WIDGET_SET_FLAGS(tmp, GTK_CAN_DEFAULT);
	gtk_container_add(GTK_CONTAINER(child->bbox), tmp);
	gtk_button_box_set_child_secondary(GTK_BUTTON_BOX(child->bbox),
					   tmp, TRUE);
	g_signal_connect(G_OBJECT(tmp), "clicked", 
			 G_CALLBACK(context_help), 
			 GINT_TO_POINTER(hcode));
	gtk_widget_show(tmp);
    }

    gtk_widget_show(child->dlg);
}





