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

/* calculator.c for gretl */

#define NTESTS 6
#define NPTESTS 2
#define NPVAL 7
#define NLOOKUPS 7
#define NGRAPHS 4
#define NRAND 8
#define NTESTENTRY 7
#define NLOOKUPENTRY 4

#include "gretl.h"
#include "calculator.h"
#include "dlgutils.h"
#include "gpt_control.h"
#include "lib_private.h"
#include "cmdstack.h"

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
    CALC_NORMAL,
    CALC_STUDENT,
    CALC_CHISQ,
    CALC_SNEDECOR,
    CALC_GAMMA,
    CALC_BINOMIAL,
    CALC_POISSON
};

enum {
    RAND_UNIFORM,
    RAND_NORMAL,
    RAND_STUDENT,
    RAND_CHISQ,
    RAND_SNEDECOR,
    RAND_GAMMA,
    RAND_BINOMIAL,
    RAND_POISSON
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

enum {
    C_INT,
    C_DBL,
    C_POS_DBL
};

/* getval: get a numerical value from a text entry box

   t == C_DBL: parse the entry as a double
   t == C_POS_DBL: parse as a positive double
   t == C_INT: parse as positive int

   (In this context, when an int is wanted, it's always a 
   positive int: df, sample size, number of trials.)

   If t is C_DBL or C_POS_DBL and something is wrong with the
   input, we flag this by returning NADBL.  In the case of
   C_INT, we flag an error by returning -1.

   Possible refinement: should we accept names of variables
   in these entry fields?
*/

static double getval (GtkWidget *w, int t)
{
    const gchar *s = gtk_entry_get_text(GTK_ENTRY(w));
    double x = NADBL;
    int k;

    if (s == NULL || *s == '\0') {
	errbox(_("Incomplete entry"));
    } else if (t == C_INT) {
	if (check_atoi(s)) {
	    errbox(gretl_errmsg_get());
	} else {
	    k = atoi(s);
	    if (k <= 0) {
		errbox(_("Invalid entry"));
		x = -1;
	    } else {
		x = k;
	    }
	}
    } else {
	if (check_atof(s)) {
	    errbox(gretl_errmsg_get());
	} else {
	    x = atof(s);
	    if (t == C_POS_DBL && !na(x) && x <= 0.0) {
		errbox(_("Invalid entry"));
		x = NADBL;
	    }
	}
    }

    return x;
}

static int check_prob (double p)
{
    int err = 0;

    if (na(p)) {
	err = 1;
    } else if (p >= 1.0) {
	errbox(_("Invalid probability"));
	err = 1;
    }

    return err;
}

/* call plugin function to look up part of the table of
   critical values for the Durbin-Watson statistic
*/

static void dw_lookup (lookup_t *tab)
{
    void *handle = NULL;
    void (*dw)(int, PRN *) = NULL;
    PRN *prn;
    int n;

    n = getval(tab->entry[0], C_INT);
    if (n < 0) {
	return;
    }

    dw = gui_get_plugin_function("dw_lookup", &handle);
    if (dw == NULL) {
	return;
    }

    if (bufopen(&prn)) {
	close_plugin(handle);
	return;
    }  

    (*dw)(n, prn);
    close_plugin(handle);

    view_buffer(prn, 77, 300, _("gretl: critical values"), 
		STAT_TABLE, NULL);
}

static void get_critical (GtkWidget *w, CalcChild *child)
{
    lookup_t *tab, **tabs = child->calcp;
    double c = NADBL;
    double x[4];
    char st = 0;
    int d, j = 0;
    PRN *prn;

    d = gtk_notebook_get_current_page(GTK_NOTEBOOK(child->book));
    tab = tabs[d];
    
    if (d == DW_DIST) {
	/* special: just a table look-up */
	dw_lookup(tab);
	return;
    }

    switch (d) {
    case NORMAL_DIST:
	st = 'z';
	break;
    case T_DIST:
    case CHISQ_DIST:
	st = (d == T_DIST)? 't' : 'X';
	x[j] = getval(tab->entry[j], C_INT); /* df */
	if (x[j++] < 0) return;
	break;
    case F_DIST:
	st = 'F';
	x[j] = getval(tab->entry[j], C_INT); /* dfn */
	if (x[j++] < 0) return;
	x[j] = getval(tab->entry[j], C_INT); /* dfd */
	if (x[j++] < 0) return;
	break;
    case BINOMIAL_DIST:
	st = 'B';
	x[j] = getval(tab->entry[j], C_POS_DBL); /* prob */
	if (check_prob(x[j++])) return;
	x[j] = getval(tab->entry[j], C_INT); /* n */
	if (x[j++] < 0) return;
	break;
    case POISSON_DIST:
	st = 'P';
	x[j] = getval(tab->entry[j], C_POS_DBL); /* mean */
	if (na(x[j++])) return;
	break;
    default:
	break;
    }

    /* right-tail probability */
    x[j] = getval(tab->entry[j], C_POS_DBL);
    if (check_prob(x[j++])) return;

    c = gretl_get_critval(st, x);
    if (na(c)) {
	errbox(_("Failed to compute critical value"));
	return;
    }

    if (bufopen(&prn)) {
	return;
    }   

    x[j] = c;

    print_critval(st, x, prn);
    view_buffer(prn, 60, 200, _("gretl: critical values"), 
		STAT_TABLE, NULL);
}

static void get_pvalue (GtkWidget *w, CalcChild *child)
{
    lookup_t *tab, **tabs = child->calcp;
    double pv, x[3];
    char st = 0;
    int i, j = 0;
    PRN *prn;

    i = gtk_notebook_get_current_page(GTK_NOTEBOOK(child->book));
    tab = tabs[i];

    switch (i) {

    case CALC_NORMAL:
	st = 'z';
	x[j] = getval(tab->entry[j], C_DBL); /* mean */
	if (na(x[j++])) return;
	x[j] = getval(tab->entry[j], C_POS_DBL); /* s.d. */
	if (na(x[j++])) return;
	x[j] = getval(tab->entry[j], C_DBL); /* val */
	if (na(x[j])) return;
	x[0] = (x[2] - x[0]) / x[1]; /* z-score */
	break;

    case CALC_STUDENT: 
	st = 't';
	x[j] = getval(tab->entry[j], C_INT); /* df */
	if (x[j++] < 0) return;
	x[j] = getval(tab->entry[j], C_DBL); /* val */
	if (na(x[j])) return;
	break;

    case CALC_CHISQ:
	st = 'X';
	x[j] = getval(tab->entry[j], C_INT); /* df */
	if (x[j++] < 0) return;
	x[j] = getval(tab->entry[j], C_POS_DBL); /* val */
	if (na(x[j])) return;
	break;

    case CALC_POISSON: 
	st = 'P';
	x[j] = getval(tab->entry[j], C_POS_DBL); /* mean */
	if (na(x[j++])) return;
	x[j] = getval(tab->entry[j], C_INT); /* val */
	if (x[j] < 0) return;
	break;

    case CALC_SNEDECOR:
	st = 'F';
	x[j] = getval(tab->entry[j], C_INT); /* dfn */
	if (x[j++] < 0) return;
	x[j] = getval(tab->entry[j], C_INT); /* dfd */
	if (x[j++] < 0) return;
	x[j] = getval(tab->entry[j], C_POS_DBL); /* val */
	if (na(x[j])) return;
	break;

    case CALC_GAMMA: 
	st = 'G';
	for (j=0; j<3; j++) {
	    x[j] = getval(tab->entry[j], C_POS_DBL);
	    if (na(x[j])) return;
	}
	break;

    case CALC_BINOMIAL: 
	st = 'B';
	x[j] = getval(tab->entry[j], C_POS_DBL); /* prob */
	if (check_prob(x[j++])) return;
	x[j] = getval(tab->entry[j], C_INT); /* trials */
	if (x[j++] < 0) return;
	x[j] = getval(tab->entry[j], C_INT); /* val */
	if (x[j] < 0) return;
	break;

    default:
	errbox(_("Failed to compute p-value"));
	return;
    }

    if (bufopen(&prn)) return;

    pv = gretl_get_pvalue(st, x);

    if (na(pv)) {
	errbox(_("Failed to compute p-value"));
    } else {
	print_pvalue(st, x, pv, prn);
	view_buffer(prn, 78, 200, _("gretl: p-value"), PVALUE, NULL);
    }
}

static int calc_finish_genr (void)
{
    char *cmdline = get_lib_cmdline();
    PRN *prn;
    int err = 0;

    if (bufopen(&prn)) {
	return 1;
    }

    err = generate(cmdline, &Z, datainfo, OPT_NONE, prn); 

    if (err) {
	errbox(gretl_print_get_buffer(prn));
	delete_last_command();
    } else {
	populate_varlist();
	mark_dataset_as_modified();
    }

    gretl_print_destroy(prn);

    return err;
}

static void get_random (GtkWidget *w, CalcChild *child)
{
    lookup_t *tab, **tabs = child->calcp;
    const char *vname;
    double x[2] = {0};
    int i, j = 0;

    i = gtk_notebook_get_current_page(GTK_NOTEBOOK(child->book));
    tab = tabs[i];

    switch (i) {

    case RAND_UNIFORM:
	x[j] = getval(tab->entry[j], C_DBL); /* min */
	if (na(x[j++])) return;
	x[j] = getval(tab->entry[j], C_DBL); /* max */
	if (na(x[j++])) return;
	break;

    case RAND_NORMAL:
	x[j] = getval(tab->entry[j], C_DBL); /* mean */
	if (na(x[j++])) return;
	x[j] = getval(tab->entry[j], C_POS_DBL); /* s.d. */
	if (na(x[j++])) return;
	break;

    case RAND_GAMMA:
	x[j] = getval(tab->entry[j], C_POS_DBL); /* shape */
	if (na(x[j++])) return;
	x[j] = getval(tab->entry[j], C_POS_DBL); /* scale */
	if (na(x[j++])) return;
	break;

    case RAND_STUDENT: 
    case RAND_CHISQ:
	x[j] = getval(tab->entry[j], C_INT); /* df */
	if (x[j++] < 0) return;
	break;

    case RAND_POISSON: 
	x[j] = getval(tab->entry[j], C_POS_DBL); /* mean */
	if (na(x[j++])) return;
	break;

    case RAND_SNEDECOR:
	x[j] = getval(tab->entry[j], C_INT); /* dfn */
	if (x[j++] < 0) return;
	x[j] = getval(tab->entry[j], C_INT); /* dfd */
	if (x[j++] < 0) return;
	break;

    case RAND_BINOMIAL: 
	x[j] = getval(tab->entry[j], C_POS_DBL); /* prob */
	if (check_prob(x[j++])) return;
	x[j] = getval(tab->entry[j], C_INT); /* trials */
	if (x[j++] < 0) return;
	break;

    default:
	return;
    }

    vname = gtk_entry_get_text(GTK_ENTRY(tab->entry[j]));
    if (vname == NULL || *vname == '\0') {
	errbox(_("You must give a name for the variable"));
	return;
    } else if (validate_varname(vname)) {
	return;
    }

    if (i == RAND_UNIFORM && x[0] >= x[1]) {
	errbox(_("Range is non-positive!"));
	return;
    }

    switch (i) {

    case RAND_UNIFORM:
	gretl_command_sprintf("genr %s = uniform(%g,%g)", vname,
			      x[0], x[1]);
	break;

    case RAND_NORMAL:
	gretl_command_sprintf("genr %s = normal(%g,%g)", vname,
			      x[0], x[1]);
	break;

    case RAND_STUDENT: 
	gretl_command_sprintf("genr %s = student(%g)", vname, x[0]);
	break;

    case RAND_CHISQ:
	gretl_command_sprintf("genr %s = chisq(%g)", vname, x[0]);
	break;

    case RAND_POISSON: 
	/* FIXME allow variable as param? */
	gretl_command_sprintf("genr %s = poisson(%g)", vname, x[0]);
	break;

    case RAND_SNEDECOR:
	gretl_command_sprintf("genr %s = randF(%g,%g)", vname, 
			      x[0], x[1]);
	break;

    case RAND_GAMMA:
	gretl_command_sprintf("genr %s = rgamma(%g,%g)", vname, 
			      x[0], x[1]);
	break;

    case RAND_BINOMIAL: 
	gretl_command_sprintf("genr %s = binomial(%g,%g)", vname, 
			      x[1], x[0]);
	break;
    }

    if (check_and_record_command()) {
	return;
    }

    calc_finish_genr();
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
		snedecor_critval(a, df1, df2);
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

/* non-parametric test */

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

static void do_h_test (test_t *test, double *x, int n1, int n2)
{
    double se, ts, pv, z;
    int j, grf;
    PRN *prn;

    if (bufopen(&prn)) {
	return;
    }

    grf = GTK_TOGGLE_BUTTON(test->extra)->active;

    switch (test->code) {

    case ONE_MEAN:
	se = x[1] / sqrt((double) n1);
	ts = (x[0] - x[2]) / se;

	pprintf(prn, _("Null hypothesis: population mean = %g\n"), x[2]);
	pprintf(prn, _("Sample size: n = %d\n"), n1);
	pprintf(prn, _("Sample mean = %g, std. deviation = %g\n"), 
		x[0], x[1]);

	if (GTK_TOGGLE_BUTTON(test->check)->active) {
	    pprintf(prn, _("Test statistic: z = (%g - %g)/%g = %g\n"), 
		    x[0], x[2], se, ts);
	    pv = normal_pvalue_2(ts);
	    print_pv(prn, pv, pv / 2.0);
	    if (grf) {
		htest_graph(0, ts, 0, 0);
	    }
	} else {
	    pprintf(prn, _("Test statistic: t(%d) = (%g - %g)/%g = %g\n"), n1-1,
		    x[0], x[2], se, ts);
	    pv = student_pvalue_2(ts, n1 - 1);
	    print_pv(prn, pv, 0.5 * pv);
	    if (grf) {
		htest_graph(1, ts, n1-1, 0);
	    }
	}
	break;

    case ONE_VARIANCE:
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

    case ONE_PROPN:
	se = sqrt(x[1] * (1.0 - x[1]) / n1);
	ts = (x[0] - x[1]) / se;

	pprintf(prn, _("Null hypothesis: population proportion = %g\n"), x[1]);
	pprintf(prn, _("Sample size: n = %d\n"), n1);
	pprintf(prn, _("Sample proportion = %g\n"), x[0]);
	pprintf(prn, _("Test statistic: z = (%g - %g)/%g = %g\n"), 
		x[0], x[1], se, ts);

	pv = normal_pvalue_2(ts);
	print_pv(prn, pv, pv / 2.0);
	if (grf) {
	    htest_graph(0, ts, 0, 0);
	}
	break;

    case TWO_MEANS:
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
	    se = sqrt(ts / n1 + ts / n2);
	} else {
	    double v1 = x[1] * x[1] / n1;
	    double v2 = x[3] * x[3] / n2;

	    se = sqrt(v1 + v2);
	}

	ts = (x[0] - x[2] - x[4]) / se;

	if (j == 1) {
	    pprintf(prn, _("Test statistic: t(%d) = (%g - %g)/%g = %g\n"),
		    n1 + n2 - 2, x[0], x[2], se, ts);
	    if (ts > 0) {
		pv = student_pvalue_2(ts, n1 + n2 - 2);
	    } else {
		pv = student_pvalue_2(-ts, n1 + n2 - 2);
	    }
	    print_pv(prn, pv, 0.5 * pv);
	    if (grf) {
		htest_graph(1, ts, n1 + n2 - 2, 0);
	    }
	} else {
	    if (x[4] > 0.0) {
		pprintf(prn, _("Test statistic: z = (%g - %g - %g)/%g = %g\n"),
			x[0], x[2], x[4], se, ts);
	    } else if (x[4] < 0.0) {
		pprintf(prn, _("Test statistic: z = [(%g - %g) - (%g)]/%g = %g\n"),
			x[0], x[2], x[4], se, ts);
	    } else {
		pprintf(prn, _("Test statistic: z = (%g - %g)/%g = %g\n"),
			x[0], x[2], se, ts);
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
	pprintf(prn, _("Null hypothesis: The population variances are "
		       "equal\n"));
	pprintf(prn, _("Sample 1:\n n = %d, variance = %g\n"), n1, x[0]);
	pprintf(prn, _("Sample 2:\n n = %d, variance = %g\n"), n2, x[1]);

	if (x[0] > x[1]) {
	    ts = x[0] / x[1];
	    pprintf(prn, _("Test statistic: F(%d, %d) = %g\n"), 
		    n1 - 1, n2 - 1, ts);
	    pv = snedecor_cdf_comp(ts, n1 - 1, n2 - 1);
	} else {
	    ts = x[1] / x[0];
	    pprintf(prn, _("Test statistic: F(%d, %d) = %g\n"), 
		    n2 - 1, n1 - 1, ts);
	    pv = snedecor_cdf_comp(ts, n2 - 1, n1 - 1);
	}

	print_pv(prn, 2.0 * pv, pv);
	if (grf) {
	    htest_graph(3, ts, n1 - 1, n2 - 1);
	}
	break;

    case TWO_PROPNS:
	pprintf(prn, _("Null hypothesis: the population proportions are "
		       "equal\n"));
	pputc(prn, '\n');
	pprintf(prn, _("Sample 1:\n n = %d, proportion = %g\n"), n1, x[0]);
	pputc(prn, '\n');
	pprintf(prn, _("Sample 2:\n n = %d, proportion = %g\n"), n2, x[1]);
	pputc(prn, '\n');

	x[2] = (n1*x[0] + n2*x[1]) / (n1 + n2);
	se = sqrt((x[2] * (1.0-x[2])) * (1.0/n1 + 1.0/n2));
	ts = (x[0] - x[1]) / se;

	pprintf(prn, _("Test statistic: z = (%g - %g) / %g = %g\n"),
		x[0], x[1], se, ts);

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

/* FIXME : should we record a relevant command when a test is done
   using dataset variables? (Two means, or two variances)
*/

static void h_test (GtkWidget *w, test_t *test)
{
    int j, k, n1 = 0, n2 = 0;
    double x[5] = {0};

    j = k = 0;

    switch (test->code) {

    case ONE_MEAN:
	x[j] = getval(test->entry[k++], C_DBL); /* mean */
	if (na(x[j++])) return;
	x[j] = getval(test->entry[k++], C_POS_DBL); /* s.d. */
	if (na(x[1])) return;
	n1 = getval(test->entry[k++], C_INT); /* sample */
	if (n1 < 0) return;
	x[j++] = getval(test->entry[k], C_DBL); /* val */
	if (na(x[j])) return;
	break;

    case ONE_VARIANCE:
	x[j] = getval(test->entry[k++], C_POS_DBL);
	if (na(x[j++])) return;
	n1 = getval(test->entry[k++], C_INT);
	if (n1 < 0) return;
	x[j] = getval(test->entry[k], C_POS_DBL);
	if (na(x[j])) return;
	break;

    case ONE_PROPN:
	x[j] = getval(test->entry[k++], C_POS_DBL); /* propn */
	if (na(x[j++])) return;
	n1 = getval(test->entry[k++], C_INT);
	if (n1 < 0) return;
	x[j] = getval(test->entry[k], C_POS_DBL);
	if (na(x[j])) return;

	if (n1 * x[1] < 5.0 || n1 * (1.0 - x[1]) < 5.0) {
	    infobox(_("The assumption of a normal sampling distribution\n"
		      "is not justified here.  Abandoning the test."));
	    return;
	}
	break;

    case TWO_MEANS:
	x[j] = getval(test->entry[k++], C_DBL); /* mean1 */
	if (na(x[j++])) return;
	x[j] = getval(test->entry[k++], C_POS_DBL); /* sd1 */
	if (na(x[j++])) return;
	n1 = getval(test->entry[k++], C_INT);
	if (n1 < 0) return; 

	x[j] = getval(test->entry[k++], C_DBL); /* mean2 */
	if (na(x[j++])) return;
	x[j] = getval(test->entry[k++], C_POS_DBL); /* sd2 */
	if (na(x[j++])) return;
	n2 = getval(test->entry[k++], C_INT);
	if (n2 < 0) return; 

	x[j] = getval(test->entry[k], C_DBL);
	if (na(x[j])) return;
	break;

    case TWO_VARIANCES:
	x[j] = getval(test->entry[k++], C_POS_DBL);
	if (na(x[j++])) return;
	n1 = getval(test->entry[k++], C_INT);
	if (n1 < 0) return;

	x[1] = getval(test->entry[k++], C_POS_DBL);
	if (na(x[1])) return;
	n2 = getval(test->entry[k], C_INT);
	if (n2 < 0) return;
	break;

    case TWO_PROPNS:
	x[j] = getval(test->entry[k++], C_POS_DBL);
	if (na(x[j++])) return;
	n1 = getval(test->entry[k++], C_INT);
	if (n1 < 0) return;

	x[j] = getval(test->entry[k++], C_POS_DBL);
	if (na(x[j])) return;
	n2 = getval(test->entry[k], C_INT);
	if (n2 < 0) return;
	break;

    default:
	break;
    }

    do_h_test(test, x, n1, n2);
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
    lookup_t *look, **looks = child->calcp;
    int m = 0, n = 0;
    int d;

    d = gtk_notebook_get_current_page(GTK_NOTEBOOK(child->book));
    look = looks[d];

    switch (d) {
    case NORMAL_DIST:
	break;
    case T_DIST:
    case CHISQ_DIST:
	m = getval(look->entry[0], C_INT);
	if (m < 0) return;
	break;
    case F_DIST:
	m = getval(look->entry[0], C_INT);
	n = getval(look->entry[1], C_INT);
	if (m < 0 || n < 0) return;
	break;
    }

    if (child->plot != NULL) {
	revise_distribution_plotspec(child->plot, d, m, n);
    } else {
	htest_graph(d, NADBL, m, n);
    } 
}

static void add_calc_entry (GtkWidget *tbl, gint *rows, 
			    const gchar *label, CalcChild *child,
			    int i)
{
    lookup_t **look = child->calcp;
    int c = child->code;
    GtkWidget *tmp;

    *rows += 1;

    /* label */
    gtk_table_resize(GTK_TABLE(tbl), *rows, 2);
    tmp = gtk_label_new(_(label));
    gtk_misc_set_alignment(GTK_MISC(tmp), 1, 0.5);
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      tmp, 0, 1, *rows - 1, *rows);
    gtk_widget_show(tmp);

    /* entry box */
    tmp = gtk_entry_new();
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      tmp, 1, 2, *rows - 1, *rows);
    gtk_widget_show(tmp);
    look[i]->entry[*rows-2] = tmp;

    g_signal_connect(G_OBJECT(tmp), "activate", 
		     (c == CALC_PVAL)? G_CALLBACK(get_pvalue) : 
		     (c == CALC_RAND)? G_CALLBACK(get_random) :
		     (c == CALC_GRAPH || c == CALC_GRAPH_ADD)? 
		     G_CALLBACK(dist_graph) :
		     G_CALLBACK(get_critical),
		     child);
}

static void make_lookup_tab (CalcChild *child, int d)
{
    GtkWidget *tmp, *box, *tbl;
    gint rows;
    const gchar *titles[] = {
	N_("normal"), 
	N_(" t "), 
	N_("chi-square"), 
	N_(" F "), 
	N_("binomial"),
	N_("poisson"),
	N_(" DW "),
    };
   
    box = gtk_vbox_new(FALSE, 0);
    gtk_container_set_border_width(GTK_CONTAINER(box), 10);
    gtk_widget_show(box);

    tmp = gtk_label_new(_(titles[d]));
    gtk_widget_show(tmp);
    gtk_notebook_append_page(GTK_NOTEBOOK(child->book), box, tmp);   

    rows = 1;
    tbl = gtk_table_new(rows, 2, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE (tbl), 5);
    gtk_table_set_col_spacings(GTK_TABLE (tbl), 5);
    gtk_box_pack_start(GTK_BOX(box), tbl, FALSE, FALSE, 0);
    gtk_widget_show(tbl);

    switch (d) {
    case NORMAL_DIST:
	break;
    case T_DIST:
	add_calc_entry(tbl, &rows, N_("df"), child, d);
	break;
    case CHISQ_DIST:
	add_calc_entry(tbl, &rows, N_("df"), child, d);
	break;
    case F_DIST:
	add_calc_entry(tbl, &rows, N_("dfn"), child, d);
	add_calc_entry(tbl, &rows, N_("dfd"), child, d);
	break;	
    case BINOMIAL_DIST:
	add_calc_entry(tbl, &rows, N_("Prob"), child, d);
	add_calc_entry(tbl, &rows, N_("trials"), child, d);
	break;
    case POISSON_DIST:
	add_calc_entry(tbl, &rows, N_("mean"), child, d);
	break;
    case DW_DIST:
	add_calc_entry(tbl, &rows, N_("n"), child, d);
	break;
    default:
	break;
    } 

    if (child->code == CALC_DIST && d != DW_DIST) {
	add_calc_entry(tbl, &rows, N_("right-tail probability"), child, d);
    }
}

/* make a tab (notebook) page, for p-value lookup */

static void make_pval_tab (CalcChild *child, int d) 
{
    lookup_t **tab = child->calcp;
    GtkWidget *tempwid, *box, *tbl;
    gint rows;
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

    rows = 1;
    tbl = gtk_table_new(rows, 2, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(tbl), 5);
    gtk_table_set_col_spacings(GTK_TABLE(tbl), 5);
    gtk_box_pack_start(GTK_BOX(box), tbl, FALSE, FALSE, 0);
    gtk_widget_show(tbl);
   
    switch (d) {

    case CALC_NORMAL: 
	add_calc_entry(tbl, &rows, N_("mean"), child, d);
	gtk_entry_set_text(GTK_ENTRY(tab[d]->entry[0]), "0");
	add_calc_entry(tbl, &rows, N_("std. deviation"), child, d);
	gtk_entry_set_text(GTK_ENTRY(tab[d]->entry[1]), "1");
	add_calc_entry(tbl, &rows, N_("value"), child, d);
	break;

    case CALC_STUDENT:
	add_calc_entry(tbl, &rows, N_("df"), child, d);
	add_calc_entry(tbl, &rows, N_("value"), child, d);
	break;

    case CALC_CHISQ:
	add_calc_entry(tbl, &rows, N_("df"), child, d);
	add_calc_entry(tbl, &rows, N_("value"), child, d);
	break;

    case CALC_SNEDECOR:
	add_calc_entry(tbl, &rows, N_("dfn"), child, d);
	add_calc_entry(tbl, &rows, N_("dfd"), child, d);
	add_calc_entry(tbl, &rows, N_("value"), child, d);
	break;

    case CALC_GAMMA:
	add_calc_entry(tbl, &rows, N_("shape"), child, d);
	add_calc_entry(tbl, &rows, N_("scale"), child, d);
	add_calc_entry(tbl, &rows, N_("value"), child, d);
	break;

    case CALC_BINOMIAL:
	add_calc_entry(tbl, &rows, N_("Prob"), child, d);
	add_calc_entry(tbl, &rows, N_("trials"), child, d);
	add_calc_entry(tbl, &rows, N_("value"), child, d);
	break;

    case CALC_POISSON:
	add_calc_entry(tbl, &rows, N_("mean"), child, d);
	add_calc_entry(tbl, &rows, N_("value"), child, d);
	break;

    default:
	break;
    } 
}

/* make a tab (notebook) page, for r.v. generation */

static void make_rand_tab (CalcChild *child, int d) 
{
    lookup_t **tab = child->calcp;
    GtkWidget *tempwid, *box, *tbl;
    gint rows;
    const gchar *titles[] = {
	N_("uniform"),
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

    rows = 1;
    tbl = gtk_table_new(rows, 2, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(tbl), 5);
    gtk_table_set_col_spacings(GTK_TABLE(tbl), 5);
    gtk_box_pack_start(GTK_BOX(box), tbl, FALSE, FALSE, 0);
    gtk_widget_show(tbl);

    switch (d) {

    case RAND_UNIFORM:
	add_calc_entry(tbl, &rows, N_("minimum"), child, d);
	gtk_entry_set_text(GTK_ENTRY(tab[d]->entry[0]), "0");
	add_calc_entry(tbl, &rows, N_("maximum"), child, d);
	gtk_entry_set_text(GTK_ENTRY(tab[d]->entry[1]), "1");
	break;

    case RAND_NORMAL: 
	add_calc_entry(tbl, &rows, N_("mean"), child, d);
	gtk_entry_set_text(GTK_ENTRY(tab[d]->entry[0]), "0");
	add_calc_entry(tbl, &rows, N_("std. deviation"), child, d);
	gtk_entry_set_text(GTK_ENTRY(tab[d]->entry[1]), "1");
	break;

    case RAND_STUDENT:
	add_calc_entry(tbl, &rows, N_("df"), child, d);
	break;

    case RAND_CHISQ:
	add_calc_entry(tbl, &rows, N_("df"), child, d);
	break;

    case RAND_SNEDECOR:
	add_calc_entry(tbl, &rows, N_("dfn"), child, d);
	add_calc_entry(tbl, &rows, N_("dfd"), child, d);
	break;

    case RAND_GAMMA:
	add_calc_entry(tbl, &rows, N_("shape"), child, d);
	add_calc_entry(tbl, &rows, N_("scale"), child, d);
	break;

    case RAND_BINOMIAL:
	add_calc_entry(tbl, &rows, N_("Prob"), child, d);
	add_calc_entry(tbl, &rows, N_("trials"), child, d);
	break;

    case RAND_POISSON:
	add_calc_entry(tbl, &rows, N_("mean"), child, d);
	break;

    default:
	break;
    } 

    add_calc_entry(tbl, &rows, N_("name"), child, d);
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

static void add_test_var_selector (GtkWidget *tbl, gint *rows, 
				   test_t *test, int pos,
				   int labelit)
{
    GtkWidget *label, *tmp;
    gchar **pbuf;

    *rows += 1;
    gtk_table_resize(GTK_TABLE(tbl), *rows, 2);
    if (labelit) {
	gchar *tmp = g_strdup_printf(_("Variable %d"), pos + 1);
	label = gtk_label_new(tmp);
	g_free(tmp);
    } else {
	label = gtk_label_new(_("Variable"));
    }
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      label, 0, 1, *rows - 1, *rows);
    gtk_widget_show(label);

    tmp = gtk_combo_new();
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      tmp, 1, 2, *rows - 1, *rows);
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

static void add_test_combo (GtkWidget *tbl, gint *rows, 
			    test_t *test, int pos)
{
    GtkWidget *button, *tmp;
    gchar **pbuf;

    *rows += 1;
    gtk_table_resize(GTK_TABLE(tbl), *rows, 2);
    button = gtk_check_button_new_with_label(_("Use variable from dataset"));
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      button, 0, 1, *rows - 1, *rows);
    gtk_widget_show(button);

    tmp = gtk_combo_new();
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      tmp, 1, 2, *rows - 1, *rows);
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

static void add_test_entry (GtkWidget *tbl, gint *rows, 
			    const gchar *label, test_t *test, 
			    int i)
{
    GtkWidget *tmp;

    *rows += 1;
    gtk_table_resize(GTK_TABLE(tbl), *rows, 2);
    tmp = gtk_label_new(label);
    gtk_misc_set_alignment(GTK_MISC(tmp), 1, 0.5);
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      tmp, 0, 1, *rows - 1, *rows);
    gtk_widget_show(tmp);
    tmp = gtk_entry_new();
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      tmp, 1, 2, *rows - 1, *rows);
    gtk_widget_show(tmp);
    test->entry[i] = tmp;

    g_signal_connect(G_OBJECT(tmp), "activate", 
		     G_CALLBACK(h_test), test);
}

static void add_test_label (GtkWidget *tbl, gint *rows, 
			    const gchar *label)
{
    GtkWidget *tmp;

    *rows += 1;
    gtk_table_resize(GTK_TABLE(tbl), *rows, 2);
    tmp = gtk_label_new(label);
    gtk_misc_set_alignment(GTK_MISC(tmp), 0, 0.5);
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      tmp, 0, 2, *rows - 1, *rows);
    gtk_widget_show(tmp);
}

static void add_test_check (GtkWidget *tbl, gint *rows, 
			    const gchar *label, test_t *test,
			    gboolean val)
{
    GtkWidget *tmp;

    *rows += 1;
    gtk_table_resize(GTK_TABLE(tbl), *rows, 2);
    tmp = gtk_check_button_new_with_label(label);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), val);
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      tmp, 0, 2, *rows - 1, *rows);
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
    gint i, rows;
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

    rows = 1;
    tbl = gtk_table_new(rows, 2, FALSE);
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
	add_test_var_selector(tbl, &rows, test, 0, 1);
	add_test_var_selector(tbl, &rows, test, 1, 1);

	/* option radios */
	rows += 3;
	gtk_table_resize(GTK_TABLE(tbl), rows, 2);

	test->radio[0] = gtk_radio_button_new_with_label(NULL, _("Sign test"));
	gtk_table_attach_defaults(GTK_TABLE(tbl), test->radio[0], 0, 2, 
				  rows - 3, rows - 2);
	gtk_widget_show(test->radio[0]);

	group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(test->radio[0]));
	test->radio[1] = gtk_radio_button_new_with_label(group, _("Wilcoxon rank sum test"));
	gtk_table_attach_defaults(GTK_TABLE(tbl), test->radio[1], 0, 2, 
				  rows - 2, rows - 1);
	gtk_widget_show(test->radio[1]);

	group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(test->radio[1]));
	test->radio[2] = gtk_radio_button_new_with_label(group, _("Wilcoxon signed rank test"));
	gtk_table_attach_defaults(GTK_TABLE(tbl), test->radio[2], 0, 2, 
				  rows - 1, rows);
	gtk_widget_show(test->radio[2]);
	break;

    case NP_RUNS: 
	add_test_var_selector(tbl, &rows, test, 0, 0);
	break;

    default:
	break;
    } 

    /* check box for extra option */
    rows += 1;
    gtk_table_resize(GTK_TABLE(tbl), rows, 2);
    if (idx == NP_DIFF) {
	tmp = gtk_check_button_new_with_label(_("Show details"));
    } else {
	tmp = gtk_check_button_new_with_label(_("Use first difference"));
    }
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 2, 
			      rows - 1, rows);
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

/* make tab (notebook page) for hypothesis test */

static void make_test_tab (CalcChild *child, int idx) 
{
    test_t **tests = child->calcp;
    test_t *test = tests[idx];
    GtkWidget *tempwid, *box, *tbl;
    int nv = 0;
    gint i, rows;
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

    rows = 1;
    tbl = gtk_table_new(rows, 2, FALSE);
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
	add_test_combo(tbl, &rows, test, 0);
    }
   
    switch (idx) {

    case ONE_MEAN: 
	add_test_entry(tbl, &rows, _("sample mean"), test, 0);
	add_test_entry(tbl, &rows, _("std. deviation"), test, 1);
	add_test_entry(tbl, &rows, _("sample size"), test, 2);
	add_test_entry(tbl, &rows, _("H0: mean ="), test, 3);
	add_test_check(tbl, &rows, _("Assume standard deviation is "
		       "population value"), test, FALSE);
	break;

    case ONE_VARIANCE: 
	add_test_entry(tbl, &rows, _("sample variance"), test, 0);
	add_test_entry(tbl, &rows, _("sample size"), test, 1);
	add_test_entry(tbl, &rows, _("H0: variance ="), test, 2);
	break;

    case ONE_PROPN: /* proportion */
	add_test_entry(tbl, &rows, _("sample proportion"), test, 0);
	add_test_entry(tbl, &rows, _("sample size"), test, 1);
	add_test_entry(tbl, &rows, _("H0: proportion ="), test, 2);
	break;

    case TWO_MEANS:
	add_test_entry(tbl, &rows, _("mean of sample 1"), test, 0);
	add_test_entry(tbl, &rows, _("std. deviation, sample 1"), test, 1);
	add_test_entry(tbl, &rows, _("size of sample 1"), test, 2);
	if (nv > 0) {
	    add_test_combo(tbl, &rows, test, 3);
	}
	add_test_entry(tbl, &rows, _("mean of sample 2"), test, 3);
	add_test_entry(tbl, &rows, _("std. deviation, sample 2"), test, 4);
	add_test_entry(tbl, &rows, _("size of sample 2"), test, 5);
	add_test_entry(tbl, &rows, _("H0: Difference of means ="), test, 6);
	gtk_entry_set_text(GTK_ENTRY(test->entry[6]), "0");
	add_test_check(tbl, &rows, _("Assume common population standard "
		       "deviation"), test, TRUE);
	break;

    case TWO_VARIANCES:
	add_test_entry(tbl, &rows, _("variance of sample 1"), test, 0);
	add_test_entry(tbl, &rows, _("size of sample 1"), test, 1);
	if (nv > 0) {
	    add_test_combo(tbl, &rows, test, 2);
	}	
	add_test_entry(tbl, &rows, _("variance of sample 2"), test, 2);
	add_test_entry(tbl, &rows, _("size of sample 2"), test, 3);
	add_test_label(tbl, &rows, _("H0: Ratio of variances = 1"));
	break;

    case TWO_PROPNS:
	add_test_entry(tbl, &rows, _("proportion, sample 1"), test, 0);
	add_test_entry(tbl, &rows, _("size of sample 1"), test, 1);
	if (nv > 0) {
	    add_test_combo(tbl, &rows, test, 2);
	}	
	add_test_entry(tbl, &rows, _("proportion, sample 2"), test, 2);
	add_test_entry(tbl, &rows, _("size of sample 2"), test, 3);
	add_test_label(tbl, &rows, _("H0: Difference of proportions = 0"));
	break;

    default:
	break;
    } 

    /* add check box for showing graph of sampling dist. */
    rows += 1;
    gtk_table_resize(GTK_TABLE(tbl), rows, 2);
    tempwid = gtk_check_button_new_with_label(_("Show graph of sampling "
						"distribution"));
    gtk_table_attach_defaults(GTK_TABLE(tbl), tempwid, 0, 2, 
			      rows - 1, rows);
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
	    (c == CALC_PVAL)? NPVAL : 
	    (c == CALC_PVAL)? NRAND : 
	    NGRAPHS;

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
	    (c == CALC_RAND)? NRAND :
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

static int calc_help_code (int c)
{
    int hc = 0;

    if (c == CALC_TEST) {
	hc = HTEST;
    } else if (c == CALC_NPTEST) {
	hc = HTESTNP;
    } else if (c == CALC_RAND) {
	hc = GENR_RANDOM;
    }

    return hc;
}

void stats_calculator (gpointer data, guint code, GtkWidget *widget) 
{
    GtkWidget *tmp = NULL;
    static GtkWidget *winptr[6];
    GtkWidget *oldwin;
    CalcChild *child;
    const char *calc_titles[] = {
	N_("gretl: p-value finder"),
	N_("gretl: critical values"),
	N_("gretl: test calculator"),
	N_("gretl: nonparametric tests"),
	N_("gretl: distribution graphs"),
	N_("gretl: add distribution graph"),
	N_("gretl: add random variable")
    };
    int i, hcode, nv = 0;

    g_return_if_fail(code == CALC_PVAL || 
		     code == CALC_DIST || 
		     code == CALC_TEST ||
		     code == CALC_NPTEST ||
		     code == CALC_GRAPH ||
		     code == CALC_GRAPH_ADD ||
		     code == CALC_RAND);

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
	    make_pval_tab(child, i);
	}	
    } else if (code == CALC_DIST) {
	for (i=0; i<NLOOKUPS; i++) {
	    make_lookup_tab(child, i);
	}
    } else if (code == CALC_RAND) {
	for (i=0; i<NRAND; i++) {
	    make_rand_tab(child, i);
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
		     (code == CALC_RAND)? G_CALLBACK(get_random) :
		     (code == CALC_DIST)? G_CALLBACK(get_critical) :
		     (code == CALC_NPTEST)? G_CALLBACK(np_test_global) :
		     (code == CALC_GRAPH || code == CALC_GRAPH_ADD)? 
		     G_CALLBACK(dist_graph) :
		     G_CALLBACK(h_test_global),
		     child);
    gtk_widget_show(tmp);

    /* Help button? */
    hcode = calc_help_code(code);
    if (hcode) { 
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
