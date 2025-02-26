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

#include "gretl.h"
#include "calculator.h"
#include "dlgutils.h"
#include "gpt_control.h"
#include "lib_private.h"
#include "cmdstack.h"
#include "winstack.h"
#include "gui_utils.h"

#define NTESTS 6
#define NPTESTS 3
#define NPVAL 8
#define NDISTS 9
#define NGRAPHS 8
#define NRAND 9
#define NTESTENTRY 7
#define NDISTENTRY 4

typedef struct CalcChild_ CalcChild;
typedef struct test_t_ test_t;
typedef struct dist_t_ dist_t;

struct CalcChild_ {
    int code;
    int n_pages;
    GtkWidget *dlg;
    GtkWidget *vbox;
    GtkWidget *bbox;
    GtkWidget *book;
    gpointer calcp;
    gpointer winp;
    png_plot *plot;
    GCallback callback;
};

struct test_t_ {
    int category;
    int code;
    GtkWidget *entry[NTESTENTRY];
    GtkWidget *combo[2];
    GtkWidget *check;
    GtkWidget *radio[3];
    GtkWidget *extra;
};

struct dist_t_ {
    int flags;
    GtkWidget *entry[NDISTENTRY];
    GtkWidget *check;
};

enum {
    UNIFORM_DIST,
    NORMAL_DIST,
    T_DIST,
    CHISQ_DIST,
    F_DIST,
    GAMMA_DIST,
    BINOMIAL_DIST,
    POISSON_DIST,
    WEIBULL_DIST,
    DW_DIST,
    MAX_DIST
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
    NP_RUNS,
    NP_CORR
};

static void plot_cdf (GtkWidget *parent);

/* functions relating to distribution graphics */

int use_alt_form (int d, double *parms)
{
    if (d == CHISQ_DIST && parms[0] > 69) {
	return 1;
    } else if (d == BINOMIAL_DIST) {
	double p = parms[0];
	int n = parms[1];

	return (n*p > 8 && n*(1-p) > 8);
    } else {
	return 0;
    }
}

static int get_dist_and_params (const char *s, int *d, double *x)
{
    int ret = 1;

    if (!strncmp(s, "# standard", 10)) {
	*d = NORMAL_DIST;
	x[0] = 0;
	x[1] = 1;
    } else if (sscanf(s, "# N(%lf,%lf)", &x[0], &x[1]) == 2) {
	*d = NORMAL_DIST;
    } else if (sscanf(s, "# t(%lf)", &x[0]) == 1) {
	*d = T_DIST;
    } else if (sscanf(s, "# chi-square(%lf)", &x[0]) == 1) {
	*d = CHISQ_DIST;
    } else if (sscanf(s, "# F(%lf,%lf)", &x[0], &x[1]) == 2) {
	*d = F_DIST;
    } else if (sscanf(s, "# Binomial(%lf,%lf)", &x[1], &x[0]) == 2) {
	*d = BINOMIAL_DIST;
    } else if (sscanf(s, "# Poisson(%lf)", &x[0]) == 1) {
	*d = POISSON_DIST;
    } else if (sscanf(s, "# Weibull(%lf,%lf)", &x[0], &x[1]) == 2) {
	*d = WEIBULL_DIST;
    } else if (sscanf(s, "# Gamma(%lf,%lf)", &x[0], &x[1]) == 2) {
	*d = GAMMA_DIST;
    } else {
	ret = 0;
    }

    return ret;
}

static int current_graph_dist (png_plot *plot)
{
    GPT_SPEC *spec = plot_get_spec(plot);
    const char *s;
    int i, d = -1;

    for (i=0; i<spec->n_literal && d<0; i++) {
	s = spec->literal[i];
	if (!strncmp(s, "# standard", 10)) {
	    d = NORMAL_DIST;
	} else if (!strncmp(s, "# N(", 4)) {
	    d = NORMAL_DIST;
	} else if (!strncmp(s, "# t(", 4)) {
	    d = T_DIST;
	} else if (!strncmp(s, "# chi-square(", 13)) {
	    d = CHISQ_DIST;
	} else if (!strncmp(s, "# F(", 4)) {
	    d = F_DIST;
	} else if (!strncmp(s, "# Binomial(", 11)) {
	    d = BINOMIAL_DIST;
	} else if (!strncmp(s, "# Poisson(", 10)) {
	    d = POISSON_DIST;
	} else if (!strncmp(s, "# Weibull(", 10)) {
	    d = WEIBULL_DIST;
	} else if (!strncmp(s, "# Gamma(", 8)) {
	    d = GAMMA_DIST;
	}
    }

    return d;
}

static gchar *
htest_graph_title (int dist, double x, double *parms)
{
    gchar *s = NULL;

    if (dist == NORMAL_DIST) {
	s = g_strdup(_("Gaussian sampling distribution"));
    } else if (dist == T_DIST) {
	s = g_strdup_printf(_("t(%d) sampling distribution"), (int) parms[0]);
    } else if (dist == CHISQ_DIST) {
	s = g_strdup_printf(_("Chi-square(%d) sampling distribution"),
			    (int) parms[0]);
    } else if (dist == F_DIST) {
	s = g_strdup_printf(_("F(%d, %d) sampling distribution"),
			    (int) parms[0], (int) parms[1]);
    }

    return s;
}

static gchar *
dist_graph_title (int dist, double *parms)
{
    gchar *s = NULL;

    if (dist == NORMAL_DIST) {
	s = g_strdup_printf(_("N(%g, %g)"), parms[0], parms[1] * parms[1]);
    } else if (dist == T_DIST) {
	s = g_strdup_printf(_("t(%d)"), (int) parms[0]);
    } else if (dist == CHISQ_DIST) {
	s = g_strdup_printf(_("Chi-square(%d)"), (int) parms[0]);
    } else if (dist == F_DIST) {
	s = g_strdup_printf(_("F(%d, %d)"), (int) parms[0], (int) parms[1]);
    } else if (dist == BINOMIAL_DIST) {
	s = g_strdup_printf(_("Binomial(%d, %g)"), (int) parms[1], parms[0]);
    } else if (dist == POISSON_DIST) {
	s = g_strdup_printf(_("Poisson(%g)"), parms[0]);
    } else if (dist == WEIBULL_DIST) {
	s = g_strdup_printf(_("Weibull(%g, %g)"), parms[0], parms[1]);
    } else if (dist == GAMMA_DIST) {
	s = g_strdup_printf(_("Gamma(%g, %g)"), parms[0], parms[1]);
    }

    return s;
}

static gchar *dist_comment_line (int dist, double *parms)
{
    gchar *s = NULL;

    if (dist == NORMAL_DIST) {
	s = g_strdup_printf("# N(%g,%g)", parms[0], parms[1]);
    } else if (dist == T_DIST) {
	s = g_strdup_printf("# t(%d)", (int) parms[0]);
    } else if (dist == CHISQ_DIST) {
	s = g_strdup_printf("# chi-square(%d)", (int) parms[0]);
    } else if (dist == F_DIST) {
	s = g_strdup_printf("# F(%d,%d)", (int) parms[0], (int) parms[1]);
    } else if (dist == BINOMIAL_DIST) {
	s = g_strdup_printf("# Binomial(%d,%g)", (int) parms[1], parms[0]);
    } else if (dist == POISSON_DIST) {
	s = g_strdup_printf("# Poisson(%g)", parms[0]);
    } else if (dist == WEIBULL_DIST) {
	s = g_strdup_printf("# Weibull(%g,%g)", parms[0], parms[1]);
    } else if (dist == GAMMA_DIST) {
	s = g_strdup_printf("# Gamma(%g,%g)", parms[0], parms[1]);
    }

    return s;
}

enum {
    ID_D,
    ID_N,
    ID_P,
    ID_L,
    ID_M,
    ID_S,
    ID_SHP,
    ID_SCL,
    ID_MAX
};

#define F_STDN   "normal(x)=1/(sqrt(2*pi))*exp(-(x)**2/2)"
#define F_NORM   "normal(x,mu,s)=1/(s*sqrt(2*pi))*exp(-(x-mu)**2/(2*s*s))"
#define F_BINV   "Binv(p,q)=exp(lgamma(p+q)-lgamma(p)-lgamma(q))"
#define F_CHI    "chi(x,m)=x<0?0.0/0.0:x**(0.5*m-1.0)*exp(-0.5*x)/gamma(0.5*m)/2**(0.5*m)"
#define F_ALTCHI "bigchi(x,m)=x<0?0.0/0.0:exp((0.5*m-1.0)*log(x)-0.5*x-lgamma(0.5*m)-m*0.5*log(2.0))"
#define F_F      "f(x,m,n)=x<0?0.0/0.0:Binv(0.5*m,0.5*n)*(m/n)**(0.5*m)*" \
                 "x**(0.5*m-1.0)/(1.0+m/n*x)**(0.5*(m+n))"
#define F_STUD   "stud(x,m)=Binv(0.5*m,0.5)/sqrt(m)*(1.0+(x*x)/m)**(-0.5*(m+1.0))"
#define F_POIS   "poisson(z,k)=k<0?0.0/0.0:exp(-z)*(z**k)/(int(k))!"
#define F_WEIB   "weibull(x,shp,scl)=x<0?0.0/0.0:(shp/scl)*(x/scl)**(shp-1.0)*exp(-(x/scl)**shp)"
#define F_GAMMA  "gdens(x,shp,scl)=x<0?0.0/0.0:1.0/(gamma(shp)*(scl**shp))*(x**(shp-1.0))*exp(-(x/scl))"
#define F_BINOM  "binom(k,n,p)=k<0||k>n?0.0:exp(lgamma(n+1)-lgamma(n-k+1)-lgamma(k+1)" \
                 "+k*log(p)+(n-k)*log(1.0-p))"
#define F_ALTBIN "bigbin(x,mu,s)=x<0?0.0/0.0:1/(s*sqrt(2*pi))*exp(-(x-mu)**2/(2*s*s))"

static void
dist_xmin_xmax (int d, double *parm, double *xmin, double *xmax, int alt)
{
    double arg;
    int dcode = D_NONE;

    if (d == NORMAL_DIST) {
	*xmin = parm[0] - 4.5 * parm[1];
	*xmax = parm[0] + 4.5 * parm[1];
    } else if (d == T_DIST) {
	*xmin = -4.5;
	*xmax = 4.5;
    } else if (d == CHISQ_DIST) {
	dcode = D_CHISQ;
	*xmin = 0;
	arg = 0.005;
    } else if (d == F_DIST) {
	dcode = D_SNEDECOR;
	*xmin = 0;
	if (parm[0] + parm[1] < 16) {
	    arg = 0.009;
	} else {
	    arg = 0.005;
	}
    } else if (d == BINOMIAL_DIST) {
	if (alt) {
	    int n = parm[1];
	    double p = parm[0];
	    double m = n * p;
	    double s = sqrt(m * (1 - p));

	    *xmin = m - 4 * s;
	    *xmax = m + 4 * s;
	    if (*xmin < 0) {
		*xmin = 0;
	    }
	} else {
	    dcode = D_BINOMIAL;
	    *xmin = 0;
	    arg = 0.001;
	}
    } else if (d == POISSON_DIST) {
	dcode = D_POISSON;
	*xmin = 0;
	arg = 0.0015;
    } else if (d == WEIBULL_DIST) {
	dcode = D_WEIBULL;
	*xmin = 0;
	arg = 0.0004;
    } else if (d == GAMMA_DIST) {
	dcode = D_GAMMA;
	*xmin = 0;
	arg = 0.001;
    }

    if (dcode != D_NONE) {
	*xmax = gretl_get_critval(dcode, parm, arg);
    }
}

static double dist_xmax (int d, double *parm)
{
    double arg = NADBL;
    int dcode = D_NONE;

    switch (d) {
    case CHISQ_DIST:
	dcode = D_CHISQ;
	arg = 0.005;
	break;
    case F_DIST:
	dcode = D_SNEDECOR;
	if (parm[0] + parm[1] < 16) {
	    arg = 0.009;
	} else {
	    arg = 0.005;
	}
	break;
    case BINOMIAL_DIST:
	dcode = D_BINOMIAL;
	arg = 0.001;
	break;
    case POISSON_DIST:
	dcode = D_POISSON;
	arg = 0.0015;
	break;
    case WEIBULL_DIST:
	dcode = D_WEIBULL;
	arg = 0.0004;
	break;
    case GAMMA_DIST:
	dcode = D_GAMMA;
	arg = 0.0002;
	break;
    }

    return gretl_get_critval(dcode, parm, arg);
}

static int n_literal_lines (int d, int ptype)
{
    int n = 0;

    switch (d) {
    case NORMAL_DIST:
	n = (ptype == PLOT_H_TEST)? 2 : 4;
	break;
    case T_DIST:
    case WEIBULL_DIST:
    case GAMMA_DIST:
	n = 4;
	break;
    case CHISQ_DIST:
    case POISSON_DIST:
	n = 3;
	break;
    case F_DIST:
	n = 5;
	break;
    case BINOMIAL_DIST:
	n = 4;
	break;
    }

    return n;
}

static double normal_pdf_height (double s)
{
    return 1 / (s * sqrt(M_2PI));
}

static int tic_step (int t)
{
    return 1 + t / 20;
}

static void
range_from_test_stat (int d, double x, double *parms, double *spike,
		      FILE *fp)
{
    double xx, x1;

    if (d == NORMAL_DIST || d == T_DIST) {
	xx = fabs(x);
	x1 = ((xx > 3.5)? xx + .5 : 3.5);
	*spike = .25;
	fprintf(fp, "set xrange [%.3f:%.3f]\n", -x1, x1);
	fprintf(fp, "set yrange [0:.50]\n");
	fprintf(fp, "set xlabel \"%s\"\n", _("Standard errors"));
    } else {
	x1 = dist_xmax(d, parms);
	if (x > x1) {
	    x1 = 1.1 * x;
	}
	*spike = 1.0 / x1;
	fprintf(fp, "set xrange [0:%.3f]\n", x1);
    }
}

static void
range_from_dist (int d, double *parms, int alt, FILE *fp)
{
    double x, tmin = 0, tmax = 0;

    dist_xmin_xmax(d, parms, &tmin, &tmax, alt);

    fprintf(stderr, "range_from_dist: got %g:%g\n", tmin, tmax);

    fprintf(fp, "set trange [%g:%g]\n", tmin, tmax);

    if (d == NORMAL_DIST) {
	x = normal_pdf_height(parms[1]);
	fprintf(fp, "set yrange [0:%g]\n", x * 1.1);
    } else if (d == T_DIST) {
	fputs("set yrange [0:.50]\n", fp);
    } else if (d == BINOMIAL_DIST || d == POISSON_DIST) {
	fprintf(fp, "set xtics %d\n", tic_step(tmax - tmin));
    }
}

static void
adjust_range_from_dist (int d, double *parms, int alt, GPT_SPEC *spec)
{
    double *t_range = spec->range[GP_T_RANGE];
    double tmin = 0, tmax = 0;

    dist_xmin_xmax(d, parms, &tmin, &tmax, alt);

    if (t_range[0] == 0 && tmin < 0) {
	; /* don't adjust? */
    } else if (tmin < t_range[0]) {
	t_range[0] = tmin;
    }

    if (tmax > t_range[1]) {
	t_range[1] = tmax;
    }

    if (d == NORMAL_DIST && !na(spec->range[GP_Y_RANGE][1])) {
	double ymax = normal_pdf_height(parms[1]);

	if (spec->range[GP_Y_RANGE][1] < ymax) {
	    spec->range[GP_Y_RANGE][1] = 1.05 * ymax;
	}
    }

    if (d != NORMAL_DIST && d != T_DIST) {
	spec->range[GP_Y_RANGE][0] = spec->range[GP_Y_RANGE][1] = NADBL;
    }

    if (d == BINOMIAL_DIST && alt) {
	sprintf(spec->xtics, "%d", tic_step(t_range[1] - t_range[0]));
    }

    spec->samples = 200;
}

static gchar *make_plot_line (int d, int alt, const int *ids)
{
    gchar *s = NULL;
    int k, j;

    switch (d) {
    case NORMAL_DIST:
	k = ids[ID_M] + 1;
	j = ids[ID_S] + 1;
	s = g_strdup_printf("t,normal(t,mu%d,s%d)", k, j);
	break;
    case T_DIST:
	k = ids[ID_D] + 1;
	s = g_strdup_printf("t,stud(t,df%d)", k);
	break;
    case CHISQ_DIST:
	k = ids[ID_D] + 1;
	s = g_strdup_printf("t,%s(t,df%d)", (alt)? "bigchi" : "chi", k);
	break;
    case F_DIST:
	k = ids[ID_D] + 1;
	s = g_strdup_printf("t,f(t,df%d,df%d)", k, k + 1);
	break;
    case BINOMIAL_DIST:
	k = ids[ID_N] + 1;
	j = ids[ID_P] + 1;
	if (alt) {
	    s = g_strdup_printf("int(t),bigbin(int(t)+.5,n%d*p%d,sqrt(n%d*p%d*(1-p%d)))",
		    k, j, k, j, j);
	} else {
	    s = g_strdup_printf("int(t),binom(int(t),n%d,p%d)", k, j);
	}
	break;
    case POISSON_DIST:
	k = ids[ID_L] + 1;
	s = g_strdup_printf("int(t),poisson(lambda%d,int(t))", k);
	break;
    case WEIBULL_DIST:
	k = ids[ID_SHP] + 1;
	j = ids[ID_SCL] + 1;
	s = g_strdup_printf("t,weibull(t,shp%d,scl%d)", k, j);
	break;
    case GAMMA_DIST:
	k = ids[ID_SHP] + 1;
	j = ids[ID_SCL] + 1;
	s = g_strdup_printf("t,gdens(t,shp%d,scl%d)", k, j);
	break;
    }

    return s;
}

static void htest_graph (int d, double x, double *parms)
{
    double spike = 0.0;
    gchar *title = NULL;
    FILE *fp;
    int alt = 0, err = 0;

    fp = open_plot_input_file(PLOT_H_TEST, 0, &err);
    if (err) {
	return;
    }

    alt = use_alt_form(d, parms);

    print_keypos_string(GP_KEY_RIGHT_TOP, fp);

    gretl_push_c_numeric_locale();

    range_from_test_stat(d, x, parms, &spike, fp);

    /* header */
    fprintf(fp, "# literal lines = %d\n", n_literal_lines(d, PLOT_H_TEST));
    title = dist_comment_line(d, parms);
    fprintf(fp, "%s\n", title);
    g_free(title);

    /* required variables and formulae */
    switch (d) {
    case NORMAL_DIST:
	fprintf(fp, "%s\n", F_STDN);
	break;
    case T_DIST:
	fprintf(fp, "df=%.1f\n", parms[0]);
	fprintf(fp, "%s\n", F_BINV);
	fprintf(fp, "%s\n", F_STUD);
	break;
    case CHISQ_DIST:
	fprintf(fp, "df=%.1f\n", parms[0]);
	fprintf(fp, "%s\n", (alt)? F_ALTCHI : F_CHI);
	break;
    case F_DIST:
	fprintf(fp, "dfn=%.1f\n", parms[0]);
	fprintf(fp, "dfd=%.1f\n", parms[1]);
	fprintf(fp, "%s\n", F_BINV);
	fprintf(fp, "%s\n", F_F);
	break;
    }

    fputs("plot \\\n", fp);

    switch (d) {
    case NORMAL_DIST:
	fputs("normal(x)", fp);
	break;
    case T_DIST:
	fputs("stud(x,df)", fp);
	break;
    case CHISQ_DIST:
	fprintf(fp, "%s(x,df)", (alt)? "bigchi" : "chi");
	break;
    case F_DIST:
	fputs("f(x,dfn,dfd)", fp);
	break;
    }

    title = htest_graph_title(d, x, parms);
    fprintf(fp, " title \"%s\" w lines , \\\n", title);
    fprintf(fp, "'-' using 1:2 title \"%s\" w impulses\n",
	    _("test statistic"));
    fprintf(fp, "%g %g\n", x, spike);
    fputs("e\n", fp);
    g_free(title);

    gretl_pop_c_numeric_locale();

    err = finalize_plot_input_file(fp);
    gui_graph_handler(err);
}

static void dist_graph (int d, double *parms)
{
    int ids[ID_MAX] = {0};
    gchar *pline = NULL;
    gchar *comment = NULL;
    gchar *title = NULL;
    FILE *fp;
    int alt = 0, err = 0;

    fp = open_plot_input_file(PLOT_PROB_DIST, 0, &err);
    if (err) {
	return;
    }

    alt = use_alt_form(d, parms);

    print_keypos_string(GP_KEY_RIGHT_TOP, fp);
    fputs("set parametric\n", fp);

    title = dist_graph_title(d, parms);

    gretl_push_c_numeric_locale();

    range_from_dist(d, parms, alt, fp);

    /* header */
    fprintf(fp, "# literal lines = %d\n", n_literal_lines(d, PLOT_PROB_DIST));
    comment = dist_comment_line(d, parms);
    fprintf(fp, "%s\n", comment);
    g_free(comment);

    /* required variables and formulae */
    switch (d) {
    case NORMAL_DIST:
	fprintf(fp, "mu1=%g\n", parms[0]);
	fprintf(fp, "s1=%g\n", parms[1]);
	fprintf(fp, "%s\n", F_NORM);
	break;
    case T_DIST:
	fprintf(fp, "df1=%.1f\n", parms[0]);
	fprintf(fp, "%s\n", F_BINV);
	fprintf(fp, "%s\n", F_STUD);
	break;
    case CHISQ_DIST:
	fprintf(fp, "df1=%.1f\n", parms[0]);
	fprintf(fp, "%s\n", (alt)? F_ALTCHI : F_CHI);
	break;
    case F_DIST:
	fprintf(fp, "df1=%.1f\n", parms[0]);
	fprintf(fp, "df2=%.1f\n", parms[1]);
	fprintf(fp, "%s\n", F_BINV);
	fprintf(fp, "%s\n", F_F);
	break;
    case BINOMIAL_DIST:
	fprintf(fp, "n1=%d\n", (int) parms[1]);
	fprintf(fp, "p1=%g\n", parms[0]);
	fprintf(fp, "%s\n", (alt)? F_ALTBIN : F_BINOM);
	break;
    case POISSON_DIST:
	fprintf(fp, "lambda1=%g\n", parms[0]);
	fprintf(fp, "%s\n", F_POIS);
	break;
    case WEIBULL_DIST:
	fprintf(fp, "shp1=%f\n", parms[0]);
	fprintf(fp, "scl1=%f\n", parms[1]);
	fprintf(fp, "%s\n", F_WEIB);
	break;
    case GAMMA_DIST:
	fprintf(fp, "shp1=%f\n", parms[0]);
	fprintf(fp, "scl1=%f\n", parms[1]);
	fprintf(fp, "%s\n", F_GAMMA);
	break;
    }

    fprintf(fp, "plot \\\n");

    pline = make_plot_line(d, alt, ids);

    if (d == BINOMIAL_DIST || d == POISSON_DIST) {
	fprintf(fp, "%s title \"%s\" w linespoints\n", pline, title);
    } else {
	fprintf(fp, "%s title \"%s\" w lines\n", pline, title);
    }

    g_free(pline);
    g_free(title);

    gretl_pop_c_numeric_locale();

    err = finalize_plot_input_file(fp);
    gui_graph_handler(err);
}

#define ALT_CHI MAX_DIST
#define ALT_BIN (ALT_CHI + 1)

static void revise_distribution_plot (png_plot *plot, int d, double *parms)
{
    GPT_SPEC *spec = plot_get_spec(plot);
    char *title = NULL;
    const char *f1 = NULL, *f2 = NULL;
    char v1[32] = {0};
    char v2[32] = {0};
    char got[MAX_DIST+2] = {0};
    int ids[ID_MAX] = {0};
    double x[2] = {0};
    int i, k, alt = 0;
    int err = 0;

    gretl_push_c_numeric_locale();

    /* check what we already have in the plot spec */

    for (i=0; i<spec->n_literal; i++) {
	const char *s = spec->literal[i];
	int id, prevd;

	if (*s == '#') {
	    if (get_dist_and_params(s, &prevd, x)) {
		if (prevd == d && x[0] == parms[0] && x[1] == parms[1]) {
		    /* no-op: line is already present */
		    return;
		}
		if (prevd == CHISQ_DIST && use_alt_form(prevd, x)) {
		    got[ALT_CHI] = 1;
		} else if (prevd == BINOMIAL_DIST && use_alt_form(prevd, x)) {
		    got[ALT_BIN] = 1;
		} else {
		    got[prevd] = 1;
		}
	    }
	} else if (sscanf(s, "df%d=", &id) == 1) {
	    if (id > ids[ID_D]) {
		ids[ID_D] = id;
	    }
	} else if (sscanf(s, "n%d=", &id) == 1) {
	    if (id > ids[ID_N]) {
		ids[ID_N] = id;
	    }
	} else if (sscanf(s, "p%d=", &id) == 1) {
	    if (id > ids[ID_P]) {
		ids[ID_P] = id;
	    }
	} else if (sscanf(s, "lambda%d=", &id) == 1) {
	    if (id > ids[ID_L]) {
		ids[ID_L] = id;
	    }
	} else if (sscanf(s, "mu%d=", &id) == 1) {
	    if (id > ids[ID_M]) {
		ids[ID_M] = id;
	    }
	} else if (sscanf(s, "s%d=", &id) == 1) {
	    if (id > ids[ID_S]) {
		ids[ID_S] = id;
	    }
	} else if (sscanf(s, "shp%d=", &id) == 1) {
	    if (id > ids[ID_SHP]) {
		ids[ID_SHP] = id;
	    }
	} else if (sscanf(s, "scl%d=", &id) == 1) {
	    if (id > ids[ID_SCL]) {
		ids[ID_SCL] = id;
	    }
	}
    }

    /* alternate forms for some plots */
    alt = use_alt_form(d, parms);

    /* maybe adjust plot range */
    adjust_range_from_dist(d, parms, alt, spec);

    /* add a comment line for the new plot */
    title = dist_comment_line(d, parms);
    err = strings_array_add(&spec->literal, &spec->n_literal, title);
    g_free(title);

    if (err) {
	gui_errmsg(err);
	goto bailout;
    }

    /* add parameter line(s) if needed */

    switch (d) {
    case NORMAL_DIST:
	k = ids[ID_M] + 1;
	sprintf(v1, "mu%d=%g", k, parms[0]);
	k = ids[ID_S] + 1;
	sprintf(v2, "s%d=%g", k, parms[1]);
	break;
    case T_DIST:
    case CHISQ_DIST:
    case F_DIST:
	k = ids[ID_D] + 1;
	sprintf(v1, "df%d=%.1f", k, parms[0]);
	if (d == F_DIST) {
	    sprintf(v2, "df%d=%.1f", k + 1, parms[1]);
	}
	break;
    case BINOMIAL_DIST:
	k = ids[ID_N] + 1;
	sprintf(v1, "n%d=%d", k, (int) parms[1]);
	k = ids[ID_P] + 1;
	sprintf(v2, "p%d=%g", k, parms[0]);
	break;
    case POISSON_DIST:
	k = ids[ID_L] + 1;
	sprintf(v1, "lambda%d=%g", k, parms[0]);
	break;
    case WEIBULL_DIST:
    case GAMMA_DIST:
	k = ids[ID_SHP] + 1;
	sprintf(v1, "shp%d=%f", k, parms[0]);
	k = ids[ID_SCL] + 1;
	sprintf(v2, "scl%d=%f", k, parms[1]);
	break;
    }

    if (*v1) {
	err = strings_array_add(&spec->literal, &spec->n_literal, v1);
    }
    if (!err && *v2) {
	err = strings_array_add(&spec->literal, &spec->n_literal, v2);
    }

    if (err) {
	gui_errmsg(err);
	return;
    }

    /* add any required formula lines */

    switch (d) {
    case NORMAL_DIST:
	if (!got[NORMAL_DIST]) f1 = F_NORM;
	break;
    case T_DIST:
	if (!got[T_DIST] && !got[F_DIST]) {
	    f1 = F_BINV;
	}
	if (!got[T_DIST]) f2 = F_STUD;
	break;
    case CHISQ_DIST:
	if (alt && !got[ALT_CHI]) {
	    f1 = F_ALTCHI;
	} else if (!alt && !got[CHISQ_DIST]) {
	    f1 = F_CHI;
	}
	break;
    case F_DIST:
	if (!got[F_DIST] && !got[T_DIST]) {
	    f1 = F_BINV;
	}
	if (!got[F_DIST]) f2 = F_F;
	break;
    case BINOMIAL_DIST:
	if (alt && !got[ALT_BIN]) {
	    f1 = F_ALTBIN;
	} else if (!alt && !got[BINOMIAL_DIST]) {
	    f1 = F_BINOM;
	}
	break;
    case POISSON_DIST:
	if (!got[POISSON_DIST]) f1 = F_POIS;
	break;
    case WEIBULL_DIST:
	if (!got[WEIBULL_DIST]) f1 = F_WEIB;
	break;
    case GAMMA_DIST:
	if (!got[GAMMA_DIST]) f1 = F_GAMMA;
	break;
    }

    if (f1 != NULL) {
	err = strings_array_add(&spec->literal, &spec->n_literal, f1);
    }
    if (!err && f2 != NULL) {
	err = strings_array_add(&spec->literal, &spec->n_literal, f2);
    }

    if (!err) {
	/* add new plot line */
	err = plotspec_add_line(spec);
    }

    if (err) {
	gui_errmsg(err);
	goto bailout;
    }

    i = spec->n_lines - 1;

    if (d == BINOMIAL_DIST || d == POISSON_DIST) {
	spec->lines[i].style = GP_STYLE_LINESPOINTS;
    } else {
	spec->lines[i].style = GP_STYLE_LINES;
    }

    g_free(spec->lines[i].title);
    spec->lines[i].title = dist_graph_title(d, parms);

    g_free(spec->lines[i].formula);
    spec->lines[i].formula = make_plot_line(d, alt, ids);

    redisplay_edited_plot(plot);

 bailout:

    gretl_pop_c_numeric_locale();
}

/* end of graphics functions */

static double get_real_const (const char *s, EntryValType t)
{
    if (t == C_DBL || t == C_POS_DBL) {
	if (!strcmp(s, "e")) {
	    return 2.71828182845904523536;
	} else if (!strcmp(s, "pi")) {
	    return M_PI;
	}
    }

    return NADBL;
}

/* entry_get_numeric_value: get a numerical value from a text entry box

   t == C_DBL: parse the entry as a double
   t == C_POS_DBL: parse as a positive double
   t == C_INT: parse as non-negative int
   t == C_POS_INT: parse as positive int
   t == C_FRAC: parse as 0 < p < 1

   (In this context, when an int is wanted, it's always non-negative
   and usually positive: df, sample size, number of trials.)

   If t is C_DBL, C_POS_DBL or C_FRAC and something is wrong with
   the input, we flag this by returning NADBL.  In the case of
   C_INT and C_POS_INT, we flag an error by returning -1.
*/

double entry_get_numeric_value (GtkWidget *w, EntryValType t)
{
    const gchar *text = gtk_entry_get_text(GTK_ENTRY(w));
    char s[32];
    double x = NADBL;
    int sub = 0;
    int k, bad = 0;

    gretl_error_clear();

    if (text == NULL || *text == '\0') {
	warnbox(_("Incomplete entry"));
	return (t == C_INT || t == C_POS_INT)? -1 : NADBL;
    }

    *s = '\0';

    while (isspace(*text)) text++;
    strncat(s, text, 31);
    tailstrip(s);

    /* try for constants: e, pi */
    x = get_real_const(s, t);
    if (!na(x)) {
	return x;
    }

    if (get_local_decpoint() != '.') {
	gretl_push_c_numeric_locale();
	gretl_charsub(s, ',', '.');
	sub = 1;
    }

    /* a formula? */
    if (*s == '=') {
	x = generate_scalar(s+1, NULL, &bad);
	if (bad || na(x)) {
	    bad = 1;
	    goto do_message;
	}
    }

    if (t == C_INT || t == C_POS_INT) {
	if (na(x) && check_atoi(s)) {
	    bad = 1;
	} else {
	    if (na(x)) {
		k = atoi(s);
	    } else {
		k = (int) x;
	    }
	    if ((t == C_INT && k < 0) ||
		(t == C_POS_INT && k <= 0)) {
		bad = 1;
	    } else {
		x = k;
	    }
	}
    } else {
	if (na(x) && check_atof(s)) {
	    bad = 1;
	} else {
	    if (na(x)) {
		x = atof(s);
	    }
	    if (t == C_POS_DBL) {
		if (!na(x) && x <= 0.0) {
		    bad = 1;
		}
	    } else if (t == C_FRAC) {
		if (!na(x) && (x <= 0.0 || x >= 1.0)) {
		    bad = 1;
		}
	    }
	}
    }

 do_message:

    if (bad) {
	const char *msg = gretl_errmsg_get();

	warnbox((*msg != '\0')? msg : _("Invalid entry"));
	gtk_editable_select_region(GTK_EDITABLE(w), 0, -1);
	gtk_widget_grab_focus(w);
	x = (t == C_INT || t == C_POS_INT)? -1 : NADBL;
    }

    if (sub) {
	gretl_pop_c_numeric_locale();
    }

    gretl_error_clear();

    return x;
}

#define getval(w,t) entry_get_numeric_value(w,t)

/* call plugin function to look up the table of critical values for
   the Durbin-Watson statistic
*/

static void dw_lookup_call (dist_t *tab)
{
    int (*dw)(int, int, gretl_matrix **) = NULL;
    gretl_vector *v = NULL;
    PRN *prn;
    int n, k, err = 0;

    n = getval(tab->entry[0], C_POS_INT);
    if (n < 0) {
	return;
    }

    k = getval(tab->entry[1], C_POS_INT);
    if (k < 0) {
	return;
    }

    dw = gui_get_plugin_function("dw_lookup");
    if (dw == NULL) {
	return;
    }

    if (bufopen(&prn)) {
	return;
    }

    err = (*dw)(n, k, &v);

    if (!err) {
	pprintf(prn, "%s, n = %d, k = %d\n\n",
		/* xgettext:no-c-format */
		_("5% critical values for Durbin-Watson statistic"),
		(int) gretl_vector_get(v, 2),
		(int) gretl_vector_get(v, 3));
	pprintf(prn, "  dL = %6.4f\n", gretl_vector_get(v, 0));
	pprintf(prn, "  dU = %6.4f\n", gretl_vector_get(v, 1));
	gretl_vector_free(v);
	view_buffer(prn, 77, 300, _("gretl: critical values"),
		    STAT_TABLE, NULL);
    } else {
	gui_errmsg(err);
    }
}

static int dist_from_page (int code, int page)
{
    if (code == CALC_RAND) {
	return page;
    } else if (code == CALC_PVAL) {
	return page + 1;
    } else {
	switch (page) {
	case 0:
	    return NORMAL_DIST;
	case 1:
	    return T_DIST;
	case 2:
	    return CHISQ_DIST;
	case 3:
	    return F_DIST;
	case 4:
	    return BINOMIAL_DIST;
	case 5:
	    return POISSON_DIST;
	case 6:
	    return WEIBULL_DIST;
	case 7:
	    return GAMMA_DIST;
	case 8:
	    return DW_DIST;
	}
    }

    return 0;
}

/* translate from the subset-encoding of distributions
   used here to the coding used in pvalues.c
*/

static int d_to_pdist (int d)
{
    if (d == UNIFORM_DIST) {
	return D_UNIFORM;
    } else if (d == NORMAL_DIST) {
	return D_NORMAL;
    } else if (d == T_DIST) {
	return D_STUDENT;
    } else if (d == CHISQ_DIST) {
	return D_CHISQ;
    } else if (d == F_DIST) {
	return D_SNEDECOR;
    } else if (d == GAMMA_DIST) {
	return D_GAMMA;
    } else if (d == BINOMIAL_DIST) {
	return D_BINOMIAL;
    } else if (d == POISSON_DIST) {
	return D_POISSON;
    } else if (d == WEIBULL_DIST) {
	return D_WEIBULL;
    } else if (d == DW_DIST) {
	return D_DW;
    } else {
	return D_NONE;
    }
}

static int
get_dist_entry_vector (int code, dist_t *tab, int d, double *x,
		       int *pj)
{
    int j = 0;

    switch (d) {
    case UNIFORM_DIST:
	x[j] = getval(tab->entry[j], C_DBL); /* min */
	if (na(x[j])) return 1; else j++;
	x[j] = getval(tab->entry[j], C_DBL); /* max */
	if (na(x[j])) return 1; else j++;
	break;
    case NORMAL_DIST:
	x[j] = getval(tab->entry[j], C_DBL); /* mean */
	if (na(x[j])) return 1; else j++;
	x[j] = getval(tab->entry[j], C_POS_DBL); /* s.d. */
	if (na(x[j])) return 1; else j++;
	break;
    case T_DIST:
    case CHISQ_DIST:
	x[j] = getval(tab->entry[j], C_POS_INT); /* df */
	if (x[j++] < 0) return 1;
	break;
    case F_DIST:
	x[j] = getval(tab->entry[j], C_POS_INT); /* dfn */
	if (x[j++] < 0) return 1;
	x[j] = getval(tab->entry[j], C_POS_INT); /* dfd */
	if (x[j++] < 0) return 1;
	break;
    case GAMMA_DIST:
    case WEIBULL_DIST:
	x[j] = getval(tab->entry[j], C_POS_DBL); /* shape */
	if (na(x[j])) return 1; else j++;
	x[j] = getval(tab->entry[j], C_POS_DBL); /* scale */
	if (na(x[j])) return 1; else j++;
	break;
    case BINOMIAL_DIST:
	x[j] = getval(tab->entry[j], C_FRAC); /* prob */
	if (na(x[j])) return 1; else j++;
	x[j] = getval(tab->entry[j], C_POS_INT); /* n */
	if (x[j++] < 0) return 1;
	break;
    case POISSON_DIST:
	x[j] = getval(tab->entry[j], C_POS_DBL); /* mean */
	if (na(x[j])) return 1; else j++;
	break;
    }

    if (pj != NULL) {
	*pj = j;
    }

    return 0;
}

static void print_normal_critval (const double *parm, double a,
				  double c, PRN *prn)
{
    double mu = parm[0];
    double s = parm[1];

    if (mu == 0 && s == 1) {
	pprintf(prn, "%s", _("Standard normal distribution"));
    } else {
	pprintf(prn, "N(%g, %g)", mu, s * s);
	c = c * s + mu;
    }

    pputs(prn, "\n ");
    pprintf(prn, _("right-tail probability = %g"), a);
    pputs(prn, "\n ");
    pprintf(prn, _("complementary probability = %g"), 1.0 - a);
    if (a < 0.5) {
	pputs(prn, "\n ");
	pprintf(prn, _("two-tailed probability = %g"), 2.0 * a);
    }
    pputs(prn, "\n\n ");
    pprintf(prn, _("Critical value = %g"), c);
    pputc(prn, '\n');
}

static void get_critical (GtkWidget *w, CalcChild *child)
{
    dist_t **tabs = child->calcp;
    double c = NADBL;
    double a, parm[2];
    int i, d, j = 0;
    int pdist;
    PRN *prn;

    i = gtk_notebook_get_current_page(GTK_NOTEBOOK(child->book));
    d = dist_from_page(child->code, i);

    if (d == DW_DIST) {
	/* special: just a table look-up */
	dw_lookup_call(tabs[i]);
	return;
    }

    if (get_dist_entry_vector(child->code, tabs[i], d, parm, &j)) {
	return;
    }

    /* right-tail probability */
    a = getval(tabs[i]->entry[j], C_FRAC);
    if (na(a)) return;

    pdist = d_to_pdist(d);

    c = gretl_get_critval(pdist, parm, a);
    if (na(c)) {
	errbox(_("Failed to compute critical value"));
	return;
    }

    if (bufopen(&prn)) {
	return;
    }

    if (d == NORMAL_DIST) {
	print_normal_critval(parm, a, c, prn);
    } else {
	print_critval(pdist, parm, a, c, prn);
    }

    view_buffer(prn, 60, 200, _("gretl: critical values"),
		STAT_TABLE, NULL);
}

static void get_pvalue (GtkWidget *w, CalcChild *child)
{
    dist_t **tabs = child->calcp;
    double pv, x = 0, parm[2];
    int pdist = 0;
    int i, d, j = 0;
    PRN *prn;

    i = gtk_notebook_get_current_page(GTK_NOTEBOOK(child->book));
    d = dist_from_page(child->code, i);

    if (get_dist_entry_vector(child->code, tabs[i], d, parm, &j)) {
	return;
    }

    /* x-value for which to get p-value */

    switch (d) {
    case NORMAL_DIST:
    case T_DIST:
	x = getval(tabs[i]->entry[j], C_DBL);
	if (na(x)) return;
	break;
    case CHISQ_DIST:
    case F_DIST:
    case GAMMA_DIST:
    case WEIBULL_DIST:
	x = getval(tabs[i]->entry[j], C_POS_DBL);
	if (na(x)) return;
	break;
    case BINOMIAL_DIST:
    case POISSON_DIST:
	x = getval(tabs[i]->entry[j], C_INT);
	if (x < 0) return;
	break;
    };

    pdist = d_to_pdist(d);

    if (d == NORMAL_DIST) {
	/* transform to z-score */
	x = (x - parm[0]) / parm[1];
    }

    if (bufopen(&prn)) return;

    pv = gretl_get_pvalue(pdist, parm, x);

    if (na(pv)) {
	errbox(_("Failed to compute p-value"));
    } else {
	print_pvalue(pdist, parm, x, pv, prn);
	view_buffer(prn, 78, 200, _("gretl: p-value"), PVAL, NULL);
    }
}

static int calc_finish_genr (void)
{
    PRN *prn;
    int err = 0;

    if (bufopen(&prn)) {
	return 1;
    }

    err = gui_run_genr(get_lib_cmdline(), dataset, OPT_NONE, prn);

    if (err) {
	errbox(gretl_print_get_buffer(prn));
    } else {
	record_command_verbatim();
	populate_varlist();
	mark_dataset_as_modified();
    }

    gretl_print_destroy(prn);

    return err;
}

static void get_random (GtkWidget *w, CalcChild *child)
{
    dist_t **tabs = child->calcp;
    const char *vname;
    double x[2] = {0};
    int i, d, j = 0;
    int err;

    i = gtk_notebook_get_current_page(GTK_NOTEBOOK(child->book));
    d = dist_from_page(child->code, i);

    if (get_dist_entry_vector(child->code, tabs[i], d, x, &j)) {
	return;
    }

    vname = gtk_entry_get_text(GTK_ENTRY(tabs[i]->entry[j]));

    if (vname == NULL || *vname == '\0') {
	warnbox(_("You must give a name for the variable"));
	return;
    } else if (gui_validate_varname(vname,
				    GRETL_TYPE_SERIES,
				    child->dlg)) {
	return;
    }

    if (d == UNIFORM_DIST && x[0] >= x[1]) {
	warnbox(_("Range is non-positive!"));
	return;
    }

    gretl_push_c_numeric_locale();

    switch (d) {
    case UNIFORM_DIST:
	lib_command_sprintf("series %s = randgen(u,%g,%g)", vname,
			    x[0], x[1]);
	break;
    case NORMAL_DIST:
	lib_command_sprintf("series %s = randgen(N,%g,%g)", vname,
			    x[0], x[1]);
	break;
    case T_DIST:
	lib_command_sprintf("series %s = randgen(t,%g)", vname, x[0]);
	break;
    case CHISQ_DIST:
	lib_command_sprintf("series %s = randgen(X,%g)", vname, x[0]);
	break;
    case F_DIST:
	lib_command_sprintf("series %s = randgen(F,%g,%g)", vname,
			    x[0], x[1]);
	break;
    case GAMMA_DIST:
	lib_command_sprintf("series %s = randgen(G,%g,%g)", vname,
			    x[0], x[1]);
	break;
    case WEIBULL_DIST:
	lib_command_sprintf("series %s = randgen(W,%g,%g)", vname,
			    x[0], x[1]);
	break;
    case BINOMIAL_DIST:
	lib_command_sprintf("series %s = randgen(B,%g,%g)", vname,
			    x[0], x[1]);
	break;
    case POISSON_DIST:
	/* FIXME allow variable as param? */
	lib_command_sprintf("series %s = randgen(P,%g)", vname, x[0]);
	break;
    }

    gretl_pop_c_numeric_locale();

    err = calc_finish_genr();

    if (!err) {
	int v = current_series_index(dataset, vname);

	gtk_widget_hide(child->dlg);
	infobox_printf(_("Generated series %s (ID %d)"), vname, v);
	gtk_widget_destroy(child->dlg);
    }
}

static void print_pv (PRN *prn, double p1, double p2)
{
    pprintf(prn, _("Two-tailed p-value = %.4g\n(one-tailed = %.4g)\n"),
	    p1, p2);
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
    v1 = series_index(dataset, var1);

    if (v1 == dataset->v) {
	gui_errmsg(E_UNKVAR);
	return;
    }

    if (test->code == NP_DIFF || test->code == NP_CORR) {
	var2 = gtk_entry_get_text(GTK_ENTRY(test->entry[1]));
	v2 = series_index(dataset, var2);
	if (v2 == dataset->v) {
	    gui_errmsg(E_UNKVAR);
	    return;
	}
    }

    if (bufopen(&prn)) {
	return;
    }

    if (test->code == NP_DIFF) {
	int list[3] = { 2, v1, v2 };

	if (test->extra != NULL && button_is_active(test->extra)) {
	    opt |= OPT_V;
	}

	if (button_is_active(test->radio[0])) {
	    opt |= OPT_G;
	} else if (button_is_active(test->radio[1])) {
	    opt |= OPT_R;
	} else if (button_is_active(test->radio[2])) {
	    opt |= OPT_I;
	}

	err = diff_test(list, dataset, opt, prn);
    } else if (test->code == NP_RUNS) {
	if (test->extra != NULL && button_is_active(test->extra)) {
	    opt |= OPT_D;
	}
	if (test->check != NULL && button_is_active(test->check)) {
	    opt |= OPT_E;
	}
	err = runs_test(v1, dataset, opt, prn);
    } else if (test->code == NP_CORR) {
	int list[3] = { 2, v1, v2 };

	if (button_is_active(test->radio[0])) {
	    /* Kendall */
	    opt |= OPT_K;
	} else if (button_is_active(test->radio[1])) {
	   /* Spearman */
	    opt |= OPT_S;
	}

	if (test->extra != NULL && button_is_active(test->extra)) {
	    opt |= OPT_V;
	}
	if (opt & OPT_K) {
	    err = kendall_tau(list, dataset, opt, prn);
	} else {
	    err = spearman_rho(list, dataset, opt, prn);
	}
    }

    if (err) {
	gui_errmsg(err);
    } else {
	gchar *title = g_strdup_printf("gretl: %s", _("nonparametric test"));

	view_buffer(prn, 78, 380, title, PRINT, NULL);
	g_free(title);
	/* record successful command */
	if (test->code == NP_DIFF) {
	    lib_command_sprintf("difftest %s %s%s",
				dataset->varname[v1],
				dataset->varname[v2],
				print_flags(opt, DIFFTEST));
	    record_command_verbatim();
	} else if (test->code == NP_RUNS) {
	    lib_command_sprintf("runs %s%s", dataset->varname[v1],
				print_flags(opt, RUNS));
	    record_command_verbatim();
	} else if (test->code == NP_CORR) {
	    lib_command_sprintf("corr %s %s%s",
				dataset->varname[v1],
				dataset->varname[v2],
				print_flags(opt, CORR));
	    record_command_verbatim();
	}
    }
}

static void do_two_means_test (double d0,
			       int n1, double xbar1, double s1,
			       int n2, double xbar2, double s2,
			       int common_variance,
			       double *ptest,
			       PRN *prn)
{
    double z, se, test, pv;
    double v1 = s1 * s1;
    double v2 = s2 * s2;
    int df;

    pprintf(prn, _("Null hypothesis: Difference of means = %g\n"), d0);
    pputc(prn, '\n');

    /* sample 1 info */
    pprintf(prn, _("Sample 1:\n n = %d, mean = %g, s.d. = %g\n"),
	    n1, xbar1, s1);
    z = s1 / sqrt((double) n1);
    pprintf(prn, _(" standard error of mean = %g\n"), z);
    z *= tcrit95(n1 - 1);
    pprintf(prn, _(" 95%% confidence interval for mean: %g to %g\n"),
	    xbar1 - z, xbar1 + z);
    pputc(prn, '\n');

    /* sample 2 info */
    pprintf(prn, _("Sample 2:\n n = %d, mean = %g, s.d. = %g\n"),
	    n2, xbar2, s2);
    z = s2 / sqrt((double) n2);
    pprintf(prn, _(" standard error of mean = %g\n"), z);
    z *= tcrit95(n2 - 1);
    pprintf(prn, _(" 95%% confidence interval for mean: %g to %g\n"),
	    xbar2 - z, xbar2 + z);
    pputc(prn, '\n');

    if (common_variance) {
	double v;

	v = ((n1-1) * v1 + (n2-1) * v2) / (n1 + n2 - 2);
	se = sqrt(v / n1 + v / n2);
	df = n1 + n2 - 2;
    } else {
	se = sqrt(v1/n1 + v2/n2);
	df = satterthwaite_df(v1, n1, v2, n2);
    }

    test = (xbar1 - xbar2 - d0) / se;

    if (d0 > 0.0) {
	pprintf(prn, _("Test statistic: t(%d) = (%g - %g - %g)/%g = %g\n"),
		df, xbar1, xbar2, d0, se, test);
    } else if (d0 < 0.0) {
	pprintf(prn, _("Test statistic: t(%d) = [(%g - %g) - (%g)]/%g = %g\n"),
		df, xbar1, xbar2, d0, se, test);
    } else {
	pprintf(prn, _("Test statistic: t(%d) = (%g - %g)/%g = %g\n"),
		df, xbar1, xbar2, se, test);
    }

    if (test > 0) {
	pv = student_pvalue_2(df, test);
    } else {
	pv = student_pvalue_2(df, -test);
    }
    print_pv(prn, pv, pv / 2);

    *ptest = test;
}

static void do_h_test (test_t *test, double *x, int n1, int n2)
{
    double se, ts, pv;
    double gparm[2] = {0};
    int common, grf;
    PRN *prn;

    if (bufopen(&prn)) {
	return;
    }

    grf = button_is_active(test->extra);

    switch (test->code) {

    case ONE_MEAN:
	se = x[1] / sqrt((double) n1);
	ts = (x[0] - x[2]) / se;

	pprintf(prn, _("Null hypothesis: population mean = %g\n"), x[2]);
	pprintf(prn, _("Sample size: n = %d\n"), n1);
	pprintf(prn, _("Sample mean = %g, std. deviation = %g\n"),
		x[0], x[1]);

	if (button_is_active(test->check)) {
	    pprintf(prn, _("Test statistic: z = (%g - %g)/%g = %g\n"),
		    x[0], x[2], se, ts);
	    pv = normal_pvalue_2(ts);
	    print_pv(prn, pv, pv / 2.0);
	    if (grf) {
		gparm[1] = 1;
		htest_graph(NORMAL_DIST, ts, gparm);
	    }
	} else {
	    pprintf(prn, _("Test statistic: t(%d) = (%g - %g)/%g = %g\n"), n1-1,
		    x[0], x[2], se, ts);
	    pv = student_pvalue_2(n1 - 1, ts);
	    print_pv(prn, pv, 0.5 * pv);
	    if (grf) {
		gparm[0] = n1 - 1;
		htest_graph(T_DIST, ts, gparm);
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

	if (x[0] >= x[1]) {
	    pv = chisq_cdf_comp(n1 - 1, ts);
	} else {
	    pv = chisq_cdf(n1 - 1, ts);
	}
	print_pv(prn, 2.0 * pv, pv);
	if (grf) {
	    gparm[0] = n1 - 1;
	    htest_graph(CHISQ_DIST, ts, gparm);
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
	    gparm[1] = 1;
	    htest_graph(NORMAL_DIST, ts, gparm);
	}
	break;

    case TWO_MEANS:
	common = button_is_active(test->check);
	do_two_means_test(x[4], n1, x[0], x[1], n2, x[2], x[3],
			  common, &ts, prn);
	if (grf) {
	    if (common) {
		gparm[0] = n1 + n2 - 2;
		htest_graph(T_DIST, ts, gparm);
	    } else {
		gparm[1] = 1;
		htest_graph(NORMAL_DIST, ts, gparm);
	    }
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
	    pv = snedecor_cdf_comp(n1 - 1, n2 - 1, ts);
	} else {
	    ts = x[1] / x[0];
	    pprintf(prn, _("Test statistic: F(%d, %d) = %g\n"),
		    n2 - 1, n1 - 1, ts);
	    pv = snedecor_cdf_comp(n2 - 1, n1 - 1, ts);
	}

	print_pv(prn, 2.0 * pv, pv);
	if (grf) {
	    gparm[0] = n1 - 1;
	    gparm[1] = n2 - 1;
	    htest_graph(F_DIST, ts, gparm);
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
	se = sqrt((x[2] * (1 - x[2])) * (1.0/n1 + 1.0/n2));
	ts = (x[0] - x[1]) / se;

	pprintf(prn, _("Test statistic: z = (%g - %g) / %g = %g\n"),
		x[0], x[1], se, ts);

	pv = normal_pvalue_2(ts);
	print_pv(prn, pv, pv / 2.0);
	if (grf) {
	    gparm[1] = 1;
	    htest_graph(NORMAL_DIST, ts, gparm);
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
    int k = 0, n1 = 0, n2 = 0;
    double x[5] = {0};

    switch (test->code) {

    case ONE_MEAN:
	x[0] = getval(test->entry[k++], C_DBL); /* mean */
	if (na(x[0])) return;
	x[1] = getval(test->entry[k++], C_POS_DBL); /* s.d. */
	if (na(x[1])) return;
	n1 = getval(test->entry[k++], C_POS_INT); /* n */
	if (n1 < 0) return;
	x[2] = getval(test->entry[k], C_DBL); /* val */
	if (na(x[2])) return;
	break;

    case ONE_VARIANCE:
	x[0] = getval(test->entry[k++], C_POS_DBL);
	if (na(x[0])) return;
	n1 = getval(test->entry[k++], C_POS_INT);
	if (n1 < 0) return;
	x[1] = getval(test->entry[k], C_POS_DBL);
	if (na(x[1])) return;
	break;

    case ONE_PROPN:
	x[0] = getval(test->entry[k++], C_FRAC); /* propn */
	if (na(x[0])) return;
	n1 = getval(test->entry[k++], C_POS_INT);
	if (n1 < 0) return;
	x[1] = getval(test->entry[k], C_FRAC); /* H0 propn */
	if (na(x[1])) return;

	if (n1 * x[1] < 5.0 || n1 * (1.0 - x[1]) < 5.0) {
	    infobox(_("The assumption of a normal sampling distribution\n"
		      "is not justified here.  Abandoning the test."));
	    return;
	}
	break;

    case TWO_MEANS:
	x[0] = getval(test->entry[k++], C_DBL); /* mean1 */
	if (na(x[0])) return;
	x[1] = getval(test->entry[k++], C_POS_DBL); /* sd1 */
	if (na(x[1])) return;
	n1 = getval(test->entry[k++], C_POS_INT);
	if (n1 < 0) return;

	x[2] = getval(test->entry[k++], C_DBL); /* mean2 */
	if (na(x[2])) return;
	x[3] = getval(test->entry[k++], C_POS_DBL); /* sd2 */
	if (na(x[3])) return;
	n2 = getval(test->entry[k++], C_POS_INT);
	if (n2 < 0) return;

	x[4] = getval(test->entry[k], C_DBL);
	if (na(x[4])) return;
	break;

    case TWO_VARIANCES:
	x[0] = getval(test->entry[k++], C_POS_DBL);
	if (na(x[0])) return;
	n1 = getval(test->entry[k++], C_POS_INT);
	if (n1 < 0) return;

	x[1] = getval(test->entry[k++], C_POS_DBL);
	if (na(x[1])) return;
	n2 = getval(test->entry[k], C_POS_INT);
	if (n2 < 0) return;
	break;

    case TWO_PROPNS:
	x[0] = getval(test->entry[k++], C_FRAC);
	if (na(x[0])) return;
	n1 = getval(test->entry[k++], C_POS_INT);
	if (n1 < 0) return;

	x[1] = getval(test->entry[k++], C_FRAC);
	if (na(x[1])) return;
	n2 = getval(test->entry[k], C_POS_INT);
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

static void dist_graph_callback (GtkWidget *w, CalcChild *child)
{
    dist_t **tabs = child->calcp;
    double x[2] = {0};
    int i, d;

    i = gtk_notebook_get_current_page(GTK_NOTEBOOK(child->book));
    d = dist_from_page(child->code, i);

    if (d == NORMAL_DIST && tabs[i]->flags) {
	plot_cdf(child->dlg);
	return;
    }

    if (get_dist_entry_vector(child->code, tabs[i], d, x, NULL)) {
	return;
    }

    if (child->plot != NULL) {
	revise_distribution_plot(child->plot, d, x);
    } else {
	dist_graph(d, x);
    }
}

static void toggle_dist_flag (GtkToggleButton *b,
			      dist_t *dist)
{
    int i;

    dist->flags = gtk_toggle_button_get_active(b);

    for (i=0; i<NDISTENTRY; i++) {
	if (dist->entry[i] == NULL) {
	    break;
	}
	gtk_widget_set_sensitive(dist->entry[i],
				 dist->flags == 0);
    }
}

static void
calc_checkbox (GtkWidget *tbl, gint *rows, const gchar *label,
	       CalcChild *child, int i)
{
    dist_t **dist = child->calcp;
    GtkWidget *w;

    *rows += 1;

    gtk_table_resize(GTK_TABLE(tbl), *rows, 2);
    w = gtk_check_button_new_with_label(_(label));
    gtk_table_attach_defaults(GTK_TABLE(tbl),
			      w, 0, 1, *rows - 1, *rows);
    gtk_widget_show(w);
    dist[i]->check = w;

    g_signal_connect(G_OBJECT(w), "toggled",
		     G_CALLBACK(toggle_dist_flag),
		     dist[i]);
}

static void
calc_entry_with_default (GtkWidget *tbl, gint *rows,
			 const gchar *label, CalcChild *child,
			 int i, const char *deflt)
{
    dist_t **dist = child->calcp;
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
    dist[i]->entry[*rows-2] = tmp;

    /* default text */
    if (deflt != NULL) {
	gtk_entry_set_text(GTK_ENTRY(tmp), deflt);
    }

    g_signal_connect(G_OBJECT(tmp), "activate", child->callback,
		     child);
}

static void calc_entry (GtkWidget *tbl, gint *rows,
			const gchar *label, CalcChild *child,
			int i)
{
    calc_entry_with_default(tbl, rows, label, child, i, NULL);
}

/* make a tab (notebook page) for a given distribution */

static void make_dist_tab (CalcChild *child, int i)
{
    GtkWidget *tmp, *box, *tbl;
    gint rows = 1;
    const gchar *titles[] = {
	N_("uniform"),
	N_("normal"),
	N_(" t "),
	N_("chi-square"),
	N_(" F "),
	N_("gamma"),
	N_("binomial"),
	N_("poisson"),
	N_("weibull"),
	N_(" DW "),
    };
    int d = dist_from_page(child->code, i);

    box = gtk_vbox_new(FALSE, 0);
    gtk_container_set_border_width(GTK_CONTAINER(box), 10);
    gtk_widget_show(box);

    tmp = gtk_label_new(_(titles[d]));
    gtk_widget_show(tmp);
    gtk_notebook_append_page(GTK_NOTEBOOK(child->book), box, tmp);

    tbl = gtk_table_new(rows, 2, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(tbl), 5);
    gtk_table_set_col_spacings(GTK_TABLE(tbl), 5);
    gtk_box_pack_start(GTK_BOX(box), tbl, FALSE, FALSE, 0);
    gtk_widget_show(tbl);

    switch (d) {
    case UNIFORM_DIST:
	calc_entry_with_default(tbl, &rows, N_("minimum"),
				child, i, "0");
	calc_entry_with_default(tbl, &rows, N_("maximum"),
				child, i, "1");
	break;
    case NORMAL_DIST:
	calc_entry_with_default(tbl, &rows, N_("mean"),
				child, i, "0");
	calc_entry_with_default(tbl, &rows, N_("std. deviation"),
				child, i, "1");
	break;
    case T_DIST:
	calc_entry(tbl, &rows, N_("df"), child, i);
	break;
    case CHISQ_DIST:
	calc_entry(tbl, &rows, N_("df"), child, i);
	break;
    case F_DIST:
	calc_entry(tbl, &rows, N_("dfn"), child, i);
	calc_entry(tbl, &rows, N_("dfd"), child, i);
	break;
    case GAMMA_DIST:
    case WEIBULL_DIST:
	calc_entry(tbl, &rows, N_("shape"), child, i);
	calc_entry(tbl, &rows, N_("scale"), child, i);
	break;
    case BINOMIAL_DIST:
	calc_entry(tbl, &rows, N_("Prob"), child, i);
	calc_entry(tbl, &rows, N_("trials"), child, i);
	break;
    case POISSON_DIST:
	calc_entry(tbl, &rows, N_("mean"), child, i);
	break;
    case DW_DIST:
	calc_entry(tbl, &rows, N_("sample size, n"), child, i);
	calc_entry(tbl, &rows, N_("number of regressors\n(excluding the constant)"),
		   child, i);
	break;
    default:
	break;
    }

    if (child->code == CALC_DIST) {
	if (d != DW_DIST) {
	    calc_entry(tbl, &rows, N_("right-tail probability"), child, i);
	}
    } else if (child->code == CALC_PVAL) {
	calc_entry(tbl, &rows, N_("value"), child, i);
    } else if (child->code == CALC_RAND) {
	calc_entry(tbl, &rows, N_("name"), child, i);
    } else if (child->code == CALC_GRAPH) {
	if (d == NORMAL_DIST) {
	    calc_checkbox(tbl, &rows, N_("CDF instead of density"), child, i);
	}
    }
}

static int get_restriction_vxy (const char *s, int *vx, int *vy,
				GretlOp *yop, double *yval)
{
    char test[VNAMELEN];
    char *p, *q = NULL;
    char *str = g_strdup(s);
    int err = 0;

    /* FIXME : use genr here */

    if (str == NULL) {
	return 1;
    }

    p = strchr(str, '(');
    *p = 0;
    p++;

    if (gretl_scan_varname(str, test) != 1) {
	err = 1;
    } else {
	*vx = series_index(dataset, test);
	if (*vx >= dataset->v) {
	    gui_errmsg(E_UNKVAR);
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
            } else if (!strncmp(q, "==", 2)) {
                *yop = OP_EQ;
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
	if (gretl_scan_varname(p, test) != 1) {
	    err = 1;
	} else {
	    *vy = series_index(dataset, test);
	    if (*vy >= dataset->v) {
		gui_errmsg(E_UNKVAR);
		err = 1;
	    }
	}
    }

    if (!err) {
        q += strspn(q, " ");
        if (!strncmp(q, "TRUE", 4)) {
            *yval = 1;
        } else if (!strncmp(q, "FALSE", 5)) {
            *yval = 0;
        } else if (sscanf(q, "%lf", yval) != 1) {
	    err = 1;
	}
    }

    g_free(str);

    return err;
}

static void entry_set_float (test_t *test, int i, double x)
{
    char numstr[32];

    sprintf(numstr, "%.10g", x);
    gtk_entry_set_text(GTK_ENTRY(test->entry[i]), numstr);
}

static void entry_set_int (test_t *test, int i, int k)
{
    char numstr[16];

    sprintf(numstr, "%d", k);
    gtk_entry_set_text(GTK_ENTRY(test->entry[i]), numstr);
}

static void entry_set_blank (test_t *test, int i)
{
    gtk_entry_set_text(GTK_ENTRY(test->entry[i]), "");
}

/* fill out the sample statistics boxes based on the user's
   choice of variable (or variable plus restriction) */

static void populate_stats (GtkWidget *w, gpointer p)
{
    test_t *test = g_object_get_data(G_OBJECT(p), "test");
    int pos = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(p), "pos"));
    gchar **pbuf = g_object_get_data(G_OBJECT(p), "pbuf");
    int t, n = sample_size(dataset);
    int vx = -1, vy = -1;
    GretlOp yop = 0;
    gchar *buf = NULL;
    double x1, x2, yval;
    int err = 0;

    g_return_if_fail(GTK_IS_COMBO_BOX(p));
    if (!gtk_widget_is_sensitive(p)) {
	return;
    }

    buf = combo_box_get_active_text(p);
    if (buf == NULL || *buf == '\0') {
	g_free(buf);
	return;
    }

    if (pbuf != NULL) {
	if (*pbuf != NULL && !strcmp(buf, *pbuf)) {
	    /* no real change */
	    g_free(buf);
	    return;
	}
	g_free(*pbuf);
	*pbuf = g_strdup(buf);
    }

    if (strchr(buf, '(') != NULL) {
	/* e.g. "cholest (gender = 1)" */
	err = get_restriction_vxy(buf, &vx, &vy, &yop, &yval);
    } else {
	vx = series_index(dataset, buf);
	if (vx >= dataset->v) {
	    err = 1;
	}
    }

    g_free(buf);

    if (err) {
	GtkWidget *entry;

	/* scrub any existing stats entries */
	entry_set_blank(test, pos);
	entry_set_blank(test, pos + 1);
	if (test->code == ONE_MEAN || test->code == TWO_MEANS) {
	    entry_set_blank(test, pos + 2);
	}
	/* highlight error region */
	entry = gtk_bin_get_child(GTK_BIN(p));
	gtk_editable_select_region(GTK_EDITABLE(entry), 0, -1);
	return;
    }

    for (t=dataset->t1; t<=dataset->t2; t++) {
	if (na(dataset->Z[vx][t]) ||
	    (vy > 0 && !eval_ytest(dataset->Z[vy][t], yop, yval))) {
	    n--;
	}
    }

    if (n == 0) {
	errbox_printf(_("Data missing for variable '%s'"), dataset->varname[vx]);
	return;
    }

    if (test->code == ONE_MEAN || test->code == TWO_MEANS) {
	if (vy < 0) {
	    x1 = gretl_mean(dataset->t1, dataset->t2, dataset->Z[vx]);
	    x2 = gretl_stddev(dataset->t1, dataset->t2, dataset->Z[vx]);
	} else {
	    x1 = gretl_restricted_mean(dataset->t1, dataset->t2,
				       dataset->Z[vx], dataset->Z[vy],
				       yop, yval);
	    x2 = gretl_restricted_stddev(dataset->t1, dataset->t2,
					 dataset->Z[vx], dataset->Z[vy],
					 yop, yval);
	}
	entry_set_float(test, pos, x1);
	entry_set_float(test, pos + 1, x2);
	entry_set_int(test, pos + 2, n);
    } else if (test->code == ONE_VARIANCE || test->code == TWO_VARIANCES) {
	if (vy < 0) {
	    x1 = gretl_variance(dataset->t1, dataset->t2, dataset->Z[vx]);
	} else {
	    x1 = gretl_restricted_variance(dataset->t1, dataset->t2,
					   dataset->Z[vx], dataset->Z[vy],
					   yop, yval);
	}
	entry_set_float(test, pos, x1);
	entry_set_int(test, pos + 1, n);
    } else if (test->code == ONE_PROPN || test->code == TWO_PROPNS) {
	if (vy < 0) {
	    x1 = gretl_mean(dataset->t1, dataset->t2, dataset->Z[vx]);
	} else {
	    x1 = gretl_restricted_mean(dataset->t1, dataset->t2,
				       dataset->Z[vx], dataset->Z[vy],
				       yop, yval);
	}
	entry_set_float(test, pos, x1);
	entry_set_int(test, pos + 1, n);
    }
}

static int var_is_ok (int i, test_t *test)
{
    int ret = 1;

    if (series_is_hidden(dataset, i)) {
	ret = 0;
    } else {
	int need_dummy = test->category == CALC_TEST &&
	    (test->code == ONE_PROPN || test->code == TWO_PROPNS);

	if (need_dummy) {
	    ret = gretl_isdummy(dataset->t1, dataset->t2, dataset->Z[i]);
	}
    }

    return ret;
}

static void add_vars_to_combo (GtkWidget *box, test_t *test, int pos)
{
    int i, vmin = (pos > 0)? 2 : 1;

    for (i=vmin; i<dataset->v; i++) {
	if (var_is_ok(i, test)) {
	    combo_box_append_text(box, dataset->varname[i]);
	}
    }

    if (pos > 0) {
	/* add first variable at the end of the list */
	for (i=1; i<dataset->v; i++) {
	    if (var_is_ok(i, test)) {
		combo_box_append_text(box, dataset->varname[i]);
		break;
	    }
	}
    }

    gtk_combo_box_set_active(GTK_COMBO_BOX(box), 0);
}

static void toggle_combo_ok (GtkWidget *toggle, gpointer p)
{
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(toggle))) {
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
	    g_free(*pbuf);
	}
	free(pbuf);
	g_object_set_data(G_OBJECT(w), "pbuf", NULL);
    }
}

static void select_child_callback (GtkComboBox *b, gpointer p)
{
    if (gtk_combo_box_get_active(b) >= 0) {
	populate_stats(NULL, p);
    }
}

static void add_test_var_selector (GtkWidget *tbl, gint *row,
				   test_t *test, int i,
				   int labelit)
{
    GtkWidget *label, *tmp;
    gchar **pbuf;

    *row += 1;
    gtk_table_resize(GTK_TABLE(tbl), *row, 2);
    if (labelit) {
	gchar *tmp = g_strdup_printf(_("Variable %d"), i + 1);

	label = gtk_label_new(tmp);
	g_free(tmp);
    } else {
	label = gtk_label_new(_("Variable"));
    }
    gtk_table_attach_defaults(GTK_TABLE(tbl),
			      label, 0, 1, *row - 1, *row);
    gtk_widget_show(label);

    tmp = combo_box_text_new_with_entry();
    gtk_table_attach_defaults(GTK_TABLE(tbl),
			      tmp, 1, 2, *row - 1, *row);
    gtk_widget_show(tmp);
    g_object_set_data(G_OBJECT(tmp), "test", test);
    g_object_set_data(G_OBJECT(tmp), "pos", GINT_TO_POINTER(i));
    test->entry[i] = gtk_bin_get_child(GTK_BIN(tmp));

    pbuf = malloc(sizeof *pbuf);
    if (pbuf != NULL) {
	*pbuf = NULL;
	g_object_set_data(G_OBJECT(tmp), "pbuf", pbuf);
	g_signal_connect(G_OBJECT(tmp), "destroy", G_CALLBACK(free_pbuf), NULL);
    }

    if (i > 0) {
	test->combo[1] = tmp;
    } else {
	test->combo[0] = tmp;
    }

    add_vars_to_combo(tmp, test, i);
}

static void add_test_combo (GtkWidget *tbl, gint *rows,
			    test_t *test, int pos)
{
    GtkWidget *button, *tmp;
    GtkWidget *entry;
    gchar **pbuf;

    *rows += 1;
    gtk_table_resize(GTK_TABLE(tbl), *rows, 2);
    button = gtk_check_button_new_with_label(_("Use variable from dataset"));
    gtk_table_attach_defaults(GTK_TABLE(tbl),
			      button, 0, 1, *rows - 1, *rows);
    gtk_widget_show(button);

    tmp = combo_box_text_new_with_entry();
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

    add_vars_to_combo(tmp, test, pos);
    gtk_widget_set_sensitive(tmp, FALSE);

    entry = gtk_bin_get_child(GTK_BIN(tmp));
    g_signal_connect(G_OBJECT(GTK_ENTRY(entry)), "key-press-event",
		     G_CALLBACK(catch_combo_key), tmp);
    g_signal_connect(G_OBJECT(GTK_COMBO_BOX(tmp)), "changed",
		     G_CALLBACK(select_child_callback), tmp);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(toggle_combo_ok), tmp);
}

static void test_entry (GtkWidget *tbl, gint *row,
			const gchar *label, test_t *test,
			int i)
{
    GtkWidget *tmp;

    *row += 1;
    gtk_table_resize(GTK_TABLE(tbl), *row, 2);
    tmp = gtk_label_new(label);
    gtk_misc_set_alignment(GTK_MISC(tmp), 1, 0.5);
    gtk_table_attach_defaults(GTK_TABLE(tbl),
			      tmp, 0, 1, *row - 1, *row);
    gtk_widget_show(tmp);
    tmp = gtk_entry_new();
    gtk_table_attach_defaults(GTK_TABLE(tbl),
			      tmp, 1, 2, *row - 1, *row);
    gtk_widget_show(tmp);
    test->entry[i] = tmp;

    g_signal_connect(G_OBJECT(tmp), "activate",
		     G_CALLBACK(h_test), test);
}

static void add_test_label (GtkWidget *tbl, gint *row,
			    const gchar *label)
{
    GtkWidget *tmp;

    *row += 1;
    gtk_table_resize(GTK_TABLE(tbl), *row, 2);
    tmp = gtk_label_new(label);
    gtk_misc_set_alignment(GTK_MISC(tmp), 0, 0.5);
    gtk_table_attach_defaults(GTK_TABLE(tbl),
			      tmp, 0, 2, *row - 1, *row);
    gtk_widget_show(tmp);
}

static void add_test_check (GtkWidget *tbl, gint *row,
			    const gchar *label, test_t *test,
			    gboolean val)
{
    GtkWidget *tmp;

    *row += 1;
    gtk_table_resize(GTK_TABLE(tbl), *row, 2);
    tmp = gtk_check_button_new_with_label(label);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), val);
    gtk_table_attach_defaults(GTK_TABLE(tbl),
			      tmp, 0, 2, *row - 1, *row);
    gtk_widget_show(tmp);
    test->check = tmp;
}

static int n_ok_series (void)
{
    int i, nv = 0;

    if (dataset != NULL) {
	for (i=1; i<dataset->v; i++) {
	    if (!series_is_hidden(dataset, i)) {
		nv++;
	    }
	}
    }

    return nv;
}

static int n_ok_dummies (void)
{
    int i, nv = 0;

    if (dataset != NULL) {
	for (i=1; i<dataset->v; i++) {
	    if (!series_is_hidden(dataset, i) &&
		gretl_isdummy(dataset->t1, dataset->t2, dataset->Z[i])) {
		nv++;
	    }
	}
    }

    return nv;
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
	N_("Runs test"),
	N_("Correlation")
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

    case NP_CORR:
	add_test_var_selector(tbl, &rows, test, 0, 1);
	add_test_var_selector(tbl, &rows, test, 1, 1);

	/* option radios */
	rows += 2;
	gtk_table_resize(GTK_TABLE(tbl), rows, 2);

	test->radio[0] = gtk_radio_button_new_with_label(NULL, _("Kendall's tau"));
	gtk_table_attach_defaults(GTK_TABLE(tbl), test->radio[0], 0, 2,
				  rows - 2, rows - 1);
	gtk_widget_show(test->radio[0]);

	group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(test->radio[0]));
	test->radio[1] = gtk_radio_button_new_with_label(group, _("Spearman's rho"));
	gtk_table_attach_defaults(GTK_TABLE(tbl), test->radio[1], 0, 2,
				  rows - 1, rows);
	gtk_widget_show(test->radio[1]);
	break;

    default:
	break;
    }

    /* check box for extra option */
    rows += 1;
    gtk_table_resize(GTK_TABLE(tbl), rows, 2);
    if (idx == NP_DIFF || idx == NP_CORR) {
	tmp = gtk_check_button_new_with_label(_("Show details"));
    } else {
	tmp = gtk_check_button_new_with_label(_("Use first difference"));
    }
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 2,
			      rows - 1, rows);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), FALSE);

    gtk_widget_show(tmp);
    test->extra = tmp;

    if (idx == NP_RUNS) {
	rows += 1;
	gtk_table_resize(GTK_TABLE(tbl), rows, 2);
	tmp = gtk_check_button_new_with_label(_("Assume positive and negative are equiprobable"));
	gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 2,
				  rows - 1, rows);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), FALSE);
	gtk_widget_show(tmp);
	test->check = tmp;
    }

    if (idx == NP_DIFF) {
	gtk_widget_set_sensitive(tmp, FALSE);
	desensitize_conditional_on(tmp, test->radio[0]);
    }
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
	test_entry(tbl, &rows, _("sample mean"), test, 0);
	test_entry(tbl, &rows, _("std. deviation"), test, 1);
	test_entry(tbl, &rows, _("sample size"), test, 2);
	test_entry(tbl, &rows, _("H0: mean ="), test, 3);
	add_test_check(tbl, &rows, _("Assume standard deviation is "
				     "population value"), test, FALSE);
	break;

    case ONE_VARIANCE:
	test_entry(tbl, &rows, _("sample variance"), test, 0);
	test_entry(tbl, &rows, _("sample size"), test, 1);
	test_entry(tbl, &rows, _("H0: variance ="), test, 2);
	break;

    case ONE_PROPN: /* proportion */
	test_entry(tbl, &rows, _("sample proportion"), test, 0);
	test_entry(tbl, &rows, _("sample size"), test, 1);
	test_entry(tbl, &rows, _("H0: proportion ="), test, 2);
	break;

    case TWO_MEANS:
	test_entry(tbl, &rows, _("mean of sample 1"), test, 0);
	test_entry(tbl, &rows, _("std. deviation, sample 1"), test, 1);
	test_entry(tbl, &rows, _("size of sample 1"), test, 2);
	if (nv > 0) {
	    add_test_combo(tbl, &rows, test, 3);
	}
	test_entry(tbl, &rows, _("mean of sample 2"), test, 3);
	test_entry(tbl, &rows, _("std. deviation, sample 2"), test, 4);
	test_entry(tbl, &rows, _("size of sample 2"), test, 5);
	test_entry(tbl, &rows, _("H0: Difference of means ="), test, 6);
	gtk_entry_set_text(GTK_ENTRY(test->entry[6]), "0");
	add_test_check(tbl, &rows, _("Assume common population standard "
				     "deviation"), test, TRUE);
	break;

    case TWO_VARIANCES:
	test_entry(tbl, &rows, _("variance of sample 1"), test, 0);
	test_entry(tbl, &rows, _("size of sample 1"), test, 1);
	if (nv > 0) {
	    add_test_combo(tbl, &rows, test, 2);
	}
	test_entry(tbl, &rows, _("variance of sample 2"), test, 2);
	test_entry(tbl, &rows, _("size of sample 2"), test, 3);
	add_test_label(tbl, &rows, _("H0: Ratio of variances = 1"));
	break;

    case TWO_PROPNS:
	test_entry(tbl, &rows, _("proportion, sample 1"), test, 0);
	test_entry(tbl, &rows, _("size of sample 1"), test, 1);
	if (nv > 0) {
	    add_test_combo(tbl, &rows, test, 2);
	}
	test_entry(tbl, &rows, _("proportion, sample 2"), test, 2);
	test_entry(tbl, &rows, _("size of sample 2"), test, 3);
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
    int i;

    *wp = NULL;

    if (c == CALC_TEST || c == CALC_NPTEST) {
	test_t **test = child->calcp;

	for (i=0; i<child->n_pages; i++) {
	    free(test[i]);
	}
	free(test);
    } else {
	dist_t **dist = child->calcp;

	for (i=0; i<child->n_pages; i++) {
	    free(dist[i]);
	}
	free(dist);
    }

    free(child);
}

static test_t *test_holder_new (int c, int i)
{
    test_t *test = mymalloc(sizeof *test);

    if (test != NULL) {
	int j;

	test->category = c;
	test->code = i;
	test->check = NULL;
	test->extra = NULL;

	test->combo[0] = test->combo[1] = NULL;
	test->radio[0] = test->radio[1] = test->radio[2] = NULL;

	for (j=0; j<NTESTENTRY; j++) {
	    test->entry[j] = NULL;
	}
    }

    return test;
}

static dist_t *dist_holder_new (void)
{
    dist_t *dist = mymalloc(sizeof *dist);

    if (dist != NULL) {
	int j;

	dist->flags = 0;
	dist->check = NULL;

	for (j=0; j<NDISTENTRY; j++) {
	    dist->entry[j] = NULL;
	}
    }

    return dist;
}

static int child_allocate_calcp (CalcChild *child)
{
    int c = child->code;
    int n = child->n_pages;
    int i, err = 0;

    child->calcp = NULL;

    if (c == CALC_TEST || c == CALC_NPTEST) {
	test_t **test = mymalloc(n * sizeof *test);

	if (test != NULL) {
	    child->calcp = test;
	    for (i=0; i<n && !err; i++) {
		test[i] = test_holder_new(c, i);
		if (test[i] == NULL) {
		    err = E_ALLOC;
		}
	    }
	}
    } else {
	dist_t **dist = mymalloc(n * sizeof *dist);

	if (dist != NULL) {
	    child->calcp = dist;
	    for (i=0; i<n && !err; i++) {
		dist[i] = dist_holder_new();
		if (dist[i] == NULL) {
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

    if (code == CALC_TEST) {
	child->n_pages = NTESTS;
	child->callback = G_CALLBACK(h_test_global);
    } else if (code == CALC_NPTEST) {
	child->n_pages = NPTESTS;
	child->callback = G_CALLBACK(np_test_global);
    } else if (code == CALC_PVAL) {
	child->n_pages = NPVAL;
	child->callback = G_CALLBACK(get_pvalue);
    } else if (code == CALC_DIST) {
	child->n_pages = NDISTS;
	child->callback = G_CALLBACK(get_critical);
    } else if (code == CALC_RAND) {
	child->n_pages = NRAND;
	child->callback = G_CALLBACK(get_random);
    } else {
	child->n_pages = NGRAPHS;
	child->callback = G_CALLBACK(dist_graph_callback);
    }

    child->plot = (code == CALC_GRAPH_ADD)? p : NULL;

    if (child_allocate_calcp(child)) {
	free(child);
	return NULL;
    }

    child->dlg = gretl_gtk_window();

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
    gtk_box_set_spacing(GTK_BOX(child->bbox), 10);

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

static void calc_disable_page (CalcChild *child, int i)
{
    GtkWidget *p = gtk_notebook_get_nth_page(GTK_NOTEBOOK(child->book), i);

    gtk_widget_set_sensitive(p, FALSE);
    p = gtk_notebook_get_tab_label(GTK_NOTEBOOK(child->book), p);
    gtk_widget_set_sensitive(p, FALSE);
}

static void
configure_graph_add_tabs (CalcChild *child, png_plot *plot)
{
    int i = current_graph_dist(plot);

    gtk_notebook_set_current_page(GTK_NOTEBOOK(child->book), i);
}

static void switch_child_role (GtkWidget *win, png_plot *plot)
{
    CalcChild *child;
    GtkNotebook *book;
    dist_t **dists, *dist;
    int i, d;

    child = g_object_get_data(G_OBJECT(win), "gchild");
    book = GTK_NOTEBOOK(child->book);
    child->plot = plot;
    child->code = CALC_GRAPH_ADD;

    gtk_window_set_title(GTK_WINDOW(win),
			 _("gretl: add distribution graph"));

    make_graph_window_transient(win, plot);

    i = gtk_notebook_get_current_page(book);
    d = dist_from_page(child->code, i);
    dists = child->calcp;
    dist = dists[i];

    if (d == NORMAL_DIST) {
	gtk_widget_set_sensitive(dist->check, FALSE);
    } else if (d == T_DIST || d == CHISQ_DIST || d == POISSON_DIST) {
	gtk_editable_select_region(GTK_EDITABLE(dist->entry[0]), 0, -1);
	gtk_widget_grab_focus(dist->entry[0]);
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

static void real_stats_calculator (int code, gpointer data)
{
    GtkWidget *tmp = NULL;
    static GtkWidget *winptr[CALC_RAND + 1];
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
	    warnbox(_("No suitable data are available"));
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

    for (i=0; i<child->n_pages; i++) {
	if (code == CALC_TEST) {
	    make_test_tab(child, i);
	} else if (code == CALC_NPTEST) {
	    make_nptest_tab(child, i);
	} else {
	    make_dist_tab(child, i);
	}
    }

    if (code == CALC_GRAPH_ADD) {
	configure_graph_add_tabs(child, data);
    }

    if (code == CALC_NPTEST && nv < 2) {
	calc_disable_page(child, NP_DIFF);
	calc_disable_page(child, NP_CORR);
	gtk_notebook_set_current_page(GTK_NOTEBOOK(child->book), NP_RUNS);
    }

    /* Close button */
    tmp = gtk_button_new_from_stock(GTK_STOCK_CLOSE);
    gtk_widget_set_can_default(tmp, TRUE);
    gtk_container_add(GTK_CONTAINER(child->bbox), tmp);
    g_signal_connect_swapped(G_OBJECT(tmp), "clicked",
			     G_CALLBACK(gtk_widget_destroy),
			     child->dlg);

    /* OK button */
    tmp = gtk_button_new_from_stock(GTK_STOCK_OK);
    gtk_widget_set_can_default(tmp, TRUE);
    gtk_container_add(GTK_CONTAINER(child->bbox), tmp);
    g_signal_connect(G_OBJECT(tmp), "clicked", child->callback, child);

    /* Help button? */
    hcode = calc_help_code(code);
    if (hcode) {
	context_help_button(child->bbox, hcode);
    }

    gtk_widget_show_all(child->dlg);
    window_list_add(child->dlg, STAT_TABLE);
}

static gchar *maybe_trim_y_equals (gchar *s)
{
    gchar *p = s;

    while (isspace(*p)) p++;

    if (*p == 'y' || *p == 'Y') {
	p++;
	while (isspace(*p)) p++;
	if (*p == '=') {
	    s = p + 1;
	}
    }

    return s;
}

/* for gnuplot: convert '^' to '**' for exponentiation */

static gchar *formula_mod (const gchar *s)
{
    const gchar *p = s;
    int n = strlen(s) + 1;
    gchar *q, *ret;

    while (*p) {
	if (*p == '^') n++;
	p++;
    }

    ret = g_malloc(n);
    q = ret;

    while (*s) {
	if (*s == '^') {
	    *q++ = '*';
	    *q++ = '*';
	    s++;
	} else {
	    *q++ = *s++;
	}
    }

    *q = '\0';

    return ret;
}

struct curve_plotter {
    GtkWidget *dlg;
    GtkWidget *entry;
    gchar *formula;
    double xmin;
    double xmax;
};

static void do_plot_curve (GtkWidget *w, struct curve_plotter *p)
{
    gchar *s1, *s0 = get_genr_string(p->entry, NULL);
    FILE *fp = NULL;
    int err = 0;

    if (s0 == NULL || *s0 == '\0') {
	return;
    }

    /* It's "natural", but not accepted on a gnuplot 'plot'
       line, to type something like "y = x**2"; here we just
       want the bit to the right of the equals sign.
    */
    s1 = maybe_trim_y_equals(s0);

    g_free(p->formula);

    if (strchr(s1, '^')) {
	p->formula = formula_mod(s1);
    } else {
	p->formula = g_strdup(s1);
    }

    g_free(s0);

    fp = open_plot_input_file(PLOT_CURVE, 0, &err);
    if (err) {
	return;
    }

    print_keypos_string(GP_KEY_RIGHT_TOP, fp);

    gretl_push_c_numeric_locale();

    fprintf(fp, "set xrange [%g:%g]\n", p->xmin, p->xmax);
    fprintf(fp, "plot \\\n");
    if (strstr(p->formula, " with ") ||
	strstr(p->formula, " w ") ||
	strstr(p->formula, "title ")) {
	fprintf(fp, "%s\n", p->formula);
    } else {
	fprintf(fp, "%s notitle w lines\n", p->formula);
    }

    gretl_pop_c_numeric_locale();

    err = finalize_plot_input_file(fp);
    gui_graph_handler(err);

    if (!err) {
	gtk_widget_destroy(p->dlg);
    }
}

/* plot a curve specified via formula (no data required) */

static void plot_a_curve (void)
{
    static struct curve_plotter plotter;
    GtkWidget *dialog, *hbox, *vbox;
    GtkWidget *button, *tmp;
    GtkAdjustment *adj;

    if (plotter.dlg != NULL) {
	gtk_window_present(GTK_WINDOW(plotter.dlg));
	return;
    }

    if (plotter.xmin == 0 && plotter.xmax == 0) {
	plotter.xmax = 10;
    }

    plotter.dlg = dialog =
	gretl_dialog_new(_("gretl: plot a curve"), mdata->main, 0);

    g_signal_connect(G_OBJECT(dialog), "destroy",
		     G_CALLBACK(gtk_widget_destroyed), &plotter.dlg);

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));

    /* gnuplot formula entry box */
    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new(_("formula"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    plotter.entry = tmp = gtk_entry_new();
    gtk_entry_set_width_chars(GTK_ENTRY(tmp), 32);
    if (plotter.formula != NULL) {
	gtk_entry_set_text(GTK_ENTRY(tmp), plotter.formula);
    }
    gtk_entry_set_activates_default(GTK_ENTRY(tmp), TRUE);
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    /* x-axis spin buttons (min and range) */
    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new(_("x minimum"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    adj = (GtkAdjustment *) gtk_adjustment_new(plotter.xmin, -100, 100, 1, 0, 0);
    tmp = gtk_spin_button_new(adj, 1, 0);
    g_signal_connect(tmp, "value-changed",
		     G_CALLBACK(set_double_from_spin), &plotter.xmin);
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 0);
    tmp = gtk_label_new("  ");
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    tmp = gtk_label_new(_("x maximum"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
    adj = (GtkAdjustment *) gtk_adjustment_new(plotter.xmax, -100, 1000, 1, 0, 0);
    tmp = gtk_spin_button_new(adj, 1, 0);
    g_signal_connect(tmp, "value-changed",
		     G_CALLBACK(set_double_from_spin), &plotter.xmax);
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));

    /* "Cancel" button */
    cancel_delete_button(hbox, dialog);

    /* "OK" button */
    button = ok_button(hbox);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(do_plot_curve), &plotter);
    gtk_widget_grab_default(button);

    /* Help button */
    context_help_button(hbox, GPT_CURVE);

    gtk_widget_show_all(dialog);
}

static void do_plot_cdf (GtkWidget *w, GtkWidget *dlg)
{
    const char *formulae[] = {
	"normcdf(x)=0.5+0.5*erf(x/sqrt(2.0))",
	"logcdf(x)=1.0/(1+exp(-x))"
    };
    const char *titles[] = {
	N_("normal CDF"),
	N_("logistic CDF")
    };
    FILE *fp = NULL;
    double xmax = 4.0;
    int opt, err = 0;

    fp = open_plot_input_file(PLOT_CURVE, 0, &err);
    if (err) {
	return;
    }

    opt = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(dlg), "opt"));

    print_keypos_string(GP_KEY_LEFT_TOP, fp);

    if (opt > 0) {
	xmax = 6.0;
    }

    gretl_push_c_numeric_locale();

    if (opt == 0 || opt == 1) {
	fputs("# literal lines = 2\n", fp);
	fprintf(fp, "%s\n", formulae[opt]);
	fputs("set zeroaxis\n", fp);
    } else {
	fputs("# literal lines = 3\n", fp);
	fprintf(fp, "%s\n", formulae[0]);
	fprintf(fp, "%s\n", formulae[1]);
	fputs("set zeroaxis\n", fp);
    }

    fprintf(fp, "set xrange [%g:%g]\n", -xmax, xmax);
    fputs("plot \\\n", fp);

    if (opt == 0) {
	fprintf(fp, "normcdf(x) title \"%s\" w lines\n", _(titles[opt]));
    } else if (opt == 1) {
	fprintf(fp, "logcdf(x) title \"%s\" w lines\n", _(titles[opt]));
    } else {
	fprintf(fp, "normcdf(x) title \"%s\" w lines , \\\n", _(titles[0]));
	fprintf(fp, "logcdf(x) title \"%s\" w lines\n", _(titles[1]));
    }

    gretl_pop_c_numeric_locale();

    err = finalize_plot_input_file(fp);
    gui_graph_handler(err);

    gtk_widget_destroy(dlg);
}

static void set_cdf_opt (GtkWidget *button, GtkWidget *dlg)
{
    gpointer p = g_object_get_data(G_OBJECT(button), "opt");

    g_object_set_data(G_OBJECT(dlg), "opt", p);
}

static void plot_cdf (GtkWidget *parent)
{
    static GtkWidget *dialog;
    GtkWidget *hbox, *vbox;
    GtkWidget *button;
    GSList *group;

    if (dialog) {
	gtk_window_present(GTK_WINDOW(dialog));
	return;
    }

    dialog = gretl_dialog_new(_("gretl: plot CDF"), parent, 0);

    g_signal_connect(G_OBJECT(dialog), "destroy",
		     G_CALLBACK(gtk_widget_destroyed), &dialog);

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));

    /* normal */
    button = gtk_radio_button_new_with_label(NULL, _("standard normal"));
    pack_in_hbox(button, vbox, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON (button), TRUE);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_cdf_opt), dialog);

    /* logistic */
    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
    button = gtk_radio_button_new_with_label(group, _("logistic"));
    g_object_set_data(G_OBJECT(button), "opt", GINT_TO_POINTER(1));
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_cdf_opt), dialog);
    pack_in_hbox(button, vbox, 0);

    /* both */
    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));
    button = gtk_radio_button_new_with_label(group, _("both CDFs"));
    g_object_set_data(G_OBJECT(button), "opt", GINT_TO_POINTER(2));
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_cdf_opt), dialog);
    pack_in_hbox(button, vbox, 0);

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));

    /* "Cancel" button */
    cancel_delete_button(hbox, dialog);

    /* "OK" button */
    button = ok_button(hbox);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(do_plot_cdf), dialog);
    gtk_widget_grab_default(button);

    gtk_widget_show_all(dialog);
}

static int stats_calculator_code (GtkAction *action)
{
    const gchar *s = gtk_action_get_name(action);

    if (!strcmp(s, "PValues"))
	return CALC_PVAL;
    else if (!strcmp(s, "StatsTables"))
	return CALC_DIST;
    else if (!strcmp(s, "TestStats"))
	return CALC_TEST;
    else if (!strcmp(s, "NonparamTests"))
	return CALC_NPTEST;
    else if (!strcmp(s, "DistGraphs"))
	return CALC_GRAPH;
    else if (!strcmp(s, "AddRandom"))
	return CALC_RAND;
    else if (!strcmp(s, "PlotCurve"))
	return CALC_PLOT;
    else
	return 0;
}

void stats_calculator (GtkAction *action, gpointer data)
{
    int code = stats_calculator_code(action);

    g_return_if_fail(code == CALC_PVAL ||
		     code == CALC_DIST ||
		     code == CALC_TEST ||
		     code == CALC_NPTEST ||
		     code == CALC_GRAPH ||
		     code == CALC_RAND ||
		     code == CALC_PLOT);

    if (code == CALC_PLOT) {
	plot_a_curve();
    } else {
	real_stats_calculator(code, data);
    }
}

void dist_graph_add (gpointer p)
{
    real_stats_calculator(CALC_GRAPH_ADD, p);
}
