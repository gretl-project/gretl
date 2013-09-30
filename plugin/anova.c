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
#include "version.h"

#include <glib.h>

/* provides for one-way and two-way ANOVA */

struct anova {
    int n;          /* observations used */
    int nt;         /* number of treatment 'levels' */
    int nb;         /* number of 'blocks' (if applicable) */
    double SST;     /* total sum of squares */
    double SSTr;    /* treatment sum of squares */
    double SSB;     /* 'block' sum of squares (if applicable) */
    double SSE;     /* error sum of squares */
    double F;       /* F-test */
    double pval;    /* p-value for F-test */
    double *cmeans; /* column means */
    double *rmeans; /* row means */
    int *ccount;    /* column counts */
    int *rcount;    /* row counts */
    double *tvec;   /* workspace follows */
    double *bvec;
    gretl_matrix *tvals;
    gretl_matrix *bvals;
};

static void anova_init (struct anova *v)
{
    v->n = v->nt = v->nb = 0;
    v->SST = v->SSTr = v->SSB = v->SSE = 0.0;
    v->F = v->pval = NADBL;
    v->cmeans = v->rmeans = NULL;
    v->ccount = v->rcount = NULL;
    v->tvec = v->bvec = NULL;
    v->tvals = v->bvals = NULL;
}

static void anova_free (struct anova *v)
{
    free(v->cmeans);
    free(v->rmeans);
    free(v->ccount);
    free(v->rcount);
    free(v->tvec);
    free(v->bvec);

    gretl_matrix_free(v->tvals);
    gretl_matrix_free(v->bvals);
}

static int print_anova (struct anova *v, PRN *prn)
{
    int dftotal, dftreat, dfblock, dferr;
    double mst, msr, mse;
    int n, c1, c2, c3;

    dftotal = v->n - 1;
    dftreat = v->nt - 1;
    dfblock = (v->nb > 0)? v->nb - 1 : 0;
    dferr = dftotal - dftreat - dfblock;

    pputs(prn, "\n\n");

    c1 = g_utf8_strlen(_("Sum of squares"), -1);
    c2 = g_utf8_strlen(_("df"), -1);
    c3 = g_utf8_strlen(_("Mean square"), -1);

    c1 = (c1 < 35)? 35 : c1;
    c2 = (c2 > 8)? c2 + 1 : (c2 < 8)? 8 : c2;
    c3 = (c3 > 16)? c3 + 1 : (c3 < 16)? 16 : c3;

    /* header strings are right-aligned */
    n = g_utf8_strlen(_("Sum of squares"), -1);
    bufspace(c1 - n, prn);
    pputs(prn, _("Sum of squares"));
    n = g_utf8_strlen(_("df"), -1);
    bufspace(c2 + 1 - n, prn);
    pputs(prn, _("df"));
    n = g_utf8_strlen(_("Mean square"), -1);
    bufspace(c3 + 1 - n, prn);
    pputs(prn, _("Mean square"));
    pputs(prn, "\n\n");
    c1 = 16;

    /* Mean Square, treatment */
    msr = v->SSTr / dftreat;
    /* string left-aligned with initial offset of 2 */
    n = g_utf8_strlen(_("Treatment"), -1);
    bufspace(2, prn);
    pputs(prn, _("Treatment"));
    bufspace(16 - n, prn);	
    pprintf(prn, " %*g %*d %*g\n", c1, v->SSTr, c2, dftreat, c3, msr);

    if (dfblock > 0) {
	/* Mean Square, block */
	double msb = v->SSB / dfblock;

	/* string left-aligned with initial offset of 2 */
	n = g_utf8_strlen(_("Block"), -1);
	bufspace(2, prn);
	pputs(prn, _("Block"));
	bufspace(16 - n, prn);	
	pprintf(prn, " %*g %*d %*g\n", c1, v->SSB, c2, dfblock, c3, msb);
    }    

    /* Mean Square, errors */
    mse = v->SSE / dferr;
    /* string left-aligned with initial offset of 2 */
    n = g_utf8_strlen(_("Residual"), -1);
    bufspace(2, prn);
    pputs(prn, _("Residual"));
    bufspace(16 - n, prn);
    pprintf(prn, " %*g %*d %*g\n", c1, v->SSE, c2, dferr, c3, mse);

    /* Mean Square, total */
    mst = v->SST / dftotal;
    /* string left-aligned with initial offset of 2 */
    n = g_utf8_strlen(_("Total"), -1);
    bufspace(2, prn);
    pputs(prn, _("Total"));
    bufspace(16 - n, prn);
    pprintf(prn, " %*g %*d %*g\n", c1, v->SST, c2, dftotal, c3, mst);

    pputc(prn, '\n');

    if (na(v->F)) {
	pprintf(prn, "  F(%d, %d) = %g / %g (%s)\n\n", dftreat, dferr, 
		msr, mse, _("undefined"));
    } else {
	pprintf(prn, "  F(%d, %d) = %g / %g = %g ", 
		dftreat, dferr, msr, mse, v->F);
	if (v->pval < .0001) {
	    pprintf(prn, "[%s %.3g]\n\n", _("p-value"), v->pval);
	} else if (!na(v->pval)) {
	    pprintf(prn, "[%s %.4f]\n\n", _("p-value"), v->pval); 
	}
    }

    return 0;
}

/* one-way anova only: print the mean and standard deviation
   of the response at each level of treatment */

static int anova_print_means (struct anova *v, const double *xt,
			      const double *y, double ybar,
			      int t1, int t2, PRN *prn)
{
    double d, *csd = malloc(v->nt * sizeof *csd);
    int c1, c2, c3, c4;
    int i, t;

    if (csd == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<v->nt; i++) {
	csd[i] = 0.0;
    }

    for (t=t1; t<=t2; t++) {
	if (!na(xt[t]) && !na(y[t])) {
	    for (i=0; i<v->nt; i++) {
		if (xt[t] == v->tvals->val[i]) {
		    d = y[t] - v->cmeans[i];
		    csd[i] += d * d;
		    break;
		}
	    }
	}
    }

    c1 = g_utf8_strlen(_("Level"), -1);
    c2 = g_utf8_strlen(_("n"), -1);
    c3 = g_utf8_strlen(_("mean"), -1);
    c4 = g_utf8_strlen(_("std. dev"), -1);

    c1 = (c1 < 8)? 8 : c1;
    c2 = (c2 > 6)? c2 + 1 : (c2 < 6)? 6 : c2;
    c3 = (c3 > 10)? c3 + 1 : (c3 < 10)? 10 : c3;
    c4 = (c4 > 12)? c4 + 1 : (c4 < 12)? 12 : c4;

    pprintf(prn, "  %-*s %*s %*s %*s\n\n", c1, _("Level"), c2, _("n"), 
	    c3, _("mean"), c4, _("std. dev"));

    for (i=0; i<v->nt; i++) {
	if (v->ccount[i] > 1) {
	    csd[i] /= v->ccount[i] - 1;
	    csd[i] = sqrt(csd[i]);
	    pprintf(prn, "  %-*g %*d %*g %#*.5g\n", c1, v->tvals->val[i], 
		    c2, v->ccount[i], c3, v->cmeans[i], c4, csd[i]);
	} else {
	    pprintf(prn, "  %-*g %*d %*g %*s\n", c1, v->tvals->val[i], 
		    c2, v->ccount[i], c3, v->cmeans[i], c4, "NA");
	}
    } 

    pprintf(prn, "\n  %s = %g\n\n", _("Grand mean"), ybar);

    free(csd);

    return 0;
}

/* allocate arrays to hold the valid, in-sample values of the
   treatment and block variables */

static int anova_make_arrays (const double *xb, struct anova *v)
{
    int err = 0;

    v->tvec = malloc(v->n * sizeof *v->tvec);

    if (v->tvec == NULL) {
	err = E_ALLOC;
    } else if (xb != NULL) {
	v->bvec = malloc(v->n * sizeof *v->bvec);
	if (v->bvec == NULL) {
	    free(v->tvec);
	    v->tvec = NULL;
	    err = E_ALLOC;
	}
    }	

    return err;
}

/* construct vectors holding the distinct values of the 
   treatment and block variables */

static int anova_make_value_vecs (struct anova *v)
{
    int err = 0;

    v->tvals = gretl_matrix_values(v->tvec, v->n, OPT_S, &err);

    if (!err && v->tvals->rows < 2) {
	gretl_errmsg_set("Insufficient observations");
	err = E_DATA;
    }

    if (!err && v->bvec != NULL) {
	v->bvals = gretl_matrix_values(v->bvec, v->n, OPT_S, &err);
	if (!err && v->bvals->rows < 2) {
	    gretl_errmsg_set("Insufficient observations");
	    err = E_DATA;
	}
    }
    
    if (!err) {
	v->nt = v->tvals->rows;
	if (v->bvals != NULL) {
	    v->nb = v->bvals->rows;
	}
    }

    return err;
}

/* allocate and initialize arrays to calculate the column
   and row means of the ANOVA table */

static int anova_accounting_arrays (struct anova *v)
{
    int i;

    v->cmeans = malloc(v->nt * sizeof *v->cmeans);
    v->ccount = malloc(v->nt * sizeof *v->ccount);

    if (v->cmeans == NULL || v->ccount == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<v->nt; i++) {
	v->cmeans[i] = 0.0;
	v->ccount[i] = 0;
    }

    if (v->nb > 0) {
	v->rmeans = malloc(v->nb * sizeof *v->rmeans);
	v->rcount = malloc(v->nb * sizeof *v->rcount);

	if (v->rmeans == NULL || v->rcount == NULL) {
	    return E_ALLOC;
	}

	for (i=0; i<v->nb; i++) {
	    v->rmeans[i] = 0.0;
	    v->rcount[i] = 0;
	}
    }	

    return 0;
}

static void anova_add_F_stat (struct anova *v)
{
    int dfn = v->nt - 1;
    int dfd = v->n - 1 - dfn;
    double MSE, MSTr = v->SSTr / dfn;

    if (v->nb > 0) {
	dfd -= v->nb - 1;
    }

    MSE = v->SSE / dfd;

    if (MSE > 0) {
	v->F = MSTr / MSE;
	v->pval = snedecor_cdf_comp(dfn, dfd, v->F);
    }

    record_test_result(v->F, v->pval, NULL);    
}

#define anova_obs_ok(y,x,z,t) (!na(y[t]) && !na(x[t]) && \
                               (z == NULL || !na(z[t])))

/* For one-way anova @list contains response and treatment; for
   two-way it should in addition contain the block variable.
   @opt can contain OPT_Q to suppress printing.
*/

int gretl_anova (const int *list, const DATASET *dset, 
		 gretlopt opt, PRN *prn)
{
    struct anova v;
    const double *y, *xt, *xb;
    double ybar, dev;
    int i, t, t1, t2;
    int missvals = 0;
    int err = 0;

    if (list[0] < 2 || list[0] > 3) {
	return E_DATA;
    }

    anova_init(&v);

    t1 = dset->t1;
    t2 = dset->t2;

    list_adjust_sample(list, &t1, &t2, dset, &missvals);

    v.n = t2 - t1 + 1 - missvals;
    if (v.n < 2) {
	return E_TOOFEW;
    }

    y = dset->Z[list[1]];
    xt = dset->Z[list[2]];
    xb = (list[0] == 3)? dset->Z[list[3]] : NULL;

    /* check that treatment (and block, if present) are discrete */

    if (!series_is_discrete(dset, list[2]) && 
	!gretl_isdiscrete(t1, t2, xt)) {
	gretl_errmsg_set(_("anova: the treatment variable must be discrete"));
	return E_DATA;
    }

    if (xb != NULL && !series_is_discrete(dset, list[3]) && 
	!gretl_isdiscrete(t1, t2, xb)) {
	gretl_errmsg_set(_("anova: the block variable must be discrete"));
	return E_DATA;
    }

    v.n = 0;
    for (t=t1; t<=t2; t++) {
	if (anova_obs_ok(y, xt, xb, t)) {
	    v.n += 1;
	}
    }
    
    if (v.n < 2) {
	return E_TOOFEW;
    }

    err = anova_make_arrays(xb, &v);
    if (err) {
	return err;
    }

    /* fill tvec and bvec; calculate grand mean */

    ybar = 0.0;
    i = 0;
    for (t=t1; t<=t2; t++) {
	if (anova_obs_ok(y, xt, xb, t)) {
	    v.tvec[i] = xt[t];
	    ybar += y[t];
	    if (v.bvec != NULL) {
		v.bvec[i] = xb[t];
	    }
	    i++;
	}
    }

    ybar /= v.n;

    err = anova_make_value_vecs(&v);
    if (err) {
	goto bailout;
    }

    err = anova_accounting_arrays(&v);
    if (err) {
	goto bailout;
    }    

    /* find column (treatment) means */

    for (t=t1; t<=t2; t++) {
	if (anova_obs_ok(y, xt, xb, t)) {
	    dev = y[t] - ybar;
	    v.SST += dev * dev;
	    for (i=0; i<v.nt; i++) {
		if (xt[t] == v.tvals->val[i]) {
		    v.cmeans[i] += y[t];
		    v.ccount[i] += 1;
		    break;
		}
	    }
	}
    }

    for (i=0; i<v.nt; i++) {
	v.cmeans[i] /= v.ccount[i];
    }

    /* sums of squares */

    if (v.nb > 0) {
	/* two-way ANOVA */
	for (t=t1; t<=t2; t++) {
	    if (anova_obs_ok(y, xt, xb, t)) {
		for (i=0; i<v.nb; i++) {
		    if (xb[t] == v.bvals->val[i]) {
			v.rmeans[i] += y[t];
			v.rcount[i] += 1;
			break;
		    }
		}
	    }
	}

	for (i=0; i<v.nb; i++) {
	    v.rmeans[i] /= v.rcount[i];
	    dev = v.rmeans[i] - ybar;
	    v.SSB += dev * dev * v.rcount[i];
	}

	for (i=0; i<v.nt; i++) {
	    dev = v.cmeans[i] - ybar;
	    v.SSTr += dev * dev * v.ccount[i];
	}

	v.SSE = v.SST - v.SSTr - v.SSB;
    } else {
	/* one-way ANOVA */
	for (t=t1; t<=t2; t++) {
	    if (!na(xt[t]) && !na(y[t])) {
		for (i=0; i<v.nt; i++) {
		    if (xt[t] == v.tvals->val[i]) {
			dev = y[t] - v.cmeans[i];
			v.SSE += dev * dev;
			break;
		    }
		}
	    }
	}
	v.SSTr = v.SST - v.SSE;
    }

    anova_add_F_stat(&v);
	
    if (!(opt & OPT_Q)) {
	const char *yname = dset->varname[list[1]];
	const char *tname = dset->varname[list[2]];

	pputc(prn, '\n');
	pprintf(prn, _("%s, response = %s, treatment = %s:"), 
		_("Analysis of Variance"), yname, tname);

	err = print_anova(&v, prn);

	if (!err && v.nb == 0) {
	    anova_print_means(&v, xt, y, ybar, t1, t2, prn);
	}
    } 

 bailout:

    anova_free(&v);

    return err;
}
