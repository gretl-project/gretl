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
#include "matrix_extra.h"

#define ADEBUG 0

enum {
    DPD_TWOSTEP  = 1 << 0,
    DPD_ORTHDEV  = 1 << 1,
    DPD_TIMEDUM  = 1 << 2,
    DPD_WINCORR  = 1 << 3,
    DPD_NEWSTYLE = 1 << 4,
    DPD_SYSTEM   = 1 << 5,
    DPD_DPDSTYLE = 1 << 6
};

typedef struct dpdinfo_ dpdinfo;
typedef struct unit_info_ unit_info; /* old-style struct */

struct unit_info_ {
    int t1;      /* first usable obs for unit */
    int t2;      /* last usable obs */
    int nobs;    /* number of usable observations */
    char *skip;  /* mask for obs to be skipped, if any (1 = skip, 0 = OK) */ 
};

struct diag_info {
    int v;       /* ID number of variable */
    int minlag;  /* minimum lag order */
    int maxlag;  /* maximum lag order */
};

struct dpdinfo_ {
    int flags;            /* option flags */
    int step;             /* what step are we on? (1 or 2) */
    int yno;              /* ID number of dependent var */
    int p;                /* lag order for dependent variable */
    int qmax;             /* longest lag of y used as instrument */
    int qmin;             /* shortest lag of y used as instrument */
    int nx;               /* number of independent variables */
    int nzr;              /* number of regular instruments */
    int nzb;              /* number of block-diagonal instruments (other than y) */
    int nz;               /* total columns in instrument matrix */
    int pc0;              /* column in Z where predet vars start */
    int xc0;              /* column in Z where exog vars start */
    int N;                /* total number of units in sample */
    int effN;             /* number of units with usable observations */
    int T;                /* total number of observations per unit */
    int minTi;            /* minumum equations (> 0) for any given unit */
    int maxTi;            /* maximum equations for any given unit */
    int k;                /* number of parameters estimated */
    int nobs;             /* total observations actually used */
    int t1;               /* initial offset into dataset */
    int t1min;            /* first usable observation, any unit */
    int ndum;             /* number of time dummies to use */
    double SSR;           /* sum of squared residuals */
    double s2;            /* residual variance */
    double AR1;           /* z statistic for AR(1) errors */
    double AR2;           /* z statistic for AR(2) errors */
    double sargan;        /* overidentification test statistic */
    double wald;          /* Wald test statistic */
    int wdf;              /* degrees of freedom for Wald test */
    int *xlist;           /* list of independent variables */
    int *ilist;           /* list of instruments */
    gretl_matrix_block *B1; /* matrix holder */
    gretl_matrix_block *B2; /* matrix holder */
    gretl_matrix *beta;   /* parameter estimates */
    gretl_matrix *vbeta;  /* parameter variance matrix */
    gretl_matrix *uhat;   /* residuals, differenced version */
    gretl_matrix *H;      /* step 1 error covariance matrix */
    gretl_matrix *A;      /* \sum Z'_i H Z_i */
    gretl_matrix *Acpy;
    gretl_matrix *V;      /* covariance matrix */
    gretl_matrix *ZT;     /* transpose of full instrument matrix */
    gretl_matrix *Zi;     /* per-unit instrument matrix */
    gretl_matrix *Y;      /* transformed dependent var */
    gretl_matrix *X;      /* lagged differences of y, indep vars, etc. */
    gretl_matrix *tmp1;   /* workspace */
    gretl_matrix *kmtmp;
    gretl_matrix *kktmp;
    gretl_matrix *den;
    gretl_matrix *L1;
    gretl_matrix *XZA;
    gretl_matrix *R1;
    gretl_matrix *XZ;
    struct diag_info *d;   /* info on block-diagonal instruments */
    unit_info *ui;         /* old-style info on panel units */

    /* The members above should not be touched, since they are needed
       to support the "arbond" command. The members below are
       specific to the new "dpanel" approach, and can be changed
       as needed, as the code develops.
    */
									
    int ndiff;             /* total differenced observations */
    int nlev;              /* total levels observations */
    int max_ni;            /* max number of (possibly stacked) obs per unit */
    char *used;            /* global record of observations used */
    gretl_matrix *ZY;      /* cross-moment matrix */
    gretl_matrix *ZZ;      /* ditto */
    int *laglist;          /* (possibly discontinuous) list of lags */
};

#define data_index(dpd,i) (i * dpd->T + dpd->t1)

static void dpdinfo_free (dpdinfo *dpd)
{
    int i;

    if (dpd == NULL) {
	return;
    }

    gretl_matrix_block_destroy(dpd->B1);
    gretl_matrix_block_destroy(dpd->B2);

    gretl_matrix_free(dpd->V);

    free(dpd->xlist);
    free(dpd->ilist);
    free(dpd->laglist);

    free(dpd->d);
    free(dpd->used);

    if (dpd->ui != NULL) {
	for (i=0; i<dpd->N; i++) {
	    free(dpd->ui[i].skip);
	}
	free(dpd->ui);
    }

    free(dpd);
}

static int dpd_allocate_matrices (dpdinfo *dpd)
{
    int T = dpd->maxTi;

    dpd->ZY = dpd->ZZ = NULL;

    if (dpd->flags & DPD_NEWSTYLE) {
	/* temporary hack */
	dpd->uhat = dpd->ZT = dpd->H = dpd->A = NULL;
	dpd->Acpy = dpd->Zi = dpd->Y = dpd->X = NULL;

	dpd->tmp1 = dpd->kmtmp = dpd->kktmp = dpd->den = NULL;
	dpd->L1 = dpd->XZA = dpd->R1 = dpd->XZ = NULL;

	dpd->B1 = gretl_matrix_block_new(&dpd->beta,  dpd->k, 1,
					 &dpd->vbeta, dpd->k, dpd->k,
					 NULL);

	if (dpd->B1 == NULL) {
	    return E_ALLOC;
	} else {
	    gretl_matrix_zero(dpd->vbeta);
	    return 0;	
	}
    }

    dpd->B1 = gretl_matrix_block_new(&dpd->beta,  dpd->k, 1,
				     &dpd->vbeta, dpd->k, dpd->k,
				     &dpd->uhat,  dpd->nobs, 1,
				     &dpd->ZT,    dpd->nz, dpd->nobs,
				     &dpd->H,     T, T,
				     &dpd->A,     dpd->nz, dpd->nz,
				     &dpd->Acpy,  dpd->nz, dpd->nz,
				     &dpd->Zi,    T, dpd->nz,
				     &dpd->Y,     dpd->nobs, 1,
				     &dpd->X,     dpd->nobs, dpd->k,
				     NULL);
    if (dpd->B1 == NULL) {
	return E_ALLOC;
    }

    dpd->B2 = gretl_matrix_block_new(&dpd->tmp1,  dpd->nz, dpd->nz,
				     &dpd->kmtmp, dpd->k, dpd->nz,
				     &dpd->kktmp, dpd->k, dpd->k,
				     &dpd->den,   dpd->k, dpd->k,
				     &dpd->L1,    1, dpd->nz,
				     &dpd->XZA,   dpd->k, dpd->nz,
				     &dpd->R1,    dpd->nz, 1,
				     &dpd->XZ,    dpd->k, dpd->nz,
				     NULL);

    if (dpd->B2 == NULL) {
	return E_ALLOC;
    }    

    return 0;
}

/* if the const has been included among the regressors but not the
   instruments, add it to the instruments */

static int maybe_add_const_to_ilist (dpdinfo *dpd)
{
    int i, addc = 0;
    int err = 0;

    if (dpd->xlist == NULL || (dpd->nzr == 0 && dpd->nzb == 0)) {
	/* no x's, or all x's treated as exogenous already */
	return 0;
    }

    for (i=1; i<=dpd->xlist[0]; i++) {
	if (dpd->xlist[i] == 0) {
	    addc = 1;
	    break;
	}
    }

    if (addc && dpd->ilist != NULL) {
	for (i=1; i<=dpd->ilist[0]; i++) {
	    if (dpd->ilist[i] == 0) {
		addc = 0;
		break;
	    }
	}
    }

    if (addc) {
	dpd->ilist = gretl_list_append_term(&dpd->ilist, 0);
	if (dpd->ilist == NULL) {
	    err = E_ALLOC;
	} else {
	    dpd->nzr += 1;
	}
    }

    return err;
}

static int dpd_make_laglist (dpdinfo *dpd)
{
    int err = 0;

    /* FIXME this should be hooked up to user input */
    dpd->laglist = gretl_list_new(dpd->p);

    if (dpd->laglist == NULL) {
	err = E_ALLOC;
    } else {
	int i;

	for (i=1; i<=dpd->p; i++) {
	    dpd->laglist[i] = i;
	}
    }

    return err;
}

static int dpd_make_lists (dpdinfo *dpd, const int *list, int xpos)
{
    int i, nz = 0, spos = 0;
    int err = 0;

    /* do we have a separator, between exog vars and instruments? */
    for (i=xpos; i<=list[0]; i++) {
	if (list[i] == LISTSEP) {
	    spos = i;
	    break;
	}
    }

#if ADEBUG
    printlist(list, "incoming list in dpd_make_lists");
    fprintf(stderr, "separator pos = %d\n", spos);
#endif

    if (spos > 0) {
	dpd->nx = spos - xpos;
	nz = list[0] - spos;
    } else {
	dpd->nx = list[0] - (xpos - 1);
    }

#if ADEBUG
    fprintf(stderr, "got intial dpd->nx = %d, nz = %d\n", dpd->nx, nz);
#endif

    if (dpd->nx > 0) {
	/* compose indep vars list */
	dpd->xlist = gretl_list_new(dpd->nx);
	if (dpd->xlist == NULL) {
	    err = E_ALLOC;
	} else {
	    for (i=0; i<dpd->nx; i++) {
		dpd->xlist[i+1] = list[xpos + i];
	    }
#if ADEBUG
	    printlist(dpd->xlist, "dpd->xlist");
#endif
	    if (nz == 0 && dpd->nzb == 0) {
		/* implicitly all x vars are exogenous */
		dpd->ilist = gretl_list_copy(dpd->xlist);
		if (dpd->ilist == NULL) {
		    err = E_ALLOC;
		} else {
		    dpd->nzr = dpd->ilist[0];
		}
	    }
	}
    }

    if (!err && nz > 0) {
	/* compose regular instruments list */
	dpd->ilist = gretl_list_new(nz);
	if (dpd->ilist == NULL) {
	    err = E_ALLOC;
	} else {
	    for (i=0; i<nz; i++) {
		dpd->ilist[i+1] = list[spos + i + 1];
	    }
	    dpd->nzr = nz;
	}
    } 

    if (!err) {
	err = maybe_add_const_to_ilist(dpd);
    }

#if ADEBUG
    printlist(dpd->ilist, "dpd->ilist");
#endif

    return err;
}

/* Process the incoming list and figure out various parameters
   including the maximum lag of the dependent variable, dpd->p,
   and the number of exogenous regressors, dpd->nx.

   (If we're going to handle discontinuous lags via
   dpd->laglist, this would be the place to set it up.)

   FIXME: qmin doesn't actually do anything below, yet.
*/

static int dpd_process_list (dpdinfo *dpd, const int *list)
{
    int xpos = 0;
    int err = 0;

    dpd->p = list[1];
    dpd->qmin = 2;

    if (list[2] == LISTSEP) {
	/* arbond p ; y ... */
	dpd->qmax = 0;
	dpd->yno = list[3];
	xpos = 4;
    } else if (list[3] == LISTSEP && list[0] >= 4) {
	/* arbond p qmax ; y ... */
	dpd->qmax = list[2];
	dpd->yno = list[4];
	xpos = 5;
    } else if (list[4] == LISTSEP && list[0] >= 5) {
	/* arbond p qmax qmin ; y ... */
	dpd->qmax = list[2];
	dpd->qmin = list[3];
	dpd->yno = list[5];
	xpos = 6;
    } else {
	err = E_PARSE;
    }

    if (!err) {
	/* FIXME: are these tests all valid? */
	if (dpd->p < 1) {
	    fprintf(stderr, "arbond lag order = %d < 1\n", dpd->p);
	    err = E_INVARG;
	} else if (dpd->qmax != 0 && dpd->qmax < dpd->p + 1) {
	    fprintf(stderr, "arbond qmax = %d < dpd->p + 1 = %d\n", 
		    dpd->qmax, dpd->p + 1);
	    err = E_INVARG;
	} else if (dpd->qmax != 0 && dpd->qmin > dpd->qmax) {
	    fprintf(stderr, "arbond qmin = %d > dpd->qmax = %d\n", 
		    dpd->qmin, dpd->qmax);
	    err = E_INVARG;
	}
    }

    if (!err && list[0] >= xpos) {
	err = dpd_make_lists(dpd, list, xpos);
    }

    if (!err && (dpd->flags & DPD_NEWSTYLE)) {
	err = dpd_make_laglist(dpd);
    }

    return err;
}

static int dpd_add_unit_info (dpdinfo *dpd)
{
    int i, err = 0;

    dpd->ui = malloc(dpd->N * sizeof *dpd->ui);

    if (dpd->ui == NULL) {
	err = E_ALLOC;
    } else {
	for (i=0; i<dpd->N; i++) {
	    dpd->ui[i].skip = NULL;
	    dpd->ui[i].t1 = 0;
	    dpd->ui[i].t2 = 0;
	    dpd->ui[i].nobs = 0;
	}
    }

    return err;
}

static int dpd_flags_from_opt (gretlopt opt)
{
    /* apply Windmeijer correction unless OPT_A is given */
    int f = (opt & OPT_A)? 0 : DPD_WINCORR;

    if (opt & OPT_D) {
	/* include time dummies */
	f |= DPD_TIMEDUM;
    }

    if (opt & OPT_H) {
	/* use orthogonal deviations instead of first diffs */
	f |= DPD_ORTHDEV;
    }

    if (opt & OPT_T) {
	/* two-step estimation */
	f |= DPD_TWOSTEP;
    }

    if (opt & OPT_B) {
	/* new calculation method */
	f |= DPD_NEWSTYLE;

	if (opt & OPT_L) {
	    /* system GMM: include levels equations */
	    f |= DPD_SYSTEM;
	}

	if (opt & OPT_X) {
	    /* compute H as per Ox/DPD */
	    f |= DPD_DPDSTYLE;
	}	
    }

    return f;
}

static dpdinfo *dpdinfo_new (const int *list, const double **Z,
			     const DATAINFO *pdinfo, gretlopt opt, 
			     struct diag_info *d, int nzb,
			     int *err)
{
    dpdinfo *dpd = NULL;
    int NT;

    if (list[0] < 3) {
	*err = E_PARSE;
	return NULL;
    }

    dpd = malloc(sizeof *dpd);
    if (dpd == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    /* set pointer members to NULL just in case */
    dpd->B1 = dpd->B2 = NULL;
    dpd->V = NULL;
    dpd->ui = NULL;
    dpd->used = NULL;
    dpd->xlist = dpd->ilist = dpd->laglist = NULL;

    dpd->flags = dpd_flags_from_opt(opt);

    dpd->d = d;
    dpd->nzb = nzb;
    dpd->step = 1;
    dpd->nx = 0;
    dpd->nzr = 0;
    dpd->t1min = 0;
    dpd->ndum = 0;

    *err = dpd_process_list(dpd, list);
    if (*err) {
	goto bailout;
    }

    NT = sample_size(pdinfo);

    dpd->t1 = pdinfo->t1;             /* start of sample range */
    dpd->T = pdinfo->pd;              /* max obs. per individual */
    dpd->effN = dpd->N = NT / dpd->T; /* included individuals */
    dpd->k = dpd->p + dpd->nx;        /* # of coeffs on lagged dep var, indep vars */
    dpd->minTi = 0;
    dpd->maxTi = 0;
    dpd->max_ni = 0;
    dpd->nz = 0;
    dpd->pc0 = 0;
    dpd->xc0 = 0;
    dpd->nobs = 0;
    dpd->ndiff = 0;
    dpd->nlev = 0;
    dpd->SSR = NADBL;
    dpd->s2 = NADBL;
    dpd->AR1 = NADBL;
    dpd->AR2 = NADBL;
    dpd->sargan = NADBL;
    dpd->wald = NADBL;
    dpd->wdf = 0;

#if ADEBUG
    fprintf(stderr, "yno = %d, p = %d, qmax = %d, qmin = %d, nx = %d, k = %d\n",
	    dpd->yno, dpd->p, dpd->qmax, dpd->qmin, dpd->nx, dpd->k);
    fprintf(stderr, "t1 = %d, T = %d, N = %d\n", dpd->t1, dpd->T, dpd->N);
#endif

    *err = dpd_add_unit_info(dpd);

 bailout:

    if (*err) {
	dpdinfo_free(dpd);
	dpd = NULL;
    }

    return dpd;
}

/* See if we have valid values for the dependent variable (in
   differenced or deviations form) plus p lags of same, and all of
   the independent variables, at the given observation, s.
   TODO: either generalize this or add a parallel function
   for handling equations in levels.
*/

static int obs_is_usable (dpdinfo *dpd, const double **Z, int s)
{
    int imax = dpd->p + 1;
    int i;

    /* FIXME orthgonal deviations: we can compute such a
       deviation for time-slot s if we have a value for
       y at s-1 (given that we shift these forward) plus
       one or more valid observations at s, s+1, ...
       If I'm thinking correctly we don't necessarily 
       require a valid value at s, or in other words
       the first block below is too conservative.
    */

    for (i=0; i<=imax; i++) {
	if (na(Z[dpd->yno][s-i])) {
	    return 0;
	}
    }

    if (dpd->xlist != NULL) {
	/* check the independent vars */
	for (i=1; i<=dpd->xlist[0]; i++) {
	    if (na(Z[dpd->xlist[i]][s])) {
		return 0;
	    }
	}
    }

    return 1;
}

static int bzcols (dpdinfo *dpd, int i)
{
    int j, k, nc = 0;
    int t = i + dpd->p + 1; /* ?? */

    for (j=0; j<dpd->nzb; j++) {
	for (k=dpd->d[j].minlag; k<=dpd->d[j].maxlag; k++) {
	    if (t - k >= 0) {
		nc++;
	    }
	}
    }

    return nc;
}

/* find the first y observation, all units */

static int dpd_find_t1min (dpdinfo *dpd, const double *y)
{
    int t1min = dpd->T - 1;
    int i, s, t;

    for (i=0; i<dpd->N; i++) {
	/* loop across units */
	if (t1min > 0) {
	    s = data_index(dpd, i);
	    for (t=0; t<dpd->T; t++) {
		if (!na(y[s+t])) {
		    if (t < t1min) {
			t1min = t;
		    }
		    break;
		}
	    }
	}
    }

    return t1min;
}

/* check the sample for unit/individual i */

static int dpd_sample_check_unit (dpdinfo *dpd, const double **Z, int i,
				  int s, int *t1imin, int *t2max)
{
    unit_info *unit = &dpd->ui[i];
    int t1i = dpd->T - 1, t2i = 0; 
    int tmin = dpd->p + 1;
    int Ti = 0;
    char *mask = NULL;
    int t;

#if ADEBUG
    fprintf(stderr, "Checking unit %d: s = %d\n", i, s);
#endif

    mask = calloc(dpd->T, 1);
    if (mask == NULL) {
	return E_ALLOC;
    }

    unit->nobs = 0;

    /* For the given unit: identify the observations at which we can
       form the required delta y terms, have the requisite independent
       variables, and can construct at least one orthogonality
       condition using a lagged level of y.
    */

    for (t=tmin; t<dpd->T; t++) {
	if (obs_is_usable(dpd, Z, s + t)) {
	    unit->nobs += 1;
	    if (t < t1i) t1i = t;
	    if (t > t2i) t2i = t;
	} else {
	    mask[t] = 1;
	}
    }

    if (unit->nobs == 0) {
	dpd->effN -= 1;
	unit->t1 = -1;
	unit->t2 = -1;
	unit->nobs = 0;
#if ADEBUG
	fprintf(stderr, "unit %d not usable\n", i);
#endif
	return 0;
    }

    Ti = t2i - t1i + 1;

    if (unit->nobs < Ti) {
	/* there were gaps for this unit: save the mask */
	unit->skip = mask;
    } else {
	/* the mask is not needed */
	free(mask);
    }

    if (unit->nobs > dpd->maxTi) {
	dpd->maxTi = unit->nobs;
    }

    dpd->nobs += unit->nobs;

    if (t1i < *t1imin) {
	*t1imin = t1i;
    }	

    if (t2i > *t2max) {
	*t2max = t2i;
    }

    unit->t1 = t1i;
    unit->t2 = t2i;

#if ADEBUG
    fprintf(stderr, "t1 = %d, t2 = %d, Ti = %d, usable obs = %d\n", 
	    t1i, t2i, Ti, unit->nobs);
#endif

    return 0;
}

/* compute the column-structure of the matrix Zi */

static void dpd_compute_Z_cols (dpdinfo *dpd, int t1min, int t2max)
{
    int tau = t2max - t1min + 1;
    int nblocks = tau - dpd->p - 1;
    int i, bcols, cols = 0;

#if ADEBUG
    fprintf(stderr, "\ntau = %d (dpd->p = %d)\n", tau, dpd->p);
#endif

    for (i=0; i<nblocks; i++) {
	/* block-diagonal lagged y values */
	cols = (dpd->p + i > dpd->qmax - 1)? dpd->qmax - 1 : dpd->p + i;
#if ADEBUG
	fprintf(stderr, "block %d: adding %d cols for y lags\n", i, cols);
#endif
	dpd->nz += cols;
	if (dpd->nzb > 0) {
	    /* other block-diagonal instruments */
	    bcols = bzcols(dpd, i);
	    dpd->nz += bcols;
#if ADEBUG
	    fprintf(stderr, " plus %d cols for z lags\n", bcols);
#endif
	}
    }

#if ADEBUG
    fprintf(stderr, "'basic' m = %d\n", dpd->nz);
#endif

    dpd->qmax = cols + 1;
    /* record the column where the exogenous vars start, in xc0 */
    dpd->xc0 = dpd->nz;
    dpd->nz += dpd->nzr;
    dpd->nz += dpd->ndum;

#if ADEBUG
    fprintf(stderr, "total m = %d (dummies = %d, exog = %d)\n", 
	    dpd->nz, dpd->ndum, dpd->nzr);
#endif
}

static int dpd_sample_check (dpdinfo *dpd, const DATAINFO *pdinfo,
			     const double **Z)
{
    int t1min, t1imin = dpd->T - 1;
    int s, t2max = 0;
    int i, err = 0;

    /* find "global" first y observation */
    t1min = dpd_find_t1min(dpd, Z[dpd->yno]);

#if ADEBUG
    fprintf(stderr, "dpd_sample_check, initial scan: "
	    "dpd->T = %d, t1min = %d\n", dpd->T, t1min);
#endif

    if (dpd->qmax == 0) {
	dpd->qmax = dpd->T;
    }

    s = pdinfo->t1;

    for (i=0; i<dpd->N && !err; i++) {
	err = dpd_sample_check_unit(dpd, Z, i, s, &t1imin, &t2max);
	s += dpd->T;
    }

    if (err) {
	return err;
    }

    /* record first usable obs, any unit */
    dpd->t1min = t1imin;

    /* figure number of time dummies, if wanted */
    if (dpd->flags & DPD_TIMEDUM) {
	dpd->ndum = t2max - t1imin;
	dpd->k += dpd->ndum;
    }

    /* is t1min actually "reachable"? */
    if (t1min < t1imin - dpd->qmax) {
	t1min = t1imin - dpd->qmax;
    }

#if ADEBUG
    fprintf(stderr, "Number of units with usable observations = %d\n", dpd->effN);
    fprintf(stderr, "Total usable observations = %d\n", dpd->nobs);
    fprintf(stderr, "Max equations for any given unit = %d\n", dpd->maxTi);
    fprintf(stderr, "Maximal relevant time-series span: %d to %d = %d\n", 
	    t1min, t2max, t2max - t1min + 1);
#endif

    if (dpd->effN == 0) {
	err = E_MISSDATA;
    } else {
	/* compute the number of columns in Zi */
	dpd_compute_Z_cols(dpd, t1min, t2max);
	dpd->max_ni = dpd->maxTi;
    }

    return err;
}

/* should a certain observation be skipped? */

#define skip_obs(u,t) ((u->skip == NULL)? 0 : u->skip[t])

/* should we skip a certain panel unit altogether? */

#define skip_unit(u) (u->t1 < 0)

/* should a certain panel unit be skipped when computing
   a z-statistic for autocorrelated errors? 
*/

static int ar_skip_unit (unit_info *u, int k)
{
    int t;

    for (t=u->t1+k; t<=u->t2; t++) {
	if (!skip_obs(u, t) && !skip_obs(u, t-k)) {
	    return 0;
	}
    }

    return 1;
}

/* See if we have sufficient data to calculate the z-statistic for
   AR(k) errors: return the length of the needed arrays, or 0 if we
   can't do it.
*/

static int ar_data_check (dpdinfo *dpd, int k)
{
    unit_info *ui;
    int i, t, T = 0;

    for (i=0; i<dpd->N; i++) {
	ui = &dpd->ui[i];
	if (skip_unit(ui)) {
	    continue;
	}
	for (t=ui->t1+k; t<=ui->t2; t++) {
	    if (!skip_obs(ui, t) && !skip_obs(ui, t-k)) {
		T++;
	    }
	}
    }

    return T;
}

static int dpd_const_pos (dpdinfo *dpd)
{
    int i;

    if (dpd->xlist == NULL) {
	return 0;
    }

    for (i=1; i<=dpd->xlist[0]; i++) {
	if (dpd->xlist[i] == 0) {
	    return i;
	}
    }

    return 0;
}

/* FIXME try to detect and omit periodic dummies? */

static int dpd_wald_test (dpdinfo *dpd)
{
    gretl_matrix_block *B;
    gretl_matrix *vcv = NULL;
    gretl_vector *b = NULL;
    double x = 0.0;
    int cpos = dpd_const_pos(dpd);
    int i, j, k, kc = dpd->p + dpd->nx;
    int ri, rj;
    int err;

    /* position of const in coeff vector? */
    if (cpos != 0) {
	k = kc - 1;
	cpos += dpd->p - 1;
    } else {
	k = kc;
	cpos = -1;
    }

    B = gretl_matrix_block_new(&b,   k, 1,
			       &vcv, k, k,
			       NULL);
    if (B == NULL) {
	return E_ALLOC;
    }

    ri = 0;
    for (i=0; i<kc; i++) {
	if (i != cpos) {
	    b->val[ri++] = dpd->beta->val[i];
	}
    }

    ri = 0;
    for (i=0; i<kc; i++) {
	if (i != cpos) {
	    rj = 0;
	    for (j=0; j<kc; j++) {
		if (j != cpos) {
		    x = gretl_matrix_get(dpd->vbeta, i, j);
		    gretl_matrix_set(vcv, ri, rj++, x);
		}
	    }
	    ri++;
	}
    } 

    err = gretl_invert_symmetric_matrix(vcv);

    if (!err) {
	x = gretl_scalar_qform(b, vcv, &err);
    }

    if (!err) {
	dpd->wald = x;
	dpd->wdf = k;
    }

#if ADEBUG
    fprintf(stderr, "Wald chi^2(%d) = %g\n", k, x);
#endif

    gretl_matrix_block_destroy(B);
    
    return err;
}

static int sargan_test (dpdinfo *dpd)
{
    gretl_matrix_block *B;
    gretl_matrix *Zu = NULL;
    gretl_matrix *m1 = NULL;
    int err = 0;

    B = gretl_matrix_block_new(&Zu, dpd->nz, 1,
			       &m1, dpd->nz, 1,
			       NULL);
    if (B == NULL) {
	return E_ALLOC;
    }

    gretl_matrix_multiply(dpd->ZT, dpd->uhat, Zu);

    if ((dpd->flags & DPD_NEWSTYLE) && dpd->step == 1) {
	; /* FIXME regularize this somehow */
    } else {
	gretl_matrix_divide_by_scalar(dpd->A, dpd->effN);
    }
    gretl_matrix_multiply(dpd->A, Zu, m1);

    dpd->sargan = gretl_matrix_dot_product(Zu, GRETL_MOD_TRANSPOSE,
					   m1, GRETL_MOD_NONE,
					   &err);

    if (dpd->step == 1) {
	/* allow for scale factor in H matrix */
	if (dpd->flags & DPD_ORTHDEV) {
	    dpd->sargan /= dpd->s2;
	} else {
	    fprintf(stderr, "before scaling, sargan = %g\n", dpd->sargan);
	    dpd->sargan *= 2.0 / dpd->s2; 
	}
    }

#if ADEBUG
    fprintf(stderr, "Sargan df = m - k = %d - %d\n",
	    dpd->nz, dpd->k);
    fprintf(stderr, "Sargan test: Chi-square(%d) = %g\n",
	    dpd->nz - dpd->k, dpd->sargan);
#endif

    gretl_matrix_block_destroy(B);

    return err;
}

/* Compute the z test for AR errors, if possible.  This should
   perhaps be rolled into dpd_variance() for the sake of
   efficiency.
*/

static int ar_test (dpdinfo *dpd, const gretl_matrix *C)
{
    gretl_matrix_block *B = NULL;
    gretl_matrix *v;
    gretl_matrix *vk;
    gretl_matrix *X;  /* this is the trimmed "X_{*}" */
    gretl_matrix *vkX; 
    gretl_matrix *tmpk;
    gretl_matrix *ui; 
    gretl_matrix *m1; 
    gretl_matrix *SZv;
    double x, num, den;
    double den2, den3;
    int s, t, q, Q;
    int sbig, k = 0;
    int i, j, err = 0;

 restart:

    Q = ar_data_check(dpd, ++k);
    if (Q == 0) {
	gretl_matrix_block_destroy(B);
	return 0;
    }

#if ADEBUG 
    fprintf(stderr, "AR(%d) test: number of usable obs = %d\n", k, Q);
#endif

    if (k == 1) {
	B = gretl_matrix_block_new(&v,    Q, 1,
				   &vk,   Q, 1,
				   &X,    Q, dpd->k,
				   &vkX,  1, dpd->k, 
				   &tmpk, 1, dpd->k,
				   &ui,   dpd->maxTi, 1,
				   &m1,   dpd->nz, 1,
				   &SZv,  dpd->nz, 1,
				   NULL);
	if (B == NULL) {
	    return E_ALLOC;
	}
    } else {
	/* k == 2 */
	gretl_matrix_reuse(v,  Q, 1);
	gretl_matrix_reuse(vk, Q, 1);
	gretl_matrix_reuse(X,  Q, -1);
	gretl_matrix_reuse(m1, dpd->nz, 1);
    }

    gretl_matrix_zero(SZv);
    den = 0.0;
    q = s = sbig = 0;

    for (i=0; i<dpd->N; i++) {
	unit_info *unit = &dpd->ui[i];
	int Ti = unit->nobs;
	double den1i = 0.0;
	int si = 0;

	if (Ti == 0) {
	    continue;
	}

	if (ar_skip_unit(unit, k)) {
	    sbig += Ti;
	    s += Ti;
	    continue;
	}

	gretl_matrix_reuse(ui, Ti, -1);
	gretl_matrix_reuse(dpd->Zi, dpd->nz, Ti);

	/* extract full-length Z'_i and u_i */

	for (t=unit->t1; t<=unit->t2; t++) {
	    if (!skip_obs(unit, t)) {
		for (j=0; j<dpd->nz; j++) {
		    x = gretl_matrix_get(dpd->ZT, j, sbig);
		    gretl_matrix_set(dpd->Zi, j, si, x);
		}
		x = dpd->uhat->val[sbig];
		gretl_vector_set(ui, si, x);
		sbig++;
		si++;
	    }
	}

	for (t=unit->t1; t<unit->t1+k; t++) {
	    /* skip any obs prior to t1 + k */
	    if (!skip_obs(unit, t)) {
		s++;
	    }
	}

	/* extract lagged residuals vk along with v_{*} and X_{*} */

	for (t=unit->t1+k; t<=unit->t2; t++) {
	    if (!skip_obs(unit, t)) {
		if (!skip_obs(unit, t-k)) {
		    v->val[q] = dpd->uhat->val[s];
		    vk->val[q] = dpd->uhat->val[s-k];
		    den1i += v->val[q] * vk->val[q];
		    for (j=0; j<dpd->k; j++) {
			x = gretl_matrix_get(dpd->X, s, j);
			gretl_matrix_set(X, q, j, x);
		    }
		    q++;
		}
		s++;
	    }
	}

	/* cumulate Z'_i u_i u_i'_* u_{i,-k} */
	gretl_matrix_multiply_by_scalar(ui, den1i);
	gretl_matrix_multiply_mod(dpd->Zi, GRETL_MOD_NONE,
				  ui, GRETL_MOD_NONE,
				  SZv, GRETL_MOD_CUMULATE);

	den += den1i * den1i;
    }

#if ADEBUG > 1
    gretl_matrix_print(vk, "lagged residuals");
    gretl_matrix_print(v, "current residuals");
#endif

    /* the numerator */
    num = 0.0;
    for (i=0; i<Q; i++) {
	num += vk->val[i] * v->val[i];
    }

    /* the denominator has three components, the first of which
       is already handled */

    /* required component "vkX" = \hat{v}'_{-k} X_{*} */
    gretl_matrix_multiply_mod(vk, GRETL_MOD_TRANSPOSE,
			      X, GRETL_MOD_NONE,
			      vkX, GRETL_MOD_NONE);
    
    /* vk' X_* (X'ZAZ'X)^{-1} X'ZA(sum stuff) */
    gretl_matrix_multiply(vkX, C, tmpk);
    gretl_matrix_reuse(m1, 1, dpd->nz);
    gretl_matrix_multiply(tmpk, dpd->XZA, m1);
    den2 = gretl_matrix_dot_product(m1, GRETL_MOD_NONE,
				    SZv, GRETL_MOD_NONE,
				    &err);
    den -= 2.0 * den2;

    /* additive term vk' X_* vbeta X_*' vk */
    gretl_matrix_multiply(vkX, dpd->vbeta, tmpk);
    den3 = gretl_matrix_dot_product(tmpk, GRETL_MOD_NONE,
				    vkX, GRETL_MOD_TRANSPOSE,
				    &err);
    den += den3;

    if (den < 0) {
	err = E_NAN;
    } else {
	double z = num / sqrt(den);

#if ADEBUG
	fprintf(stderr, "AR(%d) test: z = %g / sqrt(%g) = %.4g\n", k, 
		num, den, z);
#endif

	if (k == 1) {
	    dpd->AR1 = z;
	    goto restart;
	} else {
	    dpd->AR2 = z;
	}
    }

    gretl_matrix_block_destroy(B);

    return err;
}

/* Windmeijer, Journal of Econometrics, 126 (2005), page 33:
   finite-sample correction for step-2 variance matrix.  Note: from a
   computational point of view this calculation is expressed more
   clearly in the working paper by Bond and Windmeijer, "Finite Sample
   Inference for GMM Estimators in Linear Panel Data Models" (2002),
   in particular the fact that the matrix named "D" below must be
   built column by column, taking the derivative of the "W" (or "A")
   matrix with respect to the successive independent variables.
*/

static int windmeijer_correct (dpdinfo *dpd, const gretl_matrix *uhat1,
			       const gretl_matrix *varb1)
{
    gretl_matrix_block *B;
    gretl_matrix *aV;  /* standard asymptotic variance */
    gretl_matrix *D;   /* finite-sample factor */
    gretl_matrix *dWj; /* one component of the above */
    gretl_matrix *ui;  /* per-unit residuals */
    gretl_matrix *xij; /* per-unit X_j values */
    gretl_matrix *TT;  /* workspace follows */
    gretl_matrix *mT;  
    gretl_matrix *km;  
    gretl_matrix *k1; 
    int totobs = dpd->nobs;
    int i, j, t;
    int err = 0;

    aV = gretl_matrix_copy(dpd->vbeta);
    if (aV == NULL) {
	return E_ALLOC;
    }

    if (dpd->flags & DPD_SYSTEM) {
	/* levels included */
	totobs += dpd->ndiff;
    }    

    B = gretl_matrix_block_new(&D,   dpd->k, dpd->k,
			       &dWj, dpd->nz, dpd->nz,
			       &ui,  dpd->max_ni, 1,
			       &xij, dpd->max_ni, 1,
			       &TT,  dpd->max_ni, dpd->max_ni,
			       &mT,  dpd->nz, totobs,
			       &km,  dpd->k, dpd->nz,
			       &k1,  dpd->k, 1,
			       NULL);
    if (B == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    /* form -(1/N) * asyV * XZW^{-1} */
    gretl_matrix_multiply(aV, dpd->XZA, dpd->kmtmp);
    gretl_matrix_multiply_by_scalar(dpd->kmtmp, -1.0 / dpd->effN);

    /* form W^{-1}Z'v_2 */
    gretl_matrix_multiply(dpd->A, dpd->ZT, mT);
    gretl_matrix_multiply(mT, dpd->uhat, dpd->R1);

    for (j=0; j<dpd->k; j++) { /* loop across the X's */
	int s = 0;

	gretl_matrix_zero(dWj);

	/* form dWj = -(1/N) \sum Z_i'(X_j u' + u X_j')Z_i */

	for (i=0; i<dpd->N; i++) {
	    int ni = dpd->ui[i].nobs;

	    if (ni == 0) {
		continue;
	    }

	    gretl_matrix_reuse(ui, ni, 1);
	    gretl_matrix_reuse(xij, ni, 1);
	    gretl_matrix_reuse(dpd->Zi, ni, dpd->nz);
	    gretl_matrix_reuse(TT, ni, ni);

	    /* extract ui (first-step residuals) */
	    for (t=0; t<ni; t++) {
		ui->val[t] = uhat1->val[s++];
	    }

	    /* extract xij */
	    gretl_matrix_extract_matrix(xij, dpd->X, s - ni, j,
					GRETL_MOD_NONE);
	    gretl_matrix_multiply_mod(ui, GRETL_MOD_NONE,
				      xij, GRETL_MOD_TRANSPOSE,
				      TT, GRETL_MOD_NONE);
	    gretl_matrix_add_self_transpose(TT);

	    /* extract Zi */
	    gretl_matrix_extract_matrix(dpd->Zi, dpd->ZT, 0, s - ni,
					GRETL_MOD_TRANSPOSE);

	    gretl_matrix_qform(dpd->Zi, GRETL_MOD_TRANSPOSE,
			       TT, dWj, GRETL_MOD_CUMULATE);
	}

	gretl_matrix_multiply_by_scalar(dWj, -1.0 / dpd->effN);

	/* D[.,j] = -aV * XZW^{-1} * dWj * W^{-1}Z'v_2 */
	gretl_matrix_multiply(dpd->kmtmp, dWj, km);
	gretl_matrix_multiply(km, dpd->R1, k1);

	/* write into appropriate column of D (k x k) */
	for (i=0; i<dpd->k; i++) {
	    gretl_matrix_set(D, i, j, k1->val[i]);
	}
    }

    /* add to AsyV: D * AsyV */
    gretl_matrix_multiply_mod(D, GRETL_MOD_NONE,
			      aV, GRETL_MOD_NONE,
			      dpd->vbeta, GRETL_MOD_CUMULATE);

    /* add to AsyV: AsyV * D' */
    gretl_matrix_multiply_mod(aV, GRETL_MOD_NONE,
			      D, GRETL_MOD_TRANSPOSE,
			      dpd->vbeta, GRETL_MOD_CUMULATE);

    /* add to AsyV: D * var(\hat{\beta}_1) * D' */
    gretl_matrix_qform(D, GRETL_MOD_NONE, varb1,
		       dpd->vbeta, GRETL_MOD_CUMULATE);

 bailout:

    gretl_matrix_block_destroy(B);
    gretl_matrix_free(aV);

    return err;
}

/* second-step (asymptotic) variance */

static int dpd_variance_2 (dpdinfo *dpd,
			   gretl_matrix *u1,
			   gretl_matrix *V1)
{
    int err;

    err = gretl_matrix_qform(dpd->XZ, GRETL_MOD_NONE, dpd->V,
			     dpd->vbeta, GRETL_MOD_NONE);

    if (!err) {
	err = gretl_invert_symmetric_matrix(dpd->vbeta);
    }

    if (!err) {
	gretl_matrix_multiply_by_scalar(dpd->vbeta, dpd->effN);
    }  

    if (!err && u1 != NULL && V1 != NULL) {
	err = windmeijer_correct(dpd, u1, V1);
    } 

    return err;
}

/* 
   Compute the step-1 robust variance matrix:

   N * C^{-1} * (X'*Z*A_N*\hat{V}_N*A_N*Z'*X) * C^{-1} 

   where 

   C = X'*Z*A_N*Z'*X 
   
   and

   \hat{V}_N = N^{-1} \sum Z_i'*v_i*v_i'*Z_i,

   (v_i being the step-1 residuals).

*/

static int dpd_variance_1 (dpdinfo *dpd)
{
    gretl_matrix *kk, *V, *ui;
    int i, t, k, c;
    int err = 0;

    kk = gretl_matrix_alloc(dpd->k, dpd->k);
    V = gretl_zero_matrix_new(dpd->nz, dpd->nz);
    ui = gretl_column_vector_alloc(dpd->max_ni);

    if (kk == NULL || V == NULL || ui == NULL) {
	gretl_matrix_free(kk);
	gretl_matrix_free(V);
	gretl_matrix_free(ui);
	return E_ALLOC;
    }

    c = k = 0;

    for (i=0; i<dpd->N; i++) {
	int ni = dpd->ui[i].nobs;

	if (ni == 0) {
	    continue;
	}

	/* get per-unit instruments matrix, Zi */
	gretl_matrix_reuse(dpd->Zi, ni, dpd->nz);
	gretl_matrix_reuse(ui, ni, 1);
	gretl_matrix_extract_matrix(dpd->Zi, dpd->ZT, 0, c,
				    GRETL_MOD_TRANSPOSE);
	c += ni;
	
	/* load residuals into the ui vector */
	for (t=0; t<ni; t++) {
	    ui->val[t] = dpd->uhat->val[k++];
	}

	gretl_matrix_multiply_mod(ui, GRETL_MOD_TRANSPOSE,
				  dpd->Zi, GRETL_MOD_NONE,
				  dpd->L1, GRETL_MOD_NONE);
	gretl_matrix_multiply_mod(dpd->L1, GRETL_MOD_TRANSPOSE,
				  dpd->L1, GRETL_MOD_NONE,
				  V, GRETL_MOD_CUMULATE);
    }

    gretl_matrix_divide_by_scalar(V, dpd->effN);

    /* form X'Z A_N Z'X  */
    gretl_matrix_multiply(dpd->XZ, dpd->A, dpd->kmtmp);
    gretl_matrix_qform(dpd->kmtmp, GRETL_MOD_NONE, V,
		       kk, GRETL_MOD_NONE);

    /* pre- and post-multiply by C^{-1} */
    if (!(dpd->flags & DPD_NEWSTYLE)) {
	/* the new-style \hat{\beta} calculation has 
	   already done this inversion */
	err = gretl_invert_symmetric_matrix(dpd->den);
    }
    if (!err) {
	gretl_matrix_qform(dpd->den, GRETL_MOD_NONE, kk, 
			   dpd->vbeta, GRETL_MOD_NONE);
	gretl_matrix_multiply_by_scalar(dpd->vbeta, dpd->effN);
    }

    if (!err && (dpd->flags & DPD_TWOSTEP)) {
	/* preserve V for the second stage */
	dpd->V = V;
    } else {
	gretl_matrix_free(V);
    }

    gretl_matrix_free(ui);
    gretl_matrix_free(kk);

    return err;
}

static void dpd_residuals (dpdinfo *dpd)
{
    const double *b = dpd->beta->val;
    double x, ut;
    int i, j, t, k = 0;

    dpd->SSR = 0.0;

    for (i=0; i<dpd->N; i++) {
	int Ti = dpd->ui[i].nobs;

	for (t=0; t<Ti; t++) {
	    ut = dpd->Y->val[k]; 
	    for (j=0; j<dpd->k; j++) {
		x = gretl_matrix_get(dpd->X, k, j);
		ut -= b[j] * x;
	    }
	    dpd->SSR += ut * ut;
	    dpd->uhat->val[k++] = ut;
	}
    }

    dpd->s2 = dpd->SSR / (dpd->nobs - dpd->k);
}

static int dpd_variance (dpdinfo *dpd)
{
    gretl_matrix *u1 = NULL;
    gretl_matrix *V1 = NULL;
    int err;

    if (dpd->step == 2 && (dpd->flags & DPD_WINCORR)) {
	/* we'll need a copy of the step-1 residuals, and also
	   the step-1 var(\hat{\beta})
	*/
	u1 = gretl_matrix_copy(dpd->uhat);
	V1 = gretl_matrix_copy(dpd->vbeta);
	if (u1 == NULL || V1 == NULL) {
	    return E_ALLOC;
	}
    }	

    dpd_residuals(dpd);

    if (dpd->step == 2) {
	err = dpd_variance_2(dpd, u1, V1);
	gretl_matrix_free(u1);
	gretl_matrix_free(V1);
    } else {
	err = dpd_variance_1(dpd);
    }

    if (err) {
	return err;
    }

#if ADEBUG
    int i;

    gretl_matrix_print(dpd->vbeta, "Var(beta)");
    for (i=0; i<dpd->k; i++) {
	x = gretl_matrix_get(dpd->vbeta, i, i);
	fprintf(stderr, "se(beta[%d]) = %g\n", i, sqrt(x));
    }
    fprintf(stderr, "\nSSR = %.11g\n", dpd->SSR);
    fprintf(stderr, "sigma^2 = %.7g\n", dpd->s2);
    fprintf(stderr, "sigma = %.7g\n", sqrt(dpd->s2));
#endif

    /* if we're on the final step, carry out the various tests */

    if (dpd->step == 2 || !(dpd->flags & DPD_TWOSTEP)) {
	if (!(dpd->flags & DPD_ORTHDEV)) {
	    if (dpd->step == 2) {
		err = gretl_invert_symmetric_matrix(dpd->den);
	    }
	    if (!err) {
		ar_test(dpd, dpd->den);
	    }
	}

	sargan_test(dpd);
	dpd_wald_test(dpd);
    }

    return err;
}

static int next_obs (unit_info *ui, int j0, int n)
{
    int j;

    for (j=j0; j<n; j++) {
	if (ui->skip[j + ui->t1] == 0) {
	    return j;
	}
    }

    return 0;
}

/* construct the H matrix for first-differencing
   as applied to unit i */

static int make_first_diff_matrix (dpdinfo *dpd, int i)
{
    static int *rc;
    unit_info *unit;
    int n, m;
    double x;
    int k, j, adjacent, skip = 0;

    if (dpd == NULL) {
	/* clean-up signal */
	free(rc);
	rc = NULL;
	return 0;
    }

    if (rc == NULL) {
	rc = malloc(dpd->T * sizeof *rc);
	if (rc == NULL) {
	    return E_ALLOC;
	}
    }

    unit = &dpd->ui[i];
    n = unit->t2 - unit->t1 + 1;
    m = unit->nobs;

    if (m < n) {
	skip = 1;
	j = next_obs(unit, 0, n);
	for (k=0; k<m; k++) {
	    rc[k] = j;
	    j = next_obs(unit, j+1, n);
	}
    }

    gretl_matrix_reuse(dpd->H, m, m);

    for (j=0; j<m; j++) {
	for (k=j; k<m; k++) {
	    if (skip) {
		adjacent = (abs(rc[k] - rc[j]) == 1); 
	    } else {
		adjacent = (abs(k-j) == 1);
	    }
	    x = (k==j)? 2 : (adjacent)? -1 : 0;
	    gretl_matrix_set(dpd->H, j, k, x);
	    gretl_matrix_set(dpd->H, k, j, x);
	}
    }

    return 0;
}

/* In the "dpanel" case dpd->uhat may contain both differenced
   residuals and levels residuals (stacked per unit). In that case
   trim the vector down so that it conly contains the levels residuals
   prior to transcribing the residuals in series form.
*/

static int dpanel_adjust_uhat (dpdinfo *dpd)
{
    double *tmp;
    int i, k, s, t;

    tmp = malloc(dpd->nlev * sizeof *tmp);
    if (tmp == NULL) {
	return E_ALLOC;
    }

    k = s = 0;
    for (i=0; i<dpd->N; i++) {
	/* skip residuals in differences */
	k += dpd->ui[i].t1;
	for (t=0; t<dpd->ui[i].t2; t++) {
	    tmp[s++] = dpd->uhat->val[k++];
	}
    }

    for (t=0; t<dpd->nlev; t++) {
	dpd->uhat->val[t] = tmp[t];
    }

    gretl_matrix_reuse(dpd->uhat, dpd->nlev, 1);
    free(tmp);

    return 0;
}

static int dpd_finalize_model (MODEL *pmod, dpdinfo *dpd,
			       const int *list, const char *istr,
			       const double *y, const DATAINFO *pdinfo,
			       gretlopt opt)
{
    char tmp[32];
    char prefix;
    int i, j;
    int err = 0;

    pmod->t1 = pdinfo->t1;
    pmod->t2 = pdinfo->t2;
    pmod->dfn = dpd->k;
    pmod->dfd = dpd->nobs - dpd->k;

    pmod->list = gretl_list_copy(list);
    if (pmod->list == NULL) {
	pmod->errcode = E_ALLOC;
	return pmod->errcode;
    }

    gretl_model_set_int(pmod, "yno", dpd->yno);
    gretl_model_set_int(pmod, "n_included_units", dpd->effN);

    if (dpd->flags & DPD_NEWSTYLE) {
	pmod->ci = DPANEL;
    } else {
	pmod->ci = ARBOND;
    }

    pmod->ncoeff = dpd->k;
    pmod->nobs = dpd->nobs;
    pmod->full_n = pdinfo->n;
    pmod->ess = dpd->SSR;
    if (dpd->s2 >= 0) {
	pmod->sigma = sqrt(dpd->s2);
    }

    pmod->rsq = pmod->adjrsq = NADBL;
    pmod->fstt = pmod->chisq = NADBL;
    pmod->lnL = NADBL;
  
    gretl_model_allocate_params(pmod, dpd->k);
    if (pmod->errcode) {
	return pmod->errcode;
    }

    prefix = (dpd->flags & DPD_ORTHDEV)? 'O' : 'D';

    j = 0;
    for (i=0; i<dpd->p; i++) {
	pmod->params[j][0] = '\0';
	sprintf(tmp, "%c%.10s(-%d)", prefix, pdinfo->varname[dpd->yno], i+1);
	strncat(pmod->params[j++], tmp, 15);
    }

    for (i=0; i<dpd->nx; i++) {
	strcpy(pmod->params[j++], pdinfo->varname[dpd->xlist[i+1]]);
    }

    for (i=0; i<dpd->ndum; i++) {
	sprintf(pmod->params[j++], "T%d", i + 2);
    }

    err = gretl_model_allocate_storage(pmod);

    if (!err) {
	for (i=0; i<dpd->k; i++) {
	    pmod->coeff[i] = dpd->beta->val[i];
	}
	err = gretl_model_write_vcv(pmod, dpd->vbeta);
    }

    /* add uhat, yhat */

    if (!err && dpd->nlev > 0) {
	err = dpanel_adjust_uhat(dpd);
    }

    if (!err) {
	if (dpd->flags & DPD_NEWSTYLE) {
	    int t, s = 0;

	    for (t=0; t<pdinfo->n; t++) {
		if (dpd->used[t]) {
		    pmod->uhat[t] = dpd->uhat->val[s++];
		    pmod->yhat[t] = y[t] - pmod->uhat[t];
		}
	    }
	} else {	    
	    /* old-style */
	    int s, t, k = 0;

	    for (i=0; i<dpd->N; i++) {
		unit_info *unit = &dpd->ui[i];

		if (skip_unit(unit)) {
		    continue;
		}
		for (t=0; t<dpd->T; t++) {
		    if (t >= unit->t1 && t <= unit->t2) {
			if (!skip_obs(unit, t)) {
			    s = data_index(dpd, i) + t;
			    pmod->uhat[s] = dpd->uhat->val[k];
			    pmod->yhat[s] = y[s] - pmod->uhat[s];
			    k++;
			}
		    }
		}
	    }
	} 
    }  

    /* additional dpd-specific data */

    if (!err) {
	gretl_model_set_int(pmod, "step", dpd->step);
	if (!na(dpd->AR1)) {
	    gretl_model_set_double(pmod, "AR1", dpd->AR1);
	}
	if (!na(dpd->AR2)) {
	    gretl_model_set_double(pmod, "AR2", dpd->AR2);
	}
	if (!na(dpd->sargan)) {
	    gretl_model_set_int(pmod, "sargan_df", dpd->nz - dpd->k);
	    gretl_model_set_double(pmod, "sargan", dpd->sargan);
	}
	if (!na(dpd->wald)) {
	    gretl_model_set_int(pmod, "wald_df", dpd->wdf);
	    gretl_model_set_double(pmod, "wald", dpd->wald);
	}
	if (istr != NULL && *istr != '\0') {
	    gretl_model_set_string_as_data(pmod, "istr", gretl_strdup(istr));
	}
	if (dpd->flags & DPD_TIMEDUM) {
	    pmod->opt |= OPT_D;
	}
	if ((dpd->flags & DPD_TWOSTEP) && !(dpd->flags & DPD_WINCORR)) {
	    gretl_model_set_int(pmod, "asy", 1);
	}
	if (dpd->A != NULL) {
	    gretl_model_set_int(pmod, "ninst", dpd->A->rows);
	}
	if (pmod->ci == DPANEL) {
	    if (dpd->flags & DPD_SYSTEM) {
		pmod->opt |= OPT_L;
	    }
	    if (dpd->flags & DPD_DPDSTYLE) {
		pmod->opt |= OPT_X;
	    }
	    if (dpd->maxTi > 0 && dpd->minTi > 0 && dpd->maxTi > dpd->minTi) {
		gretl_model_set_int(pmod, "Tmin", dpd->minTi);
		gretl_model_set_int(pmod, "Tmax", dpd->maxTi);
	    }
	}	    
    }

    return err;
}

/* exclude any independent variables that are zero at all
   relevant observations */

static int dpd_zero_check (dpdinfo *dpd, const double **Z)
{
    unit_info *ui;
    const double *x;
    int drop = 0;
    int i, j, k, t;

    for (j=0; j<dpd->nx; j++) {
	int all0 = 1;

	x = Z[dpd->xlist[j+1]];
	for (i=0; i<dpd->N && all0; i++) {
	    ui = &dpd->ui[i];
	    if (skip_unit(ui)) {
		continue;
	    }
	    for (t=ui->t1; t<=ui->t2 && all0; t++) {
		k = data_index(dpd, i) + t;
		if (x[k] != 0.0 && !na(x[k])) {
		    all0 = 0;
		}
	    }
	}
	if (all0) {
	    dpd->xlist[j+1] = -1;
	    drop++;
	}
    }

    if (drop > 0) {
	for (i=1; i<=dpd->xlist[0]; i++) {
	    if (dpd->xlist[i] < 0) {
		gretl_list_delete_at_pos(dpd->xlist, i);
		i--;
	    }
	}
	dpd->nx = dpd->xlist[0];
	dpd->k -= drop;
    }

    return 0;
}

/* Based on reduction of the A matrix, trim ZT to match and
   adjust the sizes of all workspace matrices that have a 
   dimension involving dpd->nz.
*/

static void real_shrink_matrices (dpdinfo *dpd, const char *mask)
{
    fprintf(stderr, "A matrix: shrinking m from %d to %d\n", 
	    dpd->nz, dpd->A->rows);

    gretl_matrix_cut_rows(dpd->ZT, mask);

    dpd->nz = dpd->A->rows;
    gretl_matrix_reuse(dpd->Acpy,  dpd->nz, dpd->nz);
    gretl_matrix_reuse(dpd->tmp1,  dpd->nz, dpd->nz);
    gretl_matrix_reuse(dpd->kmtmp, -1, dpd->nz);
    gretl_matrix_reuse(dpd->L1,    -1, dpd->nz);
    gretl_matrix_reuse(dpd->XZA,   -1, dpd->nz);
    gretl_matrix_reuse(dpd->XZ,    -1, dpd->nz);
    gretl_matrix_reuse(dpd->R1,    dpd->nz, -1);
}

/* Remove zero rows/cols from A, as indicated by mask, and delete the
   corresponding columns from Z.  At this point we leave open the
   question of whether the reduced A matrix is positive definite.
*/

static int reduce_Z_and_A (dpdinfo *dpd, const char *mask)
{
    int err = gretl_matrix_cut_rows_cols(dpd->A, mask);

    if (!err) {
	real_shrink_matrices(dpd, mask);
    }

    return err;
}

/* We already removed any zero rows/columns from A, but it couldn't
   be inverted.  Now we try reducing the dimension of A to its rank.
   We need to read from the backup, Acpy, since dpd->A will have been
   mangled in the failed inversion attempt.  We proceed to check that
   the reduced A is invertible; if not we flag an error.
*/

static int try_alt_inverse (dpdinfo *dpd)
{
    char *mask = NULL;
    int err = 0;

    gretl_matrix_copy_values(dpd->A, dpd->Acpy); 

    mask = gretl_matrix_rank_mask(dpd->A, &err);

    if (!err) {
	err = gretl_matrix_cut_rows_cols(dpd->A, mask);
    }

    if (!err) {
	err = gretl_invert_symmetric_matrix(dpd->A);
	if (!err) {
	    real_shrink_matrices(dpd, mask);
	} else {
	    fprintf(stderr, "try_alt_inverse: error inverting\n");
	}
    }

    free(mask);

    return err;
}

static int dpd_calculate (dpdinfo *dpd)
{
    int err = 0;

    /* find Z'X, if this job is not already done */
    if (dpd->step == 1) {
	err = gretl_matrix_multiply_mod(dpd->X, GRETL_MOD_TRANSPOSE,
					dpd->ZT, GRETL_MOD_TRANSPOSE,
					dpd->XZ, GRETL_MOD_NONE);
    }

#if ADEBUG
    gretl_matrix_print(dpd->XZ, "XZ' (dpd_calculate)");
    gretl_matrix_print(dpd->A, "A (dpd_calculate)");
#endif

    /* find X'Z A */
    gretl_matrix_multiply(dpd->XZ, dpd->A, dpd->XZA);

    /* calculate "numerator", X'ZAZ'y */
    gretl_matrix_multiply(dpd->ZT, dpd->Y, dpd->R1);

#if ADEBUG
    gretl_matrix_print(dpd->R1, "Z'y (dpd_calculate)");
#endif
    gretl_matrix_multiply(dpd->XZA, dpd->R1, dpd->beta);

    /* calculate "denominator", X'ZAZ'X */
    gretl_matrix_multiply_mod(dpd->XZA, GRETL_MOD_NONE,
			      dpd->XZ, GRETL_MOD_TRANSPOSE,
			      dpd->den, GRETL_MOD_NONE);
    gretl_matrix_xtr_symmetric(dpd->den);

    gretl_matrix_copy_values(dpd->kktmp, dpd->den);
    err = gretl_cholesky_decomp_solve(dpd->kktmp, dpd->beta);

#if ADEBUG
    if (!err) {
	int i;

	fprintf(stderr, "%d-step estimates:\n\n", dpd->step);
	for (i=0; i<dpd->k; i++) {
	    fprintf(stderr, "beta[%d] = %g\n", i, dpd->beta->val[i]);
	}
	fputc('\n', stderr);
    }
#endif

    if (!err) {
	err = dpd_variance(dpd);
    }

    return err;
}

static int dpd_step_2 (dpdinfo *dpd)
{
    int err = 0;

#if ADEBUG
    gretl_matrix_print(dpd->V, "V, in dpd_step_2");
#endif

    if (gretl_matrix_rows(dpd->V) > dpd->effN) {
	/* we know this case requires special attention */
	err = gretl_SVD_invert_matrix(dpd->V);
	if (!err) {
	    gretl_matrix_xtr_symmetric(dpd->V);
	}	
    } else {
	gretl_matrix *Vcpy;

	if (dpd->Acpy != NULL) {
	    Vcpy = dpd->Acpy;
	    gretl_matrix_copy_values(Vcpy, dpd->V);
	} else {
	    Vcpy = gretl_matrix_copy(dpd->V);
	    if (Vcpy == NULL) {
		return E_ALLOC;
	    }
	}	

 	err = gretl_invert_symmetric_matrix(dpd->V);
	
	if (err) {
	    /* revert the data in dpd->V */
	    gretl_matrix_copy_values(dpd->V, Vcpy);
	    err = gretl_SVD_invert_matrix(dpd->V);
	    if (!err) {
		gretl_matrix_xtr_symmetric(dpd->V);
	    }
	}
    
	if (Vcpy != dpd->Acpy) {
	    gretl_matrix_free(Vcpy);
	}
    }

    if (!err) {
	/* A <- V^{-1} */
	gretl_matrix_copy_values(dpd->A, dpd->V);
	dpd->step = 2;
	err = dpd_calculate(dpd);
    }

    if (err) {
	fprintf(stderr, "step 2: dpd_calculate returning %d\n", err);
    }

    return err;
}

static double odev_at_lag (const double *x, int t, int lag, int pd)
{
    double ret, xbar = 0.0;
    int s, Tt, n = 0;

    t -= lag + 1;

    if (t < 0 || na(x[t])) {
	return NADBL;
    }

    Tt = pd - (t % pd) - (lag + 1);

    for (s=1; s<=Tt; s++) {
	if (!na(x[t+s]) && !na(x[t+s+lag])) {
	    xbar += x[t+s];
	    n++;
	}
    }

    if (n > 0) {
	xbar /= n;
	ret = sqrt(n / (n + 1.0)) * (x[t] - xbar);
    } else {
	ret = NADBL;
    }

    return ret;
}

static int dpd_make_y_X (dpdinfo *dpd, const double **Z, 
			 const DATAINFO *pdinfo)
{
    const double *y = Z[dpd->yno];
    unit_info *unit;
    int odev = (dpd->flags & DPD_ORTHDEV);
    int i, j, s, t, k = 0;
    double x;

    for (i=0; i<dpd->N; i++) {
	unit = &dpd->ui[i];
	if (skip_unit(unit)) {
	    continue;
	}
	for (t=unit->t1; t<=unit->t2; t++) {
	    if (skip_obs(unit, t)) {
		continue;
	    }
	    s = data_index(dpd, i) + t;
	    /* current difference (or deviation) of dependent var */
	    if (odev) {
		x = odev_at_lag(y, s, 0, pdinfo->pd);
		gretl_vector_set(dpd->Y, k, x);
	    } else {
		gretl_vector_set(dpd->Y, k, y[s] - y[s-1]);
	    }
	    for (j=0; j<dpd->p; j++) {
		/* lagged difference of dependent var */
		if (odev) {
		    x = odev_at_lag(y, s, j + 1, pdinfo->pd);
		    gretl_matrix_set(dpd->X, k, j, x);
		} else {
		    gretl_matrix_set(dpd->X, k, j, y[s-j-1] - y[s-j-2]);
		}
	    }
	    for (j=0; j<dpd->nx; j++) {
		/* independent vars */
		x = Z[dpd->xlist[j+1]][s];
		gretl_matrix_set(dpd->X, k, j + dpd->p, x);
	    }
	    for (j=0; j<dpd->ndum; j++) {
		/* time dummies */
		x = (t - dpd->t1min - 1 == j)? 1 : 0;
		gretl_matrix_set(dpd->X, k, j + dpd->p + dpd->nx, x);
	    }	    
	    k++;
	}
    }

#if ADEBUG
    gretl_matrix_print(dpd->Y, "Y (arbond)");
    gretl_matrix_print(dpd->X, "X (arbond)");
#endif

    return 0;
}

static int dpd_make_Z_and_A (dpdinfo *dpd, const double **Z)
{
    const double *y = Z[dpd->yno];
    int i, j, k, s, t, c = 0;
    int zi, zj, zk;
    double x;
    char *zmask;
#if ADEBUG
    char zstr[8];
#endif
    int err = 0;

    gretl_matrix_zero(dpd->A);
    gretl_matrix_zero(dpd->XZ);
    gretl_matrix_zero(dpd->R1);

    for (i=0; i<dpd->N && !err; i++) {
	unit_info *unit = &dpd->ui[i];
	int ycols = dpd->p;   /* intial y block width */
	int offj = 0;         /* initialize column offset */
	int Ti = unit->nobs;

	if (Ti == 0) {
	    continue;
	}

	gretl_matrix_reuse(dpd->Zi, Ti, dpd->nz);
	gretl_matrix_zero(dpd->Zi);

	for (t=dpd->p+1; t<unit->t1; t++) {
	    /* unbalanced case: compute initial column offset */
	    offj += ycols;
	    if (ycols < dpd->qmax - 1) {
		ycols++;
	    }
	    zj = 0;
	    for (zi=0; zi<dpd->nzb; zi++) {
		for (zk=dpd->d[zi].minlag; zk<=dpd->d[zi].maxlag; zk++) {
		    if (t - zk >= 0) {
			zj++;
		    }
		}
	    }
	    offj += zj;
	}	    

	k = 0;
	for (t=unit->t1; t<=unit->t2; t++) {
	    int skip = skip_obs(unit, t);
	    int offincr = ycols;

	    if (!skip) {
		/* lagged y (GMM instr) columns */
		for (j=0; j<ycols; j++) {
		    s = data_index(dpd, i) + t - (ycols + 1) + j;
		    if (!na(y[s])) {
			gretl_matrix_set(dpd->Zi, k, j + offj, y[s]);
		    }
		}
	    }
	    
	    /* additional block-diagonal columns, if required --
	       needs checking for the unbalanced case 
	    */
	    zj = 0;
	    for (zi=0; zi<dpd->nzb; zi++) {
		for (zk=dpd->d[zi].minlag; zk<=dpd->d[zi].maxlag; zk++) {
		    if (t - zk >= 0) { /* ?? */
			if (!skip) {
			    s = data_index(dpd, i) + t - zk;
			    x = Z[dpd->d[zi].v][s];
			    if (!na(x)) {
				gretl_matrix_set(dpd->Zi, k, zj + offj + ycols, x);
			    }
			}
			zj++;
		    }
		}
	    }
	    offincr += zj;

	    if (!skip) {
		/* additional full-length instrument columns */
		s = data_index(dpd, i) + t;
		for (j=0; j<dpd->nzr; j++) {
		    x = Z[dpd->ilist[j+1]][s];
		    gretl_matrix_set(dpd->Zi, k, dpd->xc0 + j, x);
		}
		/* plus time dummies, if wanted */
		for (j=0; j<dpd->ndum; j++) {
		    x = (t - dpd->t1min - 1 == j)? 1 : 0;
		    gretl_matrix_set(dpd->Zi, k, dpd->xc0 + dpd->nzr + j, x);
		}
		k++; /* increment target row */	
	    }

	    offj += offincr; /* starting column for next block */
	    if (ycols < dpd->qmax - 1) {
		/* increment y block width, if we're not already at max */
		ycols++;
	    }
	}

#if ADEBUG
	sprintf(zstr, "Z_%d", i + 1);
	gretl_matrix_print(dpd->Zi, zstr);
#endif

	/* Cumulate Z_i' H Z_i into A_N */
	if (dpd->flags & DPD_ORTHDEV) {
	    /* orthogonal deviations: "H" is identity matrix */
	    gretl_matrix_multiply_mod(dpd->Zi, GRETL_MOD_TRANSPOSE,
				      dpd->Zi, GRETL_MOD_NONE,
				      dpd->A, GRETL_MOD_CUMULATE);
	} else {
	    err = make_first_diff_matrix(dpd, i);
	    gretl_matrix_qform(dpd->Zi, GRETL_MOD_TRANSPOSE,
			       dpd->H, dpd->A, GRETL_MOD_CUMULATE);
	}

	/* Write Zi into ZT at offset 0, c */
	gretl_matrix_inscribe_matrix(dpd->ZT, dpd->Zi, 0, c, GRETL_MOD_TRANSPOSE);
	c += Ti;
    }

    if (!(dpd->flags & DPD_ORTHDEV)) {
	/* clean up */
	make_first_diff_matrix(NULL, 0);
    }

    if (!err) {
	/* mask zero rows of ZT, if required */
	zmask = gretl_matrix_zero_row_mask(dpd->ZT, &err);
	if (zmask != NULL) {
	    err = reduce_Z_and_A(dpd, zmask);
	    free(zmask);
	}
    } 

#if ADEBUG
    gretl_matrix_print(dpd->A, "\\sum Z_i' H Z_i");
#endif

    if (!err) {
	gretl_matrix_divide_by_scalar(dpd->A, dpd->N);
    }

#if ADEBUG
    gretl_matrix_print(dpd->ZT, "ZT");
    gretl_matrix_print(dpd->A, "N^{-1} * \\sum Z_i' H Z_i");
#endif

    return err;
}

static int parse_diag_info (const char *s, struct diag_info *d,
			    const DATAINFO *pdinfo)
{
    char vname[VNAMELEN];
    int v, m1, m2;
    int err = 0;

    if (s == NULL) {
	err = E_ALLOC;
    } else if (sscanf(s, "GMM(%15[^, ],%d,%d)", vname, &m1, &m2) != 3) {
	err = E_PARSE;
    } else {
	if (m2 == 0) {
	    m2 = 99;
	}
	v = series_index(pdinfo, vname);
	if (v == pdinfo->v) {
	    err = E_UNKVAR;
	} else if (m1 < 0 || m2 < m1) {
	    err = E_DATA;
	} else {
	    d->v = v;
	    d->minlag = m1;
	    d->maxlag = m2;
	}
    }

    return err;
}

/* parse requests of the form 

      GMM(xvar, minlag, maxlag)

   which call for inclusion of lags of xvar in block-diagonal
   fashion
*/

static int dpd_parse_istr (const char *istr, const DATAINFO *pdinfo,
			   struct diag_info **pd, int *ns)
{
    struct diag_info *d = NULL;
    char *s0 = gretl_strdup(istr);
    char *s, *p, *spec;
    int i, err = 0;

    if (s0 == NULL) {
	return E_ALLOC;
    }

    /* count ')'-terminated fields */
    s = s0;
    while (*s) {
	if (*s == ')') {
	    *ns += 1;
	}
	s++;
    }

    if (*ns == 0) {
	err = E_PARSE;
    } else {
	d = malloc(*ns * sizeof *d);
	if (d == NULL) {
	    err = E_ALLOC;
	}
    }

    /* parse and record individual GMM instrument specs */
    s = s0;
    i = 0;
    while (*s && !err) {
	while (*s == ' ') s++;
	p = s;
	while (*p && *p != ')') p++;
	spec = gretl_strndup(s, p - s + 1);
	err = parse_diag_info(spec, &d[i++], pdinfo);
	free(spec);
	s = p + 1;
    }

    free(s0);

    if (err) {
	free(d);
	*ns = 0;
    } else {
	*pd = d;
    }

    return err;
}

static int dpd_invert_A_N (dpdinfo *dpd)
{
    int err = 0;

    /* make a backup in case the first attempt fails */
    gretl_matrix_copy_values(dpd->Acpy, dpd->A);

    /* first try straight inversion */
    err = gretl_invert_symmetric_matrix(dpd->A);

    if (err) {
	/* try again, reducing A based on its rank */
	fprintf(stderr, "inverting dpd->A failed on first pass\n");
	err = try_alt_inverse(dpd);
    }

#if ADEBUG
    gretl_matrix_print(dpd->A, "A_N");
#endif

    return err;
}

/* public interface: driver for Arellano-Bond type estimation */

MODEL
arbond_estimate (const int *list, const char *istr, const double **Z, 
		 const DATAINFO *pdinfo, gretlopt opt,
		 PRN *prn)
{
    struct diag_info *d = NULL;
    dpdinfo *dpd = NULL;
    int nzb = 0;
    MODEL mod;
    int err = 0;

    gretl_model_init(&mod);
    gretl_model_smpl_init(&mod, pdinfo);

    /* parse special instrument info, if present */
    if (istr != NULL && *istr != '\0') {
	mod.errcode = dpd_parse_istr(istr, pdinfo, &d, &nzb);
	if (mod.errcode) {
	    fprintf(stderr, "Error %d in dpd_parse_istr\n", mod.errcode);
	    return mod;
	}
    }

    /* initialize (including some memory allocation) */
    dpd = dpdinfo_new(list, Z, pdinfo, opt, d, nzb, &mod.errcode);
    if (mod.errcode) {
	fprintf(stderr, "Error %d in dpd_init\n", mod.errcode);
	return mod;
    }

    /* see if we have a usable sample */
    err = dpd_sample_check(dpd, pdinfo, Z);
    if (err) {
	fprintf(stderr, "Error %d in dpd_sample_check\n", err);
    }

    if (!err && dpd->nx > 0) {
	/* cut out any all-zero variables */
	dpd_zero_check(dpd, Z);
    }

    if (!err) {
	/* main workspace allocation */
	err = dpd_allocate_matrices(dpd);
    }

    if (!err) {
	/* build y^* and X^* */
	err = dpd_make_y_X(dpd, Z, pdinfo);
    }

    if (!err) {
	/* build instrument matrix blocks, Z_i, and insert into
	   big Z' matrix; cumulate first-stage A_N as we go 
	*/
	err = dpd_make_Z_and_A(dpd, Z);
    }

    if (!err) {
	/* invert A_N: note that we allow two attempts */
	err = dpd_invert_A_N(dpd);
    }

    if (!err) {
	/* first-step calculation */
	err = dpd_calculate(dpd);
    }

    if (!err && (opt & OPT_T)) {
	/* second step, if wanted */
	err = dpd_step_2(dpd);
    }

    if (err && !mod.errcode) {
	mod.errcode = err;
    }

    if (!mod.errcode) {
	/* write estimation info into model struct */
	mod.errcode = dpd_finalize_model(&mod, dpd, list, istr, 
					 Z[dpd->yno], pdinfo, opt);
    }

    dpdinfo_free(dpd);

    return mod;
}

#include "dpanel.c"


