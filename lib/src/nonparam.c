/*
 *  Copyright (c) by Ramu Ramanathan and Allin Cottrell
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
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

/* non-parametric stats for gretl */

#include "libgretl.h"
#include "internal.h"

static int inverse_compare_doubles (const void *a, const void *b)
{
    const double *da = (const double *) a;
    const double *db = (const double *) b;
     
    return (*da < *db) - (*da > *db);
}

/**
 * spearman:
 * @list: list of (two) variables to process.
 * @Z: data matrix.
 * @pdinfo: information on the data set.
 * @opt: if non-zero, print both the "raw" and the ranked data.
 * @prn: gretl printing struct.
 *
 * Calculates and prints Spearman's rank correlation coefficient for the two
 * variables specified in the @list.
 * 
 * Returns: 0 on successful completion, 1 on error.
 */

int spearman (const LIST list, double **Z, const DATAINFO *pdinfo, 
	      const int opt, PRN *prn)
{
    double xx, yy, *sx, *sy, *rx, *ry, *tmp;
    double xdate, rsum, avg, z = 0;
    int i, j, vx, vy, t, t1 = pdinfo->t1, t2 = pdinfo->t2;
    int rcount;
    size_t nn;

    if (list[0] != 2) {
	pprintf(prn, _("spearman command requires two variables\n"));
	return 1;
    }

    sx = malloc((t2 - t1 + 1) * sizeof *sx);
    sy = malloc((t2 - t1 + 1) * sizeof *sy);
    rx = malloc((t2 - t1 + 1) * sizeof *rx);
    ry = malloc((t2 - t1 + 1) * sizeof *ry);
    tmp = malloc((t2 - t1 + 1) * sizeof *tmp);
    if (sx == NULL || sy == NULL || rx == NULL || ry == NULL 
	|| tmp == NULL) return E_ALLOC;

    /* get sorted versions of x, y */
    vx = list[1];
    vy = list[2];
    i = -1;
    for (t=t1; t<=t2; t++) {
	xx = Z[vx][t];
	yy = Z[vy][t];
	if (na(xx) || na(yy)) continue;
	i++;
	sx[i] = xx;
	sy[i] = yy;
    }
    nn = i + 1;
    qsort(sx, nn, sizeof *sx, inverse_compare_doubles);
    qsort(sy, nn, sizeof *sy, inverse_compare_doubles);

    /* make rankings by comparing "raw" with sorted */
    i = -1;
    for (t=t1; t<=t2; t++) {
	xx = Z[vx][t];
	yy = Z[vy][t];
	if (na(xx) || na(yy)) continue;
	i++;
	for (j=0; j<nn; j++) {
	    if (floateq(xx, sx[j])) {
		rx[i] = (double) j + 1;
		break;
	    }
	}
	for (j=0; j<nn; j++) {
	    if (floateq(yy, sy[j])) {
		ry[i] = (double) j + 1;
		break;
	    }
	}	
    }
    /* fix up duplicated ranks */
    for (t=0; t<nn; t++) tmp[t] = rx[t];
    qsort(tmp, nn, sizeof *tmp, _compare_doubles);
    for (t=0; t<nn; ) {
	rsum = xx = tmp[t];
	rcount = 1;
	for (j=t+1; j<nn; j++) {
	    if (floateq(tmp[j], xx)) {
		rsum += (xx + j - t);
		rcount++;
	    } else break;
	}
	t += rcount;
	if (rcount > 1) {
	    avg = rsum / rcount;
	    for (i=0; i<nn; i++) 
		if (floateq(rx[i], xx)) rx[i] = avg;
	}
    }

    for (t=0; t<nn; t++) tmp[t] = ry[t];
    qsort(tmp, nn, sizeof *tmp, _compare_doubles);
    for (t=0; t<nn; ) {
	rsum = yy = tmp[t];
	rcount = 1;
	for (j=t+1; j<nn; j++) {
	    if (floateq(tmp[j], yy)) {
		rsum += (yy + j - t);
		rcount++;
	    } else break;
	}
	t += rcount;
	if (rcount > 1) {
	    avg = rsum / rcount;
	    for (i=0; i<nn; i++) 
		if (floateq(ry[i], yy)) ry[i] = avg;
	}
    }	

    /* calculate and print rho */
    xx = 0.;
    for (i=0; i<nn; i++) 
	xx += (rx[i] - ry[i]) * (rx[i] - ry[i]);
    yy = 1.0 - 6 * xx / (nn * (nn * nn - 1));
    xx = sqrt(1./(nn - 1));
    pprintf(prn, _("\nFor the variables '%s' and '%s'\n"), pdinfo->varname[vx],
	    pdinfo->varname[vy]);
    pprintf(prn, _("Spearman's rank correlation coefficient (rho) = %f\n"), yy);
    pprintf(prn, _("Under the null hypothesis of no correlation, rho "
	    "follows N(0, %f)\n"), xx);
    if (nn >= 20) {
	z = fabs(yy/xx);
	pprintf(prn, _("z-score = %f, with one-tailed p-value %f\n"), z, 
		normal(z));
    } else {
	pprintf(prn, _("Sample is too small to calculate a p-value based on "
		"the normal distribution\n"));
    }

    if (opt) { /* print raw and ranked data */
	if (pdinfo->pd == 1) pprintf(prn, "\n Obs ");
	else pprintf(prn, "\n\n     Obs ");
	pprintf(prn, "%13s%13s%13s%13s\n\n", pdinfo->varname[vx], "rank",
	       pdinfo->varname[vy], "rank");
	i = 0;
	for (t=t1; t<=t2; t++)   {
	    if (pdinfo->markers) { 
		pprintf(prn, "%8s ", pdinfo->S[t]); 
	    } else {
		xdate = date(t, pdinfo->pd, pdinfo->sd0);
		if (dataset_is_daily(pdinfo)) {
		    char datestr[9];
		    
		    ntodate(datestr, t, pdinfo);
		    pprintf(prn, "%8s ", datestr);
		}
		else if (pdinfo->pd == 1) 
		    pprintf(prn, "%4d ", (int) xdate);
		else if (pdinfo->pd < 10) 
		    pprintf(prn, "%8.1f ", xdate);
		else 
		    pprintf(prn, "%8.2f ", xdate);
	    }
	    xx = Z[vx][t];
	    yy = Z[vy][t];
	    if (!(na(xx)) && !(na(yy))) {
		_printxs(xx, 15, PRINT, prn);
		_printxs(rx[i], 15, PRINT, prn);
		_printxs(yy, 15, PRINT, prn);
		_printxs(ry[i], 15, PRINT, prn);
		i++;
	    }
	    pprintf(prn, "\n");
	}
    }
    free(sx);
    free(sy);
    free(rx);
    free(ry);
    free(tmp);
    return 0;
}

/**
 * runs_test:
 * @varno: ID number of the variable to process.
 * @Z: data matrix.
 * @pdinfo: information on the data set.
 * @prn: gretl printing struct.
 *
 * Performs, and prints the results of, the runs test for randomness
 * for the variable specified by @varno.
 * 
 * Returns: 0 on successful completion, non-zero on error.
 */

int runs_test (const int varno, double **Z, const DATAINFO *pdinfo, 
	       PRN *prn)
{
    int t, t1 = pdinfo->t1, t2 = pdinfo->t2, n = pdinfo->n, runs = 1;
    int nn;
    double xx, *x, mean, sd, z;

    nn = t2 - t1 + 1;
    x = malloc(nn * sizeof *x);
    if (x == NULL) return E_ALLOC;

    nn = 0;
    for (t=t1; t<=t2; t++) {
	xx = Z[varno][t];
	if (na(xx)) continue;
	else x[nn++] = xx;
    }
    if (nn <= 1) {
	pprintf(prn, _("\nInsufficient data for runs test\n"));
	free(x);
	return 1;
    }
    for (t=1; t<nn; t++) {
	if ((x[t] > 0 && x[t-1] < 0) || (x[t] < 0 && x[t-1] > 0)) 
	    runs++;
    }
    mean = (1 + nn/2.0);
    sd = sqrt((double) n - 1)/2.0;
    z = fabs((runs - mean)/sd);
    pprintf(prn, _("\nNumber of runs (R) in the variable '%s' = %d\n"), 
	    pdinfo->varname[varno], runs);
    pprintf(prn, _("Under the null hypothesis of randomness, R "
	    "follows N(%f, %f)\n"), mean, sd);
    pprintf(prn, _("z-score = %f, with two-tailed p-value %f\n"), z, 
	    2 * normal(z));    
    free(x);
    return 0;
}



