/*
 *  Copyright (c) 2003 by Allin Cottrell
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

#include "libgretl.h"
#include "gretl_matrix.h"

/* #define PCA_DEBUG 1 */

static double *standardize (const double *x, int n)
{
    double *sx;
    double xbar, sd;
    int i, err;

    err = moments(0, n-1, x, &xbar, &sd, NULL, NULL, 1);
    if (err) return NULL;

    sx = malloc(n * sizeof *sx);
    if (sx == NULL) return NULL;

    for (i=0; i<n; i++) {
	if (na(x[i])) {
	    sx[i] = NADBL;
	} else {
	    sx[i] = (x[i] - xbar) / sd;
	}
    }

    return sx;
}

int pca_from_corrmat (CORRMAT *corrmat, double ***pZ,
		      DATAINFO *pdinfo, unsigned char oflag,
		      PRN *prn)
{
    gretl_matrix *m;
    double x, y;
    int i, j, n = corrmat->list[0];
    int idx;
    double *evals;

    m = gretl_matrix_alloc(n, n);
    if (m == NULL) return E_ALLOC;

    for (i=0; i<n; i++) {
	for (j=0; j<n; j++) {
	    idx = ijton(i+1, j+1, n);
	    x = corrmat->xpx[idx];
	    gretl_matrix_set(m, i, j, x);
	}
    }

    evals = gretl_symmetric_matrix_eigenvals(m, 1);
    if (evals == NULL) {
	gretl_matrix_free(m);
	return 1;
    }

    pputs(prn, "Principal Components Analysis\n\n");
    pputs(prn, "Eigenanalysis of the Correlation Matrix\n\n");

    pputs(prn, "Component  Eigenvalue  Proportion   Cumulative\n");

    y = 0.0;
    for (i=n-1; i>=0; i--) {
	y += evals[i] / n;
	pprintf(prn, "    %d%13.4f%13.4f%13.4f\n", n - i,
		evals[i], evals[i] / n, y);
	x += evals[i];
    }
    pputs(prn, "\n");

#ifdef PCA_DEBUG
    fprintf(stderr, "check: sum of evals = %g\n", x);
#endif

    pputs(prn, "Eigenvectors (component loadings)\n\n");

    pputs(prn, "Variable  ");
    for (i=1; i<=n; i++) {
	pprintf(prn, "%8s%d", "PC", i);
    }
    pputs(prn, "\n");
    for (i=0; i<n; i++) {
	pprintf(prn, "%-10s", pdinfo->varname[corrmat->list[i+1]]);
	for (j=n-1; j>=0; j--) {
	    pprintf(prn, "%9.3f", gretl_matrix_get(m, i, j));
	}
	pputs(prn, "\n");
    }
    pputs(prn, "\n");

    if (oflag) {
	/* add components with eigenvalues > 1 to the dataset */
	int v = pdinfo->v;
	int nc = 0, err = 0;
	double **sZ = NULL;
	int add_all = (oflag == 'a');
	int *list;

	if (add_all) {
	    nc = n;
	} else {
	    for (i=0; i<n; i++) {
		if (evals[i] > 1.0) nc++;
	    }
	}

	list = malloc((nc + 1) * sizeof *list);
	if (list == NULL) err = E_ALLOC;

	if (!err) {
	    /* build list of PCs (with eigenvals > 1?) */
	    list[0] = nc;
	    j = 1;
	    for (i=n-1; i>=0; i--) {
		if (add_all || evals[i] > 1.0) list[j++] = i;
	    }
#ifdef PCA_DEBUG
	    printlist(list, "pclist");
#endif
	    err = dataset_add_vars(nc, pZ, pdinfo);
	}

	if (!err) {
	    /* construct standardized versions of variables */
	    sZ = malloc(n * sizeof *sZ);
	    if (sZ == NULL) err = E_ALLOC;
	    else {
		for (i=0; i<n; i++) sZ[i] = NULL;
		for (i=0; i<n; i++) {
		    int oldv = corrmat->list[i+1];

#ifdef PCA_DEBUG
		    fprintf(stderr, "Getting standardized version of "
			    "var %d\n", oldv);
#endif
		    sZ[i] = standardize((const double *) (*pZ)[oldv], 
					pdinfo->n);
		    if (sZ[i] == NULL) {
			err = E_ALLOC;
			break;
		    }
		}
		if (err) {
		    for (i=0; i<n; i++) free(sZ[i]);
		    free(sZ);
		    sZ = NULL;
		}
	    }
	}

	if (!err) {
	    for (i=1; i<=list[0]; i++) {
		int newv = v + i - 1;
		int pcnum = list[i];
		int t;

		sprintf(pdinfo->varname[newv], "PC%d", i);
		sprintf(VARLABEL(pdinfo, newv), "Component with "
			"eigenvalue = %.4f", evals[pcnum]);
		for (t=0; t<pdinfo->n; t++) {
#ifdef PCA_DEBUG
		    fprintf(stderr, "Obs %d\n", t);
#endif
		    (*pZ)[newv][t] = 0.0;
		    for (j=0; j<n; j++) {
			double load = gretl_matrix_get(m, j, pcnum);
			double val = sZ[j][t];

#ifdef PCA_DEBUG
			fprintf(stderr, "j=%d,pcnum=%d,load=%g,val=%g\n",
				j,pcnum,load,val);
#endif

			if (na(val)) {
			    (*pZ)[newv][t] = NADBL;
			    break;
			} else {
			    (*pZ)[newv][t] += load * val;
			}
		    }
		} /* end loop over observations */
	    } /* end loop over components */
	} /* end !err conditional */

	free(list);
	if (sZ != NULL) {
	    for (i=0; i<n; i++) free(sZ[i]);
	    free(sZ);
	}

    } /* end oflag conditional */

    free(evals);
    gretl_matrix_free(m);

    return 0;
}
