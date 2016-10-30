/*
  nelmax: this is based closely on the nelmin function as written
  in C by John Burkardt; see
  http://people.sc.fsu.edu/~jburkardt/c_src/asa047
  It is converted to a maximizer, and modified for use with
  the gretl library.

  Here follows a portion of Burkardt's original notice:

  Purpose:

  NELMIN minimizes a function using the Nelder-Mead algorithm.

  Discussion:

  This routine seeks the minimum value of a user-specified function.

  Simplex function minimisation procedure due to Nelder+Mead(1965),
  as implemented by O'Neill(1971, Appl.Statist. 20, 338-45), with
  subsequent comments by Chambers+Ertel(1974, 23, 250-1), Benyon(1976,
  25, 97) and Hill(1978, 27, 380-2)

  The function to be minimized must be defined by a function of
  the form

  function fn ( x, f )
  double fn
  double x(*)

  and the name of this subroutine must be declared EXTERNAL in the
  calling routine and passed as the argument FN.

  This routine does not include a termination test using the
  fitting of a quadratic surface.

  Licensing:

  This code is distributed under the GNU LGPL license. 

  Modified:

  28 October 2010

  Author:

  Original FORTRAN77 version by R ONeill.
  C version by John Burkardt.
*/

static double nm_call (BFGS_CRIT_FUNC cfunc,
		       double *b, void *data,
		       int *ncalls, int minimize)
{
    double ret = cfunc(b, data);

    *ncalls += 1;

    return minimize ? ret : na(ret) ? ret : -ret;
}

static int
nelder_mead (BFGS_CRIT_FUNC cfunc, int n, double start[], double xmin[],
	     double *ynewlo, double reqmin, double step[], int konvge,
	     int maxcalls, int *ncalls, int *nresets, void *data,
	     gretlopt opt, PRN *prn)
{
    gretl_matrix *pmat = NULL;
    double ccoeff = 0.5;
    double ecoeff = 2.0;
    double rcoeff = 1.0;
    double eps = 0.001;
    double del = 1.0;
    int i, ihi, ilo;
    int j, jcount;
    int l, nn;
    double *wspace;
    double *p;
    double *p2star;
    double *pbar;
    double *pstar;
    double *y;
    double rq, x, z;
    double y2star;
    double ylo;
    double ystar;
    int outer, inner;
    int getmin;
    int err = 0;

    if (reqmin <= 0.0 || n < 1 || konvge < 1) {
	return E_INVARG;
    }

    /* use a gretl_matrix in case we want to print it */
    pmat = gretl_matrix_alloc(n, n + 1);
    wspace = malloc((4*n + 1) * sizeof *wspace);

    if (pmat == NULL || wspace == NULL) {
	gretl_matrix_free(pmat);
	free(wspace);
	return E_ALLOC;
    }

    p = pmat->val;
    pstar = wspace;
    p2star = pstar + n;
    pbar = p2star + n;
    y = pbar + n;

    *ncalls = *nresets = 0;
    jcount = konvge; 

    nn = n + 1;
    rq = reqmin * n;

    /* maximize by default, but minimize if OPT_I is given */
    getmin = (opt & OPT_I)? 1 : 0;
    
    /* Initial or restarted loop */

    for (outer=1; ; outer++) {
	for (i = 0; i < n; i++) {
	    p[i+n*n] = start[i];
	}
	y[n] = nm_call(cfunc, start, data, ncalls, getmin);

	if (opt & OPT_V) {
	    if (outer == 1) {
		pprintf(prn, "\nNelder-Mead outer iteration %d: function value = %#g\n",
			outer, y[n]);
	    } else {
		pprintf(prn, "Outer iteration %d (reset)\n", outer);
	    }
	}

	/* construct the simplex */
	for (j = 0; j < n; j++) {
	    x = start[j];
	    start[j] += step[j] * del;
	    for (i = 0; i < n; i++) {
		p[i+j*n] = start[i];
	    }
	    y[j] = nm_call(cfunc, start, data, ncalls, getmin);
	    start[j] = x;
	}

	/* find the lowest y value */
	ylo = y[0];
	ilo = 0;
	for (i = 1; i < nn; i++) {
	    if (y[i] < ylo) {
		ylo = y[i];
		ilo = i;
	    }
	}

	/* Inner loop */

	for (inner=1; *ncalls < maxcalls; inner++) {
	    *ynewlo = y[0];
	    ihi = 0;

	    for (i = 1; i < nn; i++) {
		if (*ynewlo < y[i]) {
		    *ynewlo = y[i];
		    ihi = i;
		}
	    }
	    /* Calculate pbar, the centroid of the simplex vertices
	       excepting the vertex with y-value ynewlo.
	    */
	    for (i = 0; i < n; i++) {
		z = 0.0;
		for (j = 0; j < nn; j++) {
		    z += p[i+j*n];
		}
		z -= p[i+ihi*n];
		pbar[i] = z / n;
	    }

	    /* Reflection through the centroid */
	    for (i = 0; i < n; i++) {
		pstar[i] = pbar[i] + rcoeff * (pbar[i] - p[i+ihi*n]);
	    }
	    ystar = nm_call(cfunc, pstar, data, ncalls, getmin);

	    if ((opt & OPT_V) && (inner == 1 || inner % 10 == 0)) {
		pprintf(prn, " inner iter %3d: function value %#g\n",
			inner, ystar);
	    }

	    if (ystar < ylo) {
		/* Successful reflection, so extension */
		for (i = 0; i < n; i++) {
		    p2star[i] = pbar[i] + ecoeff * (pstar[i] - pbar[i]);
		}
		y2star = nm_call(cfunc, p2star, data, ncalls, getmin);
		/* Check extension */
		if (ystar < y2star) {
		    for (i = 0; i < n; i++) {
			p[i+ihi*n] = pstar[i];
		    }
		    y[ihi] = ystar;
		} else {
		    for (i = 0; i < n; i++) {
			p[i+ihi*n] = p2star[i];
		    }
		    y[ihi] = y2star;
		}
	    } else {
		l = 0;
		for (i = 0; i < nn; i++) {
		    if (ystar < y[i]) {
			l++;
		    }
		}
		if (l > 1) {
		    for (i = 0; i < n; i++) {
			p[i+ihi*n] = pstar[i];
		    }
		    y[ihi] = ystar;
		} else if (l == 0) {
		    for (i = 0; i < n; i++) {
			p2star[i] = pbar[i] + ccoeff * (p[i+ihi*n] - pbar[i]);
		    }
		    y2star = nm_call(cfunc, p2star, data, ncalls, getmin);
		    /* Contract the whole simplex */
		    if (y[ihi] < y2star) {
			for (j = 0; j < nn; j++) {
			    for (i = 0; i < n; i++) {
				p[i+j*n] = (p[i+j*n] + p[i+ilo*n]) * 0.5;
				xmin[i] = p[i+j*n];
			    }
			    y[j] = nm_call(cfunc, xmin, data, ncalls, getmin);
			}
			ylo = y[0];
			ilo = 0;
			for (i = 1; i < nn; i++) {
			    if (y[i] < ylo) {
				ylo = y[i];
				ilo = i;
			    }
			}
			continue;
		    } else {
			for (i = 0; i < n; i++) {
			    p[i+ihi*n] = p2star[i];
			}
			y[ihi] = y2star;
		    }
		} else if (l == 1) {
		    for (i = 0; i < n; i++) {
			p2star[i] = pbar[i] + ccoeff * (pstar[i] - pbar[i]);
		    }
		    y2star = nm_call(cfunc, p2star, data, ncalls, getmin);
		    /* Retain reflection? */
		    if (y2star <= ystar) {
			for (i = 0; i < n; i++) {
			    p[i+ihi*n] = p2star[i];
			}
			y[ihi] = y2star;
		    } else {
			for (i = 0; i < n; i++) {
			    p[i+ihi*n] = pstar[i];
			}
			y[ihi] = ystar;
		    }
		}
	    }

	    /* Check if ylo improved */
	    if (y[ihi] < ylo) {
		ylo = y[ihi];
		ilo = ihi;
	    }
	    jcount--;

	    if (jcount > 0) {
		continue;
	    }
	    
	    /* Check to see if minimum reached */
	    if (*ncalls <= maxcalls) {
		jcount = konvge;
		z = 0.0;
		for (i = 0; i < nn; i++) {
		    z += y[i];
		}
		x = z / nn;
		z = 0.0;
		for (i = 0; i < nn; i++) {
		    z += (y[i] - x) * (y[i] - x);
		}
		if (z <= rq) {
		    break;
		}
	    }
	}
	
	/* Factorial tests to check that ynewlo is a
	   local minimum */
	for (i = 0; i < n; i++) {
	    xmin[i] = p[i+ilo*n];
	}
	*ynewlo = y[ilo];

	if (*ncalls > maxcalls) {
	    err = E_NOCONV;
	    break;
	}

	err = 0;

	for (i = 0; i < n; i++) {
	    double xsave = xmin[i];
	    double dx = step[i] * eps;

	    xmin[i] = xsave + dx;
	    z = nm_call(cfunc, xmin, data, ncalls, getmin);
	    if (z < *ynewlo) {
		err = E_NOCONV;
		break;
	    }
	    xmin[i] = xsave - dx;
	    z = nm_call(cfunc, xmin, data, ncalls, getmin);
	    if (z < *ynewlo) {
		err = E_NOCONV;
		break;
	    }
	    xmin[i] = xsave;
	}

	if (err == 0) {
	    if (opt & OPT_V) {
		pprintf(prn, "Found optimum %#g after %d function calls, "
			"%d resets\n\n", *ynewlo, *ncalls, *nresets);
	    }
	    break;
	}
	
	/* prepare to restart the procedure */
	for (i = 0; i < n; i++) {
	    start[i] = xmin[i];
	}
	del = eps;
	*nresets += 1;
    }
  
    gretl_matrix_free(pmat);
    free(wspace);

    return err;
}
