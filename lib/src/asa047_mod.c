/*
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

void nelmax (double fn(const double x[], void *data ), int n,
	     double start[], double xmin[], double *ynewlo,
	     double reqmin, double step[], int konvge, int kcount, 
	     int *icount, int *numres, int *ifault, void *data)
{
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
    double rq, x, z;
    double *y;
    double y2star;
    double ylo;
    double ystar;

    if (reqmin <= 0.0 || n < 1 || konvge < 1) {
	*ifault = 1;
	return;
    }

    p = malloc((n * (n + 1)) * sizeof *p);
    wspace = malloc((4*n + 1) * sizeof *wspace);

    if (p == NULL || wspace == NULL) {
	*ifault = 3;
	return;
    }

    pstar = wspace;
    p2star = pstar + n;
    pbar = p2star + n;
    y = pbar + n;

    *icount = *numres = 0;
    jcount = konvge; 

    nn = n + 1;
    rq = reqmin * n;
    
    /* Initial or restarted loop */

    for ( ; ; ) {
	for (i = 0; i < n; i++) {
	    p[i+n*n] = start[i];
	}
	y[n] = -fn(start, data);
	*icount += 1;

	for (j = 0; j < n; j++) {
	    x = start[j];
	    start[j] = start[j] + step[j] * del;
	    for (i = 0; i < n; i++) {
		p[i+j*n] = start[i];
	    }
	    y[j] = -fn(start, data);
	    *icount += 1;
	    start[j] = x;
	}
	/*                 
          The simplex construction is complete.
                    
          Find highest and lowest Y values.  YNEWLO = Y(IHI) indicates
          the vertex of the simplex to be replaced.
	*/                
	ylo = y[0];
	ilo = 0;

	for (i = 1; i < nn; i++) {
	    if (y[i] < ylo) {
		ylo = y[i];
		ilo = i;
	    }
	}
	/*
	  Inner loop
	*/
	for ( ; ; ) {
	    if (*icount >= kcount) {
		break;
	    }
	    *ynewlo = y[0];
	    ihi = 0;

	    for (i = 1; i < nn; i++) {
		if (*ynewlo < y[i]) {
		    *ynewlo = y[i];
		    ihi = i;
		}
	    }
	    /*
	      Calculate PBAR, the centroid of the simplex vertices
	      excepting the vertex with Y value YNEWLO.
	    */
	    for (i = 0; i < n; i++) {
		z = 0.0;
		for (j = 0; j < nn; j++) {
		    z = z + p[i+j*n];
		}
		z = z - p[i+ihi*n];  
		pbar[i] = z / n;
	    }
	    /*
	      Reflection through the centroid
	    */
	    for (i = 0; i < n; i++) {
		pstar[i] = pbar[i] + rcoeff * (pbar[i] - p[i+ihi*n]);
	    }
	    ystar = -fn(pstar, data);
	    *icount += 1;
	    /*
	      Successful reflection, so extension
	    */
	    if (ystar < ylo) {
		for (i = 0; i < n; i++) {
		    p2star[i] = pbar[i] + ecoeff * (pstar[i] - pbar[i]);
		}
		y2star = -fn(p2star, data);
		*icount += 1;
		/*
		  Check extension
		*/
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
		    y2star = -fn(p2star, data);
		    *icount += 1;
		    /*
		      Contract the whole simplex
		    */
		    if (y[ihi] < y2star) {
			for (j = 0; j < nn; j++) {
			    for (i = 0; i < n; i++) {
				p[i+j*n] = (p[i+j*n] + p[i+ilo*n]) * 0.5;
				xmin[i] = p[i+j*n];
			    }
			    y[j] = -fn(xmin, data);
			    *icount += 1;
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
		    y2star = -fn(p2star, data);
		    *icount += 1;
		    /*
		      Retain reflection?
		    */
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
	    /*
	      Check if YLO improved
	    */
	    if (y[ihi] < ylo) {
		ylo = y[ihi];
		ilo = ihi;
	    }
	    jcount--;

	    if (jcount > 0) {
		continue;
	    }
	    /*
	      Check to see if minimum reached
	    */
	    if (*icount <= kcount) {
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
#if 0
		    fprintf(stderr, "conv check satisfied: %g <= %g\n", z, rq);
#endif
		    break;
		}
	    }
	}
	
	/*
	  factorial tests to check that YNEWLO is a local minimum.
	*/
	for (i = 0; i < n; i++) {
	    xmin[i] = p[i+ilo*n];
	}
	*ynewlo = y[ilo];

	if (*icount > kcount) {
	    *ifault = 2;
	    break;
	}

	*ifault = 0;

	for (i = 0; i < n; i++) {
	    del = step[i] * eps;
	    xmin[i] = xmin[i] + del;
	    z = -fn(xmin, data);
	    *icount += 1;
	    if (z < *ynewlo) {
		*ifault = 2;
		break;
	    }
	    xmin[i] = xmin[i] - del - del;
	    z = -fn(xmin, data);
	    *icount += 1;
	    if (z < *ynewlo) {
		*ifault = 2;
		break;
	    }
	    xmin[i] = xmin[i] + del;
	}

	if (*ifault == 0) {
	    break;
	}
	
	/* restart the procedure */
	for (i = 0; i < n; i++) {
	    start[i] = xmin[i];
	}
	del = eps;
	*numres += 1;
    }
  
    free(p);
    free(wspace);

    return;
}
