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
#include "compat.h"

/*
    compat.c -  ASCII graphics for backward compatibility
*/

static void initpx (int n, char *pp)
{
    int i;

    pp[0] = '|';
    for (i=1; i<=n; i++) {
	pp[i] = ' ';
    }
    pp[n + 1] = '\n';
}

static void drawline (int n, PRN *prn)
{
    int t;

    pputs(prn, "       |+");
    for (t=1; t<=n; t++) {
	pputc(prn, (t % 10)? '-' : '+');
    }
    pputc(prn, '\n');
}

static void prntdate (int nt, int n, const DATASET *dset, PRN *prn)
{
    if (n != dset->t2 - dset->t1 + 1) {
	pprintf(prn, "%5d  ", nt);
    } else {
	double x = date_as_double(nt, dset->pd, dset->sd0);

	if (dset->pd == 1) {
	    pprintf(prn, "%5.0f ", x);
	} else if (dset->pd < 10) {
	    pprintf(prn, "%5g ", x);
	} else {
	    pprintf(prn, "%6.2f", x);
	}
    }
}

static void printgx (double x, PRN *prn)
{
    char word[32];
    int lw;

    sprintf(word, "%11g", x);
    lw = strlen(word);
    pputs(prn, word);
    bufspace(13 - lw, prn);
}

static int z_to_xy (int v1, int v2, double *x, double *y,
		    const DATASET *dset)
{
    int t, m = 0;

    for (t=dset->t1; t<=dset->t2; t++)  {
	if (na(dset->Z[v1][t]) || na(dset->Z[v2][t])) {
	    continue;
	}
	x[m] = dset->Z[v1][t];
	y[m] = dset->Z[v2][t];
	m++;
    }

    return m;
}

static int z_to_xyz (int v1, int v2, int v3,
		     double *x, double *y, double *z,
		     const DATASET *dset)
{
    int t, m = 0;

    for (t=dset->t1; t<=dset->t2; t++)  {
	if (na(dset->Z[v1][t]) || na(dset->Z[v2][t]) || na(dset->Z[v3][t])) {
	    continue;
	}
	x[m] = dset->Z[v1][t];
	y[m] = dset->Z[v2][t];
	z[m] = dset->Z[v3][t];
	m++;
    }

    return m;
}

/* graph one or two y variables against a given x variable */

static int graphyzx (const double *y1, const double *y2,
		     const double *x, int n,
		     const char *yname,
		     const char *y2name,
		     const char *xname,
		     gretlopt opt, PRN *prn)
{
    int ix, iy1, iy2, lx, ly, xzero, yzero;
    int nrows, nr2, nc2, ls, lw, t1, t2;
    double xx, xmin, xmax, xrange;
    double ymin, ymax, yrange;
    int ncols = 60;
    char p[41][132];
    char word[32];
    int i, j;

    t1 = 0;
    t2 = n - 1;

    if (y2 != NULL) {
	double y1min, y1max;
	double y2min, y2max;

	gretl_minmax(t1, t2, y1, &y1min, &y1max);
	gretl_minmax(t1, t2, y2, &y2min, &y2max);
	ymin = (y1min < y2min)? y1min : y2min;
	ymax = (y1max > y2max)? y1max : y2max;
    } else {
	gretl_minmax(t1, t2, y1, &ymin, &ymax);
    }

    if (na(ymin) || na(ymax)) {
	return E_MISSDATA;
    }

    yrange = ymax - ymin;
    xzero = yzero = 0;

    /* set the number of rows to be used */
    if (opt & OPT_T) {
	nrows = 40;
    } else {
	nrows = (y2 != NULL)? 16 : 18;
    }

    nr2 = nrows / 2;
    nc2 = ncols / 2;

    gretl_minmax(t1, t2, x, &xmin, &xmax);

    if (na(xmin) || na(xmax)) {
	return E_MISSDATA;
    }

    xrange = xmax - xmin;

    /* Initialize picture matrix */

    for (i=0; i<=nrows; i++) {
	p[i][0] = (i % 5 == 0)? '+' : '|';
	for (j=1; j<=ncols+1; j++) {
	    p[i][j] = ' ';
	}
    }

    if (xmin < 0 && xmax > 0) {
	xzero = 0.5 - xmin * ncols / xrange;
	for (i=0; i<=nrows; i++) {
	    p[i][xzero+1] = '|';
	}
    }

    if (ymin < 0 && ymax > 0) {
	yzero = 0.5 - ymin * nrows / yrange;
	for (j=0; j<=ncols; j++) {
	    p[yzero][j+1] = '-';
	}
    }

    /* replace blanks in PICTURE with o's that correspond to the
       scaled values of the specified variables */

    if (y2 != NULL) {
	for (i=0; i<n; i++) {
	    ix = (floatneq(xrange, 0.0))? ((x[i] - xmin)/xrange) * ncols : nc2;
	    iy1 = (floatneq(yrange, 0.0))? ((y1[i] - ymin)/yrange) * nrows : nr2;
	    iy2 = (floatneq(yrange, 0.0))? ((y2[i] - ymin)/yrange) * nrows : nr2;
	    if (iy1 != iy2) {
		p[iy1][ix+1] = 'o';
		p[iy2][ix+1] = 'x';
	    } else {
		p[iy1][ix+1] = '+';
	    }
	}
    } else {
	for (i=0; i<n; i++) {
	    ix = (floatneq(xrange, 0.0))? ((x[i] - xmin)/xrange) * ncols : nc2;
	    iy1 = (floatneq(yrange, 0.0))? ((y1[i] - ymin)/yrange) * nrows : nr2;
	    p[iy1][ix+1] = 'o';
	}
    }

    /* print out the 2-dimensional picture matrix */

    if (y2name == NULL) {
	pprintf(prn, "  %s\n", yname);
    } else {
	pputs(prn, "  ");
	pprintf(prn, _("'o' stands for %s and 'x' stands for %s (+ means they "
		       "are equal)"), yname, y2name);
	pprintf(prn, "\n\n  %s, %s\n", yname, y2name);
    }

    for (i=nrows; i>=0; i--) {
	if (i && i == yzero) {
	    pputs(prn, "        0.0  ");
	} else if (i == nrows || i % 5 == 0) {
	    xx = ymin + ((ymax - ymin) * i/nrows);
	    printgx(xx, prn);
	} else {
	    bufspace(13, prn);
	}
	for (j=0; j<=ncols+1; j++) {
	    pputc(prn, p[i][j]);
	}
	pputc(prn, '\n');
    }

    bufspace(13, prn);
    pputc(prn, '|');
    for (j=0; j<=ncols; j++) {
	pputc(prn, (j % 10)? '-' : '+');
    }

    pputc(prn, '\n');
    bufspace(14, prn);

    sprintf(word, "%g", xmin);
    lx = strlen(word);
    lw = 13 + lx;
    pputs(prn, word);

    sprintf(word, "%s", xname);
    ly = strlen(word);
    ls = 30 - lx - ly / 2;
    bufspace(ls, prn);
    pputs(prn, word);

    lw = lw + ls + ly;

    sprintf(word, "%g", xmax);
    ls = strlen(word);
    if (ls < 7) {
	bufspace(73 - lw, prn);
    } else {
	lw = lw + ls;
	bufspace(79 - lw, prn);
    }
    pprintf(prn, "%s\n\n", word);

    return 0;
}

/* OPT_O: force the two variables to be plotted on the same scale;
   otherwise they are scaled to fit
*/

static int ascii_plot (const int *list, const DATASET *dset,
		       gretlopt opt, PRN *prn)
{
    int i, nc2, vy, vz, cntrline;
    int ix = 0, iy = 0, iz = 0, n = 0;
    int ncols = 70;
    int t1 = dset->t1;
    int t2 = dset->t2;
    char word[32], px[132];
    char s1[10], s2[10];
    double xmin, xmax, xrange, ymin, ymax, yrange, xymin, xymax;
    double xyrange, xx, yy;
    double *x, *y;
    int ls, t;

    x = malloc(2 * dset->n * sizeof *x);
    if (x == NULL) {
	return E_ALLOC;
    }

    y = x + dset->n;

    pputc(prn, '\n');

    nc2 = ncols / 2;
    vy = list[1];
    strcpy(s1, dset->varname[vy]);

    if (list[0] == 1) {
	/* only one variable is to be plotted */
	n = transcribe_array(x, dset->Z[vy], dset);
	gretl_minmax(t1, t2, x, &xmin, &xmax);
	xrange = xmax - xmin;
	cntrline = (floatgt(xmax, 0) && floatlt(xmin, 0))? 1 : 0;

	/* print headings */
	pprintf(prn, _("%25cNOTE: o stands for %s\n\n%8c"), ' ', s1, ' ');

	sprintf(word, "x-min = %g", xmin);
	ls = 8 + strlen(word);
	pputs(prn, word);

	sprintf(word, "x-max = %g", xmax);
	ls = 78 - ls - strlen(word);
	bufspace(ls, prn);
	pprintf(prn, "%s\n", word);

	if (cntrline) {
	    iy = -(xmin / xrange) * ncols;
	    bufspace(iy + 7, prn);
	    pputs(prn, "0.0\n");
	}

	drawline(ncols, prn);

	gretl_push_c_numeric_locale();

	for (t=t1; t<=t2; ++t) {
	    xx = dset->Z[vy][t];
	    if (na(xx)) {
		continue;
	    }
	    prntdate(t, n, dset, prn);
	    ix = (floatneq(xrange, 0.0))? ((xx-xmin)/xrange) * ncols : nc2;
	    initpx(ncols, px);
	    if (cntrline) {
		px[iy+1] = '|';
	    }
	    px[ix+1] = 'o';
	    for (i=0; i<=ncols+1; i++) {
		pprintf(prn, "%c", px[i]);
	    }
	    if (ix == ncols) {
		pputc(prn, '\n');
	    }
	}

	gretl_pop_c_numeric_locale();

	pputs(prn, "\n\n");

	free(x);

	return 0;
    }

    /* two variables are to be plotted */
    vz = list[2];
    strcpy(s2, dset->varname[vz]);
    n = z_to_xy(vy, vz, x, y, dset);

    /* find maximum and minimum using all values from both arrays */
    gretl_minmax(t1, t2, x, &xmin, &xmax);
    xrange = xmax - xmin;
    gretl_minmax(t1, t2, y, &ymin, &ymax);
    yrange = ymax - ymin;
    xymin = (xmin <= ymin) ? xmin : ymin;
    xymax = (xmax >= ymax) ? xmax : ymax;
    xyrange = xymax - xymin;

    /* print headings for the two variables */
    pprintf(prn, _("%17cNOTE: o stands for %s,   x stands for %s\n%17c+ means %s "
	   "and %s are equal when scaled\n"), ' ', s1, s2, ' ', s1, s2);

    if (opt & OPT_O) {
	pprintf(prn, _("%20c%s and %s are plotted on same scale\n\n%8c"),
	       ' ', s1, s2, ' ');

	sprintf(word, "xy-min = %g", xymin);
	ls = 8 + strlen(word);
	pputs(prn, word);

	sprintf(word, "xy-max = %g", xymax);
	ls = 78 - ls - strlen(word);
	bufspace(ls, prn);
	pprintf(prn, "%s\n", word);
    } else {
	pputc(prn, '\n');

	sprintf(word, "        o-min = %g", ymin);
	ls = strlen(word);
	pputs(prn, word);

	sprintf(word, "o-max = %g", ymax);
	ls = 78 - ls - strlen(word);
	bufspace(ls, prn);
	pprintf(prn, "%s\n", word);

	sprintf(word, "        x-min = %g", xmin);
	ls = strlen(word);
	pputs(prn, word);

	sprintf(word, "x-max = %g", xmax);
	ls = 78 - ls - strlen(word);
	bufspace(ls, prn);
	pprintf(prn, "%s\n", word);
    }

    /*  First x and y values are scaled, then we check to see which scaled
	value is smaller: print that one first, then print the larger
	scaled value. If the scaled values are equal then print a "+" ,
	otherwise print an "x" for the first variable and an "o" for the
	second variable.
    */

    pputc(prn, '\n');

    cntrline = (floatgt(xymax, 0) && floatlt(xymin, 0))? 1 : 0;
    if (cntrline) {
	iz = (-xymin / xyrange) * ncols;
	bufspace(iz + 7, prn);
	pputs(prn, "0.0\n");
    }

    drawline(ncols, prn);

    gretl_push_c_numeric_locale();

    for (t=t1; t<=t2; ++t) {

	xx = dset->Z[vy][t];
	yy = dset->Z[vz][t];

	if (na(xx) || na(yy)) {
	    continue;
	}

	prntdate(t, n, dset, prn);

	if (opt & OPT_O) {
	    ix = (floatneq(xyrange, 0.0))? ((xx-xymin)/xyrange) * ncols : nc2;
	    iy = (floatneq(xyrange, 0.0))? ((yy-xymin)/xyrange) * ncols : nc2;
	} else {
	    ix = (floatneq(xrange, 0.0))? ((xx-xmin)/xrange) * ncols : nc2;
	    iy = (floatneq(yrange, 0.0))? ((yy-ymin)/yrange) * ncols : nc2;
	}

	initpx(ncols, px);

	if (iz) {
	    px[iz+1] = '|';
	}

	if (ix == iy) {
	    px[ix+1] = '+';
	} else {
	    px[ix+1] = 'o';
	    px[iy+1] = 'x';
	}

	for (i=0; i<=ncols+1; i++) {
	    pprintf(prn, "%c", px[i]);
	}

	if (ix == ncols || iy == ncols) {
	    pputc(prn, '\n');
	}
    }

    gretl_pop_c_numeric_locale();

    pputs(prn, "\n\n");

    free(x);

    return 0;
}

static int
ascii_scatter (const int *list, const DATASET *dset,
	       gretlopt opt, PRN *prn)
{
    int T = sample_size(dset);
    int vx, vy1, vy2 = -1;
    double *x = NULL;
    double *y1, *y2 = NULL;
    int err = 0;

    if (list[0] > 2) {
	x = malloc(3 * T * sizeof *x);
    } else {
	x = malloc(2 * T * sizeof *x);
    }

    if (x == NULL) {
	return E_ALLOC;
    }

    y1 = x + T;
    if (list[0] > 2) {
	y2 = y1 + T;
    }

    vy1 = list[1];

    /* Put values from Z array into x and y arrays */
    if (list[0] == 2) {
	vx = list[2];
	T = z_to_xy(vx, vy1, x, y1, dset);
    } else {
	vy2 = list[2];
	vx = list[3];
	T = z_to_xyz(vx, vy1, vy2, x, y1, y2, dset);
    }

    pputc(prn, '\n');
    err = graphyzx(y1, y2, x, T, dset->varname[vy1],
		   (vy2 < 0)? NULL : dset->varname[vy2],
		   dset->varname[vx], opt, prn);
    pputc(prn, '\n');

    free(x);

    return err;
}

/**
 * textplot:
 * @list: list of series to plot.
 * @dset: dataset struct.
 * @opt: see below.
 * @prn: gretl printing struct.
 *
 * Produces, using ascii graphics, either a scatter plot of
 * the first k-1 series in @list against the kth, or (if @opt
 * includes %OPT_S) a time-series plot (or more generally
 * a plot by observation).
 *
 * In the case of a scatter plot, %OPT_T (tall) can be used to
 * request the use of 40 rows rather than the default of 20.
 * In the case of a plot by observation %OPT_O (one-scale)
 * can be used to force the use of a single scale (otherwise
 * the series may be scaled to fit).
 *
 * Returns: 0 on successful completion, error code on error.
 */

int textplot (const int *list, const DATASET *dset,
	      gretlopt opt, PRN *prn)
{
    if (opt & OPT_S) {
	/* time series */
	return ascii_plot(list, dset, opt, prn);
    } else if (list[0] < 2) {
	return E_ARGS;
    } else {
	return ascii_scatter(list, dset, opt, prn);
    }
}
