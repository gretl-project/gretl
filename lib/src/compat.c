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
#include "libset.h"
#include "gretl_scalar.h"
#include "compat.h"

/*
    compat.c - older functions retained for backward compatibility
*/

/* ASCII graphics */

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

static void prntdate (int nt, int n, const DATAINFO *pdinfo, PRN *prn)
{
    if (n != pdinfo->t2 - pdinfo->t1 + 1) {
	pprintf(prn, "%5d  ", nt);
    } else {
	double xx = date(nt, pdinfo->pd, pdinfo->sd0);

	if (pdinfo->pd == 1) {
	    pprintf(prn, "%5.0f ", xx);
	} else if (pdinfo->pd < 10) {
	    pprintf(prn, "%5g ", xx);
	} else {
	    pprintf(prn, "%6.2f", xx);
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
		    const double **Z, const DATAINFO *pdinfo)
{
    int t, m = 0;

    for (t=pdinfo->t1; t<=pdinfo->t2; t++)  {
	if (na(Z[v1][t]) || na(Z[v2][t])) {
	    continue;
	}
	x[m] = Z[v1][t];
	y[m] = Z[v2][t];
	m++;
    }

    return m;
}

static int z_to_xyz (int v1, int v2, int v3, 
		     double *x, double *y, double *z,
		     const double **Z, const DATAINFO *pdinfo)
{
    int t, m = 0;

    for (t=pdinfo->t1; t<=pdinfo->t2; t++)  {
	if (na(Z[v1][t]) || na(Z[v2][t]) || na(Z[v3][t])) {
	    continue;
	}
	x[m] = Z[v1][t];
	y[m] = Z[v2][t];
	z[m] = Z[v3][t];
	m++;
    }

    return m;
}

/* graph one or two y variables against a given x variable */

static int graphyzx (const double *y1, const double *y2, const double *x, 
		     int n, const char *yname, const char *y2name,
		     const char *xname, gretlopt opt, PRN *prn)
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
    if (opt & OPT_O) {
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
	pprintf(prn, "%14s\n", yname);
    } else {
	pprintf(prn, _("%7co stands for %s and x stands for %s (+ means they "
		"are equal)\n\n%9s, %s\n"), ' ', 
		yname, y2name, yname, y2name);
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
	    pprintf(prn, "%c", p[i][j]);
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

/**
 * graphyx:
 * @y: y-axis data.
 * @x: x-axis data.
 * @n: number of observations.
 * @yname: y-axis label.
 * @xname: x-axis label.
 * @prn: gretl printing struct.
 *
 * Generates a simple ascii scatter-plot of @y against @x and 
 * prints the plot to @prn.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int graphyx (const double *y, const double *x, int n,
	     const char *yname, const char *xname, 
	     PRN *prn)
{
    return graphyzx(y, NULL, x, n, yname, NULL, xname, 
		    OPT_NONE, prn);
}

/* OPT_O: force the two variables to be plotted on the same; 
   otherwise they are scaled to fit
*/

static int ascii_plot (const int *list, const double **Z, 
		       const DATAINFO *pdinfo, gretlopt opt, 
		       PRN *prn)
{
    int i, nc2, vy, vz, cntrline;
    int ix = 0, iy = 0, iz = 0, n = 0;
    int ncols = 70;
    int ls, lineno, t, t1 = pdinfo->t1, t2 = pdinfo->t2;
    char word[32], px[132];
    char s1[10], s2[10];
    double xmin, xmax, xrange, ymin, ymax, yrange, xymin, xymax; 
    double xyrange, xx, yy;
    double *x, *y;

    int pause = gretl_get_text_pause();

    x = malloc(pdinfo->n * sizeof *x);
    y = malloc(pdinfo->n * sizeof *y);

    if (x == NULL || y == NULL) {
	return E_ALLOC;
    }

    pputc(prn, '\n');

    nc2 = ncols / 2;
    vy = list[1];
    strcpy(s1, pdinfo->varname[vy]);

    if (list[0] == 1) {
	/* only one variable is to be plotted */
	n = ztox(vy, x, Z, pdinfo);
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

	lineno = 1;
	for (t=t1; t<=t2; ++t) {
	    xx = Z[vy][t];
	    if (na(xx)) {
		continue;
	    }
	    if (pause && (lineno % PAGELINES == 0)) {
		scroll_pause();
		lineno = 1;
	    }
	    prntdate(t, n, pdinfo, prn);
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
	    lineno++;
	}

	gretl_pop_c_numeric_locale();

	pputs(prn, "\n\n");

	free(x);
	free(y);

	return 0;
    }

    /* two variables are to be plotted */
    vz = list[2];
    strcpy(s2, pdinfo->varname[vz]);
    n = z_to_xy(vy, vz, x, y, Z, pdinfo);

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
    lineno = 6;

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

    lineno = 1;

    for (t=t1; t<=t2; ++t) {
	if (pause && (lineno % PAGELINES == 0)) {
	    scroll_pause();
	    lineno = 1;
	}

	xx = Z[vy][t];
	yy = Z[vz][t];

	if (na(xx) || na(yy)) {
	    continue;
	}

	prntdate(t, n, pdinfo, prn);

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

	lineno++;
    }

    gretl_pop_c_numeric_locale();

    pputs(prn, "\n\n");

    free(x);
    free(y);

    return 0;
}

/**
 * ascii_graph:
 * @list: contains ID numbers of variables to graph.
 * @Z: data matrix.
 * @pdinfo: data information struct.
 * @opt: if includes OPT_O, use 40 rows, otherwise use 20 rows;
 * if includes OPT_T do a time-series plot.
 * @prn: gretl printing struct.
 *
 * Graph (using ascii graphics) one variable against another, as given
 * in @list: by default, the first variable will appear on the y-axis, 
 * the second on the x-axis.  But if @opt contains %OPT_T the listed
 * variables will be plotted against time (or by observation).
 *
 * Returns: 0 on successful completion, error code on error.
 */

int ascii_graph (const int *list, const double **Z, const DATAINFO *pdinfo, 
		 gretlopt opt, PRN *prn)
{
    int T = pdinfo->t2 - pdinfo->t1 + 1;
    int vx, vy1, vy2 = -1;
    double *x = NULL;
    double *y1 = NULL;
    double *y2 = NULL;
    int err = 0;

    if (opt & OPT_T) {
	return ascii_plot(list, Z, pdinfo, opt, prn);
    }

    if (list[0] < 2) {
	return E_ARGS; 
    }

    x = malloc(T * sizeof *x);
    y1 = malloc(T * sizeof *y1);

    if (x == NULL || y1 == NULL) {
	free(x);
	free(y1);
	return E_ALLOC;    
    }

    if (list[0] > 2) {
	y2 = malloc(T * sizeof *y2);
	if (y2 == NULL) {
	    free(x);
	    free(y1);
	    return E_ALLOC;    
	}
    }

    vy1 = list[1];

    /* Put values from Z array into x and y arrays */
    if (list[0] == 2) {
	vx = list[2];
	T = z_to_xy(vx, vy1, x, y1, Z, pdinfo);
    } else {
	vy2 = list[2];
	vx = list[3];
	T = z_to_xyz(vx, vy1, vy2, x, y1, y2, Z, pdinfo);
    }

    pputc(prn, '\n');

    err = graphyzx(y1, y2, x, T, pdinfo->varname[vy1], 
		   (vy2 < 0)? NULL : pdinfo->varname[vy2], 
		   pdinfo->varname[vx], opt, prn);

    pputc(prn, '\n');

    free(x); 
    free(y1); 
    free(y2);

    return err;
}

#undef RHODEBUG

/**
 * rhodiff:
 * @param: please see the gretl help on rhodiff() for syntax.
 * @list: list of variables to process.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 *
 * Generates and adds to the data set rho-differenced versions
 * of the variables given in @list.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int rhodiff (char *param, const int *list, double ***pZ, DATAINFO *pdinfo)
{
    int i, j, maxlag, p, t, t1, nv;
    int v = pdinfo->v, n = pdinfo->n;
    char parmbit[VNAMELEN];
    double *rhot;

#ifdef RHODEBUG
    fprintf(stderr, "rhodiff: param = '%s'\n", param);
#endif

    maxlag = count_fields(param);
    rhot = malloc(maxlag * sizeof *rhot);
    if (rhot == NULL) {
	return E_ALLOC;
    }

    if (maxlag > pdinfo->t1) {
	t1 = maxlag;
    } else {
	t1 = pdinfo->t1;
    }

#ifdef RHODEBUG
    fprintf(stderr, "rhodiff: maxlag = %d, t1 = %d\n", maxlag, t1);
#endif

    /* parse "param" string */
    j = strlen(param);
    p = 0;
    for (i=0; i<j; i++) {
	if ((i == 0 || param[i] == ' ') && i < (j - 1)) {
	    sscanf(param + i + (i? 1 : 0), "%15s", parmbit); 
#ifdef RHODEBUG
	    fprintf(stderr, "rhodiff: parmbit = '%s'\n", parmbit);
#endif
	    if (isalpha((unsigned char) parmbit[0])) {
		if (gretl_is_scalar(parmbit)) {
		    rhot[p] = gretl_scalar_get_value(parmbit, NULL);
		} else {
		    nv = varindex(pdinfo, parmbit);
		    if (nv == v) {
			free(rhot);
			return E_UNKVAR;
		    }
		    rhot[p] = get_xvalue(nv, (const double **) *pZ, pdinfo);
		}
	    } else {
		rhot[p] = dot_atof(parmbit);
	    }
	    p++;
	}
    }

    if (dataset_add_series(list[0], pZ, pdinfo)) {
	return E_ALLOC;
    }

    for (i=1; i<=list[0]; i++) {
	int vr = v + i - 1;
	double xx;

	j = list[i];

	*pdinfo->varname[vr] = 0;
	strncat(pdinfo->varname[vr], pdinfo->varname[j], VNAMELEN-3);
	strcat(pdinfo->varname[vr], "_r");

	sprintf(VARLABEL(pdinfo, vr), _("= rho-differenced %s"),
		pdinfo->varname[j]);

	for (t=0; t<n; t++) {
	    if (t < t1 || t > pdinfo->t2) {
		(*pZ)[vr][t] = NADBL;
		continue;
	    }
	    xx = (*pZ)[j][t];
	    if (na(xx)) {
		(*pZ)[vr][t] = NADBL;
		continue;
	    }
	    for (p=0; p<maxlag; p++) {
		if (na((*pZ)[j][t-p-1])) {
		    xx = NADBL;
		    break;
		} else {
		    xx -= rhot[p] * (*pZ)[j][t-p-1];
		}
	    }
	    (*pZ)[vr][t] = xx;
	}
    }

    free(rhot);

    return 0;
}
