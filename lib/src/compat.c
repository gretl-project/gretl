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
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

#include "libgretl.h"
#include "gretl_private.h"
#include "libset.h"
#include "compat.h"

/*
    compat.c - older functions retained for backward compatibility
*/

/* ASCII graphics */

static void initpx (int nn, char *pp)
{
    int i;

    pp[0] = '|';
    for (i=1; i<=nn; i++) {
	pp[i] = ' ';
    }
    pp[nn+1] = '\n';
}

static void drawline (int nn, PRN *prn)
{
    int t;

    pputs(prn, "       |+");

    for (t=1; t<=nn; t++) {
	if (t % 10 == 0) {
	    pputc(prn, '+');
	} else {
	    pputc(prn, '-');
	}
    }

    pputc(prn, '\n');
}

static void 
prntdate (int nt, int n, const DATAINFO *pdinfo, PRN *prn)
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

static void printgx (const double xx, PRN *prn)
{
    static char word[32];
    int lw;

    sprintf(word, "%11g", xx);
    lw = strlen(word);
    pputs(prn, word);
    bufspace(13 - lw, prn);
} 

void graphyzx (const int *list, const double *zy1, const double *zy2, 
	       const double *zx, int n, const char *yname, 
	       const char *xname, const DATAINFO *pdinfo, 
	       gretlopt oflag, PRN *prn)
/*
  if n > 0 graphs zy1 against zx, otherwise
  graphs zy1[i] and zy2[i] against zx[i] for i = 1, 2, .... n
  no of rows = 40 if oflag = OPT_O, else it is = 18 or 16
*/
{
    register int i, j;
    int ix, iy1, iy2, lx, ly, xzero, yzero, nrows, nr2, ncols, nc2,
	ls, lw, t1, t2, option = 0;
    double xmin, xmax, xrange, ymin, ymax, yrange, y1min, y1max; 
    double xx, y2min, y2max;
    char p[41][132];
    static char word[32];

    if (pdinfo != NULL) {
	t1 = pdinfo->t1;
	t2 = pdinfo->t2;
    } else {
	t1 = 0;
	t2 = (n < 0)? -n - 1 : n - 1;
    }

    if (n < 0) {
	n = -n;
	option = 1;
	gretl_minmax(t1, t2, zy1, &y1min, &y1max);
	gretl_minmax(t1, t2, zy2, &y2min, &y2max);
	ymin = (y1min < y2min)? y1min : y2min;
	ymax = (y1max > y2max)? y1max : y2max;
    } else {
	gretl_minmax(t1, t2, zy1, &ymin, &ymax);
    }

    yrange = ymax - ymin;
    xzero = yzero = 0;

    /* set the number of columns and rows to be used */
    ncols = 60;
    if (oflag & OPT_O) nrows = 40;
    else nrows = option ? 16 : 18 ;
    nr2 = nrows/2;
    nc2 = ncols/2;

    gretl_minmax(t1, t2, zx, &xmin, &xmax);
    xrange = xmax - xmin;

    /* Initialize picture matrix */
    for (i=0; i<=nrows; ++i) {
	p[i][0] = (i%5 == 0)? '+' : '|'; 
	for (j=1; j<=ncols+1; j++) p[i][j] = ' ';
    }

    if (xmin < 0 && xmax > 0) {
	xzero = 0.5 - xmin * ncols / xrange;
	for (i=0; i<=nrows; i++) p[i][xzero+1] = '|';
    }
    if (ymin < 0 && ymax > 0) {
	yzero = 0.5 - ymin * nrows / yrange;
	for (j=0; j<=ncols; j++) p[yzero][j+1] = '-';
    }

    /*  loop replaces blanks in PICTURE with o's that correspond to the
	scaled values of the specified variables */
    if (option) {
	for (i=0; i<n; ++i) {
	    ix = (floatneq(xrange, 0.0))? 
		((zx[i] - xmin)/xrange)*ncols : nc2;
	    iy1 = (floatneq(yrange, 0.0))? 
		((zy1[i] - ymin)/yrange)*nrows : nr2;
	    iy2 = (floatneq(yrange, 0.0))? 
		((zy2[i] - ymin)/yrange)*nrows : nr2;
	    if (iy1 != iy2) {
		p[iy1][ix+1] = 'o';
		p[iy2][ix+1] = 'x';
	    }
	    else p[iy1][ix+1] = '+';
	}
    } else for (i=0; i<n; ++i) {
	ix = (floatneq(xrange, 0.0))? 
	    ((zx[i] - xmin)/xrange)*ncols : nc2;
	iy1 = (floatneq(yrange, 0.0))? 
	    ((zy1[i] - ymin)/yrange)*nrows : nr2;
	p[iy1][ix+1] = 'o';
    }

    /* loop prints out the matrix PICTURE that is stored in the
       2-dimensional p matrix. */
    if (!option) {
	pprintf(prn, "%14s\n", yname);
    } else if (list) {
	pprintf(prn, _("%7co stands for %s and x stands for %s (+ means they "
		"are equal)\n\n%9s, %s\n"), ' ', 
		yname, pdinfo->varname[list[2]], yname, 
		pdinfo->varname[list[2]]);
    }

    for (i=nrows; i>=0; --i) {
	if (i && i == yzero) {
	    pputs(prn, "        0.0  ");
	} else if (i == nrows || i%5 == 0) {
	    xx = ymin + ((ymax-ymin) * i/nrows);
	    printgx(xx, prn);
	} else {
	    bufspace(13, prn);
	}
	for (j=0; j<=ncols+1; ++j) {
	    pprintf(prn, "%c", p[i][j]);
	}
	pputc(prn, '\n');
    }

    bufspace(13, prn);
    pputc(prn, '|');
    for (j=0; j<=ncols; j++) {
	if (j%10 == 0) pputc(prn, '+');
	else pputc(prn, '-');
    }
    pputc(prn, '\n');

    bufspace(14, prn);
    sprintf(word, "%g", xmin);
    lx = strlen(word);
    lw = 13 + lx;
    pputs(prn, word);
    sprintf(word, "%s", xname);
    ly = strlen(word);
    ls = 30 - lx - ly/2;
    bufspace(ls, prn);
    pputs(prn, word);
    lw = lw + ls + ly; 
    sprintf(word, "%g", xmax);

    ls = strlen(word);
    if (ls < 7) bufspace(73 - lw, prn);
    else { 
	lw = lw + ls;
	bufspace(79 - lw, prn);
    }
    pprintf(prn, "%s\n\n", word);
}

/**
 * ascii_plot:
 * @list: contains ID numbers of variables to plot.
 * @Z: data matrix.
 * @pdinfo: data information struct.
 * @oflag: if non-zero, forces two variables to be plotted on the same
 * scale (otherwise they will be scaled to fit).
 * @prn: gretl printing struct.
 *
 * Plot (using ascii graphics) either one or two variables, as given
 * in @list.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int ascii_plot (const LIST list, double **Z, const DATAINFO *pdinfo, 
		gretlopt oflag, PRN *prn)
/*
	plot var1 ;		plots var1 values
	plot var1 var2 ;	plots var1 and var2 values
	plot -o var1 var2 ;	plots var1 and var2 on same scale
				this is useful to plot observed and
				predicted values
	each row is an observation; the values are scaled horizontally
*/
{
    int  i, n, ix, iy=0, iz=0, ncols, nc2, l0, vy, vz, cntrline;
    int ls, lineno, t, t1 = pdinfo->t1, t2 = pdinfo->t2;
    char px[132];
    static char word[32];
    char s1[10], s2[10];
    double xmin, xmax, xrange, ymin, ymax, yrange, xymin, xymax; 
    double xyrange, xxx, yy;
    double *x, *y;

    int pause = gretl_get_text_pause();

    x = malloc(pdinfo->n * sizeof *x);
    y = malloc(pdinfo->n * sizeof *y);

    if (x == NULL || y == NULL) return E_ALLOC;

    l0 = list[0];
    pputc(prn, '\n');
    ncols = 70;
    nc2 = ncols/2;
    vy = list[1];
    strcpy(s1, pdinfo->varname[vy]);
    n = 0;

    if (l0 == 1) {
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
	ls = 78-ls-strlen(word);
	bufspace(ls, prn);
	pprintf(prn, "%s\n", word); 
	if (cntrline) {
	    iy = -(xmin / xrange) * ncols;
	    bufspace(iy+7, prn);
	    pputs(prn, "0.0\n"); 
	}
	drawline(ncols, prn);

#ifdef ENABLE_NLS
	setlocale(LC_NUMERIC, "C");
#endif

	lineno = 1;
	for (t=t1; t<=t2; ++t) {
	    xxx = Z[vy][t];
	    if (na(xxx)) continue;
	    if (pause && (lineno % PAGELINES == 0)) {
		takenotes(0);
		lineno = 1;
	    }
	    prntdate(t, n, pdinfo, prn);
	    ix = (floatneq(xrange, 0.0))? ((xxx-xmin)/xrange) * ncols : nc2;
	    initpx(ncols, px);
	    if (cntrline) px[iy+1] = '|';
	    px[ix+1] = 'o';
	    for (i=0; i<=ncols+1; i++) {
		pprintf(prn, "%c", px[i]); 
	    }
	    if (ix == ncols) pputc(prn, '\n');
	    lineno++;
	}

#ifdef ENABLE_NLS
	setlocale(LC_NUMERIC, "");
#endif
	pputs(prn, "\n\n");
	free(x);
	free(y);

	return 0;
    }

    /* two variables are to be plotted */
    vz = list[2];
    strcpy(s2, pdinfo->varname[vz]);
    n = z_to_xy(vy, vz, x, y, pdinfo, Z);
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

    if (oflag & OPT_O) {
	pprintf(prn, _("%20c%s and %s are plotted on same scale\n\n%8c"),
	       ' ', s1, s2, ' ');
	sprintf(word, "xy-min = %g", xymin);
	ls = 8 + strlen(word);
	pputs(prn, word);
	sprintf(word, "xy-max = %g", xymax);
	ls = 78-ls-strlen(word);
	bufspace(ls, prn);
	pprintf(prn, "%s\n", word);
    } else {
	pputc(prn, '\n');
	sprintf(word, "        o-min = %g", ymin);
	ls = strlen(word);
	pputs(prn, word);
	sprintf(word, "o-max = %g", ymax);
	ls = 78-ls-strlen(word);
	bufspace(ls, prn);
	pprintf(prn, "%s\n", word);
	sprintf(word, "        x-min = %g", xmin);
	ls = strlen(word);
	pputs(prn, word);
	sprintf(word, "x-max = %g", xmax);
	ls = 78-ls-strlen(word);
	bufspace(ls, prn);
	pprintf(prn, "%s\n", word);
    }

    /*  First x and y values are scaled, then we check to see which scaled
	value is smaller: print that one first, then print the larger
	scaled value. If the scaled values are equal then print a "+" ,
	otherwise prints an "x" for the first variable and an "o" for the
	second variable.
    */

    pputc(prn, '\n');
    cntrline = (floatgt(xymax, 0) && floatlt(xymin, 0))? 1 : 0;
    if (cntrline) {
	iz = (-xymin / xyrange) * ncols;
	bufspace(iz+7, prn);
	pputs(prn, "0.0\n");
    }
    drawline(ncols, prn);

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "C");
#endif

    lineno = 1;
    for (t=t1; t<=t2; ++t) {
	if (pause && (lineno % PAGELINES == 0)) {
	    takenotes(0);
	    lineno = 1;
	}

	xxx = Z[vy][t];
	yy = Z[vz][t];

	if (na(xxx) || na(yy)) continue;

	prntdate(t, n, pdinfo, prn);
	if (oflag & OPT_O) {
	    ix = (floatneq(xyrange, 0.0))? ((xxx-xymin)/xyrange) * ncols : nc2;
	    iy = (floatneq(xyrange, 0.0))? ((yy-xymin)/xyrange) * ncols : nc2;
	}
	else {
	    ix = (floatneq(xrange, 0.0))? ((xxx-xmin)/xrange) * ncols : nc2;
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

	if (ix == ncols || iy == ncols) pputc(prn, '\n');

	lineno++;
    }

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif

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
 * @oflag: if non-zero, use 40 rows, otherwise use 20 rows.
 * @prn: gretl printing struct.
 *
 * Graph (using ascii graphics) one variable against another, as given
 * in @list: the first variable will appear on the y-axis, the second
 * on the x-axis.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int ascii_graph (const LIST list, double **Z, const DATAINFO *pdinfo, 
		 gretlopt oflag, PRN *prn)
/*
  graph var1 var2 ;	graphs var1 (y-axis) against var2 (x-axis)
			in 20 rows and 60 columns
  graph -o var1 var2 ;	graphs var1 against var2 in 40 rows, 60 cols.
*/
{
    int m, vx, vy, vz, l0, t;
    int t1 = pdinfo->t1, t2 = pdinfo->t2;
    double *x, *y, xx, xy, xz, *uhat;

    if (list[0] < 2) return E_ARGS; 

    m = list_dups(list, GRAPH);
    if (m) {
	fprintf(stderr, _("var no. %d duplicated in command list.\n"), m);
	return 1;
    }

    pputc(prn, '\n');
    l0 = list[0];
    vy = list[1];

    x = malloc(pdinfo->n * sizeof *x);
    y = malloc(pdinfo->n * sizeof *y);
    uhat = malloc(pdinfo->n * sizeof *uhat);
    if (x == NULL || y == NULL || uhat == NULL) {
	return E_ALLOC;    
    }

    /* Put values from z matrix into x and y arrays */
    m = 0;
    if (l0 == 2) {
	vx = list[2];
	m = z_to_xy(vx, vy, x, y, pdinfo, Z);
	graphyzx(list, y, uhat, x, m, pdinfo->varname[vy], 
		 pdinfo->varname[vx], pdinfo, oflag, prn);
    }
    else {
	vz = list[2];
	vx = list[3];
	for(t=t1; t<=t2; t++) {
	    xx = Z[vx][t];
	    xy = Z[vy][t];
	    xz = Z[vz][t];
	    if (na(xx) || na(xy) || na(xz)) continue;
	    else {
		x[m] = xx;
		y[m] = xy;
		uhat[m] = xz;
		m++;
	    }
	}
	graphyzx(list, y, uhat, x, -m, pdinfo->varname[vy], 
		 pdinfo->varname[vx], pdinfo, oflag, prn);
    }
    pputc(prn, '\n');

    free(x); free(y); free(uhat);

    return 0;
}

/* end ASCII graphics */


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

int rhodiff (char *param, const LIST list, double ***pZ, DATAINFO *pdinfo)
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
	    sscanf(param + i + (i? 1 : 0), "%8s", parmbit); 
#ifdef RHODEBUG
	    fprintf(stderr, "rhodiff: parmbit = '%s'\n", parmbit);
#endif
	    if (isalpha((unsigned char) parmbit[0])) {
		nv = varindex(pdinfo, parmbit);
		if (nv == v) {
		    free(rhot);
		    return E_UNKVAR;
		}
		rhot[p] = get_xvalue(nv, (const double **) *pZ, pdinfo);
	    } else {
		rhot[p] = dot_atof(parmbit);
	    }
	    p++;
	}
    }

    if (dataset_add_vars(list[0], pZ, pdinfo)) {
	return E_ALLOC;
    }

    for (i=1; i<=list[0]; i++) {
	int vr = v + i - 1;
	double xx;

	j = list[i];

	strncat(pdinfo->varname[vr], pdinfo->varname[j], 7);
	strcat(pdinfo->varname[vr], "#");

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

/* stuff to support "sim" command */

/* ensure there's no gap betweeen a minus sign and the term
   it qualifies, in a "sim" command */

static char *close_minus (char *str)
{
    char *q, *p = str;

    while ((p = strchr(p, '-'))) {
	q = ++p;
	while (*q == ' ') q++;
	if (q != p) {
	    memmove(p, q, strlen(q) + 1);
	}
    }

    return str;
}

/* construct a descriptive label for a variable that is
   modified via the "sim" command */

static char *make_sim_label (char *label, const char *vname,
			     char **parm, int nparm)
{
    int k, neg, started = 0;
    char term[32];

    sprintf(label, "%s(t)=", vname);

    for (k=0; k<nparm; k++) {
	if (isdigit(*parm[k]) && dot_atof(parm[k]) == 0.0) {
	    continue;
	}
	if (k == 0) {
	    strcpy(term, parm[k]);
	} else {
	    neg = (*parm[k] == '-');
	    sprintf(term, "%s%s*%s(t-%d)", (neg)? "-" : 
		    (started)? "+" : "",
		    parm[k] + neg, vname, k);
	}
	if (strlen(label) + strlen(term) >= MAXLABEL - 4) {
	    if (strlen(label) < MAXLABEL - 4) {
		strcat(label, "...");
	    }
	    break;
	} else {
	    strcat(label, term);
	}
	started = 1;
    }
	
    return label;
}

int simulate (char *cmd, double ***pZ, DATAINFO *pdinfo)
{
    int f, i, t, t1, t2, tstart, m, nv = 0, pv;
    char varname[VNAMELEN], parm[16], tmpstr[MAXLEN];
    char *isconst = NULL, **toks = NULL;
    double xx, yy, *a = NULL;
    int vtok = 0, err = 0;

    *gretl_errmsg = '\0';

    close_minus(cmd);

    f = count_fields(cmd);
    m = f - 2; /* default: allow for command word varname */

    a = malloc(m * sizeof *a);
    isconst = malloc(m * sizeof *isconst);
    toks = malloc((f - 1) * sizeof *toks);

    if (a == NULL || isconst == NULL || toks == NULL) {
	err = E_ALLOC;
	goto sim_bailout;
    }

    for (i=0; i<m; i++) isconst[i] = 1;

    *tmpstr = 0;
    strncat(tmpstr, cmd, MAXLEN - 1);
    
    strtok(tmpstr, " "); /* discard the "sim" command word */
    for (i=0; i<f-1; i++) {
	toks[i] = strtok(NULL, " ");
    }

    /* allow for implicit starting and ending dates */
    if (isalpha(*toks[0])) {
	t1 = pdinfo->t1;
	t2 = pdinfo->t2;
    } else {
	m -= 2;
	vtok = 2;
	/* try getting valid obs from stobs and endobs */
	t1 = dateton(toks[0], pdinfo);
	t2 = dateton(toks[1], pdinfo);
	if (*gretl_errmsg || t1 < 0 || t1 > t2 || t2 >= pdinfo->n) {

	    if (t1 < 0 || t2 >= pdinfo->n) {
		strcpy(gretl_errmsg, _("Observation number out of bounds"));
	    } else if (t1 > t2 ) {
		strcpy(gretl_errmsg, _("Invalid null sample"));
	    }

	    err = 1;
	    goto sim_bailout;
	}
    }

    /* name of var to simulate */
    *varname = 0;
    strncat(varname, toks[vtok], 8);
    nv = varindex(pdinfo, varname);

    if (nv > 0 && nv < pdinfo->v && pdinfo->vector[nv] == 0) {
	sprintf(gretl_errmsg, _("variable %s is a scalar"), 
		pdinfo->varname[nv]);
	err = 1;
	goto sim_bailout;
    }
		
    if (nv == 0 || nv >= pdinfo->v) {
	sprintf(gretl_errmsg, (nv)? _("For 'sim', the variable must already "
				      "exist") :
		_("You can't use the constant for this purpose"));
	err = 1;
	goto sim_bailout;
    }

    /* get the parameter terms */
    for (i=0; i<m; i++) {
	int neg = 0;
	const char *p = parm;

	*parm = '\0';
	strncat(parm, toks[i + vtok + 1], sizeof parm - 1);

	if (*parm == '-') {
	    neg = 1;
	    p++;
	}

	if (isalpha((unsigned char) *p)) {
	    pv = varindex(pdinfo, p);
	    if (pv == 0 || pv >= pdinfo->v) {
		sprintf(gretl_errmsg, _("Bad varname '%s' in sim"), p);
		err = 1;
		goto sim_bailout;
	    } else {
		isconst[i] = !pdinfo->vector[pv];
		a[i] = (isconst[i])? (*pZ)[pv][0] : (double) pv;
	    }
	} else {
	    a[i] = dot_atof(p);
	}

	if (neg) a[i] = -a[i];
    }

    tstart = t1;
    if (tstart < m - 1) tstart = m - 1;

    for (t=tstart; t<=t2; t++) {
	xx = 0.;
	for (i=0; i<m; i++) {
	    if (isconst[i]) {
		if (i == 0) xx += a[i];
		else xx += a[i] * (*pZ)[nv][t-i];
	    } else {
		int neg = 0;

		pv = (int) a[i];
		if (pv < 0) {
		    neg = 1;
		    pv = -pv;
		}
		yy = (*pZ)[pv][t];
		if (na(yy)) {
		    xx = NADBL;
		    break;
		}
		if (neg) yy = -yy;
		if (i == 0) xx += yy;
		else xx += yy * (*pZ)[nv][t-i];
	    }
	}
	(*pZ)[nv][t] = xx;
    }

 sim_bailout:

    if (!err && nv > 0) {
	sprintf(gretl_msg, "%s %s %s (ID %d)", 
		_("Replaced"), _("vector"), pdinfo->varname[nv], nv);
	make_sim_label(VARLABEL(pdinfo, nv), pdinfo->varname[nv],
		       toks + vtok + 1, m);
    }

    free(a);
    free(isconst);
    free(toks);

    return err;
}

/* .......................................................... */

int gretl_multiply (char *s, int *list, char *sfx, double ***pZ,
		    DATAINFO *pdinfo)
{
    int i, t, v = 0, nv, n = pdinfo->n, lv, l0 = list[0];
    int slen;
    double m = 0;
    char tmp[VNAMELEN];

    /* parse s */
    if (isdigit((unsigned char) *s)) {
	m = dot_atof(s);
    } else {
	v = varindex(pdinfo, s);
	if (v == pdinfo->v) return E_UNKVAR; 
    }

    if (dataset_add_vars(l0, pZ, pdinfo)) {
	return E_ALLOC;
    }

    slen = strlen(sfx);

    /* fill out values */
    for (i=1; i<=l0; i++) {
	nv = pdinfo->v - l0 - 1 + i;
	lv = list[i];

	for (t=0; t<n; t++) {
	    (*pZ)[nv][t] = NADBL;
	}

	for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	    if (na((*pZ)[lv][t])) {
		(*pZ)[nv][t] = NADBL;
		continue;
	    }
	    if (v) {
		double yy = (pdinfo->vector[v])? 
		    (*pZ)[v][t] : (*pZ)[v][0];

		if (na(yy)) {
		    (*pZ)[nv][t] = NADBL;
		} else {
		    (*pZ)[nv][t] = yy * (*pZ)[lv][t];
		}
	    } else {
		(*pZ)[nv][t] = m * (*pZ)[lv][t];
	    }
	}

	/* do names and labels */
	strcpy(tmp, pdinfo->varname[lv]);
	gretl_trunc(tmp, 8 - slen);
	strcat(tmp, sfx);
	strcpy(pdinfo->varname[nv], tmp);
	if (v) {
	    sprintf(VARLABEL(pdinfo, nv), "%s = %s * %s",
		    pdinfo->varname[nv], pdinfo->varname[v], 
		    pdinfo->varname[lv]); 
	} else {
	    sprintf(VARLABEL(pdinfo, nv), "%s = %g * %s",
		    pdinfo->varname[nv], m, pdinfo->varname[lv]); 
	}
    }

    return 0;
}
