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

/* graphing.c for gretl */

#ifdef OS_WIN32
# include "../../winconfig.h"
#else
# include "../../config.h"
#endif
#include "libgretl.h"
#include "internal.h"
#include <unistd.h>

#ifdef OS_WIN32
# include <windows.h>
#endif

extern void _printstr (PRN *prn, const double xx, int *ls);
extern double _gamma_func (double x);
extern double _gammadist (double s1, double s2, double x, int control);

/* ........................................................ */

static int printv (FILE *fp, const int nt, const int v1, 
		   const int *list, double ***pZ)
{
    register int i;
    int v2 = list[0], ls = 0, miss = 0;
    double xx;
    PRN prn;

    prn.fp = fp;
    prn.buf = NULL;
    
    for (i=v1; i<=v2; i++)  {
	xx = (*pZ)[list[i]][nt];
	if (na(xx)) {
	    fprintf(fp, "? ");
	    miss = 1;
	} else 
	    _printstr(&prn, xx, &ls);
    }
    fprintf(fp, "\n");
    return miss;
}

/* ........................................................ */

static void prntdate (const int nt, const int n, 
		      const DATAINFO *pdinfo, PRN *prn)
{
    int t1 = pdinfo->t1, t2 = pdinfo->t2;
    double xx;

    if (n != t2 - t1 + 1) pprintf(prn, "%5d  ", nt);
    else {
	xx = date(nt, pdinfo->pd, pdinfo->sd0);
	if (pdinfo->pd == 1) pprintf(prn, "%5.0f ", xx);
	else {
	    if (pdinfo->pd < 10) pprintf(prn, "%5g ", xx);
	    else pprintf(prn, "%6.2f", xx);
	}
    }
}

/* ........................................................... */

static int get_timevar (DATAINFO *pdinfo, char *timevar)
{
    switch (pdinfo->pd) {
    case 1: 
	strcpy(timevar, "annual"); return 0;
    case 4: 
	strcpy(timevar, "qtrs"); break;
    case 12: 
	strcpy(timevar, "months"); break;
    case 24: 
	strcpy(timevar, "hours"); break;
    default:
	return -1;
    }
    return 1;
}

/* ........................................................ */

int _ztoxy (const int v1, const int v2, double *px, double *py, 
           const DATAINFO *pdinfo, double **Z)
{
    int m = 0, t, t1 = pdinfo->t1, t2 = pdinfo->t2;
    double x1, x2;

    for (t=t1; t<=t2; t++)  {
	x1 = Z[v1][t];
	x2 = Z[v2][t];
	if (na(x1) || na(x2)) continue;
	px[m] = x1;
	py[m++] = x2;
    }
    return m;
}

/* ........................................................ */

static void initpx (const int nn, char *pp)
{
    int i;

    pp[0] = '|';
    for (i=1; i<=nn; i++) pp[i] = ' ';
    pp[nn+1] = '\n';
}

/* ........................................................ */

static void drawline (const int nn, PRN *prn)
{
    int t;

    pprintf(prn, "       |+");
    for (t=1; t<=nn; ++t) {
	if(t%10 == 0) pprintf(prn, "+");
	else pprintf(prn, "-");
    }
    pprintf(prn, "\n");
}

/**
 * plot:
 * @list: contains ID numbers of variables to plot.
 * @Z: data matrix.
 * @pdinfo: data information struct.
 * @oflag: if non-zero, forces two variables to be plotted on the same
 * scale (otherwise they will be scaled to fit).
 * @pause: if non-zero, pause after showing each screen of data.
 * @prn: gretl printing struct.
 *
 * Plot (using ascii graphics) either one or two variables, as given
 * in @list.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int plot (const LIST list, double **Z, const DATAINFO *pdinfo, 
	  int oflag, int pause, PRN *prn)
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

    x = malloc(pdinfo->n * sizeof *x);
    y = malloc(pdinfo->n * sizeof *y);
    if (x == NULL || y == NULL) return E_ALLOC;

    l0 = list[0];
    pprintf(prn, "\n");
    ncols = 70;
    nc2 = ncols/2;
    vy = list[1];
    strcpy(s1, pdinfo->varname[vy]);
    n = 0;

    if (l0 == 1) {
	/* only one variable is to be plotted */
	n = ztox(vy, x, Z, pdinfo);
	_minmax(t1, t2, x, &xmin, &xmax);
	xrange = xmax - xmin;
	cntrline = (floatgt(xmax, 0) && floatlt(xmin, 0))? 1 : 0;
	/* print headings */
	pprintf(prn, _("%25cNOTE: o stands for %s\n\n%8c"), ' ', s1, ' ');
	sprintf(word, "x-min = %g", xmin);
	ls = 8 + strlen(word);
	pprintf(prn, "%s", word);
	sprintf(word, "x-max = %g", xmax);
	ls = 78-ls-strlen(word);
	_bufspace(ls, prn);
	pprintf(prn, "%s\n", word); 
	if (cntrline) {
	    iy = (-xmin/xrange)*ncols;
	    _bufspace(iy+7, prn);
	    pprintf(prn, "0.0\n"); 
	}
	drawline(ncols, prn);
	/* plot values */
	lineno = 3;
	for (t=t1; t<=t2; ++t) {
	    xxx = Z[vy][t];
	    if (na(xxx)) continue;
	    if (pause) page_break(1, &lineno, 0);
	    lineno++;
	    prntdate(t, n, pdinfo, prn);
	    ix = (floatneq(xrange, 0.0))? ((xxx-xmin)/xrange) * ncols : nc2;
	    initpx(ncols, px);
	    if (cntrline) px[iy+1] = '|';
	    px[ix+1] = 'o';
	    for (i=0; i<=ncols+1; i++) 
		pprintf(prn, "%c", px[i]); 
	    if (ix == ncols) pprintf(prn, "\n");
	}
	pprintf(prn, "\n\n");
	free(x);
	free(y);
	return 0;
    }

    /* two variables are to be plotted */
    vz = list[2];
    strcpy(s2, pdinfo->varname[vz]);
    n = _ztoxy(vy, vz, x, y, pdinfo, Z);
    /* find maximum and minimum using all values from both arrays */
    _minmax(t1, t2, x, &xmin, &xmax);
    xrange = xmax - xmin;
    _minmax(t1, t2, y, &ymin, &ymax);
    yrange = ymax - ymin;
    xymin = (xmin <= ymin) ? xmin : ymin;
    xymax = (xmax >= ymax) ? xmax : ymax;
    xyrange = xymax - xymin;
    /* print headings for the two variables */
    pprintf(prn, _("%17cNOTE: o stands for %s,   x stands for %s\n%17c+ means %s "
	   "and %s are equal when scaled\n"), ' ', s1, s2, ' ', s1, s2);
    lineno = 6;
    if (oflag == OPT_O) {
	pprintf(prn, _("%20c%s and %s are plotted on same scale\n\n%8c"),
	       ' ', s1, s2, ' ');
	sprintf(word, "xy-min = %g", xymin);
	ls = 8 + strlen(word);
	pprintf(prn, "%s", word);
	sprintf(word, "xy-max = %g", xymax);
	ls = 78-ls-strlen(word);
	_bufspace(ls, prn);
	pprintf(prn, "%s\n", word);
    }
    else {
	pprintf(prn, "\n");
	sprintf(word, "        o-min = %g", ymin);
	ls = strlen(word);
	pprintf(prn, "%s", word);
	sprintf(word, "o-max = %g", ymax);
	ls = 78-ls-strlen(word);
	_bufspace(ls, prn);
	pprintf(prn, "%s\n", word);
	sprintf(word, "        x-min = %g", xmin);
	ls = strlen(word);
	pprintf(prn, "%s", word);
	sprintf(word, "x-max = %g", xmax);
	ls = 78-ls-strlen(word);
	_bufspace(ls, prn);
	pprintf(prn, "%s\n", word);
    }
    /*  First x and y values are scaled, then it checks to see which scaled
	value is smaller and prints that one first, then prints the larger
	scaled value. If the scaled values are equal then it prints a "+" ,
	otherwise it prints an "x" for the first variable and an "o" for the
	second variable.
    */
    pprintf(prn, "\n");
    cntrline = (floatgt(xymax, 0) && floatlt(xymin, 0))? 1 : 0;
    if (cntrline) {
	iz = (-xymin/xyrange)*ncols;
	_bufspace(iz+7, prn);
	pprintf(prn, "0.0\n");
    }
    drawline(ncols, prn);
    for (t=t1; t<=t2; ++t) {
	if (pause) page_break(1, &lineno, 0);
	lineno++;
	xxx = Z[vy][t];
	yy = Z[vz][t];
	if (na(xxx) || na(yy)) continue;
	prntdate(t, n, pdinfo, prn);
	if (oflag == OPT_O) {
	    ix = (floatneq(xyrange, 0.0))? ((xxx-xymin)/xyrange)*ncols : nc2;
	    iy = (floatneq(xyrange, 0.0))? ((yy-xymin)/xyrange)*ncols : nc2;
	}
	else {
	    ix = (floatneq(xrange, 0.0))? ((xxx-xmin)/xrange)*ncols : nc2;
	    iy = (floatneq(yrange, 0.0))? ((yy-ymin)/yrange)*ncols : nc2;
	}
	initpx(ncols, px);
	if (iz) px[iz+1] = '|';
	if (ix == iy) px[ix+1] = '+';
	else {
	    px[ix+1] = 'o';
	    px[iy+1] = 'x';
	}
	for (i=0; i<=ncols+1; i++) pprintf(prn, "%c", px[i]);
	if (ix == ncols || iy == ncols) pprintf(prn, "\n");
    }
    pprintf(prn, "\n\n");
    free(x);
    free(y);
    return 0;
}

/**
 * graph:
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

int graph (const LIST list, double **Z, const DATAINFO *pdinfo, 
	   const int oflag, PRN *prn)
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

    m = _list_dups(list, GRAPH);
    if (m) {
	fprintf(stderr, _("var no. %d duplicated in command list.\n"), m);
	return 1;
    }

    pprintf(prn, "\n");
    l0 = list[0];
    vy = list[1];

    x = malloc(pdinfo->n * sizeof *x);
    y = malloc(pdinfo->n * sizeof *y);
    uhat = malloc(pdinfo->n * sizeof *uhat);
    if (x == NULL || y == NULL || uhat == NULL) return E_ALLOC;    

    /* Put values from z matrix into x and y arrays */
    m = 0;
    if (l0 == 2) {
	vx = list[2];
	m = _ztoxy(vx, vy, x, y, pdinfo, Z);
	_graphyzx(list, y, uhat, x, m, pdinfo->varname[vy], 
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
	_graphyzx(list, y, uhat, x, -m, pdinfo->varname[vy], 
		  pdinfo->varname[vx], pdinfo, oflag, prn);
    }
    pprintf(prn, "\n");
    free(x); free(y); free(uhat);
    return 0;
}

/* ........................................................ */

static int factorized_vars (double ***pZ, 
			    const int t1, const int t2,
			    double **y1, double **y2,
			    const int ynum, const int dum)
{
    int i = 0, fn, t;
    double xx;

    fn = t2 - t1 + 1;

    *y1 = malloc(fn * sizeof **y1);
    *y2 = malloc(fn * sizeof **y2);
    if (*y1 == NULL || *y2 == NULL) return 1;

    for (t=t1; t<=t2; t++) {
	if (na((*pZ)[ynum][t])) {
	    (*y1)[i] = NADBL;
	    (*y2)[i] = NADBL;
	} else {
	    xx = (*pZ)[dum][t];
	    if (floateq(xx, 1.)) {
		(*y1)[i] = (*pZ)[ynum][t];
		(*y2)[i] = NADBL;
	    } else {
		(*y1)[i] = NADBL;
		(*y2)[i] = (*pZ)[ynum][t];
	    }
	}
	i++;
    }
    return 0;
}

/**
 * gnuplot_init:
 * @ppaths: pointer to path information struct.
 * @fpp: pointer to stream to be opened.
 *
 * If GNUPLOT_PNG is defined and we're in GUI mode, writes a unique
 * temporary filename into the plotfile member of @ppaths; opens
 * plotfile for writing as @fpp; If GNUPLOT_PNG is defined and we're in 
 * GUI mode, writes PNG terminal type into plotfile.
 *
 * Returns: 0 on success, 1 on failure
 */

int gnuplot_init (PATHS *ppaths, FILE **fpp)
{
#ifdef GNUPLOT_PNG
    if (GRETL_GUI(ppaths)) {
	sprintf(ppaths->plotfile, "%sgpttmp.XXXXXX", ppaths->userdir);
	if (mktemp(ppaths->plotfile) == NULL) return 1;
    }
    *fpp = fopen(ppaths->plotfile, "w");
    if (*fpp == NULL) return 1;
    if (GRETL_GUI(ppaths)) {
	fprintf(*fpp, "set term png color\n");
	fprintf(*fpp, "set output 'gretltmp.png'\n");
    }
#else
    *fpp = fopen(ppaths->plotfile, "w");
    if (*fpp == NULL) return 1;
#endif
    return 0;
}

/**
 * gnuplot_display:
 * @ppaths: pointer to path information struct.
 *
 * Executes gnuplot, passing as an argument ppaths->plotfile.
 *
 * Returns: the return value from the system command.
 */

int gnuplot_display (const PATHS *ppaths)
{
    int err = 0;
    char plotcmd[MAXLEN];

#ifdef OS_WIN32
    sprintf(plotcmd, "\"%s\" \"%s\"", ppaths->gnuplot, ppaths->plotfile);
    if (WinExec(plotcmd, SW_SHOWNORMAL) < 32) err = 1;
#else
# ifdef GNUPLOT_PNG
    sprintf(plotcmd, "%s%s \"%s\"", ppaths->gnuplot, 
	    (GRETL_GUI(ppaths))? "" : " -persist", ppaths->plotfile);
# else
    sprintf(plotcmd, "%s -persist \"%s\"", ppaths->gnuplot, ppaths->plotfile);
# endif /* GNUPLOT_PNG */
    if (system(plotcmd)) err = 1;
#endif /* OS_WIN32 */
    return err;
}

/**
 * gnuplot:
 * @list: list of variables to plot, by ID number.
 * @lines: vector of 1s and 0s to indicate whether variables should
 * be represented by lines or not (or NULL).
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 * @ppaths: path information struct.
 * @plot_count: pointer to count of graphs drawn so far.
 * @batch: if non-zero, the plot commands will be saved to file instead
 * of being sent to gnuplot.
 * @gui: should be non-zero if called from GUI client program.
 * @opt:
 *
 * Writes a gnuplot plot file to display the values of the
 * variables in @list and calls gnuplot to make the graph.
 *
 * Returns: 0 on successful completion, -1 if the gnuplot system
 * command fails, or 1 if there are missing data values.
 */

int gnuplot (LIST list, const int *lines, 
	     double ***pZ, DATAINFO *pdinfo,
	     PATHS *ppaths, int *plot_count, 
	     const int batch, const int gui, const int opt)
{
    FILE *fq = NULL;
    int t, t1 = pdinfo->t1, t2 = pdinfo->t2, lo = list[0];
    int i, j, oddman = 0;
    char s1[9], s2[9], xlabel[12], withstring[8];
    int tscale = 0;   /* time series scaling needed? */
    int ts_plot = 1;  /* plotting against time on x-axis? */
    int pdist = 0;    /* plotting probability dist. */
    double a = 0, b = 0;
    double *yvar1 = NULL, *yvar2 = NULL;
    int xvar, miss = 0, ols_ok = 0, tmplist[4];

    if (opt == OPT_M || lines == NULL) {
	strcpy(withstring, "w i");
	pdist = 1;
    }

    if (batch) {  
	*plot_count += 1; 
	sprintf(ppaths->plotfile, "%sgpttmp%02d.plt", ppaths->userdir, 
		*plot_count);
	fq = fopen(ppaths->plotfile, "w");
	if (fq == NULL) return E_FOPEN;
    } else
	if (gnuplot_init(ppaths, &fq)) return E_FOPEN;

    if (strcmp(pdinfo->varname[list[lo]], "time") == 0) {
	if (get_timevar(pdinfo, s2) >= 0) {
	    plotvar(pZ, pdinfo, s2);
	    list[lo] = varindex(pdinfo, s2);
	}
	strcpy(xlabel, _("Observation"));
	if (lo > 2 && lo < 7) tscale = 1;
    } else {
	if (opt == OPT_Z || opt == OPT_RESIDZ)
	    strcpy(xlabel, pdinfo->varname[list[2]]); 
	else
	    strcpy(xlabel, pdinfo->varname[list[lo]]);
	ts_plot = 0;
    }
    if (strcmp(pdinfo->varname[list[lo]], "qtrs") == 0 ||
	strcmp(pdinfo->varname[list[lo]], "months") == 0) {
	ts_plot = 1;
	strcpy(xlabel, _("period"));
    }

    /* add a simple regression line if appropriate */
    if (!pdist && lo == 2 && ts_plot == 0) {
	MODEL plotmod;

	tmplist[0] = 3;
	tmplist[1] = list[1];
	tmplist[2] = list[2];	
	tmplist[3] = 0;	
	_init_model(&plotmod, pdinfo);
	plotmod = lsq(tmplist, pZ, pdinfo, OLS, 0, 0.0);
	if (!plotmod.errcode) {
	    /* is the fit significant? */
	    b = plotmod.coeff[1];
	    if (tprob(b / plotmod.sderr[1], plotmod.dfd) < .10) {
		ols_ok = 1;
		a = plotmod.coeff[2];
	    }
	}
	clear_model(&plotmod, NULL, NULL, pdinfo);
    }

    _adjust_t1t2(NULL, list, &t1, &t2, *pZ, NULL);
    /* if resulting sample range is empty, complain */
    if (t2 == t1) return -999;

    if (opt == OPT_Z || opt == OPT_RESIDZ) { /* separation by dummy variable */
	if (lo != 3) return -1;
	if (factorized_vars(pZ, t1, t2, &yvar1, &yvar2, list[1], list[3])) {
	    fclose(fq);
	    return -1;
	}
    }    

    if (ts_plot) {
	fprintf(fq, "# timeseries %d\n", pdinfo->pd);
	if (pdinfo->pd == 4) {
	    if ((t2 - t1) / 4 < 8) {
		fputs("set xtics nomirror 0,1\n", fq); 
		fputs("set mxtics 4\n", fq);
	    }
	}
	if (pdinfo->pd == 12) {
	    if ((t2 - t1) / 12 < 8) {
		fputs("set xtics nomirror 0,1\n", fq); 
		fputs("set mxtics 12\n", fq);
	    }
	}
    }

    /* titling and so on */
    fprintf(fq, "set xlabel '%s'\n", xlabel);
    fprintf(fq, "set xzeroaxis\n"); 
    fprintf(fq, "set missing \"?\"\n");
    if (lo == 2) {
	if (ols_ok) 
	    fprintf(fq, "set title '%s versus %s (with least squares fit)\n",
		    pdinfo->varname[list[1]], xlabel);
	if (opt == OPT_RESID) {
	    fprintf(fq, "set title 'Regression residuals (= observed - "
		    "fitted %s)'\n", pdinfo->varname[list[1]]);
	    fprintf(fq, "set ylabel 'residual'\n");
	    fprintf(fq, "set nokey\n");
	} else {
	    fprintf(fq, "set ylabel '%s'\n", pdinfo->varname[list[1]]);
	    fprintf(fq, "set nokey\n");
	}
    } else if (opt == OPT_RESIDZ) {
	fprintf(fq, "set title 'Regression residuals (= observed - "
		"fitted %s)'\n", pdinfo->varname[list[1]]);
	fprintf(fq, "set ylabel 'residual'\n");
	fprintf(fq, "set key left top\n");
    } else if (opt == OPT_FA) {
	if (list[3] == pdinfo->v - 1) /* x var is just time or index */
	    fprintf(fq, "set title 'Actual and fitted %s\n",
		    pdinfo->varname[list[2]]);
	else
	    fprintf(fq, "set title 'Actual and fitted %s versus %s\n",
		    pdinfo->varname[list[2]], pdinfo->varname[list[3]]);
	fprintf(fq, "set ylabel '%s'\n", pdinfo->varname[list[2]]);
	fprintf(fq, "set key left top\n");	
    } else
	fprintf(fq, "set key left top\n");

    xvar = (opt == OPT_Z || opt == OPT_RESIDZ)? list[lo - 1] : list[lo];
    if (isdummy(xvar, t1, t2, *pZ)) {
	fputs("set xrange[-1:2]\n", fq);	
	fputs("set xtics (\"0\" 0, \"1\" 1)\n", fq);
    } else {
	double xmin, xmax, xrange;

	_minmax(t1, t2, (*pZ)[xvar], &xmin, &xmax);
	xrange = xmax - xmin;
	xmin -= xrange * .025;
	xmax += xrange * .025;
	fprintf(fq, "set xrange [%g:%g]\n", xmin, xmax);
    }

    if (tscale) { /* two or more vars plotted against time */
	double ymin[6], ymax[6];
	int oddcount;

	/* find minima, maxima of the vars */
	for (i=1; i<lo; i++) 
	    _minmax(t1, t2, (*pZ)[list[i]], &(ymin[i]), &(ymax[i]));
	tscale = 0;
	for (i=1; i<lo; i++) {
	    oddcount = 0;
	    for (j=1; j<lo; j++) {
		if (j == i) continue;
		if (ymax[i] > 5*ymax[j] || ymax[j] > 5*ymax[i]) {
		    tscale = 1;
		    oddcount++;
		}
	    }
	    if (oddcount == lo - 2) {
		oddman = i;
		break;
	    }
	}
    }

    if (tscale) {
	fprintf(fq, "set ytics nomirror\n");
	fprintf(fq, "set y2tics\n");
	fputs("plot \\\n", fq);
	for (i=1; i<lo; i++) {
	    if (i != oddman) 
		fprintf(fq, "'-' using 1:($2) axes x1y1 title '%s (left)' "
			"w lines", 
			pdinfo->varname[list[i]]);
	    else 
		fprintf(fq, "'-' using 1:($2) axes x1y2 title '%s (right)' "
			"w lines", 
			pdinfo->varname[list[i]]);
	    if (i == lo - 1) fprintf(fq, "\n");
	    else fprintf(fq, " , \\\n");
	}
    } else if (opt == OPT_Z || opt == OPT_RESIDZ) { 
	/* FIXME OPT_Z with time series? */
	fputs("plot \\\n", fq);
	if (opt == OPT_Z) 
	    strcpy(s1, pdinfo->varname[list[1]]);
	else 
	    strcpy(s1, "residual");
	strcpy(s2, pdinfo->varname[list[3]]);
	fprintf(fq, " '-' using 1:($2) title '%s (%s=1)', \\\n", s1, s2);
	fprintf(fq, " '-' using 1:($2) title '%s (%s=0)'\n", s1, s2);
    } else {
	fputs("plot \\\n", fq);
	for (i=1; i<lo; i++)  {
	    if (opt == OPT_FA) {
		if (i == 1) strcpy(s1, _("fitted"));
		else strcpy(s1, _("actual"));
	    } else
		strcpy(s1, pdinfo->varname[list[i]]);
	    if (!pdist) { 
		withstring[0] = '\0';
		if ((gui)? lines[i-1] : lines[0]) strcpy(withstring, "w lines");
	    }
	    fprintf(fq, " '-' using 1:($2) title '%s' %s", 
		    s1, withstring);
	    if (i < (lo - 1) || ols_ok) fputs(" , \\\n", fq); 
	    else fputc('\n', fq);
	}
    } 
    if (ols_ok) 
	fprintf(fq, _("%f + %f*x title 'least squares fit' w lines\n"),
		a, b);

    /* supply the data to gnuplot inline 
       -- is something wrong here?? */
    if (opt == OPT_Z || opt == OPT_RESIDZ) {
	double xx, yy;

	for (i=0; i<2; i++) {
	    for (t=t1; t<=t2; t++) {
		xx = (*pZ)[list[2]][t];
		if (na(xx)) continue;
		yy = (i)? yvar2[t-t1] : yvar1[t-t1];
		if (na(yy))
		    fprintf(fq, "%f ?\n", xx);
		else
		    fprintf(fq, "%f %f\n", xx, yy);
	    }
	    fprintf(fq, "e\n");
	}
	free(yvar1);
	free(yvar2);
    } else {
	tmplist[0] = 2;
	tmplist[1] = list[lo];
	for (i=1; i<lo; i++)  {
	    tmplist[2] = list[i];
	    for (t=t1; t<=t2; t++) {
		if (gui && miss == 0) 
		    miss = printv(fq, t, 1, tmplist, pZ); 
		else
		    printv(fq, t, 1, tmplist, pZ);
	    }
	    fprintf(fq, "e\n");
	}
    }

#ifdef OS_WIN32
    fprintf(fq, "pause -1\n");
#endif
    fclose(fq);

    if (!batch) {
	if (gnuplot_display(ppaths)) miss = -1;
    }
    return miss;
}

/**
 * multi_scatters:
 * @list: list of variables to plot, by ID number.
 * @pos: 
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 * @ppaths: path information struct.
 *
 * Writes a gnuplot plot file to display up to 6 small X-Y graphs.
 * variables in @list and calls gnuplot to make the graph.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int multi_scatters (const LIST list, const int pos, double ***pZ, 
		    const DATAINFO *pdinfo, PATHS *ppaths)
{
    int i, t, err = 0, xvar, yvar, *plotlist;
    int nplots, m;
    FILE *fp = NULL;
    double xx;

    if (pos > 2) { /* plot several yvars against one xvar */
	yvar = 0;
	plotlist = malloc(pos * sizeof *plotlist);
	xvar = list[list[0]];
    } else {       /* plot one yvar against several xvars */
	yvar = list[1];
	plotlist = malloc((list[0] + 1 - pos) * sizeof *plotlist);
	xvar = 0;
    }
    if (plotlist == NULL)
	return E_ALLOC;

    if (yvar) {
	plotlist[0] = list[0] - pos;
	for (i=1; i<=plotlist[0]; i++)
	   plotlist[i] = list[i+pos]; 
    } else {
	plotlist[0] = pos - 1;
	for (i=1; i<pos; i++)
	   plotlist[i] = list[i]; 
    }

    /* max 6 plots */
    if (plotlist[0] > 6) plotlist[0] = 6;
    nplots = plotlist[0];

    if (gnuplot_init(ppaths, &fp)) return E_FOPEN;

    fprintf(fp, "# multiple scatterplots\n");
    fprintf(fp, "set size 1.0,1.0\nset origin 0.0,0.0\n"
	    "set multiplot\n");
    fputs("set nokey\n", fp);
    fputs("set noxtics\nset noytics\n", fp);
    for (i=0; i<nplots; i++) {  
	if (nplots <= 4) {
	    fprintf(fp, "set size 0.45,0.5\n");
	    fprintf(fp, "set origin ");
	    if (i == 0) fprintf(fp, "0.0,0.5\n");
	    else if (i == 1) fprintf(fp, "0.5,0.5\n");
	    else if (i == 2) fprintf(fp, "0.0,0.0\n");
	    else if (i == 3) fprintf(fp, "0.5,0.0\n");
	} else {
	    fprintf(fp, "set size 0.31,0.45\n");
	    fprintf(fp, "set origin ");
	    if (i == 0) fprintf(fp, "0.0,0.5\n");
	    else if (i == 1) fprintf(fp, "0.32,0.5\n");
	    else if (i == 2) fprintf(fp, "0.64,0.5\n");
	    else if (i == 3) fprintf(fp, "0.0,0.0\n");
	    else if (i == 4) fprintf(fp, "0.32,0.0\n");
	    else if (i == 5) fprintf(fp, "0.64,0.0\n");
	}
	fprintf(fp, "set xlabel '%s'\n",
		(yvar)? pdinfo->varname[plotlist[i+1]] :
		pdinfo->varname[xvar]);
	fprintf(fp, "set ylabel '%s'\n", 
		(yvar)? pdinfo->varname[yvar] :
		pdinfo->varname[plotlist[i+1]]);
	fprintf(fp, "plot '-' using 1:2\n");
	for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	    m = (yvar)? plotlist[i+1] : xvar;
	    xx = (*pZ)[m][t];
	    if (na(xx)) fprintf(fp, "? ");
	    else fprintf(fp, "%f ", xx);
	    m = (yvar)? yvar : plotlist[i+1];
	    xx = (*pZ)[m][t];
	    if (na(xx)) fprintf(fp, "?\n");
	    else fprintf(fp, "%f\n", xx);
	}
	fprintf(fp, "e\n");

    } 
    fprintf(fp, "set nomultiplot\n");
#ifdef OS_WIN32
    fprintf(fp, "\npause -1\n");
#endif
    fclose(fp);
    err = gnuplot_display(ppaths);
    free(plotlist);
    return err;
}

/**
 * plot_freq:
 * @freq: frequency distribution struct.
 * @ppaths: path information struct.
 * @dist: distribution code (see #dist_codes).
 *
 * Plot the actual frequency distribution for a variable versus a
 * theoretical distribution, Gaussian or gamma (or none).
 *
 * Returns: 0 on successful completion, error code on error.
 */

int plot_freq (FREQDIST *freq, PATHS *ppaths, int dist)
{
    double alpha = 0.0, beta = 0.0, lambda = 1.0;
    FILE *fp = NULL;
    int i, K = freq->numbins;

    if (gnuplot_init(ppaths, &fp)) return E_FOPEN;

    fprintf(fp, "# frequency plot\n");

    if (dist) {
	double propn, plotmin = 0.0, plotmax = 0.0;

	/* find the endpts that straddle the mean... */
	for (i=0; i<K ; i++) 
	    if (freq->endpt[i] > freq->xbar) break;

	/* OK, they are k-1 and k: now find the proportion of the 
	   theoretical distribution they enclose, and calculate a
	   height adjustment factor for the impulses */

	if (dist == NORMAL) {
	    propn = normal((freq->endpt[i-1] - freq->xbar)/freq->sdx) -
		normal((freq->endpt[i] - freq->xbar)/freq->sdx);
	    lambda = 1.0 / (propn * freq->n * sqrt(2 * M_PI) * freq->sdx);
	    fprintf(fp, "sigma = %f\n", freq->sdx);
	    fprintf(fp, "mu = %f\n", freq->xbar);
	    plotmin = freq->xbar - 3.3 * freq->sdx;
	    if (freq->midpt[0] < plotmin) plotmin = freq->midpt[0];
	    plotmax = freq->xbar + 3.3 * freq->sdx;
	    fprintf(fp, _("set label 'Test statistic for normality:'"
		    " at graph .05, graph .9\n"));
	    fprintf(fp, _("set label 'Chi-squared(2) = %.3f, pvalue %.5f'"
		    " at graph .05, graph .85\n"), 
		    freq->chisqu, chisq(freq->chisqu, 2));	
	}
	else if (dist == GAMMA) {
	    double xx, height, var = freq->sdx * freq->sdx;
	
	    /* scale param = variance/mean */
	    beta = var / freq->xbar;
	    /* shape param = mean/scale */
	    alpha = freq->xbar / beta;

	    propn = _gammadist(freq->xbar, var, freq->endpt[i], 2) -
		_gammadist(freq->xbar, var, freq->endpt[i-1], 2);
	    xx = (freq->endpt[i] + freq->endpt[i-1])/2.0;
	    height = pow(xx, alpha - 1.0) * exp(-xx / beta) /
		(_gamma_func(alpha) * pow(beta, alpha));
	    lambda = height/(freq->n * propn);
	    fprintf(fp, "beta = %f\n", beta);
	    fprintf(fp, "alpha = %f\n", alpha);
	    plotmin = 0.0;
	    plotmax = freq->xbar + 4.0 * freq->sdx;
	}

	/* adjust max if needed */
	if (freq->midpt[K-1] > plotmax) plotmax = freq->midpt[K-1];

	fprintf(fp, "set xrange [%.3f:%.3f]\n", plotmin, plotmax);
	fprintf(fp, "set key right top\n");
	fprintf(fp, "plot \\\n");

    } else { /* not dist */
	lambda = 1.0 / freq->n;
	fprintf(fp, "set nokey\n");
	fprintf(fp, "set xlabel 'frequency distribution for %s'\n", 
		freq->varname);	
    }

    /* plot instructions */
    if (!dist) {
	fprintf(fp, "plot '-' using 1:($2) w impulses\n");
    } else if (dist == NORMAL) {
	fprintf(fp, "(1/(sqrt(2*pi)*sigma)*exp(-(x-mu)**2/(2*sigma**2))) "
		"title 'N(%.4f,%.4f)' w lines , \\\n"
		"'-' using 1:($2) title '%s' w impulses\n",
		freq->xbar, freq->sdx, freq->varname);
    }
    else if (dist == GAMMA) {
	fprintf(fp, "x**(alpha-1.0)*exp(-x/beta)/(gamma(alpha)*(beta**alpha)) "
		"title 'gamma(%.4f,%.4f)' w lines , \\\n"
		"'-' using 1:($2) title '%s' w impulses\n",
		alpha, beta, freq->varname); 
    }

    /* send sample data inline */
    for (i=0; i<K; i++) 
	fprintf(fp, "%f %f\n", freq->midpt[i], lambda * freq->f[i]);
    fprintf(fp, "e\n");

#ifdef OS_WIN32
    fprintf(fp, "pause -1\n");
#endif
    if (fp) fclose(fp);
    return gnuplot_display(ppaths);
}

/* ......................................................... */ 

int plot_fcast_errs (const int n, const double *obs, 
		     const double *depvar, const double *yhat, 
		     const double *maxerr, const char *varname, 
		     PATHS *ppaths)
{
    FILE *fp = NULL;
    int t;

    if (gnuplot_init(ppaths, &fp)) return E_FOPEN;

    fprintf(fp, "# forecasts with 95 pc conf. interval\n");
    fprintf(fp, "set key left top\n"
	    "plot \\\n'-' using 1:2 title '%s' w lines , \\\n"
	    "'-' using 1:2 title 'fitted' w lines , \\\n"
	    "'-' using 1:2:3 title '95%% confidence interval' "
	    "w errorbars\n", varname);
    /* send data inline */
    for (t=0; t<n; t++)
	fprintf(fp, "%f %f\n", obs[t], depvar[t]);
    fprintf(fp, "e\n");
    for (t=0; t<n; t++)
	fprintf(fp, "%f %f\n", obs[t], yhat[t]);
    fprintf(fp, "e\n");
    for (t=0; t<n; t++)
	fprintf(fp, "%f %f %f\n", obs[t], yhat[t], maxerr[t]);
    fprintf(fp, "e\n");

#ifdef OS_WIN32
    fprintf(fp, "pause -1\n");
#endif
    fclose(fp);
    return gnuplot_display(ppaths);
}

/* ........................................................... */

void free_plot (GPT_SPEC *plot)
{
    int i;

    if (plot->lines) free(plot->lines);
    if (plot->data) free(plot->data);
    if (plot->literal[0]) {
	for (i=0; i<4; i++)
	    free(plot->literal[i]);
    }
}

/* ........................................................... */

int open_gnuplot_pipe (const PATHS *ppaths, GPT_SPEC *plot)
     /* add file or pipe to plot struct */
{
    FILE *fp;
#ifndef OS_WIN32 
    char gnuplot_pipe[MAXLEN]; 
#endif

#ifdef OS_WIN32
    fp = fopen(ppaths->plotfile, "w");
#else
    sprintf(gnuplot_pipe, "gnuplot"); /* should be user-specified name? */
    fp = popen(gnuplot_pipe, "w");
#endif
    plot->edit = 1;
    if (fp == NULL) return 1;
    plot->fp = fp;
    return 0;
}

/* ........................................................... */

int termtype_to_termstr (char *termtype, char *termstr)
{
    int cmds = 0;

    if (!strcmp(termtype, "postscript")) 
	strcpy(termstr, "postscript eps");
    else if (!strcmp(termtype, "fig")) 
	strcpy(termstr, "fig");
    else if (!strcmp(termtype, "latex")) 
	strcpy(termstr, "latex");
    else if (!strcmp(termtype, "png")) 
	strcpy(termstr, "png small color");
    else if (!strcmp(termtype, "plot commands")) 
	cmds = 1;
    else strcpy(termstr, termtype);
    return cmds;
}

/* ........................................................... */

int go_gnuplot (GPT_SPEC *plot, char *fname, PATHS *ppaths)
     /* ship out a plot struct, to gnuplot or file.  
	N.B. under unix plot->fp will be a pipe */
{
    FILE *fp;
    int i, k, t, dump = 0, plotn, lo = plot->list[0];
    int err = 0, miss, datlines;
    char termstr[72];
    double xx;

    dump = termtype_to_termstr(plot->termtype, termstr);
    if (dump) {  /* dump of gnuplot commands to named file */
	if (fname == NULL) return 1;  /* impossible */
	fp = fopen(fname, "w");
	if (fp == NULL) return 1;
    } else {     /* output to gnuplot, for screen or other "term" */
	fp = plot->fp;
	if (fname != NULL) { /* not a screen display */
	    fprintf(fp, "set term %s\n", termstr);
	    fprintf(fp, "set output '%s'\n", fname);
	}
    }

    fprintf(fp, "set title '%s'\n", plot->titles[0]);
    fprintf(fp, "set xlabel '%s'\n", plot->titles[1]);
    fprintf(fp, "set ylabel '%s'\n", plot->titles[2]);
    if (plot->y2axis)
	fprintf(fp, "set y2label '%s'\n", plot->titles[3]);

    fprintf(fp, "set xzeroaxis\n");
    fprintf(fp, "set missing \"?\"\n");
    if (strcmp(plot->keyspec, "none") == 0)
	fprintf(fp, "set nokey\n");
    else
	fprintf(fp, "set key %s\n", plot->keyspec);

    k = (plot->y2axis)? 3: 2;
    for (i=0; i<k; i++) {
	fprintf(fp, "set %srange [%s:%s]\n",
		(i==0)? "x" : (i==1)? "y" : "y2",
		plot->range[i][0], plot->range[i][1]);
    }

    /* customized xtics? */
    if (strlen(plot->xtics))
	fprintf(fp, "set xtics %s\n", plot->xtics);
    if (strlen(plot->mxtics))
	fprintf(fp, "set mxtics %s\n", plot->mxtics);

    /* using two y axes? */
    if (plot->y2axis) {
	fprintf(fp, "set ytics nomirror\n");
	fprintf(fp, "set y2tics\n");
    }

    if (plot->code == FREQ || plot->code == PERGM) { 
	if (plot->code == FREQ)
	    fprintf(fp, "# frequency plot\n");
	else
	    fprintf(fp, "# periodogram\n");
	for (i=0; i<4; i++)
	    fprintf(fp, "%s\n", plot->literal[i]);
    } 
    fputs("plot \\\n", fp);

    datlines = lo - 1;
    for (i=1; i<lo; i++) {
	if (strcmp(plot->lines[i-1].scale, "NA")) {
	    fprintf(fp, "'-' using 1:($2*%s) ", 
		    plot->lines[i-1].scale);
	} else {
	    fprintf(fp, "%s ", plot->lines[i-1].formula); 
	    datlines--;
	}
	fprintf(fp, "axes x1y%d title '%s' w %s", 
		plot->lines[i-1].yaxis,
		plot->lines[i-1].title,
		plot->lines[i-1].style);
	if (i == lo - 1) fprintf(fp, "\n");
	else fprintf(fp, ", \\\n");
    } 

    /* supply the data to gnuplot inline */
    miss = 0;
    plotn = plot->t2 - plot->t1 + 1;
    for (i=1; i<=datlines; i++) {  
	for (t=plot->t1; t<=plot->t2; t++) {
	    xx = plot->data[t - plot->t1];
	    if (na(xx)) {
		fprintf(fp, "? ");
		miss = 1;
	    } else fprintf(fp, "%f ", xx);
	    xx = plot->data[plotn * i + t - plot->t1];
	    if (na(xx)) {
		fprintf(fp, "?\n");
		miss = 1;
	    } else fprintf(fp, "%f\n", xx);
	}
	fprintf(fp, "e\n");
    }

    fflush(fp);
    if (dump) fclose(fp);
    
#ifdef OS_WIN32
    if (!dump) {
	char plotcmd[MAXLEN];

	fprintf(fp, "pause -1\n");
	fclose(fp);
	sprintf(plotcmd, "\"%s\" < \"%s\"", ppaths->pgnuplot, ppaths->plotfile); 
	if (system(plotcmd)) err = 1;
    }
#endif
    if (miss) err = 2;
    return err;
}
