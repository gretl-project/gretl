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

#include "libgretl.h"
#include "internal.h"
#include "../../cephes/libprob.h"

#include <unistd.h>

#ifdef OS_WIN32
# include <windows.h>
#endif

/* pvalues.c */
extern double gamma_dist (double s1, double s2, double x, int control);

static char gnuplot_path[MAXLEN];

static const char *auto_ols_string = "# plot includes automatic OLS line\n";

static const char *get_gnuplot_pallette (int i, int plottype);

/* ........................................................ */

static int printvars (FILE *fp, int t, const int *list, double **Z,
		      const char *label, double offset)
{
    int i, miss = 0;
    double xx;

    for (i=1; i<=list[0]; i++)  {
	xx = Z[list[i]][t];
	if (na(xx)) {
	    fputs("? ", fp);
	    miss = 1;
	} else {
	    if (i == 1) { /* the x variable */
		xx += offset;
	    }
	    fprintf(fp, "%.8g ", xx);
	}
    }

    if (label != NULL) {
	fprintf(fp, "# %s", label);
    }

    fputc('\n', fp);

    return miss;
}

/* ........................................................ */

static void prntdate (int nt, int n, 
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
	strcpy(timevar, "hrs"); break;
    default:
	return -1;
    }
    return 1;
}

/* ........................................................ */

int _ztoxy (int v1, int v2, double *px, double *py, 
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

static void initpx (int nn, char *pp)
{
    int i;

    pp[0] = '|';
    for (i=1; i<=nn; i++) pp[i] = ' ';
    pp[nn+1] = '\n';
}

/* ........................................................ */

static void drawline (int nn, PRN *prn)
{
    int t;

    pputs(prn, "       |+");
    for (t=1; t<=nn; ++t) {
	if(t%10 == 0) pputs(prn, "+");
	else pputs(prn, "-");
    }
    pputs(prn, "\n");
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
	  unsigned char oflag, int pause, PRN *prn)
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
    pputs(prn, "\n");
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
	pputs(prn, word);
	sprintf(word, "x-max = %g", xmax);
	ls = 78-ls-strlen(word);
	_bufspace(ls, prn);
	pprintf(prn, "%s\n", word); 
	if (cntrline) {
	    iy = (-xmin/xrange)*ncols;
	    _bufspace(iy+7, prn);
	    pputs(prn, "0.0\n"); 
	}
	drawline(ncols, prn);
	/* plot values */
	lineno = 3;
#ifdef ENABLE_NLS
	setlocale(LC_NUMERIC, "C");
#endif
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
	    if (ix == ncols) pputs(prn, "\n");
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
    if (oflag == 'o') {
	pprintf(prn, _("%20c%s and %s are plotted on same scale\n\n%8c"),
	       ' ', s1, s2, ' ');
	sprintf(word, "xy-min = %g", xymin);
	ls = 8 + strlen(word);
	pputs(prn, word);
	sprintf(word, "xy-max = %g", xymax);
	ls = 78-ls-strlen(word);
	_bufspace(ls, prn);
	pprintf(prn, "%s\n", word);
    }
    else {
	pputs(prn, "\n");
	sprintf(word, "        o-min = %g", ymin);
	ls = strlen(word);
	pputs(prn, word);
	sprintf(word, "o-max = %g", ymax);
	ls = 78-ls-strlen(word);
	_bufspace(ls, prn);
	pprintf(prn, "%s\n", word);
	sprintf(word, "        x-min = %g", xmin);
	ls = strlen(word);
	pputs(prn, word);
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
    pputs(prn, "\n");
    cntrline = (floatgt(xymax, 0) && floatlt(xymin, 0))? 1 : 0;
    if (cntrline) {
	iz = (-xymin/xyrange)*ncols;
	_bufspace(iz+7, prn);
	pputs(prn, "0.0\n");
    }
    drawline(ncols, prn);
#ifdef ENABLE_NLS
	setlocale(LC_NUMERIC, "C");
#endif
    for (t=t1; t<=t2; ++t) {
	if (pause) page_break(1, &lineno, 0);
	lineno++;
	xxx = Z[vy][t];
	yy = Z[vz][t];
	if (na(xxx) || na(yy)) continue;
	prntdate(t, n, pdinfo, prn);
	if (oflag == 'o') {
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
	if (ix == ncols || iy == ncols) pputs(prn, "\n");
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
	   unsigned char oflag, PRN *prn)
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

    pputs(prn, "\n");
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
    pputs(prn, "\n");
    free(x); free(y); free(uhat);
    return 0;
}

/* ........................................................ */

static int factorized_vars (double ***pZ, 
			    int t1, int t2,
			    double **y1, double **y2,
			    int ynum, int dum)
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

#ifdef OS_WIN32

int gnuplot_has_ttf (void)
{
    /* we know the gnuplot supplied with gretl for win32
       does TrueType fonts */
    return 1;
}

static int gnuplot_has_filledcurve (void)
{
    /* and we know it does filledcurve */
    return 1;
}

#else

static int old_gnuplot_png (void)
{
    /* "color" is wanted for gnuplot 3.7, but not 3.8 */
    static int c = -1; 

    if (c == -1) {
	char cmd[512];

	sprintf(cmd, "echo \"set term png color\" | %s 2>/dev/null",
		(*gnuplot_path == 0)? "gnuplot" : gnuplot_path);
	c = system(cmd);
    }
    return !c;
}

int gnuplot_has_ttf (void)
{
    static int t = -1; 

    if (t == -1) {
	char cmd[512];

	sprintf(cmd, "echo \"set term png font arial 8\" | %s 2>/dev/null",
		(*gnuplot_path == 0)? "gnuplot" : gnuplot_path);
	t = system(cmd);
    }
    return !t;
}

static int gnuplot_has_filledcurve (void)
{
    static int t = -1; 

    if (t == -1) {
	char cmd[512];

	sprintf(cmd, "echo \"set data style filledcurve\" | %s 2>/dev/null",
		(*gnuplot_path == 0)? "gnuplot" : gnuplot_path);
	t = system(cmd);
    }
    return !t;
}

#endif /* OS_WIN32 */

/**
 * get_gretl_png_term_line:
 *
 * Constructs a suitable line for sending to gnuplot to invoke
 * the PNG "terminal".  Checks the environment for setting of 
 * GRETL_PNG_GRAPH_FONT. Also appends a color-specification string 
 * if the gnuplot PNG driver supports this.
 *
 * Returns: the term string, "set term png ..."
 */

const char *get_gretl_png_term_line (const PATHS *ppaths, int plottype)
{
    static char png_term_line[256];
    char font_string[128];
    char color_string[64];
    int oldgp = 0, gpttf = 1;
    const char *grfont = NULL;

    *font_string = 0;
    *color_string = 0;
#ifndef OS_WIN32
    oldgp = old_gnuplot_png();
    gpttf = gnuplot_has_ttf();
#endif

    /* plot font setup */
    if (gpttf) {
	if (*ppaths->pngfont != 0) {
	    grfont = ppaths->pngfont;
	} else {
	    grfont = getenv("GRETL_PNG_GRAPH_FONT");
	}
	if (grfont != NULL && *grfont != 0) {
	    sprintf(font_string, " font %s", grfont);
	}
    }

    /* plot color setup */
    if (oldgp) {
	strcpy(color_string, " color");
    } else {
	int i;

	strcpy(color_string, " xffffff x000000 x202020");
	for (i=0; i<3; i++) {
	    strcat(color_string, " ");
	    strcat(color_string, get_gnuplot_pallette(i, plottype));
	}
    }	

    sprintf(png_term_line, "set term png%s%s",
	    font_string, color_string);

#ifdef GNUPLOT_DEBUG
    fprintf(stderr, "png term line:\n'%s'\n", png_term_line);
#endif

    return png_term_line;
}

/**
 * gnuplot_init:
 * @ppaths: pointer to path information struct.
 * @plottype: code for the type of plot.
 * @fpp: pointer to stream to be opened.
 *
 * If GNUPLOT_PNG is defined and we're in GUI mode, writes a unique
 * temporary filename into the plotfile member of @ppaths; opens
 * plotfile for writing as @fpp; If GNUPLOT_PNG is defined and we're in 
 * GUI mode, writes PNG terminal type into plotfile.
 *
 * Returns: 0 on success, 1 on failure
 */

int gnuplot_init (PATHS *ppaths, int plottype, FILE **fpp)
{
#ifdef GNUPLOT_PNG
    if (*gnuplot_path == 0) {
	strcpy(gnuplot_path, ppaths->gnuplot);
    }

    if (GRETL_GUI(ppaths)) {
	sprintf(ppaths->plotfile, "%sgpttmp.XXXXXX", ppaths->userdir);
	if (mktemp(ppaths->plotfile) == NULL) return 1;
    } else {
	sprintf(ppaths->plotfile, "%sgpttmp.plt", ppaths->userdir);
    }

    *fpp = fopen(ppaths->plotfile, "w");
    if (*fpp == NULL) return 1;

    if (GRETL_GUI(ppaths)) {
	fprintf(*fpp, "%s\n", get_gretl_png_term_line(ppaths, plottype));
	fprintf(*fpp, "set output '%sgretltmp.png'\n", ppaths->userdir);
    }
#else /* not GNUPLOT_PNG */
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

#ifdef GNUPLOT_PNG

int gnuplot_display (const PATHS *ppaths)
{
    int err = 0;
    char plotcmd[MAXLEN];

# ifdef OS_WIN32
    sprintf(plotcmd, "\"%s\" \"%s\"", ppaths->gnuplot, ppaths->plotfile);
    err = winfork(plotcmd, NULL, SW_SHOWMINIMIZED, 0);
# else
    sprintf(plotcmd, "%s%s \"%s\"", ppaths->gnuplot, 
	    (GRETL_GUI(ppaths))? "" : " -persist", ppaths->plotfile);
    if (system(plotcmd)) err = 1;
# endif /* OS_WIN32 */

    return err;
}

#else /* following: old non-GNUPLOT_PNG version */

int gnuplot_display (const PATHS *ppaths)
{
    int err = 0;
    char plotcmd[MAXLEN];

# ifdef OS_WIN32
    sprintf(plotcmd, "\"%s\" \"%s\"", ppaths->gnuplot, ppaths->plotfile);
    if (WinExec(plotcmd, SW_SHOWNORMAL) < 32) err = 1;
# else
    sprintf(plotcmd, "%s -persist \"%s\"", ppaths->gnuplot, ppaths->plotfile);
    if (system(plotcmd)) err = 1;
# endif /* OS_WIN32 */

    return err;
}

#endif /* GNUPLOT_PNG alternation */

enum {
    GTITLE_VLS,
    GTITLE_RESID,
    GTITLE_AF,
    GTITLE_AFV
} graph_titles;

static void make_gtitle (FILE *fp, int code, const char *n1, const char *n2)
{
    char title[128];

    switch (code) {
    case GTITLE_VLS:
	sprintf(title, I_("%s versus %s (with least squares fit)"),
		n1, n2);
	break;
    case GTITLE_RESID:
	sprintf(title, I_("Regression residuals (= observed - fitted %s)"), n1);
	break;
    case GTITLE_AF:
	sprintf(title, I_("Actual and fitted %s"), n1);
	break;
    case GTITLE_AFV:
	sprintf(title, I_("Actual and fitted %s versus %s"), n1, n2);
	break;
    default:
	*title = 0;
	break;
    }

    if (*title) fprintf(fp, "set title '%s'\n", title);
}

static const char *front_strip (const char *s)
{
    while (*s) {
	if (isspace(*s) || *s == '{') s++;
	else break;
    }	
    return s;
}

static void line_out (const char *s, int len, FILE *fp)
{
    char *p = malloc(len + 1);

    if (p != NULL) {
	*p = 0;
	strncat(p, s, len);
	fprintf(fp, "%s\n", front_strip(p));
	free(p);
    }
}

static void print_gnuplot_literal_lines (const char *s, FILE *fp)
{
    const char *p;

    p = s = front_strip(s);

    while (*s && *s != '}') {
	if (*s == ';') {
	    line_out(p, s - p, fp);
	    p = s + 1;
	}
	s++;
    }
}

static const char *get_series_name (const DATAINFO *pdinfo, int v)
{
    if (pdinfo->varinfo != NULL && pdinfo->varinfo[v] != NULL
	&& *DISPLAYNAME(pdinfo, v)) {
	return DISPLAYNAME(pdinfo, v);
    } else {
	return pdinfo->varname[v];
    }
}

/**
 * gnuplot:
 * @list: list of variables to plot, by ID number.
 * @lines: vector of 1s and 0s to indicate whether variables should
 * be represented by lines or not (or NULL).
 * @literal: commands to be passed to gnuplot.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 * @ppaths: path information struct.
 * @plot_count: pointer to count of graphs drawn so far.
 * @flags: bitwise OR of zero or more options from #gnuplot_flags
 *
 * Writes a gnuplot plot file to display the values of the
 * variables in @list and calls gnuplot to make the graph.
 *
 * Returns: 0 on successful completion, -1 if the gnuplot system
 * command fails, or 1 if there are missing data values.
 */

int gnuplot (LIST list, const int *lines, const char *literal,
	     double ***pZ, DATAINFO *pdinfo, PATHS *ppaths, 
	     int *plot_count, unsigned char flags)
{
    FILE *fq = NULL;
    int t, t1 = pdinfo->t1, t2 = pdinfo->t2, lo = list[0];
    int i, j, oddman = 0;
    char s1[MAXDISP], s2[MAXDISP], xlabel[MAXDISP], withstring[8];
    char depvar[9];
    int yscale = 0;   /* two y axis scales needed? */
    int ts_plot = 1;  /* plotting against time on x-axis? */
    int pdist = 0;    /* plotting probability dist. */
    double a = 0, b = 0, offset = 0, xrange = 0;
    double *yvar1 = NULL, *yvar2 = NULL;
    int xvar, miss = 0, ols_ok = 0, tmplist[4];
    int npoints;

    *withstring = 0;

    if ((flags & GP_IMPULSES) || lines == NULL) {
	if (flags & ~GP_OLS_OMIT) {
	    strcpy(withstring, "w i");
	}
	pdist = 1;
    }

    *depvar = 0;
    if (flags & GP_RESIDS) {
	/* a hack to get the name of the dependent variable into
	   the graph */
	strcpy(depvar, pdinfo->varname[list[lo]]);
	lo--;
    }

    /* organize the output */
    if ((flags & GP_FILE) && *ppaths->plotfile) {
	fq = fopen(ppaths->plotfile, "w");
    } else if ((flags & GP_BATCH) && plot_count != NULL) {  
	if (*ppaths->plotfile == 0) {
	    *plot_count += 1; 
	    sprintf(ppaths->plotfile, "%sgpttmp%02d.plt", ppaths->userdir, 
		    *plot_count);
	}
	fq = fopen(ppaths->plotfile, "w");
	if (fq == NULL) return E_FOPEN;
    } else {
	if (gnuplot_init(ppaths, PLOT_REGULAR, &fq)) return E_FOPEN;
    }

    /* setting yscale at this point is not a commitment, we will
       do some further tests below */
    if (lo > 2 && lo < 7 && 
	(flags & ~GP_RESIDS) && (flags & ~GP_FA)
	&& (flags & ~GP_DUMMY)) {
	yscale = 1;
    }

    if (strcmp(pdinfo->varname[list[lo]], "time") == 0) {
	if (get_timevar(pdinfo, s2) >= 0) {
	    int pv = plotvar(pZ, pdinfo, s2);

	    if (pv > 0) list[lo] = pv;
	    else {
		if (fq != NULL) fclose(fq);
		return E_ALLOC;
	    }
	}
	/* strcpy(xlabel, I_("Observation")); */
	*xlabel = 0;
    } else {
	if (flags & GP_DUMMY) {
	    strcpy(xlabel, get_series_name(pdinfo, list[2])); 
	} else {
	    strcpy(xlabel, get_series_name(pdinfo, list[lo]));
	}
	ts_plot = 0;
    }

    if (strcmp(pdinfo->varname[list[lo]], "qtrs") == 0 ||
	strcmp(pdinfo->varname[list[lo]], "months") == 0) {
	ts_plot = 1;
	*xlabel = 0;
	/* strcpy(xlabel, I_("period")); */
    }

    /* add a simple regression line if appropriate */
    if (!pdist && !(flags & GP_OLS_OMIT) && lo == 2 && ts_plot == 0) {
	MODEL plotmod;

	tmplist[0] = 3;
	tmplist[1] = list[1];
	tmplist[2] = 0;
	tmplist[3] = list[2];	
	_init_model(&plotmod, pdinfo);
	plotmod = lsq(tmplist, pZ, pdinfo, OLS, 0, 0.0);
	if (!plotmod.errcode) {
	    /* is the fit significant? */
	    b = plotmod.coeff[1];
	    if (tprob(b / plotmod.sderr[1], plotmod.dfd) < .10) {
		ols_ok = 1;
		a = plotmod.coeff[0];
	    }
	}
	clear_model(&plotmod, pdinfo);
    }

    _adjust_t1t2(NULL, list, &t1, &t2, *pZ, NULL);
    /* if resulting sample range is empty, complain */
    if (t2 == t1) return -999;
    npoints = t2 - t1 + 1;

    if (flags & GP_DUMMY) { /* separation by dummy variable */
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

    fprintf(fq, "set xlabel '%s'\n", xlabel);
    fputs("set xzeroaxis\n", fq); 
    fputs("set missing \"?\"\n", fq);

    if (lo == 2) {
	/* only two variables */
	if (ols_ok) {
	    fputs(auto_ols_string, fq);
	    if (flags & GP_FA) {
		make_gtitle(fq, GTITLE_AFV, get_series_name(pdinfo, list[1]), 
			    get_series_name(pdinfo, list[2]));
	    } else {
		make_gtitle(fq, GTITLE_VLS, get_series_name(pdinfo, list[1]), 
			    xlabel);
	    }
	}
	if (flags & GP_RESIDS && !(flags & GP_DUMMY)) { 
	    make_gtitle(fq, GTITLE_RESID, depvar, NULL);
	    fprintf(fq, "set ylabel '%s'\n", I_("residual"));
	    fputs("set nokey\n", fq);
	} else {
	    fprintf(fq, "set ylabel '%s'\n", get_series_name(pdinfo, list[1]));
	    fputs("set nokey\n", fq);
	}
    } else if ((flags & GP_RESIDS) && (flags & GP_DUMMY)) { 
	make_gtitle(fq, GTITLE_RESID, depvar, NULL);
	fprintf(fq, "set ylabel '%s'\n", I_("residual"));
	fputs("set key left top height 1 width 1 box\n", fq);
    } else if (flags & GP_FA) {
	if (list[3] == pdinfo->v - 1) { 
	    /* x var is just time or index: is this always right? */
	    make_gtitle(fq, GTITLE_AF, get_series_name(pdinfo, list[2]), NULL);
	} else {
	    make_gtitle(fq, GTITLE_AFV, get_series_name(pdinfo, list[2]), 
			get_series_name(pdinfo, list[3]));
	}
	fprintf(fq, "set ylabel '%s'\n", get_series_name(pdinfo, list[2]));
	fputs("set key left top height 1 width 1 box\n", fq);	
    } else {
	fputs("set key left top height 1 width 1 box\n", fq);
    }

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "C");
#endif

    xvar = (flags & GP_DUMMY)? list[lo - 1] : list[lo];

    if (isdummy((*pZ)[xvar], t1, t2)) {
	fputs("set xrange[-1:2]\n", fq);	
	fputs("set xtics (\"0\" 0, \"1\" 1)\n", fq);
	xrange = 3;
    } else {
	double xmin, xmax;

	_minmax(t1, t2, (*pZ)[xvar], &xmin, &xmax);
	xrange = xmax - xmin;
	xmin -= xrange * .025;
	xmax += xrange * .025;
	fprintf(fq, "set xrange [%.8g:%.8g]\n", xmin, xmax);
	xrange = xmax - xmin;
    }

    if (yscale) { 
	/* two or more y vars plotted against some x: test to
	   see if we want to use two y axes */
	double ymin[6], ymax[6];
	double ratio;
	int oddcount;

	/* find minima, maxima of the vars */
	for (i=1; i<lo; i++) {
	    _minmax(t1, t2, (*pZ)[list[i]], &(ymin[i]), &(ymax[i]));
	}
	yscale = 0;
	for (i=1; i<lo; i++) {
	    oddcount = 0;
	    for (j=1; j<lo; j++) {
		if (j == i) continue;
		ratio = ymax[i] / ymax[j];
		if (ratio > 5.0 || ratio < 0.2) {
		    yscale = 1;
		    oddcount++;
		}
	    }
	    if (oddcount == lo - 2) {
		oddman = i;
		break;
	    }
	}
    }

    if (yscale) {
	fputs("set ytics nomirror\n", fq);
	fputs("set y2tics\n", fq);
    }

    if (literal != NULL && *literal != 0) {
	print_gnuplot_literal_lines(literal, fq);
    }

    if (yscale) {
	fputs("plot \\\n", fq);
	for (i=1; i<lo; i++) {
	    if (i != oddman) {
		fprintf(fq, "'-' using 1:($2) axes x1y1 title '%s (%s)' %s",
			get_series_name(pdinfo, list[i]), I_("left"),
			(pdist)? "w impulses" : 
			(ts_plot)? "w lines" : "w points");
	    } else {
		fprintf(fq, "'-' using 1:($2) axes x1y2 title '%s (%s)' %s",
			get_series_name(pdinfo, list[i]), I_("right"),
			(pdist)? "w impulses" : 
			(ts_plot)? "w lines" : "w points");
	    }
	    if (i == lo - 1) fputc('\n', fq);
	    else fputs(" , \\\n", fq);
	}
    } else if (flags & GP_DUMMY) { 
	fputs("plot \\\n", fq);
	if (!(flags & GP_RESIDS)) {
	    strcpy(s1, get_series_name(pdinfo, list[1]));
	} else {
	    strcpy(s1, I_("residual"));
	}
	strcpy(s2, get_series_name(pdinfo, list[3]));
	fprintf(fq, " '-' using 1:($2) title '%s (%s=1)', \\\n", s1, s2);
	fprintf(fq, " '-' using 1:($2) title '%s (%s=0)'\n", s1, s2);
    } else {
	fputs("plot \\\n", fq);
	for (i=1; i<lo; i++)  {
	    if (flags & GP_FA) {
		if (i == 1) strcpy(s1, I_("fitted"));
		else strcpy(s1, I_("actual"));
	    } else {
		strcpy(s1, get_series_name(pdinfo, list[i]));
	    }
	    if (!pdist) { 
		if (flags & GP_GUI) {
		    if (lines[i-1]) {
			strcpy(withstring, "w lines");
		    } else {
			strcpy(withstring, "w points");
		    }
		} else {
		    if (lines[0]) {
			strcpy(withstring, "w lines");
		    } else {
			strcpy(withstring, "w points");
		    }
		}
	    }
	    fprintf(fq, " '-' using 1:($2) title '%s' %s", 
		    s1, withstring);
	    if (i < (lo - 1) || ols_ok) fputs(" , \\\n", fq); 
	    else fputc('\n', fq);
	}
    } 

    if (ols_ok) {
	fprintf(fq, "%g + %g*x title '%s' w lines\n", a, b, 
		I_("least squares fit"));
    }

    /* multi impulse plot? calculate offset for lines */
    if ((flags & GP_IMPULSES) && list[lo] > 2) {
	offset = 0.10 * xrange / npoints;
    }

    /* supply the data to gnuplot inline */
    if (flags & GP_DUMMY) {
	double xx, yy;

	for (i=0; i<2; i++) {
	    for (t=t1; t<=t2; t++) {
		xx = (*pZ)[list[2]][t];
		if (na(xx)) continue;
		yy = (i)? yvar2[t-t1] : yvar1[t-t1];
		if (na(yy)) {
		    fprintf(fq, "%.8g ?\n", xx);
		} else {
		    fprintf(fq, "%.8g %.8g", xx, yy);
		    if (!ts_plot && pdinfo->markers) {
			fprintf(fq, " # %s", pdinfo->S[t]);
		    }
		    fputc('\n', fq);
		}
	    }
	    fputs("e\n", fq);
	}
	free(yvar1);
	free(yvar2);
    } else {
	tmplist[0] = 2;
	tmplist[1] = list[lo];
	for (i=1; i<lo; i++) {
	    double xoff = offset * (i - 1);

	    tmplist[2] = list[i];
	    for (t=t1; t<=t2; t++) {
		int t_miss;
		const char *label = NULL;

		if (!ts_plot && pdinfo->markers && i == 1) {
		    label = pdinfo->S[t];
		}
		t_miss = printvars(fq, t, tmplist, *pZ, label, xoff);
		if ((flags & GP_GUI) && miss == 0) {
		    miss = t_miss;
		}
	    }
	    fputs("e\n", fq);
	}
    }

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif

#if defined(OS_WIN32) && !defined(GNUPLOT_PNG)
    fputs("pause -1\n", fq);
#endif
    fclose(fq);

    if (!(flags & GP_BATCH)) {
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

int multi_scatters (const LIST list, int pos, double ***pZ, 
		    const DATAINFO *pdinfo, PATHS *ppaths,
		    int *plot_count, unsigned char flags)
{
    int i, t, err = 0, xvar, yvar, *plotlist;
    int nplots, m;
    FILE *fp = NULL;
    double xx;

    if (pos > 2) { 
	/* plot several yvars against one xvar */
	yvar = 0;
	plotlist = malloc(pos * sizeof *plotlist);
	xvar = list[list[0]];
    } else {       
	/* plot one yvar against several xvars */
	yvar = list[1];
	plotlist = malloc((list[0] + 1 - pos) * sizeof *plotlist);
	xvar = 0;
    }

    if (plotlist == NULL) return E_ALLOC;

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

    /* organize the output */
    if ((flags & GP_FILE) && *ppaths->plotfile) {
	fp = fopen(ppaths->plotfile, "w");
    } else if ((flags & GP_BATCH) && plot_count != NULL) {  
	if (*ppaths->plotfile == 0) {
	    *plot_count += 1; 
	    sprintf(ppaths->plotfile, "%sgpttmp%02d.plt", ppaths->userdir, 
		    *plot_count);
	}
	fp = fopen(ppaths->plotfile, "w");
	if (fp == NULL) return E_FOPEN;
    } else {
	if (gnuplot_init(ppaths, PLOT_MULTI_SCATTER, &fp)) return E_FOPEN;
    }

    fputs("# multiple scatterplots\n", fp);
    fputs("set size 1.0,1.0\nset origin 0.0,0.0\n"
	  "set multiplot\n", fp);
    fputs("set nokey\n", fp);
    fputs("set noxtics\nset noytics\n", fp);
    for (i=0; i<nplots; i++) {  
	if (nplots <= 4) {
	    fputs("set size 0.45,0.5\n", fp);
	    fputs("set origin ", fp);
	    if (i == 0) fputs("0.0,0.5\n", fp);
	    else if (i == 1) fputs("0.5,0.5\n", fp);
	    else if (i == 2) fputs("0.0,0.0\n", fp);
	    else if (i == 3) fputs("0.5,0.0\n", fp);
	} else {
	    fputs("set size 0.31,0.45\n", fp);
	    fputs("set origin ", fp);
	    if (i == 0) fputs("0.0,0.5\n", fp);
	    else if (i == 1) fputs("0.32,0.5\n", fp);
	    else if (i == 2) fputs("0.64,0.5\n", fp);
	    else if (i == 3) fputs("0.0,0.0\n", fp);
	    else if (i == 4) fputs("0.32,0.0\n", fp);
	    else if (i == 5) fputs("0.64,0.0\n", fp);
	}
	fprintf(fp, "set xlabel '%s'\n",
		(yvar)? pdinfo->varname[plotlist[i+1]] :
		pdinfo->varname[xvar]);
	fprintf(fp, "set ylabel '%s'\n", 
		(yvar)? pdinfo->varname[yvar] :
		pdinfo->varname[plotlist[i+1]]);
	fputs("plot '-' using 1:2\n", fp);

#ifdef ENABLE_NLS
	setlocale(LC_NUMERIC, "C");
#endif
	for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	    m = (yvar)? plotlist[i+1] : xvar;
	    xx = (*pZ)[m][t];
	    if (na(xx)) fputs("? ", fp);
	    else fprintf(fp, "%.8g ", xx);
	    m = (yvar)? yvar : plotlist[i+1];
	    xx = (*pZ)[m][t];
	    if (na(xx)) fputs("?\n", fp);
	    else fprintf(fp, "%.8g\n", xx);
	}
	fputs("e\n", fp);
#ifdef ENABLE_NLS
	setlocale(LC_NUMERIC, "");
#endif
    } 
    fputs("set nomultiplot\n", fp);
#if defined(OS_WIN32) && !defined(GNUPLOT_PNG)
    fputs("\npause -1\n", fp);
#endif
    fclose(fp);

    if (!(flags & GP_BATCH)) {
	err = gnuplot_display(ppaths);
    }

    free(plotlist);

    return err;
}

/**
 * gnuplot_3d:
 * @list: list of variables to plot, by ID number: Y, X, Z
 * @literal: literal command(s) to pass to gnuplot (or NULL)
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 * @ppaths: path information struct.
 * @plot_count: pointer to counter variable for plots (unused)
 * @flags: unused at present.
 *
 * Writes a gnuplot plot file to display a 3D plot (Z on
 * the vertical axis, X and Y on base plane).
 *
 * Returns: 0 on successful completion, error code on error.
 */

int gnuplot_3d (LIST list, const char *literal,
		double ***pZ, DATAINFO *pdinfo, PATHS *ppaths, 
		int *plot_count, unsigned char flags)
{
    FILE *fq = NULL;
    int t, t1 = pdinfo->t1, t2 = pdinfo->t2, lo = list[0];
    int tmplist[5];
    char surface[64];

    if (lo != 3) {
	fprintf(stderr, "gnuplot_3d needs three variables (only)\n");
	return -1;
    }

    sprintf(ppaths->plotfile, "%sgpttmp.plt", ppaths->userdir);
    fq = fopen(ppaths->plotfile, "w");
    if (fq == NULL) {
	return E_FOPEN;
    }

    _adjust_t1t2(NULL, list, &t1, &t2, *pZ, NULL);
    /* if resulting sample range is empty, complain */
    if (t2 == t1) return -999;

    *surface = 0;

    if (1) {
	MODEL pmod;
	double umin, umax, vmin, vmax;

	tmplist[0] = 4;
	tmplist[1] = list[3];
	tmplist[2] = 0;
	tmplist[3] = list[2];
	tmplist[4] = list[1];

	_minmax(t1, t2, (*pZ)[list[2]], &umin, &umax);
	_minmax(t1, t2, (*pZ)[list[1]], &vmin, &vmax);

	_init_model(&pmod, pdinfo);
	pmod = lsq(tmplist, pZ, pdinfo, OLS, 0, 0.0);
	if (!pmod.errcode && !na(pmod.fstt) &&
	    fdist(pmod.fstt, pmod.dfn, pmod.dfd) < .10) {
	    double uadj = (umax - umin) * 0.02;
	    double vadj = (vmax - vmin) * 0.02;

	    sprintf(surface, "[u=%g:%g] [v=%g:%g] "
		    "%g+(%g)*u+(%g)*v title '', ", 
		    umin - uadj, umax + uadj, 
		    vmin - vadj, vmax + vadj,
		    pmod.coeff[0], pmod.coeff[1],
		    pmod.coeff[2]);
	} 
	clear_model(&pmod, pdinfo);
    }

    fprintf(fq, "set xlabel '%s'\n", get_series_name(pdinfo, list[2]));
    fprintf(fq, "set ylabel '%s'\n", get_series_name(pdinfo, list[1]));
    fprintf(fq, "set zlabel '%s'\n", get_series_name(pdinfo, list[3]));
    fputs("set missing \"?\"\n", fq);

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "C");
#endif

    if (literal != NULL && *literal != 0) {
	print_gnuplot_literal_lines(literal, fq);
    }

#ifdef OS_WIN32
    fprintf(fq, "splot %s'-' title '' w p lt 3\n", surface);
#else
    fprintf(fq, "splot %s'-' title ''\n", surface);
#endif

    /* supply the data to gnuplot inline */
    tmplist[0] = 3;
    tmplist[1] = list[2];
    tmplist[2] = list[1];
    tmplist[3] = list[3];
    for (t=t1; t<=t2; t++) {
	const char *label = NULL;

	if (pdinfo->markers) label = pdinfo->S[t];
	printvars(fq, t, tmplist, *pZ, label, 0.0);
    }	
    fputs("e\n", fq);

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif

    fclose(fq);

    return 0;
}

/* ........................................................... */

static int gnuplot_bars_sane (FREQDIST *freq, int dist, double lambda,
			      double plotmin, double plotmax, 
			      double barwidth, double barskip)
{
    int i;

    for (i=0; i<freq->numbins; i++) { 
	double x1, x2;

	x1 = (i == 0)? freq->endpt[i+1] - barwidth : freq->endpt[i];
	x2 = (i == freq->numbins - 1)? 
	    freq->endpt[i] + barwidth : freq->endpt[i+1];

	if (dist) {
	    if (x1 < plotmin) x1 = plotmin;
	    if (x2 > plotmax) x2 = plotmax;
	}

	if (x1 + barskip >= x2 - barskip) return 0;
    }

    return 1;
}

#if 0
static double tallest_spike (FREQDIST *freq, double d)
{
    int i;
    double x, ret = 0;

    for (i=0; i<freq->numbins; i++) {
	x = d * freq->f[i];
	if (x > ret) ret = x;
    }

    return ret;
}
#endif

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
    char withstring[32];
    double plotmin = 0.0, plotmax = 0.0;
    double barwidth = freq->endpt[K-1] - freq->endpt[K-2];
    double barskip = 0.005 * (freq->endpt[K] - freq->endpt[0]);
    int plottype = PLOT_FREQ_SIMPLE;
    int use_bars = gnuplot_has_filledcurve();

    if (freq->numbins > 16) barskip /= 2.0;

    if (dist == NORMAL) plottype = PLOT_FREQ_NORMAL;
    else if (dist == GAMMA) plottype = PLOT_FREQ_GAMMA;

    if (gnuplot_init(ppaths, plottype, &fp)) return E_FOPEN;

    *withstring = 0;

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "C");
#endif

    fputs("# frequency plot ", fp);

    if (dist) {
	double propn;

	/* find the endpts that straddle the mean... */
	for (i=0; i<K ; i++) 
	    if (freq->endpt[i] > freq->xbar) break;

	/* OK, they are k-1 and k: now find the proportion of the 
	   theoretical distribution they enclose, and calculate a
	   height adjustment factor for the impulses */

	if (dist == NORMAL) {
	    char chilbl[80];

	    fputs("(against normal)\n", fp);

	    propn = normal((freq->endpt[i-1] - freq->xbar)/freq->sdx) -
		normal((freq->endpt[i] - freq->xbar)/freq->sdx);
	    lambda = 1.0 / (propn * freq->n * sqrt(2 * M_PI) * freq->sdx);
	    fprintf(fp, "sigma = %g\n", freq->sdx);
	    fprintf(fp, "mu = %g\n", freq->xbar);

	    plotmin = freq->endpt[1] - barwidth;
	    if (plotmin > freq->xbar - 3.3 * freq->sdx) {
		plotmin = freq->xbar - 3.3 * freq->sdx;
	    }
	    plotmax = freq->endpt[K-1] + barwidth;
	    if (plotmax < freq->xbar + 3.3 * freq->sdx) {
		plotmax = freq->xbar + 3.3 * freq->sdx;
	    }

	    if (!na(freq->chisqu)) {
		fprintf(fp, "set label '%s:' at graph .03, graph .97\n",
			I_("Test statistic for normality"));
		sprintf(chilbl, I_("Chi-squared(2) = %.3f, pvalue %.5f"), 
			freq->chisqu, chisq(freq->chisqu, 2));
		fprintf(fp, "set label '%s' at graph .03, graph .93\n", chilbl);
	    }	
	}
	else if (dist == GAMMA) {
	    double xx, height, var = freq->sdx * freq->sdx;

	    fputs("(against gamma)\n", fp);
	
	    /* scale param = variance/mean */
	    beta = var / freq->xbar;
	    /* shape param = mean/scale */
	    alpha = freq->xbar / beta;

	    propn = gamma_dist(freq->xbar, var, freq->endpt[i], 2) -
		gamma_dist(freq->xbar, var, freq->endpt[i-1], 2);
	    xx = (freq->endpt[i] + freq->endpt[i-1])/2.0;
	    height = pow(xx, alpha - 1.0) * exp(-xx / beta) /
		(cephes_gamma(alpha) * pow(beta, alpha));
	    lambda = height/(freq->n * propn);
	    fprintf(fp, "beta = %g\n", beta);
	    fprintf(fp, "alpha = %g\n", alpha);
	    plotmin = 0.0;
	    plotmax = freq->xbar + 4.0 * freq->sdx;
	}

	/* adjust min, max if needed */
	if (freq->midpt[0] < plotmin) plotmin = freq->midpt[0];
	if (freq->midpt[K-1] > plotmax) plotmax = freq->midpt[K-1];

	fprintf(fp, "set xrange [%.8g:%.8g]\n", plotmin, plotmax);
	fputs("set key right top\n", fp);
    } else { /* plain frequency plot */
	fputs("(simple)\n", fp);

	lambda = 1.0 / freq->n;
	fputs("set nokey\n", fp);
	fprintf(fp, "set xlabel '%s %s'\n", 
		I_("Frequency distribution for"), freq->varname);	
    }

    if (isnan(lambda)) {
	if (fp) fclose(fp);
	return 1;
    }

#if 0
    /* need to allow for the function plotted to go higher */
    fprintf(fp, "set yrange [0:%.8g]\n", 
	    1.2 * tallest_spike(freq, lambda));
#endif

    if (use_bars) {
	/* this won't work for "weird" distributions */
	if (!gnuplot_bars_sane(freq, dist, lambda, plotmin, plotmax,
			       barwidth, barskip)) {
	    use_bars = 0;
	}
    }

    /* plot instructions */
    if (use_bars) {
	strcat(withstring, "w filledcurve");
    } else {
	strcat(withstring, "w impulses");
    }

    if (!dist) {
	fprintf(fp, "plot '-' using 1:($2) %s\n", withstring);
    } else if (dist == NORMAL) {
	fputs("plot \\\n", fp);
	fprintf(fp, "'-' using 1:($2) title '%s' %s , \\\n"
		"(1/(sqrt(2*pi)*sigma)*exp(-(x-mu)**2/(2*sigma**2))) "
		"title 'N(%.4f,%.4f)' w lines\n",
		freq->varname, withstring, freq->xbar, freq->sdx);
    }
    else if (dist == GAMMA) {
	fputs("plot \\\n", fp);
	fprintf(fp, "'-' using 1:($2) title '%s' %s ,\\\n"
		"x**(alpha-1.0)*exp(-x/beta)/(gamma(alpha)*(beta**alpha)) "
		"title 'gamma(%.4f,%.4f)' w lines\n",
		freq->varname, withstring, alpha, beta); 
    }

    /* send sample data inline */

    if (use_bars) {
	/* construct solid histogram bars */
	for (i=0; i<K; i++) { 
	    double y = lambda * freq->f[i];
	    double x1, x2;

	    x1 = (i == 0)? freq->endpt[i+1] - barwidth : freq->endpt[i];
	    x2 = (i == K-1)? freq->endpt[i] + barwidth : freq->endpt[i+1];

	    if (dist) {
		if (x1 < plotmin) x1 = plotmin;
		if (x2 > plotmax) x2 = plotmax;
	    }

	    fprintf(fp, "%.8g 0.0\n", x1 + barskip);
	    fprintf(fp, "%.8g %.8g\n", x1 + barskip, y);
	    fprintf(fp, "%.8g %.8g\n", x2 - barskip, y);
	    fprintf(fp, "%.8g 0\n", x2 - barskip);
	}
    } else {
	/* plot with impulses */
	for (i=0; i<K; i++) { 
	    fprintf(fp, "%.8g %.8g\n", freq->midpt[i], lambda * freq->f[i]);
	}
    }

    fputs("e\n", fp);

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif

#if defined(OS_WIN32) && !defined(GNUPLOT_PNG)
    fputs("pause -1\n", fp);
#endif
    if (fp) fclose(fp);
    return gnuplot_display(ppaths);
}

/* ......................................................... */ 

int plot_fcast_errs (int n, const double *obs, 
		     const double *depvar, const double *yhat, 
		     const double *maxerr, const char *varname, 
		     PATHS *ppaths)
{
    FILE *fp = NULL;
    double xmin, xmax, xrange;
    int t;

    if (gnuplot_init(ppaths, PLOT_FORECAST, &fp)) return E_FOPEN;

    fputs("# forecasts with 95 pc conf. interval\n", fp);

    _minmax(0, n - 1, obs, &xmin, &xmax);
    xrange = xmax - xmin;
    xmin -= xrange * .025;
    xmax += xrange * .025;
    fprintf(fp, "set xrange [%.8g:%.8g]\n", xmin, xmax);

    fprintf(fp, "set key left top\n"
	    "plot \\\n'-' using 1:2 title '%s' w lines , \\\n"
	    "'-' using 1:2 title '%s' w lines , \\\n"
	    "'-' using 1:2:3 title '%s' w errorbars\n", 
	    varname, I_("fitted"), I_("95 percent confidence interval"));

    /* send data inline */
#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "C");
#endif
    for (t=0; t<n; t++) {
	fprintf(fp, "%.8g %.8g\n", obs[t], depvar[t]);
    }
    fputs("e\n", fp);
    for (t=0; t<n; t++) {
	fprintf(fp, "%.8g %.8g\n", obs[t], yhat[t]);
    }
    fputs("e\n", fp);
    for (t=0; t<n; t++) {
	fprintf(fp, "%.8g %.8g %.8g\n", obs[t], yhat[t], maxerr[t]);
    }
    fputs("e\n", fp);
#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif

#if defined(OS_WIN32) && !defined(GNUPLOT_PNG)
    fputs("pause -1\n", fp);
#endif
    fclose(fp);

    return gnuplot_display(ppaths);
}

/* ........................................................... */

void free_plotspec (GPT_SPEC *spec)
{
    int i;

    if (spec == NULL) return;

    if (spec->lines != NULL) free(spec->lines);
    if (spec->data != NULL) free(spec->data);

    for (i=0; i<4; i++) {
	if (spec->literal[i] != NULL) {
	    free(spec->literal[i]);
	}
    }

    if (spec->labels != NULL) {
	for (i=0; i<spec->nlabels; i++) {
	    free(spec->labels[i]);
	}
	free(spec->labels);
    }

    free(spec);
}

/* ........................................................... */

int termtype_to_termstr (const char *termtype, char *termstr,
			 const PATHS *ppaths)
{
    int cmds = 0;

    if (!strcmp(termtype, "postscript color")) 
	strcpy(termstr, "postscript eps color"); 
    else if (!strcmp(termtype, "postscript")) 
	strcpy(termstr, "postscript eps"); 
    else if (!strcmp(termtype, "fig")) 
	strcpy(termstr, "fig");
    else if (!strcmp(termtype, "latex")) 
	strcpy(termstr, "latex");
    else if (!strcmp(termtype, "png")) { 
	const char *png_str = 
	    get_gretl_png_term_line(ppaths, PLOT_REGULAR);

	strcpy(termstr, png_str + 9);
    }
    else if (!strcmp(termtype, "plot commands")) 
	cmds = 1;
    else strcpy(termstr, termtype);
    return cmds;
}

/* ........................................................... */

static char *escape_quotes (const char *s)
{
    if (strchr(s, '"') == NULL) return NULL;
    else {
	int qcount = 0;
	char *ret, *r;
	const char *p = s;

	while (*p) {
	    if (*p == '"') qcount++;
	    p++;
	}
	ret = malloc(strlen(s) + 1 + qcount);
	if (ret == NULL) return NULL;

	r = ret;
	while (*s) {
	    if (*s == '"') {
		*r++ = '\\';
		*r++ = '"';
	    } else {
		*r++ = *s;
	    }
	    s++;
	}

	*r = 0;
	return ret;
    }
}

/* ........................................................... */

int print_plotspec_details (const GPT_SPEC *spec, FILE *fp)
{
    int i, k, t, datlines;
    int plotn, nlines = spec->nlines;
    int miss = 0;
    double xx;

    if (!string_is_blank(spec->titles[0])) {
	if ((spec->flags & GPTSPEC_OLS_HIDDEN) && 
	    is_auto_ols_string(spec->titles[0])) ;
	else fprintf(fp, "set title '%s'\n", spec->titles[0]);
    }
    if (!string_is_blank(spec->titles[1])) {
	fprintf(fp, "set xlabel '%s'\n", spec->titles[1]);
    }
    if (!string_is_blank(spec->titles[2])) {
	fprintf(fp, "set ylabel '%s'\n", spec->titles[2]);
    }
    if ((spec->flags & GPTSPEC_Y2AXIS) && !string_is_blank(spec->titles[3])) {
	fprintf(fp, "set y2label '%s'\n", spec->titles[3]);
    }

    for (i=0; i<MAX_PLOT_LABELS; i++) {
	if (!string_is_blank(spec->text_labels[i].text)) {
	    char *label = escape_quotes(spec->text_labels[i].text);

	    fprintf(fp, "set label \"%s\" at %s %s\n", 
		    (label != NULL)? label : spec->text_labels[i].text,
		    spec->text_labels[i].pos,
		    spec->text_labels[i].just);
	    if (label != NULL) free(label);
	}
    }

    fputs("set xzeroaxis\n", fp);
    fputs("set missing \"?\"\n", fp);

    if (strcmp(spec->keyspec, "none") == 0) {
	fputs("set nokey\n", fp);
    } else {
	fprintf(fp, "set key %s\n", spec->keyspec);
    }

    k = (spec->flags & GPTSPEC_Y2AXIS)? 3: 2;
    for (i=0; i<k; i++) {
	fprintf(fp, "set %srange [%s:%s]\n",
		(i==0)? "x" : (i==1)? "y" : "y2",
		spec->range[i][0], spec->range[i][1]);
    }

    /* customized xtics? */
    if (!string_is_blank(spec->xtics)) {
	fprintf(fp, "set xtics %s\n", spec->xtics);
    }
    if (!string_is_blank(spec->mxtics)) {
	fprintf(fp, "set mxtics %s\n", spec->mxtics);
    }

    /* using two y axes? */
    if (spec->flags & GPTSPEC_Y2AXIS) {
	fputs("set ytics nomirror\n", fp);
	fputs("set y2tics\n", fp);
    }

    /* in case of plots that are editable (see gui client), it is
       important to re-write the comment string that identifies the
       sort of graph, so that it will be recognized by type when
       it is redisplayed */

    if (spec->code == PLOT_FORECAST) {
	fputs("# forecasts with 95 pc conf. interval\n", fp);
    }
    else if (spec->code == PLOT_CORRELOGRAM) {
	fputs("# correlogram\n", fp); 
    }
    else if (spec->code == PLOT_FREQ_SIMPLE) {
	fputs("# frequency plot (simple)\n", fp); 
    }
    else if (spec->code == PLOT_FREQ_NORMAL || 
	     spec->code == PLOT_FREQ_GAMMA ||
	     spec->code == PLOT_PERIODOGRAM) { 
	if (spec->code == PLOT_FREQ_NORMAL) {
	    fputs("# frequency plot (against normal)\n", fp); 
	} else if (spec->code == PLOT_FREQ_GAMMA) {
	    fputs("# frequency plot (against gamma)\n", fp); 
	} else {
	    fputs("# periodogram\n", fp);
	}
	for (i=0; i<4; i++) {
	    fprintf(fp, "%s\n", spec->literal[i]);
	}
    } 

    if (spec->flags & GPTSPEC_AUTO_OLS) {
	fputs(auto_ols_string, fp);
	if ((spec->flags & GPTSPEC_OLS_HIDDEN) && nlines > 1) {
	    nlines--;
	}
    }

    fputs("plot \\\n", fp);

    datlines = nlines; /* ? */ 

    for (i=0; i<nlines; i++) {
	if (strcmp(spec->lines[i].scale, "NA")) {
	    fprintf(fp, "'-' using 1:($2*%s) ", 
		    spec->lines[i].scale);
	} else {
	    fprintf(fp, "%s ", spec->lines[i].formula); 
	    datlines--;
	}
	if (spec->lines[i].yaxis != 1) {
	    fprintf(fp, "axes x1y%d ", spec->lines[i].yaxis);
	}
	fprintf(fp, "title '%s' w %s", 
		spec->lines[i].title,
		spec->lines[i].style);
	if (i == nlines - 1) {
	    fputc('\n', fp);
	} else {
	    fputs(", \\\n", fp);
	}
    } 

    /* supply the data to gnuplot inline */
#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "C");
#endif
    miss = 0;
    plotn = spec->t2 - spec->t1 + 1;
    for (i=1; i<=datlines; i++) {  
	for (t=spec->t1; t<=spec->t2; t++) {
	    xx = spec->data[t - spec->t1];
	    if (na(xx)) {
		fputs("? ", fp);
		miss = 1;
	    } else {
		fprintf(fp, "%.8g ", xx);
	    }
	    xx = spec->data[plotn * i + t - spec->t1];
	    if (na(xx)) {
		fputc('?', fp);
		miss = 1;
	    } else {
		fprintf(fp, "%.8g", xx);
	    }
	    if (spec->labels != NULL && i == 1) {
		fprintf(fp, " # %s", spec->labels[t]);
	    }
	    fputc('\n', fp);
	}
	fputs("e\n", fp);
    }
#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif

    return miss;
}

/* ........................................................... */

int go_gnuplot (GPT_SPEC *spec, char *fname, PATHS *ppaths)
     /* ship out a plot struct, to gnuplot or file */
{
    FILE *fp = NULL;
    int dump = 0;
    int err = 0, miss;
    char termstr[72];

    dump = termtype_to_termstr(spec->termtype, termstr, ppaths);

    if (dump) {  
	/* dump of gnuplot commands to named file */
	if (fname == NULL) return 1;  /* impossible */
	fp = fopen(fname, "w");
	if (fp == NULL) return 1;
    } else {     
	/* output to gnuplot, for screen or other "term" */
#ifdef GNUPLOT_PIPE
	fp = spec->fp; /* pipe */
#else
	if (spec->fp == NULL) {
	    fp = fopen(ppaths->plotfile, "w");
	}
	if (fp == NULL) return 1;
#endif /* GNUPLOT_PIPE */
	if (fname != NULL) { 
	    /* file, not screen display */
	    fprintf(fp, "set term %s\n", termstr);
#ifdef ENABLE_NLS
	    if (strstr(termstr, "postscript")) {
		fputs("set encoding iso_8859_1\n", fp);
	    }
#endif /* ENABLE_NLS */
	    fprintf(fp, "set output '%s'\n", fname);
	}
    }

    miss = print_plotspec_details(spec, fp);
    fflush(fp);

    if (dump) {
	/* we're finished */
	fclose(fp);
    }
    
#ifndef GNUPLOT_PIPE
    if (!dump) {
	char plotcmd[MAXLEN];
# ifdef OS_WIN32
	int winshow = 0;

	if (fname == NULL) { /* sending plot to screen */
	    fputs("pause -1\n", fp);
	    winshow = 1;
	} 
# endif
	fclose(fp);
	spec->fp = NULL;
	sprintf(plotcmd, "\"%s\" \"%s\"", ppaths->gnuplot, ppaths->plotfile);
# ifdef OS_WIN32
	if (winshow) {
	    err = (WinExec(plotcmd, SW_SHOWNORMAL) < 32);
	} else {
	    err = winfork(plotcmd, NULL, SW_SHOWMINIMIZED, 0);
	}
# else
	if (system(plotcmd)) err = 1;
# endif 
    }
#endif /* GNUPLOT_PIPE */
    if (miss) err = 2;
    return err;
}

/* ........................................................... */

int rmplot (const LIST list, double **Z, DATAINFO *pdinfo, PRN *prn,
	    PATHS *ppaths)
{
    int err;
    void *handle;
    int (*range_mean_graph) (int, double **, const DATAINFO *, 
                             PRN *, PATHS *);

    if (open_plugin("range-mean", &handle)) return 1;

    range_mean_graph = get_plugin_function("range_mean_graph", handle);

    if (range_mean_graph == NULL) {
        pputs(prn, _("Couldn't load plugin function\n"));
        close_plugin(handle);
        return 1;
    }

    err = range_mean_graph (list[1], Z, pdinfo, prn, ppaths);

    close_plugin(handle);

    if (!err) {
        return gnuplot_display(ppaths);
    } else {
	return err;
    }
}

/* ........................................................... */

int is_auto_ols_string (const char *s)
{
    if (strstr(s, "automatic OLS")) return 1;
    if (strstr(s, I_("with least squares fit"))) return 1;
    return 0;
}

/* ........................................................... */

static char gnuplot_pallette[3][8] = {
    "xff0000", 
    "x0000ff", 
    "x00cc00"  /* full-intensity green is not very legible */
};

static const char *get_gnuplot_pallette (int i, int ptype)
{
#ifdef GNUPLOT_DEBUG
    fprintf(stderr, "get_gnuplot_pallette: i=%d, ptype=%d\n",
	    i, ptype);
#endif
    if (i == 0 && (ptype == PLOT_FREQ_SIMPLE ||
		   ptype == PLOT_FREQ_NORMAL || 
		   ptype == PLOT_FREQ_GAMMA)) {
	return "xaabbcc";
    }
    else if (i >= 0 && i < 3) {
	return gnuplot_pallette[i];
    } else {
	return "";
    }
}

static int colstr_is_valid (const char *colstr)
{
    int i;
    const char *ok = "0123456789abcdef";

    if (*colstr != 'x') return 0;
    if (strlen(colstr) != 7) return 0;
    for (i=1; i<7; i++) {
	if (strchr(ok, colstr[i]) == NULL) return 0;
    }

    return 1;
}

void set_gnuplot_pallette (int i, const char *colstr)
{
    if (i >= 0 && i < 3 && colstr_is_valid(colstr)) {
	strcpy(gnuplot_pallette[i], colstr);
    } else {
	fprintf(stderr, "Invalid color spec, '%s'\n", colstr);
    }
}
