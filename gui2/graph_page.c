/*
 *  Copyright (c) 2004 by Allin Cottrell
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

#include "gretl.h"
#include "graph_page.h"

enum {
    PS_OUTPUT,
    PDF_OUTPUT
} output_types;

typedef struct _graphpage graphpage;

struct _graphpage {
    int output;
    int color;
    int ngraphs;
    char **fnames;
};

const char *gpage_base = "gretl_graphpage";
const char *gpage_tex = "gretl_graphpage.tex";
const char *gpage_plt = "gretl_graphpage.plt";

static graphpage gpage;

static void gpage_errmsg (char *msg, int gui)
{
    if (gui) {
	errbox(msg);
    } else {
	gretl_errmsg_set(msg);
    }
}

static void graph_page_init (void)
{
    gpage.output = PS_OUTPUT;
    gpage.color = 0;
    gpage.ngraphs = 0;
    gpage.fnames = NULL;
}

static void doctop (int otype, FILE *fp)
{
    fputs("\\documentclass{article}\n", fp);
    fprintf(fp, "\\usepackage[%s]{graphicx}\n", 
	    (otype == PS_OUTPUT)? "dvips" : "pdftex");
}

static int geomline (int ng, FILE *fp)
{
    double width = 6.5, height = 9.0;
    double tmarg = 1.0, lmarg = 1.0;

    if (ng == 3) {
	height = 9.5;
	tmarg = 0.8;
    } else if (ng < 7) {
	width = 7.0;
	lmarg = 0.75;
    } else {
	fprintf(stderr, "ng (%d) out of range\n", ng);
	return 1;
    }

    fprintf(fp, "\\usepackage[body={%gin,%gin},"
	    "top=%gin,left=%gin,nohead]{geometry}\n\n",
	    width, height, tmarg, lmarg);

    return 0;
}

static void common_setup (FILE *fp)
{
    fputs("\\begin{document}\n\n"
	  "\\thispagestyle{empty}\n\n"
	  "\\begin{center}\n", fp);
}

static int tex_graph_setup (int ng, FILE *fp)
{
    char fname[FILENAME_MAX];
    double scale = 1.0;
    double vspace = 1.0;
    int i;

    if (ng > 6) {
	fprintf(stderr, "ng (%d) out of range\n", ng);
	return 1;
    }

    if (ng == 1) {
	fputs("\n\\vspace{2in}\n\n", fp);
	strcpy(fname, "gretl_gpage_1");
	fprintf(fp, "\\includegraphics{%s}\n", fname);
    }

    else if (ng == 2) {
	sprintf(fname, "gretl_gpage_%d", 1);
	fprintf(fp, "\\includegraphics{%s}\n\n", fname);
	fprintf(fp, "\\vspace{%gin}\n\n", vspace);
	sprintf(fname, "gretl_gpage_%d", 2);
	fprintf(fp, "\\includegraphics{%s}\n\n", fname);
    }

    else if (ng == 3) {
	scale = 0.8;
	vspace = 0.25;
	for (i=0; i<3; i++) {
	    sprintf(fname, "gretl_gpage_%d", i + 1);
	    fprintf(fp, "\\includegraphics[scale=%g]{%s}\n\n",
		    scale, fname);
	    fprintf(fp, "\\vspace{%gin}\n", vspace);
	}
    }    

    else {
	scale = 0.9;
	vspace = 0.25;
	fputs("\\begin{tabular}{cc}\n", fp);
	for (i=0; i<ng; i++) {
	    sprintf(fname, "gretl_gpage_%d", i + 1);
	    fprintf(fp, "\\includegraphics[scale=%g]{%s}",
		    scale, fname);
	    if (i % 2 == 0) {
		fputs(" &\n  ", fp);
	    } else {
		fprintf(fp, " \\\\ [%gin]\n", vspace);
	    }
	}
	fputs("\\end{tabular}\n", fp);
    } 

    return 0;
}

static void common_end (FILE *fp)
{
    fputs("\\end{center}\n\n"
	  "\\end{document}\n", fp);
}

static int make_graphpage_tex (void)
{
    FILE *fp;
    int err = 0;

    fp = fopen(gpage_tex, "w");
    if (fp == NULL) return 1;

    doctop(gpage.output, fp);

    err = geomline(gpage.ngraphs, fp);

    if (!err) {
	common_setup(fp);
	err = tex_graph_setup(gpage.ngraphs, fp);
    }

    if (!err) common_end(fp);

    fclose(fp);

    return err;
}

int graph_page_add_file (const char *fname)
{
    char **fnames;
    int ng = gpage.ngraphs + 1;

    fnames = myrealloc(gpage.fnames, ng * sizeof *fnames);
    if (fnames == NULL) return 1;

    fnames[ng - 1] = g_strdup(fname);
    gpage.fnames = fnames;
    gpage.ngraphs = ng;

    return 0;
}

static int gnuplot_compile (const char *fname)
{
    char cmd[128];

    sprintf(cmd, "gnuplot %s", fname);
    return system(cmd);
}

static int gp_make_outfile (const char *fname, int i, double scale,
			    int output, int color)
{
    char outfile[128], line[128];
    FILE *fp, *fq;
    int err = 0;

    fp = fopen(fname, "r");
    if (fp == NULL) return 1;

    fq = fopen(gpage_plt, "w");
    if (fq == NULL) {
	fclose(fp);
	return 1;
    }
    
    if (output == PDF_OUTPUT) {
	fprintf(fq, "set term png%s font verdana 6 size %g,%g\n", 
		((color)? "" : " mono"),
		480.0 * scale, 360.0 * scale);
	sprintf(outfile, "gretl_gpage_%d.png", i);
    } else {
	fprintf(fq, "set term postscript eps%s\n", (color)? " color" : "");
	sprintf(outfile, "gretl_gpage_%d.ps", i);
	if (scale != 1.0) {
	    fprintf(fq, "set size %g,%g\n", scale, scale);
	}
    }

    fprintf(fq, "set output '%s'\n", outfile);

    while (fgets(line, sizeof line, fp)) {
	if (!strncmp(line, "set out", 7)) continue;
	if (!strncmp(line, "set term", 8)) continue;
	if (!strncmp(line, "set size", 8)) continue;
	fputs(line, fq);
    }

    fclose(fp);
    fclose(fq);

    err = gnuplot_compile(gpage_plt);

    return err;
}

static int latex_compile_graph_page (void)
{
    char cmd[128];
    int err;

    if (gpage.output == PDF_OUTPUT) {
	sprintf(cmd, "pdflatex %s", gpage_base);
	err = system(cmd);
    } else {
	sprintf(cmd, "latex %s", gpage_base);
	err = system(cmd);
	if (!err) {
	    sprintf(cmd, "dvips %s", gpage_base);
	    err = system(cmd);
	}
    }

    return err;
}

static int make_gp_output (void)
{
    double scale = 1.0;
    int i;
    int err = 0;

    if (gpage.ngraphs == 3) {
	scale = 0.8;
    } else if (gpage.ngraphs > 3) {
	scale = 0.75;
    }

    for (i=0; i<gpage.ngraphs && !err; i++) {
	err = gp_make_outfile(gpage.fnames[i], i + 1, scale, 
			      gpage.output, gpage.color);
    }

    remove(gpage_plt); 

    return err;
}

static int real_display_gpage (void)
{
    char cmd[128];

    if (gpage.output == PDF_OUTPUT) {
	sprintf(cmd, "acroread %s.pdf", gpage_base);
    } else {
	sprintf(cmd, "gv %s.ps", gpage_base);
    }

    return system(cmd);
}

static void gpage_cleanup (void)
{
    char fname[128];
    int i;

    for (i=0; i<gpage.ngraphs; i++) {
	sprintf(fname, "gretl_gpage_%d.ps", i + 1);
	remove(fname);
    }

    sprintf(fname, "%s.tex", gpage_base);
    remove(fname);
    sprintf(fname, "%s.dvi", gpage_base);
    remove(fname);
    sprintf(fname, "%s.log", gpage_base);
    remove(fname);
    sprintf(fname, "%s.aux", gpage_base);
    remove(fname);
}

int display_graph_page (void)
{
    int err = 0;

    if (gpage.ngraphs == 0) {
	gpage_errmsg(_("The graph page is empty"), 1);
	return 1;
    }

    /* write the LaTeX driver file */
    err = make_graphpage_tex();

    if (!err) {
	/* transform individual plot files and compile 
	   using gnuplot */
	err = make_gp_output();
    }

    if (!err) {
	err = latex_compile_graph_page();
    }

    if (!err) {
	/* compile LaTeX and display output */
	real_display_gpage();
    }

    gpage_cleanup();

    return err;
}

void clear_graph_page (void)
{
    int i;

    for (i=0; i<gpage.ngraphs; i++) {
	free(gpage.fnames[i]);
    }

    free(gpage.fnames);

    graph_page_init();
}


