#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

static void graphpage_init (graphpage *gpage)
{
    gpage->output = PS_OUTPUT;
    gpage->color = 0;
    gpage->ngraphs = 0;
    gpage->fnames = NULL;
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

static int make_graphpage_tex (graphpage *gpage)
{
    FILE *fp;
    int err = 0;

    fp = fopen(gpage_tex, "w");
    if (fp == NULL) return 1;

    doctop(gpage->output, fp);

    err = geomline(gpage->ngraphs, fp);

    if (!err) {
	common_setup(fp);
	err = tex_graph_setup(gpage->ngraphs, fp);
    }

    if (!err) common_end(fp);

    fclose(fp);

    return err;
}

static char *my_strdup (const char *s)
{
    char *ret = malloc(strlen(s) + 1);

    if (ret != NULL) strcpy(ret, s);
    return ret;
}

static int graphpage_add_file (graphpage *gpage, const char *fname)
{
    char **fnames;
    int ng = gpage->ngraphs + 1;

    fnames = realloc(gpage->fnames, ng * sizeof *fnames);
    if (fnames == NULL) return 1;

    fnames[ng - 1] = my_strdup(fname);
    gpage->fnames = fnames;
    gpage->ngraphs = ng;

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
	fprintf(stderr, "Gnuplot: doing PNG\n");
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

static int latex_compile (graphpage *gpage)
{
    char cmd[128];
    int err;

    if (gpage->output == PDF_OUTPUT) {
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

static int make_gp_output (graphpage *gpage)
{
    double scale = 1.0;
    int i;
    int err = 0;

    if (gpage->ngraphs == 3) {
	scale = 0.8;
    } else if (gpage->ngraphs > 3) {
	scale = 0.75;
    }

    for (i=0; i<gpage->ngraphs && !err; i++) {
	err = gp_make_outfile(gpage->fnames[i], i + 1, scale, 
			      gpage->output, gpage->color);
    }

    remove(gpage_plt); 

    return err;
}

static int real_display_gpage (graphpage *gpage)
{
    char cmd[128];

    if (gpage->output == PDF_OUTPUT) {
	sprintf(cmd, "acroread %s.pdf", gpage_base);
    } else {
	sprintf(cmd, "gv %s.ps", gpage_base);
    }

    return system(cmd);
}

static void gpage_cleanup (graphpage *gpage)
{
    char fname[128];
    int i;

    for (i=0; i<gpage->ngraphs; i++) {
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

static int display_graphpage (graphpage *gpage)
{
    int err = 0;

    /* write the LaTeX driver file */
    err = make_graphpage_tex(gpage);

    if (!err) {
	/* transform individual plot files and compile 
	   using gnuplot */
	err = make_gp_output(gpage);
    }

    if (!err) {
	err = latex_compile(gpage);
    }

    if (!err) {
	/* compile LaTeX and display output */
	real_display_gpage(gpage);
    }

    gpage_cleanup(gpage);

    return err;
}

int main (int argc, char **argv)
{
    graphpage gpage;
    const char *gnames[] = {
	"test.gp", "test.gp", "test.gp", "test.gp",
	"test.gp", "test.gp", NULL
    };
    int i, ng = 4;

    if (argc > 1) {
	ng = atoi(argv[1]);
	if (ng == 0) exit(EXIT_FAILURE);
    }

    graphpage_init(&gpage);

    if (argc == 3) {
	if (!strcmp(argv[2], "pdf")) {
	    gpage.output = PDF_OUTPUT;
	}
	else if (!strcmp(argv[2], "color")) {
	    gpage.color = 1;
	}
    }

    for (i=0; i<ng; i++) {
	graphpage_add_file(&gpage, gnames[i]);
    }

    display_graphpage(&gpage);

    return 0;
}
