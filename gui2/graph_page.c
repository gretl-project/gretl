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

#include "gretl.h"
#include "graph_page.h"
#include "session.h"
#include "gpt_control.h"
#include "texprint.h"
#include "guiprint.h"

#ifdef G_OS_WIN32 
# include <io.h>
# include "gretlwin32.h"
#else
# include <unistd.h>
# include <sys/stat.h>
#endif

#define GRAPHS_MAX 8

typedef struct _graphpage graphpage;

struct _graphpage {
    int term;
    int mono;
    int ngraphs;
    char **fnames;
};

static graphpage gpage;
static char gpage_base[FILENAME_MAX];
static char gpage_tex_base[FILENAME_MAX];

static void gpage_filenames_init (const char *base)
{
    if (base == NULL) {
	strcpy(gpage_base, gretl_dotdir());
	strcat(gpage_base, "gretl_graphpage");
	strcpy(gpage_tex_base, "gretl_graphpage");
    } else {
	const char *p;

	strcpy(gpage_base, base);
	if (has_suffix(gpage_base, ".tex") || 
	    has_suffix(gpage_base, ".pdf") ||
	    has_suffix(gpage_base, ".eps")) {
	    gpage_base[strlen(gpage_base) - 4] = '\0';
	} else if (has_suffix(gpage_base, ".ps")) {
	    gpage_base[strlen(gpage_base) - 3] = '\0';
	}
	p = strrchr(gpage_base, SLASH);
	if (p != NULL) {
	    strcpy(gpage_tex_base, p + 1);
	} else {
	    strcpy(gpage_tex_base, gpage_base);
	}
    }
}

static char *gpage_fname (const char *ext, int i)
{
    static char fname[MAXLEN];

    strcpy(fname, gpage_base);

    if (i > 0) {
	char num[6];

	sprintf(num, "_%d", i);
	strcat(fname, num);
    }

    if (ext != NULL) {
	strcat(fname, ext);
    }

    return fname;
}

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
    gpage.mono = 0;
    gpage.ngraphs = 0;
    gpage.fnames = NULL;
}

static void doctop (FILE *fp)
{
    const char *paper = in_usa()? "letterpaper" : "a4paper";
    const char *driver = "pdftex";

    if (gpage.term == GP_TERM_EPS) {
	driver = "dvips";
    }

    fprintf(fp, "\\documentclass[%s]{article}\n", paper);
    fprintf(fp, "\\usepackage[%s]{graphicx}\n", driver);
}

/* A4 is 210mm * 297mm */

static int geomline (int ng, FILE *fp)
{
    double width = 7.0, height = 10.0;
    double tmarg = 0.5, lmarg = 0.75;
    char unit[3];

    if (in_usa()) {
	strcpy(unit, "in");
    } else {
	strcpy(unit, "mm");
	width = 190;
	height = 277.0;
	tmarg = lmarg = 10.0;
    }

    gretl_push_c_numeric_locale();

    fprintf(fp, "\\usepackage[body={%g%s,%g%s},"
	    "top=%g%s,left=%g%s,nohead]{geometry}\n\n",
	    width, unit, height, unit, 
	    tmarg, unit, lmarg, unit);

    gretl_pop_c_numeric_locale();

    return 0;
}

static void common_setup (FILE *fp)
{
    fputs("\\begin{document}\n\n"
	  "\\thispagestyle{empty}\n\n"
	  "\\vspace*{\\stretch{1}}\n\n"
	  "\\begin{center}\n", fp);
}

static int oddgraph (int ng, int i)
{
    return (ng % 2) && (i == ng - 1);
}

enum {
    SCALE_LARGE,
    SCALE_REGULAR,
    SCALE_MEDIUM,
    SCALE_SMALL
};

static int tex_graph_setup (int ng, FILE *fp)
{
    char fname[FILENAME_MAX];
    double scale[] = { 1.2, 1.0, 0.9, 0.85 };
    double s, vspace = 1.0;
    int i;

    if (ng > GRAPHS_MAX) {
	fprintf(stderr, "ng (%d) out of range\n", ng);
	return 1;
    }

    gretl_push_c_numeric_locale();

    if (ng == 1) {
	s = scale[SCALE_LARGE];
	sprintf(fname, "%s_1", gpage_tex_base);
	fprintf(fp, "\\includegraphics[scale=%g]{%s}\n", s, fname);
    } else if (ng == 2) {
	s = scale[SCALE_REGULAR];
	sprintf(fname, "%s_%d", gpage_tex_base, 1);
	fprintf(fp, "\\includegraphics[scale=%g]{%s}\n\n", s, fname);
	fprintf(fp, "\\vspace{%gin}\n\n", vspace);
	sprintf(fname, "%s_%d", gpage_tex_base, 2);
	fprintf(fp, "\\includegraphics{%s}\n\n", fname);
    } else if (ng == 3) {
	s = scale[SCALE_MEDIUM];
	vspace = 0.25;
	for (i=0; i<3; i++) {
	    sprintf(fname, "%s_%d", gpage_tex_base, i + 1);
	    fprintf(fp, "\\includegraphics[scale=%g]{%s}\n\n", s, fname);
	    fprintf(fp, "\\vspace{%gin}\n", vspace);
	}
    } else {
	if (ng > 6) {
	    s = scale[SCALE_SMALL];
	    vspace = 0.20;
	} else {
	    s = scale[SCALE_MEDIUM];
	    vspace = 0.25;
	}	    
	fputs("\\begin{tabular}{cc}\n", fp);
	for (i=0; i<ng; i++) {
	    sprintf(fname, "%s_%d", gpage_tex_base, i + 1);
	    if (oddgraph(ng, i)) {
		fprintf(fp, "\\multicolumn{2}{c}{\\includegraphics[scale=%g]{%s}}",
			s, fname);
	    } else {
		fprintf(fp, "\\includegraphics[scale=%g]{%s}", s, fname);
		if (i % 2 == 0) {
		    fputs(" &\n  ", fp);
		} else if (i < ng - 1) {
		    fprintf(fp, " \\\\ [%gin]\n", vspace);
		}
	    }
	}
	fputs("\\end{tabular}\n", fp);
    }

    gretl_pop_c_numeric_locale();

    return 0;
}

static void common_end (FILE *fp)
{
    fputs("\\end{center}\n\n"
	  "\\vspace*{\\stretch{2}}\n\n"
	  "\\end{document}\n", fp);
}

static int make_graphpage_tex (void)
{
    char *fname;
    FILE *fp;
    int err = 0;

    fname = gpage_fname(".tex", 0);

    fp = gretl_fopen(fname, "w");
    if (fp == NULL) {
	return 1;
    }

    doctop(fp);

    err = geomline(gpage.ngraphs, fp);

    if (!err) {
	common_setup(fp);
	err = tex_graph_setup(gpage.ngraphs, fp);
    }

    if (!err) {
	common_end(fp);
    }

    fclose(fp);

    return err;
}

int graph_page_add_file (const char *fname)
{
    char **fnames;
    int ng = gpage.ngraphs + 1;

    if (ng > GRAPHS_MAX) {
	gpage_errmsg(_("The graph page is full"), 1);
	return 1;
    }

    fnames = myrealloc(gpage.fnames, ng * sizeof *fnames);
    if (fnames == NULL) {
	return E_ALLOC;
    }

    fnames[ng - 1] = g_strdup(fname);
    gpage.fnames = fnames;
    gpage.ngraphs = ng;

    mark_session_changed();

    return 0;
}

int in_graph_page (const char *fname)
{
    int i;

    for (i=0; i<gpage.ngraphs; i++) {
	if (!strcmp(fname, gpage.fnames[i])) {
	    return 1;
	}
    }

    return 0;
}

static int graph_page_add_last_graph (void)
{
    const char *fname = last_session_graph_name();
    int err = 0;

    if (fname != NULL && *fname != '\0') {
	if (!in_graph_page(fname)) {
	    err = graph_page_add_file(fname);
	}
    }

    return err;
}

static int gnuplot_compile (const char *fname)
{
    const char *gnuplot = gretl_gnuplot_path();
    char plotcmd[MAXLEN];
    int err = 0;

#ifdef G_OS_WIN32
    sprintf(plotcmd, "\"%s\" \"%s\"", gnuplot, fname);
#else
    sprintf(plotcmd, "%s \"%s\"", gnuplot, fname);
#endif

    err = gretl_spawn(plotcmd);

    return err;
}

static double gp_fontscale = 1.0;

static int graph_page_set_font_scale (const char *s)
{
    double x;

    s += strspn(s, " ");
    x = dot_atof(s);

    if (x > 0.0 && x < 4.0) {
	gp_fontscale = x;
	return 0;
    } else {
	gretl_errmsg_sprintf("'%s': invalid fontscale", s);
	return E_DATA;
    }
}

static int gp_pdfcairo_fontsize (void)
{
#ifdef G_OS_WIN32
    int basesize = 10;
#else
    static int basesize;
#endif
    int fontsize;

#ifndef G_OS_WIN32
    if (basesize == 0) {
	double gpver = gnuplot_version();

	basesize = (gpver > 4.4)? 10 : 5;
    }
#endif

    if (gp_fontscale != 1.0) {
	fontsize = basesize * gp_fontscale;
	if (fontsize < 4) {
	    fontsize = 4;
	}
    } else {
	fontsize = basesize;
    }

    return fontsize;
}

static int gp_make_outfile (const char *gfname, int i, double scale)
{
    int latin = 0;
    int pdfterm = 0;
    char *fname;
    FILE *fp, *fq;
    int err = 0;

    fp = gretl_fopen(gfname, "r");
    if (fp == NULL) {
	file_read_errbox(gfname);
	return E_FOPEN;
    }

    fname = gpage_fname(".plt", 0);

    fq = gretl_fopen(fname, "w");
    if (fq == NULL) {
	file_write_errbox(fname);
	fclose(fp);
	return E_FOPEN;
    }

    if (gpage.term == GP_TERM_PDF) {
	pdfterm = gnuplot_pdf_terminal();
    }

    if (gpage.term == GP_TERM_EPS || pdfterm != GP_PDF_CAIRO) {
	latin = iso_latin_version();
	fprintf(fq, "set encoding iso_8859_%d\n", latin);
    }

    gretl_push_c_numeric_locale();
    
    if (gpage.term == GP_TERM_PDF) {
	if (pdfterm == GP_PDF_CAIRO) {
	    int fontsize = gp_pdfcairo_fontsize();

	    fprintf(fq, "set term pdfcairo font \"sans,%d\"%s", fontsize,
		    (gpage.mono)? " monochrome dashed" : " ");
	} else {
	    fprintf(fq, "set term pdf%s", (gpage.mono)? " monochrome dashed" : " color");
	}
	fname = gpage_fname(".pdf", i);
    } else {
	fprintf(fq, "set term postscript eps%s", (gpage.mono)? " monochrome" : " color");
	fname = gpage_fname(".ps", i);
    }

    if (scale != 1.0) {
	fprintf(fq, " size %g,%g\n", scale * 5.0, scale * 3.5);
    } else {
	fputc('\n', fq);
    }    

    gretl_pop_c_numeric_locale();

    fprintf(fq, "set output '%s'\n", fname);

    filter_gnuplot_file(gpage.term, latin, gpage.mono, fp, fq);

    fclose(fp);
    fclose(fq);

    fname = gpage_fname(".plt", 0);
    err = gnuplot_compile(fname);

    return err;
}

#if defined(G_OS_WIN32)

static int get_dvips_path (char *path)
{
    int ret;
    char *p;

    ret = SearchPath(NULL, "dvips.exe", NULL, MAXLEN, path, &p);

    return (ret == 0);
}

#else

static int spawn_dvips (char *texsrc)
{
    GError *error = NULL;
    gchar *sout = NULL;
    gchar *argv[5];
    char outfile[MAXLEN]; 
    int ok, status;
    int ret = 0;

    sprintf(outfile, "%s.ps", texsrc);

    argv[0] = "dvips";
    argv[1] = "-o";
    argv[2] = outfile;
    argv[3] = texsrc;
    argv[4] = NULL;

    ok = g_spawn_sync(gretl_dotdir(),
		      argv,
		      NULL,    /* envp */
		      G_SPAWN_SEARCH_PATH |
		      G_SPAWN_STDERR_TO_DEV_NULL,
		      NULL,    /* child_setup */
		      NULL,    /* user_data */
		      &sout,   /* standard output */
		      NULL,    /* standard error */
		      &status, /* exit status */
		      &error);

    if (!ok) {
	errbox(error->message);
	g_error_free(error);
	ret = LATEX_EXEC_FAILED;
    } else if (status != 0) {
	gchar *errmsg;

	errmsg = g_strdup_printf("%s\n%s", 
				 _("Failed to process TeX file"),
				 sout);
	errbox(errmsg);
	g_free(errmsg);
	ret = 1;
    }

    if (sout != NULL) g_free(sout);

    return ret;
}

#endif /* Windows or not */

int dvips_compile (char *texshort)
{
#ifdef G_OS_WIN32
    static char dvips_path[MAXLEN];
    char tmp[MAXLEN];
#endif
    int err = 0;

#if defined(G_OS_WIN32)
    if (*dvips_path == 0 && get_dvips_path(dvips_path)) {
	win_show_last_error();
	return 1;
    }

    sprintf(tmp, "\"%s\" -o %s.ps %s", dvips_path, texshort, texshort);
    if (win_run_sync(tmp, gretl_dotdir())) {
	return 1;
    }
#else
    err = spawn_dvips(texshort);
#endif 

    return err;
}

static int latex_compile_graph_page (void)
{
    char *fname;
    int err;

    /* ensure we don't get stale output */
    if (gpage.term == GP_TERM_EPS) {
	fname = gpage_fname(".ps", 0);
    } else {
	fname = gpage_fname(".pdf", 0);
    }

    gretl_remove(fname);

    err = latex_compile(gpage_tex_base);

    if (err == LATEX_ERROR) {
	char *fname = gpage_fname(".log", 0);

	view_file(fname, 0, 1, 78, 350, VIEW_FILE);
    }

    if (gpage.term == GP_TERM_EPS && !err) {
	err = dvips_compile(gpage_base);
    }    

    return err;
}

static int make_gp_output (void)
{
    char *fname;
    double scale = 1.0;
    int i;
    int err = 0;

    if (gpage.ngraphs == 3) {
	scale = 0.8;
    } else if (gpage.ngraphs > 3) {
	scale = 0.75;
    }

    for (i=0; i<gpage.ngraphs && !err; i++) {
	err = gp_make_outfile(gpage.fnames[i], i + 1, scale);
    }

    fname = gpage_fname(".plt", 0);
    gretl_remove(fname);

    return err;
}

static int real_display_gpage (void)
{
#ifndef G_OS_WIN32
    const char *viewer;
#endif
    char *fname;
    int err = 0;

    if (gpage.term == GP_TERM_PDF) {
	fname = gpage_fname(".pdf", 0);
    } else {
	fname = gpage_fname(".ps", 0);
    }

#if defined(G_OS_WIN32)
    err = win32_open_file(fname);
#elif defined(OS_OSX)
    err = osx_open_file(fname);
#else
    viewer = (gpage.term == GP_TERM_PDF)? "viewpdf" : "viewps";
    err = gretl_fork(viewer, fname);
#endif

    return err;
}

static void gpage_cleanup (void)
{
    char *fname;
    int i;

    for (i=0; i<gpage.ngraphs; i++) {
	if (gpage.term == GP_TERM_PDF) {
	    fname = gpage_fname(".pdf", i + 1);
	} else {
	    fname = gpage_fname(".ps", i + 1);
	}
	gretl_remove(fname);
    }

    fname = gpage_fname(".tex", 0);
    gretl_remove(fname);
    fname = gpage_fname(".dvi", 0);
    gretl_remove(fname);
    fname = gpage_fname(".log", 0);
    gretl_remove(fname);
    fname = gpage_fname(".aux", 0);
    gretl_remove(fname);
}

static gchar *gpage_switch_compiler (int term)
{
    gchar *orig = g_strdup(latex);
    char *p, tmp[MAXSTR];
    char test[4];
    int len0, have_pdf;

    strcpy(tmp, latex);
    p = strrchr(tmp, SLASH);
    if (p == NULL) {
	len0 = 0;
	p = tmp;
    } else {
	len0 = p - tmp + 1;
	p++;
    }

    *test = '\0';
    strncat(test, p, 3);
    gretl_lower(test);

    have_pdf = (strcmp(test, "pdf") == 0);

    if (term == GP_TERM_PDF && !have_pdf) {
	/* switch to pdflatex */
	*latex = '\0';
	strncat(latex, tmp, len0);
	strncat(latex, "pdf", 3);
	strcat(latex, p);
    } else if (term != GP_TERM_PDF && have_pdf) {
	/* switch to plain latex */
	*latex = '\0';
	strncat(latex, tmp, len0);
	strcat(latex, p + 3);
    }

    return orig;
}

static void gpage_revert_compiler (gchar *orig)
{
    if (orig != NULL) {
	strcpy(latex, orig);
	g_free(orig);
    }
}

int display_graph_page (GtkWidget *parent)
{
    const char *opts[] = {
	N_("color"),
	N_("monochrome")
    };
    const char *sdir = get_session_dirname();
    gchar *latex_orig = NULL;
    int resp, err = 0;

    if (gpage.ngraphs == 0) {
	gpage_errmsg(_("The graph page is empty"), 1);
	return 1;
    }

    resp = radio_dialog(_("graph page options"), NULL, opts, 
			2, 0, 0, parent);
    if (resp < 0) {
	return 0;
    }

    gpage.mono = (resp == 1);

    gretl_chdir(sdir);

    gpage_filenames_init(NULL);

    if (get_tex_use_pdf()) {
	if (gnuplot_pdf_terminal()) {
	    gpage.term = GP_TERM_PDF;
	} else {
	    gpage.term = GP_TERM_EPS;
	    latex_orig = gpage_switch_compiler(gpage.term);
	}
    } else {
	gpage.term = GP_TERM_EPS;
    }

    /* write the LaTeX driver file */
    err = make_graphpage_tex();

    if (!err) {
	/* transform individual plot files and compile 
	   using gnuplot */
	err = make_gp_output();
    }

    if (!err) {
	/* compile LaTeX */
	err = latex_compile_graph_page();
    }

    if (!err) {
	/* display output file */
	err = real_display_gpage();
    }

    gpage_revert_compiler(latex_orig);
    gpage_cleanup();

    if (err) {
	gui_errmsg(err);
    }

    return err;
}

void clear_graph_page (int on_exit)
{
    int i;

    if (!on_exit && gpage.ngraphs > 0) {
	mark_session_changed();
    }

    for (i=0; i<gpage.ngraphs; i++) {
	free(gpage.fnames[i]);
    }
    free(gpage.fnames);

    graph_page_init();
}

int graph_page_get_n_graphs (void)
{
    return gpage.ngraphs;
}

int save_graph_page (const char *fname)
{    
    const char *sdir = get_session_dirname();
    char *latex_orig = NULL;
    int err = 0;

    gretl_chdir(sdir);

    gpage_filenames_init(fname);

    if (get_tex_use_pdf()) {
	if (gnuplot_pdf_terminal()) {
	    gpage.term = GP_TERM_PDF;
	} else {
	    gpage.term = GP_TERM_EPS;
	    latex_orig = gpage_switch_compiler(gpage.term);
	}
    } else {
	gpage.term = GP_TERM_EPS;
    }

    /* write the LaTeX driver file */
    err = make_graphpage_tex();

    if (!err) {
	/* transform individual plot files and compile 
	   using gnuplot */
	err = make_gp_output();
    }

    gpage_revert_compiler(latex_orig);

    return err;
}

static int print_graph_page_direct (const char *fname, 
				    gretlopt opt)
{
    char thisdir[MAXLEN];
    char *latex_orig = NULL;
    int err = 0;

    if (gpage.ngraphs == 0) {
	gpage_errmsg(_("The graph page is empty"), 1);
	return 1;
    }

    if (getcwd(thisdir, MAXLEN - 1) == NULL) {
	*thisdir = '\0';
    }

    gpage.mono = (opt & OPT_M) ? 1 : 0;

    gretl_chdir(get_session_dirname());

    gpage_filenames_init(NULL);

    if (has_suffix(fname, ".pdf")) {
	gpage.term = GP_TERM_PDF;
	if (!get_tex_use_pdf()) {
	    latex_orig = gpage_switch_compiler(gpage.term);
	}	
    } else {
	gpage.term = GP_TERM_EPS;
	if (get_tex_use_pdf()) {
	    latex_orig = gpage_switch_compiler(gpage.term);
	}
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
	char *output;

	if (gpage.term == GP_TERM_PDF) {
	    output = gpage_fname(".pdf", 0);
	} else {
	    output = gpage_fname(".ps", 0);
	}
	if (*thisdir != '\0') {
	    /* come back out of dotdir */
	    chdir(thisdir);
	}
	fname = gretl_maybe_switch_dir(fname);
	err = gretl_copy_file(output, fname);
	if (!err) {
	    fprintf(stderr, "graphpg: wrote %s\n", fname);
	}
	remove(output);
    }

    gpage_revert_compiler(latex_orig);
    gpage_cleanup();

    return err;
}

int graph_page_parse_line (const char *line, gretlopt opt)
{
    const char *outfile = NULL;
    char cmdword[16];
    int gotcmd = 0;
    int err = 0;

    if (sscanf(line, "%*s %15s", cmdword) == 1) {
	gotcmd = 1;
    }

    if (gotcmd && (opt & OPT_O)) {
	/* the --output option rules out the various
	   command words */
	return E_PARSE;
    } else if (opt & OPT_O) {
	/* --output="filename" */
	outfile = get_optval_string(GRAPHPG, OPT_O);
	if (outfile == NULL) {
	    return E_PARSE;
	} else {
	    return print_graph_page_direct(outfile, opt);
	}
    } else if (!gotcmd) {
	return E_PARSE;
    }

    if (!strcmp(cmdword, "add")) {
	err = graph_page_add_last_graph();
    } else if (!strcmp(cmdword, "show")) {
	err = display_graph_page(NULL);
    } else if (!strcmp(cmdword, "free")) {
	if (gpage.ngraphs > 0) {
	    clear_graph_page(0);
	}
    } else if (!strcmp(cmdword, "fontscale")) {
	const char *s = strstr(line, "fontscale");

	err = graph_page_set_font_scale(s + 9);
    }

    return err;
}
