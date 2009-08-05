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
#include "gretl_foreign.h"

#include <glib.h>

#ifdef USE_RLIB
# include "libset.h"
# include <Rinternals.h> /* for SEXP and friends */
#endif

#ifdef G_OS_WIN32
# include <windows.h>
#else
# include <signal.h>
#endif

#define FDEBUG 0

static char **foreign_lines;
static int foreign_started;
static int foreign_n_lines; 
static int foreign_lang;

/* foreign_opt may include OPT_D to send data, OPT_Q to operate
   quietly (don't display output from foreign program)
*/
static gretlopt foreign_opt;

/* "dotdir" filenames */
static gchar *gretl_dotdir;
static gchar *gretl_Rprofile;
static gchar *gretl_Rsrc;
static gchar *gretl_Rout;
static gchar *gretl_Oxprog;

enum {
    LANG_R = 1,
    LANG_OX
};

static void destroy_foreign (void)
{
    free_strings_array(foreign_lines, foreign_n_lines);
    foreign_lines = NULL;
    foreign_started = 0;
    foreign_n_lines = 0;
    foreign_opt = OPT_NONE;
}

static int set_foreign_lang (const char *lang, PRN *prn)
{
    int err = 0;

    if (g_ascii_strcasecmp(lang, "R") == 0) {
	foreign_lang = LANG_R;
    } else if (g_ascii_strcasecmp(lang, "ox") == 0) {
#ifdef USE_OX
	foreign_lang = LANG_OX;
#else
	pprintf(prn, "%s: not supported\n", lang);
	err = E_DATA;
#endif
    } else {
	pprintf(prn, "%s: unknown language\n", lang);
	err = E_DATA;
    }

    return err;
}

static const gchar *gretl_ox_filename (void)
{
    if (gretl_Oxprog == NULL) {
	const char *dotdir = gretl_dot_dir();

	gretl_Oxprog = g_strdup_printf("%sgretltmp.ox", dotdir);
    }

    return gretl_Oxprog;
}

static void make_gretl_R_names (void)
{
    static int done;

    if (!done) {
	gretl_dotdir = g_strdup(gretl_dot_dir());
#ifdef G_OS_WIN32
	slash_convert(gretl_dotdir, FROM_BACKSLASH);
#endif
	gretl_Rprofile = g_strdup_printf("%sgretl.Rprofile", gretl_dotdir);
	gretl_Rsrc = g_strdup_printf("%sRsrc", gretl_dotdir);
	gretl_Rout = g_strdup_printf("%sR.out", gretl_dotdir);
	done = 1;
    }
}

#ifdef G_OS_WIN32

static int lib_run_R_sync (gretlopt opt, PRN *prn)
{
    char Rterm[MAX_PATH];
    gchar *cmd;
    int err = 0;

    err = R_path_from_registry(Rterm, RTERM);
    if (err) {
	return E_EXTERNAL;
    }

    cmd = g_strdup_printf("\"%s\" --no-save --no-init-file --no-restore-data "
			  "--slave", Rterm);

    err = winfork(cmd, NULL, SW_SHOWMINIMIZED, CREATE_NEW_CONSOLE);

    if (!err && !(opt & OPT_Q)) {
	FILE *fp = gretl_fopen(gretl_Rout, "r");

	if (fp == NULL) {
	    err = E_FOPEN;
	} else {
	    char line[1024];

	    while (fgets(line, sizeof line, fp)) {
		pputs(prn, line);
	    }
	    fclose(fp);
	}
    }

    g_free(cmd);

    return err;
}

static int lib_run_ox_sync (gretlopt opt, PRN *prn)
{
    const char *path = gretl_ox_path();
    const char *fname = gretl_ox_filename();
    gchar *cmd, *sout = NULL;
    int err;

    cmd = g_strdup_printf("\"%s\" \"%s\"", path, fname);
    err = gretl_win32_grab_output(cmd, &sout);

    if (sout != NULL && *sout != '\0') {
	pputs(prn, sout);
    }

    g_free(sout);
    g_free(cmd);

    return err;
}

#else /* !G_OS_WIN32 */

static int lib_run_prog_sync (char **argv, gretlopt opt, PRN *prn)
{
    gchar *sout = NULL;
    gchar *errout = NULL;
    gint status = 0;
    GError *gerr = NULL;
    int err = 0;

    signal(SIGCHLD, SIG_DFL);

    g_spawn_sync(NULL, argv, NULL, G_SPAWN_SEARCH_PATH,
		 NULL, NULL, &sout, &errout,
		 &status, &gerr);

    if (gerr != NULL) {
	pprintf(prn, "%s\n", gerr->message); 
	g_error_free(gerr);
	err = 1;
    } else if (status != 0) {
	pprintf(prn, "%s exited with status %d", argv[0], status);
	if (sout != NULL && *sout != '\0') {
	    pputs(prn, sout);
	    pputc(prn, '\n');
	} 
	if (errout != NULL && *errout != '\0') {
	    pputs(prn, errout);
	    pputc(prn, '\n');
	}	
	err = 1;
    } else if (sout != NULL) {
	if (!(opt & OPT_Q)) {
	    /* with OPT_Q, don't print non-error output */
	    pputs(prn, sout);
	    pputc(prn, '\n');
	}
    } else {
	pprintf(prn, "%s: %s\n", argv[0], "Got no output");
	err = 1;
    }

    g_free(sout);
    g_free(errout);

    return err;
}

static int lib_run_R_sync (gretlopt opt, PRN *prn)
{
    char *argv[] = {
	"R", 
	"--no-save",
	"--no-init-file",
	"--no-restore-data",
	"--slave",
	NULL
    };

    return lib_run_prog_sync(argv, opt, prn);
}

static int lib_run_ox_sync (gretlopt opt, PRN *prn)
{
    char *argv[3];
    int err;

    argv[0] = (char *) gretl_ox_path();
    argv[1] = (char *) gretl_ox_filename();
    argv[2] = NULL;

    err = lib_run_prog_sync(argv, opt, prn);

    return err;
}

#endif /* switch on MS Windows or not */

static int write_ox_gretl_io_file (void)
{
    static int written;

    if (!written) {
	const char *dotdir = gretl_dot_dir();
	gchar *fname;
	FILE *fp;

	fname = g_strdup_printf("%sgretl_io.ox", dotdir);
	fp = gretl_fopen(fname, "w");
	g_free(fname);

	if (fp == NULL) {
	    return E_FOPEN;
	} else {
	    fputs("gretl_export (const X, const str)\n{\n", fp);
#ifdef G_OS_WIN32
            gchar *dotcpy = g_strdup(dotdir);
            charsub(dotcpy, '\\', '/');
            fprintf(fp, "  decl fname = \"%s\" ~ str;\n", dotcpy);
#else
	    fprintf(fp, "  decl fname = \"%s\" ~ str;\n", dotdir);
#endif
	    fputs("  decl fp = fopen(fname, \"w\");\n", fp);
	    fputs("  fprint(fp, \"%d %d\", rows(X), columns(X));\n", fp);
	    fputs("  fprint(fp, \"%.15g\", X);\n", fp);
	    fputs("  fclose(fp);\n}\n\n", fp);

	    fputs("gretl_loadmat (const str)\n{\n", fp);
#ifdef G_OS_WIN32
            fprintf(fp, "  decl fname = \"%s\" ~ str;\n", dotcpy);
            g_free(dotcpy);
#else
	    fprintf(fp, "  decl fname = \"%s\" ~ str;\n", dotdir);
#endif
	    fputs("  decl X = loadmat(fname);\n", fp);
	    fputs("  return X;\n}\n", fp);

	    fclose(fp);
	    written = 1;
	}
    }

    return 0;
}

static void add_gretl_include (FILE *fp)
{
    const char *dotdir = gretl_dot_dir();

    if (strchr(dotdir, ' ')) {
	fprintf(fp, "#include \"%sgretl_io.ox\"\n", dotdir);
    } else {
	fprintf(fp, "#include <%sgretl_io.ox>\n", dotdir);
    }
}

/**
 * write_gretl_ox_file:
 * @buf: text buffer containing Ox code.
 * @opt: should contain %OPT_G for use from GUI.
 * @fname: location to receive name of file written, or %NULL.
 *
 * Writes the content of @buf into a file in the gretl user's
 * "dotdir".
 *
 * Returns: 0 on success, non-zero on error.
 */

int write_gretl_ox_file (const char *buf, gretlopt opt, const char **pfname)
{
    const gchar *fname = gretl_ox_filename();
    FILE *fp;
    int err;

    err = write_ox_gretl_io_file();

    fp = gretl_fopen(fname, "w");

    if (fp == NULL) {
	return E_FOPEN;
    } else {
	if (buf != NULL) {
	    /* pass on the material supplied in the 'buf' argument */
	    char line[1024];

	    bufgets_init(buf);
	    while (bufgets(line, sizeof line, buf)) {
		fputs(line, fp);
		if (!err && strstr(line, "oxstd.h")) {
		    add_gretl_include(fp);
		}
	    }
	    bufgets_finalize(buf);
	} else {
	    /* put out the stored 'foreign' lines */
	    int i;

	    for (i=0; i<foreign_n_lines; i++) { 
		fprintf(fp, "%s\n", foreign_lines[i]);
		if (!err && strstr(foreign_lines[i], "oxstd.h")) {
		    add_gretl_include(fp);
		}
	    }
	}
	fclose(fp);
	if (pfname != NULL) {
	    *pfname = fname;
	}
    }

    return 0;
}

/* write out current dataset in R format, and, if this succeeds,
   write appropriate R commands to @fp to source the data
*/

static int write_data_for_R (const double **Z, 
			     const DATAINFO *pdinfo,
			     FILE *fp)
{
    gchar *Rdata, *Rline;
    int err;

    Rdata = g_strdup_printf("%sRdata.tmp", gretl_dotdir);
    Rline = g_strdup_printf("store \"%s\" -r", Rdata);
    g_free(Rline);

    err = write_data(Rdata, NULL, Z, pdinfo, OPT_R, NULL);
    if (err) {
	g_free(Rdata);
	return err;
    }

    fputs("# load data from gretl\n", fp);
    fprintf(fp, "gretldata <- read.table(\"%s\", header=TRUE)\n", Rdata);
    g_free(Rdata);

    if (dataset_is_time_series(pdinfo)) {
	char *p, datestr[OBSLEN];
	int subper = 1;
	    
	ntodate(datestr, pdinfo->t1, pdinfo);
	p = strchr(datestr, ':');
	if (p != NULL) {
	    subper = atoi(p + 1);
	}
	    
	fprintf(fp, "gretldata <- ts(gretldata, start=c(%d, %d), frequency = %d)\n", 
		atoi(datestr), subper, pdinfo->pd);
    } else {
	fputs("attach(gretldata)\n", fp);
    }

    return err;
}

/* define an R function for passing data back to gretl */

static void write_R_export_func (FILE *fp) 
{
    fprintf(fp, "gretl.dotdir <- \"%s\"\n", gretl_dotdir);
    fputs("gretl.export <- function(x) {\n", fp);
    fprintf(fp, "  prefix <- \"%s\"\n", gretl_dotdir);
    fputs("  sx <- as.character(substitute(x))\n", fp);
    fputs("  if (is.ts(x)) {\n", fp);
    fputs("    fname <- paste(prefix, sx, \".csv\", sep=\"\")\n", fp);
    fputs("    dfx <- data.frame(x)\n", fp);
    fputs("    if (ncol(dfx) == 1) {\n", fp);
    fputs("      colnames(dfx) <- sx;\n", fp);
    fputs("    }\n", fp);
    fputs("    write.csv(dfx, file=fname, row.names=F)\n", fp);
    fputs("  } else if (is.data.frame(x)) {\n", fp);
    fputs("    fname <- paste(prefix, sx, \".csv\", sep=\"\")\n", fp);
    fputs("    write.csv(x, file=fname, row.names=F)\n", fp);
    fputs("  } else if (is.matrix(x)) {\n", fp);
    fputs("    fname <- paste(prefix, sx, \".mat\", sep=\"\")\n", fp);
    fputs("    write(dim(x), fname)\n", fp);
    fputs("    write(t(x), file=fname, ncolumns=ncol(x), append=TRUE)\n", fp);
    fputs("  }\n", fp);
    fputs("}\n", fp);
}

/* basic content which can either go into gretl.Rprofile or into
   Rsrc for sourcing */

static void put_startup_content (FILE *fp)
{
    fputs("vnum <- as.double(R.version$major) + (as.double(R.version$minor) / 10.0)\n", 
	  fp);
    fputs("if (vnum > 2.41) library(utils)\n", fp);
    fputs("library(stats)\n", fp);
    fputs("if (vnum <= 1.89) library(ts)\n", fp);
    write_R_export_func(fp);
}

/* Set up a gretl-specific R profile, and put notice of its existence
   into the environment.  Used when exec'ing the R binary (only) */

static int write_gretl_R_profile (gretlopt opt)
{
    FILE *fp;
    int err;

#if FDEBUG
    printf("writing R profile: starting\n");
#endif

    err = gretl_setenv("R_PROFILE", gretl_Rprofile);
    if (err) {
	return err;
    }     

    fp = gretl_fopen(gretl_Rprofile, "w");

    if (fp == NULL) {
	err = E_FOPEN;
    } else {
	put_startup_content(fp);
	fprintf(fp, "source(\"%s\", %s = TRUE)\n", 
		gretl_Rsrc, (opt & OPT_V)? "echo" : "print.eval");
	fclose(fp);
    }

#if FDEBUG
    printf("writing R profile: returning %d\n", err);
#endif

    return err;
}

/* Write an R command file to be sourced by R.  @buf may contain R
   commands assembled via the GUI; if it is NULL the current "foreign"
   block is used as input.

   @opt may contain the following:

   OPT_I: indicates that we're in the context of an interactive R
   session.

   OPT_D: indicates that the current gretl dataset should be sent
   to R.

   OPT_G: we're being called via the gretl GUI.

   OPT_L: indicates that the source file is intended for use
   via the R shared library.
*/

static int write_R_source_file (const char *buf,
				const double **Z, 
				const DATAINFO *pdinfo,
				gretlopt opt)
{
    FILE *fp = gretl_fopen(gretl_Rsrc, "w");
    int i, err = 0;

#if FDEBUG
    printf("write R source file: starting\n");
#endif

    if (fp == NULL) {
	err = E_FOPEN;
    } else {
#ifdef G_OS_WIN32
	if (!(opt & OPT_I)) {
	    /* non-interactive */
	    fprintf(fp, "sink(\"%s\")\n", gretl_Rout);
	}
#endif

	if (opt & OPT_L) {
	    /* we're using the R shared library */
	    static int startup_done;

	    if (!startup_done) {
#if FDEBUG
		printf("Rlib: writing 'startup' material\n");
#endif
		put_startup_content(fp);
		startup_done = 1;
	    }
	    fprintf(fp, "sink(\"%s\")\n", gretl_Rout);
	}

	if (opt & OPT_D) {
	    /* send data */
	    err = write_data_for_R(Z, pdinfo, fp);
	}

	if (buf != NULL) {
	    /* pass on the script supplied in the 'buf' argument */
	    fputs("# load script from gretl\n", fp);
	    fputs(buf, fp);
	} else if (!(opt & OPT_G)) {
	    /* non-GUI */
	    for (i=0; i<foreign_n_lines; i++) { 
		fprintf(fp, "%s\n", foreign_lines[i]);
	    }
	}

	if (opt & OPT_L) {
	    fputs("sink()\n", fp);
	}

#ifdef G_OS_WIN32
	if (!(opt & OPT_I) && !(opt & OPT_L)) {
	    /* Rterm on Windows won't exit without this? */
	    fputs("q()\n", fp);
	}
#endif

	fclose(fp);
    }

#if FDEBUG
    printf("write R source file: returning %d\n", err);
#endif

    return err;
}

/* Write files to be read by R: profile to be read on startup and
   command source file.  This is called when we're exec'ing the R
   binary.  OPT_G in @opt indicates that this function is being called
   from the GUI program; @buf may contain R commands taken from a GUI
   window.
*/

int write_gretl_R_files (const char *buf,
			 const double **Z, 
			 const DATAINFO *pdinfo,
			 gretlopt opt)
{ 
    int err = 0;

#if FDEBUG
    printf("write_gretl_R_files: starting\n");
#endif

    make_gretl_R_names();

    /* write a temporary R profile so R knows what to do */
    err = write_gretl_R_profile(opt);
    if (err) {
	fprintf(stderr, "error writing gretl.Rprofile\n");
    } 

    if (!err) {
	/* write commands and/or data to file, to be sourced in R */
	err = write_R_source_file(buf, Z, pdinfo, opt);
	if (err) {
	    fprintf(stderr, "error writing gretl's Rsrc\n");
	} 	
    }

#if FDEBUG
    printf("write_gretl_R_files: returning %d\n", err);
#endif

    return err;
}

void delete_gretl_R_files (void)
{
#if FDEBUG
    printf("deleting gretl R files...\n");
#endif

    if (gretl_Rprofile != NULL) {
	gretl_remove(gretl_Rprofile);
    }
    if (gretl_Rsrc != NULL) {
	gretl_remove(gretl_Rsrc);
    }
}

void delete_gretl_ox_file (void)
{
    if (gretl_Oxprog != NULL) {
	gretl_remove(gretl_Oxprog);
    }
}

/* The following code block is used if we're implementing
   gretl's R support by dlopening the R shared library
   (as opposed to executing the R binary).
*/

#ifdef USE_RLIB

static void *Rhandle;  /* handle to the R library */
static int Rlib_err;   /* initialization error record */
static int Rinit;      /* are we initialized or not? */

static SEXP current_arg;
static SEXP current_call;

/* ponters to, and renamed versions of, the R global variables
   we'll need */

SEXP *PR_GlobalEnv;
SEXP *PR_NilValue;
SEXP *PR_UnboundValue;

SEXP VR_GlobalEnv;
SEXP VR_NilValue;
SEXP VR_UnboundValue;

/* renamed, pointerized versions of the R functions we need */

static double *(*R_REAL) (SEXP);

static SEXP (*R_CDR) (SEXP);
static SEXP (*R_allocList) (int);
static SEXP (*R_allocMatrix) (SEXPTYPE, int, int);
static SEXP (*R_allocVector) (SEXPTYPE, R_len_t);
static SEXP (*R_findFun) (SEXP, SEXP);
static SEXP (*R_findVar) (SEXP, SEXP);
static SEXP (*R_SETCAR) (SEXP, SEXP);
static SEXP (*R_protect) (SEXP);
static SEXP (*R_ScalarReal) (double);
static SEXP (*R_catch) (SEXP, SEXP, int *);
static SEXP (*R_install) (const char *);
static SEXP (*R_mkString) (const char *);

static Rboolean (*R_isMatrix) (SEXP);
static Rboolean (*R_isLogical) (SEXP);
static Rboolean (*R_isInteger) (SEXP);
static Rboolean (*R_isReal) (SEXP);

static int (*R_initEmbeddedR) (int, char **);
static int (*R_ncols) (SEXP);
static int (*R_nrows) (SEXP);
static int (*R_TYPEOF) (SEXP);

static void (*R_endEmbeddedR) (int);
static void (*R_unprotect) (int);
static void (*R_PrintValue) (SEXP);
static void (*R_SET_TYPEOF) (SEXP, int);
static void (*R_SET_TAG) (SEXP, SEXP);

static int *(*R_LOGICAL) (SEXP);

#ifdef WIN32
static char *(*R_get_HOME) (void);
#endif

/* utility function to cumulate errors from dlsym */

static void *dlget (void *handle, const char *name, int *err)
{
    void *p = gretl_dlsym(handle, name);
    
    if (p == NULL) {
	fprintf(stderr, "dlget: couldn't find '%s'\n", name);
	*err += 1;
    }

    return p;
}

/* dlopen the R library and grab all the symbols we need:
   several function pointers and a few global variables
*/

static int load_R_symbols (void)
{
    const char *libpath = gretl_rlib_path();
    int err = 0;

#if FDEBUG
    printf("Loading libR symbols\n");
#endif

    Rhandle = gretl_dlopen(libpath, 1);
    if (Rhandle == NULL) {
	err = E_EXTERNAL;
	goto bailout;
    } 

    R_CDR           = dlget(Rhandle, "CDR", &err);
    R_REAL          = dlget(Rhandle, "REAL", &err);
    R_allocList     = dlget(Rhandle, "Rf_allocList", &err);
    R_allocMatrix   = dlget(Rhandle, "Rf_allocMatrix", &err);
    R_allocVector   = dlget(Rhandle, "Rf_allocVector", &err);
    R_endEmbeddedR  = dlget(Rhandle, "Rf_endEmbeddedR", &err);
    R_findFun       = dlget(Rhandle, "Rf_findFun", &err);
    R_findVar       = dlget(Rhandle, "Rf_findVar", &err);
    R_initEmbeddedR = dlget(Rhandle, "Rf_initEmbeddedR", &err);
    R_install       = dlget(Rhandle, "Rf_install", &err);
    R_isMatrix      = dlget(Rhandle, "Rf_isMatrix", &err);
    R_isLogical     = dlget(Rhandle, "Rf_isLogical", &err);
    R_isInteger     = dlget(Rhandle, "Rf_isInteger", &err);
    R_isReal        = dlget(Rhandle, "Rf_isReal", &err);
    R_mkString      = dlget(Rhandle, "Rf_mkString", &err);
    R_ncols         = dlget(Rhandle, "Rf_ncols", &err);
    R_nrows         = dlget(Rhandle, "Rf_nrows", &err);
    R_PrintValue    = dlget(Rhandle, "Rf_PrintValue", &err);
    R_protect       = dlget(Rhandle, "Rf_protect", &err);
    R_ScalarReal    = dlget(Rhandle, "Rf_ScalarReal", &err);
    R_unprotect     = dlget(Rhandle, "Rf_unprotect", &err);
    R_catch         = dlget(Rhandle, "R_tryEval", &err);
    R_SETCAR        = dlget(Rhandle, "SETCAR", &err);
    R_SET_TYPEOF    = dlget(Rhandle, "SET_TYPEOF", &err); 
    R_TYPEOF        = dlget(Rhandle, "TYPEOF", &err);
    R_SET_TAG       = dlget(Rhandle, "SET_TAG", &err);
    R_LOGICAL       = dlget(Rhandle, "LOGICAL", &err);

#ifdef WIN32
    R_get_HOME = dlget(Rhandle, "get_R_HOME", &err);
#endif

    if (!err) {
	PR_GlobalEnv    = (SEXP *) dlget(Rhandle, "R_GlobalEnv", &err);
	PR_NilValue     = (SEXP *) dlget(Rhandle, "R_NilValue", &err);
	PR_UnboundValue = (SEXP *) dlget(Rhandle, "R_UnboundValue", &err);
    }

    if (err) {
	close_plugin(Rhandle);
	Rhandle = NULL;
	err = E_EXTERNAL;
    }

 bailout:

#if FDEBUG
    printf("load_R_symbols: returning %d\n", err);
#endif

    return err;
}

void gretl_R_cleanup (void)
{
#if FDEBUG
    printf("gretl_R_cleanup: Rinit = %d\n", Rinit);
#endif

    if (Rinit) {
	R_endEmbeddedR(0);
	close_plugin(Rhandle);
	Rhandle = NULL;
    }
}

/* called from gretl_paths.c on revising the Rlib path:
   allow for the possibility that the path was wrong but is
   now OK
*/

void gretl_R_reset_error (void)
{
    Rlib_err = 0;
}

/* Initialize the R library for use with gretl.  Note that we only
   need do this once per gretl session.  We need to check that the
   environment is set to R's liking first, otherwise initialization
   will fail -- and abort gretl too!
*/

static int gretl_Rlib_init (void)
{
    char *Rhome;
    int err = 0;

#if FDEBUG
    printf("gretl_Rlib_init: starting\n");
#endif

#ifndef WIN32
    Rhome = getenv("R_HOME");
    if (Rhome == NULL) {
	fprintf(stderr, "To use Rlib, the variable R_HOME must be set\n");
	err = E_EXTERNAL;
	goto bailout;
    }
#endif

    err = load_R_symbols();
    if (err) {
	fprintf(stderr, "gretl_Rlib_init: failed to load R functions\n");
	goto bailout;
    }

#ifdef WIN32
    Rhome = R_get_HOME();
    if (Rhome == NULL) {
	fprintf(stderr, "To use Rlib, the variable R_HOME must be set\n");
	err = E_EXTERNAL;
	goto bailout;
    }
#endif

    /* ensure common filenames are in place */
    make_gretl_R_names();

    /* and ensure that gretl.Rprofile doesn't get in the way */
    gretl_remove(gretl_Rprofile);

    if (!err) {
	char *argv[] = { 
	    "gretl", 
	    "--no-save", 
	    "--silent", 
	};
	int ok, argc = 3;

	ok = R_initEmbeddedR(argc, argv);
	if (ok) {
	    VR_GlobalEnv = *PR_GlobalEnv;
	    VR_NilValue = *PR_NilValue;
	    VR_UnboundValue = *PR_UnboundValue;
	    Rinit = 1;
	} else {
	    close_plugin(Rhandle);
	    Rhandle = NULL;
	    err = Rlib_err = E_EXTERNAL;
	}
    }

 bailout:

#if FDEBUG
    printf("gretl_Rlib_init: returning %d\n", err);
#endif

    return err;
}

/* run R's source() function on an R command file written by
   gretl -- shared library version */

static int lib_run_Rlib_sync (gretlopt opt, PRN *prn) 
{
    int err = 0;

#if FDEBUG
    printf("lib_run_Rlib_sync: starting\n");
#endif

    if (!Rinit) {
	err = gretl_Rlib_init();
    }

    if (!err) {
	SEXP expr, p;

	/* make echo/print.eval argument */
	R_protect(p = R_allocVector(LGLSXP, 1));
	R_LOGICAL(p)[0] = TRUE;

	/* expression source(f, print.eval=p) */
	R_protect(expr = R_allocVector(LANGSXP, 3));
	R_SETCAR(expr, R_install("source")); 
	R_SETCAR(R_CDR(expr), R_mkString(gretl_Rsrc));
	R_SETCAR(R_CDR(R_CDR(expr)), p);
	R_SET_TAG(R_CDR(R_CDR(expr)), 
		  R_install((opt & OPT_V)? "echo" : "print.eval"));

	R_catch(expr, NULL, &err);
	R_unprotect(2);
    }

    if (!err && prn != NULL) {
	FILE *fp = gretl_fopen(gretl_Rout, "r");

	if (fp != NULL) {
	    char line[512];

	    while (fgets(line, sizeof line, fp)) {
		pputs(prn, line);
	    }
	    fclose(fp);
	    gretl_remove(gretl_Rout);
	}
    }

#if FDEBUG
    printf("lib_run_Rlib_sync: returning %d\n", err);
#endif

    return (err)? E_EXTERNAL : 0;
}

static SEXP find_R_function (const char *name)
{
    SEXP fun;
    SEXPTYPE t;

    fun = R_findVar(R_install(name), VR_GlobalEnv);
    t = R_TYPEOF(fun);

    if (t == PROMSXP) {
	/* eval promise if need be */
	int err = 1;

	fun = R_catch(fun, VR_GlobalEnv, &err);
	if (!err) {
	    t = R_TYPEOF(fun);
	}
    }

    if (t != CLOSXP && t != BUILTINSXP &&
	t != BUILTINSXP && t != SPECIALSXP) {
	return VR_UnboundValue;
    }

    return fun;
}

/* Check if we should be using the R shared library for executing the
   code in a "foreign" block.  This is disabled if the user has not
   done "set R_lib on", and can be prohibited by the environment
   variable GRETL_NO_RLIB.  It may also be blocked if we already tried
   and failed to initialize the library for gretl's use.  (The
   fallback will be to call the R binary.)
*/

static int gretl_use_Rlib (void)
{
    int ret = 0;

#if FDEBUG
    printf("gretl_use_Rlib: starting\n");
#endif

    if (!Rlib_err && libset_get_bool(R_LIB) && !getenv("GRETL_NO_RLIB")) {
	/* use of library is not blocked */
	if (Rinit) {
	    /* already opened, fine */
	    ret = 1;
	} else {
	    /* try opening library */
	    Rlib_err = gretl_Rlib_init();
	    ret = !Rlib_err;
	}
    }

#if FDEBUG
    printf("gretl_use_Rlib: using %s\n", (ret)? "library" : "executable");
#endif    

    return ret;
}

/* Used in "genr", to see if @name denotes an R function,
   either built-in or possibly user-defined.  The lookup
   is conditional on the user's doing "set R_functions on".
*/

int get_R_function_by_name (const char *name) 
{
    int ret = 0;
    
    if (libset_get_bool(R_FUNCTIONS) && gretl_use_Rlib()) {
	SEXP fun = find_R_function(name);

	ret = (fun == VR_UnboundValue)? 0 : 1;
    } 

    return ret;
}

/* gretl_R_function_add... : these functions are used to convert from
   gretl types to R constructs for passing to R functions
*/

int gretl_R_function_add_scalar (double x) 
{
    current_arg = R_CDR(current_arg);
    R_SETCAR(current_arg, R_ScalarReal(x));
    return 0;
}

int gretl_R_function_add_vector (const double *x, int t1, int t2) 
{
    SEXP res = R_allocVector(REALSXP, t2 - t1 + 1);
    int i;

    if (res == NULL) {
	return E_ALLOC;
    }

    current_arg = R_CDR(current_arg);

    for (i=t1; i<=t2; i++) {
    	R_REAL(res)[i-t1] = x[i];
    }
    
    R_SETCAR(current_arg, res);
    return 0;
}

int gretl_R_function_add_matrix (const gretl_matrix *m) 
{
    int nr = gretl_matrix_rows(m);
    int nc = gretl_matrix_cols(m);
    SEXP res;
    int i, j;

    current_arg = R_CDR(current_arg);
    res = R_allocMatrix(REALSXP, nr, nc);
    if (res == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<nr; i++) {
	for (j=0; j<nc; j++) {
	    R_REAL(res)[i + j * nr] = gretl_matrix_get(m, i, j);
    	}
    }
    
    R_SETCAR(current_arg, res);

    return 0;
}

/* called from geneval.c only, and should be pre-checked */

int gretl_R_get_call (const char *name, int argc) 
{
    SEXP call, e;

    call = R_findFun(R_install(name), VR_GlobalEnv);

    if (call == VR_NilValue) {
	fprintf(stderr, "gretl_R_get_call: no definition for function %s\n", 
		name);
	R_unprotect(1); /* is this OK? */
	return E_EXTERNAL;
    } 

    R_protect(e = R_allocList(argc + 1));
    R_SET_TYPEOF(e, LANGSXP);
    R_SETCAR(e, R_install(name));
    current_call = current_arg = e;
 
    return 0;
}

static int R_type_to_gretl_type (SEXP s)
{
    if (R_isMatrix(s)) {
	return GRETL_TYPE_MATRIX;
    } else if (R_isLogical(s)) {
	return GRETL_TYPE_BOOL;
    } else if (R_isInteger(s)) {
	return GRETL_TYPE_INT;
    } else if (R_isReal(s)) {
	return GRETL_TYPE_DOUBLE;
    } else {
	return GRETL_TYPE_NONE;
    }
}

/* execute an R function and try to convert the value returned
   into a gretl type
*/

int gretl_R_function_exec (const char *name, int *rtype, void **ret) 
{
    SEXP res;
    int err = 0;

    if (gretl_messages_on()) {
	R_PrintValue(current_call);
    }

    res = R_catch(current_call, VR_GlobalEnv, &err);
    if (err) {
	return E_EXTERNAL;
    }

    *rtype = R_type_to_gretl_type(res);

    if (*rtype == GRETL_TYPE_MATRIX) {
	gretl_matrix *pm;
	int nc = R_ncols(res);
	int nr = R_nrows(res);
	int i, j;

	pm = gretl_matrix_alloc(nr, nc);
	if (pm == NULL) {
	    return E_ALLOC;
	}

	for (i=0; i<nr; i++) {
	    for (j=0; j<nc; j++) {
		gretl_matrix_set(pm, i, j, R_REAL(res)[i + j * nr]);
	    }
	}
	R_unprotect(1);
	*ret = pm;
    } else if (gretl_scalar_type(*rtype)) {
	double *realres = R_REAL(res);
	double *dret = *ret;

	*dret = *realres;

    	R_unprotect(1);
    } else {
	err = E_EXTERNAL;
    }

    return err;
}

#endif /* USE_RLIB */

static int foreign_block_init (const char *line, gretlopt opt, PRN *prn)
{
    int err = 0;

    foreign_opt = OPT_NONE;

    if (!strncmp(line, "foreign ", 8)) {
	char lang[16];

	line += 8;
	line += strspn(line, " ");
	if (!strncmp(line, "language", 8)) {
	    line += 8;
	    line += strspn(line, " =");
	    if (sscanf(line, "%15s", lang) == 1) {
		err = set_foreign_lang(lang, prn);
	    } else {
		err = E_PARSE;
	    }
	} else {
	    err = E_PARSE;
	}
    } else {
	/* we'll default to R for now */
	foreign_lang = LANG_R;
    }

    if (!err) {
	foreign_opt = opt;
    }

    return err;
}

/**
 * foreign_append_line:
 * @line: command line.
 * @opt: may include %OPT_V for verbose operation
 * @prn: struct for printing output.
 *
 * Appends @line to an internally stored block of "foreign"
 * commands, or starts a new block if no such block is
 * currently defined.
 * 
 * Returns: 0 on success, non-zero on error.
 */

int foreign_append_line (const char *line, gretlopt opt, PRN *prn)
{
    int err = 0;

    if (string_is_blank(line)) {
	return 0;
    }

    if (!foreign_started) {
	/* starting from scratch */
	err = foreign_block_init(line, opt, prn);
	if (!err) {
	    foreign_started = 1;
	}
    } else {
	/* appending */
	err = strings_array_add(&foreign_lines, &foreign_n_lines, line);
	if (err) {
	    destroy_foreign();
	}
    }	

    return err;
}

static int run_R_binary (const double **Z, const DATAINFO *pdinfo, 
			 gretlopt opt, PRN *prn)
{
    int err;

    /* write both profile and Rsrc files */

    err = write_gretl_R_files(NULL, Z, pdinfo, opt);
    if (err) {
	delete_gretl_R_files();
    } else {
	err = lib_run_R_sync(opt, prn);
    }

    return err;
}

static int run_R_lib (const double **Z, const DATAINFO *pdinfo, 
		      gretlopt opt, PRN *prn)
{
    int err;

    /* we don't want gretl.Rprofile in the way */
    gretl_remove(gretl_Rprofile);

    /* by passing OPT_L below we indicate that we're
       using the library */
    err = write_R_source_file(NULL, Z, pdinfo, opt | OPT_L);
    if (!err) {
	err = lib_run_Rlib_sync(opt, prn);
    }

    return err;
}

/**
 * foreign_execute:
 * @Z: data array.
 * @pdinfo: dataset information.
 * @opt: may include %OPT_V for verbose operation
 * @prn: struct for printing output.
 *
 * Executes a block of commands previously established via
 * calls to foreign_append_line().
 * 
 * Returns: 0 on success, non-zero on error.
 */

int foreign_execute (const double **Z, const DATAINFO *pdinfo, 
		     gretlopt opt, PRN *prn)
{
    int i, err = 0;

    if (foreign_lang == LANG_R) {
	make_gretl_R_names();
    }

    if (opt & OPT_V) {
	/* verbose: echo the stored commands */
	for (i=0; i<foreign_n_lines; i++) {
	    pprintf(prn, "> %s\n", foreign_lines[i]);
	}
    }

    foreign_opt |= opt;

    if (foreign_lang == LANG_R) {
#ifndef USE_RLIB
	err = run_R_binary(Z, pdinfo, foreign_opt, prn);	
#else
	if (gretl_use_Rlib()) {
	    err = run_R_lib(Z, pdinfo, foreign_opt, prn);
	} else {
	    err = run_R_binary(Z, pdinfo, foreign_opt, prn);
	}
#endif	    
    } else if (foreign_lang == LANG_OX) {
	err = write_gretl_ox_file(NULL, foreign_opt, NULL);
	if (err) {
	    delete_gretl_ox_file();
	} else {
	    err = lib_run_ox_sync(foreign_opt, prn);
	}
    } else {
	/* "can't happen" */
	err = E_DATA;
    }
    
    destroy_foreign();

    return err;
}

