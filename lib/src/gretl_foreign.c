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
# define STRICT_R_HEADERS
# include <R.h>
# include <Rinternals.h>
# include <Rembedded.h> 
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
static gretlopt foreign_opt;

/* foreign_opt may include OPT_D to send data, OPT_Q to operate
   quietly (don't display output from foreign program)
*/

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

    if (!strcmp(lang, "R")) {
	foreign_lang = LANG_R;
    } else if (!strcmp(lang, "ox")) {
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

static gchar *gretl_ox_filename (void)
{
    const char *dotdir = gretl_dot_dir();
    
    return g_strdup_printf("%sgretltmp.ox", dotdir);
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
	gchar *Rout = g_strdup_printf("%sR.out", gretl_dot_dir());
	FILE *fp;

	fp = gretl_fopen(Rout, "r");
	if (fp == NULL) {
	    err = E_FOPEN;
	} else {
	    char line[1024];

	    while (fgets(line, sizeof line, fp)) {
		pputs(prn, line);
	    }
	    fclose(fp);
	}
	g_free(Rout);
    }

    g_free(cmd);

    return err;
}

static int lib_run_ox_sync (gretlopt opt, PRN *prn)
{
    return E_EXTERNAL;
}

#else

static int lib_run_prog_sync (char **argv, gretlopt opt, PRN *prn)
{
    gchar *out = NULL;
    gchar *errout = NULL;
    gint status = 0;
    GError *gerr = NULL;
    int err = 0;

    signal(SIGCHLD, SIG_DFL);

    g_spawn_sync(NULL, argv, NULL, G_SPAWN_SEARCH_PATH,
		 NULL, NULL, &out, &errout,
		 &status, &gerr);

    if (gerr != NULL) {
	pprintf(prn, "%s\n", gerr->message); 
	g_error_free(gerr);
	err = 1;
    } else if (status != 0) {
	if (errout != NULL) {
	    if (*errout == '\0') {
		pprintf(prn, "%s exited with status %d", argv[0], status);
	    } else {
		pputs(prn, errout);
		pputc(prn, '\n');
	    } 
	}
	err = 1;
    } else if (out != NULL) {
	if (!(opt & OPT_Q)) {
	    /* with OPT_Q, don't print non-error output */
	    pputs(prn, out);
	    pputc(prn, '\n');
	}
    } else {
	pprintf(prn, "%s\n", "Got no output");
	err = 1;
    }

    g_free(out);
    g_free(errout);

    return err;
}

static int lib_run_R_sync (gretlopt opt, PRN *prn)
{
    gchar *argv[] = {
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
    gchar *fname;
    gchar *argv[] = {
	"oxl", 
	NULL,
	NULL
    };
    int err;

    fname = gretl_ox_filename();
    argv[1] = fname;
    err = lib_run_prog_sync(argv, opt, prn);
    g_free(fname);

    return err;
}

#endif

int write_gretl_ox_file (const char *buf, gretlopt opt)
{
    gchar *fname = gretl_ox_filename();
    FILE *fp;
    int i, err = 0;

    fp = gretl_fopen(fname, "w");

    if (fp == NULL) {
	err = E_FOPEN;
    } else {
	if (buf != NULL) {
	    /* pass on the material supplied in the 'buf' argument */
	    fputs(buf, fp);
	} else if (!(opt & OPT_G)) {
	    /* non-GUI */
	    for (i=0; i<foreign_n_lines; i++) { 
		fprintf(fp, "%s\n", foreign_lines[i]);
	    }
	}
	fclose(fp);
    }

    g_free(fname);

    return err;
}

/* write out current dataset in R format, and, if this succeeds,
   write appropriate R commands to @fp to source the data
*/

static int write_data_for_R (const double **Z, 
			     const DATAINFO *pdinfo,
			     const gchar *dotdir,
			     FILE *fp)
{
    gchar *Rline = NULL, *Rdata = NULL;
    int err;

    Rdata = g_strdup_printf("%sRdata.tmp", dotdir);
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
	    
	ntodate_full(datestr, pdinfo->t1, pdinfo);
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

/* define an R function for passing data or matrices back to gretl */

static void write_R_export_func (const gchar *dotdir, FILE *fp) 
{
    fprintf(fp, "gretl.dotdir <- \"%s\"\n", dotdir);
    fputs("gretl.export <- function(x) {\n", fp);
    fprintf(fp, "  prefix <- \"%s\"\n", dotdir);
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

/* set up for a temporary R profile, so R knows what to do */

static int write_gretl_R_profile (const char *Rprofile, const char *Rsrc,
				  const gchar *dotdir, gretlopt opt)
{
    FILE *fp;
    int err;

#ifdef G_OS_WIN32
    err = !SetEnvironmentVariable("R_PROFILE", Rprofile);
#else
    err = setenv("R_PROFILE", Rprofile, 1);
#endif
    if (err) {
	return err;
    }     

    fp = gretl_fopen(Rprofile, "w");
    if (fp == NULL) {
	return E_FOPEN;
    }

    /* profile preamble */
    fputs("vnum <- as.double(R.version$major) + (as.double(R.version$minor) / 10.0)\n", 
	  fp);
    fputs("if (vnum > 2.41) library(utils)\n", fp);
    fputs("library(stats)\n", fp);
    fputs("if (vnum <= 1.89) library(ts)\n", fp);

    write_R_export_func(dotdir, fp);

    fprintf(fp, "source(\"%s\", echo=TRUE)\n", Rsrc);

    fclose(fp);

    return 0;
}

static int write_R_source_file (const char *Rsrc, const char *buf,
				const double **Z, const DATAINFO *pdinfo,
				const gchar *dotdir, gretlopt opt)
{
    FILE *fp = gretl_fopen(Rsrc, "w");
    int i, err = 0;

    if (fp == NULL) {
	err = E_FOPEN;
    } else {
#ifdef G_OS_WIN32
	if (!(opt & OPT_I)) {
	    /* non-interactive */
	    gchar *Rout = g_strdup_printf("%sR.out", dotdir);

	    fprintf(fp, "sink(\"%s\")\n", Rout);
	    g_free(Rout);
	}
#endif

	if (opt & OPT_D) {
	    /* send data */
	    err = write_data_for_R(Z, pdinfo, dotdir, fp);
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

#ifdef G_OS_WIN32
	if (!(opt & OPT_I)) {
	    /* Rterm on Windows won't exit without this? */
	    fputs("q()\n", fp);
	}
#endif

	fclose(fp);
    }

    return err;
}

/* opt should contain OPT_G (only) when called from the
   GUI program */

int write_gretl_R_files (const char *buf,
			 const double **Z, const DATAINFO *pdinfo,
			 gretlopt opt)
{ 
    const char *dotdir = gretl_dot_dir();
    gchar *dotcpy = g_strdup(dotdir);
    gchar *Rprofile, *Rsrc;
    int err = 0;

#ifdef G_OS_WIN32
    slash_convert(dotcpy, FROM_BACKSLASH);
#endif

    /* organize filenames */
    Rprofile = g_strdup_printf("%sgretl.Rprofile", dotcpy);
    Rsrc = g_strdup_printf("%sRsrc", dotcpy);
    
    /* first write a temporary R profile, so R knows what to do */
    err = write_gretl_R_profile(Rprofile, Rsrc, dotcpy, opt);
    if (err) {
	fprintf(stderr, "error writing gretl Rprofile\n");
	return err;
    }

    /* then write commands and/or data to file, to be sourced
       by the R profile */
    err = write_R_source_file(Rsrc, buf, Z, pdinfo, dotcpy, opt);

    g_free(dotcpy);
    g_free(Rprofile);
    g_free(Rsrc);

    return err;
}

void delete_gretl_R_files (void)
{
    const char *dotdir = gretl_dot_dir();
    gchar *Rprofile, *Rsrc;
    
    Rprofile = g_strdup_printf("%sgretl.Rprofile", dotdir);
    Rsrc = g_strdup_printf("%sRsrc", dotdir);

    gretl_remove(Rprofile);
    gretl_remove(Rsrc);

    g_free(Rprofile);
    g_free(Rsrc);
}

void delete_gretl_ox_file (void)
{
    gchar *fname = gretl_ox_filename();

    gretl_remove(fname);
    g_free(fname);
}

#ifdef USE_RLIB

static void *Rhandle;

SEXP *PVR_GlobalEnv;
SEXP *PVR_NilValue;
SEXP *PVR_UnboundValue;

SEXP VR_GlobalEnv;
SEXP VR_NilValue;
SEXP VR_UnboundValue;

static double *(*VREAL) (SEXP);

static SEXP (*VCDR) (SEXP);
static SEXP (*VRf_allocList) (int);
static SEXP (*VRf_allocMatrix) (SEXPTYPE, int, int);
static SEXP (*VRf_allocVector) (SEXPTYPE, R_len_t);
static SEXP (*VRf_findFun) (SEXP, SEXP);
static SEXP (*VRf_findVar) (SEXP, SEXP);
static SEXP (*VSETCAR) (SEXP, SEXP);
static SEXP (*VRf_lang2) (SEXP, SEXP);
static SEXP (*VRf_protect) (SEXP);
static SEXP (*VRf_ScalarReal) (double);
static SEXP (*VR_tryEval) (SEXP, SEXP, int *);
static SEXP (*VRf_install) (const char *);
static SEXP (*VRf_mkString) (const char *);

static Rboolean (*VRf_isMatrix) (SEXP);
static Rboolean (*VRf_isLogical) (SEXP);
static Rboolean (*VRf_isInteger) (SEXP);
static Rboolean (*VRf_isReal) (SEXP);

static int (*VRf_initEmbeddedR) (int, char **);
static int (*VRf_ncols) (SEXP);
static int (*VRf_nrows) (SEXP);
static int (*VTYPEOF) (SEXP);

static void (*VRf_endEmbeddedR) (int);
static void (*VRf_unprotect) (int);
static void (*VRf_PrintValue) (SEXP);
static void (*VSET_TYPEOF) (SEXP, int);

#ifdef WIN32
static char *(*Vget_R_HOME) (void);
#endif

static void *dlget (void *handle, const char *name, int *err)
{
    void *p = gretl_dlsym(handle, name);
    
    if (p == NULL) {
	fprintf(stderr, "dlget: couldn't find '%s'\n", name);
	*err += 1;
    }

    return p;
}

static int load_R_symbols (void)
{
    const char *libpath = gretl_rlib_path();
    int err = 0;

    Rhandle = gretl_dlopen(libpath, 1);
    if (Rhandle == NULL) {
	return E_EXTERNAL;
    } 

    VCDR              = dlget(Rhandle, "CDR", &err);
    VREAL             = dlget(Rhandle, "REAL", &err);
    VRf_allocList     = dlget(Rhandle, "Rf_allocList", &err);
    VRf_allocMatrix   = dlget(Rhandle, "Rf_allocMatrix", &err);
    VRf_allocVector   = dlget(Rhandle, "Rf_allocVector", &err);
    VRf_endEmbeddedR  = dlget(Rhandle, "Rf_endEmbeddedR", &err);
    VRf_findFun       = dlget(Rhandle, "Rf_findFun", &err);
    VRf_findVar       = dlget(Rhandle, "Rf_findVar", &err);
    VRf_initEmbeddedR = dlget(Rhandle, "Rf_initEmbeddedR", &err);
    VRf_install       = dlget(Rhandle, "Rf_install", &err);
    VRf_isMatrix      = dlget(Rhandle, "Rf_isMatrix", &err);
    VRf_isLogical     = dlget(Rhandle, "Rf_isLogical", &err);
    VRf_isInteger     = dlget(Rhandle, "Rf_isInteger", &err);
    VRf_isReal        = dlget(Rhandle, "Rf_isReal", &err);
    VRf_lang2         = dlget(Rhandle, "Rf_lang2", &err);
    VRf_mkString      = dlget(Rhandle, "Rf_mkString", &err);
    VRf_ncols         = dlget(Rhandle, "Rf_ncols", &err);
    VRf_nrows         = dlget(Rhandle, "Rf_nrows", &err);
    VRf_PrintValue    = dlget(Rhandle, "Rf_PrintValue", &err);
    VRf_protect       = dlget(Rhandle, "Rf_protect", &err);
    VRf_ScalarReal    = dlget(Rhandle, "Rf_ScalarReal", &err);
    VRf_unprotect     = dlget(Rhandle, "Rf_unprotect", &err);
    VR_tryEval        = dlget(Rhandle, "R_tryEval", &err);
    VSETCAR           = dlget(Rhandle, "SETCAR", &err);
    VSET_TYPEOF       = dlget(Rhandle, "SET_TYPEOF", &err); 
    VTYPEOF           = dlget(Rhandle, "TYPEOF", &err);

#ifdef WIN32
    Vget_R_HOME = dlget(Rhandle, "get_R_HOME", &err);
#endif

    if (!err) {
	PVR_GlobalEnv    = (SEXP *) dlget(Rhandle, "R_GlobalEnv", &err);
	PVR_NilValue     = (SEXP *) dlget(Rhandle, "R_NilValue", &err);
	PVR_UnboundValue = (SEXP *) dlget(Rhandle, "R_UnboundValue", &err);
    }

    if (err) {
	close_plugin(Rhandle);
	Rhandle = NULL;
	err = E_EXTERNAL;
    }

    return err;
}

static int Rinit;

static SEXP current_arg;
static SEXP current_call;

void gretl_R_cleanup (void)
{
    if (Rinit) {
	VRf_endEmbeddedR(0);
	close_plugin(Rhandle);
    }
}

static int gretl_Rlib_init (void)
{
    char *Rhome;
    int err = 0;

#ifndef WIN32
    Rhome = getenv("R_HOME");
    if (Rhome == NULL) {
	fprintf(stderr, "To use Rlib, the variable R_HOME must be set\n");
	return E_EXTERNAL;
    }
#endif

    err = load_R_symbols();
    if (err) {
	fprintf(stderr, "gretl_Rlib_init: failed to load R functions\n");
	return err;
    }

#ifdef WIN32
    Rhome = Vget_R_HOME();
    if (Rhome == NULL) {
	fprintf(stderr, "To use Rlib, the variable R_HOME must be set\n");
	return E_EXTERNAL;
    }
#endif

    if (!err) {
	char *argv[] = { "R", "--silent" };

	VRf_initEmbeddedR(2, argv);
	VR_GlobalEnv = *PVR_GlobalEnv;
	VR_NilValue = *PVR_NilValue;
	VR_UnboundValue = *PVR_UnboundValue;
	Rinit = 1;
    }

    return err;
}

/* run R's source() function on the gretl-written R command file */

static int lib_run_Rlib_sync (gretlopt opt, PRN *prn) 
{
    gchar *Rsrc;
    SEXP e;
    int err = 0;

    if (!Rinit) {
	if ((err = gretl_Rlib_init())) {
	    fprintf(stderr, "lib_run_Rlib_sync: failed on gretl_Rlib_init\n");
	    return err;
	}
    }
    
    Rsrc = g_strdup_printf("%sRsrc", gretl_dot_dir());

    VRf_protect(e = VRf_lang2(VRf_install("source"), VRf_mkString(Rsrc)));
    VR_tryEval(e, VR_GlobalEnv, NULL);
    VRf_unprotect(1);

    g_free(Rsrc);

    return err;
}

static SEXP find_R_function (const char *name)
{
    SEXP fun;
    SEXPTYPE t;

    fun = VRf_findVar(VRf_install(name), VR_GlobalEnv);
    t = VTYPEOF(fun);

    if (t == PROMSXP) {
	/* eval promise if need be */
	int err = 1;

	fun = VR_tryEval(fun, VR_GlobalEnv, &err);
	if (!err) {
	    t = VTYPEOF(fun);
	}
    }

    if (t != CLOSXP && t != BUILTINSXP &&
	t != BUILTINSXP && t != SPECIALSXP) {
	return VR_UnboundValue;
    }

    return fun;
}

static int Rlib_err;

/* called from gretl_paths.c on revising the Rlib path:
   allow for the possibility that the path was wrong but is
   now OK
*/

void gretl_R_reset_error (void)
{
    Rlib_err = 0;
}

static int gretl_use_Rlib (void)
{
    int ret = 0;

    if (!Rlib_err && !getenv("GRETL_NO_RLIB")) {
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
    fprintf(stderr, "R: using %s\n", (ret)? "library" : "executable");
#endif    

    return ret;
}

/* used in "genr", to see if @name denotes an R function,
   either built-in or possibly user-defined
*/

int get_R_function_by_name (const char *name) 
{
    int ret = 0;
    
    if (gretl_use_Rlib()) {
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
    current_arg = VCDR(current_arg);
    VSETCAR(current_arg, VRf_ScalarReal(x));
    return 0;
}

int gretl_R_function_add_vector (const double *x, int t1, int t2) 
{
    SEXP res = VRf_allocVector(REALSXP, t2 - t1 + 1);
    int i;

    if (res == NULL) {
	return E_ALLOC;
    }

    current_arg = VCDR(current_arg);

    for (i=t1; i<=t2; i++) {
    	VREAL(res)[i-t1] = x[i];
    }
    
    VSETCAR(current_arg, res);
    return 0;
}

int gretl_R_function_add_matrix (const gretl_matrix *m) 
{
    int nr = gretl_matrix_rows(m);
    int nc = gretl_matrix_cols(m);
    SEXP res;
    int i, j;

    current_arg = VCDR(current_arg);
    res = VRf_allocMatrix(REALSXP, nr, nc);
    if (res == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<nr; i++) {
	for (j=0; j<nc; j++) {
	    VREAL(res)[i + j * nr] = gretl_matrix_get(m, i, j);
    	}
    }
    
    VSETCAR(current_arg, res);

    return 0;
}

/* called from geneval.c only, and should be pre-checked */

int gretl_R_get_call (const char *name, int argc) 
{
    SEXP call, e;

    call = VRf_findFun(VRf_install(name), VR_GlobalEnv);

    if (call == VR_NilValue) {
	fprintf(stderr, "gretl_R_get_call: no definition for function %s\n", 
		name);
	VRf_unprotect(1); /* is this OK? */
	return E_EXTERNAL;
    } 

    VRf_protect(e = VRf_allocList(argc + 1));
    VSET_TYPEOF(e, LANGSXP);
    VSETCAR(e, VRf_install(name));
    current_call = current_arg = e;
 
    return 0;
}

static int R_type_to_gretl_type (SEXP s)
{
    if (VRf_isMatrix(s)) {
	return GRETL_TYPE_MATRIX;
    } else if (VRf_isLogical(s)) {
	return GRETL_TYPE_BOOL;
    } else if (VRf_isInteger(s)) {
	return GRETL_TYPE_INT;
    } else if (VRf_isReal(s)) {
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
	VRf_PrintValue(current_call);
    }

    res = VR_tryEval(current_call, VR_GlobalEnv, &err);
    if (err) {
	return E_EXTERNAL;
    }

    *rtype = R_type_to_gretl_type(res);

    if (*rtype == GRETL_TYPE_MATRIX) {
	gretl_matrix *pm;
	int nc = VRf_ncols(res);
	int nr = VRf_nrows(res);
	int i, j;

	pm = gretl_matrix_alloc(nr, nc);
	if (pm == NULL) {
	    return E_ALLOC;
	}

	for (i=0; i<nr; i++) {
	    for (j=0; j<nc; j++) {
		gretl_matrix_set(pm, i, j, VREAL(res)[i + j * nr]);
	    }
	}
	VRf_unprotect(1);
	*ret = pm;
    } else if (gretl_scalar_type(*rtype)) {
	double *realres = VREAL(res);
	double *dret = *ret;

	*dret = *realres;

    	VRf_unprotect(1);
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

/* starting a "foreign" block from scratch, or adding a line
   to an existing block */

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

/* respond to "end foreign" [opts] */

int foreign_execute (const double **Z, const DATAINFO *pdinfo, 
		     gretlopt opt, PRN *prn)
{
    int i, err = 0;

    if (opt & OPT_V) {
	for (i=0; i<foreign_n_lines; i++) {
	    pprintf(prn, "> %s\n", foreign_lines[i]);
	}
    }

    foreign_opt |= opt;

    if (foreign_lang == LANG_R) {
	err = write_gretl_R_files(NULL, Z, pdinfo, foreign_opt);
	if (err) {
	    delete_gretl_R_files();
	} else {
#ifdef USE_RLIB
	    if (gretl_use_Rlib()) {
		lib_run_Rlib_sync(foreign_opt, prn);
	    } else {
		lib_run_R_sync(foreign_opt, prn);
	    }
#else
	    lib_run_R_sync(foreign_opt, prn);
#endif
	}
    } else if (foreign_lang == LANG_OX) {
	err = write_gretl_ox_file(NULL, foreign_opt);
	if (err) {
	    delete_gretl_ox_file();
	} else {
	    lib_run_ox_sync(foreign_opt, prn);
	}
    } 
    
    destroy_foreign();

    return err;
}

