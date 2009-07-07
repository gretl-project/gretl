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
    LANG_OX,
    LANG_RLIB
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
    } else if (!strcmp(lang, "RLib")) {
	foreign_lang = LANG_RLIB;
    } else if (!strcmp(lang, "ox")) {
	foreign_lang = LANG_OX;
    } else {
	pprintf(prn, "%s: unknown language\n");
	err = E_DATA;
    }

    return err;
}

#ifdef G_OS_WIN32

static int lib_run_R_sync (gretlopt opt, PRN *prn)
{
    gchar *Rterm, *cmd;
    int err = 0;

    Rterm = R_path_from_registry();
    if (Rterm == NULL) {
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
    g_free(Rterm);

    return err;
}

#else

static int lib_run_R_sync (gretlopt opt, PRN *prn)
{
    gchar *argv[6] = {
	"R", 
	"--no-save",
	"--no-init-file",
	"--no-restore-data",
	"--slave",
	NULL
    };
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
		pprintf(prn, "R exited with status %d", status);
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

#endif

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

#ifdef USE_RLIB

static int Rinit;

static SEXP current_arg;
static SEXP current_call;

void gretl_R_cleanup (void)
{
    if (Rinit) {
	Rf_endEmbeddedR(0);
    }
}

static int gretl_R_init (void)
{
    int err = 0;

    if (getenv("R_HOME") == NULL) {
	fprintf(stderr, "To use Rlib, the variable R_HOME must be set\n");
	err = E_EXTERNAL;
    } else {
	char *argv[] = { "R", "--silent" };

	Rf_initEmbeddedR(2, argv);
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
	if ((err = gretl_R_init())) {
	    return err;
	}
    }
    
    Rsrc = g_strdup_printf("%sRsrc", gretl_dot_dir());

    PROTECT(e = lang2(install("source"), mkString(Rsrc)));
    R_tryEval(e, R_GlobalEnv, NULL);
    UNPROTECT(1);

    g_free(Rsrc);

    return err;
}

static SEXP find_R_function (const char *name)
{
    SEXP fun;
    SEXPTYPE t;

    fun = findVar(install(name), R_GlobalEnv);
    t = TYPEOF(fun);

    if (t == PROMSXP) {
	/* eval promise if need be */
	int err = 1;

	fun = R_tryEval(fun, R_GlobalEnv, &err);
	if (!err) {
	    t = TYPEOF(fun);
	}
    }

    if (t != CLOSXP && t != BUILTINSXP &&
	t != BUILTINSXP && t != SPECIALSXP) {
	return R_UnboundValue;
    }

    return fun;
}

/* used in "genr" */

int get_R_function_by_name (const char *name) 
{
    SEXP fun;

    if (!Rinit) {
	int err = gretl_R_init();

	if (err) {
	    return err;
	}
    } 

    fun = find_R_function(name);

    return (fun == R_UnboundValue)? 0 : 1;
}

int gretl_R_function_add_scalar (double x) 
{
    current_arg = CDR(current_arg);
    SETCAR(current_arg, ScalarReal(x));
    return 0;
}

int gretl_R_function_add_vector (const double *x, int t1, int t2) 
{
    SEXP res = allocVector(REALSXP, t2 - t1 + 1);
    int i;

    if (res == NULL) {
	return E_ALLOC;
    }

    current_arg = CDR(current_arg);

    for (i=t1; i<=t2; i++) {
    	REAL(res)[i-t1] = x[i];
    }
    
    SETCAR(current_arg, res);
    return 0;
}

int gretl_R_function_add_matrix (const gretl_matrix *m) 
{
    int nr = gretl_matrix_rows(m);
    int nc = gretl_matrix_cols(m);
    SEXP res;
    int i, j;

    current_arg = CDR(current_arg);
    res = allocMatrix(REALSXP, nr, nc);
    if (res == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<nr; i++) {
	for (j=0; j<nc; j++) {
	    REAL(res)[i + j * nr] = gretl_matrix_get(m, i, j);
    	}
    }
    
    SETCAR(current_arg, res);
    return 0;
}

int gretl_R_get_call (const char *name, int argc) 
{
    SEXP call, e;

    call = findFun(install(name), R_GlobalEnv);

    if (call == R_NilValue) {
	fprintf(stderr, "gretl_R_get_call: no definition for function %s\n", 
		name);
	UNPROTECT(1); /* is this OK? */
	return E_EXTERNAL;
    } 

    PROTECT(e = allocList(argc + 1));
    SET_TYPEOF(e, LANGSXP);
    SETCAR(e, install(name));
    current_call = current_arg = e;
 
    return 0;
}

int gretl_R_function_exec (const char *name, int *rtype, void **ret) 
{
    SEXP res;
    int err = 0;

    /* FIXME make this optional? */
    PrintValue(current_call);

    res = R_tryEval(current_call, R_GlobalEnv, NULL);

    if (isMatrix(res)) {
	gretl_matrix *pm;
	int nc = ncols(res);
	int nr = nrows(res);
	int i, j;

	pm = gretl_matrix_alloc(nr, nc);
	if (pm == NULL) {
	    return E_ALLOC;
	}

	for (i=0; i<nr; i++) {
	    for (j=0; j<nc; j++) {
		gretl_matrix_set(pm, i, j,REAL(res)[i + j * nr]);
	    }
	}
	UNPROTECT(1);
	*ret = pm;
	*rtype = GRETL_TYPE_MATRIX;
    } else if (isVectorAtomic(res)) {
	double *realres = REAL(res);
	double *dret = *ret;

	*dret = *realres;

    	UNPROTECT(1);
	*rtype = GRETL_TYPE_DOUBLE;
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

    if (foreign_lang == LANG_RLIB) {
#ifdef USE_RLIB
	err = write_gretl_R_files(NULL, Z, pdinfo, foreign_opt);
	if (err) {
	    delete_gretl_R_files();
	} else {
	    lib_run_Rlib_sync(foreign_opt, prn);
	}
#else
	pputs(prn, "language=Rlib: not supported\n");
	destroy_foreign();
	return E_EXTERNAL;
#endif
    } else if (foreign_lang == LANG_R) {
	err = write_gretl_R_files(NULL, Z, pdinfo, foreign_opt);
	if (err) {
	    delete_gretl_R_files();
	} else {
	    lib_run_R_sync(foreign_opt, prn);
	}
    } else if (foreign_lang == LANG_OX) {
	pputs(prn, "language=ox: not supported yet\n");
	err = E_NOTIMP;
    } 
    
    destroy_foreign();

    return err;
}

