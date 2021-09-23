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

/*  Apparatus for interfacing with the lpsolve library for
    solution of linear programming problems.
*/

#include "libgretl.h"
#include "version.h"

#if defined(PKGBUILD)
# define PRELINKED
#elif defined(WIN32)
# include <windows.h>
#else
# include <dlfcn.h>
#endif

/* The bits of the lpsolve API that we need; we prefer not to include
   lp_lib.h since it contains some defines that conflict with GLib.
*/

typedef struct _lprec lprec;
typedef double REAL;

enum { LE = 1, GE, EQ };

/* reporting levels */
#define NEUTRAL    0
#define CRITICAL   1
#define SEVERE     2
#define IMPORTANT  3
#define NORMAL     4
#define DETAILED   5
#define FULL       6

/* solver status values */
#define UNKNOWNERROR  -5
#define DATAIGNORED   -4
#define NOBFP         -3
#define NOMEMORY      -2
#define NOTRUN        -1
#define OPTIMAL        0
#define SUBOPTIMAL     1
#define INFEASIBLE     2
#define UNBOUNDED      3
#define DEGENERATE     4
#define NUMFAILURE     5
#define USERABORT      6
#define TIMEOUT        7
#define RUNNING        8
#define PRESOLVED      9

#ifdef PRELINKED

/* In our macOS and Windows packages we ease the user's path
   by including the lpsolve shared library and linking this
   plugin against it. Otherwise we look to load the required
   symbols on demand at run-time. On Linux this should be
   easy.
*/

lprec *make_lp (int rows, int columns);
lprec *read_lp (FILE *fp, int verbose, char *lp_name);
unsigned char set_add_rowmode (lprec *lp, unsigned char s);
void delete_lp (lprec *lp);
void set_verbose (lprec *lp, int verbose);
void set_maxim (lprec *lp);
void set_minim (lprec *lp);
unsigned char set_lp_name (lprec *lp, char *s);
unsigned char set_obj_fn (lprec *lp, REAL *row);
unsigned char add_constraint (lprec *lp, REAL *row,
			      int constr_type, REAL rh);
unsigned char set_col_name (lprec *lp, int col, char *name);
unsigned char set_row_name (lprec *lp, int row, char *name);
char *get_col_name (lprec *lp, int col);
char *get_row_name (lprec *lp, int row);
unsigned char set_int (lprec *lp, int col, unsigned char s);
int get_Nrows (lprec *lp);
int get_Ncolumns (lprec *lp);
int solve (lprec *lp);
void print_objective (lprec *lp);
void print_solution (lprec *lp, int columns);
void print_constraints (lprec *lp, int columns);
void print_duals (lprec *lp);
REAL get_objective (lprec *lp);
REAL get_accuracy (lprec *lp);
unsigned char get_primal_solution (lprec *lp, REAL *pv);
unsigned char get_dual_solution (lprec *lp, REAL *duals);
unsigned char get_sensitivity_rhs (lprec *lp, REAL *duals,
				   REAL *from, REAL *till);
void set_outputstream (lprec *lp, FILE *fp);

#else

static lprec *(*make_lp) (int rows, int columns);
static lprec *(*read_lp) (FILE *fp, int verbose, char *lp_name);
static unsigned char (*set_add_rowmode) (lprec *lp, unsigned char s);
static void (*delete_lp) (lprec *lp);
static void (*set_verbose) (lprec *lp, int verbose);
static void (*set_maxim) (lprec *lp);
static void (*set_minim) (lprec *lp);
static unsigned char (*set_lp_name) (lprec *lp, char *s);
static unsigned char (*set_obj_fn) (lprec *lp, REAL *row);
static unsigned char (*add_constraint) (lprec *lp, REAL *row,
					int constr_type, REAL rh);
static unsigned char (*set_col_name) (lprec *lp, int col, char *name);
static unsigned char (*set_row_name) (lprec *lp, int row, char *name);
static char *(*get_col_name) (lprec *lp, int col);
static char *(*get_row_name) (lprec *lp, int row);
static unsigned char (*set_int) (lprec *lp, int col, unsigned char s);
static int (*get_Nrows) (lprec *lp);
static int (*get_Ncolumns) (lprec *lp);
static int (*solve) (lprec *lp);
static void (*print_objective) (lprec *lp);
static void (*print_solution) (lprec *lp, int columns);
static void (*print_constraints) (lprec *lp, int columns);
static void (*print_duals) (lprec *lp);
static REAL (*get_objective) (lprec *lp);
static REAL (*get_accuracy) (lprec *lp);
static unsigned char (*get_primal_solution) (lprec *lp, REAL *pv);
static unsigned char (*get_dual_solution) (lprec *lp, REAL *duals);
static unsigned char (*get_sensitivity_rhs) (lprec *lp, REAL *duals,
					     REAL *from, REAL *till);
static void (*set_outputstream) (lprec *lp, FILE *fp);

static void *lphandle;            /* handle to the lpsolve library */
static int gretl_lpsolve_err;     /* initialization error record */
static int gretl_lpsolve_initted; /* are we initialized or not? */

static void *lpget (void *handle, const char *name, int *err)
{
#ifdef WIN32
    void *p = GetProcAddress(handle, name);
#else
    void *p = dlsym(handle, name);
#endif

    if (p == NULL && err != NULL) {
	printf("lpget: couldn't find '%s'\n", name);
	*err += 1;
    }

    return p;
}

static int gretl_lpsolve_init (void)
{
    int err = 0;

    if (gretl_lpsolve_initted) {
	return gretl_lpsolve_err;
    }

#ifdef WIN32
    lphandle = LoadLibrary(gretl_lpsolve_path());
#else
    lphandle = dlopen(gretl_lpsolve_path(), RTLD_NOW);
#endif

    if (lphandle == NULL) {
	err = E_EXTERNAL;
    } else {
	make_lp             = lpget(lphandle, "make_lp", &err);
	read_lp             = lpget(lphandle, "read_lp", &err);
	set_add_rowmode     = lpget(lphandle, "set_add_rowmode", &err);
	set_lp_name         = lpget(lphandle, "set_lp_name", &err);
	delete_lp           = lpget(lphandle, "delete_lp", &err);
	set_verbose         = lpget(lphandle, "set_verbose", &err);
	set_maxim           = lpget(lphandle, "set_maxim", &err);
	set_minim           = lpget(lphandle, "set_minim", &err);
	set_obj_fn          = lpget(lphandle, "set_obj_fn", &err);
	add_constraint      = lpget(lphandle, "add_constraint", &err);
	set_col_name        = lpget(lphandle, "set_col_name", &err);
	set_row_name        = lpget(lphandle, "set_row_name", &err);
	get_col_name        = lpget(lphandle, "get_col_name", &err);
	get_row_name        = lpget(lphandle, "get_row_name", &err);
	set_int             = lpget(lphandle, "set_int", &err);
	get_Nrows           = lpget(lphandle, "get_Nrows", &err);
	get_Ncolumns        = lpget(lphandle, "get_Ncolumns", &err);
	solve               = lpget(lphandle, "solve", &err);
	print_objective     = lpget(lphandle, "print_objective", &err);
	print_solution      = lpget(lphandle, "print_solution", &err);
	print_constraints   = lpget(lphandle, "print_constraints", &err);
	print_duals         = lpget(lphandle, "print_duals", &err);
	get_objective       = lpget(lphandle, "get_objective", &err);
	get_primal_solution = lpget(lphandle, "get_primal_solution", &err);
	get_dual_solution   = lpget(lphandle, "get_dual_solution", &err);
	get_sensitivity_rhs = lpget(lphandle, "get_sensitivity_rhs", &err);
	set_outputstream    = lpget(lphandle, "set_outputstream", &err);
	get_accuracy        = lpget(lphandle, "get_accuracy", NULL);

	if (err) {
	    close_plugin(lphandle);
	    lphandle = NULL;
	    err = E_EXTERNAL;
	}
    }

    gretl_lpsolve_err = err;
    gretl_lpsolve_initted = 1;

    return err;
}

#endif /* not PRELINKED to lpsolve library */

static void lp_row_from_mrow (double *targ, int nv,
			      const gretl_matrix *m,
			      int i)
{
    int j;

    for (j=1; j<=nv; j++) {
	targ[j] = gretl_matrix_get(m, i, j-1);
    }
}

static int lp_ctype_from_string (const char *s)
{
    if (!strcmp(s, "=")) return EQ;
    else if (!strcmp(s, "<=")) return LE;
    else if (!strcmp(s, ">=")) return GE;
    else return -1;
}

/* We'll accept constraint-type specifiers either as strings
   or as integer values */

static int *get_constraint_types (gretl_bundle *b,
				  int nc, int *err)
{
    gretl_array *S = NULL;
    gretl_matrix *m = NULL;
    int *ret = NULL;
    GretlType t;
    void *ptr;

    ptr = gretl_bundle_get_data(b, "ctypes", &t, NULL, err);

    if (!*err) {
	if (t == GRETL_TYPE_ARRAY) {
	    S = ptr;
	    if (gretl_array_get_type(S) != GRETL_TYPE_STRINGS ||
		gretl_array_get_length(S) != nc) {
		*err = E_DATA;
	    }
	} else if (t == GRETL_TYPE_MATRIX) {
	    m = ptr;
	    if (gretl_vector_get_length(m) != nc) {
		*err = E_DATA;
	    }
	} else {
	    *err = E_DATA;
	}
    }

    if (!*err) {
	int i;

	ret = malloc(nc * sizeof *ret);
	if (ret == NULL) {
	    *err = E_ALLOC;
	}
	for (i=0; i<nc && !*err; i++) {
	    if (S != NULL) {
		ret[i] = lp_ctype_from_string(gretl_array_get_data(S, i));
	    } else {
		ret[i] = m->val[i];
	    }
	    if (ret[i] < LE || ret[i] > EQ) {
		*err = E_DATA;
	    }
	}
	if (*err) {
	    free(ret);
	    ret = NULL;
	}
    }

    return ret;
}

static lprec *lp_model_from_bundle (gretl_bundle *b,
				    const char ***pcnames,
				    const char ***prnames,
				    gretlopt opt,
				    int *err)
{
    lprec *lp = NULL;
    const gretl_matrix *O = NULL; /* objective coeffs */
    const gretl_matrix *R = NULL; /* constraints */
    const gretl_matrix *K = NULL; /* integer spec */
    int *ctypes = NULL;           /* constraint types */
    int nv = 0;  /* number of variables */
    int nc = 0;  /* number of constraints */
    int ni = 0;  /* number of integer variables */

    /* required model data */
    O = gretl_bundle_get_matrix(b, "objective", err);
    if (!*err) {
	R = gretl_bundle_get_matrix(b, "constraints", err);
    }

    if (!*err) {
	/* checks on @O and @R */
	nv = gretl_vector_get_length(O);
	if (nv == 0) {
	    *err = E_DATA;
	} else if (R->cols != nv + 1) {
	    *err = E_DATA;
	} else if (R->rows < 1) {
	    *err = E_DATA;
	} else {
	    nc = R->rows;
	}
	if (!*err) {
	    ctypes = get_constraint_types(b, nc, err);
	}
    }

    if (!*err) {
	/* check optional setting of integer variables */
	K = gretl_bundle_get_matrix(b, "intvars", NULL);
	if (K != NULL) {
	    ni = gretl_vector_get_length(K);
	    if (ni == 0 || ni > nv) {
		*err = E_DATA;
	    }
	}
    }

    if (!*err) {
	const char **cnames = gretl_matrix_get_colnames(O);
	const char **rnames = gretl_matrix_get_rownames(R);
	double rhs, *row = calloc(nv+1, sizeof *row);
	int col, i, j;

	lp = make_lp(0, nv);

	if (lp == NULL) {
	    *err = E_ALLOC;
	} else {
	    /* set lpsolve verbosity level */
	    if (opt & OPT_V) {
		set_verbose(lp, NORMAL);
	    } else {
		set_verbose(lp, CRITICAL);
	    }
	    /* set rowmode: this has to be placed carefully! */
	    set_add_rowmode(lp, 1);
	    for (j=0; j<nv; j++) {
		row[j+1] = O->val[j];
		if (cnames != NULL) {
		    set_col_name(lp, j+1, (char *) cnames[j]);
		}
	    }
	    set_obj_fn(lp, row);
	    /* integer variable specifiers? */
	    for (j=0; j<ni; j++) {
		col = K->val[j];
		if (col > 0 && col <= nv) {
		    set_int(lp, col, 1);
		}
	    }
	    /* constraints */
	    for (i=0; i<nc; i++) {
		lp_row_from_mrow(row, nv, R, i);
		rhs = gretl_matrix_get(R, i, nv);
		add_constraint(lp, row, ctypes[i], rhs);
		if (rnames != NULL) {
		    set_row_name(lp, i+1, (char *) rnames[i]);
		}
	    }
	    /* turn off rowmode or lpsolve will likely crash */
	    set_add_rowmode(lp, 0);
	    *pcnames = cnames;
	    *prnames = rnames;
	}
	free(row);
    }

    free(ctypes);

    return lp;
}

static lprec *lp_model_from_file (const char *fname,
				  const char *buf,
				  PRN *prn, int *err)
{
    lprec *lp = NULL;
    gchar *tmp = NULL;
    FILE *fp = NULL;

    if (fname != NULL) {
	fp = gretl_fopen(fname, "r");
    } else {
	tmp = gretl_make_dotpath("tmp.lp");
	fp = gretl_fopen(tmp, "wb");
	if (fp != NULL) {
	    fputs(buf, fp);
	    fclose(fp);
	    fp = gretl_fopen(tmp, "rb");
	}
    }

    if (fp == NULL) {
	*err = E_FOPEN;
    } else {
	if (fname != NULL && prn != NULL) {
	    pprintf(prn, "Reading input from '%s'\n", fname);
	}
	lp = read_lp(fp, NORMAL, NULL);
	fclose(fp); /* check that lpsolve doesn't close close this? */
    }

    if (tmp != NULL) {
	gretl_remove(tmp);
	g_free(tmp);
    }

    return lp;
}

static int get_name_len (const char **S, int ns)
{
    int i, len, maxlen = 0;

    for (i=0; i<ns; i++) {
	len = strlen(S[i]);
	if (len > maxlen) {
	    maxlen = len;
	}
    }

    return maxlen;
}

static void print_lpsolve_output (lprec *lp,
				  gretl_matrix *VV,
				  gretl_matrix *VC,
				  PRN *prn)
{
    gchar *fname = gretl_make_dotpath("lptmp.txt");
    FILE *fp = gretl_fopen(fname, "wb");
    const char **SVV = gretl_matrix_get_rownames(VV);
    const char **SVC = gretl_matrix_get_rownames(VC);
    int vlen = 0, clen = 0;

    if (SVV != NULL) {
	vlen = get_name_len(SVV, VV->rows);
    }
    if (SVC != NULL) {
	clen = get_name_len(SVC, VC->rows);
    }

    if (fp != NULL) {
	gchar *buf = NULL;
	int i;

	set_outputstream(lp, fp);
	print_objective(lp);
	fputs("\nValues of the variables:\n", fp);
	for (i=0; i<VV->rows; i++) {
	    if (vlen > 0) {
		fprintf(fp, "%*s  %#g\n", vlen, SVV[i], VV->val[i]);
	    } else {
		fprintf(fp, "C%d  %#g\n", i+1, VV->val[i]);
	    }
	}
	fputs("\nValues of the constraints:\n", fp);
	for (i=0; i<VC->rows; i++) {
	    if (clen > 0) {
		fprintf(fp, "%*s  %g\n", clen, SVC[i], VC->val[i]);
	    } else {
		fprintf(fp, "R%d  %g\n", i+1, VC->val[i]);
	    }
	}
	print_duals(lp);
	fputc('\n', fp);
	fflush(fp);
	set_outputstream(lp, stdout);
	if (g_file_get_contents(fname, &buf, NULL, NULL)) {
	    pputs(prn, buf);
	    g_free(buf);
	}
	gretl_remove(fname);
    }

    g_free(fname);
}

static char **get_lp_colnames (lprec *lp, int n)
{
    char **S = strings_array_new(n);
    int i;

    for (i=0; i<n; i++) {
	S[i] = gretl_strdup(get_col_name(lp, i+1));
    }
    return S;
}

static char **get_lp_rownames (lprec *lp, int n)
{
    char **S = strings_array_new(n);
    int i;

    for (i=0; i<n; i++) {
	S[i] = gretl_strdup(get_row_name(lp, i+1));
    }
    return S;
}

static int get_lp_model_data (lprec *lp, gretl_bundle *ret,
			      const char **cnames,
			      const char **rnames,
			      gretlopt opt, PRN *prn)
{
    int nr = get_Nrows(lp);
    int nc = get_Ncolumns(lp);
    int i, k, psize = 1 + nr + nc;
    double *prim = malloc(psize * sizeof *prim);
    double *dual = malloc(psize * sizeof *dual);
    gretl_matrix *VC = gretl_matrix_alloc(nr, 1);
    gretl_matrix *VV = gretl_matrix_alloc(nc, 1);
    gretl_matrix *SP = gretl_matrix_alloc(nr, 1);
    char **S = NULL;

    if (prim == NULL || dual == NULL ||
	VC == NULL || VV == NULL || SP == NULL) {
	return E_ALLOC;
    }

    if (cnames != NULL) {
	S = strings_array_dup((char **) cnames, VV->rows);
	gretl_matrix_set_rownames(VV, S);
    } else {
	S = get_lp_colnames(lp, VV->rows);
	gretl_matrix_set_rownames(VV, S);
    }

    if (rnames != NULL) {
	S = strings_array_dup((char **) rnames, VC->rows);
	gretl_matrix_set_rownames(VC, S);
	S = strings_array_dup((char **) rnames, VC->rows);
	gretl_matrix_set_rownames(SP, S);
    } else {
	S = get_lp_rownames(lp, VC->rows);
	gretl_matrix_set_rownames(VC, S);
	S = get_lp_rownames(lp, VC->rows);
	gretl_matrix_set_rownames(SP, S);
    }

    gretl_bundle_set_scalar(ret, "objective", get_objective(lp));

    /* note: get_accuracy() requires lpsolve >= 5.5.2.11 */
#ifdef PRELINKED
    gretl_bundle_set_scalar(ret, "accuracy", get_accuracy(lp));
#else
    if (get_accuracy != NULL) {
	gretl_bundle_set_scalar(ret, "accuracy", get_accuracy(lp));
    }
#endif

    if (get_primal_solution(lp, prim)) {
	k = 1;
	for (i=0; i<nr; i++) {
	    VC->val[i] = prim[k++];
	}
	for (i=0; i<nc; i++) {
	    VV->val[i] = prim[k++];
	}
    }

    if (get_dual_solution(lp, dual)) {
	for (i=0; i<nr; i++) {
	    SP->val[i] = dual[i+1];
	}
    }

    if (opt & OPT_S) {
	int ns, n = psize - 1;
	gretl_matrix *SE = gretl_matrix_alloc(n, 3);
	double *x = SE->val;

	if (!get_sensitivity_rhs(lp, x, x + n, x + 2*n)) {
	    gretl_matrix_free(SE);
	} else {
	    S = gretl_string_split("dual lower upper", &ns, NULL);
	    gretl_matrix_set_colnames(SE, S);
	    gretl_bundle_donate_data(ret, "sensitivity", SE, GRETL_TYPE_MATRIX, 0);

	}
    }

    if (prn != NULL) {
	print_lpsolve_output(lp, VV, VC, prn);
    }

    gretl_bundle_donate_data(ret, "constraint_values", VC, GRETL_TYPE_MATRIX, 0);
    gretl_bundle_donate_data(ret, "variable_values", VV, GRETL_TYPE_MATRIX, 0);
    gretl_bundle_donate_data(ret, "shadow_prices", SP, GRETL_TYPE_MATRIX, 0);

    free(prim);
    free(dual);

    return 0;
}

static gretlopt lp_options_from_bundle (gretl_bundle *b,
					const char **pfname,
					const char **pbuf,
					gretl_bundle **pb,
					const char **pname)
{
    gretlopt opt = OPT_NONE;

    if (gretl_bundle_get_bool(b, "minimize", 0)) {
	opt |= OPT_I;
    }
    if (gretl_bundle_get_bool(b, "verbose", 0)) {
	opt |= OPT_V;
    }
    if (gretl_bundle_get_bool(b, "sensitivity", 0)) {
	opt |= OPT_S;
    }

    if (gretl_bundle_has_key(b, "model_name")) {
	*pname = gretl_bundle_get_string(b, "model_name", NULL);
    }

    if (gretl_bundle_has_key(b, "lp_filename")) {
	*pfname = gretl_bundle_get_string(b, "lp_filename", NULL);
    } else if (gretl_bundle_has_key(b, "lp_buffer")) {
	*pbuf = gretl_bundle_get_string(b, "lp_buffer", NULL);
    } else if (gretl_bundle_has_key(b, "lp_bundle")) {
	*pb = gretl_bundle_get_bundle(b, "lp_bundle", NULL);
    }

    return opt;
}

/* Run the lpsolve solve() function, with provision to
   catch output that would by default go to stdout in
   case of error, or when verbose output is requested.
*/

static int maybe_catch_solve (lprec *lp, gretlopt opt,
			      PRN *prn)
{
    gchar *fname = gretl_make_dotpath("solve.txt");
    FILE *fp = gretl_fopen(fname, "wb");
    int retval;

    if (fp != NULL) {
	gchar *buf = NULL;

	set_outputstream(lp, fp);
	retval = solve(lp);
	set_outputstream(lp, stdout);
	if ((opt & OPT_V) || retval != OPTIMAL) {
	    if (g_file_get_contents(fname, &buf, NULL, NULL)) {
		pputs(prn, buf);
		g_free(buf);
	    }
	}
	gretl_remove(fname);
    } else {
	retval = solve(lp);
    }

    g_free(fname);

    return retval;
}

/* driver function */

gretl_bundle *gretl_lpsolve (gretl_bundle *b, PRN *prn, int *err)
{
    gretl_bundle *ret = NULL;
    gretl_bundle *lpb = NULL;
    const char **cnames = NULL;
    const char **rnames = NULL;
    const char *lpfname = NULL;
    const char *lpbuf = NULL;
    const char *lpname = NULL;
    lprec *lp = NULL;
    PRN *vprn = NULL;
    int msg_set = 0;
    gretlopt opt;

#ifndef PRELINKED
    if (!gretl_lpsolve_initted) {
	gretl_lpsolve_init();
    }

    if (gretl_lpsolve_err) {
	gretl_errmsg_sprintf("lpsolve: failed to load %s", gretl_lpsolve_path());
	*err = gretl_lpsolve_err;
	return NULL;
    }
#endif

    opt = lp_options_from_bundle(b, &lpfname, &lpbuf, &lpb, &lpname);
    vprn = (opt & OPT_V)? prn : NULL;

    if (lpfname != NULL || lpbuf != NULL) {
	lp = lp_model_from_file(lpfname, lpbuf, vprn, err);
    } else if (lpb != NULL) {
	lp = lp_model_from_bundle(lpb, &cnames, &rnames, opt, err);
    } else {
	gretl_errmsg_set("lpsolve: didn't find a model specification");
	msg_set = 1;
	*err = E_ARGS;
    }

    if (*err) {
	if (!msg_set) {
	    gretl_errmsg_set("lpsolve: failed to build model");
	}
    } else {
	if (lpname != NULL) {
	    set_lp_name(lp, (char *) lpname);
	}
	if (opt & OPT_I) {
	    set_minim(lp);
	} else {
	    set_maxim(lp);
	}
	*err = maybe_catch_solve(lp, opt, prn);
	if (*err) {
	    gretl_errmsg_set("lpsolve: solution failed");
	} else {
	    ret = gretl_bundle_new();
	    *err = get_lp_model_data(lp, ret, cnames, rnames,
				     opt, vprn);
	}
    }

    if (lp != NULL) {
	delete_lp(lp);
    }

    return ret;
}
