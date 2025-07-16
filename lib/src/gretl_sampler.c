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
#include "gretl_func.h"
#include "gretl_typemap.h"
#include "uservar.h"
#include "libset.h"
#include "gretl_sampler.h"

/* first: apparatus pertaining to gibbs_via_bundles() */

#define BDEBUG 0

static gretl_matrix *run_extractor_func (const char *fname,
                                         gretl_bundle *b,
                                         PRN *prn)
{
    gretl_matrix *ret = NULL;
    fncall *fcall = NULL;
    ufunc *func;
    int err = 0;

    func = get_user_function_by_name(fname);

    if (func != NULL) {
        fcall = fncall_new(func, 0);
        if (fcall != NULL) {
            err = push_anon_function_arg(fcall, GRETL_TYPE_BUNDLE, b);
        }
        if (!err) {
            err = gretl_function_exec(fcall, GRETL_TYPE_MATRIX, NULL,
                                      &ret, prn);
            if (err) {
                fprintf(stderr, "run_extractor_func: err %d\n", err);
            }
        }
    }

    return ret;
}

static gretl_matrix *extract_param (const char *s, int i,
                                    gretl_bundle *allInfo,
                                    PRN *prn, int *err)
{
    gretl_matrix *ret = NULL;
    gretl_array *dists;
    gchar *tmp = g_strdup(s);
    int n, k;

    dists = gretl_bundle_get_array(allInfo, "dists", err);

    g_strstrip(tmp);
    n = strlen(tmp);

    if (*tmp == '[') {
        if (tmp[n-1] == ']') {
            /* matrix from private array */
            gretl_bundle *b = gretl_array_get_data(dists, i);
            gretl_array *a = gretl_bundle_get_array(b, "private", err);

            sscanf(tmp+1, "%d]", &k);
            ret = gretl_array_get_data(a, k-1);
        } else {
            *err = E_INVARG;
        }
    } else if (*tmp == '(') {
        if (tmp[n-1] == ')') {
            /* @value from specified dists member */
            gretl_bundle *b;

            sscanf(tmp+1, "%d)", &k);
            b = gretl_array_get_data(dists, k-1);
            ret = gretl_bundle_get_matrix(b, "value", err);
        } else {
            *err = E_INVARG;
        }
    } else if (!strncmp(tmp, "__", 2)) {
        /* should be a user-written function */
        ret = run_extractor_func(tmp + 2, allInfo, prn);
    } else {
        *err = E_INVARG;
    }

#if BDEBUG
    fprintf(stderr, "extract_param: i=%d, s='%s'\n", i, s);
#endif

    g_free(tmp);

    return ret;
}

static gretl_matrix *make_N_matrix (gretl_matrix *m,
                                    gretl_matrix *S,
                                    int dim, int id)
{
    gretl_matrix *C = NULL;
    gretl_matrix *Mn = NULL;
    gretl_matrix *v = NULL;
    int err = 0;

    /* C = cholesky(S) */
    C = gretl_matrix_copy(S);
    err = gretl_matrix_cholesky_decomp(C);
    if (err) {
        fprintf(stderr, "make_N_matrix: cholesky(S) failed\n");
        return NULL;
    }

    if (id == D_NORMAL2) {
        /* v = S * m */
        v = gretl_matrix_alloc(S->rows, m->cols);
        gretl_matrix_multiply(S, m, v);
    } else {
        v = gretl_matrix_copy(m);
    }

    /* Mn = mnormal(dim,1) */
    Mn = gretl_random_matrix_new(dim, 1, D_NORMAL);

    /* v += C * Mn */
    gretl_matrix_multiply_mod(C,  GRETL_MOD_NONE,
                              Mn, GRETL_MOD_NONE,
                              v,  GRETL_MOD_CUMULATE);
    gretl_matrix_transpose_in_place(v);
#if BDEBUG
    gretl_matrix_print(v, "v' in make_N_matrix");
#endif

    gretl_matrix_free(C);
    gretl_matrix_free(Mn);

    return v;
}

#define mapstr(m,i) gretl_array_get_data(m,i)

static double gen_one (gretl_bundle *B, int which, PRN *prn)
{
    gretl_array *dists;
    gretl_bundle *this;
    gretl_array *map;
    const char *s;
    int i, id;
    int err = 0;

    dists = gretl_bundle_get_array(B, "dists", &err);
    this = gretl_array_get_data(dists, which);
    id = gretl_bundle_get_int(this, "id", &err);
    map = gretl_bundle_get_array(this, "map", &err);

    if (id == D_NORMAL || id == D_NORMAL2) {
        gretl_matrix *md, *m, *S;
        gretl_matrix **MM[3] = {&md, &m, &S};
        int mfree[3];
        gretl_matrix *V;
        int dim;

        for (i=0; i<3 && !err; i++) {
            s = gretl_array_get_data(map, i);
            mfree[i] = (*s == '_');
            *(MM[i]) = extract_param(s, which, B, prn, &err);
        }
        dim = (int) md->val[0];
        V = make_N_matrix(m, S, dim, id);
        for (i=0; i<3; i++) {
            if (mfree[i]) gretl_matrix_free(*(MM[i]));
        }
        gretl_bundle_donate_data(this, "value", V, GRETL_TYPE_MATRIX, 0);
    } else if (id == D_UNIFORM) {
        gretl_matrix *ma, *mb;
        gretl_matrix **MM[2] = {&ma, &mb};
        int mfree[2];
        double d, parm[2];

        for (i=0; i<2 && !err; i++) {
            s = gretl_array_get_data(map, i);
            mfree[i] = (*s == '_');
            *(MM[i]) = extract_param(s, which, B, prn, &err);
        }
        parm[0] = ma->val[0];
        parm[1] = mb->val[0];
        for (i=0; i<2; i++) {
            if (mfree[i]) gretl_matrix_free(*(MM[i]));
        }
        d = gretl_get_random_scalar(D_UNIFORM, parm, &err);
        gretl_bundle_set_scalar(this, "value", d);
    } else if (id == D_BINOMIAL) {
        gretl_matrix *mp;
        int mfree = 0;
        double parm[2] = {0, 1};
        double p, d;

        s = gretl_array_get_data(map, 0);
        mfree = (*s == '_');
        mp = extract_param(s, which, B, prn, &err);
        p = mp->val[0];
        if (mfree) gretl_matrix_free(mp);
        d = gretl_get_random_scalar(D_UNIFORM, parm, &err);
        gretl_bundle_set_scalar(this, "value", d < p);
    } else if (id == D_GAMMA || id == D_IGAMMA) {
        gretl_matrix *mp, *ma;
        gretl_matrix **MM[2] = {&mp, &ma};
        int mfree[2];
        double d, parm[2];

        for (i=0; i<2 && !err; i++) {
            s = gretl_array_get_data(map, i);
            mfree[i] = (*s == '_');
            *(MM[i]) = extract_param(s, which, B, prn, &err);
        }
        parm[0] = mp->val[0];
        parm[1] = 1.0 / ma->val[0];
        for (i=0; i<2; i++) {
            if (mfree[i]) gretl_matrix_free(*(MM[i]));
        }
        d = gretl_get_random_scalar(D_GAMMA, parm, &err);
        if (id == D_GAMMA) {
            gretl_bundle_set_scalar(this, "value", d);
        } else {
            gretl_bundle_set_scalar(this, "value", 1/d);
        }
    }

    return err;
}

/* initialises the array; returns a int array with the
   dimension of each element */

static int *gibbs_init (gretl_bundle *allInfo, PRN *prn, int *err)
{
    gretl_array *dists;
    gretl_bundle *this;
    const char *idstr;
    int id;
    int *ret;
    int n, i;

    dists = gretl_bundle_get_array(allInfo, "dists", err);
    if (*err) {
        return NULL;
    }

    n = gretl_array_get_length(dists);
    ret = calloc(n+1, sizeof *ret);

    for (i=0; i<n && !*err; i++) {
        ret[i] = 0;
        this = gretl_array_get_data(dists, i);
        idstr = gretl_bundle_get_string(this, "dist", err);
        id = dist_code_from_string(idstr);
        gretl_bundle_set_int(this, "id", id);
        if (id == D_NORMAL || id == D_NORMAL2) {
            gretl_array *map;
            gretl_matrix *md;
            gretl_matrix *Mn;
            int dim;

            map = gretl_bundle_get_array(this, "map", err);
            md = extract_param(mapstr(map, 0), i, allInfo, prn, err);
            dim = md->val[0];
            Mn = gretl_random_matrix_new(dim, 1, D_NORMAL);
            gretl_bundle_donate_data(this, "value", Mn, GRETL_TYPE_MATRIX, 0);
            ret[i] = dim;
        } else if (id == D_UNIFORM) {
            gretl_bundle_set_scalar(this, "value", 0.5);
            ret[i] = 1;
        } else if (id == D_BINOMIAL) {
            double parm[2] = {0, 1};
            double d;

            d = gretl_get_random_scalar(D_UNIFORM, parm, err);
            gretl_bundle_set_scalar(this, "value", d > 0.5);
            ret[i] = 1;
        } else if (id == D_GAMMA || id == D_IGAMMA) {
            gretl_bundle_set_scalar(this, "value", 1);
            ret[i] = 1;
        } else {
            fprintf(stderr, "unknown distribution code '%s'\n", idstr);
        }
        ret[n] += ret[i];
    }

    return ret;
}

gretl_matrix *gibbs_via_bundles (gretl_bundle *B, int T,
                                 PRN *prn, int *err)
{
    gretl_matrix *ret = NULL;
    gretl_array *dists;
    gretl_bundle *dist;
    gretl_matrix *m;
    GretlType gt;
    double x;
    int *nparams;
    int i, j, k, t;
    int n;

    /* the array of distributions to sample */
    dists = gretl_bundle_get_array(B, "dists", err);
    if (*err) {
        return NULL;
    }

    /* the number of distributions */
    n = gretl_array_get_length(dists);

    /* initialize and count params */
    nparams = gibbs_init(B, prn, err);

    if (!*err) {
        /* allocate output matrix */
        ret = gretl_zero_matrix_new(T, nparams[n]);
    }

    gretl_iteration_push();

    for (t=0; t<T && !*err; t++) {
        k = 0;
        gt = 0;
        for (i=0; i<n && !*err; i++) {
            *err = gen_one(B, i, prn);
#if BDEBUG
            fprintf(stderr, "*** t = %d, i = %d ***\n", t, i);
            fprintf(stderr, " gen_one err %d\n", *err);
#endif
            if (!*err) {
                dist = gretl_array_get_data(dists, i);
                gt = gretl_bundle_get_member_type(dist, "value", NULL);
                if (gt == GRETL_TYPE_MATRIX) {
                    m = gretl_bundle_get_matrix(dist, "value", err);
                } else {
                    x = gretl_bundle_get_scalar(dist, "value", err);
                }
            }
            if (!*err && gt == GRETL_TYPE_MATRIX) {
                for (j=0; j<nparams[i]; j++) {
                    gretl_matrix_set(ret, t, k++, m->val[j]);
                }
            } else if (!*err) {
                gretl_matrix_set(ret, t, k++, x);
            }
        }
    }

    gretl_iteration_pop();

    free(nparams);
    if (*err) {
        gretl_matrix_free(ret);
        ret = NULL;
    }

    return ret;
}

/* End of apparatus pertaining to gibbs_via_bundles().
   Start of apparatus for "gibbs" block command.
*/

static GENERATOR **allocate_generators (int nv, int *err)
{
    GENERATOR **genrs = malloc(nv * sizeof *genrs);
    int i;

    if (genrs == NULL) {
        *err = E_ALLOC;
    } else {
        for (i=0; i<nv; i++) {
            genrs[i] = NULL;
        }
    }

    return genrs;
}

static int get_lhs_info (const char **psrc,
                         GretlType *pt,
                         char *vname)
{
    const char *src = *psrc;
    char s[VNAMELEN];
    GretlType t;

    sscanf(src, "%31[^= ]", s);
    t = gretl_type_from_string(s);
    if (t == GRETL_TYPE_NONE) {
        t = GRETL_TYPE_ANY;
        strcpy(vname, s);
    } else if (t != GRETL_TYPE_DOUBLE &&
               t != GRETL_TYPE_MATRIX) {
        return E_TYPES;
    } else {
        /* got 'scalar' or 'matrix': skip forward */
        src += strspn(src, " ");
        src += strlen(s);
        src += strspn(src, " ");
        sscanf(src, "%31[^= ]", vname);
    }

    *pt = t;
    *psrc = src;

    return 0;
}

static int gibbs_record_result (GENERATOR *genr,
                                gretl_matrix *H,
                                int t, int *pc)
{
    GretlType gt = genr_get_output_type(genr);
    int c = *pc;
    int err = 0;

    if (gt == GRETL_TYPE_DOUBLE) {
        double gx = genr_get_output_scalar(genr);

        gretl_matrix_set(H, t, c++, gx);
    } else {
        gretl_matrix *gm = genr_get_output_matrix(genr);
        int j, nvals;

        if (gm == NULL) {
            err = E_TYPES;
        } else {
            nvals = gm->rows * gm->cols;
            for (j=0; j<nvals; j++) {
                gretl_matrix_set(H, t, c++, gm->val[j]);
            }
        }
    }

    *pc = c;

    return err;
}

static int should_record (char *vname, char **S, int n)
{
    int i;

    for (i=0; i<n; i++) {
        if (!strcmp(vname, S[i])) {
            return 1;
        }
    }

    return 0;
}

static gretl_matrix *gibbs_via_genrs (char **init, int ni,
                                      char **iterate, int ng,
                                      int burnin, int N,
                                      guint8 *record,
                                      char **recnames, int nr,
                                      PRN *prn, int *err)
{
    gretl_matrix *ret = NULL;
    gretl_matrix *gm = NULL;
    GENERATOR **genrs = NULL;
    char vname[VNAMELEN];
    GretlType gt;
    const char *s;
    int ncols = 0;
    int iters;
    int i, j, k, t, c;

    for (i=0; i<ni && !*err; i++) {
        s = init[i];
        *err = get_lhs_info(&s, &gt, vname);
        if (!*err) {
            *err = generate(s, NULL, gt, OPT_NONE, prn);
        }
        if (*err) {
            pprintf(prn, "generate failed on init[%d]\n", i+1);
            pprintf(prn, " > %s\n", init[i]);
        } else {
            user_var *uv = get_user_var_by_name(vname);
            int add_cols = should_record(vname, recnames, nr);

            gt = user_var_get_type(uv);
            if (gt == GRETL_TYPE_DOUBLE) {
                if (add_cols) {
                    ncols++;
                }
            } else if (gt == GRETL_TYPE_MATRIX) {
                if (add_cols) {
                    gm = user_var_get_value(uv);
                    ncols += gm->rows * gm->cols;
                }
            } else {
                pprintf(prn, "bad type from init[%d]\n", i+1);
                pprintf(prn, " > %s\n", init[i]);
                *err = E_TYPES;
            }
        }
    }

    if (*err) {
        return ret;
    }

    genrs = allocate_generators(ng, err);
    if (!*err) {
        ret = gretl_matrix_alloc(N, ncols);
        if (ret == NULL) {
            *err = E_ALLOC;
        }
    }

    iters = burnin + N;
    gretl_iteration_push();

    /* the main iteration */
    for (t=0; t<iters && !*err; t++) {
        k = t - burnin;
        c = 0;
        for (i=0; i<ng && !*err; i++) {
            if (t == 0) {
                s = iterate[i];
                *err = get_lhs_info(&s, &gt, vname);
                if (!*err) {
                    genrs[i] = genr_compile(s, NULL, gt, OPT_NONE,
                                            prn, err);
                }
            } else {
                /* already compiled */
                *err = execute_genr(genrs[i], NULL, prn);
            }
            if (*err) {
                pprintf(prn, "genr[%d] failed at iter %d\n> %s\n",
                        i, t, iterate[i]);
            } else if (k >= 0 && record[i]) {
                *err = gibbs_record_result(genrs[i], ret, k, &c);
            }
        }
    }

    gretl_iteration_pop();

    for (j=0; j<ng; j++) {
        destroy_genr(genrs[j]);
    }
    free(genrs);

    return ret;
}

static char **gibbs_lines;
static int gibbs_started;
static int gibbs_n_lines;
static int gibbs_burnin;
static int gibbs_N;
static char *gibbs_output;

static void gibbs_destroy (void)
{
    if (gibbs_lines != NULL) {
        strings_array_free(gibbs_lines, gibbs_n_lines);
        gibbs_lines = NULL;
    }
    if (gibbs_output != NULL) {
        free(gibbs_output);
        gibbs_output = NULL;
    }
    gibbs_started = 0;
    gibbs_n_lines = 0;
    gibbs_burnin = 0;
    gibbs_N = 0;
}

static int parse_gibbs_params (const char *s)
{
    const char *targ[] = {
        "burnin", "N", "output"
    };
    char p[VNAMELEN] = {0};
    int i, err = 0;

    for (i=0; i<3; i++) {
        s += strspn(s, " ");
        if (*s == '\0') {
            err = E_PARSE;
            break;
        }
        if (!sscanf(s, "%8[^ =]", p) || strcmp(p, targ[i])) {
            err = E_PARSE;
            break;
        }
        s += strlen(p);
        s += strspn(s, " ");
        if (*s != '=') {
            err = E_PARSE;
            break;
        }
        s += 1 + strspn(s, " ");
        p[0] = '\0';
        if (!sscanf(s, "%31[^ ]", p)) {
            err = E_PARSE;
            break;
        }
        if (i == 0) {
            gibbs_burnin = gretl_int_from_string(p, &err);
        } else if (i == 1) {
            gibbs_N = gretl_int_from_string(p, &err);
        } else {
            err = check_identifier(p);
            if (!err) {
                gibbs_output = gretl_strdup(p);
            }
        }
        s += strlen(p);
    }

    if (!err) {
        s += strspn(s, " ");
        if (*s != '\0') {
            err = E_PARSE;
        }
    }
    if (err) {
        fprintf(stderr, "> %s\n", s);
        gibbs_destroy();
    }

    return err;
}

int gibbs_block_start (const char *line, PRN *prn)
{
    int err = 0;

    if (gibbs_started) {
        gretl_errmsg_sprintf(_("%s: a block is already started"),
                             gretl_command_word(GIBBS));
        return E_DATA;
    }

    /* parse out burnin, N and output */
    err = parse_gibbs_params(line);

    if (!err) {
        gibbs_started = 1;
    }

    return err;
}

int gibbs_block_append (const char *line)
{
    int err = 0;

    if (!gibbs_started) {
        gretl_errmsg_sprintf(_("%s: no block is in progress"),
                             gretl_command_word(GIBBS));
        err = E_DATA;
    } else if (!string_is_blank(line)) {
        err = strings_array_add(&gibbs_lines, &gibbs_n_lines, line);
        if (err) {
            gibbs_destroy();
        }
    }

    return err;
}

int gibbs_execute (gretlopt opt, PRN *prn)
{
    gretl_matrix *H = NULL;
    char **init = NULL;
    char **iter = NULL;
    guint8 *record = NULL;
    char vname[VNAMELEN];
    char **recnames = NULL;
    char *s;
    int verbose;
    int ni = 0;
    int ng = 0;
    int nr = 0;
    int i, iniprev = 0;
    int err = 0;

    verbose = (opt & OPT_V)? 1 : 0;

    if (verbose) {
        pprintf(prn, "gibbs: burnin = %d, N = %d, output %s\n",
                gibbs_burnin, gibbs_N, gibbs_output);
    }

    for (i=0; i<gibbs_n_lines && !err; i++) {
        s = gibbs_lines[i];
        if (!strncmp(s, "init: ", 6)) {
            if (i > 0 && !iniprev) {
                pprintf(prn, "Initializers must come first\n");
                pprintf(prn, "> %s\n", s);
                err = E_INVARG;
            } else {
                shift_string_left(s, 6);
                iniprev = 1;
                ni++;
            }
        } else {
            if (ng == 0) {
                ng = gibbs_n_lines - ni;
                record = calloc(ng, 1);
            }
            if (!strncmp(s, "record: ", 8)) {
                record[i-ni] = 1;
                shift_string_left(s, 8);
                s += strspn(s, " ");
                sscanf(s, "%31[^ =]", vname);
                strings_array_add(&recnames, &nr, vname);
            }
            iniprev = 0;
        }
    }

    if (!err) {
        if (ni > 0) {
            init = gibbs_lines;
        }
        iter = gibbs_lines + ni;
    }

    if (!err) {
        H = gibbs_via_genrs(init, ni, iter, ng,
                            gibbs_burnin, gibbs_N,
                            record, recnames, nr,
                            prn, &err);
        if (H != NULL) {
            err = user_var_add_or_replace(gibbs_output,
                                          GRETL_TYPE_MATRIX,
                                          H);
        }
        if (!err && verbose) {
            pprintf(prn, "%s is %d x %d\n", gibbs_output, H->rows, H->cols);
        }
    }

    strings_array_free(recnames, nr);
    free(record);
    gibbs_destroy();

    return err;
}
