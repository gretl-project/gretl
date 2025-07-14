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

#if 0
    fprintf(stderr, "get_lhs_info: vname '%s', type '%s'\n",
            vname, gretl_type_get_name(t));
#endif

    return 0;
}

static int is_inivar (const char *s, char **S, int ns)
{
    int i;

    for (i=0; i<ns; i++) {
        if (!strcmp(s, S[i])) {
            return 1;
        }
    }

    return 0;
}

gretl_matrix *gibbs_via_genrs (char **init, int ni,
                               char **iterate, int ng,
                               int burnin, int N,
                               PRN *prn, int *err)
{
    gretl_matrix *ret = NULL;
    gretl_matrix *gm = NULL;
    GENERATOR **genrs = NULL;
    GretlType *gtypes = NULL;
    GretlType gt;
    char **inivars = NULL;
    char vname[VNAMELEN];
    const char *s;
    gint8 *record = NULL;
    int ncols = 0;
    int nvals;
    int iters;
    double gx;
    int i, j, k, t, c;

    gtypes = calloc(ni, sizeof *gtypes);
    record = calloc(ng, sizeof *record);
    inivars = strings_array_new(ni);

    for (i=0; i<ni && !*err; i++) {
        s = init[i];
        *err = get_lhs_info(&s, &gt, vname);
        if (!*err) {
            gtypes[i] = gt;
            inivars[i] = gretl_strdup(vname);
            *err = generate(s, NULL, gt, OPT_NONE, prn);
        }
        if (*err) {
            fprintf(stderr, "generate failed on init[%d]\n", i);
            fprintf(stderr, " > %s\n", init[i]);
        } else {
            user_var *uv = get_user_var_by_name(vname);

            gtypes[i] = user_var_get_type(uv);
            if (gtypes[i] == GRETL_TYPE_DOUBLE) {
                ncols++;
            } else if (gtypes[i] == GRETL_TYPE_MATRIX) {
                gm = user_var_get_value(uv);
                ncols += gm->rows * gm->cols;
            } else {
                fprintf(stderr, "E_TYPES on init[%d]\n", i);
                *err = E_TYPES;
            }
        }
    }

    if (*err) {
        free(record);
        free(gtypes);
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
        k = c = 0;
        for (i=0; i<ng && !*err; i++) {
            if (t == 0) {
                s = iterate[i];
                *err = get_lhs_info(&s, &gt, vname);
                if (!*err) {
                    genrs[i] = genr_compile(s, NULL, gt, OPT_NONE,
                                            prn, err);
                }
                if (!*err && is_inivar(vname, inivars, ni)) {
                    record[i] = 1;
                }
            } else {
                /* already compiled */
                *err = execute_genr(genrs[i], NULL, prn);
            }
            if (*err) {
                fprintf(stderr, "genr[%d] failed at iter %d\n> %s\n",
                        i, t, iterate[i]);
            } else if (t >= burnin && record[i]) {
                if (gtypes[k] == GRETL_TYPE_DOUBLE) {
                    gx = genr_get_output_scalar(genrs[i]);
                    gretl_matrix_set(ret, t, c++, gx);
                } else {
                    gm = genr_get_output_matrix(genrs[i]);
                    if (gm == NULL) {
                        fprintf(stderr, "genr[%d] didn't produce a matrix\n", i);
                        *err = 1;
                    } else {
                        nvals = gm->rows * gm->cols;
                        for (j=0; j<nvals; j++) {
                            gretl_matrix_set(ret, t, c++, gm->val[j]);
                        }
                    }
                }
                k++;
            }
        }
    }

    gretl_iteration_pop();

    for (j=0; j<ng; j++) {
        destroy_genr(genrs[j]);
    }
    free(genrs);
    free(record);
    free(gtypes);
    strings_array_free(inivars, ni);

    return ret;
}

/* Below: apparatus pertaining to gibbs_via_bundles() */

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

    free(nparams);
    if (*err) {
        gretl_matrix_free(ret);
        ret = NULL;
    }

    return ret;
}
