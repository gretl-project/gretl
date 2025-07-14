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

static gretl_matrix *extract_param (const char *s, int i,
                                    gretl_bundle *allInfo)
{
    gretl_array *dists;
    gchar *tmp = g_strdup(s);
    gretl_matrix *ret = NULL;
    int n, k, err = 0;

    dists = gretl_bundle_get_array(allInfo, "dists", &err);

    g_strstrip(tmp);
    n = strlen(tmp);
    if (*tmp == '[') {
        if (tmp[n-1] == ']') {
            gretl_bundle *b = gretl_array_get_data(dists, i);
            gretl_array *a = gretl_bundle_get_array(b, "private", &err);

            sscanf(tmp+1, "%d]", &k);
            ret = gretl_array_get_data(a, k-1);
        }
    } else if (*tmp == '(') {
        if (tmp[n-1] == ')') {
            gretl_bundle *b;

            sscanf(tmp+1, "%d)", &k);
            b = gretl_array_get_data(dists, k-1);
            ret = gretl_bundle_get_matrix(b, "value", &err); /* copy? */
        }
    } else if (!strncmp(tmp, "__", 2)) {
        /* a user-written function? */
#if 0
        gchar *call = g_strdup_printf("matrix XXX");
        err = generate(line, dset, GRETL_TYPE_MATRIX, OPT_NONE,
                       prn);
#else        
        gchar *call = g_strdup_printf("%s(allInfo)", tmp + 2);

        ret = generate_matrix(call, NULL, &err);
        fprintf(stderr, "HERE '%s' genr err %d\n", call, err);
        g_free(call);
    }

    if (ret == NULL) {
        ret = gretl_matrix_alloc(1,1);
        ret->val[0] = NADBL;
    }

    g_free(tmp);

    return ret;
}

#define mapstr(m,i) gretl_array_get_data(m,i)

static double gen_one (gretl_bundle *B, int which)
{
    gretl_array *dists;
    gretl_bundle *this;
    gretl_array *map;
    const char *id;
    int err = 0;

    dists = gretl_bundle_get_array(B, "dists", &err);
    this = gretl_array_get_data(dists, which);
    id = gretl_bundle_get_string(this, "dist", &err);
    map = gretl_bundle_get_array(this, "map", &err);

    if (!strcmp(id, "N") || !strcmp(id, "N2")) {
        int N2 = strcmp(id, "N2") == 0;
        gretl_matrix *md;
        gretl_matrix *m;
        gretl_matrix *S;
        gretl_matrix *C;
        gretl_matrix *V;
        gretl_matrix *Mn;
        int dim;

        md = extract_param(mapstr(map, 0), which, B);
        m  = extract_param(mapstr(map, 1), which, B);
        S  = extract_param(mapstr(map, 2), which, B);
        dim = (int) md->val[0];
        if (N2) {
            gretl_matrix *Sm;

            Sm = gretl_matrix_alloc(S->rows, m->cols);
            gretl_matrix_multiply(S, m, Sm);
            V = gretl_matrix_vectorize_new(Sm);
            gretl_matrix_free(Sm);
        } else {
            V = gretl_matrix_vectorize_new(m);
        }
        C = gretl_matrix_copy(S);
        err = gretl_matrix_cholesky_decomp(C);
        Mn = gretl_random_matrix_new(dim, 1, D_NORMAL);

        gretl_matrix_multiply_mod(C,  GRETL_MOD_NONE,
                                  Mn, GRETL_MOD_NONE,
                                  V,  GRETL_MOD_CUMULATE);
        /* FIXME efficiency! */
        gretl_matrix_transpose_in_place(V);
        gretl_bundle_donate_data(this, "value", V, GRETL_TYPE_MATRIX, 0);
    } else if (!strcmp(id, "U")) {
        gretl_matrix *ma, *mb;
        double d, parm[2];

        ma = extract_param(mapstr(map, 0), which, B);
        mb = extract_param(mapstr(map, 1), which, B);
        parm[0] = ma->val[0]; parm[1] = mb->val[0];
        d = gretl_get_random_scalar(D_UNIFORM, parm, &err);
        gretl_bundle_set_scalar(this, "value", d);
        /* free stuff? */
    } else if (!strcmp(id, "B")) {
        gretl_matrix *mp;
        double parm[2] = {0, 1};
        double p, d;

        mp = extract_param(mapstr(map, 0), which, B);
        p = mp->val[0];
        d = gretl_get_random_scalar(D_UNIFORM, parm, &err);
        gretl_bundle_set_scalar(this, "value", d < p);
    } else if (!strcmp(id, "G") || !strcmp(id, "IG")) {
        gretl_matrix *mp, *ma;
        double d, parm[2];

        mp = extract_param(mapstr(map, 0), which, B);
        ma = extract_param(mapstr(map, 1), which, B);
        parm[0] = mp->val[0];
        parm[1] = ma->val[0];
        d = gretl_get_random_scalar(D_GAMMA, parm, &err);
        if (!strcmp(id, "G")) {
            gretl_bundle_set_scalar(this, "value", d);
        } else {
            gretl_bundle_set_scalar(this, "value", 1/d);
        }
    }

    return err;
}

/* initialises the array; returns a int array with the
   dimension of each element */

static int *gibbs_init (gretl_bundle *allInfo, int *err)
{
    gretl_array *dists;
    gretl_bundle *this;
    const char *id;
    int *ret;
    int n, i;

    dists = gretl_bundle_get_array(allInfo, "dists", err);
    if (*err) {
        return NULL;
    }

    n = gretl_array_get_length(dists);
    ret = malloc(n * sizeof *ret);

    for (i=0; i<n && !*err; i++) {
        ret[i] = 0;
        this = gretl_array_get_data(dists, i);
        id = gretl_bundle_get_string(this, "dist", err);
        if (!strcmp(id, "N") || !strcmp(id, "N2")) {
            gretl_array *map;
            gretl_matrix *md;
            gretl_matrix *Mn;
            int dim;

            map = gretl_bundle_get_array(this, "map", err);
            md = extract_param(mapstr(map, 0), i, allInfo);
            fprintf(stderr, "md: %p\n", (void *) md);
            gretl_matrix_print(md, "md");
            dim = md->val[0];
            Mn = gretl_random_matrix_new(dim, 1, D_NORMAL);
            gretl_bundle_donate_data(this, "value", Mn, GRETL_TYPE_MATRIX, 0);
            ret[i] = dim;
        } else if (!strcmp(id, "U")) {
            gretl_bundle_set_scalar(this, "value", 0.5);
            ret[i] = 1;
        } else if (!strcmp(id, "B")) {
            double parm[2] = {0, 1};
            double d;

            d = gretl_get_random_scalar(D_UNIFORM, parm, err);
            gretl_bundle_set_scalar(this, "value", d > 0.5);
            ret[i] = 1;
        } else if (!strcmp(id, "G") || !strcmp(id, "IG")) {
            gretl_bundle_set_scalar(this, "value", 1);
            ret[i] = 1;
        }
        fprintf(stderr, "ret[%d] = %d\n", i, ret[i]);
    }

    return ret;
}

gretl_matrix *gibbs_via_bundles (gretl_bundle *B, int T,
                                 PRN *prn, int *err)
{
    gretl_matrix *ret = NULL;
    gretl_array *dists;
    gretl_bundle *localB;
    gretl_bundle *dist;
    gretl_matrix *m;
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
    pprintf(prn, "gibbs_via_bundles: n=%d, T=%d\n", n, T);


    localB = gretl_bundle_copy(B, err);
    if (!*err) {
        nparams = gibbs_init(localB, err);
    }

    if (!*err) {
        int c = 0;

        for (i=0; i<n; i++) {
            c += nparams[i];
        }
        ret = gretl_zero_matrix_new(T, c);
        pprintf(prn, " output cols=%d\n", c);
    }

    for (t=0; t<T && !*err; t++) {
        pprintf(prn, " t = %d\n", t);
        k = 0;
        for (i=0; i<n && !*err; i++) {
            pprintf(prn, "   i = %d\n", i);
            *err = gen_one(localB, i);
            dist = gretl_array_get_data(dists, i);
            m = gretl_bundle_get_matrix(dist, "value", err);
            for (j=0; j<nparams[i]; j++) {
                gretl_matrix_set(ret, t, k++, m->val[j]);
            }
        }
    }

    pprintf(prn, "gibbs_via_bundles: *err = %d\n", *err);

    free(nparams);
    gretl_bundle_destroy(localB);

    return ret;
}
