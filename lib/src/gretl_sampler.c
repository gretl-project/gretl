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

/* Apparatus to support the experimental gibbs block command */

#include "libgretl.h"
#include "gretl_typemap.h"
#include "uservar.h"
#include "libset.h"
#include "matrix_extra.h"
#include "gretl_sampler.h"

/* file-scope globals */
static char **gibbs_lines;
static char **gibbs_temps;
static int gibbs_started;
static int gibbs_n_lines;
static int gibbs_n_temps;
static int gibbs_burnin;
static int gibbs_N;
static int gibbs_thin;
static char *gibbs_output;

typedef struct gibbs_var_info_ {
    char *name;
    GretlType gt;
    int startcol;
    int ncols;
    int discrete;
} gibbs_var_info;

static gibbs_var_info *gibbs_info_alloc (int nr)
{
    gibbs_var_info *gvi = malloc(nr * sizeof *gvi);
    int i;

    for (i=0; i<nr; i++) {
        gvi[i].name = NULL;
        gvi[i].startcol = 0;
        gvi[i].ncols = 0;
        gvi[i].discrete = 0;
    }

    return gvi;
}

static void gibbs_info_destroy (gibbs_var_info *gvi, int nr)
{
    if (gvi != NULL) {
        int i;

        for (i=0; i<nr; i++) {
            free(gvi[i].name);
        }
        free(gvi);
    }
}

static void gibbs_delete_temps (void)
{
    int i;

    for (i=0; i<gibbs_n_temps; i++) {
        user_var_delete_by_name(gibbs_temps[i], NULL);
    }
}

static void gibbs_destroy (void)
{
    if (gibbs_lines != NULL) {
        strings_array_free(gibbs_lines, gibbs_n_lines);
        gibbs_lines = NULL;
    }
    if (gibbs_temps != NULL) {
        gibbs_delete_temps();
        strings_array_free(gibbs_temps, gibbs_n_temps);
        gibbs_temps = NULL;
    }
    if (gibbs_output != NULL) {
        free(gibbs_output);
        gibbs_output = NULL;
    }
    gibbs_started = 0;
    gibbs_n_lines = 0;
    gibbs_n_temps = 0;
    gibbs_burnin = 0;
    gibbs_N = 0;
    gibbs_thin = 0;
}

static void gibbs_mark_as_temp (const char *s)
{
    strings_array_add(&gibbs_temps, &gibbs_n_temps, s);
}

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

/* For a statement within a gibbs block: extract a leading type name
   ('scalar' or 'matrix') if present and thereby set the content of
   @pt. Also extract the identifier of the left-hand side object and
   write it to @vname. Pass back in @psrc the portion of the statement
   starting with @vname.
*/

static int gibbs_get_lhs_info (const char **psrc,
                               GretlType *pt,
                               char *vname)
{
    const char *src = *psrc;
    char s[VNAMELEN] = {0};
    GretlType t;

    sscanf(src, "%31[^= ]", s);
    t = gretl_type_from_string(s);
    if (t == GRETL_TYPE_NONE) {
        if (strlen(s) == gretl_namechar_spn(s)) {
            strcpy(vname, s);
        }
        t = GRETL_TYPE_ANY;
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

/* Write the result from the compiled @genr into row @k of the matrix
   @H, starting at column c = *pc. Pass back in @pc the column that was
   reached.
*/

static int gibbs_record_result (GENERATOR *genr,
                                int recsize,
                                gretl_matrix *H,
                                int k, int *pc)
{
    GretlType gt = genr_get_output_type(genr);
    int c = *pc;
    int err = 0;

    if (gt == GRETL_TYPE_DOUBLE) {
        double gx = genr_get_output_scalar(genr);

        gretl_matrix_set(H, k, c++, gx);
    } else {
        gretl_matrix *gm = genr_get_output_matrix(genr);
        int j, nvals;

        if (gm == NULL) {
            err = E_TYPES;
        } else {
            nvals = gm->rows * gm->cols;
            if (nvals != recsize) {
                err = E_DATA;
            } else {
                for (j=0; j<nvals; j++) {
                    gretl_matrix_set(H, k, c++, gm->val[j]);
                }
            }
        }
    }

    *pc = c;

    return err;
}

static void matrix_dset_set_names (DATASET *mdset, int *list)
{
    int i;

    for (i=1; i<=list[0]; i++) {
        sprintf(mdset->varname[i], "col%d", list[i]);
    }
}

static void gibbs_attach_param_stats (gretl_bundle *GB,
                                      gretl_matrix *m,
                                      const char *pname)
{
    gretl_bundle *b;

    b = gretl_bundle_get_bundle(GB, pname, NULL);
    if (b != NULL) {
        gretl_bundle_donate_data(b, "stats", m, GRETL_TYPE_MATRIX, 0);
    }
}

static void gibbs_param_info (gibbs_var_info *gvi, int nr,
                              gretl_matrix *M,
                              gretl_bundle *GB,
                              PRN *prn)
{
    DATASET *mdset;
    gretl_matrix *m;
    GretlType gt;
    gretlopt opt = OPT_S;
    int *collist;
    int i, c1;
    int err = 0;

    if (prn == NULL) {
        opt |= OPT_Q;
    } else {
        pprintf(prn, "variables recorded in %s:\n\n", gibbs_output);
    }

    for (i=0; i<nr; i++) {
        c1 = gvi[i].startcol + 1;
        if (gvi[i].gt == GRETL_TYPE_DOUBLE) {
            /* a scalar variable */
            int list[2] = {1, c1};

            if (prn != NULL) {
                pprintf(prn, "scalar %s", gvi[i].name);
            }
            mdset = gretl_dataset_from_matrix(M, list, OPT_B, &err);
            if (!err) {
                matrix_dset_set_names(mdset, list);
                collist = gretl_consecutive_list_new(1, mdset->v - 1);
                err = list_summary(collist, 0, mdset, opt | OPT_M, prn);
                if (!err) {
                    m = get_last_result_data(&gt, &err);
                    if (m != NULL) {
                        gibbs_attach_param_stats(GB, m, gvi[i].name);
                    }
                }
                free(collist);
            }
            destroy_dataset(mdset);
        } else {
            /* a matrix variable */
            int c2 = c1 + gvi[i].ncols - 1;
            int *list = gretl_consecutive_list_new(c1, c2);

            if (prn != NULL) {
                pprintf(prn, "matrix %s", gvi[i].name);
            }
            mdset = gretl_dataset_from_matrix(M, list, OPT_B, &err);
            if (!err) {
                matrix_dset_set_names(mdset, list);
                collist = gretl_consecutive_list_new(1, mdset->v - 1);
                err = list_summary(collist, 0, mdset, opt, prn);
                if (!err) {
                    m = get_last_result_data(&gt, &err);
                    if (m != NULL) {
                        gibbs_attach_param_stats(GB, m, gvi[i].name);
                    }
                }
                free(collist);
            }
            destroy_dataset(mdset);
            free(list);
        }
    }
    if (prn != NULL) {
        pputc(prn, '\n');
    }
}

static void gibbs_record_info (gibbs_var_info *gvi,
                               const char *vname,
                               GretlType gt,
                               int c, int recsize)
{
    gvi->name = gretl_strdup(vname);
    gvi->gt = gt;
    gvi->startcol = c;
    gvi->ncols = recsize;
}

static gretl_matrix *do_run_sampler (char **init, int ni,
                                     char **iter, int ng,
                                     guint8 *record,
                                     gibbs_var_info *gvi,
                                     gretlopt opt,
                                     DATASET *dset,
                                     PRN *prn,
                                     int *err)
{
    gretl_matrix *ret = NULL;
    gretl_matrix *gm = NULL;
    GENERATOR **genrs = NULL;
    char vname[VNAMELEN];
    GretlType gt;
    const char *str;
    int cleanup = !(opt & OPT_K);
    int quiet = (opt & OPT_Q);
    int nrec = gibbs_N;
    int ncols = 0;
    int msize = 0;
    int iters;
    int i, j, k, s, t, c;

    if (gibbs_thin > 0) {
	nrec /= gibbs_thin;
    }

    /* run the @init generators and ensure that they all produce
       a scalar or matrix result
    */
    if (!quiet && ni > 0) {
        pputs(prn, "initializing\n");
    }
    for (i=0; i<ni && !*err; i++) {
        str = init[i];
        *err = gibbs_get_lhs_info(&str, &gt, vname);
        if (!*err) {
            if (cleanup && !gretl_is_user_var(vname)) {
                gibbs_mark_as_temp(vname);
            }
            *err = generate(str, dset, gt, OPT_NONE, prn);
        }
        if (*err) {
            pprintf(prn, "generate failed on init[%d]\n", i+1);
            pprintf(prn, " > %s\n", init[i]);
        } else {
            user_var *uv = get_user_var_by_name(vname);

            gt = user_var_get_type(uv);
            if (gt != GRETL_TYPE_DOUBLE && gt != GRETL_TYPE_MATRIX) {
                pprintf(prn, "bad type from init[%d]\n", i+1);
                pprintf(prn, " > %s\n", init[i]);
                *err = E_TYPES;
            }
        }
    }

    if (!*err) {
        genrs = allocate_generators(ng, err);
    }
    if (*err) {
        goto bailout;
    }

    /* Compile and run the @iter generators. We don't record the results
       here, but we check that the genrs work, and for those that are
       tagged for recording we store the sizes of the objects generated.
    */
    if (!quiet) {
        pputs(prn, "checking iteration statements\n");
    }
    for (i=0, j=0; i<ng && !*err; i++) {
        str = iter[i];
        *err = gibbs_get_lhs_info(&str, &gt, vname);
        if (!*err) {
            if (cleanup && !gretl_is_user_var(vname)) {
                gibbs_mark_as_temp(vname);
            }
            genrs[i] = genr_compile(str, dset, gt, OPT_NONE,
                                    prn, err);
        }
        if (*err) {
            pprintf(prn, "generate failed on iter[%d]\n", i+1);
            pprintf(prn, " > %s\n", iter[i]);
        } else {
            user_var *uv = get_user_var_by_name(vname);
            int colnum = ncols;

            gt = user_var_get_type(uv);
            if (gt != GRETL_TYPE_DOUBLE && gt != GRETL_TYPE_MATRIX) {
                pprintf(prn, "bad type from iter[%d]\n", i+1);
                pprintf(prn, " > %s\n", iter[i]);
                *err = E_TYPES;
            } else if (record[i]) {
                if (gt == GRETL_TYPE_DOUBLE) {
                    ncols++;
                } else {
                    gm = user_var_get_value(uv);
                    msize = gm->rows * gm->cols;
                    ncols += msize;
                    record[i] = msize;
                }
                gibbs_record_info(&gvi[j++], vname, gt, colnum, record[i]);
            }
        }
    }

    if (*err) {
        goto bailout;
    }

    ret = gretl_matrix_alloc(nrec, ncols);
    if (ret == NULL) {
        *err = E_ALLOC;
        goto bailout;
    }

    iters = gibbs_burnin + gibbs_N;
    if (!quiet) {
        pprintf(prn, "starting %d iterations\n", iters);
    }
    gretl_iteration_push();

    /* the main iteration, using compiled genrs */
    for (t=0, k=0; t<iters && !*err; t++) {
        int discard = 0;

        s = t - gibbs_burnin;
	if (gibbs_thin > 0 && s >= 0) {
	    discard = s % gibbs_thin > 0;
	}
        c = 0;
        for (i=0; i<ng && !*err; i++) {
            *err = execute_genr(genrs[i], dset, prn);
            if (*err) {
                pprintf(prn, "genr[%d] failed at iter %d\n> %s\n",
                        i, t, iter[i]);
            } else if (s >= 0 && record[i] && !discard) {
                *err = gibbs_record_result(genrs[i], record[i],
                                           ret, k, &c);
            }
        }
        if (s >= 0 && !discard) {
            k++;
        }
    }

    gretl_iteration_pop();
    if (!quiet && !*err) {
        pputs(prn, "iteration completed\n");
    }

 bailout:

    if (genrs != NULL) {
        for (j=0; j<ng; j++) {
            destroy_genr(genrs[j]);
        }
        free(genrs);
    }

    return ret;
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
    } else if (line == NULL || *line == '\0') {
        return E_ARGS;
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

static char *gibbs_type_name (GretlType gt)
{
    if (gt == GRETL_TYPE_DOUBLE) {
        return gretl_strdup("scalar");
    } else {
        return gretl_strdup("matrix");
    }
}

static gretl_bundle *make_gibbs_bundle (gretl_matrix *H,
                                        gibbs_var_info *gvi,
                                        int nr, int *err)
{
    gretl_bundle *b = gretl_bundle_new();
    gretl_bundle *bi = NULL;
    int i;

    if (b == NULL) {
        *err = E_ALLOC;
        return NULL;
    }

    gretl_bundle_donate_data(b, "H", H, GRETL_TYPE_MATRIX, 0);
    /* general metadata */
    gretl_bundle_set_int(b, "burnin", gibbs_burnin);
    gretl_bundle_set_int(b, "iterations", gibbs_N);
    gretl_bundle_set_int(b, "thinning", gibbs_thin);

    for (i=0; i<nr; i++) {
        /* per-variable metadata */
        bi = gretl_bundle_new();
        gretl_bundle_donate_data(bi, "type", gibbs_type_name(gvi[i].gt),
                                 GRETL_TYPE_STRING, 0);
        gretl_bundle_set_int(bi, "startcol", gvi[i].startcol + 1);
        gretl_bundle_set_int(bi, "ncols", gvi[i].ncols);
        gretl_bundle_set_int(bi, "discrete", gvi[i].discrete);
        gretl_bundle_donate_data(b, gvi[i].name, bi,
                                 GRETL_TYPE_BUNDLE, 0);
    }

    return b;
}

int gibbs_execute (gretlopt opt, DATASET *dset, PRN *prn)
{
    gibbs_var_info *gvi = NULL;
    gretl_bundle *GB = NULL;
    gretl_matrix *H = NULL;
    char **init = NULL;
    char **iter = NULL;
    guint8 *record = NULL;
    char *str;
    int quiet;
    int ni = 0;
    int ng = 0;
    int nr = 0;
    int i, iniprev = 0;
    int err = 0;

    if (opt & OPT_T) {
	int thin = get_optval_int(GIBBS, OPT_T, &err);

	if (err) {
	    gibbs_destroy();
	    return err;
	} else {
	    gibbs_thin = thin;
	}
    }

    quiet = (opt & OPT_Q)? 1 : 0;

    if (!quiet) {
        pprintf(prn, "gibbs: burnin = %d, N = %d, thinning = %d, output %s\n",
                gibbs_burnin, gibbs_N, gibbs_thin, gibbs_output);
    }

    for (i=0; i<gibbs_n_lines && !err; i++) {
        str = gibbs_lines[i];
        if (!strncmp(str, "init ", 5)) {
            if (i > 0 && !iniprev) {
                pprintf(prn, "Initializers must come first\n");
                pprintf(prn, "> %s\n", str);
                err = E_INVARG;
            } else {
                shift_string_left(str, 5);
                iniprev = 1;
                ni++;
            }
        } else {
            if (ng == 0) {
                ng = gibbs_n_lines - ni;
                record = calloc(ng, 1);
            }
            if (!strncmp(str, "record ", 7)) {
                record[i-ni] = 1;
                shift_string_left(str, 7);
                nr++;
            }
            iniprev = 0;
        }
    }

    if (!err && ng == 0) {
        gretl_errmsg_sprintf("%s: no statements to be iterated", "gibbs");
        err = E_ARGS;
    }
    if (!err && (nr == 0 || gibbs_N < gibbs_thin)) {
        gretl_errmsg_sprintf("%s: nothing to be recorded", "gibbs");
        err = E_ARGS;
    }
    if (!err) {
        if (ni > 0) {
            init = gibbs_lines;
        }
        iter = gibbs_lines + ni;
        gvi = gibbs_info_alloc(nr);
    }

    if (!err) {
        H = do_run_sampler(init, ni, iter, ng,
                           record, gvi, opt, dset,
                           prn, &err);
        if (H != NULL) {
            GB = make_gibbs_bundle(H, gvi, nr, &err);
        }
        if (GB != NULL) {
            err = user_var_add_or_replace(gibbs_output,
                                          GRETL_TYPE_BUNDLE,
                                          GB);
        }
        if (!err) {
            if (quiet) {
                gibbs_param_info(gvi, nr, H, GB, NULL);
            } else {
                pprintf(prn, "results matrix is %d x %d\n", H->rows, H->cols);
                gibbs_param_info(gvi, nr, H, GB, prn);
            }
        }
    }

    free(record);
    gibbs_destroy();
    gibbs_info_destroy(gvi, nr);

    return err;
}
