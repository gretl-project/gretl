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
#include "gretl_sampler.h"

/* file-scope globals */
static char **gibbs_lines;
static char **gibbs_temps;
static int gibbs_started;
static int gibbs_n_lines;
static int gibbs_n_temps;
static int gibbs_burnin;
static int gibbs_N;
static char *gibbs_output;

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

/* Write the result from the compiled @genr into row @s of the matrix
   @H, starting at column c = *pc. Pass back in @pc the column that was
   reached.
*/

static int gibbs_record_result (GENERATOR *genr,
                                int recsize,
                                gretl_matrix *H,
                                int s, int *pc)
{
    GretlType gt = genr_get_output_type(genr);
    int c = *pc;
    int err = 0;

    if (gt == GRETL_TYPE_DOUBLE) {
        double gx = genr_get_output_scalar(genr);

        gretl_matrix_set(H, s, c++, gx);
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
                    gretl_matrix_set(H, s, c++, gm->val[j]);
                }
            }
        }
    }

    *pc = c;

    return err;
}

static gretl_matrix *do_run_sampler (char **init, int ni,
                                     char **iter, int ng,
                                     int burnin, int N,
                                     guint8 *record,
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
    int cleanup = (opt & OPT_C);
    int verbose = (opt & OPT_V);
    int ncols = 0;
    int iters;
    int i, j, s, t, c;

    /* run the @init generators and ensure that they all produce
       a scalar or matrix result
    */
    if (verbose && ni > 0) {
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

    /* Compile and run the @iter generators. We don't record the results
       here, but we check that the genrs work, and for those that are
       tagged for recording we store the sizes of the objects generated.
    */
    if (verbose) {
        pputs(prn, "checking iteration statements\n");
    }
    for (i=0; i<ng && !*err; i++) {
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

            gt = user_var_get_type(uv);
            if (gt == GRETL_TYPE_DOUBLE) {
                if (record[i]) {
                    ncols++;
                }
            } else if (gt == GRETL_TYPE_MATRIX) {
                if (record[i]) {
                    gm = user_var_get_value(uv);
                    ncols += gm->rows * gm->cols;
                    record[i] = ncols;
                }
            } else {
                pprintf(prn, "bad type from iter[%d]\n", i+1);
                pprintf(prn, " > %s\n", iter[i]);
                *err = E_TYPES;
            }
        }
    }

    if (*err) {
        goto bailout;
    }

    ret = gretl_matrix_alloc(N, ncols);
    if (ret == NULL) {
        *err = E_ALLOC;
        goto bailout;
    }

    iters = burnin + N;
    if (verbose) {
        pputs(prn, "starting iteration\n");
    }
    gretl_iteration_push();

    /* the main iteration, using compiled genrs */
    for (t=0; t<iters && !*err; t++) {
        s = t - burnin;
        c = 0;
        for (i=0; i<ng && !*err; i++) {
            *err = execute_genr(genrs[i], dset, prn);
            if (*err) {
                pprintf(prn, "genr[%d] failed at iter %d\n> %s\n",
                        i, t, iter[i]);
            } else if (s >= 0 && record[i]) {
                *err = gibbs_record_result(genrs[i], record[i],
                                           ret, s, &c);
            }
        }
    }

    gretl_iteration_pop();
    if (verbose && !*err) {
        pputs(prn, "iteration completed\n");
    }

 bailout:

    for (j=0; j<ng; j++) {
        destroy_genr(genrs[j]);
    }
    free(genrs);

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

int gibbs_execute (gretlopt opt, DATASET *dset, PRN *prn)
{
    gretl_matrix *H = NULL;
    char **init = NULL;
    char **iter = NULL;
    guint8 *record = NULL;
    char *str;
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
        gretl_errmsg_sprintf("%s: no statements to be interated", "gibbs");
        err = E_ARGS;
    }
    if (!err && nr == 0) {
        gretl_errmsg_sprintf("%s: nothing to be recorded", "gibbs");
        err = E_ARGS;
    }
    if (!err) {
        if (ni > 0) {
            init = gibbs_lines;
        }
        iter = gibbs_lines + ni;
    }

    if (!err) {
        H = do_run_sampler(init, ni, iter, ng,
                           gibbs_burnin, gibbs_N,
                           record, opt, dset, prn, &err);
        if (H != NULL) {
            err = user_var_add_or_replace(gibbs_output,
                                          GRETL_TYPE_MATRIX,
                                          H);
        }
        if (!err && verbose) {
            pprintf(prn, "%s is %d x %d\n", gibbs_output, H->rows, H->cols);
        }
    }

    free(record);
    gibbs_destroy();

    return err;
}
