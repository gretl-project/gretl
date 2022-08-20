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

/* prog_loop.c - code to support "progressive" loops */

static void loop_model_free (LOOP_MODEL *lmod)
{
    int i, n;

#if LOOP_DEBUG > 1
    fprintf(stderr, "loop_model_free: lmod at %p, model0 at %p\n",
            (void *) lmod, (void *) lmod->model0);
#endif

    n = 4 * lmod->model0->ncoeff;

    for (i=0; i<n; i++) {
        mpf_clear(lmod->bigarray[i]);
    }

    free(lmod->bigarray);
    free(lmod->cbak);
    free(lmod->cdiff);

    gretl_model_free(lmod->model0);
}

/* Reset the loop model */

static void loop_model_zero (LOOP_MODEL *lmod, int started)
{
    int i, bnc = 4 * lmod->nc;

#if LOOP_DEBUG > 1
    fprintf(stderr, "loop_model_zero: %p\n", (void *) lmod);
#endif

    for (i=0; i<bnc; i++) {
        if (started) {
            mpf_set_d(lmod->bigarray[i], 0.0);
        } else {
            mpf_init(lmod->bigarray[i]);
        }
    }

    for (i=0; i<lmod->nc; i++) {
        lmod->cbak[i] = lmod->sbak[i] = NADBL;
        lmod->cdiff[i] = lmod->sdiff[i] = 0;
    }

    lmod->n = 0;
}

/* Set everything in lmod to 0/null in case of failure */

static void loop_model_init (LOOP_MODEL *lmod, int lno)
{
    lmod->lineno = lno;
    lmod->nc = 0;
    lmod->model0 = NULL;
    lmod->bigarray = NULL;
    lmod->cbak = NULL;
    lmod->cdiff = NULL;
}

/* Start up a LOOP_MODEL struct: copy @pmod into place and
   allocate storage */

static int loop_model_start (LOOP_MODEL *lmod, MODEL *pmod)
{
    int nc = pmod->ncoeff;
    int err = 0;

#if LOOP_DEBUG > 1
    fprintf(stderr, "init: copying model at %p\n", (void *) pmod);
#endif

    lmod->model0 = gretl_model_copy(pmod);
    if (lmod->model0 == NULL) {
        return E_ALLOC;
    }

    lmod->nc = nc;

    lmod->bigarray = malloc(nc * 4 * sizeof *lmod->bigarray);
    if (lmod->bigarray == NULL) {
        return E_ALLOC;
    }

    lmod->sum_coeff = lmod->bigarray;
    lmod->ssq_coeff = lmod->sum_coeff + nc;
    lmod->sum_sderr = lmod->ssq_coeff + nc;
    lmod->ssq_sderr = lmod->sum_sderr + nc;

    lmod->cbak = malloc(nc * 2 * sizeof *lmod->cbak);
    if (lmod->cbak == NULL) {
        err = E_ALLOC;
    } else {
        lmod->sbak = lmod->cbak + nc;
    }

    if (!err) {
        lmod->cdiff = malloc(nc * 2 * sizeof *lmod->cdiff);
        if (lmod->cdiff == NULL) {
            err = E_ALLOC;
        } else {
            lmod->sdiff = lmod->cdiff + nc;
        }
    }

    if (!err) {
        loop_model_zero(lmod, 0);
#if LOOP_DEBUG > 1
        fprintf(stderr, " model copied to %p, returning 0\n",
                (void *) lmod->model0);
#endif
    }

    if (err) {
        free(lmod->bigarray);
        free(lmod->cbak);
        free(lmod->cdiff);
    }

    return err;
}

static void loop_print_free (LOOP_PRINT *lprn)
{
    int i;

    for (i=0; i<lprn->nvars; i++) {
        mpf_clear(lprn->sum[i]);
        mpf_clear(lprn->ssq[i]);
    }

    strings_array_free(lprn->names, lprn->nvars);

    free(lprn->sum);
    free(lprn->ssq);
    free(lprn->xbak);
    free(lprn->diff);
    free(lprn->na);
}

static void loop_print_zero (LOOP_PRINT *lprn, int started)
{
    int i;

    lprn->n = 0;

    for (i=0; i<lprn->nvars; i++) {
        if (started) {
            mpf_set_d(lprn->sum[i], 0.0);
            mpf_set_d(lprn->ssq[i], 0.0);
        } else {
            mpf_init(lprn->sum[i]);
            mpf_init(lprn->ssq[i]);
        }
        lprn->xbak[i] = NADBL;
        lprn->diff[i] = 0;
        lprn->na[i] = 0;
    }
}

/* allocate and initialize @lprn, based on the number of
   elements in @namestr */

static int loop_print_start (LOOP_PRINT *lprn, const char *namestr)
{
    int i, nv;

    if (namestr == NULL || *namestr == '\0') {
        gretl_errmsg_set("'print' list is empty");
        return E_DATA;
    }

    lprn->names = gretl_string_split(namestr, &lprn->nvars, NULL);
    if (lprn->names == NULL) {
        return E_ALLOC;
    }

    nv = lprn->nvars;

    for (i=0; i<nv; i++) {
        if (!gretl_is_scalar(lprn->names[i])) {
            gretl_errmsg_sprintf(_("'%s': not a scalar"), lprn->names[i]);
            strings_array_free(lprn->names, lprn->nvars);
            lprn->names = NULL;
            lprn->nvars = 0;
            return E_DATA;
        }
    }

    lprn->sum = malloc(nv * sizeof *lprn->sum);
    if (lprn->sum == NULL) goto cleanup;

    lprn->ssq = malloc(nv * sizeof *lprn->ssq);
    if (lprn->ssq == NULL) goto cleanup;

    lprn->xbak = malloc(nv * sizeof *lprn->xbak);
    if (lprn->xbak == NULL) goto cleanup;

    lprn->diff = malloc(nv * sizeof *lprn->diff);
    if (lprn->diff == NULL) goto cleanup;

    lprn->na = malloc(nv);
    if (lprn->na == NULL) goto cleanup;

    loop_print_zero(lprn, 0);

    return 0;

 cleanup:

    strings_array_free(lprn->names, lprn->nvars);
    lprn->names = NULL;
    lprn->nvars = 0;

    free(lprn->sum);
    free(lprn->ssq);
    free(lprn->xbak);
    free(lprn->diff);
    free(lprn->na);

    lprn->sum = NULL;
    lprn->ssq = NULL;
    lprn->xbak = NULL;
    lprn->diff = NULL;
    lprn->na = NULL;

    return E_ALLOC;
}

static void loop_print_init (LOOP_PRINT *lprn, int lno)
{
    lprn->lineno = lno;
    lprn->nvars = 0;
    lprn->names = NULL;
    lprn->sum = NULL;
    lprn->ssq = NULL;
    lprn->xbak = NULL;
    lprn->diff = NULL;
    lprn->na = NULL;
}

static LOOP_PRINT *get_loop_print_by_line (LOOPSET *loop, int lno, int *err)
{
    LOOP_PRINT *prns;
    int i, np = loop->n_prints;

    for (i=0; i<np; i++) {
        if (loop->prns[i].lineno == lno) {
            return &loop->prns[i];
        }
    }

    prns = realloc(loop->prns, (np + 1) * sizeof *prns);
    if (prns == NULL) {
        *err = E_ALLOC;
        return NULL;
    } else {
        loop->prns = prns;
    }

    loop_print_init(&loop->prns[np], lno);
    loop->n_prints += 1;

    return &loop->prns[np];
}

static void loop_store_free (LOOP_STORE *lstore)
{
    destroy_dataset(lstore->dset);
    lstore->dset = NULL;

    strings_array_free(lstore->names, lstore->nvars);
    lstore->nvars = 0;
    lstore->names = NULL;

    free(lstore->fname);
    lstore->fname = NULL;

    lstore->lineno = -1;
    lstore->n = 0;
    lstore->opt = OPT_NONE;
}

static int loop_store_set_filename (LOOP_STORE *lstore,
                                    const char *fname,
                                    gretlopt opt)
{
    if (fname == NULL || *fname == '\0') {
        return E_ARGS;
    }

    lstore->fname = gretl_strdup(fname);
    if (lstore->fname == NULL) {
        return E_ALLOC;
    }

    lstore->opt = opt;

    return 0;
}

static void loop_store_init (LOOP_STORE *lstore)
{
    lstore->lineno = -1;
    lstore->n = 0;
    lstore->nvars = 0;
    lstore->names = NULL;
    lstore->fname = NULL;
    lstore->opt = OPT_NONE;
    lstore->dset = NULL;
}

/* check, allocate and initialize loop data storage */

static int loop_store_start (LOOPSET *loop, const char *names,
                             const char *fname, gretlopt opt)
{
    LOOP_STORE *lstore = &loop->store;
    int i, n, err = 0;

    if (names == NULL || *names == '\0') {
        gretl_errmsg_set("'store' list is empty");
        return E_DATA;
    }

    lstore->names = gretl_string_split(names, &lstore->nvars, NULL);
    if (lstore->names == NULL) {
        return E_ALLOC;
    }

    err = loop_store_set_filename(lstore, fname, opt);
    if (err) {
        return err;
    }

    n = (loop->itermax > 0)? loop->itermax : DEFAULT_NOBS;

    lstore->dset = create_auxiliary_dataset(lstore->nvars + 1, n, 0);
    if (lstore->dset == NULL) {
        return E_ALLOC;
    }

#if LOOP_DEBUG > 1
    fprintf(stderr, "loop_store_init: created sZ, v = %d, n = %d\n",
            lstore->dset->v, lstore->dset->n);
#endif

    for (i=0; i<lstore->nvars && !err; i++) {
        const char *s = lstore->names[i];

        if (!gretl_is_scalar(s)) {
            gretl_errmsg_sprintf(_("'%s': not a scalar"), s);
            err = E_DATA;
        } else {
            strcpy(lstore->dset->varname[i+1], s);
        }
    }

    return err;
}

static int loop_store_update (LOOPSET *loop, int j,
                              const char *names,
                              const char *fname,
                              gretlopt opt)
{
    LOOP_STORE *lstore = &loop->store;
    int i, t, err = 0;

    if (lstore->lineno >= 0 && lstore->lineno != j) {
        gretl_errmsg_set("Only one 'store' command is allowed in a "
                         "progressive loop");
        return E_DATA;
    }

    if (lstore->dset == NULL) {
        /* not started yet */
        err = loop_store_start(loop, names, fname, opt);
        if (err) {
            return err;
        }
        lstore->lineno = j;
        loop->cmds[j].flags |= LOOP_CMD_PDONE;
    }

    t = lstore->n;

    if (t >= lstore->dset->n) {
        if (extend_loop_dataset(lstore)) {
            err = E_ALLOC;
        }
    }

    for (i=0; i<lstore->nvars && !err; i++) {
        lstore->dset->Z[i+1][t] =
            gretl_scalar_get_value(lstore->names[i], &err);
    }

    if (!err) {
        lstore->n += 1;
    }

    return err;
}

/* See if we already have a LOOP_MODEL in place for the command
   on line @lno of the loop.  If so, return it, else create
   a new LOOP_MODEL and return it.
*/

static LOOP_MODEL *
get_loop_model_by_line (LOOPSET *loop, int lno, int *err)
{
    LOOP_MODEL *lmods;
    int n = loop->n_loop_models;
    int i;

#if LOOP_DEBUG > 1
    fprintf(stderr, "get_loop_model_by_line: loop->n_loop_models = %d\n",
            loop->n_loop_models);
#endif

    for (i=0; i<n; i++) {
        if (loop->lmodels[i].lineno == lno) {
            return &loop->lmodels[i];
        }
    }

    lmods = realloc(loop->lmodels, (n + 1) * sizeof *loop->lmodels);
    if (lmods == NULL) {
        *err = E_ALLOC;
        return NULL;
    }

    loop->lmodels = lmods;
    loop_model_init(&loop->lmodels[n], lno);
    loop->n_loop_models += 1;

    return &loop->lmodels[n];
}

#define realdiff(x,y) (fabs((x)-(y)) > 2.0e-13)

/* Update the info stored in LOOP_MODEL based on the results in pmod.
   If this is the first use we have to do some allocation first.
*/

static int loop_model_update (LOOP_MODEL *lmod, MODEL *pmod)
{
    mpf_t m;
    int j, err = 0;

#if LOOP_DEBUG > 1
    fprintf(stderr, "loop_model_update: lmod = %p, pmod = %p\n",
            (void *) lmod, (void *) pmod);
#endif

    if (lmod == NULL) {
        fprintf(stderr, "loop_model_update: got NULL loop model\n");
        return E_DATA;
    }

    if (lmod->nc == 0) {
        /* not started yet */
        err = loop_model_start(lmod, pmod);
        if (err) {
            return err;
        }
    } else if (pmod->ncoeff != lmod->nc) {
        gretl_errmsg_set(_("progressive loop: model must be of constant size"));
        return E_DATA;
    }

    mpf_init(m);

    for (j=0; j<pmod->ncoeff; j++) {
        mpf_set_d(m, pmod->coeff[j]);
        mpf_add(lmod->sum_coeff[j], lmod->sum_coeff[j], m);
        mpf_mul(m, m, m);
        mpf_add(lmod->ssq_coeff[j], lmod->ssq_coeff[j], m);

        mpf_set_d(m, pmod->sderr[j]);
        mpf_add(lmod->sum_sderr[j], lmod->sum_sderr[j], m);
        mpf_mul(m, m, m);
        mpf_add(lmod->ssq_sderr[j], lmod->ssq_sderr[j], m);
        if (!na(lmod->cbak[j]) && realdiff(pmod->coeff[j], lmod->cbak[j])) {
            lmod->cdiff[j] = 1;
        }
        if (!na(lmod->sbak[j]) && realdiff(pmod->sderr[j], lmod->sbak[j])) {
            lmod->sdiff[j] = 1;
        }
        lmod->cbak[j] = pmod->coeff[j];
        lmod->sbak[j] = pmod->sderr[j];
    }

    mpf_clear(m);

    lmod->n += 1;

#if LOOP_DEBUG > 1
    fprintf(stderr, "loop_model_update: returning %d\n", err);
#endif

    return err;
}

/* Update the LOOP_PRINT struct @lprn using the current values of the
   specified variables. If this is the first use we need to do some
   allocation first.
*/

static int loop_print_update (LOOPSET *loop, int j, const char *names)
{
    LOOP_PRINT *lprn;
    int err = 0;

    lprn = get_loop_print_by_line(loop, j, &err);

    if (!err && lprn->names == NULL) {
        /* not started yet */
        err = loop_print_start(lprn, names);
        if (!err) {
            loop->cmds[j].flags |= LOOP_CMD_PDONE;
        }
    }

    if (!err) {
        mpf_t m;
        double x;
        int i;

        mpf_init(m);

        for (i=0; i<lprn->nvars; i++) {
            if (lprn->na[i]) {
                continue;
            }
            x = gretl_scalar_get_value(lprn->names[i], &err);
            if (err) {
                break;
            }
            if (na(x)) {
                lprn->na[i] = 1;
                continue;
            }
            mpf_set_d(m, x);
            mpf_add(lprn->sum[i], lprn->sum[i], m);
            mpf_mul(m, m, m);
            mpf_add(lprn->ssq[i], lprn->ssq[i], m);
            if (!na(lprn->xbak[i]) && realdiff(x, lprn->xbak[i])) {
                lprn->diff[i] = 1;
            }
            lprn->xbak[i] = x;
        }

        mpf_clear(m);

        lprn->n += 1;
    }

    return err;
}

static void print_loop_coeff (const DATASET *dset,
                              const LOOP_MODEL *lmod,
                              int i, PRN *prn)
{
    char pname[VNAMELEN];
    char tmp[NAMETRUNC];
    mpf_t c1, c2, m, sd1, sd2;
    unsigned long ln = lmod->n;

    mpf_init(c1);
    mpf_init(c2);
    mpf_init(m);
    mpf_init(sd1);
    mpf_init(sd2);

    mpf_div_ui(c1, lmod->sum_coeff[i], ln);
    if (lmod->cdiff[i] == 0) {
        mpf_set_d(sd1, 0.0);
    } else {
        mpf_mul(m, c1, c1);
        mpf_mul_ui(m, m, ln);
        mpf_sub(m, lmod->ssq_coeff[i], m);
        mpf_div_ui(sd1, m, ln);
        if (mpf_cmp_d(sd1, 0.0) > 0) {
            mpf_sqrt(sd1, sd1);
        } else {
            mpf_set_d(sd1, 0.0);
        }
    }

    mpf_div_ui(c2, lmod->sum_sderr[i], ln);
    if (lmod->sdiff[i] == 0) {
        mpf_set_d(sd2, 0.0);
    } else {
        mpf_mul(m, c2, c2);
        mpf_mul_ui(m, m, ln);
        mpf_sub(m, lmod->ssq_sderr[i], m);
        mpf_div_ui(sd2, m, ln);
        if (mpf_cmp_d(sd2, 0.0) > 0) {
            mpf_sqrt(sd2, sd2);
        } else {
            mpf_set_d(sd2, 0.0);
        }
    }

    gretl_model_get_param_name(lmod->model0, dset, i, pname);
    maybe_trim_varname(tmp, pname);
    pprintf(prn, "%*s", 15, tmp); /* FIXME length */
    pprintf(prn, "%#14g %#14g %#14g %#14g\n", mpf_get_d(c1), mpf_get_d(sd1),
            mpf_get_d(c2), mpf_get_d(sd2));

    mpf_clear(c1);
    mpf_clear(c2);
    mpf_clear(m);
    mpf_clear(sd1);
    mpf_clear(sd2);
}

static void loop_model_print (LOOP_MODEL *lmod, const DATASET *dset,
                              PRN *prn)
{
    char startdate[OBSLEN], enddate[OBSLEN];
    int i;

    ntolabel(startdate, lmod->model0->t1, dset);
    ntolabel(enddate, lmod->model0->t2, dset);

    pputc(prn, '\n');
    pprintf(prn, _("%s estimates using the %d observations %s-%s\n"),
            _(estimator_string(lmod->model0, prn)), lmod->model0->nobs,
            startdate, enddate);
    print_model_vcv_info(lmod->model0, dset, prn);
    pprintf(prn, _("Statistics for %d repetitions\n"), lmod->n);
    pprintf(prn, _("Dependent variable: %s\n\n"),
            gretl_model_get_depvar_name(lmod->model0, dset));

    pputs(prn, _("                     mean of      std. dev. of     mean of"
                 "     std. dev. of\n"
                 "                    estimated      estimated"
                 "      estimated      estimated\n"
                 "      Variable     coefficients   coefficients   std. errors"
                 "    std. errors\n\n"));

    for (i=0; i<lmod->model0->ncoeff; i++) {
        print_loop_coeff(dset, lmod, i, prn);
    }

    pputc(prn, '\n');
}

static void loop_print_print (LOOP_PRINT *lprn, PRN *prn)
{
    bigval mean, m, sd;
    int len, maxlen = 7;
    int i, n;
    const char *s;

    if (lprn == NULL) {
        return;
    }

    n = lprn->n;

    mpf_init(mean);
    mpf_init(m);
    mpf_init(sd);

    for (i=0; i<lprn->nvars; i++) {
        len = strlen(lprn->names[i]);
        if (len > maxlen) {
            maxlen = len;
        }
    }

    pprintf(prn, _("Statistics for %d repetitions\n"), n);
    pputc(prn, '\n');
    bufspace(maxlen + 1, prn);

    len = get_utf_width(_("mean"), 14);
    pprintf(prn, "%*s ", len, _("mean"));

    len = get_utf_width(_("std. dev"), 14);
    pprintf(prn, "%*s\n", len, _("std. dev"));

    for (i=0; i<lprn->nvars; i++) {
        s = lprn->names[i];
        if (lprn->na[i]) {
            pprintf(prn, "%*s", maxlen + 1, s);
            pprintf(prn, "%14s %14s\n", "NA   ", "NA   ");
            continue;
        }
        mpf_div_ui(mean, lprn->sum[i], (unsigned long) n);
        if (lprn->diff[i] == 0) {
            mpf_set_d(sd, 0.0);
        } else {
            mpf_mul(m, mean, mean);
            mpf_mul_ui(m, m, (unsigned long) n);
            mpf_sub(sd, lprn->ssq[i], m);
            mpf_div_ui(sd, sd, (unsigned long) n);
            if (mpf_cmp_d(sd, 0.0) > 0) {
                mpf_sqrt(sd, sd);
            } else {
                mpf_set_d(sd, 0.0);
            }
        }
        pprintf(prn, "%*s", maxlen + 1, s);
        pprintf(prn, "%#14g %#14g\n", mpf_get_d(mean), mpf_get_d(sd));
    }

    mpf_clear(mean);
    mpf_clear(m);
    mpf_clear(sd);

    pputc(prn, '\n');
}

static int loop_store_save (LOOP_STORE *lstore, PRN *prn)
{
    int *list;
    int err = 0;

    list = gretl_consecutive_list_new(1, lstore->dset->v - 1);
    if (list == NULL) {
        return E_ALLOC;
    }

    lstore->dset->t2 = lstore->n - 1;
    pprintf(prn, _("store: using filename %s\n"), lstore->fname);
    err = write_data(lstore->fname, list, lstore->dset, lstore->opt, prn);

    if (err) {
        pprintf(prn, _("write of data file failed\n"));
    }

    free(list);

    return err;
}

static int extend_loop_dataset (LOOP_STORE *lstore)
{
    double *x;
    int oldn = lstore->dset->n;
    int n = oldn + DEFAULT_NOBS;
    int i, t;

    for (i=0; i<lstore->dset->v; i++) {
        x = realloc(lstore->dset->Z[i], n * sizeof *x);
        if (x == NULL) {
            return E_ALLOC;
        }
        lstore->dset->Z[i] = x;
        for (t=oldn; t<n; t++) {
            lstore->dset->Z[i][t] = (i == 0)? 1.0 : NADBL;
        }
    }

    lstore->dset->n = n;
    lstore->dset->t2 = n - 1;

    ntolabel(lstore->dset->endobs, n - 1, lstore->dset);

    return 0;
}

static void progressive_loop_zero (LOOPSET *loop)
{
    int i;

    /* What we're doing here is debatable: could we get
       away with just "zeroing" the relevant structures
       in an appropriate way, rather than destroying
       them? Maybe, but so long as we're destroying them
       we have to remove the "started" flags from
       associated "print" and "store" commands, or else
       things will go awry on the second execution of
       a nested progressive loop.
    */

    if (loop->cmds != NULL) {
        for (i=0; i<loop->n_cmds; i++) {
            if (loop->cmds[i].ci == PRINT ||
                loop->cmds[i].ci == STORE) {
                /* reset */
                loop->cmds[i].flags &= ~LOOP_CMD_PDONE;
            }
        }
    }

    for (i=0; i<loop->n_loop_models; i++) {
        loop_model_free(&loop->lmodels[i]);
    }

    loop->lmodels = NULL;
    loop->n_loop_models = 0;

    for (i=0; i<loop->n_prints; i++) {
        loop_print_free(&loop->prns[i]);
    }

    loop->prns = NULL;
    loop->n_prints = 0;

    loop_store_free(&loop->store);
}

static int model_command_post_process (ExecState *s,
                                       DATASET *dset,
                                       LOOPSET *loop,
                                       int j)
{
    int prog = loop_is_progressive(loop);
    int moderr = check_gretl_errno();
    int err = 0;

    if (moderr) {
        if (prog || model_print_deferred(s->cmd->opt)) {
            err = moderr;
        } else {
            errmsg(moderr, s->prn);
        }
    } else if (prog && !(s->cmd->opt & OPT_Q)) {
        LOOP_MODEL *lmod = get_loop_model_by_line(loop, j, &err);

        if (!err) {
            err = loop_model_update(lmod, s->model);
            set_as_last_model(s->model, GRETL_OBJ_EQN);
        }
    } else if (model_print_deferred(s->cmd->opt)) {
        MODEL *pmod = get_model_record_by_line(loop, j, &err);

        if (!err) {
            swap_models(s->model, pmod);
            pmod->ID = j + 1;
            set_as_last_model(pmod, GRETL_OBJ_EQN);
            model_count_minus(NULL);
        }
    } else {
        loop_print_save_model(s->model, dset, s->prn, s);
    }

    return err;
}

#define loop_literal(lc) (lc->flags & LOOP_CMD_LIT)

static void progressive_loop_finalize (LOOPSET *loop,
                                       const DATASET *dset,
                                       PRN *prn)
{
    int i, j = 0, k = 0;

    for (i=0; i<loop->n_cmds; i++) {
	loop_command *lc = &loop->cmds[i];

        if (plain_model_ci(lc->ci) && !(lc->opt & OPT_Q)) {
            loop_model_print(&loop->lmodels[j], dset, prn);
            loop_model_zero(&loop->lmodels[j], 1);
            j++;
        } else if (lc->ci == PRINT && !loop_literal(lc)) {
            loop_print_print(&loop->prns[k], prn);
            loop_print_zero(&loop->prns[k], 1);
            k++;
        } else if (lc->ci == STORE) {
            loop_store_save(&loop->store, prn);
        }
    }
}

#define prog_cmd_started(l,j) (l->cmds[j].flags & LOOP_CMD_PDONE)

#define not_ok_in_progloop(c) (NEEDS_MODEL_CHECK(c) || \
                               c == NLS ||  \
                               c == MLE ||  \
                               c == GMM)

static int handle_prog_command (LOOPSET *loop, int j,
                                CMD *cmd, int *err)
{
    loop_command *lc = &loop->cmds[j];
    int handled = 0;

    if (cmd->ci == PRINT && !loop_literal(lc)) {
        if (prog_cmd_started(loop, j)) {
            *err = loop_print_update(loop, j, NULL);
        } else {
            *err = loop_print_update(loop, j, cmd->parm2);
        }
        handled = 1;
    } else if (cmd->ci == STORE) {
        if (prog_cmd_started(loop, j)) {
            *err = loop_store_update(loop, j, NULL, NULL, 0);
        } else {
            *err = loop_store_update(loop, j, cmd->parm2, cmd->param,
                                     cmd->opt);
        }
        handled = 1;
    } else if (not_ok_in_progloop(cmd->ci)) {
        gretl_errmsg_sprintf(_("%s: not implemented in 'progressive' loops"),
                             gretl_command_word(cmd->ci));
        *err = 1;
        handled = 1;
    }

    return handled;
}
