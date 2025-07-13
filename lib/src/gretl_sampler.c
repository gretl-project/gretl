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

static gretl_matrix *do_sampler (char **init, int ni,
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

#if 0

/* the code below is not ready to compile yet */

static gretl_matrix *extract_param (const char *s, int i,
                                    const gretl_bundle *allInfo)
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
            gretl_bundle *b = gretl_array_get_data(dists, i-1);
            gretl_array *a = gretl_bundle_get_array(a, "private", &err);
            
            sscanf(tmp+1, "%d]", &k);
            ret = a[k];
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
        ret = feval(s+2, allInfo);
    }

    if (ret == NULL) {
        ret = gretl_matrix_alloc(1,1);
        ret->val[0] = NADBL;
    }

    g_free(tmp);

    return ret;
}

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
        gretl_matrix *mval;
        gretl_matrix *m;
        gretl_matrix *S;
        int dim;

        dim = extract_param(map[0], which, B);
        m = extract_param(map[1], which, B);
        S = extract_param(map[2], which, B);
        if (!strcmp(id, "N")) {
            mval = vec(m) + cholesky(S) * mnormal(dim,1); /* transp */
        } else {
            mval = vec(S*m) + cholesky(S) * mnormal(dim,1); /* transp */
        }
        gretl_bundle_donate_data(this, "value", mval, GRETL_TYPE_MATRIX, 0);    
    } else if (!strcmp(id, "U")) {
        double a = extract_param(map[0], which, B);
        double b = extract_param(map[1], which, B);
        
        this.value = randgen1(u, a, b);
    } else if (!strcmp(id, "B")) {
        double p = extract_param(map[0], which, B);
        double d = randgen1(u, 0, 1) < p;

        gretl_bundle_set_scalar(this, "value", d);
    } else if (!strcmp(id, "G") || !strcmp(id, "IG")) {
        double p = extract_param(map[0], which, B);
        double a = extract_param(map[1], which, B);
        double d = randgen1(g, p, 1/a);

        if (!strcmp(dist, "G")) {
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
    ret = calloc(n, sizeof *ret);

    for (i=0; i<n && !err; i++) {
        this = gretl_array_get_data(dists, i);
        id = gretl_bundle_get_string(this, "dist", &err);
        if (!strcmp(id, "N") || !strcmp(id, "N2")) {
            gretl_array *map;
            gretl_matrix *m;
            int dim;

            map = gretl_bundle_get_array(this, "map", &err);
            dim = extract_param(gretl_array_get_data(map, 0), i, allInfo);
            m = gretl_random_matrix_new(dim, 1, D_NORMAL);
            gretl_bundle_donate_data(this, "value", m, GRETL_TYPE_MATRIX, 0);
            ret[i] = dim;
        } else if (!strcmp(id, "U")) {
            gretl_bundle_set_scalar(this, "value", 0.5);
            ret[i] = 1;
        } else if (!strcmp(id, "B")) {
            double ru = randgen1(u, 0, 1);
            
            gretl_bundle_set_scalar(this, "value", ru > 0.5);
            ret[i] = 1;
        } else if (!strcmp(id, "G") || !strcmp(id, "IG")) {
            gretl_bundle_set_scalar(this, "value", 1);
            ret[i] = 1;
        }
    }

    return ret;
}

gretl_matrix *gibbs (const gretl_bundle *B, int T)
{
    gretl_matrix *ret;
    gretl_array *dists;
    gretl_bundle *localB;
    int *nparams;
    int i, ini, fin;
    int n, err = 0;

    /* the array of distributions to sample */
    dists = gretl_bundle_get_array(B, "dists", &err);
    if (err) {
        return NULL;
    }

    /* the number of distributions */
    n = gretl_array_get_length(dists);

    localB = gretl_bundle_copy(B, &err);
    if (!err) {
        nparams = gibbs_init(localB, &err);
    }
    
    if (!err) {
        int m = 0;
        
        for (i=0; i<n; i++) {
            m += nparams[i];
        }
        ret = gretl_zero_matrix_new(T, m);
    }

    for (t=0; t<T && !err; t++) {
        fin = 0;
        for (i=0; i<n && !err; i++) {
            err = gen_one(localB, i);
            ini = fin + 1;
            fin += nparams[i];
        }
    }

    free(nparams);
    gretl_bundle_destroy(localB);
            
    return ret;
}

#endif /* not ready */
            

