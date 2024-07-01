#include "gretl_func.h"
#include "gretl_string_table.h" /* for csvdata */
#include "csvdata.h"
#include "gretl_join.h"

/* call hansl code from dbnomics.gfn to get a series bundle */

static gretl_bundle *get_dbn_series_bundle (const char *datacode,
					    int *err)
{
    gretl_bundle *b = NULL;
    fncall *fc;

    fc = get_pkg_function_call("dbnomics_get_series", "dbnomics", NULL);
    if (fc == NULL) {
	*err = E_DATA;
    } else {
	*err = push_anon_function_arg(fc, GRETL_TYPE_STRING,
				      (void *) datacode);
	if (!*err) {
	    *err = gretl_function_exec(fc, GRETL_TYPE_BUNDLE, NULL,
				       &b, NULL);
	}
	if (b != NULL) {
	    int dberr = gretl_bundle_get_int(b, "error", NULL);

	    if (dberr) {
		const char *msg =
		    gretl_bundle_get_string(b, "errmsg", NULL);

		*err = E_DATA;
		if (msg != NULL) {
		    gretl_errmsg_set(msg);
		} else {
		    gretl_errmsg_sprintf(_("%s: no data found"), datacode);
		}
		gretl_bundle_destroy(b);
		b = NULL;
	    }
	} else if (!*err) {
	    gretl_errmsg_sprintf(_("%s: no data found"), datacode);
	    *err = E_DATA;
	}
    }

    return b;
}

static int write_dbn_csv (gretl_bundle *b,
			  const char *srcname,
			  const char *fname,
			  int *pd, char **stobs)
{
    FILE *fp = NULL;
    gretl_matrix *v = NULL;
    gretl_array *P = NULL;
    const char *Pt;
    int t, T;
    int err = 0;

    T = gretl_bundle_get_int(b, "T", &err);
    P = gretl_bundle_get_array(b, "period", &err);
    v = gretl_bundle_get_matrix(b, "value", &err);
    *pd = gretl_bundle_get_int(b, "frequency", &err);

    if (!err && (T <= 0 || *pd <= 0 || P == NULL || v == NULL)) {
	return E_DATA;
    }

    *stobs = gretl_array_get_data(P, 0);

    fp = fopen(fname, "wb");
    if (fp == NULL) {
	return E_FOPEN;
    }

    gretl_push_c_numeric_locale();

    fprintf(fp, "obs %s\n", srcname);
    for (t=0; t<T; t++) {
	Pt = gretl_array_get_data(P, t);
	if (t >= v->rows || na(v->val[t])) {
	    fprintf(fp, "%s NA\n", Pt);
	} else {
	    fprintf(fp, "%s %.12g\n", Pt, v->val[t]);
	}
    }

    gretl_pop_c_numeric_locale();
    fclose(fp);

    return 0;
}

/* convert from dbread's CompactMethod spec to join's
   AggrType
*/

static int dbn_get_aggr (CompactMethod cmethod,
			 int *seqval, int *err)
{
    AggrType aggr = AGGR_NONE;

    if (cmethod == COMPACT_SUM) {
	aggr = AGGR_SUM;
    } else if (cmethod == COMPACT_AVG) {
	aggr = AGGR_AVG;
    } else if (cmethod == COMPACT_SPREAD) {
	aggr = AGGR_MIDAS;
    } else if (cmethod == COMPACT_SOP) {
	aggr = AGGR_SEQ;
	*seqval = 1;
    } else if (cmethod == COMPACT_EOP) {
	aggr = AGGR_SEQ;
	*seqval = -1;
    } else {
	*err = E_INVARG;
    }

    return aggr;
}

static int get_one_dbnomics_series (const char *datacode,
				    const char *altname,
				    DATASET *dset,
				    CompactMethod cmethod,
				    PRN *prn)
{
    gchar *tmpname = NULL;
    char srcname[VNAMELEN] = {0};
    char *stobs = NULL;
    AggrType aggr = AGGR_NONE;
    int pd, seqval = 0;
    gretl_bundle *b = NULL;
    int err = 0;

    if (strchr(datacode, '/') == NULL) {
	err = E_INVARG;
    } else if (dset->v > 0 && cmethod != COMPACT_UNSET) {
	aggr = dbn_get_aggr(cmethod, &seqval, &err);
    }
    if (!err) {
	b = get_dbn_series_bundle(datacode, &err);
    }

    if (!err) {
	const char *rawname;
	gchar *basenam;
	double u;

	rawname = strrchr(datacode, '/') + 1;
	gretl_normalize_varname(srcname, rawname, 0, 0);
	gretl_rand_uniform_int_minmax(&u, 0, 0, 0, 9999, 0);
	basenam = g_strdup_printf("tmp%04d.csv", (int) u);
	tmpname = gretl_make_dotpath(basenam);
	err = write_dbn_csv(b, srcname, tmpname, &pd, &stobs);
	g_free(basenam);
    }

    if (!err && dset->v == 0) {
	/* empty dataset: just import */
	err = import_csv(tmpname, dset, OPT_NONE, NULL);
	if (!err && *altname != '\0') {
            int list[2] = {1,1};
            
	    rename_series(dset, list, altname, OPT_NONE);
	}
    } else if (!err) {
	/* otherwise use join */
	gretlopt jopt = OPT_NONE;
	const char *vnames[2] = {NULL, NULL};
	const char *rhname = NULL;
	char *okey = NULL;

	if (*altname != '\0') {
	    vnames[0] = altname;
	    rhname = srcname;
	} else {
	    vnames[0] = srcname;
	}

	if (pd == 7) {
	    jopt = OPT_K;
	    okey = "obs,%Y-%m-%d";
	} else if (strstr(stobs, "-Q")) {
	    jopt = OPT_K;
	    okey = "obs,%Y-Q%q";
	}
	err = gretl_join_data(tmpname, vnames, 1,
			      dset, NULL, okey, NULL,
			      rhname, aggr, seqval,
			      NULL, NULL, NULL,
			      0, jopt, prn);
    }

    if (!err) {
	/* retrieve and attach description */
	const char *label;
	const char *vname;
	int v;

	label = gretl_bundle_get_string(b, "series_name", NULL);
	if (label != NULL) {
	    vname = *altname ? altname : srcname;
	    v = current_series_index(dset, vname);
	    if (v > 0) {
		series_set_label(dset, v, label);
	    }
	}
    }

    gretl_bundle_destroy(b);

    if (tmpname != NULL) {
	gretl_remove(tmpname);
	g_free(tmpname);
    }

    return err;
}
