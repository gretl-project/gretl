
#define MAP_DEBUG 1

/* Given the current dataset and the $mapfile name recorded
   on it: get the content of $mapfile as a bundle then
   revise the bundle (a) to include only the features in the
   current sample and (b) to reflect any changes in the dataset
   (series added, deleted or modified). Return the modified
   bundle.
*/

gretl_bundle *get_current_map (const DATASET *dset, int *err)
{
    const char *sj, *id, *fname;
    gretl_bundle *fi, *pp, *jb = NULL;
    gretl_array *features = NULL;
    int n, fmax = 0;
    int i, j, dsi, fidx;
    double xj;

    fname = dataset_get_mapfile(dset);

    if (fname == NULL) {
	gretl_errmsg_set("no mapfile is present");
	*err = E_DATA;
    } else if (dataset_is_resampled(dset)) {
	/* most unlikely */
	gretl_errmsg_set("dataset is resampled!");
	*err = E_DATA;
    }

    if (!*err) {
	jb = gretl_bundle_read_from_file(fname, 0, err);
	if (!*err) {
	    features = gretl_bundle_get_array(jb, "features", err);
	}
    }

    if (!*err) {
	fmax = gretl_array_get_length(features);
	if (dset->submask != NULL) {
	    n = get_full_length_n();
	} else {
	    n = dset->n;
	}
	if (fmax != n) {
	    /* Although it may be sub-sampled the full dataset must
	       have a number of observations equal to the number of
	       features in the existing map.
	    */
	    gretl_errmsg_set("map and dataset are out of sync!");
	    *err = E_DATA;
	}
	n = fmax;
    }

    if (*err) {
	gretl_bundle_destroy(jb);
	return NULL;
    }

    dsi = -1;

    for (i=0, fidx=0; i<fmax; i++) {
	int skip = 0;
#if MAP_DEBUG
	fprintf(stderr, "get_current_map: i=%d\n", i);
#endif
	if (dset->submask != NULL) {
	    if (dset->submask[i] == 0) {
#if MAP_DEBUG
		fprintf(stderr, "  skip masked feature %d\n", i);
#endif
		skip = 1;
	    } else {
		dsi++;
	    }
	} else {
	    dsi = i;
	}
	if (dsi < dset->t1) {
	    skip = 1;
	} else if (dsi > dset->t2) {
	    /* we've got everything we need */
	    break;
	}
	if (skip) {
#if MAP_DEBUG
	    fprintf(stderr, "  skip: delete feature %d\n", i);
#endif
	    gretl_array_delete_element(features, fidx);
	    fmax--;
	} else {
	    fi = gretl_array_get_element(features, fidx, NULL, err);
	    id = gretl_bundle_get_string(fi, "id", err);
#if MAP_DEBUG
	    fprintf(stderr, "  include feature %d (%s)\n", i, id);
#endif
	    pp = gretl_bundle_get_bundle(fi, "properties", err);
	    /* clear the existing properties bundle */
	    gretl_bundle_void_content(pp);
	    /* and refill it from the dataset */
	    for (j=1; j<dset->v; j++) {
		id = dset->varname[j];
		if (is_string_valued(dset, j)) {
		    sj = series_get_string_for_obs(dset, j, dsi);
		    gretl_bundle_set_string(pp, id, sj);
		} else {
		    xj = dset->Z[j][dsi];
		    gretl_bundle_set_scalar(pp, id, xj);
		}
	    }
	    fidx++;
	}
    }

    /* delete any unwanted trailing features */
    for (j=fidx; j<n; j++) {
	gretl_array_delete_element(features, fidx);
    }

    return jb;
}
