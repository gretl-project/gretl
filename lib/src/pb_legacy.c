
/* Here: functions supporting the legacy plot band syntax that was in
   force prior to the gretl 2023c release.
*/

/* Handle the case where we get to the band-plot code from a command
   in which the data to be plotted (and hence also the band
   specification) are given in matrix form. By this point the
   plot-data have been converted to (temporary) DATASET form; here we
   retrieve the band-spec matrix, check it for conformability, and
   stick the two extra columns onto the dataset (borrowing pointers
   into the matrix content).
*/

static int parse_band_matrix_option (band_info *bi,
                                     const char *spec,
                                     gnuplot_info *gi,
                                     DATASET *dset)
{
    const char *mname = NULL;
    gchar **S;
    int i = 0;
    int err = 0;

    S = g_strsplit(spec, ",", -1);

    while (S != NULL && S[i] != NULL && !err) {
        if (i == 0) {
            mname = S[i];
        } else if (i == 1) {
            /* spec for width multiplier: optional */
            if (numeric_string(S[i])) {
                bi->factor = dot_atof(S[i]);
            } else if (gretl_is_scalar(S[i])) {
                bi->factor = gretl_scalar_get_value(S[i], &err);
            } else {
                err = invalid_field_error(S[i]);
            }
        } else {
            /* we got too many comma-separated terms */
            err = invalid_field_error(S[i]);
        }
        i++;
    }

    if (!err && (bi->factor <= 0 || na(bi->factor))) {
        err = E_INVARG;
    }
    if (!err) {
        gretl_matrix *m = get_matrix_by_name(mname);

        if (m == NULL) {
            err = E_INVARG;
        } else {
            err = process_band_matrix(m, bi, gi, dset);
        }
    }

    g_strfreev(S);

    return err;
}

static int handle_recession_bars (band_info *bi,
                                  const char *s,
                                  const DATASET *dset)
{
    int v = current_series_index(dset, s);
    int err = 0;

    if (v >= 0 && v < dset->v) {
        err = check_assign_bdummy(bi, v, dset);
    } else {
        err = E_INVARG;
    }

    return err;
}

/* Handle the band plus-minus option for all cases apart from the one
   handled just above. Here we require two comma-separated series
   identifiers for center and width.
*/

static int parse_band_pm_option (band_info *bi,
                                 const char *spec,
                                 gnuplot_info *gi,
                                 const DATASET *dset,
                                 gretlopt opt)
{
    gchar **S;
    int v, i = 0;
    int err = 0;

    if (strchr(spec, ',') == NULL) {
        return handle_recession_bars(bi, spec, dset);
    }

    S = g_strsplit(spec, ",", -1);

    while (S != NULL && S[i] != NULL && !err) {
        if (i < 2) {
            /* specs for the "center" and "width" series: required */
            if (integer_string(S[i])) {
                /* var ID number? */
                v = atoi(S[i]);
            } else {
                /* varname? */
                v = current_series_index(dset, S[i]);
            }
            if (v >= 0 && v < dset->v) {
                do_center_or_width(bi, gi, i, v);
            } else {
                err = invalid_field_error(S[i]);
            }
        } else if (i == 2) {
            /* spec for width multiplier: optional */
            if (numeric_string(S[i])) {
                bi->factor = dot_atof(S[i]);
            } else if (gretl_is_scalar(S[i])) {
                bi->factor = gretl_scalar_get_value(S[i], &err);
            } else {
                err = invalid_field_error(S[i]);
            }
        } else {
            /* we got too many comma-separated terms */
            err = invalid_field_error(S[i]);
        }
        i++;
    }

    g_strfreev(S);

#if PB_DEBUG
    fprintf(stderr, "parse_band_pm_option: err = %d\n", err);
    fprintf(stderr, "center = %d\n", bi->center);
    fprintf(stderr, "width = %d\n", bi->width);
    fprintf(stderr, "factor = %g\n", bi->factor);
#endif

    if (!err) {
        if (bi->center < 0 || bi->width < 0 ||
            bi->factor < 0 || na(bi->factor)) {
            err = E_INVARG;
        }
    }

    return err;
}

/* We're looking here for any one of three patterns:

   <style>
   <style>,<color>
   <color>

   where <style> should be "fill", "dash" or "line" (the default) and
   <color> should be a hex string such as "#00ff00" or "0x00ff00".
*/

static int parse_band_style_option (band_info *bi)
{
    const char *s = get_optval_string(GNUPLOT, OPT_J);
    int err = 0;

    if (s != NULL) {
        const char *p = strchr(s, ',');
        int do_color = 1;
        guint32 u = 0;

        if (bi->bdummy && *s != ',') {
            /* must be just a color */
            u = numeric_color_from_string(s, &err);
        } else if (*s == ',') {
            /* skipping field 1, going straight to color */
            u = numeric_color_from_string(s + 1, &err);
        } else if (p == NULL) {
            /* just got field 1, style spec */
            bi->style = style_from_string(s, &err);
            do_color = 0;
        } else {
            /* embedded comma: style + color */
            if (strlen(s) >= 8 && s[4] == ',') {
                gchar *tmp = g_strndup(s, 4);

                bi->style = style_from_string(tmp, &err);
                g_free(tmp);
            } else {
                err = invalid_field_error(s);
            }
            if (!err) {
                u = numeric_color_from_string(s + 5, &err);
            }
            if (!err && do_color) {
                sprintf(bi->rgb, "#%x", u);
            }
        }
    }

    return err;
}

static band_info **legacy_get_band_info (const char *spec,
                                         int matrix_mode,
                                         gnuplot_info *gi,
                                         DATASET *dset,
                                         gretlopt opt,
                                         int *err)
{
    band_info **bi = NULL;

    bi = malloc(sizeof *bi);
    bi[0] = band_info_new();
    if (bi[0] == NULL) {
        *err = E_ALLOC;
        return NULL;
    }

    if (matrix_mode) {
        /* In this case the band should be given in the form of a
           named matrix with two columns holding center and width,
           respectively.
        */
        *err = parse_band_matrix_option(bi[0], spec, gi, dset);
    } else {
        *err = parse_band_pm_option(bi[0], spec, gi, dset, opt);
    }
    if (!*err && (opt & OPT_J)) {
        *err = parse_band_style_option(bi[0]);
    }

    if (*err) {
        free_bbi(bi, 1);
        bi = NULL;
    }

    return bi;
}
