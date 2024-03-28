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

#include "mapinfo.h"
#include "matrix_extra.h"

#define GEODEBUG 0

typedef struct HSV_ {
    double h; /* hue, as degrees in [0,359] */
    double s; /* saturation, as fraction in [0,1] */
    double v; /* value, as fraction in [0,1] */
} HSV;

typedef struct RGB_ {
    guint8 r; /* red channel, in [0,255] */
    guint8 g; /* green channel, in [0,255] */
    guint8 b; /* blue channel, in [0,255] */
} RGB;

static void mapinfo_init (mapinfo *mi,
                          const char *fname,
                          gretl_bundle *map,
                          gretl_bundle *opts)
{
    mi->mapfile = fname;
    mi->map = map;
    mi->opts = opts;
    mi->flags = MAP_DISPLAY;
    mi->na_action = NA_OUTLINE;
    mi->proj = PRJ0;
}

/* check for the canonical case for a discrete payload:
   minimum value = 1 and all values are consecutive
*/

static int canonical_discrete (const double *v, int nv)
{
    int i, ret = 1;

    if (v[0] != 1) {
	ret = 0;
    } else {
	for (i=1; i<nv; i++) {
	    if (v[i] != v[i-1] + 1) {
		ret = 0;
		break;
	    }
	}
    }

    return ret;
}

/* check the payload series for the coded property, for
   whether it is string-valued, and for whether it is
   "canonical" in these respects.
*/

static int inspect_payload (mapinfo *mi, const DATASET *dset, int v)
{
    int ns = 0;
    int err = 0;

    mi->zname = dset->varname[v];

    if (gretl_isdummy(dset->t1, dset->t2, dset->Z[v])) {
	mi->zvals = gretl_matrix_alloc(2, 1);
	mi->zvals->val[0] = 0;
	mi->zvals->val[1] = 1;
	mi->n_codes = 2;
	mi->flags |= MAP_DUMMY;
	return 0;
    } else {
        char **strvals = series_get_string_vals(dset, v, &ns, 1);
        int coded = series_is_coded(dset, v);

        if (coded || ns > 0) {
            int nz = gretl_vector_get_length(mi->zvec);
            const double *z = mi->zvec->val;
            gretl_matrix *vals;
            int nv = 0;

            vals = gretl_matrix_values(z, nz, OPT_S, &err);
            if (err) {
                return err;
            }
            nv = gretl_vector_get_length(vals);
            if (canonical_discrete(vals->val, nv)) {
                mi->zvals = vals;
                mi->n_codes = nv;
            } else {
                gretl_matrix_free(vals);
            }
            if (ns == mi->n_codes) {
                mi->zlabels = strvals;
            }
        }
    }

    return err;
}

/* Driver function for calling the geoplot plugin to produce
   a map. To obtain the map polygons we need EITHER the name
   of the source file (GeoJSON or Shapefile), via @fname,
   OR the map info in the form of a gretl_bundle, via @map.

   The "payload" (if any) is given as a series, via @plx,
   and if it's a named series its ID number will be in @plv.
   Plotting options (if any) are provided via @opts.
*/

int geoplot_driver (const char *fname,
                    gretl_bundle *map,
		    int plv,
                    const double *plx,
                    const DATASET *dset,
                    gretl_bundle *opts)
{
    mapinfo mi = {0};
    int (*mapfunc) (mapinfo *);
    const char *mapfile = NULL;
    int free_map = 0;
    int err = 0;

    if (fname != NULL && map != NULL) {
        gretl_errmsg_set("geoplot: cannot give both filename and map bundle");
        return E_DATA;
    }

    if (map == NULL) {
        mapfile = dataset_get_mapfile(dset);
    }

    if (fname == NULL && map == NULL) {
        fname = mapfile;
        if (fname == NULL) {
            gretl_errmsg_set("geoplot: no map was specified");
            return E_DATA;
        }
    }

    mapinfo_init(&mi, fname, map, opts);

    if (plx != NULL) {
        /* convert payload series to vector for convenience in plugin */
        mi.zvec = gretl_vector_from_series(plx, dset->t1, dset->t2);
        if (mi.zvec == NULL) {
            err = E_ALLOC;
        } else if (plv > 0) {
	    /* check the payload series for the discrete property,
	       and for whether it is string-valued
	    */
	    err = inspect_payload(&mi, dset, plv);
	}
    }

#if GEODEBUG
    fprintf(stderr, "geoplot_driver: map=%p, mapfile=%p, fname=%p\n",
            (void *) map, (void *) mapfile, (void *) fname);
#endif

    /* In the case where we got @fname, do we want to produce a
       map bundle in which the actual map data are synced with
       the dataset? Probably so if @fname is just $mapfile
       (metadata loaded as dataset), and presumably not if @fname
       is an "external" reference.
    */
    if (map == NULL && mapfile != NULL) {
        if (fname == mapfile || !strcmp(fname, mapfile)) {
#if GEODEBUG
            fprintf(stderr, "geoplot_driver: calling get_current_map()\n");
#endif
            mi.map = get_current_map(dset, NULL, &err);
            free_map = 1;
            fname = NULL;
        }
    }

    if (!err) {
	/* call the geoplot plugin, which will in turn call
	   gnuplot if all goes well
	*/
        mapfunc = get_plugin_function("geoplot");
        if (mapfunc == NULL) {
            err = E_FOPEN;
        } else {
            err = mapfunc(&mi);
        }
    }

    gretl_matrix_free(mi.zvec);
    if (free_map) {
        gretl_bundle_destroy(map);
    }

    return err;
}

/* We do this only if the gnuplot version is lower than 5.4;
   for 5.4 and higher we use the "keyentry" mechanism instead.
*/

static void print_labeled_colorbox (mapinfo *mi,
                                    char **labels,
                                    int n, FILE *fp)
{
    double zmin = mi->zrange->val[0];
    double incr = (n - 1.0) / n;
    double x = zmin + incr / 2.0;
    int i;

    fputs("set cbtics (", fp);
    for (i=0; i<n; i++) {
	if (labels != NULL) {
	    fprintf(fp, "'%s' %g", labels[i], x);
	} else {
	    fprintf(fp, "'%s=%g' %g", mi->zname, zmin + i, x);
	}
        if (i < n-1) {
            fputs(", ", fp);
            x += incr;
        }
    }
    fputs(") scale 0\n", fp);
}

static RGB hsv_to_rgb256 (HSV *in)
{
    double hh, p, q, t, ff;
    double r, g, b;
    gint32 i;
    RGB out;

    if (in->s <= 0.0) {
	out.r = out.g = out.b = in->v * 255.0;
        return out;
    }
    hh = in->h;
    if (hh >= 360.0) {
	hh = 0.0;
    }
    hh /= 60.0;
    i = (gint32) hh;
    ff = hh - i;
    p = in->v * (1.0 - in->s);
    q = in->v * (1.0 - (in->s * ff));
    t = in->v * (1.0 - (in->s * (1.0 - ff)));

    switch(i) {
    case 0:
        r = in->v;
        g = t;
        b = p;
        break;
    case 1:
        r = q;
        g = in->v;
        b = p;
        break;
    case 2:
        r = p;
        g = in->v;
        b = t;
        break;
    case 3:
        r = p;
        g = q;
        b = in->v;
        break;
    case 4:
        r = t;
        g = p;
        b = in->v;
        break;
    case 5:
    default:
        r = in->v;
        g = p;
        b = q;
        break;
    }

    out.r = r * 255.0;
    out.g = g * 255.0;
    out.b = b * 255.0;

    return out;
}

/* This variant picks RGB colors that vary by Hue but have common
   Saturation and Value attributes.
*/

static int print_discrete_colors_hsv (int n, FILE *fp)
{
    char color[9];
    RGB rgb = {0};
    HSV hsv = {0, 0.5, 0.8};
    double delta;
    int i;

    delta = 360.0 / n;
    hsv.h = 0.5 * delta;

    fprintf(fp, "set palette maxcolors %d\n", n);
    fputs("set palette defined (", fp);
    for (i=0; i<n; i++) {
        rgb = hsv_to_rgb256(&hsv);
        sprintf(color, "#%02x%02x%02x", rgb.r, rgb.g, rgb.b);
        fprintf(fp, "%d '%s'", i, color);
        fputs((i < n-1)? ", " : ")\n", fp);
        hsv.h += delta;
    }

    fprintf(fp, "set cbrange [1:%d]\n", n);

    return 0;
}

/* This variant picks RGB colors that vary in all three of the
   HSV dimensions.
*/

static int print_discrete_colors_rgb (int n, FILE *fp)
{
    gretl_matrix *H = NULL;
    char color[9];
    guint8 r, g, b;
    int i, err = 0;

    H = halton_matrix(3, n, 20, &err);
    if (H == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<3*n; i++) {
        H->val[i] = floor(255 * H->val[i]);
    }

    fprintf(fp, "set palette maxcolors %d\n", n);
    fputs("set palette defined (", fp);
    for (i=0; i<n; i++) {
        r = (guint8) gretl_matrix_get(H, 0, i);
        g = (guint8) gretl_matrix_get(H, 1, i);
        b = (guint8) gretl_matrix_get(H, 2, i);
        sprintf(color, "#%02x%02x%02x", r, g, b);
        fprintf(fp, "%d '%s'", i, color);
        fputs((i < n-1)? ", " : ")\n", fp);
    }

    fprintf(fp, "set cbrange [1:%d]\n", n);
    gretl_matrix_free(H);

    return 0;
}

static int print_discrete_multicolors (mapinfo *mi, FILE *fp)
{
    if (mi->n_codes < 13) {
	return print_discrete_colors_hsv(mi->n_codes, fp);
    } else {
	return print_discrete_colors_rgb(mi->n_codes, fp);
    }
}

static int print_discrete_autocolors (mapinfo *mi, FILE *fp)
{
    if (mi->flags & MAP_DUMMY) {
        fputs("set palette maxcolors 2\n", fp);
        fputs("set palette defined (0 '#D22532', 1 '#244999')\n", fp);
        fputs("set cbrange [0:1]\n", fp);
    } else {
        print_discrete_multicolors(mi, fp);
    }

    return 0;
}

static int discrete_array_error (const char *s)
{
    gretl_errmsg_sprintf("Invalid discrete palette argument '%s'", s);
    return E_INVARG;
}

static int print_discrete_map_palette (mapinfo *mi,
				       const char *s,
                                       double gpver,
				       FILE *fp)
{
    gchar **anames = NULL;
    int n = 0, err = 0;

    anames = g_strsplit(s, ",", 2);
    if (anames == NULL) {
	gretl_errmsg_sprintf("Invalid discrete palette argument '%s'", s);
	err = E_INVARG;
    } else {
	n = g_strv_length(anames);
	if (n < 1 || n > 2) {
	    gretl_errmsg_sprintf("Invalid discrete palette argument '%s'", s);
	    err = E_INVARG;
	}
    }

    if (!err) {
        char **colors = NULL;
        char **labels = NULL;
        gretl_array *a;
        int i, nc = 0, nl = 0;
        int colors_done = 0;

        if (mi->n_codes > 0 && !strcmp(anames[0], "auto")) {
            nc = mi->n_codes;
            print_discrete_autocolors(mi, fp);
            colors_done = 1;
        } else {
            a = get_array_by_name(anames[0]);
            colors = gretl_array_get_strings(a, &nc);
            if (colors == NULL || nc < 2) {
                err = discrete_array_error(anames[0]);
            }
        }

        if (!err && n == 2) {
            a = get_array_by_name(anames[1]);
            labels = gretl_array_get_strings(a, &nl);
            if (labels == NULL || nl != nc) {
                err = discrete_array_error(anames[1]);
            }
        }

	if (!err) {
            if (!colors_done) {
                if (mi->flags & MAP_DUMMY) {
                    fputs("set cbrange [0:1]\n", fp);
                } else {
                    fprintf(fp, "set cbrange [1:%d]\n", nc);
                }
                fprintf(fp, "set palette maxcolors %d\n", nc);
                fputs("set palette defined (", fp);
                for (i=0; i<nc; i++) {
                    fprintf(fp, "%d '%s'", i, colors[i]);
                    fputs((i < nc-1)? ", " : ")\n", fp);
                }
            }
	    if (gpver < 5.4) {
		if (labels != NULL) {
		    print_labeled_colorbox(mi, labels, nl, fp);
		} else if (mi->zlabels != NULL) {
		    print_labeled_colorbox(mi, mi->zlabels,
                                           mi->n_codes, fp);
		} else {
		    fputs("unset colorbox\n", fp);
		}
	    } else {
		if (labels != NULL) {
		    mi->zlabels = labels;
		}
		fputs("unset colorbox\n", fp);
	    }
	}
    }

    if (anames != NULL) {
	g_strfreev(anames);
    }

    return err;
}

static void simple_print_map_palette (const char *p, FILE *fp)
{
    if (p == NULL || *p == '\0') {
	return;
    } else if (!strcmp(p, "blues")) {
        fputs("set palette defined (0 '#D4E4F2', 1 'steelblue')\n", fp);
    } else if (!strcmp(p, "oranges")) {
        fputs("set palette defined (0 '#E9D9B5', 1 'dark-orange')\n", fp);
    } else if (!strcmp(p, "green-to-red")) {
	fputs("set palette defined (0 '#58996E', 1 '#E1D99A', 2 '#C0414C')\n", fp);
    } else {
	fprintf(fp, "%s\n", p);
    }
}

static void tricky_print_map_palette (const char *p,
				      const double *zlim,
				      FILE *fp)
{
    const char *colors[3][3] = {
        { "#D4E4F2", "steelblue", NULL },   /* "blues" */
        { "#E9D9B5", "dark-orange", NULL }, /* "oranges" */
        { "#58996E", "#E1D99A", "#C0414C" } /* "green-to-red" */
    };
    int i = 3;

    if (p == NULL) {
	i = 3;
    } else if (!strcmp(p, "blues")) {
	i = 0;
    } else if (!strcmp(p, "oranges")) {
	i = 1;
    } else if (!strcmp(p, "green-to-red")) {
	i = 2;
    }

    /* FIXME: maybe allow specification of NA fill color? */

    fprintf(fp, "set palette defined (%.8g 'gray', ", zlim[0] - 0.002);

    if (i == 3) {
	const char *hc[] = {
            "#000000", "#7202F3", "#A11096",
            "#C63700", "#E48300", "#FFFF00"
	};
	double step = (zlim[1] - zlim[0] + 0.001) / 5;
	double z = zlim[0] - 0.001;
	int j;

	for (j=0; j<6; j++) {
	    fprintf(fp, "%.8g '%s'", z, hc[j]);
	    fputs(j == 1 ? ", \\\n" : j < 5 ? ", " : ")\n", fp);
	    z += step;
	}
    } else {
	fprintf(fp, "%.8g '%s', ", zlim[0] - 0.001, colors[i][0]);
	if (i < 2) {
	    fprintf(fp, "%.8g '%s')\n", zlim[1], colors[i][1]);
	} else {
	    fprintf(fp, "%.8g '%s', %.8g '%s')\n", (zlim[1] - zlim[0]) / 2,
		    colors[i][1], zlim[1], colors[i][2]);
	}
    }

    /* for this to work, cbrange has to be set using zlim */
    fprintf(fp, "set cbrange [%.8g:%.8g]\n", zlim[0] - .001, zlim[1]);
}

int print_map_palette (mapinfo *mi, double gpver, FILE *fp)
{
    const double *zlim = mi->zrange->val;
    const char *p = NULL;

    if (mi->opts != NULL) {
        /* try for a caller-specified palette */
        p = gretl_bundle_get_string(mi->opts, "palette", NULL);
    }

    if (mi->n_codes > 0) {
	/* the payload series is an encoding */
	if (p != NULL && *p != '\0') {
            return print_discrete_map_palette(mi, p, gpver, fp);
	} else {
            print_discrete_autocolors(mi, fp);
	    if (gpver < 5.4) {
		if (mi->zlabels != NULL || mi->zname != NULL) {
		    print_labeled_colorbox(mi, mi->zlabels,
                                           mi->n_codes, fp);
		} else {
		    fputs("unset colorbox\n", fp);
		}
	    }
	}
    } else {
	/* payload seems to be quantitative */
	if (mi->na_action == NA_FILL) {
	    tricky_print_map_palette(p, zlim, fp);
	    /* note: cbrange handled */
	} else {
	    simple_print_map_palette(p, fp);
	    /* FIXME is this always right (gnuplot version?) */
	    fprintf(fp, "set cbrange [%g:%g]\n", zlim[0], zlim[1]);
	}
    }

    return 0;
}

/* stretch_limits(): allow a little extra space in the X and Y
   dimensions so that the map doesn't entirely fill the plot area; the
   range is scaled by the factor (1 + 2*@margin).
*/

static void stretch_limits (double *targ, const gretl_matrix *minmax,
			    int col, double margin)
{
    double lo = gretl_matrix_get(minmax, 0, col);
    double hi = gretl_matrix_get(minmax, 1, col);
    double mid = 0.5 * (lo + hi);
    double hlf = 0.5 * (hi - lo) * (1 + 2*margin);

    targ[0] = mid - hlf;
    targ[1] = mid + hlf;
}

static void set_limits_from_opts (mapinfo *mi,
				  double *xlim,
				  double *ylim,
				  double margin)
{
    const gretl_matrix *mxy, *mx, *my;

    mxy = gretl_bundle_get_matrix(mi->opts, "mxy__", NULL);
    if (mxy != NULL) {
	xlim[0] = mxy->val[0];
	xlim[1] = mxy->val[1];
	ylim[0] = mxy->val[2];
	ylim[1] = mxy->val[3];
	gretl_bundle_delete_data(mi->opts, "mxy__");
	return;
    }

    mx = gretl_bundle_get_matrix(mi->opts, "xrange", NULL);
    if (mx != NULL && gretl_vector_get_length(mx) == 2) {
	xlim[0] = mx->val[0];
	xlim[1] = mx->val[1];
    } else {
	stretch_limits(xlim, mi->bbox, 0, margin);
    }

    my = gretl_bundle_get_matrix(mi->opts, "yrange", NULL);
    if (my != NULL && gretl_vector_get_length(my) == 2) {
	ylim[0] = my->val[0];
	ylim[1] = my->val[1];
    } else {
	stretch_limits(ylim, mi->bbox, 1, margin);
    }
}

void set_map_plot_limits (mapinfo *mi,
			  double *xlim,
			  double *ylim,
			  double margin)
{
    if (mi->opts != NULL) {
	set_limits_from_opts(mi, xlim, ylim, margin);
    } else {
	stretch_limits(xlim, mi->bbox, 0, margin);
	stretch_limits(ylim, mi->bbox, 1, margin);
    }
}
