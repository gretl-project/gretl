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

/* check for the canonical case, where min = 1 and all values
   are consecutive
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

/* check the payload series for the discrete property, for
   whether it is string-valued, and for whether it is
   "canonical" in these respects.
*/

static int inspect_payload (mapinfo *mi, const DATASET *dset, int v)
{
    char **strvals = NULL;
    gretl_matrix *vals = NULL;
    int discrete = 0;
    int canon = 0;
    int ns = 0;
    int nv = 0;
    int err = 0;

    if (is_string_valued(dset, v)) {
	strvals = series_get_string_vals(dset, v, &ns, 1);
	discrete = 1;
    } else if (series_is_discrete(dset, v)) {
	discrete = 1;
    }

    if (discrete) {
	int nz = gretl_vector_get_length(mi->zvec);
	const double *z = mi->zvec->val;

	vals = gretl_matrix_values(z, nz, OPT_S, &err);
	if (!err) {
	    nv = gretl_vector_get_length(vals);
	    canon = canonical_discrete(vals->val, nv);
	}
    }

    if (discrete && !err) {
	if (canon) {
	    mi->zvals = vals;
	    mi->flags |= MAP_DISCRETE;
	    if (ns == nv) {
		mi->zlabels = strvals;
	    }
	} else {
	    fprintf(stderr, "discrete not canonical: up to the caller?\n");
	    gretl_matrix_free(vals); /* don't leak memory */
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
            map = get_current_map(dset, NULL, &err);
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

/* [example] palette = "set palette maxcolors 4; \
   set palette defined (0 '#D65E5E', 1 '#8594E1', \
   2 '#85E1C3', 3 '#E1C385'); unset colorbox"
*/

static int print_discrete_map_palette (const char *s,
				       const gretl_matrix *zrange,
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
	gretl_array *S[2] = {NULL};
	int i, j, m[2] = {0};

	for (i=0; i<n; i++) {
	    S[i] = get_array_by_name(anames[i]);
	    if (gretl_array_get_type(S[i]) != GRETL_TYPE_STRINGS ||
		(m[i] = gretl_array_get_length(S[i])) < 2) {
		gretl_errmsg_sprintf("Invalid discrete palette argument '%s'", anames[i]);
		err = E_INVARG;
	    }
	}
	if (!err && n == 2 && m[1] != m[0]) {
	    gretl_errmsg_sprintf("Invalid discrete palette argument '%s'", anames[i]);
	    err = E_INVARG;
	}
	if (!err) {
	    fprintf(fp, "set palette maxcolors %d\n", m[0]);
	    fputs("set palette defined (", fp);
	    for (j=0; j<m[0]; j++) {
		s = gretl_array_get_data(S[0], j);
		fprintf(fp, "%d '%s'", j, s);
		if (j < m[0]-1) {
		    fputs(", ", fp);
		}
	    }
	    fputs(")\n", fp);
	    if (n == 2) {
		double zmin = zrange->val[0];
		double incr = (m[0] - 1.0) / m[0];
		double x = zmin + incr / 2.0;

		fputs("set cbtics (", fp);
		for (j=0; j<m[0]; j++) {
		    s = gretl_array_get_data(S[1], j);
		    fprintf(fp, "'%s' %g", s, x);
		    x += incr;
		    if (j < m[0]-1) {
			fputs(", ", fp);
		    }
		}
		fputs(") scale 0\n", fp);
	    } else {
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
    if (!strcmp(p, "blues")) {
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

int print_map_palette (mapinfo *mi, FILE *fp)
{
    const double *zlim = mi->zrange->val;
    const char *p;

    p = gretl_bundle_get_string(mi->opts, "palette", NULL);

    if (p == NULL || *p == '\0') {
	return 0;
    } else if (!strncmp(p, "discrete,", 9)) {
        return print_discrete_map_palette(p + 9, mi->zrange, fp);
    } else if (mi->na_action == NA_FILL) {
	tricky_print_map_palette(p, zlim, fp);
	/* cbrange handled */
	return 0;
    } else {
	simple_print_map_palette(p, fp);
    }

    fprintf(fp, "set cbrange [%g:%g]\n", zlim[0], zlim[1]);

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
