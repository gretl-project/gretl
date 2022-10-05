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

    if (plv > 0) {
	/* check the payload series for the discrete property,
	   and for whether it is string-valued
	*/
	; /* TODO */
    }

    if (plx != NULL) {
        /* convert payload series to vector for convenience in plugin */
        mi.zvec = gretl_vector_from_series(plx, dset->t1, dset->t2);
        if (mi.zvec == NULL) {
            err = E_ALLOC;
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
