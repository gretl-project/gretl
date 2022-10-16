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

/* mapinfo.h for gretl: the role of this file is to facilitate
   communication between the various phases of producing a
   "geoplot" map, the code for which is found in in mapinfo.c,
   plugin/geoplot.c and graphing.c.
*/

#ifndef MAPINFO_H
#define MAPINFO_H

#include "libgretl.h"

typedef enum {
    PRJ0,      /* geoplot default: quasi-Mercator */
    WGS84,     /* null projection */
    EPSG3857,  /* "proper" Mercator */
    EPSG2163,  /* US National Atlas Equal Area */
    EPSG3035,  /* Europe Equal Area */
} MapProj;

typedef enum {
    MAP_NON_STD  = 1 << 0,
    MAP_IS_IMAGE = 1 << 1,
    MAP_DISPLAY  = 1 << 2,
    MAP_DUMMY    = 1 << 3
} MapFlags;

typedef enum {
    NA_SKIP,     /* exclude regions with missing payload */
    NA_OUTLINE,  /* show outlines of such areas */
    NA_FILL      /* give such regions a specific fill color */
} NaAction;

typedef struct mapinfo_ mapinfo;

struct mapinfo_ {
    const char *mapfile;   /* name of file containing polygons */
    gretl_bundle *map;     /* bundle containing polygons */
    gretl_matrix *zvec;    /* "payload" data to plot, or NULL */
    char *zname;           /* name of payload series */
    char *plotfile;        /* plot filename specified by caller */
    char *datfile;         /* auxiliary file: gnuplot data */
    gretl_matrix *bbox;    /* plot bounding box */
    gretl_matrix *zrange;  /* range of the payload data */
    gretl_matrix *zvals;   /* ordered unique values (discrete data) */
    char **zlabels;        /* labels for @zvals (discrete data) */
    int n_codes;           /* number of coded payload values */
    gretl_bundle *opts;    /* options specified by caller, or NULL */
    MapFlags flags;        /* state flags */
    NaAction na_action;    /* method for handling missing values */
    MapProj proj;          /* choice of projection */
};

int geoplot_driver (const char *fname,
                    gretl_bundle *map,
		    int plv,
                    const double *plx,
                    const DATASET *dset,
                    gretl_bundle *opts);

int print_map_palette (mapinfo *mi, double gpver, FILE *fp);

void set_map_plot_limits (mapinfo *mi,
			  double *xlim,
			  double *ylim,
			  double margin);

#endif /* MAPINFO_H */
