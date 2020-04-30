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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "shapefile.h"
#include "libgretl.h"
#include "gretl_typemap.h"

enum { DBF, SHP, GEO };

#define GEOHUGE 1.0e100

static char *get_fullpath (char *fname)
{
    if (!g_path_is_absolute(fname)) {
	gretl_addpath(fname, 0);
    }

    return fname;
}

static int matrix_is_payload (const gretl_matrix *mat)
{
    int i;

    for (i=0; i<mat->rows; i++) {
	if (!na(mat->val[i]) && mat->val[i] != 0) {
	    return 1;
	}
    }
    return 0;
}

static int skip_object (int i, const gretl_matrix *m, double *pz)
{
    if (m == NULL) {
	return 0;
    } else {
	*pz = m->val[i];
	return na(*pz);
    }
}

static void record_extrema (double x, double y,
			    double *fmin, double *gmax)
{
    if (x < fmin[0]) {
	fmin[0] = x;
    }
    if (x > gmax[0]) {
	gmax[0] = x;
    }
    if (y < fmin[1]) {
	fmin[1] = y;
    }
    if (y > gmax[1]) {
	gmax[1] = y;
    }
}

static gretl_matrix *ring2matrix (gretl_array *ring)
{
    int i, n = gretl_array_get_length(ring);
    gretl_matrix *ret = gretl_matrix_alloc(n, 2);
    const char *sx, *sy;
    gretl_array *ri;

    for (i=0; i<n; i++) {
	ri = gretl_array_get_data(ring, i);
	sx = gretl_array_get_data(ri, 0);
	sy = gretl_array_get_data(ri, 1);
	gretl_matrix_set(ret, i, 0, atof(sx));
	gretl_matrix_set(ret, i, 1, atof(sy));
    }

    return ret;
}

static gretl_array *geojson_get_features (const char *fname,
					  int *err)
{
    gretl_bundle *(*jfunc) (const char *, const char *, int *);
    gretl_array *a = NULL;
    GError *gerr = NULL;
    gchar *JSON = NULL;
    gsize len = 0;
    gboolean ok;

    ok = g_file_get_contents(fname, &JSON, &len, &gerr);

    if (ok) {
	gretl_bundle *b = NULL;
	GretlType type = 0;

	jfunc = get_plugin_function("json_get_bundle");
	if (jfunc == NULL) {
	    *err = E_FOPEN;
	} else {
	    b = jfunc(JSON, NULL, err);
	    if (!*err) {
		a = gretl_bundle_steal_data(b, "features", &type, NULL, err);
	    }
	}
	gretl_bundle_destroy(b);
	g_free(JSON);
    } else if (gerr != NULL) {
	gretl_errmsg_set(gerr->message);
	g_error_free(gerr);
    }

    return a;
}

static gretl_matrix *geo2dat (const char *geoname,
			      const char *datname,
			      const gretl_matrix *zvec)
{
    gretl_array *features, *AC;
    gretl_array *ACj, *ACjk;
    gretl_matrix *X, *bbox = NULL;
    gretl_bundle *fi, *geom;
    double gmin[2] = {GEOHUGE, GEOHUGE};
    double gmax[2] = {-GEOHUGE, -GEOHUGE};
    FILE *fp;
    const char *gtype;
    int have_payload = 0;
    int nf, mp, nac, ncj;
    int i, j, k, p;
    int err = 0;

    features = geojson_get_features(geoname, &err);
    if (features == NULL) {
	return NULL;
    }

    if (zvec != NULL) {
	if (matrix_is_payload(zvec)) {
	    have_payload = 1;
	}
    }

    fp = gretl_fopen(datname, "wb");
    if (fp == NULL) {
	gretl_array_destroy(features);
	return NULL;
    }

    nf = gretl_array_get_length(features);

    for (i=0; i<nf; i++) {
	double x, y, z = 0;

	if (skip_object(i, zvec, &z)) {
	    continue;
	}

        if (!na(z)) {
	    fi = gretl_array_get_data(features, i);
	    geom = gretl_bundle_get_bundle(fi, "geometry", NULL);
	    gtype = gretl_bundle_get_string(geom, "type", NULL);
	    if (!strcmp(gtype, "Polygon")) {
                mp = 0;
	    } else if (!strcmp(gtype, "MultiPolygon")) {
                mp = 1;
            } else {
                gretl_errmsg_sprintf("can't handle geometry type '%s'", gtype);
		err = E_DATA;
		break;
	    }
	    AC = gretl_bundle_get_array(geom, "coordinates", NULL);
            nac = gretl_array_get_length(AC);
            if (mp == 0) {
		/* got Polygon */
                for (j=0; j<nac; j++) {
		    ACj = gretl_array_get_data(AC, j);
                    X = ring2matrix(ACj);
		    for (k=0; k<X->rows; k++) {
			x = gretl_matrix_get(X, k, 0);
			y = gretl_matrix_get(X, k, 1);
			if (have_payload) {
                            fprintf(fp, "%.8g %.8g %.8g\n", x, y, z);
			} else {
			    fprintf(fp, "%.8g %.8g\n", x, y);
			}
			record_extrema(x, y, gmin, gmax);
		    }
		    gretl_matrix_free(X);
                    if (j < nac-1) {
                        fputc('\n', fp);
		    }
		}
	    } else {
		/* got MultiPolygon */
		for (j=0; j<nac; j++) {
		    ACj = gretl_array_get_data(AC, j);
		    ncj = gretl_array_get_length(ACj);
		    for (k=0; k<ncj; k++) {
			ACjk = gretl_array_get_data(ACj, k);
                        X = ring2matrix(ACjk);
			for (p=0; p<X->rows; p++) {
			    x = gretl_matrix_get(X, p, 0);
			    y = gretl_matrix_get(X, p, 1);
			    if (have_payload) {
				fprintf(fp, "%.8g %.8g %.8g\n", x, y, z);
			    } else {
				fprintf(fp, "%.8g %.8g\n", x, y);
			    }
			    record_extrema(x, y, gmin, gmax);
			}
			gretl_matrix_free(X);
                        if (k < ncj-1) {
                            fputc('\n', fp);
                        }
		    }
		    if (j < nac-1) {
                        fputc('\n', fp);
                    }
                }
            }
            if (i < nf-1) {
                fputs("\n\n", fp); /* end of entity block */
            }
        }
    }

    fputc('\n', fp);
    fclose(fp);

    gretl_array_destroy(features);

    if (!err) {
	bbox = gretl_matrix_alloc(2, 2);
	if (bbox != NULL) {
	    gretl_matrix_set(bbox, 0, 0, gmin[0]);
	    gretl_matrix_set(bbox, 0, 1, gmin[1]);
	    gretl_matrix_set(bbox, 1, 0, gmax[0]);
	    gretl_matrix_set(bbox, 1, 1, gmax[1]);
	}
    }

    return bbox;
}

static void output_dbf_string (const char *s, FILE *fp)
{
    fputc('"', fp);
    while (*s) {
	if (*s != 0x0d && *s != 0x0a) {
	    fputc(*s, fp);
	}
	s++;
    }
    fputc('"', fp);
}

/* dbf2csv: Written by Allin Cottrell, 2020-04-13,
   based on dbfdump by Frank Warmerdam.

   Outputs the content of the .dbf component of a shapefile
   (metadata) as CSV.
*/

int dbf2csv (const char *dbfname,
	     const char *csvname,
	     gretlopt opt)
{
    DBFHandle DBF;
    FILE *fp;
    int width, decimals;
    int fcount, rcount;
    DBFFieldType etype;
    int header = 0;
    char title[32];
    int i, j;

    if (opt & OPT_H) {
	header = 1;
    }

    DBF = DBFOpen(dbfname, "rb");
    if (DBF == NULL) {
	gretl_errmsg_sprintf("DBFOpen(%s) failed", dbfname);
	return E_FOPEN;
    }

    fcount = DBFGetFieldCount(DBF);
    if (fcount == 0) {
	DBFClose(DBF);
	gretl_errmsg_set("There are no fields in this DBF table!");
	return E_DATA;
    }

    rcount = DBFGetRecordCount(DBF);
    if (rcount == 0) {
	DBFClose(DBF);
	gretl_errmsg_set("There are no records in this DBF table!");
	return E_DATA;
    }

    fp = gretl_fopen(csvname, "wb");
    if (fp == NULL) {
	DBFClose(DBF);
	return E_FOPEN;
    }

    if (header) {
        for (i=0; i<fcount; i++) {
            const char *typename = NULL;
            char native_type;

            native_type = DBFGetNativeFieldType(DBF, i);

            etype = DBFGetFieldInfo(DBF, i, title, &width, &decimals);
            if (etype == FTString)
                typename = "String";
            else if (etype == FTInteger)
                typename = "Integer";
            else if (etype == FTDouble)
                typename = "Double";
            else if (etype == FTInvalid)
                typename = "Invalid";

            fprintf(fp, "# Field %d: Type=%c/%s, Title=`%s', Width=%d, Decimals=%d\n",
		    i, native_type, typename, title, width, decimals);
        }
	fprintf(fp, "# Number of records in table: %d\n", rcount);
    }

    /* print column headings */
    for (i=0; i<fcount; i++) {
	etype = DBFGetFieldInfo(DBF, i, title, &width, &decimals);
	fputs(title, fp);
	if (i < fcount - 1) {
	    fputc(',', fp);
	}
    }
    fputc('\n', fp);

    /* Read all the records */
    for (j=0; j<rcount; j++) {
	for (i=0; i<fcount; i++) {
            etype = DBFGetFieldInfo(DBF, i, title, &width, &decimals);

	    if (DBFIsAttributeNULL(DBF, j, i)) {
		fputs("(NULL)", fp);
	    } else {
		const char *s;

		switch (etype) {
		case FTString:
		    s = DBFReadStringAttribute(DBF, j, i);
		    output_dbf_string(s, fp);
		    break;
		case FTInteger:
		    fprintf(fp, "%d", DBFReadIntegerAttribute(DBF, j, i));
		    break;
		case FTDouble:
		    fprintf(fp, "%.8g", DBFReadDoubleAttribute(DBF, j, i));
		    break;
		default:
		    break;
		}
	    }
	    if (i < fcount - 1) {
		fputc(',', fp);
	    }
	    fflush(fp);
	}

        if (DBFIsRecordDeleted(DBF, j)) {
            fputs("(DELETED)", fp);
	}
	fputc('\n', fp);
    }

    fclose(fp);
    DBFClose(DBF);

    return 0;
}

static void mercatorize (double lat, double lon,
			 int width, int height,
			 double *px, double *py)
{
    double x = (lon + 180) * (width / 360.0);
    double rad = lat * M_PI/180.0;
    double merc = log(tan(M_PI/4.0 + rad/2.0));
    double y = height/2.0 - height*merc/(2*M_PI);

    *px = x;
    *py = -y;
}

static gretl_matrix *shp2dat (const char *shpname,
			      const char *datname,
			      const gretl_matrix *zvec)
{
    gretl_matrix *bbox = NULL;
    SHPHandle SHP;
    FILE *fp;
    int n_shapetype, n_entities, i, part;
    double gmin[4], gmax[4];
    int have_payload = 0;
    int pxwidth = 0;
    int pxheight = 0;
    int mercator = 0;
    char delim = ' ';
    int nskip = 0;
    int prec = 8;
    int err = 0;

#if 0 /* not yet */
    /* extra args?? mercvec, prec, opt */
    if (opt & OPT_M) {
	mercator = 1;
    }

    if (prec == 0) {
	prec = 8;
    }

    if (mercvec != NULL) {
	/* FIXME */
	mercator = 1;
	pxwidth = mercvec->val[0];
	pxheight = mercvec->val[1];
    }
#endif

    SHP = SHPOpen(shpname, "rb");
    if (SHP == NULL) {
	return NULL;
    }

    SHPGetInfo(SHP, &n_entities, &n_shapetype, gmin, gmax);

    if (zvec != NULL) {
	if (zvec->rows != n_entities) {
	    fprintf(stderr, "data vector: expected %d rows but got %d\n",
		    n_entities, zvec->rows);
	    SHPClose(SHP);
	    return NULL;
	}
	for (i=0; i<zvec->rows; i++) {
	    if (na(zvec->val[i])) {
		nskip++;
		break;
	    }
	}
	if (nskip > 0) {
	    for (i=0; i<2; i++) {
		gmin[i] = GEOHUGE;
		gmax[i] = -GEOHUGE;
	    }
	}
	if (matrix_is_payload(zvec)) {
	    have_payload = 1;
	}
    }

    fp = gretl_fopen(datname, "wb");
    if (fp == NULL) {
	SHPClose(SHP);
	return NULL;
    }

    SHPSetFastModeReadObject(SHP, TRUE);

    for (i=0; i<n_entities; i++) {
        SHPObject *obj;
	double x, y, z = 0;
	int j;

	if (skip_object(i, zvec, &z)) {
	    continue;
	}

	obj = SHPReadObject(SHP, i);

        if (obj == NULL) {
	    fprintf(stderr, "Unable to read shape %d, terminating.\n", i);
	    err = E_DATA;
        } else if (obj->nParts > 0 && obj->PartStart[0] != 0) {
            fprintf(stderr, "PartStart[0] = %d, not zero as expected.\n",
		    obj->PartStart[0]);
	    err = E_DATA;
        }

        for (j=0, part=1; j<obj->nVertices && !err; j++) {
	    if (part < obj->nParts && obj->PartStart[part] == j) {
		part++;
		fputc('\n', fp);
            }
	    if (mercator) {
		mercatorize(obj->fY[j], obj->fX[j],
			    pxwidth, pxheight, &x, &y);
	    } else {
		x = obj->fX[j];
		y = obj->fY[j];
	    }
	    if (have_payload) {
		fprintf(fp, "%.*g%c%.*g%c%.*g\n", prec, x, delim, prec, y, delim, prec, z);
	    } else {
		fprintf(fp, "%.*g%c%.*g\n", prec, x, delim, prec, y);
	    }
	    if (nskip > 0) {
		record_extrema(x, y, gmin, gmax);
	    }
        }

	if (i < n_entities - 1) {
	    fputs("\n\n", fp);
	}

        SHPDestroyObject(obj);
    }

    fputc('\n', fp);
    fclose(fp);
    SHPClose(SHP);

    if (!err) {
	bbox = gretl_matrix_alloc(2, 2);
	if (bbox != NULL) {
	    gretl_matrix_set(bbox, 0, 0, gmin[0]);
	    gretl_matrix_set(bbox, 0, 1, gmin[1]);
	    gretl_matrix_set(bbox, 1, 0, gmax[0]);
	    gretl_matrix_set(bbox, 1, 1, gmax[1]);
	}
    }

    return bbox;
}

static char *put_ext (char *fname, const char *ext)
{
    char *p = strrchr(fname, '.');

    *p = '\0';
    strcat(p, ext);
    return fname;
}

static int do_geojson (const char *fname,
		       const char *csvname,
		       char **mapname)
{
    GError *gerr = NULL;
    gchar *JSON = NULL;
    gsize len = 0;
    gboolean ok;
    int err = 0;

    ok = g_file_get_contents(fname, &JSON, &len, &gerr);

    if (!ok) {
	if (gerr != NULL) {
	    gretl_errmsg_set(gerr->message);
	    g_error_free(gerr);
	} else {
	    fprintf(stderr, "g_file_get_contents failed for '%s'\n", fname);
	}
	err = E_DATA;
    } else {
	gretl_bundle *(*jfunc) (const char *, const char *, int *);
	FILE *fp = NULL;
	void *ptr = NULL;
	GretlType type;
	gretl_bundle *jb, *pp, *fi;
	gretl_array *features = NULL;
	gretl_array *keys = NULL;
	const char *key;
	int nf = 0, nk = 0;
	int i, j;

	jfunc = get_plugin_function("json_get_bundle");
	if (jfunc == NULL) {
	    return E_DATA;
	}
	fp = gretl_fopen(csvname, "wb");
	if (fp == NULL) {
	    gretl_errmsg_sprintf(_("Couldn't open %s for writing"), csvname);
	    return E_FOPEN;
	}
	jb = jfunc(JSON, NULL, &err);
	if (jb == NULL) {
	    gretl_errmsg_sprintf(_("Couldn't find function %s"), "json_get_bundle");
	    err = E_DATA;
	}
	if (!err) {
	    features = gretl_bundle_get_array(jb, "features", &err);
	    if (err) {
		gretl_errmsg_sprintf(_("Couldn't read '%s'"), "features");
	    }
	}
	if (!err) {
	    nf = gretl_array_get_length(features);
	    fi = gretl_array_get_element(features, 0, NULL, &err);
	}
	if (!err) {
	    pp = gretl_bundle_get_bundle(fi, "properties", &err);
	}
	if (!err) {
	    keys = gretl_bundle_get_keys(pp, &err);
	}
	if (!err) {
	    nk = gretl_array_get_length(keys);
	}
	for (j=0; j<nk && !err; j++) {
	    key = gretl_array_get_data(keys, j);
	    fprintf(fp, "%s%c", key, j < nk-1 ? ',' : '\n');
	}
	for (i=0; i<nf && !err; i++) {
	    fi = gretl_array_get_element(features, i, NULL, &err);
	    pp = gretl_bundle_get_bundle(fi, "properties", &err);
	    for (j=0; j<nk && !err; j++) {
		key = gretl_array_get_data(keys, j);
		ptr = gretl_bundle_get_data(pp, key, &type, NULL, &err);
		if (err) {
		    fprintf(stderr, "error at feature %d, propkey %d\n", i, j);
		    break;
		}
		if (type == GRETL_TYPE_STRING) {
		    fprintf(fp, "\"%s\"", (char *) ptr);
		} else if (type == GRETL_TYPE_INT) {
		    fprintf(fp, "%d", *(int *) ptr);
		} else if (type == GRETL_TYPE_DOUBLE) {
		    fprintf(fp, "%g", *(double *) ptr);
		} else {
		    fprintf(stderr, "Got property type %s\n", gretl_type_get_name(type));
		    fprintf(fp, "\"\"");
		}
		fputc(j < nk-1 ? ',' : '\n', fp);
	    }
	}

	gretl_array_destroy(keys);
	gretl_bundle_destroy(jb);
	fclose(fp);
    }

    g_free(JSON);

    if (!err) {
	*mapname = gretl_strdup(fname);
    }

    return err;
}

static int do_shapefile (const char *fname,
			 const char *csvname,
			 int ftype,
			 char **mapname)
{
    char *dbfname = gretl_strdup(fname);
    char *shpname = gretl_strdup(fname);
    char *shxname = gretl_strdup(fname);
    int err;

    if (ftype == SHP) {
	put_ext(dbfname, ".dbf");
    } else if (ftype == DBF) {
	put_ext(shpname, ".shp");
    }
    put_ext(shxname, ".shx");

    if (gretl_stat(dbfname, NULL) != 0) {
	gretl_errmsg_sprintf(_("Couldn't open '%s'"), dbfname);
	err = E_FOPEN;
	goto bailout;
    }
    if (gretl_stat(shpname, NULL) != 0) {
	gretl_errmsg_sprintf(_("Couldn't open '%s'"), shpname);
	err = E_FOPEN;
	goto bailout;
    }
    if (gretl_stat(shxname, NULL) != 0) {
	gretl_errmsg_sprintf(_("Couldn't open '%s'"), shxname);
	err = E_FOPEN;
	goto bailout;
    }

    err = dbf2csv(dbfname, csvname, OPT_NONE);

    if (!err) {
	*mapname = shpname;
	shpname = NULL;
    }

 bailout:

    free(dbfname);
    free(shpname);
    free(shxname);

    return err;
}

gretl_matrix *map2dat (const char *mapname,
		       const char *datname,
		       const gretl_matrix *zvec)
{
    char infile[MAXLEN];

    strcpy(infile, mapname);
    get_fullpath(infile);

    if (has_suffix(mapname, ".shp")) {
	return shp2dat(infile, datname, zvec);
    } else {
	return geo2dat(infile, datname, zvec);
    }
}

int map_get_data (const char *fname, DATASET *dset,
		  gretlopt opt, PRN *prn)
{
    gchar *csvname = NULL;
    gchar *base = NULL;
    char *mapname = NULL;
    int ftype;
    int err = 0;

    /* FIXME: should appending be allowed, or just "open"? */

    if (has_suffix(fname, ".dbf")) {
	ftype = DBF;
    } else if (has_suffix(fname, ".shp")) {
	ftype = SHP;
    } else {
	ftype = GEO;
    }

    base = g_path_get_basename(fname);
    csvname = gretl_make_dotpath(base);
    put_ext(csvname, ".csv");

    if (ftype == GEO) {
	err = do_geojson(fname, csvname, &mapname);
    } else {
	err = do_shapefile(fname, csvname, ftype, &mapname);
    }

    if (!err) {
	err = import_csv(csvname, dset, opt, prn);
	if (!err) {
	    dset->mapfile = mapname;
	    mapname = NULL;
	}
    }

    g_free(base);
    g_free(csvname);
    free(mapname);

    return err;
}
