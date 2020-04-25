#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "shapefile.h"
#include "libgretl.h"
#include "gretl_typemap.h"

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
	printf("DBFOpen(%s, \"r\") failed\n", dbfname);
	return E_FOPEN;
    }

    fcount = DBFGetFieldCount(DBF);
    if (fcount == 0) {
	DBFClose(DBF);
	printf("There are no fields in this table!\n");
	return E_DATA;
    }

    rcount = DBFGetRecordCount(DBF);
    if (rcount == 0) {
	DBFClose(DBF);
	printf("There are no records in this table!\n");
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

/* shp2dat: Written by Allin Cottrell, 2020-04-13,
   based on shpdump by Frank Warmerdam.

   Outputs the content of the .shp component of a shapefile
   in a form suitable for use by gnuplot.
*/

int shp2dat (const char *shpname,
	     const char *datname,
	     const gretl_matrix *mat)
{
    SHPHandle SHP;
    FILE *fp;
    int n_shapetype, n_entities, i, part;
    double adfmin[4], adfmax[4];
    int pxwidth = 0;
    int pxheight = 0;
    int mercator = 0;
    char delim = ' ';
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
        fprintf(stderr, "Unable to open:%s\n", shpname);
	return E_FOPEN;
    }

    SHPGetInfo(SHP, &n_entities, &n_shapetype, adfmin, adfmax);

    if (mat != NULL) {
	if (mat->rows != n_entities) {
	    fprintf(stderr, "data vector: expected %d rows but got %d\n",
		    n_entities, mat->rows);
	    SHPClose(SHP);
	    return E_DATA;
	}
    }

    fp = gretl_fopen(datname, "wb");
    if (fp == NULL) {
	SHPClose(SHP);
	return E_FOPEN;
    }

    for (i=0; i<n_entities; i++) {
        SHPObject *obj = SHPReadObject(SHP, i);
	double x, y, z = 0;
	int j;

	if (mat != NULL) {
	    z = gretl_matrix_get(mat, i, 0);
	    if (na(z)) {
		SHPDestroyObject(obj);
		continue;
	    }
	}

        if (obj == NULL) {
	    fprintf(stderr, "Unable to read shape %d, terminating.\n", i);
	    err = E_DATA;
            break;
        }

        if (obj->nParts > 0 && obj->panPartStart[0] != 0) {
            fprintf(stderr, "panPartStart[0] = %d, not zero as expected.\n",
		    obj->panPartStart[0]);
	    err = E_DATA;
            break;
        }

        for (j=0, part=1; j<obj->nVertices; j++) {
	    if (part < obj->nParts && obj->panPartStart[part] == j) {
		part++;
		fputc('\n', fp);
            }
	    if (mercator) {
		mercatorize(obj->padfY[j], obj->padfX[j],
			    pxwidth, pxheight, &x, &y);
	    } else {
		x = obj->padfX[j];
		y = obj->padfY[j];
	    }
	    if (mat != NULL) {
		fprintf(fp, "%.*g%c%.*g%c%.*g\n", prec, x, delim, prec, y, delim, prec, z);
	    } else {
		fprintf(fp, "%.*g%c%.*g\n", prec, x, delim, prec, y);
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

    return err;
}

static gchar *put_ext (gchar *fname, const char *ext)
{
    gchar *p = strrchr(fname, '.');

    strcat(p, ext);
    return fname;
}

static int do_geojson (const char *fname,
		       const char *csvname)
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
	gretl_array *features;
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

    return err;
}

int map_get_data (const char *fname, DATASET *dset,
		  gretlopt opt, PRN *prn)
{
    enum { DBF, SHP, GEO };
    gchar *csvname = NULL;
    gchar *base = NULL;
    int ftype;
    int err = 0;

    /* FIXME: allow appending, or just "open"? */

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
	err = do_geojson(fname, csvname);
    } else {
	gchar *dbfname = g_strdup(fname);

	if (ftype == SHP) {
	    put_ext(dbfname, ".dbf");
	}
	err = dbf2csv(dbfname, csvname, OPT_NONE);
	g_free(dbfname);
    }

    if (!err) {
	err = import_csv(csvname, dset, opt, prn);
    }

    /* FIXME: use dset->descrip to record map data location? */

    g_free(base);
    g_free(csvname);

    return err;
}
