#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "shapefile.h"
#include "libgretl.h"

static int is_int_string (const char *s)
{
    const char *d = "0123456789";

    return strlen(s) == strspn(s, d);
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
		    if (is_int_string(s)) {
			fputs(s, fp);
		    } else {
			fprintf(fp, "\"%s\"", s);
		    }
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

    fclose(fp);
    SHPClose(SHP);

    return err;
}
