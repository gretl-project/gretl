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

#include "libgretl.h"
#include "dbwrite.h"

/**
 * SECTION:dbwrite
 * @short_description: writing to a gretl database
 * @title: DB write
 * @include: gretl/libgretl.h, gretl/dbwrite.h
 *
 * Functionality for writing series to a native-format gretl database.
 */

#define DB_DEBUG 0

static void dotify (char *s)
{
    while (*s) {
	if (*s == ':') *s = '.';
	s++;
    }
}

static char pd_char (const DATASET *dset)
{
    if (dset->pd == 4) {
	return 'Q';
    } else if (dset->pd == 12) {
	return 'M';
    } else if (dataset_is_time_series(dset)) {
	return 'A';
    } else {
	return 'U';
    }
}

static int get_db_series_names (const char *idxname, char ***pnames,
				int *pnv)
{
    char line[256];
    char **vnames = NULL;
    FILE *fp;
    int i, j, nv;
    int err = 0;

    fp = gretl_fopen(idxname, "r");
    if (fp == NULL) {
	return E_FOPEN;
    }

#if DB_DEBUG
    fprintf(stderr, "get_db_series_names: opened %s\n", idxname);
#endif

    /* first pass: count the number of vars */
    i = nv = 0;
    while (fgets(line, sizeof line, fp)) {
	if (*line == '#' || string_is_blank(line)) {
	    continue;
	}
	i++;
	if (i % 2) {
	    /* odd-numbered lines hold varnames */
	    nv++;
	}
    }

#if DB_DEBUG
    fprintf(stderr, " found %d varnames\n", nv);
#endif

    if (nv == 0) {
	err = E_DATA;
    } else {
	vnames = strings_array_new_with_length(nv, VNAMELEN);
	if (vnames == NULL) {
	    err = E_ALLOC;
	}
    }

    if (err) {
	fclose(fp);
	return err;
    }

    rewind(fp);

    /* second pass: grab all the varnames */
    i = j = 0;
    while (fgets(line, sizeof line, fp) && !err) {
	if (*line == '#' || string_is_blank(line)) {
	    continue;
	}
	i++;
	if (i % 2) {
	    if (gretl_scan_varname(line, vnames[j]) != 1) {
		err = E_DATA;
	    }
	    j++;
	}
    }    

    fclose(fp);

    if (err) {
	strings_array_free(vnames, nv);
    } else {
	*pnames = vnames;
	*pnv = nv;
    }

    return err;
}

/* Given a list of variables to be appended to a gretl database,
   check that there is not already a variable in the database
   with the same name as any of those in the list.  Return 0
   if no duplicates, -1 on failure, or the positive number 
   of duplicated variables.
*/

static int 
check_for_db_duplicates (const int *list, const DATASET *dset,
			 const char *idxname, int *err)
{
    char **snames = NULL;
    int i, j, v, oldv = 0;
    int ret = 0;

#if DB_DEBUG   
    printlist(list, "check_for_db_duplicates: input save list");
#endif

    *err = get_db_series_names(idxname, &snames, &oldv);
    if (*err) {
	return -1;
    }

    for (i=1; i<=list[0]; i++) {
	v = list[i];
	for (j=0; j < oldv; j++) {
	    if (!strcmp(dset->varname[v], snames[j])) {
		ret++;
		break;
	    }
	}
    }

    strings_array_free(snames, oldv);

    return ret;
}

static int output_db_var (int v, const DATASET *dset,
			  FILE *fidx, FILE *fbin) 
{
    char stobs[OBSLEN], endobs[OBSLEN];
    int t1 = dset->t1;
    int t2 = dset->t2;
    int t, nobs;
    float val;

    if (dataset_is_time_series(dset)) {
	for (t=dset->t1; t<=dset->t2; t++) {
	    if (na(dset->Z[v][t])) t1++;
	    else break;
	}
	for (t=dset->t2; t>=t1; t--) {
	    if (na(dset->Z[v][t])) t2--;
	    else break;
	}
    }

    nobs = t2 - t1 + 1;
    if (nobs <= 0) {
	return 0;
    }

    ntodate(stobs, t1, dset);
    ntodate(endobs, t2, dset);
    dotify(stobs);
    dotify(endobs);	

    fprintf(fidx, "%s  %s\n", dset->varname[v], series_get_label(dset, v));
    fprintf(fidx, "%c  %s - %s  n = %d\n", pd_char(dset),
	    stobs, endobs, nobs);

    for (t=t1; t<=t2; t++) {
	if (na(dset->Z[v][t])) {
	    val = DBNA;
	} else {
	    val = dset->Z[v][t];
	}
	fwrite(&val, sizeof val, 1, fbin);
    }

    return 0;
}

static int write_old_bin_chunk (long offset, int nvals, 
				FILE *fin, FILE *fout)
{
    int i, err = 0;
    float val;

    fseek(fin, offset, SEEK_SET);

    for (i=0; i<nvals && !err; i++) {
	if (fread(&val, sizeof val, 1, fin) != 1) {
	    err = 1;
	} 
	if (!err) {
	    if (fwrite(&val, sizeof val, 1, fout) != 1) {
		err = 1;
	    }
	}
    } 

    return err;
}

static void list_delete_element (int *list, int m)
{
    int i;

    for (i=1; i<=list[0]; i++) {
	if (list[i] == m) {
	    gretl_list_delete_at_pos(list, i);
	    break;
	}
    }
}

/* writing to a previously existing database, replacing any existing
   variables with the same name as "new" ones, but otherwise
   preserving the existing content 
*/

static int 
append_db_data_with_replacement (const char *idxname, 
				 const char *binname,
				 int *list, 
				 const DATASET *dset) 
{
    FILE *fidx = NULL, *fbin = NULL;
    char **oldnames = NULL;
    char *mask = NULL;
    int *newlist = NULL;
    int i, j, v, oldv;
    int nrep, err = 0;

    err = get_db_series_names(idxname, &oldnames, &oldv);
    if (err) {
	return err;
    }

    mask = calloc(oldv, 1);
    if (mask == NULL) {
	err = E_ALLOC;
	goto bailout;
    }	

    newlist = gretl_list_copy(list);
    if (newlist == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    nrep = 0;
    for (i=1; i<=list[0]; i++) {
	v = list[i];
	for (j=0; j<oldv; j++) {
	    if (!strcmp(oldnames[j], dset->varname[v])) {
		/* match: remove var v from "newlist" and flag that it
		   is a replacement in "mask" */
		list_delete_element(newlist, v);
#if DB_DEBUG
		fprintf(stderr, "match: var %d and old db var %d\n", v, j);
#endif    
		mask[j] = 1;
		nrep++;
		break;
	    }
	}
    }

#if DB_DEBUG
    fprintf(stderr, "write_db_data_with_replacement: replicated vars = %d\n", 
	    nrep);
    printlist(list, "full var list");
    printlist(newlist, "new var list");
#endif    

    /* handle replacement variables first */
    if (nrep > 0) {
	char idxcpy[FILENAME_MAX];
	char bincpy[FILENAME_MAX];

	strcpy(idxcpy, idxname);
	strcat(idxcpy, ".cpy");

	strcpy(bincpy, binname);
	strcat(bincpy, ".cpy");

	err = gretl_copy_file(idxname, idxcpy);

	if (!err) {
	    err = gretl_copy_file(binname, bincpy);
	}

	if (!err) {
	    FILE *fp = NULL, *fq = NULL;
	    char line1[256], line2[256];
	    long offset = 0L;
	    int nobs;

	    fp = fq = NULL;
	    
	    fp = gretl_fopen(idxcpy, "r");
	    if (fp == NULL) {
		err = E_FOPEN;
	    }

	    if (!err) {
		fidx = gretl_fopen(idxname, "w");
		if (fidx == NULL) {
		    err = E_FOPEN;
		}
	    }

	    if (!err) {
		fq = gretl_fopen(bincpy, "rb");
		if (fq == NULL) {
		    err = E_FOPEN;
		}
	    }

	    if (!err) {
		fbin = gretl_fopen(binname, "wb");
		if (fbin == NULL) {
		    err = E_FOPEN;
		}
	    }	    

	    i = 0;
	    while (fgets(line1, sizeof line1, fp) && !err) {
		if (*line1 == '#' || string_is_blank(line1)) {
		    if (*line1 == '#') {
			fputs(line1, fidx);
		    }	
		    continue;
		}
		if (fgets(line2, sizeof line2, fp) == NULL) {
		    /* db index lines must be in pairs */
		    err = 1;
		    break;
		}
		if (sscanf(line2, "%*s  %*s - %*s  n = %d\n", &nobs) != 1) {
		    err = 1;
		    break;
		}

#if DB_DEBUG
		fprintf(stderr, "old db, var %d, nobs = %d\n", i, nobs);
#endif
		if (mask[i]) {
		    v = series_index(dset, oldnames[i]);
#if DB_DEBUG
		    fprintf(stderr, "replacing this with var %d\n", v);
#endif
		    output_db_var(v, dset, fidx, fbin);
		} else {
#if DB_DEBUG
		    fprintf(stderr, "passing through old var\n");
#endif
		    fputs(line1, fidx);
		    fputs(line2, fidx);
		    write_old_bin_chunk(offset, nobs, fq, fbin);
		}
		i++;
		offset += nobs * sizeof(float);
	    }

	    if (fp != NULL) fclose(fp);
	    if (fq != NULL) fclose(fq);
	}

	gretl_remove(idxcpy);
	gretl_remove(bincpy);
    } else {
	/* no variables to be replaced */
	fidx = gretl_fopen(idxname, "a");
	if (fidx == NULL) {
	    err = E_FOPEN;
	}
	if (!err) {
	    fbin = gretl_fopen(binname, "ab");
	    if (fbin == NULL) {
		err = E_FOPEN;
	    }
	}
    }

    if (!err) {
	/* do any newly added variables */
	for (i=1; i<=newlist[0]; i++) {
#if DB_DEBUG
	    fprintf(stderr, "adding new var, %d\n", newlist[i]);
#endif
	    output_db_var(newlist[i], dset, fidx, fbin);
	}
    }

 bailout:

    if (fidx != NULL) fclose(fidx);
    if (fbin != NULL) fclose(fbin);

    strings_array_free(oldnames, oldv);
    free(mask);
    free(newlist);

    return err;
}

static int 
open_db_files (const char *fname, char *idxname, char *binname,
	       FILE **fidx, FILE **fbin, int *append)
{
    FILE *fp;
    char base[FILENAME_MAX];
    char imode[3] = "w";
    char bmode[3] = "wb";
    char *p;

    strcpy(base, fname);
    p = strchr(base, '.');
    if (p != NULL) {
	*p = 0;
    }

    strcpy(idxname, base);
    strcat(idxname, ".idx");

    fprintf(stderr, "open_db_files: doing test open on '%s'\n", idxname);

    fp = gretl_fopen(idxname, "r");
    if (fp != NULL) {
	*append = 1;
	strcpy(imode, "a");
	strcpy(bmode, "ab");
	fclose(fp);
    }

    *fidx = gretl_fopen(idxname, imode);
    if (*fidx == NULL) {
	gretl_errmsg_sprintf(_("Couldn't open %s for writing"), idxname);
	return 1;
    }

    strcpy(binname, base);
    strcat(binname, ".bin");
    
    *fbin = gretl_fopen(binname, bmode);
    if (*fbin == NULL) {
	gretl_errmsg_sprintf(_("Couldn't open %s for writing"), binname);
	fclose(*fidx);
	if (*append == 0) {
	    gretl_remove(idxname);
	}
	return 1;
    }

    fprintf(stderr, "Opened database index '%s' in mode '%s'\n", 
	    idxname, imode);
    fprintf(stderr, "Opened database binary '%s' in mode '%s'\n", 
	    binname, bmode);

    return 0;
}

/* screen out any empty series, after discounting missing 
   obsservations
*/

static int *make_db_save_list (const int *list, const DATASET *dset,
			       int *err)
{
    int *dlist = gretl_list_new(list[0]);
    int i, t;

    if (dlist == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    dlist[0] = 0;

    for (i=1; i<=list[0]; i++) {
	int v = list[i];
	int gotobs = 0;

	for (t=dset->t1; t<=dset->t2; t++) {
	    if (!na(dset->Z[v][t])) {
		gotobs = 1;
		break;
	    }
	}

	if (!gotobs) {
	    continue;
	}

	dlist[0] += 1;
	dlist[dlist[0]] = v;
    }

    if (dlist[0] == 0) {
	*err = E_MISSDATA;
	free(dlist);
	dlist = NULL;
    }

    return dlist;
}

/**
 * write_db_data:
 * @fname: name of target database file (e.g. "foo.bin").
 * @list: list of series ID numbers.
 * @opt: option flag.
 * @dset: dataset struct.
 *
 * Writes the listed series from @dset to a gretl database. If @opt
 * includes OPT_F (force, overwrite), then in case any variables
 * in the database have the same names as some of those in @list,
 * replace the ones in the database.  Otherwise, in case of replicated
 * variables, set an error message and return E_DB_DUP.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int write_db_data (const char *fname, const int *list, gretlopt opt,
		   const DATASET *dset) 
{
    char idxname[FILENAME_MAX];
    char binname[FILENAME_MAX];
    FILE *fbin = NULL, *fidx = NULL;
    const int *mylist = list;
    int *dlist = NULL;
    int append = 0;
    int force = (opt & OPT_F);
    int i, err = 0;

    if (dataset_is_time_series(dset)) {
	if (dset->pd != 1 && dset->pd != 4 && dset->pd != 12) {
	    return E_PDWRONG;
	}
    } else if (dset->pd != 1) {
	return E_PDWRONG;
    }

    if (open_db_files(fname, idxname, binname, 
		      &fidx, &fbin, &append)) {
	return 1;
    }

    if (append) {
#if DB_DEBUG
	fprintf(stderr, "Appending to existing db\n");
#endif
	dlist = make_db_save_list(list, dset, &err);
	if (err) {
	    goto bailout;
	}

	if (force) {
#if DB_DEBUG
	    fprintf(stderr, "Got force flag, overwriting\n");
#endif
	    fclose(fidx);
	    fclose(fbin);
	    return append_db_data_with_replacement(idxname, binname, dlist,
						   dset);
	} else {
	    int dups = check_for_db_duplicates(dlist, dset, idxname, &err);

#if DB_DEBUG
	    fprintf(stderr, "No force flag, checking for dups\n");
#endif
	    if (err) {
		fputs("check_for_db_duplicates failed\n", stderr);
	    } else if (dups > 0) {
		gretl_errmsg_sprintf(_("Of the variables to be saved, %d were already "
				       "present in the database."), dups);
		/* FIXME add message for command line use, about the
		   --overwrite option */
		err = E_DB_DUP;
	    }
	    if (err) {
		goto bailout;
	    }
	    mylist = dlist;
	} 
    } 

    if (!append) {
	const char *s = get_optval_string(STORE, OPT_E);

	if (s != NULL && *s != '\0') {
	    fprintf(fidx, "# %s\n", s);
	} else {
	    fputs("# Description goes here\n", fidx);
	}
    }

    for (i=1; i<=mylist[0]; i++) {
	int v = mylist[i];

	output_db_var(v, dset, fidx, fbin);
    }

 bailout:

    if (fidx != NULL) fclose(fidx);
    if (fbin != NULL) fclose(fbin);

    if (dlist != NULL) {
	free(dlist);
    }

    return err;
}

#ifdef notyet

/* apparatus for writing XML databases */

static void
xml_write_tagged (const char *tag, const char *s, gzFile fz)
{
    gzprintf(fz, "<%s>\n", tag);
    gzputs(s, gz);
    gzprintf(fz, "</%s>\n", tag);
}

static void 
xml_write_attrib (const char *attr, const char *s, gzFile fz)
{
    gzprintf(fz, "%s=\"", attr);
    gzputs(s, gz);
    gzputs("\"", gz);
}

static void 
xml_write_int_attrib (const char *attr, int val, gzFile fz)
{
    gzprintf(fz, "%s=\"%d\" ", attr, val);
}

static void xml_write_db_header (const char *name, gzFile fz)
{
    double gretl_db_version = 1.0;

    gzputs("<?xml version=\"1.0\"?>\n"
	   "<!DOCTYPE gretldb SYSTEM \"gretldb.dtd\">\n\n",
	   fz);
    gzprintf(fz, "<gretldb name=\"%s\" version=\"%.1f\">\n", 
	     name, gretl_db_version);
}

static void xml_open_tag (const char *tag, gzFile fz)
{
    gzprintf(fz, "<%s ", tag);
}

static void xml_close_tag (const char *tag, gzFile fz)
{
    gzprintf(fz, "\n</%s>\n", tag);
}

static void xml_write_db_footer (gzFile fz)
{
    gzputs("\n</gretldb>\n", fz);
}

int write_gretl_xml_db (gretl_db *db, const char *fname)
{
    gzFile *fz = Z_NULL;
    int i, t, err = 0;

    fz = gretl_gzopen(fname, "wb");
    if (fz == Z_NULL) {
	return E_FOPEN;
    }

    xml_write_db_header(db->name, fz);

    if (db->source != NULL) {
	xml_write_tagged("source", db->source, fz);
    }
    
    if (db->descrip != NULL) {
	xml_write_tagged("description", db->descrip, fz);
    }

    if (db->codebook != NULL) {
	xml_write_tagged("codebook", db->codebook, fz);
    }

    for (i=0; i<db->nchaps; i++) {
	xml_open_tag("chapter", fz);
	xml_write_attrib("title", db->chapters[i]->title);
	if (db->chapters[i]->descrip != NULL) {
	    xml_write_tagged("description", db->chapters[i]->descrip);
	}
	xml_close_tag("chapter", fz);
    }

    for (i=0; i<db->nseries; i++) {
	xml_open_tag("series", fz);
	xml_write_attrib("name", db->series[i]->name);
	if (db->series[i]->chapter > 0) {
	    xml_write_int_attrib("chapter", db->series[i]->chapter);
	}
	if (db->series[i]->label != NULL) {
	    xml_write_attrib("label", db->series[i]->label);
	}
	if (db->series[i]->displayname != NULL) {
	    xml_write_attrib("displayname", db->series[i]->displayname);
	}
	if (db->series[i]->frequency > 0) {
	    xml_write_int_attrib("frequency", db->series[i]->frequency);
	}
	xml_write_attrib("startobs", db->series[i]->startobs);
	xml_write_attrib("endobs", db->series[i]->endobs);
	/* FIXME type, compact-method */

	xml_open_tag("observations", fz);
	xml_write_int_attrib("count", db->series[i]->nobs);

	if (db->series[i]->markers) {
	    xml_write_attrib("labels", "true");
	    for (t=0; i<db->series[i]->nobs; t++) {
		gzprintf(fz, "<obs label=\"%s\">%.8g</obs>\n", 
			 db->series[i]->S[t], db->series[i]->x[t]);
	    }
	} else {
	    for (t=0; i<db->series[i]->nobs; t++) {
		gzprintf(fz, "<obs>%.8g</obs>\n", db->series[i]->x[t]);
	    }
	}	    
	
	xml_close_tag("observations", fz);
	xml_close_tag("series", fz);
    }    
	
    xml_write_db_footer(fz);

    gzclose(fz);

    return err;
}

#endif


