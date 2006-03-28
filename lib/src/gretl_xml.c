/*
 *  Copyright (c) by Allin Cottrell
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

#include "libgretl.h"
#include "gretl_xml.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#undef XML_DEBUG

#ifdef WIN32
# include <glib.h>
static xmlDocPtr gretl_xmlParseFile (const char *fname)
{
    xmlDocPtr ptr = NULL;
    FILE *fp = fopen(fname, "r");

    if (fp != NULL) {
	fclose(fp);
	ptr = xmlParseFile(fname);
    } else {
	int save_errno = errno;
	gchar *fconv;
	gsize wrote;

	fconv = g_locale_from_utf8(fname, -1, NULL, &wrote, NULL);
	if (fconv != NULL) {
	    ptr = xmlParseFile(fconv);
	    g_free(fconv);
	}
	errno = save_errno;
    }

    return ptr;
}
#else
# define gretl_xmlParseFile(f) xmlParseFile(f)
#endif

int gretl_xml_open_doc_root (const char *fname,
			     const char *rootname,
			     xmlDocPtr *pdoc, 
			     xmlNodePtr *pnode)
{
    xmlDocPtr doc;
    xmlNodePtr node;
    int err = 0;

    *pdoc = NULL;
    *pnode = NULL;

    doc = gretl_xmlParseFile(fname);
    if (doc == NULL) {
	sprintf(gretl_errmsg, _("xmlParseFile failed on %s"), fname);
	err = 1;
    }

    if (!err) {
	node = xmlDocGetRootElement(doc);
	if (node == NULL) {
	    sprintf(gretl_errmsg, _("%s: empty document"), fname);
	    xmlFreeDoc(doc);
	    err = 1;
	}
    }

    if (!err) {
	if (xmlStrcmp(node->name, (XUC) rootname)) {
	    sprintf(gretl_errmsg, _("File of the wrong type, root node not %s"),
		    rootname);
	    xmlFreeDoc(doc);
	    err = 1;
	}
    }    

    if (!err) {
	*pdoc = doc;
	*pnode = node;
    }

    return err;
}

static char *compact_method_to_string (int method)
{
    if (method == COMPACT_SUM) return "COMPACT_SUM";
    else if (method == COMPACT_AVG) return "COMPACT_AVG";
    else if (method == COMPACT_SOP) return "COMPACT_SOP";
    else if (method == COMPACT_EOP) return "COMPACT_EOP";
    else return "COMPACT_NONE";
}

static int compact_string_to_int (const char *str)
{
    if (!strcmp(str, "COMPACT_SUM")) return COMPACT_SUM;
    else if (!strcmp(str, "COMPACT_AVG")) return COMPACT_AVG;
    else if (!strcmp(str, "COMPACT_SOP")) return COMPACT_SOP;
    else if (!strcmp(str, "COMPACT_EOP")) return COMPACT_EOP;
    else return COMPACT_NONE;
}

/* given a full filename in @src, write to @dest a "simple"
   counterpart without leading path or extension
*/

static char *simple_fname (char *dest, const char *src)
{
    char *p;
    const char *s;

    s = strrchr(src, SLASH);

    /* take last part of src filename */
    if (s != NULL) {
        strcpy(dest, s + 1);
    } else {
        strcpy(dest, src);
    }

    /* trash any extension */
    p = strrchr(dest, '.');
    if (p != NULL && strlen(dest) > 3) {
	*p = '\0';
    }

    return dest;
}

static int alt_puts (const char *s, FILE *fp, gzFile *fz)
{
    int ret = 0;

    if (fp != NULL) {
	ret = fputs(s, fp);
    } else if (fz != NULL) {
	ret = gzputs(fz, s);
    } 

    return ret;
}

static const char *data_structure_string (int s)
{
    switch (s) {
    case TIME_SERIES:
    case SPECIAL_TIME_SERIES:
	return "time-series";
    case STACKED_TIME_SERIES:
	return "stacked-time-series";
    case STACKED_CROSS_SECTION:
	return "stacked-cross-section";
    default:
	return "cross-section";
    }
}

static int savenum (const int *list, int i)
{
    if (list != NULL) {
	return list[i];
    } else {
	return i;
    }
}

/**
 * gretl_xml_put_double:
 * @tag: name to give value.
 * @x: value to put.
 * @fp: file to which to write.
 * 
 * Writes to @fp a string of the form "%s=%.15g" if the value of
 * @x is valid, otherwise "%s=NA".
 */

void gretl_xml_put_double (const char *tag, double x, FILE *fp)
{
    if (na(x)) {
	fprintf(fp, "%s=\"NA\" ", tag);
    } else {
	fprintf(fp, "%s=\"%.15g\" ", tag, x);
    }
}

/**
 * gretl_xml_put_double_array:
 * @tag: name to give array.
 * @x: values to put.
 * @n: number of values in @x.
 * @fp: file to which to write.
 * 
 */

void gretl_xml_put_double_array (const char *tag, double *x, int n,
				 FILE *fp)
{
    int i;

    fprintf(fp, "<%s count=\"%d\">\n", tag, n);
    for (i=0; i<n; i++) {
	if (na(x[i])) {
	    fputs("NA ", fp);
	} else {
	    fprintf(fp, "%.15g ", x[i]);
	}
    }
    fprintf(fp, "</%s>\n", tag);    
}

/**
 * gretl_xml_put_strings_array:
 * @tag: name to give array.
 * @strs: array of strings to put.
 * @n: number of strings in @strs.
 * @fp: file to which to write.
 * 
 */

void gretl_xml_put_strings_array (const char *tag, const char **strs, int n,
				  FILE *fp)
{
    int i;

    fprintf(fp, "<%s count=\"%d\">\n", tag, n);
    for (i=0; i<n; i++) {
	fprintf(fp, "%s ", strs[i]);
    }
    fprintf(fp, "</%s>\n", tag); 
}

/**
 * gretl_xml_put_list:
 * @tag: name to give list.
 * @list: list of integers to be written.
 * @fp: file to which to write.
 * 
 */

void gretl_xml_put_list (const char *tag, const int *list, FILE *fp)
{
    int i;

    fprintf(fp, "<%s>\n", tag);
    for (i=0; i<=list[0]; i++) {
	fprintf(fp, "%d ", list[i]);
    }
    fprintf(fp, "</%s>\n", tag); 
}

/**
 * gretl_xml_get_prop_as_int:
 * @node: XML node pointer.
 * @tag: name by which integer property is known.
 * @targ: location to write int value.
 * 
 * Returns: 1 if an int is found and read successfully, 0
 * otherwise.
 */

int gretl_xml_get_prop_as_int (xmlNodePtr node, const char *tag,
			       int *targ)
{
    xmlChar *tmp = xmlGetProp(node, (XUC) tag);
    int ret = 0;

    if (tmp != NULL) {
	*targ = atoi((const char *) tmp);
	free(tmp);
	ret = 1;
    }

    return ret;
}

/**
 * gretl_xml_get_prop_as_uchar:
 * @node: XML node pointer.
 * @tag: name by which unsigned character property is known.
 * @targ: location to write value.
 * 
 * Returns: 1 if an unsigned char is found and read successfully, 0
 * otherwise.
 */

int gretl_xml_get_prop_as_uchar (xmlNodePtr node, const char *tag,
				 unsigned char *targ)
{
    xmlChar *tmp = xmlGetProp(node, (XUC) tag);
    int ret = 0;

    if (tmp != NULL) {
	*targ = (unsigned char) atoi((const char *) tmp);
	free(tmp);
	ret = 1;
    }

    return ret;
}

/**
 * gretl_xml_get_prop_as_double:
 * @node: XML node pointer.
 * @tag: name by which floating-point property is known.
 * @targ: location to write double value.
 * 
 * Returns: 1 if a double is found and read successfully, 0
 * otherwise.
 */

int gretl_xml_get_prop_as_double (xmlNodePtr node, const char *tag,
				  double *targ)
{
    xmlChar *tmp = xmlGetProp(node, (XUC) tag);
    int ret = 0;

    if (tmp != NULL) {
	*targ = atof((const char *) tmp);
	free(tmp);
	ret = 1;
    }

    return ret;
}

/**
 * gretl_xml_get_prop_as_string:
 * @node: XML node pointer.
 * @tag: name by which string property is known.
 * 
 * Returns: allocated string, if found, else %NULL.
 */

char *gretl_xml_get_prop_as_string (xmlNodePtr node, const char *tag)
{
    xmlChar *tmp = xmlGetProp(node, (XUC) tag);

    return (char *) tmp;
}

/**
 * gretl_xml_node_get_list:
 * @node: XML node pointer.
 * @doc: XML document pointer.
 * @err: location to receive error code.
 * 
 * Returns: allocated list read from @node, or %NULL on
 * failure.
 */

int *gretl_xml_node_get_list (xmlNodePtr node, xmlDocPtr doc, int *err)
{
    xmlChar *tmp;
    const char *p;
    int *list = NULL;
    int i, n;

    tmp = xmlNodeListGetString(doc, node->xmlChildrenNode, 1);

    if (tmp == NULL) {
	*err = E_DATA;
    } else {
	p = (const char *) tmp;
	p += strspn(p, " \r\n");
	if (sscanf(p, "%d", &n) != 1) {
	    *err = E_DATA;
	} else if (n <= 0) {
	    *err = E_DATA;
	} else {
	    p += strcspn(p, " \r\n");
	    list = gretl_list_new(n);
	    if (list == NULL) {
		*err = E_ALLOC;
	    }
	}
	if (list != NULL && !*err) {
	    for (i=1; i<=n && !*err; i++) {
		if (sscanf(p, "%d", &list[i]) != 1) {
		    *err = E_DATA;
		}
		p += strspn(p, " \r\n");
		p += strcspn(p, " \r\n");
	    }
	}
	free(tmp);
    }

    if (list != NULL && *err) {
	free(list);
	list = NULL;
    }

    return list;
}

/**
 * gretl_xml_get_doubles_array:
 * @node: XML node pointer.
 * @doc: XML document pointer.
 * @err: location to receive error code.
 * 
 * Returns: allocated array of doubles read from @node, or %NULL on
 * failure.
 */

double *gretl_xml_get_doubles_array (xmlNodePtr node, xmlDocPtr doc,
				     int *err)
{
    xmlChar *tmp = xmlGetProp(node, (XUC) "count");
    const char *p;
    double *x = NULL;
    int i, n;

    if (tmp != NULL) {
	n = atoi((const char *) tmp);
	free(tmp);
	if (n > 0) {
	    x = malloc(n * sizeof *x);
	    if (x == NULL) {
		*err = E_ALLOC;
	    } else {
		tmp = xmlNodeListGetString(doc, node->xmlChildrenNode, 1);
		if (tmp == NULL) {
		    *err = E_DATA;
		} else {
		    p = (const char *) tmp;
		    p += strspn(p, " \r\n");
		    for (i=0; i<n && !*err; i++) {
			if (sscanf(p, "%lf", &x[i]) != 1) {
			    *err = E_DATA;
			}
			/* skip to end of number */
			p += strspn(p, " \r\n");
			p += strcspn(p, " \r\n");
		    }
		    free(tmp);
		}
	    }
	}
    }

    if (x != NULL && *err) {
	free(x);
	x = NULL;
    }

    return x;
}

static const char *get_gretl_encoding (void)
{
#ifdef USE_GTK2
    const char *enc = "UTF-8";
#else
    const char *enc = get_gretl_charset();

    if (enc == NULL) {
	enc = "ISO-8859-1";
    }
#endif

    return enc;
}

void gretl_xml_header (FILE *fp)
{
    const char *enc = get_gretl_encoding();

    fprintf(fp, "<?xml version=\"1.0\" encoding=\"%s\"?>\n", enc);
}

/**
 * gretl_write_gdt:
 * @fname: name of file to write.
 * @list: list of variables to write (or %NULL to write all).
 * @Z: data matrix.
 * @pdinfo: data information struct.
 * @fmt: if %GRETL_DATA_GZIPPED write gzipped data, else uncompressed.
 * @ppaths: pointer to paths information (or NULL).
 * 
 * Write out in xml a data file containing the values of the given set
 * of variables.
 * 
 * Returns: 0 on successful completion, non-zero on error.
 */

int gretl_write_gdt (const char *fname, const int *list, 
		     const double **Z, const DATAINFO *pdinfo, 
		     GretlDataFormat fmt, PATHS *ppaths)
{
    FILE *fp = NULL;
    gzFile *fz = Z_NULL;
    const char *enc;
    int gz = (fmt == GRETL_DATA_GZIPPED);
    int tsamp = pdinfo->t2 - pdinfo->t1 + 1;
    int *pmax = NULL;
    char startdate[OBSLEN], enddate[OBSLEN];
    char datname[MAXLEN], freqstr[32];
    char numstr[128];
    char *xmlbuf = NULL;
    void *handle = NULL;
    int (*show_progress) (long, long, int) = NULL;
    long sz = 0L;
    int i, t, v, nvars;
    int err = 0;

    if (gz) {
	fz = gretl_gzopen(fname, "wb");
	if (fz == Z_NULL) err = 1;
    } else {
	fp = gretl_fopen(fname, "wb");
	if (fp == NULL) err = 1;
    }

    if (err) {
	sprintf(gretl_errmsg, _("Couldn't open %s for writing"), fname);
	return 1;
    }

    if (list != NULL) {
	nvars = list[0];
    } else {
	nvars = pdinfo->v - 1;
    }

    pmax = malloc(nvars * sizeof *pmax);
    if (pmax == NULL) {
	sprintf(gretl_errmsg, _("Out of memory"));
	err = 1;
	goto cleanup;
    } 

    sz = (tsamp * nvars * sizeof(double));
    if (sz > 100000) {
	fprintf(stderr, I_("Writing %ld Kbytes of data\n"), sz / 1024);
	if (ppaths == NULL) {
	    sz = 0L;
	}
    } else {
	sz = 0L;
    }

    if (sz) {
	show_progress = get_plugin_function("show_progress", &handle);
	if (show_progress == NULL) {
	    sz = 0L;
	}
    }

    if (sz) (*show_progress)(0, sz, SP_SAVE_INIT); 

    for (i=1; i<=nvars; i++) {
	v = savenum(list, i);
	if (pdinfo->vector[v]) {
	    pmax[i-1] = get_precision(&Z[v][pdinfo->t1], tsamp, 10);
	} else {
	    pmax[i-1] = GRETL_SCALAR_DIGITS;
	}
    }

    ntodate_full(startdate, pdinfo->t1, pdinfo);
    ntodate_full(enddate, pdinfo->t2, pdinfo);

    simple_fname(datname, fname);
    xmlbuf = gretl_xml_encode(datname);
    if (xmlbuf == NULL) {
	err = 1;
	goto cleanup;
    }

    if (custom_time_series(pdinfo)) {
	sprintf(freqstr, "special:%d", pdinfo->pd);
    } else {
	sprintf(freqstr, "%d", pdinfo->pd);
    }

    enc = get_gretl_encoding();

    if (gz) {
	gzprintf(fz, "<?xml version=\"1.0\" encoding=\"%s\"?>\n"
		 "<!DOCTYPE gretldata SYSTEM \"gretldata.dtd\">\n\n"
		 "<gretldata name=\"%s\" frequency=\"%s\" "
		 "startobs=\"%s\" endobs=\"%s\" ", 
		 enc, datname, freqstr, startdate, enddate);
    } else {
	fprintf(fp, "<?xml version=\"1.0\" encoding=\"%s\"?>\n"
		"<!DOCTYPE gretldata SYSTEM \"gretldata.dtd\">\n\n"
		"<gretldata name=\"%s\" frequency=\"%s\" "
		"startobs=\"%s\" endobs=\"%s\" ", 
		enc, datname, freqstr, startdate, enddate);
    }

    free(xmlbuf);

    if (gz) {
	gzprintf(fz, "type=\"%s\">\n", data_structure_string(pdinfo->structure));
    } else {
	fprintf(fp, "type=\"%s\">\n", data_structure_string(pdinfo->structure));
    }

    /* deal with description, if any */
    if (pdinfo->descrip != NULL) {
	xmlbuf = gretl_xml_encode(pdinfo->descrip);
	if (xmlbuf == NULL) {
	    err = 1;
	    goto cleanup;
	} else {
	    if (gz) {
		gzputs(fz, "<description>");
		gzputs(fz, xmlbuf);
		gzputs(fz, "</description>\n");
	    } else {
		fprintf(fp, "<description>%s</description>\n", xmlbuf);
	    }
	    free(xmlbuf);
#ifdef XML_DEBUG
	    fprintf(stderr, "xmlbuf encoded buffer freed\n");
#endif
	}
    }

    gretl_push_c_numeric_locale();

    /* then listing of variable names and labels */
    if (gz) {
	gzprintf(fz, "<variables count=\"%d\">\n", nvars);
    } else {
	fprintf(fp, "<variables count=\"%d\">\n", nvars);
    }

    for (i=1; i<=nvars; i++) {
	v = savenum(list, i);
	xmlbuf = gretl_xml_encode(pdinfo->varname[v]);

	if (xmlbuf == NULL) {
	    err = 1;
	    goto cleanup;
	} else {
	    if (gz) {
		gzprintf(fz, "<variable name=\"%s\"", xmlbuf);
	    } else {
		fprintf(fp, "<variable name=\"%s\"", xmlbuf);
	    }
	    free(xmlbuf);
	}

	if (!pdinfo->vector[v] && !na(Z[v][0])) {
	    if (pmax[i-1] == PMAX_NOT_AVAILABLE) {
		sprintf(numstr, "\n role=\"scalar\" value=\"%.12g\"",
			Z[v][0]);
	    } else {
		sprintf(numstr, "\n role=\"scalar\" value=\"%.*f\"",
			pmax[i-1], Z[v][0]);
	    }
	    alt_puts(numstr, fp, fz);
	}

	if (*VARLABEL(pdinfo, v)) {
	    xmlbuf = gretl_xml_encode(VARLABEL(pdinfo, v));
	    if (xmlbuf == NULL) {
		err = 1;
		goto cleanup;
	    } else {
		if (gz) {
		    gzprintf(fz, "\n label=\"%s\"", xmlbuf);
		} else {
		    fprintf(fp, "\n label=\"%s\"", xmlbuf);
		}
		free(xmlbuf);
	    }
	} 

	if (*DISPLAYNAME(pdinfo, v)) {
	    xmlbuf = gretl_xml_encode(DISPLAYNAME(pdinfo, v));
	    if (xmlbuf == NULL) {
		err = 1;
		goto cleanup;
	    } else {
		if (gz) {
		    gzprintf(fz, "\n displayname=\"%s\"", xmlbuf);
		} else {
		    fprintf(fp, "\n displayname=\"%s\"", xmlbuf);
		}
		free(xmlbuf);
	    }
	} 

	if (COMPACT_METHOD(pdinfo, v) != COMPACT_NONE) {
	    const char *meth = compact_method_to_string(COMPACT_METHOD(pdinfo, v));

	    if (gz) {
		gzprintf(fz, "\n compact-method=\"%s\"", meth);
	    } else {
		fprintf(fp, "\n compact-method=\"%s\"", meth);
	    }
	} 

	alt_puts("\n/>\n", fp, fz);
    }

    alt_puts("</variables>\n", fp, fz);

    /* then listing of observations */
    if (gz) {
	gzprintf(fz, "<observations count=\"%d\" labels=\"%s\">\n",
		tsamp, (pdinfo->markers && pdinfo->S != NULL)? "true" : "false");
    } else {
	fprintf(fp, "<observations count=\"%d\" labels=\"%s\">\n",
		tsamp, (pdinfo->markers && pdinfo->S != NULL)? "true" : "false");
    }

    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	if (pdinfo->markers && pdinfo->S != NULL) {
	    if (gz) {
		gzprintf(fz, "<obs label=\"%s\">", pdinfo->S[t]);
	    } else {
		fprintf(fp, "<obs label=\"%s\">", pdinfo->S[t]);
	    }
	} else {
	    alt_puts("<obs>", fp, fz);
	}
	for (i=1; i<=nvars; i++) {
	    v = savenum(list, i);
	    if (!pdinfo->vector[v]) {
		continue;
	    }
	    if (na(Z[v][t])) {
		strcpy(numstr, "NA ");
	    } else if (pmax[i-1] == PMAX_NOT_AVAILABLE) {
		sprintf(numstr, "%.12g ", Z[v][t]);
	    } else {
		sprintf(numstr, "%.*f ", pmax[i-1], Z[v][t]);
	    }
	    alt_puts(numstr, fp, fz);
	}

	alt_puts("</obs>\n", fp, fz);

	if (sz && t && ((t - pdinfo->t1) % 50 == 0)) { 
	    (*show_progress) (50, tsamp, SP_NONE);
	}
    }

    alt_puts("</observations>\n</gretldata>\n", fp, fz);

 cleanup: 

    gretl_pop_c_numeric_locale();

    if (sz) {
	(*show_progress)(0, pdinfo->t2 - pdinfo->t1 + 1, SP_FINISH);
	close_plugin(handle);
    } 

    if (pmax) free(pmax);
    if (fp != NULL) fclose(fp);
    if (fz != Z_NULL) gzclose(fz);

    return err;
}

static int transcribe_string (char *targ, const char *src, int maxlen,
			      int convert)
{
    *targ = '\0';

#ifndef USE_GTK2
    if (convert) {
	char tmp[128] = {0};

	if (maxlen > 128) {
	    maxlen = 128;
	}
	strncat(tmp, src, maxlen - 1);
	utf8_to_iso_latin_1((unsigned char *) targ, maxlen, 
			    (unsigned char *) tmp, maxlen);
    } else {
	strncat(targ, src, maxlen - 1);
    }
#else
    strncat(targ, src, maxlen - 1);
#endif

    return 0;
}

static int process_varlist (xmlNodePtr node, DATAINFO *pdinfo, double ***pZ,
			    int to_iso_latin)
{
    xmlNodePtr cur;
    xmlChar *tmp = xmlGetProp(node, (XUC) "count");
    int i, err = 0;

    if (tmp != NULL) {
	int v;

	if (sscanf((char *) tmp, "%d", &v) == 1) {
	    pdinfo->v = v + 1;
	} else {
	    sprintf(gretl_errmsg, _("Failed to parse count of variables"));
	    err = 1;
	}
	if (!err && dataset_allocate_varnames(pdinfo)) {
	    sprintf(gretl_errmsg, _("Out of memory reading data file"));
	    err = 1;
	}
	if (!err) {
	    *pZ = malloc(pdinfo->v * sizeof **pZ);
	    if (*pZ == NULL) {
		sprintf(gretl_errmsg, _("Out of memory reading data file"));
		err = 1;
	    }
	}		
	free(tmp);
    } else {
	sprintf(gretl_errmsg, _("Got no variables"));
	err = 1;
    }

    if (err) return 1;

    /* now get individual variable info: names and labels */
    cur = node->xmlChildrenNode;
    while (cur && xmlIsBlankNode(cur)) {
	cur = cur->next;
    }

    if (cur == 0) {
	sprintf(gretl_errmsg, _("Got no variables"));
	return 1;
    }

    i = 1;
    while (cur != NULL) {
        if (!xmlStrcmp(cur->name, (XUC) "variable")) {
	    tmp = xmlGetProp(cur, (XUC) "name");
	    if (tmp != NULL) {
		transcribe_string(pdinfo->varname[i], (char *) tmp, VNAMELEN,
				  to_iso_latin); 
		free(tmp);
	    } else {
		sprintf(gretl_errmsg, _("Variable %d has no name"), i);
		return 1;
	    }
	    tmp = xmlGetProp(cur, (XUC) "label");
	    if (tmp != NULL) {
		transcribe_string(VARLABEL(pdinfo, i), (char *) tmp, MAXLABEL,
				  to_iso_latin);
		free(tmp);
	    }
	    tmp = xmlGetProp(cur, (XUC) "displayname");
	    if (tmp != NULL) {
		transcribe_string(DISPLAYNAME(pdinfo, i), (char *) tmp, MAXDISP,
				  to_iso_latin);
		free(tmp);
	    }
	    tmp = xmlGetProp(cur, (XUC) "compact-method");
	    if (tmp != NULL) {
		COMPACT_METHOD(pdinfo, i) = compact_string_to_int((char *) tmp);
		free(tmp);
	    }
	    tmp = xmlGetProp(cur, (XUC) "role");
	    if (tmp != NULL) {
		if (!strcmp((char *) tmp, "scalar")) {
		    char *val = (char *) xmlGetProp(cur, (XUC) "value");
		    
		    if (val) {
			double xx = atof(val);

			free(val);
			(*pZ)[i] = malloc(sizeof ***pZ);
			(*pZ)[i][0] = xx;
			pdinfo->vector[i] = 0;
		    }
		}
		free(tmp);
	    }
	    i++;
	}	    
	cur = cur->next;
    }
   
    if (i != pdinfo->v) {
	sprintf(gretl_errmsg, _("Number of variables does not match declaration"));
	err = 1;
    } 

    return err;
}

static int process_values (double **Z, DATAINFO *pdinfo, int t, char *s)
{
    char valstr[32];
    double x;
    int i, err = 0;

    *gretl_errmsg = '\0';

    for (i=1; i<pdinfo->v && !err; i++) {
	if (!pdinfo->vector[i]) {
	    continue;
	}
	s = strpbrk(s, "01234567890+-NA");
	if (s == NULL) {
	    fprintf(stderr, "i = %d: s == NULL in process_values()\n", i);
	    err = 1;
	} else {
	    if (*s == '\0' || sscanf(s, "%31s", valstr) != 1) {
		fputs("s is blank in process_values()\n", stderr);
		err = 1;
	    } else {
		if (!strcmp(valstr, "NA")) {
		    x = NADBL;
		} else if (check_atof(valstr)) {
		    err = 1;
		} else {
		    sscanf(valstr, "%lf", &x);
		}
	    }
	}
	if (!err) {
	    if (t < pdinfo->n) {
		Z[i][t] = x;
	    }
	    s = strpbrk(s, " \t\n\r");
	}
    }

    if (err && *gretl_errmsg == '\0') {
	sprintf(gretl_errmsg, _("Failed to parse data values at obs %d"), t+1);
    }

    return err;
}

static int process_observations (xmlDocPtr doc, xmlNodePtr node, 
				 double ***pZ, DATAINFO *pdinfo,
				 long progress, int to_iso_latin)
{
    xmlNodePtr cur;
    xmlChar *tmp;
    int n, i, t;
    void *handle;
    int (*show_progress) (long, long, int) = NULL;

    tmp = xmlGetProp(node, (XUC) "count");
    if (tmp == NULL) {
	return 1;
    } 

    if (sscanf((char *) tmp, "%d", &n) == 1) {
	pdinfo->n = n;
	free(tmp);
    } else {
	sprintf(gretl_errmsg, _("Failed to parse number of observations"));
	free(tmp);
	return 1;
    }

    if (progress > 0) {
	show_progress = get_plugin_function("show_progress", &handle);
	if (show_progress == NULL) {
	    progress = 0L;
	}
    }

    tmp = xmlGetProp(node, (XUC) "labels");
    if (tmp) {
	if (!strcmp((char *) tmp, "true")) {
	    if (dataset_allocate_obs_markers(pdinfo)) {
		sprintf(gretl_errmsg, "Out of memory");
		return 1;
	    }
	} else if (strcmp((char *) tmp, "false")) {
	    sprintf(gretl_errmsg, _("labels attribute for observations must be "
		    "'true' or 'false'"));
	    return 1;
	}
	free(tmp);
    } else {
	return 1;
    }

    if (pdinfo->endobs[0] == '\0') {
	sprintf(pdinfo->endobs, "%d", pdinfo->n);
    }

    pdinfo->t2 = pdinfo->n - 1;

    for (i=0; i<pdinfo->v; i++) {
	if (!pdinfo->vector[i]) {
	    continue;
	}
	(*pZ)[i] = malloc(pdinfo->n * sizeof ***pZ);
	if ((*pZ)[i] == NULL) {
	    return 1;
	}
    }

    for (t=0; t<pdinfo->n; t++) {
	(*pZ)[0][t] = 1.0;
    }

    /* now get individual obs info: labels and values */
    cur = node->xmlChildrenNode;
    while (cur && xmlIsBlankNode(cur)) {
	cur = cur->next;
    }

    if (cur == NULL) {
	sprintf(gretl_errmsg, _("Got no observations\n"));
	return 1;
    }

    if (progress) {
	(*show_progress)(0L, progress, SP_LOAD_INIT);
    }

    t = 0;
    while (cur != NULL) {
        if (!xmlStrcmp(cur->name, (XUC) "obs")) {
	    if (pdinfo->markers) {
		tmp = xmlGetProp(cur, (XUC) "label");
		if (tmp) {
		    transcribe_string(pdinfo->S[t], (char *) tmp, OBSLEN,
				      to_iso_latin); 
		    free(tmp);
		} else {
		    sprintf(gretl_errmsg, _("Case marker missing at obs %d"), t+1);
		    return 1;
		}
	    }

	    tmp = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);

	    if (tmp) {
		if (process_values(*pZ, pdinfo, t, (char *) tmp)) {
		    return 1;
		}
		free(tmp);
		t++;
	    } else {
		sprintf(gretl_errmsg, _("Values missing at observation %d"), t+1);
		return 1;
	    }
	}	    
	cur = cur->next;
	if (progress && t > 0 && t % 50 == 0) {
	    (*show_progress) (50L, (long) pdinfo->n, SP_NONE);
	}
    }

    if (progress) {
	(*show_progress)(0L, (long) pdinfo->n, SP_FINISH);
	close_plugin(handle);
    }

    if (t != pdinfo->n) {
	sprintf(gretl_errmsg, _("Number of observations does not match declaration"));
	return 1;
    }

    else return 0;
}

static long get_filesize (const char *fname)
{
    struct stat buf;

    if (stat(fname, &buf) == 0) {
        return buf.st_size;
    } else {
        return -1;
    }
}

static int xml_get_data_structure (xmlNodePtr node, int *dtype)
{
    xmlChar *tmp = xmlGetProp(node, (XUC) "type");
    int err = 0;

    if (tmp == NULL) {
	sprintf(gretl_errmsg, 
		_("Required attribute 'type' is missing from data file"));
	err = 1;
    } else {
	if (!strcmp((char *) tmp, "cross-section")) {
	    *dtype = CROSS_SECTION;
	} else if (!strcmp((char *) tmp, "time-series")) {
	    *dtype = TIME_SERIES;
	} else if (!strcmp((char *) tmp, "stacked-time-series")) {
	    *dtype = STACKED_TIME_SERIES;
	} else if (!strcmp((char *) tmp, "stacked-cross-section")) {
	    *dtype = STACKED_CROSS_SECTION;
	} else {
	    sprintf(gretl_errmsg, _("Unrecognized type attribute for data file"));
	    err = 1;
	}
	free(tmp);
    }

    return err;
}

static int xml_get_data_frequency (xmlNodePtr node, int *pd, int *dtype)
{
    xmlChar *tmp = xmlGetProp(node, (XUC) "frequency");
    int err = 0;

    *pd = 1;

    if (tmp != NULL) {
	if (!strncmp((char *) tmp, "special", 7)) {
	    *dtype = SPECIAL_TIME_SERIES;
	    if (sscanf((char *) tmp + 7, ":%d", pd) == 1) {
		fprintf(stderr, "custom time series, frequency %d\n", *pd);
	    } else {
		fprintf(stderr, "custom time series, using frequency 1\n");
	    }
	} else if (sscanf((char *) tmp, "%d", pd) != 1) {
	    strcpy(gretl_errmsg, _("Failed to parse data frequency"));
	    err = 1;
	}
	free(tmp);
    }

    return err;
}

static int xml_get_startobs (xmlNodePtr node, double *sd0, char *stobs,
			     int caldata)
{
    xmlChar *tmp = xmlGetProp(node, (XUC) "startobs");
    int err = 0;

    if (tmp != NULL) {
	char obstr[16];

	obstr[0] = '\0';
	strncat(obstr, (char *) tmp, 15);
	charsub(obstr, ':', '.');
	
	if (strchr(obstr, '/') != NULL && caldata) {
	    long ed = get_epoch_day((char *) tmp);

	    if (ed < 0) {
		err = 1;
	    } else {
		*sd0 = ed;
	    }
	} else {
	    double x;

	    if (sscanf(obstr, "%lf", &x) != 1) {
		err = 1;
	    } else {
		*sd0 = x;
	    }
	}

	if (err) {
	    strcpy(gretl_errmsg, _("Failed to parse startobs"));
	} else {
	    stobs[0] = '\0';
	    strncat(stobs, (char *) tmp, OBSLEN - 1);
	    colonize_obs(stobs);
	}

	free(tmp);
    }

    return err;
}

static int xml_get_endobs (xmlNodePtr node, char *endobs, int caldata)
{
    xmlChar *tmp = xmlGetProp(node, (XUC) "endobs");
    int err = 0;

    if (tmp != NULL) {
	if (caldata) {
	    long ed = get_epoch_day((char *) tmp);

	    if (ed < 0) err = 1;
	} else {
	    double x;

	    if (sscanf((char *) tmp, "%lf", &x) != 1) {
		err = 1;
	    }
	} 

	if (err) {
	    strcpy(gretl_errmsg, _("Failed to parse endobs"));
	} else {
	    endobs[0] = '\0';
	    strncat(endobs, (char *) tmp, OBSLEN - 1);
	    colonize_obs(endobs);
	}

	free(tmp);
    }

    return err;
}

static void data_read_message (const char *fname, DATAINFO *pdinfo, PRN *prn)
{
    pprintf(prn, M_("\nRead datafile %s\n"), fname);
    pprintf(prn, M_("periodicity: %d, maxobs: %d,\n"
		    "observations range: %s-%s\n"), 
	    (custom_time_series(pdinfo))? 1 : pdinfo->pd, 
	    pdinfo->n, pdinfo->stobs, pdinfo->endobs);
    pputc(prn, '\n');
}

/**
 * gretl_read_gdt:
 * @pZ: pointer to data set.
 * @ppdinfo: pointer to data information struct.
 * @fname: name of file to try.
 * @ppaths: path information struct.
 * @data_status: DATA_NONE: no datafile currently open; DATA_CLEAR: datafile
 * is open, should be cleared; DATA_APPEND: add to current dataset.
 * @prn: where messages should be written.
 * @gui: should = 1 if the function is launched from the GUI, else 0.
 * 
 * Read data from file into gretl's work space, allocating space as
 * required.
 * 
 * Returns: 0 on successful completion, non-zero otherwise.
 */

int gretl_read_gdt (double ***pZ, DATAINFO **ppdinfo, char *fname,
		    PATHS *ppaths, int data_status, PRN *prn, int gui) 
{
    DATAINFO *tmpdinfo;
    double **tmpZ = NULL;
    xmlDocPtr doc = NULL;
    xmlNodePtr cur;
    int gotvars = 0, gotobs = 0, err = 0;
    int caldata = 0;
    int to_iso_latin = 0;
    long fsz, progress = 0L;

    *gretl_errmsg = '\0';

    check_for_console(prn);

    tmpdinfo = datainfo_new();
    if (tmpdinfo == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    /* COMPAT: Do not generate nodes for formatting spaces */
    LIBXML_TEST_VERSION
	xmlKeepBlanksDefault(0);

    fsz = get_filesize(fname);
    if (fsz > 100000) {
	fprintf(stderr, "%s %ld bytes %s...\n", 
		(is_gzipped(fname))? I_("Uncompressing") : I_("Reading"),
		fsz, I_("of data"));
	if (gui) progress = fsz;
    }

    doc = gretl_xmlParseFile(fname);
    if (doc == NULL) {
	sprintf(gretl_errmsg, _("xmlParseFile failed on %s"), fname);
	err = 1;
	goto bailout;
    }

#ifndef USE_GTK2
    if (doc->encoding != NULL && strstr((char *) doc->encoding, "UTF")) {
	to_iso_latin = 1;
    }
#endif

    cur = xmlDocGetRootElement(doc);
    if (cur == NULL) {
        sprintf(gretl_errmsg, _("%s: empty document"), fname);
	err = 1;
	goto bailout;
    }

    if (xmlStrcmp(cur->name, (XUC) "gretldata")) {
        sprintf(gretl_errmsg, _("File of the wrong type, root node not gretldata"));
	err = 1;
	goto bailout;
    }

    /* set some datainfo parameters */

    err = xml_get_data_structure(cur, &tmpdinfo->structure);
    if (err) {
	goto bailout;
    } 

    err = xml_get_data_frequency(cur, &tmpdinfo->pd, &tmpdinfo->structure);
    if (err) {
	goto bailout;
    }   

    gretl_push_c_numeric_locale();

    strcpy(tmpdinfo->stobs, "1");
    caldata = dataset_is_daily(tmpdinfo) || dataset_is_weekly(tmpdinfo);

    err = xml_get_startobs(cur, &tmpdinfo->sd0, tmpdinfo->stobs, caldata);
    if (err) {
	gretl_pop_c_numeric_locale();
	goto bailout;
    }     

    *tmpdinfo->endobs = '\0';
    caldata = calendar_data(tmpdinfo);

    err = xml_get_endobs(cur, tmpdinfo->endobs, caldata);
    if (err) {
	gretl_pop_c_numeric_locale();
	goto bailout;
    }     

    /* Now walk the tree */
    cur = cur->xmlChildrenNode;
    while (cur != NULL && !err) {
        if (!xmlStrcmp(cur->name, (XUC) "description")) {
	    tmpdinfo->descrip = (char *) 
		xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
        } else if (!xmlStrcmp(cur->name, (XUC) "variables")) {
	    if (process_varlist(cur, tmpdinfo, &tmpZ, to_iso_latin)) {
		err = 1;
	    } else {
		gotvars = 1;
	    }
	} else if (!xmlStrcmp(cur->name, (XUC) "observations")) {
	    if (!gotvars) {
		sprintf(gretl_errmsg, _("Variables information is missing"));
		err = 1;
	    }
	    if (process_observations(doc, cur, &tmpZ, tmpdinfo, 
				     progress, to_iso_latin)) {
		err = 1;
	    } else {
		gotobs = 1;
	    }
	}
	if (!err) cur = cur->next;
    }

    gretl_pop_c_numeric_locale();

    if (err) {
	goto bailout;
    }

    if (!gotvars) {
	sprintf(gretl_errmsg, _("Variables information is missing"));
	err = 1;
	goto bailout;
    }

    if (!gotobs) {
	sprintf(gretl_errmsg, _("No observations were found"));
	err = 1;
	goto bailout;
    }

    if (ppaths != NULL && fname != ppaths->datfile) {
	strcpy(ppaths->datfile, fname);
    }

    data_read_message(fname, tmpdinfo, prn);

    if (data_status == DATA_APPEND) {
	err = merge_data(pZ, *ppdinfo, tmpZ, tmpdinfo, prn);
	if (err) {
	    tmpZ = NULL;
	    free(tmpdinfo);
	    tmpdinfo = NULL;
	}
    } else {
	free_Z(*pZ, *ppdinfo);
	if (data_status == DATA_CLEAR) {
	    clear_datainfo(*ppdinfo, CLEAR_FULL);
	}
	free(*ppdinfo);
	*ppdinfo = tmpdinfo;
	*pZ = tmpZ;
    }

 bailout:

    if (doc != NULL) {
	xmlFreeDoc(doc);
	xmlCleanupParser();
    }

    if (err) {
	free_Z(tmpZ, tmpdinfo);
	clear_datainfo(tmpdinfo, CLEAR_FULL);
	free(tmpdinfo);
    }

    console_off();

    return err;
}

/**
 * gretl_get_gdt_description:
 * @fname: name of file to try.
 * 
 * Read data description for gretl xml (.gdt) data file.
 * 
 * Returns: buffer containing description, or NULL on failure.
 */

char *gretl_get_gdt_description (const char *fname)
{
    xmlDocPtr doc;
    xmlNodePtr cur;
    xmlChar *buf = NULL;

    *gretl_errmsg = '\0';

    LIBXML_TEST_VERSION
	xmlKeepBlanksDefault(0);

    doc = gretl_xmlParseFile(fname);
    if (doc == NULL) {
	sprintf(gretl_errmsg, _("xmlParseFile failed on %s"), fname);
	return NULL;
    }

    cur = xmlDocGetRootElement(doc);
    if (cur == NULL) {
        sprintf(gretl_errmsg, _("%s: empty document"), fname);
	xmlFreeDoc(doc);
	return NULL;
    }

    if (xmlStrcmp(cur->name, (XUC) "gretldata")) {
        sprintf(gretl_errmsg, _("File of the wrong type, root node not gretldata"));
	xmlFreeDoc(doc);
	return NULL;
    }

    cur = cur->xmlChildrenNode;
    while (cur != NULL) {
        if (!xmlStrcmp(cur->name, (XUC) "description")) {
	    buf = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
	    break;
        }
	cur = cur->next;
    }

    xmlFreeDoc(doc);
    xmlCleanupParser();

    return (char *) buf;
}

