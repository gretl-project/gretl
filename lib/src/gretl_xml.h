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

#ifndef GRETL_XML_H
#define GRETL_XML_H

#ifdef FULL_XML_HEADERS

#define XUC const xmlChar *

void gretl_xml_header (PRN *prn);

int gretl_xml_open_doc_root (const char *fname,
			     const char *rootname,
			     xmlDocPtr *pdoc,
			     xmlNodePtr *pnode);

xmlNodePtr gretl_xml_get_gfn (const char *fname,
                              xmlDocPtr *pdoc,
                              int *err);

int gretl_xml_read_buffer (const char *buf,
                           const char *rootname,
                           xmlDocPtr *pdoc,
                           xmlNodePtr *pnode);

void gretl_xml_put_int (const char *tag, int i, PRN *prn);

void gretl_xml_put_unsigned (const char *tag, unsigned int u, PRN *prn);

void gretl_xml_put_double (const char *tag, double x, PRN *prn);

void gretl_xml_put_double_array (const char *tag, double *x, int n,
				 PRN *prn);

void gretl_xml_put_strings_array (const char *tag, const char **strs,
				  int n, PRN *prn);

void gretl_xml_put_strings_array_quoted (const char *tag,
					 const char **strs, int n,
					 PRN *prn);

void gretl_xml_put_tagged_list (const char *tag, const int *list,
				PRN *prn);

int gretl_xml_put_string (const char *str, PRN *prn);

int gretl_xml_get_prop_as_int (xmlNodePtr node, const char *tag,
			       int *i);

int gretl_xml_get_prop_as_opt (xmlNodePtr node, const char *tag,
			       gretlopt *o);

unsigned int
gretl_xml_get_prop_as_unsigned_int (xmlNodePtr node,
				    const char *tag,
				    unsigned int *u);

int gretl_xml_get_prop_as_uchar (xmlNodePtr node, const char *tag,
				 unsigned char *u);

int gretl_xml_get_prop_as_double (xmlNodePtr node, const char *tag,
				  double *x);

int gretl_xml_get_prop_as_string (xmlNodePtr node, const char *tag,
				  char **pstr);

int gretl_xml_get_prop_as_bool (xmlNodePtr node, const char *tag);

GretlType gretl_xml_get_type_property (xmlNodePtr node);

int gretl_xml_node_get_int (xmlNodePtr node, xmlDocPtr doc, int *i);

int gretl_xml_node_get_unsigned (xmlNodePtr node, xmlDocPtr doc,
				 unsigned int *u);

int gretl_xml_node_get_double (xmlNodePtr node, xmlDocPtr doc,
			       double *x);

int *gretl_xml_get_list (xmlNodePtr node, xmlDocPtr doc, int *err);

int gretl_xml_node_get_string (xmlNodePtr node, xmlDocPtr doc,
			       char **pstr);

char *gretl_xml_get_string (xmlNodePtr node, xmlDocPtr doc);

int gretl_xml_node_get_trimmed_string (xmlNodePtr node, xmlDocPtr doc,
				       char **pstr);

int gretl_xml_child_get_string (xmlNodePtr node, xmlDocPtr doc,
				const char *name, char **pstr);

int *gretl_xml_get_int_array (xmlNodePtr node, xmlDocPtr doc,
			      int *nelem, int *err);

double *gretl_xml_get_double_array (xmlNodePtr node, xmlDocPtr doc,
				    int *nelem, int *err);

cmplx *gretl_xml_get_cmplx_array (xmlNodePtr node, xmlDocPtr doc,
				  int *nelem, int *err);

char **gretl_xml_get_strings_array (xmlNodePtr node, xmlDocPtr doc,
				    int *nelem, int slop, int *err);

int gretl_xml_child_get_strings_array (xmlNodePtr node, xmlDocPtr doc,
				       const char *name, char ***pstrs,
				       int *nstrs);

gretl_matrix *gretl_xml_get_matrix (xmlNodePtr node, xmlDocPtr doc,
				    int *err);

int gretl_xml_get_submask (xmlNodePtr node, xmlDocPtr doc, char **pmask);

#endif /* FULL_XML_HEADERS */

void gretl_xml_put_tagged_string (const char *tag,
				  const char *str,
				  PRN *prn);

void gretl_xml_put_tagged_string_plus (const char *tag,
				       const char *str,
				       const char *attrib,
				       const char *attval,
				       PRN *prn);

void gretl_xml_put_tagged_list (const char *tag,
				const int *list,
				PRN *prn);

void gretl_matrix_serialize (const gretl_matrix *m,
			     const char *name,
			     PRN *prn);

void gretl_list_serialize (const int *list,
			   const char *name,
			   PRN *prn);

void gretl_finite_scalar_serialize (double x,
				    const char *name,
				    PRN *prn);

void gretl_string_serialize (const char *s,
			     const char *name,
			     PRN *prn);

int gretl_matrix_write_as_gdt (const char *fname,
			       const gretl_matrix *X,
			       const char **varnames,
			       const char **labels);

int gretl_write_gdt (const char *fname, const int *list,
		     const DATASET *dset, gretlopt opt,
		     int progress);

int gretl_write_gdt_to_prn (PRN *prn, const int *list,
			    const DATASET *dset);

int gretl_read_gdt (const char *fname, DATASET *dset,
		    gretlopt opt, PRN *prn);

int gretl_read_gdt_subset (const char *fname, DATASET *dset,
			   int *vlist, gretlopt opt);

int gretl_read_gdt_varnames (const char *fname,
			     char ***vnames,
			     int *nvars);

char *gretl_get_gdt_description (const char *fname, int *err);

int load_XML_functions_file (const char *fname, gretlopt opt, PRN *prn);

char *gretl_XML_escape (const char *s);

#ifndef __GTK_DOC_IGNORE__

void gretl_xml_init (void);

void gretl_xml_cleanup (void);

#endif /* __GTK_DOC_IGNORE__ */

#endif /* GRETL_XML_H */
