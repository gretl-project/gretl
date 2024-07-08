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

#ifndef GRETL_ARRAY_H_
#define GRETL_ARRAY_H_

void gretl_array_destroy (gretl_array *A);

void gretl_array_void_content (gretl_array *A);

void gretl_array_nullify_content (gretl_array *A);

void gretl_array_nullify_elements (gretl_array *A);

gretl_array *gretl_array_new (GretlType type, int n, int *err);

gretl_array *gretl_array_from_strings (char **S, int n,
				       int copy, int *err);

gretl_array *gretl_matrix_array_sized (int n, int r, int c,
				       int *err);

gretl_array *gretl_singleton_array (void *ptr, GretlType atype,
				    int copy, int *err);

void *gretl_array_get_element (gretl_array *A, int i,
			       GretlType *type,
			       int *err);

int gretl_array_set_element (gretl_array *A, int i,
			     void *ptr, GretlType type,
			     int copy);

int gretl_array_delete_element (gretl_array *A, int i);

void *gretl_array_get_data (gretl_array *A, int i);

int gretl_array_set_data (gretl_array *A, int i, void *ptr);

int gretl_array_set_type (gretl_array *A, GretlType type);

void *gretl_array_get_all_data (gretl_array *A);

char **gretl_array_get_strings (gretl_array *A, int *ns);

char **gretl_array_steal_strings (gretl_array *A, int *ns);

char **gretl_array_get_stringify_strings (gretl_array *A,
					  int nreq, int *ns,
					  int *err);

char *gretl_strings_array_flatten (gretl_array *A,
                                   const char *sep,
                                   int *err);

gretl_matrix *gretl_strings_array_pos (gretl_array *A,
				       const char *s,
				       int *err);

int gretl_strings_array_includes (gretl_array *A,
				  const char *s);

int gretl_array_drop_string (gretl_array *A, const char *s);

int gretl_array_drop_null (gretl_array *A);

GretlType gretl_array_get_type (gretl_array *A);

GretlType gretl_array_get_content_type (gretl_array *A);

int gretl_array_get_length (gretl_array *A);

int gretl_array_get_next_index (gretl_array *A);

int gretl_array_set_string (gretl_array *A, int i,
			    char *s, int copy);

int gretl_array_append_string (gretl_array *A,
			       char *s,
			       int copy);

int gretl_array_set_matrix (gretl_array *A, int i,
			    gretl_matrix *m,
			    int copy);

int gretl_array_append_matrix (gretl_array *A,
			       gretl_matrix *m,
			       int copy);

int gretl_array_set_bundle (gretl_array *A, int i,
			    gretl_bundle *b,
			    int copy);

int gretl_array_append_bundle (gretl_array *A,
			       gretl_bundle *b,
			       int copy);

int gretl_array_set_array (gretl_array *A, int i,
			   gretl_array *ai,
			   int copy);

int gretl_array_append_array (gretl_array *A,
			      gretl_array *a,
			      int copy);

gretl_bundle *gretl_array_get_bundle (gretl_array *A,
				      int i);

gretl_matrix *gretl_matrix_array_flatten (gretl_array *A,
					  int vcat,
					  int *err);

gretl_array *gretl_matrix_split_by (const gretl_matrix *X,
				    const gretl_matrix *v,
				    int colwise, int chunks,
                                    int *err);

int gretl_array_set_list (gretl_array *A, int i,
			  int *L, int copy);

int gretl_array_append_list (gretl_array *A,
			     int *L, int copy);

int gretl_array_append_object (gretl_array *A,
			       void *ptr,
			       int copy);

int gretl_array_copy_into (gretl_array *A1,
			   const gretl_array *A2);

gretl_array *gretl_arrays_join (gretl_array *A,
				gretl_array *B,
				int *err);

gretl_array *gretl_arrays_union (gretl_array *A,
				 gretl_array *B,
				 int *err);

gretl_array *gretl_arrays_intersection (gretl_array *A,
					gretl_array *B,
					int *err);

gretl_array *gretl_array_copy (const gretl_array *A,
			       int *err);

gretl_array *gretl_strings_sort (const gretl_array *A,
				 int descending,
				 int *err);

int gretl_array_copy_as (const char *name, const char *copyname,
			 GretlType copytype);

gretl_array *gretl_array_copy_subspec (gretl_array *A,
				       int *list,
				       int *err);

gretl_array *get_array_by_name (const char *name);

gretl_array *get_strings_array_by_name (const char *name);

gretl_array *get_strings_array_from_series (DATASET *dset,
					    int v, int *err);

gretl_array *gretl_array_pull_from_stack (const char *name,
					  int *err);

int gretl_array_print (gretl_array *A, PRN *prn);

int gretl_array_print_range (gretl_array *A,
			     int imin, int imax,
			     PRN *prn);

void gretl_array_serialize (gretl_array *A, PRN *prn);

gretl_array *gretl_array_deserialize (void *p1, void *p2,
				      int *err);

gretl_array *gretl_matrix_col_split (const gretl_matrix *m,
				     int leadcols, int maxcols,
				     int *err);

int is_strings_array_element (const char *str,
			      char *aname,
			      char *pidx);

int gretl_array_qsort (gretl_array *a, const char *fname,
		       DATASET *set, PRN *prn);

int gretl_arrays_are_equal (const gretl_array *a,
			    const gretl_array *b);

int arglist_validate (gretl_array *keys, gretl_array *args);

#endif /* GRETL_ARRAY_H_ */
