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

#ifndef GRETL_LIST_H
#define GRETL_LIST_H

#define LISTSEP (-100)

int *gretl_list_new (int nterms);

int *gretl_list_resize (int **oldlist, int nterms);

int *gretl_list_append_term (int **plist, int v);

int *gretl_list_sort (int *list);

int gretl_list_cmp (const int *list1, const int *list2);

int *gretl_null_list (void);

int *gretl_consecutive_list_new (int lmin, int lmax);

int **gretl_list_array_new (int nlists, int nterms);

void gretl_list_array_free (int **lists, int nlists);

int *gretl_list_copy (const int *src);

int *gretl_list_copy_from_pos (const int *src, int pos);

int *gretl_list_from_string (const char *str, int *err);

int *gretl_list_from_varnames (const char *str, 
			       const DATASET *dset,
			       int *err);

int *gretl_list_from_vector (const gretl_matrix *v, 
			     const DATASET *dset,
			     int *err);

int *gretl_auxlist_from_vector (const gretl_matrix *v, 
				int *err);

int *gretl_list_from_matrix (const gretl_matrix *X,
			     const char *prefix,
			     DATASET *dset,
			     int *err);

char *gretl_list_to_numeric_string (const int *list);

char *gretl_list_to_string (const int *list, 
			    const DATASET *dset,
			    int *err);

char *gretl_list_to_compact_string (const int *list,
				    const DATASET *dset,
				    int argstyle,
				    int *err);

gretl_matrix *gretl_list_to_vector (const int *list, int *err);

char *gretl_list_to_lags_string (const int *list, int *err);

char *gretl_list_get_names (const int *list, const DATASET *dset,
			    int *err);

char **gretl_list_get_names_array (const int *list, 
				   const DATASET *dset,
				   int *err);

int in_gretl_list (const int *list, int k);

int gretl_list_delete_at_pos (int *list, int pos);

int gretl_list_purge_const (int *list, const DATASET *dset);

int gretl_list_min_max (const int *list, int *lmin, int *lmax);

int *gretl_list_add (const int *orig, const int *add, int *err);

int *gretl_list_plus (const int *l1, const int *l2, int *err);

int *gretl_list_union (const int *l1, const int *l2, int *err);

int *gretl_list_intersection (const int *l1, const int *l2, int *err);

int gretl_list_append_list (int **targ, const int *src);

int gretl_list_merge_list (int **targ, const int *src);

int *gretl_list_product (const int *X, const int *Y, 
			 DATASET *dset, int *err);

int *gretl_list_omit (const int *orig, const int *omit, int minpos, int *err);

int *gretl_list_omit_last (const int *orig, int *err);

int *gretl_list_drop (const int *orig, const int *drop, int *err);

int gretl_list_diff (int *targ, const int *biglist, const int *sublist);

int *gretl_list_diff_new (const int *biglist, const int *sublist, int minpos);

int *gretl_list_build (const char *s, const DATASET *dset, int *err);

int gretl_list_insert_list (int **targ, const int *src, int pos);

int gretl_list_insert_list_minus (int **targ, const int *src, int pos);

int *gretl_list_sublist (const int *list, int pos0, int pos1);

int *gretl_list_select (const int *list, const int *sel);

int reglist_check_for_const (int *list, const DATASET *dset);

int gretl_list_const_pos (const int *list, int minpos, 
			  const DATASET *dset);

int list_members_replaced (const MODEL *pmod, const DATASET *dset);

int gretl_list_separator_position (const int *list);

int gretl_list_has_separator (const int *list);

int gretl_list_is_consecutive (const int *list);

int gretl_list_split_on_separator (const int *list, int **plist1, int **plist2);

int *gretl_lists_join_with_separator (const int *list1, const int *list2);

int gretl_list_duplicates (const int *list, GretlCmdIndex ci);

int gretl_lists_share_members (const int *list1, const int *list2);

int gretl_list_n_distinct_members (const int *list);

int *full_var_list (const DATASET *dset, int *nvars);

int gretl_is_midas_list (const int *list, const DATASET *dset);

int gretl_list_set_midas (const int *list, DATASET *dset);

void gretl_list_print (const int *list, 
		       const DATASET *dset,
		       PRN *prn);

int *varname_match_list (const DATASET *dset, 
			 const char *pattern,
			 int *err);

int *ellipsis_list (const DATASET *dset, int v1, int v2, int *err);

int varname_match_any (const DATASET *dset, const char *pattern);

#endif /* GRETL_LIST_H */
