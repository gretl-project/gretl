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

int *gretl_list_new (int nterms);

int *gretl_list_resize (int **oldlist, int nterms);

int *gretl_list_append_term (int **plist, int v);

int *gretl_list_sort (int *list);

int gretl_list_cmp (const int *list1, const int *list2);

int *gretl_null_list (void);

int *gretl_consecutive_list_new (int lmin, int lmax);

int *gretl_list_copy (const int *src);

int *gretl_list_copy_from_pos (const int *src, int pos);

int *gretl_list_from_string (const char *str, int *err);

char *gretl_list_to_string (const int *list);

char *gretl_list_to_lags_string (const int *list, int *err);

int in_gretl_list (const int *list, int k);

int gretl_list_delete_at_pos (int *list, int pos);

int gretl_list_purge_const (int *list, const double **Z,
			    const DATAINFO *pdinfo);

int *gretl_list_add (const int *orig, const int *add, int *err);

int *gretl_list_union (const int *l1, const int *l2, int *err);

int *gretl_list_intersection (const int *l1, const int *l2, int *err);

int *gretl_list_omit (const int *orig, const int *omit, int minpos, int *err);

int *gretl_list_omit_last (const int *orig, int *err);

int *gretl_list_drop (const int *orig, const int *drop, int *err);

int gretl_list_diff (int *targ, const int *biglist, const int *sublist);

int *gretl_list_diff_new (const int *biglist, const int *sublist, int minpos);

int *gretl_list_build (const char *s, const DATAINFO *pdinfo, int *err);

int gretl_list_add_list (int **targ, const int *src);

int gretl_list_insert_list (int **targ, const int *src, int pos);

int gretl_list_insert_list_minus (int **targ, const int *src, int pos);

int reglist_check_for_const (int *list, const double **Z,
			     const DATAINFO *pdinfo);

int gretl_list_const_pos (const int *list, int minpos, const double **Z,
			  const DATAINFO *pdinfo);

int list_members_replaced (const int *list, const DATAINFO *pdinfo,
			   int ref_id);

int gretl_list_position (int v, const int *list);

int gretl_list_separator_position (const int *list);

int gretl_list_has_separator (const int *list);

int gretl_list_is_consecutive (const int *list);

int gretl_list_split_on_separator (const int *list, int **plist1, int **plist2);

int *gretl_lists_join_with_separator (const int *list1, const int *list2);

int gretl_list_duplicates (const int *list, GretlCmdIndex ci);

int gretl_lists_share_members (const int *list1, const int *list2);

int *full_var_list (const DATAINFO *pdinfo, int *nvars);

int n_saved_lists (void);

int max_varno_in_saved_lists (void);

const char *get_list_name_by_index (int idx);

int *get_list_by_name (const char *name);

int append_to_list_by_name (const char *targ, const int *add);

int subtract_from_list_by_name (const char *targ, const int *sub);

int replace_list_by_name (const char *targ, const int *src);

int remember_list (const int *list, const char *name, PRN *prn);

int copy_named_list_as (const char *orig, const char *newname);

int named_list_lower_level (const char *name);

int rename_saved_list (const char *orig, const char *newname); 

int create_named_null_list (const char *name);

int delete_list_by_name (const char *name);

int destroy_saved_lists_at_level (int level);

int gretl_lists_revise (const int *dlist, int dmin);

void gretl_lists_cleanup (void);

int load_user_lists_file (const char *fname);

int gretl_serialize_lists (const char *fname);

void gretl_list_print (const char *lname, 
		       const DATAINFO *pdinfo,
		       PRN *prn);

int *varname_match_list (const DATAINFO *pdinfo, const char *pattern);

int varname_match_any (const DATAINFO *pdinfo, const char *pattern);

#endif /* GRETL_LIST_H */
