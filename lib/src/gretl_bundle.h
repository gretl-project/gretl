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

#ifndef GRETL_BUNDLE_H_
#define GRETL_BUNDLE_H_

typedef struct gretl_bundle_ gretl_bundle;
typedef struct bundled_item_ bundled_item;

#define AUTO_BUNDLE "BUNDLE_RET__"

gretl_bundle *gretl_bundle_new (void);

int gretl_is_bundle (const char *name);

int type_can_be_bundled (GretlType type);

gretl_bundle *get_gretl_bundle_by_name (const char *name);

gretl_bundle *get_gretl_bundle_by_index (int idx);

void *gretl_bundle_get_content (gretl_bundle *bundle);

void *gretl_bundle_get_data (gretl_bundle *bundle, const char *key,
			     GretlType *type, int *size, int *err);

const char *gretl_bundle_get_name (gretl_bundle *bundle);

GretlType gretl_bundle_get_type (gretl_bundle *bundle, const char *key,
				 int *err);

gretl_matrix *gretl_bundle_get_matrix (gretl_bundle *bundle,
				       const char *key,
				       int *err);

double gretl_bundle_get_scalar (gretl_bundle *bundle,
				const char *key,
				int *err);

const char *gretl_bundle_get_note (gretl_bundle *bundle, const char *key);

const char *gretl_bundle_get_creator (gretl_bundle *bundle);

void *bundled_item_get_data (bundled_item *item, GretlType *type,
			     int *size);

const char *bundled_item_get_note (bundled_item *item);

int gretl_bundle_set_data (gretl_bundle *bundle, const char *key,
			   void *ptr, GretlType type, int size);

int gretl_bundle_set_note (gretl_bundle *bundle, const char *key,
			   const char *note);

int gretl_bundle_delete_data (gretl_bundle *bundle, const char *key);

void set_bundle_add_callback (void (*callback));

int gretl_bundle_add_or_replace (gretl_bundle *bundle, const char *name);

int gretl_bundle_stack_as (gretl_bundle *bundle, const char *name);

int save_named_bundle (const char *name);

int gretl_bundle_copy_as (const char *name, const char *copyname);

gretl_bundle *gretl_bundle_copy (const gretl_bundle *bundle, int *err);

int copy_bundle_arg_as (const gretl_bundle *b, const char *newname);

int gretl_bundle_delete_by_name (const char *name, PRN *prn);

int gretl_bundle_set_name (gretl_bundle *b, const char *name);

int gretl_bundle_set_creator (gretl_bundle *b, const char *name);

int gretl_bundle_print (gretl_bundle *bundle, PRN *prn);

int data_is_bundled (void *ptr);

int gretl_bundle_is_stacked (gretl_bundle *b);

int gretl_bundle_get_n_keys (gretl_bundle *b);

gretl_bundle *gretl_bundle_pull_from_stack (const char *name,
					    int *err);

int gretl_bundle_localize (const char *origname,
			   const char *localname);

int gretl_bundle_unlocalize (const char *localname,
			     const char *origname);

gretl_bundle *gretl_bundle_union (const gretl_bundle *bundle1,
				  const gretl_bundle *bundle2,
				  int *err);

void gretl_bundle_destroy (gretl_bundle *bundle);

int destroy_saved_bundles_at_level (int level);

void destroy_user_bundles (void);

int n_user_bundles (void);

void write_bundles_to_file (FILE *fp);

int load_bundle_from_xml (void *p1, void *p2, const char *name);

#endif /* GRETL_BUNDLE_H_ */



