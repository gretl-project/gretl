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

#ifndef EDITBAR_H
#define EDITBAR_H

typedef enum {
    VIEWBAR_EDITABLE = 1 << 0,
    VIEWBAR_HAS_TEXT = 1 << 1
} ViewbarFlags;

#define winlist_item(i) (strcmp(i->icon, GRETL_STOCK_WINLIST) == 0)

#define use_toolbar_search_box(r) (r == VIEW_PKG_SAMPLE || \
				   r == VIEW_PKG_CODE || \
				   r == EDIT_PKG_CODE || \
				   r == VIEW_SCRIPT || \
				   r == SCRIPT_OUT || \
				   r == VIEW_DBSEARCH || \
				   r == X12A || r == VIEW_LOG)
extern int toolbar_icon_size;

void gretl_stock_icons_init (void);

void vwin_add_viewbar (windata_t *vwin, ViewbarFlags flags);

GtkWidget *build_text_popup (windata_t *vwin);

void gretl_stock_icons_init (void);

void gretl_tooltips_add (GtkWidget *w, const gchar *str);

GtkWidget *gretl_toolbar_new (GtkWidget *sibling);

GtkWidget *gretl_toolbar_insert (GtkWidget *tbar,
				 GretlToolItem *item,
				 GCallback func,
				 gpointer data,
				 gint pos);

GtkWidget *vwin_toolbar_insert (GretlToolItem *tool,
				GCallback func,
				GtkWidget *menu,
				windata_t *vwin,
				gint pos);

void vwin_add_tmpbar (windata_t *vwin);

#endif /* EDITBAR_H */
