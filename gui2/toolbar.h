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

#ifndef TOOLBAR_H
#define TOOLBAR_H

typedef enum {
    VIEWBAR_EDITABLE = 1 << 0,
    VIEWBAR_HAS_TEXT = 1 << 1
} ViewbarFlags;

#define winlist_item(i) (strcmp(i->icon, GRETL_STOCK_WINLIST) == 0)

void gretl_stock_icons_init (void);

void add_mainwin_toolbar (GtkWidget *vbox);

void gretl_pdf_manual (void);

void vwin_add_viewbar (windata_t *vwin, ViewbarFlags flags);

void viewbar_add_edit_items (windata_t *vwin);

GtkWidget *build_text_popup (windata_t *vwin);

void gretl_stock_icons_init (void);

void gretl_tooltips_add (GtkWidget *w, const gchar *str);

GtkWidget *gretl_toolbar_new (void);

GtkWidget *gretl_toolbar_insert (GtkWidget *tbar,
				 GretlToolItem *item,
				 GCallback func,
				 gpointer data,
				 gint pos);

void vwin_toolbar_insert_winlist (windata_t *vwin);

#endif /* TOOLBAR_H */
