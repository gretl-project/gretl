/*
 *   Copyright (c) by Allin Cottrell
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
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

#ifndef GRETLTYPES_H
#define GRETLTYPES_H

typedef struct {
    GtkWidget *dialog;
    GtkWidget *edit;
    GList *all_buttons;
    gpointer data;
    gint code;
} dialog_t;

typedef struct {
    GtkWidget *dialog;
    GtkWidget *vbox;
    GtkWidget *listbox; 
    GtkWidget *mbar;
    GtkWidget *w;
    GtkWidget *status;
    GtkWidget *popup;
    GtkItemFactory *ifac; 
    gpointer data;
    int active_var; 
    int role;
    int help_active;
    char fname[MAXLEN];
} windata_t; 

#endif /* GRETLTYPES_H */
