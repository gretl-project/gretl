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

enum extra_cmds {
    RENAME = NC,
    RELABEL,
    SMPLDUM,
    SMPLBOOL,
    MARKERS,
    STORE_MODEL,
    VAR_SUMMARY,
    GENR_NORMAL,
    GENR_UNIFORM,
    ONLINE,
    EXPORT,
    MEANTEST2,
    MODEL_GENR,
    GR_PLOT,
    GR_XY,
    GR_IMP,
    GR_DUMMY,
    GR_BOX,
    GR_NBOX,
    PANEL,    
    COMPACT,
    CONFINT,
    COVAR,
    STAT_TABLE,
    H_TEST,
    VIEW_MODEL,
    VIEW_LOG,
    VIEW_DATA,
    VIEW_SCRIPT,
    SCRIPT_OUT,
    CONSOLE,
    EDIT_BUFFER,
    EDIT_SCRIPT,
    EDIT_NOTES,
    CLI_HELP,
    CMD_LAST
};

typedef struct {
    GtkWidget *dialog;
    GtkWidget *edit;
    GList *all_buttons;
    gpointer data;
    gint code;
} dialog_t;

typedef struct {
    GtkWidget *dialog;
    GtkWidget *listbox; 
    GtkWidget *mbar;
    GtkWidget *w;
    GtkWidget *status;
    GtkWidget *popup;
    GtkItemFactory *ifac; 
    gpointer data;
    int active_var; 
    int role;
    int id;
    char fname[MAXLEN];
} windata_t; 

typedef struct {
    char *line;
    PRN *prn;
} lineprint_t;

#endif /* GRETLTYPES_H */
