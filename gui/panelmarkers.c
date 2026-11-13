/*
 *  gretl -- Gnu Regression, Econometrics and Time-series Library
 *  Copyright (C) 2001 Allin Cottrell and Riccardo "Jack" Lucchetti
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 */

/* GUI access to "genr markers = ..." for panel data: let the user
   pick a unit-identifier series (string-valued, or integer-valued)
   and a period-identifier series (integer-valued, e.g. year), show
   a live preview of the resulting markers, and then execute

     genr markers = "%s:%d", <unit>, <period>

   (or "%d:%d" if the unit series is numeric).
*/

#include "gretl.h"
#include "dlgutils.h"
#include "dialogs.h"
#include "gui_utils.h"
#include "library.h"
#include "menustate.h"
#include "panelmarkers.h"

typedef struct pm_info_ {
    GtkWidget *dlg;
    GtkWidget *ucombo;   /* unit selector */
    GtkWidget *tcombo;   /* period selector */
    GtkWidget *preview;  /* label showing sample markers */
    GtkWidget *ok;
    int *ulist;          /* candidate unit series IDs */
    int *tlist;          /* candidate period series IDs */
} pm_info;

#define PM_NPREVIEW 4   /* number of leading rows to show */

/* Non-negative integer values on every observation, and not
   constant: the same criterion as panelvars_list_ok() in datawiz.c.
   Always applied over the whole dataset, since markers apply to all
   observations whatever the current sample.
*/

static int int_index_candidate (int v)
{
    const double *x = dataset->Z[v];
    int t;

    for (t=0; t<dataset->n; t++) {
	if (na(x[t]) || x[t] < 0 || x[t] != floor(x[t])) {
	    return 0;
	}
    }

    return !gretl_isconst(0, dataset->n - 1, x);
}

/* Build the candidate lists: units may be string-valued series
   (listed first, as the most plausible) or integer-valued ones;
   periods must be integer-valued and not string-valued.
   Returns 0 if we can't offer both roles.
*/

static int get_candidates (pm_info *pmi)
{
    int *ul = gretl_list_new(dataset->v);
    int *tl = gretl_list_new(dataset->v);
    int i, pass;

    if (ul == NULL || tl == NULL) {
	free(ul); free(tl);
	return 0;
    }

    ul[0] = tl[0] = 0;

    for (pass=0; pass<2; pass++) {
	for (i=1; i<dataset->v; i++) {
	    int sv = is_string_valued(dataset, i);

	    if (pass == 0 && sv) {
		if (!gretl_isconst(0, dataset->n - 1, dataset->Z[i])) {
		    ul[++ul[0]] = i;
		}
	    } else if (pass == 1 && !sv && int_index_candidate(i)) {
		ul[++ul[0]] = i;
		tl[++tl[0]] = i;
	    }
	}
    }

    pmi->ulist = ul;
    pmi->tlist = tl;

    return ul[0] > 0 && tl[0] > 0 && (ul[0] + tl[0] > 2 || ul[1] != tl[1]);
}

/* write the marker for row @t into @targ (size OBSLEN); returns 1 if
   the result had to be truncated */

static int make_marker (int u, int p, int t, char *targ)
{
    char full[2 * VNAMELEN + 32];

    if (is_string_valued(dataset, u)) {
	const char *s = series_get_string_for_obs(dataset, u, t);

	sprintf(full, "%s:%d", s != NULL ? s : "", (int) dataset->Z[p][t]);
    } else {
	sprintf(full, "%d:%d", (int) dataset->Z[u][t], (int) dataset->Z[p][t]);
    }

    *targ = '\0';
    strncat(targ, full, OBSLEN - 1);

    return strlen(full) >= OBSLEN;
}

static int combo_series_id (GtkWidget *combo, const int *list)
{
    int k = gtk_combo_box_get_active(GTK_COMBO_BOX(combo));

    return k >= 0 ? list[k+1] : 0;
}

/* Refresh the preview and the OK button: warn about markers that
   would be truncated (OBSLEN) or non-unique.
*/

static void pm_update_preview (GtkWidget *w, pm_info *pmi)
{
    int u = combo_series_id(pmi->ucombo, pmi->ulist);
    int p = combo_series_id(pmi->tcombo, pmi->tlist);
    GString *gs = g_string_new(NULL);
    GHashTable *seen = NULL;
    char mk[OBSLEN];
    int ndup = 0, ntrunc = 0;
    int t, ok = (u > 0 && p > 0 && u != p);

    if (!ok) {
	g_string_append(gs, _("Please select two different series"));
    } else {
	seen = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, NULL);
	for (t=0; t<dataset->n; t++) {
	    ntrunc += make_marker(u, p, t, mk);
	    if (!g_hash_table_add(seen, g_strdup(mk))) {
		ndup++;
	    }
	    if (t < PM_NPREVIEW) {
		g_string_append_printf(gs, "%s\n", mk);
	    } else if (t == PM_NPREVIEW) {
		g_string_append(gs, "...\n");
	    }
	}
	if (ntrunc > 0) {
	    g_string_append_printf(gs, _("\nWarning: %d markers are too long "
					 "and will be truncated"), ntrunc);
	}
	if (ndup > 0) {
	    g_string_append_printf(gs, _("\nWarning: %d markers are not unique "
					 "(these series may not identify "
					 "unit and period)"), ndup);
	}
	g_hash_table_destroy(seen);
    }

    gtk_label_set_text(GTK_LABEL(pmi->preview), gs->str);
    gtk_widget_set_sensitive(pmi->ok, ok);
    g_string_free(gs, TRUE);
}

static void pm_apply (GtkWidget *w, pm_info *pmi)
{
    int u = combo_series_id(pmi->ucombo, pmi->ulist);
    int p = combo_series_id(pmi->tcombo, pmi->tlist);
    const char *fmt = is_string_valued(dataset, u) ? "%s:%d" : "%d:%d";
    gchar *line;
    PRN *prn;
    int err;

    if (u <= 0 || p <= 0 || u == p || bufopen(&prn)) {
	return;
    }

    line = g_strdup_printf("markers = \"%s\", %s, %s", fmt,
			   dataset->varname[u], dataset->varname[p]);
    err = gui_run_genr(line, dataset, OPT_NONE, prn);

    if (err) {
	errmsg_plus(err, gretl_print_get_buffer(prn));
    } else {
	lib_command_sprintf("genr %s", line);
	record_command_verbatim();
	mark_dataset_as_modified();
	refresh_data();
	gtk_widget_destroy(pmi->dlg);
    }

    g_free(line);
    gretl_print_destroy(prn);
}

static void pm_info_free (GtkWidget *w, pm_info *pmi)
{
    free(pmi->ulist);
    free(pmi->tlist);
    free(pmi);
}

static GtkWidget *pm_combo (GtkWidget *table, int row, const char *label,
			    const int *list, pm_info *pmi)
{
    GtkWidget *l = gtk_label_new(label);
    GtkWidget *c = gtk_combo_box_text_new();
    int i;

    gtk_misc_set_alignment(GTK_MISC(l), 1.0, 0.5);
    gtk_table_attach_defaults(GTK_TABLE(table), l, 0, 1, row, row+1);
    gtk_table_attach_defaults(GTK_TABLE(table), c, 1, 2, row, row+1);

    for (i=1; i<=list[0]; i++) {
	combo_box_append_text(c, dataset->varname[list[i]]);
    }

    return c;
}

/* A heuristic for the initial selections: prefer a period series
   named like "year"/"time"/"period"/"date" if there is one. */

static int preferred_period (const int *tlist)
{
    const char *pref[] = {"year", "time", "period", "date", "yr", NULL};
    int i, j;

    for (j=0; pref[j] != NULL; j++) {
	for (i=1; i<=tlist[0]; i++) {
	    if (!g_ascii_strcasecmp(dataset->varname[tlist[i]], pref[j])) {
		return i - 1;
	    }
	}
    }

    return 0;
}

void panel_markers_dialog (void)
{
    pm_info *pmi;
    GtkWidget *vbox, *hbox, *table, *frame, *w;

    if (!dataset_is_panel(dataset)) {
	errbox(_("This is not a panel dataset"));
	return;
    }

    pmi = calloc(1, sizeof *pmi);
    if (pmi == NULL) {
	return;
    }

    if (!get_candidates(pmi)) {
	infobox(_("Cannot find suitable series: we need one identifying the "
		  "units (ideally a string-valued series) and one "
		  "identifying the period (e.g. year)"));
	pm_info_free(NULL, pmi);
	return;
    }

    pmi->dlg = gretl_dialog_new(_("Observation markers"), mdata->main,
				GRETL_DLG_BLOCK);
    vbox = gtk_dialog_get_content_area(GTK_DIALOG(pmi->dlg));
    g_signal_connect(G_OBJECT(pmi->dlg), "destroy",
		     G_CALLBACK(pm_info_free), pmi);

    w = gtk_label_new(_("Construct observation markers of the form "
			"unit:period"));
    gtk_box_pack_start(GTK_BOX(vbox), w, FALSE, FALSE, 5);

    table = gtk_table_new(2, 2, FALSE);
    gtk_table_set_col_spacings(GTK_TABLE(table), 5);
    gtk_table_set_row_spacings(GTK_TABLE(table), 5);
    pmi->ucombo = pm_combo(table, 0, _("Unit identifier"), pmi->ulist, pmi);
    pmi->tcombo = pm_combo(table, 1, _("Period identifier"), pmi->tlist, pmi);
    gtk_box_pack_start(GTK_BOX(vbox), table, FALSE, FALSE, 5);

    frame = gtk_frame_new(_("Example"));
    pmi->preview = gtk_label_new(NULL);
    gtk_misc_set_alignment(GTK_MISC(pmi->preview), 0.0, 0.0);
    gtk_label_set_selectable(GTK_LABEL(pmi->preview), TRUE);
    gtk_container_add(GTK_CONTAINER(frame), pmi->preview);
    gtk_box_pack_start(GTK_BOX(vbox), frame, TRUE, TRUE, 5);

    /* buttons */
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(pmi->dlg));
    w = cancel_delete_button(hbox, pmi->dlg);
    pmi->ok = ok_button(hbox);
    g_signal_connect(G_OBJECT(pmi->ok), "clicked",
		     G_CALLBACK(pm_apply), pmi);
    gtk_widget_grab_default(pmi->ok);

    /* initial selection, then hook up the live preview */
    gtk_combo_box_set_active(GTK_COMBO_BOX(pmi->ucombo), 0);
    gtk_combo_box_set_active(GTK_COMBO_BOX(pmi->tcombo),
			     preferred_period(pmi->tlist));
    pm_update_preview(NULL, pmi);
    g_signal_connect(G_OBJECT(pmi->ucombo), "changed",
		     G_CALLBACK(pm_update_preview), pmi);
    g_signal_connect(G_OBJECT(pmi->tcombo), "changed",
		     G_CALLBACK(pm_update_preview), pmi);

    gtk_widget_show_all(pmi->dlg);
}
