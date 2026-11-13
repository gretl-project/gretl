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
   and a period-identifier series (e.g. year, or a string date), show
   a live preview of the resulting markers, and then execute

     genr markers = "%s:%d", <unit>, <period>

   with %s or %d for each role according to whether the series is
   string-valued or numeric.
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

/* The dataset is already a panel, with N = n/T units of T rows each.
   Use that structure to vet candidate series. Always applied over the
   whole dataset, since markers apply to all observations whatever the
   current sample.

   Unit candidate: constant within each block of T rows, with distinct
   values across blocks. If string-valued, the longest string must also
   leave room in the marker (OBSLEN - 1 chars) for ':' and at least one
   digit of the period.

   Period candidate: integer-valued, non-negative, and strictly
   increasing within each block.
*/

#define UNIT_STR_MAX (OBSLEN - 3)   /* longest usable unit string */
#define UNIT_STR_OK  (OBSLEN - 6)   /* leaves room for a 4-digit year */

static int cmp_doubles (const void *a, const void *b)
{
    double x = *(const double *) a;
    double y = *(const double *) b;

    return (x > y) - (x < y);
}

/* returns 0 if @v is not a plausible unit identifier, otherwise 1
   (usable, but long strings) or 2 (usable, comfortable fit) */

static int unit_candidate (int v)
{
    const double *x = dataset->Z[v];
    int T = dataset->pd;
    int N = dataset->n / T;
    int sv = is_string_valued(dataset, v);
    int maxlen = 0;
    double *first;
    int i, s, ret = 2;

    first = malloc(N * sizeof *first);
    if (first == NULL) {
	return 0;
    }

    for (i=0; i<N && ret; i++) {
	first[i] = x[i*T];
	if (na(first[i])) {
	    ret = 0;
	}
	for (s=1; s<T && ret; s++) {
	    if (x[i*T + s] != first[i]) {
		ret = 0; /* not constant within unit */
	    }
	}
	if (ret && sv) {
	    const char *str = series_get_string_for_obs(dataset, v, i*T);
	    int len = str != NULL ? g_utf8_strlen(str, -1) : 0;

	    if (len > maxlen) {
		maxlen = len;
	    }
	}
    }

    if (ret) {
	/* distinct across units? */
	qsort(first, N, sizeof *first, cmp_doubles);
	for (i=1; i<N; i++) {
	    if (first[i] == first[i-1]) {
		ret = 0;
		break;
	    }
	}
    }

    if (ret && sv) {
	if (maxlen > UNIT_STR_MAX) {
	    ret = 0;
	} else if (maxlen > UNIT_STR_OK) {
	    ret = 1;
	}
    }

    free(first);

    return ret;
}

/* returns 0 if @v is not a plausible period identifier, 1 if it is
   a numeric one, 2 if it is a string-valued one (such as "2008Q4"
   or "Jun1990"). For a string-valued series we cannot check ordering
   but require distinct values within each unit's block, with strings
   short enough to fit in a marker.
*/

static int period_candidate (int v)
{
    const double *x = dataset->Z[v];
    int T = dataset->pd;
    int N = dataset->n / T;
    int sv = is_string_valued(dataset, v);
    double *blk = NULL;
    int i, s, ret = sv ? 2 : 1;

    if (sv) {
	blk = malloc(T * sizeof *blk);
	if (blk == NULL) {
	    return 0;
	}
    }

    for (i=0; i<N && ret; i++) {
	for (s=0; s<T && ret; s++) {
	    double xt = x[i*T + s];

	    if (na(xt)) {
		ret = 0;
	    } else if (sv) {
		const char *str = series_get_string_for_obs(dataset, v, i*T + s);

		if (str == NULL || g_utf8_strlen(str, -1) > UNIT_STR_MAX) {
		    ret = 0;
		}
		blk[s] = xt;
	    } else if (xt < 0 || xt != floor(xt)) {
		ret = 0;
	    } else if (s > 0 && xt <= x[i*T + s - 1]) {
		ret = 0;
	    }
	}
	if (ret && sv) {
	    qsort(blk, T, sizeof *blk, cmp_doubles);
	    for (s=1; s<T; s++) {
		if (blk[s] == blk[s-1]) {
		    ret = 0;
		    break;
		}
	    }
	}
    }

    free(blk);

    return ret;
}

/* Build the candidate lists. Unit candidates are ordered with the
   comfortable string series first, then string series that will only
   just fit, then numeric ones. Returns 0 if we can't offer both roles.
*/

static int get_candidates (pm_info *pmi)
{
    int *ul = gretl_list_new(dataset->v);
    int *tl = gretl_list_new(dataset->v);
    int *rank = NULL;
    int i, r;

    if (ul == NULL || tl == NULL || dataset->pd < 2 ||
	dataset->n < 2 * dataset->pd) {
	goto bailout;
    }

    rank = calloc(dataset->v, sizeof *rank);
    if (rank == NULL) {
	goto bailout;
    }

    ul[0] = tl[0] = 0;

    for (i=1; i<dataset->v; i++) {
	int pc = period_candidate(i);

	if (pc == 1) {
	    tl[++tl[0]] = i;
	} else if (pc == 0) {
	    rank[i] = unit_candidate(i);
	}
    }
    /* string-valued periods go after numeric ones */
    for (i=1; i<dataset->v; i++) {
	if (is_string_valued(dataset, i) && period_candidate(i) == 2) {
	    tl[++tl[0]] = i;
	}
    }

    /* string-valued first (best fit first), then numeric */
    for (r=2; r>=1; r--) {
	for (i=1; i<dataset->v; i++) {
	    if (rank[i] == r && is_string_valued(dataset, i)) {
		ul[++ul[0]] = i;
	    }
	}
    }
    for (i=1; i<dataset->v; i++) {
	if (rank[i] > 0 && !is_string_valued(dataset, i)) {
	    ul[++ul[0]] = i;
	}
    }

    free(rank);
    pmi->ulist = ul;
    pmi->tlist = tl;

    return ul[0] > 0 && tl[0] > 0;

 bailout:
    free(ul);
    free(tl);
    free(rank);

    return 0;
}

/* write the marker for row @t into @targ (size OBSLEN); returns 1 if
   the result had to be truncated */

static int make_marker (int u, int p, int t, char *targ)
{
    char full[2 * (2 * VNAMELEN) + 8];
    char *q = full;
    int k;

    for (k=0; k<2; k++) {
	int v = k == 0 ? u : p;

	if (k == 1) {
	    *q++ = ':';
	}
	if (is_string_valued(dataset, v)) {
	    const char *s = series_get_string_for_obs(dataset, v, t);

	    q += sprintf(q, "%s", s != NULL ? s : "");
	} else {
	    q += sprintf(q, "%d", (int) dataset->Z[v][t]);
	}
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
    gchar *line;
    PRN *prn;
    int err;

    if (u <= 0 || p <= 0 || u == p || bufopen(&prn)) {
	return;
    }

    line = g_strdup_printf("markers = \"%%%c:%%%c\", %s, %s",
			   is_string_valued(dataset, u) ? 's' : 'd',
			   is_string_valued(dataset, p) ? 's' : 'd',
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
