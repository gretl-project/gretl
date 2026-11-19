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
    const int *ulist;    /* candidate unit series IDs */
    const int *tlist;    /* candidate period series IDs */
    int ulim;            /* max feasible unit-ID length */
} pm_info;

#define PM_NPREVIEW 4   /* number of leading rows to show */
#define PM_SEQ (-1)     /* "period" is just 1..T within each unit */

/* number of digits needed for the largest sequence number, T */

static int seq_width (void)
{
    return snprintf(NULL, 0, "%d", dataset->pd);
}

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

/* Returns 0 if @v is not a plausible unit identifier, otherwise 1
   (usable, but long strings) or 2 (usable, comfortable fit).
*/

static int is_unit_candidate (int v)
{
    const double *x = dataset->Z[v];
    int T = dataset->pd;
    int N = dataset->n / T;
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
	if (ret) {
	    const char *str = series_get_string_for_obs(dataset, v, i*T);
	    int len = str != NULL ? strlen(str) : 0;

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

    if (ret && maxlen > UNIT_STR_OK) {
	/* downgrade the return value somewhat */
	ret = 1;
    }

    free(first);

    return ret;
}

/* Returns 0 if @v is not a plausible period identifier, 1 if it is
   a numeric one, 2 if it is a string-valued one (such as "2008Q4"
   or "Jun1990"). For a string-valued series we cannot check ordering
   but require distinct values within each unit's block, with strings
   short enough to fit in a marker.
*/

static int is_period_candidate (int v)
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
	    int idx = i*T + s;
	    double xt = x[idx];

	    if (na(xt)) {
		ret = 0;
	    } else if (sv) {
		const char *str = series_get_string_for_obs(dataset, v, idx);

		if (str == NULL || strlen(str) > UNIT_STR_MAX) {
		    ret = 0;
		}
		blk[s] = xt;
	    } else if (xt < 0 || xt != floor(xt)) {
		/* we insist on non-negative integers */
		ret = 0;
	    } else if (s > 0 && xt <= x[idx - 1]) {
		/* and on monotonic increase */
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
	} else if (ret && x[0] == 1.0) {
	    /* a simple sequence 1, 2, 3,... adds no value */
	    int simple = 1;

	    for (s=1; s<T; s++) {
		if (x[s] != x[s-1] + 1) {
		    simple = 0;
		    break;
		}
	    }
	    if (simple) {
		ret = 0;
	    }
	}
    }

    free(blk);

    return ret;
}

/* Build candidate lists of series that identify panel units (@ul) and
   time (@tl), if possible, with the aim of offering apparatus to create
   series-based panel observation markers. This function is called from
   markers_callback() in library.c, where the condition for proceeding
   is that at least @ul is non-NULL.

   Unit candidates are ordered with any comfortably fitting string
   series (relative to OBSLEN) first, then string series that will only
   just fit. Period series may be numeric or string valued. If no
   period-identifying series are found we can substitute a simple
   numeric sequence.
*/

int get_candidate_lists (int **ul, int **tl)
{
    int *rank;
    int v, r;

    rank = calloc(dataset->v, sizeof *rank);
    if (rank == NULL) {
	return E_ALLOC;
    }

    for (v=1; v<dataset->v; v++) {
	int pc = is_period_candidate(v);

	if (pc == 1) {
	    gretl_list_append_term(tl, v);
	} else if (pc == 0 && is_string_valued(dataset, v)) {
	    rank[v] = is_unit_candidate(v);
	}
    }

    /* string-valued periods go after numeric ones */
    for (v=1; v<dataset->v; v++) {
	if (is_string_valued(dataset, v) && is_period_candidate(v) == 2) {
	    gretl_list_append_term(tl, v);
	}
    }

    /* units: string-valued only (best fit first) */
    for (r=2; r>0; r--) {
	for (v=1; v<dataset->v; v++) {
	    if (rank[v] == r && is_string_valued(dataset, v)) {
		gretl_list_append_term(ul, v);
	    }
	}
    }

    free(rank);

    return 0;
}

/* Return the length to which panel-unit strings must be truncated
   so that all obs markers fit into a buffer of length OBSLEN.
   In principle the result could be non-positive, which would
   indicate that the obs-markers recipe is infeasible.
*/

static int markers_dry_run (int u, int p)
{
    const char *s;
    double xt;
    int T = dataset->pd;
    int umax = 0;
    int pmax = 0;
    int i, t, len;

    /* examine the unit strings */
    for (i=0; i<dataset->n; i+=T) {
	s = series_get_string_for_obs(dataset, u, i);
	len = strlen(s);
	if (len > umax) {
	    umax = len;
	}
    }

    /* examine the period strings or numbers */
    for (t=0; t<T; t++) {
	if (p > 0) {
	    if (is_string_valued(dataset, p)) {
		s = series_get_string_for_obs(dataset, p, t);
		len = strlen(s);
	    } else {
		xt = floor(dataset->Z[p][t]);
		len = 1 + (int) floor(log10(xt));
	    }
	} else {
	    len = 1 + (int) floor(log10((double) T+1));
	}
	if (len > pmax) {
	    pmax = len;
	}
    }

    return OBSLEN - pmax - 2;
}

/* Write the marker for row @t into @targ (size OBSLEN) */

static int make_marker (int u, int p, int t,
			int ulim, char *targ)
{
    char *q = targ;
    char *tmp;
    int k;

    for (k=0; k<2; k++) {
	int v = k == 0 ? u : p;

	if (k == 1) {
	    *q++ = ':';
	}
	if (v == PM_SEQ) {
	    /* automatic numeric period ID */
	    q += sprintf(q, "%0*d", seq_width(), t % dataset->pd + 1);
	} else if (is_string_valued(dataset, v)) {
	    const char *s = series_get_string_for_obs(dataset, v, t);

	    if (s == NULL) {
		; /* "can't happen" */
	    } else if (v == u) {
		/* unit string */
		if (strlen(s) > ulim) {
		    tmp = g_strdup(s);
		    gretl_utf8_truncate_b(tmp, ulim);
		    q += sprintf(q, "%s", tmp);
		    g_free(tmp);
		} else {
		    q += snprintf(q, ulim+1, "%s", s);
		}
	    } else {
		/* period string */
		q += sprintf(q, "%s", s);
	    }
	} else {
	    /* numeric series */
	    q += sprintf(q, "%d", (int) dataset->Z[v][t]);
	}
    }

    return strlen(targ);
}

static int combo_series_id (GtkWidget *combo, const int *list)
{
    int k = gtk_combo_box_get_active(GTK_COMBO_BOX(combo));

    if (k < 0) {
	return 0;
    }

    /* the period combo has one extra entry at the end */
    return k >= list[0] ? PM_SEQ : list[k+1];
}

/* Refresh the preview and the OK button: warn about markers that
   would be truncated (OBSLEN) or non-unique.
*/

static void pm_update_preview (GtkWidget *w, pm_info *pmi)
{
    int u = combo_series_id(pmi->ucombo, pmi->ulist);
    int p = combo_series_id(pmi->tcombo, pmi->tlist);
    GString *gs = g_string_new(NULL);
    int ok = (u > 0 && p != 0 && u != p);

    if (ok) {
	GHashTable *seen = NULL;
	char mk[OBSLEN];
	int ndup = 0;
	int t;

	pmi->ulim = markers_dry_run(u, p);
	seen = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, NULL);
	for (t=0; t<dataset->n; t++) {
	    make_marker(u, p, t, pmi->ulim, mk);
	    if (!g_hash_table_add(seen, g_strdup(mk))) {
		ndup++;
	    }
	    if (t < PM_NPREVIEW) {
		g_string_append_printf(gs, "%s\n", mk);
	    } else if (t == PM_NPREVIEW) {
		g_string_append(gs, "...\n");
	    }
	}
	if (ndup > 0) {
	    g_string_append_printf(gs, _("\nWarning: %d markers are not unique"),
				   ndup);
	}
	g_hash_table_destroy(seen);
    } else {
	g_string_append(gs, _("Please select two different series"));
    }

    gtk_label_set_text(GTK_LABEL(pmi->preview), gs->str);
    gtk_widget_set_sensitive(pmi->ok, ok);
    g_string_free(gs, TRUE);
}

static void pm_record (int u, int p, int ulim)
{
    if (p == PM_SEQ) {
	/* name for temporary sequence series */
	char tmp[VNAMELEN] = "pmseq";
	int i;

	/* find a name that doesn't clash with an existing series */
	for (i=2; current_series_index(dataset, tmp) >= 0; i++) {
	    sprintf(tmp, "pmseq%d", i);
	}

	lib_command_sprintf("series %s = $obsminor", tmp);
	record_command_verbatim();
	lib_command_sprintf("genr markers = \"%%.%d%c:%%0%dd\", %s, %s",
			    ulim, is_string_valued(dataset, u) ? 's' : 'd',
			    seq_width(), dataset->varname[u], tmp);
	record_command_verbatim();
	lib_command_sprintf("delete %s", tmp);
	record_command_verbatim();
    } else {
	lib_command_sprintf("genr markers = \"%%.%d%c:%%%c\", %s, %s", ulim,
			    is_string_valued(dataset, u) ? 's' : 'd',
			    is_string_valued(dataset, p) ? 's' : 'd',
			    dataset->varname[u], dataset->varname[p]);
	record_command_verbatim();
    }
}

static void pm_apply (GtkWidget *w, pm_info *pmi)
{
    int u = combo_series_id(pmi->ucombo, pmi->ulist);
    int p = combo_series_id(pmi->tcombo, pmi->tlist);
    int i, err;

    err = dataset_allocate_obs_markers(dataset);

    if (err) {
	gui_errmsg(err);
    } else {
	for (i=0; i<dataset->n; i++) {
	    make_marker(u, p, i, pmi->ulim, dataset->S[i]);
	}
	pm_record(u, p, pmi->ulim);
	mark_dataset_as_modified();
	gtk_widget_destroy(pmi->dlg);
    }
}

static void pm_info_free (GtkWidget *w, pm_info *pmi)
{
    free(pmi);
}

static GtkWidget *pm_combo (GtkWidget *table, int row, const char *label,
			    const int *list, const char *extra,
			    pm_info *pmi)
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
    if (extra != NULL) {
	combo_box_append_text(c, extra);
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

void panel_markers_dialog (const int *ulist, const int *tlist)
{
    pm_info *pmi;
    GtkWidget *vbox, *hbox, *table, *frame, *w;

    pmi = calloc(1, sizeof *pmi);
    if (pmi == NULL) {
	return;
    }

    pmi->ulist = ulist;
    pmi->tlist = tlist;

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
    pmi->ucombo = pm_combo(table, 0, _("Unit identifier"), pmi->ulist, NULL, pmi);
    pmi->tcombo = pm_combo(table, 1, _("Period identifier"), pmi->tlist,
			  _("(sequence number 1, 2, ...)"), pmi);
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

/* Preliminary dialog to offer a choice: applicable only if we found at
   least one series suitable for identifying the panel unit.
*/

int panel_markers_choice (void)
{
    const char *opts[] = {
	N_("Construct markers from series"),
	N_("Add markers from text file")
    };

    return radio_dialog(_("Observation markers"),
			_("The dataset has no observation markers"),
			opts, 2, 0, 0 /* help code? */, NULL);
}
