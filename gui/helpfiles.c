/*
 *  Copyright (c) by Allin Cottrell
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
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

/* helpfiles.c for gretl */

#include "gretl.h"

/* helpfile stuff */
struct help_head_t {
    char *name;
    int *topics;
    int *pos;
    int ntopics;
};

static int gui_help_length, script_help_length;
static struct help_head_t **cli_heads, **gui_heads;

/* searching stuff */
static int look_for_string (const char *haystack, const char *needle, 
			    int start);
static void close_find_dialog (GtkWidget *widget, gpointer data);
static void find_in_help (GtkWidget *widget, gpointer data);
static void find_in_text (GtkWidget *widget, gpointer data);
static void find_in_clist (GtkWidget *widget, gpointer data);
static void find_string_dialog (void (*findfunc)(), gpointer data);

static GtkWidget *find_window = NULL;
static GtkWidget *find_entry;
static char *needle;

GtkItemFactoryEntry help_items[] = {
    { N_("/_Topics"), NULL, NULL, 0, "<Branch>" },    
    { N_("/_Find"), NULL, menu_find, 0, NULL },
    { NULL, NULL, NULL, 0, NULL}
};

/* ......................................................... */

struct gui_help_item {
    int code;
    char *string;
};

static struct gui_help_item gui_help_items[] = {
    { GR_PLOT,    "graphing" },
    { GR_XY,      "graphing" },
    { GR_DUMMY,   "factorized" },
    { GR_BOX,     "boxplots" },
    { GR_NBOX,    "boxplots" },
    { RANGE_MEAN, "range-mean" },
    { ONLINE,     "online" },
    { MARKERS,    "markers" },
    { EXPORT,     "export" },
    { SMPLBOOL,   "sampling" },
    { SMPLDUM,    "sampling" },
    { COMPACT,    "compact" },
    { VSETMISS,   "missing" },
    { GSETMISS,   "missing" },
    { GUI_HELP,   "dialog" },
    { 0,          NULL },
};

/* state the topic headings from the help files so they 
   can be translated */
const char *intl_topics[] = {
    N_("Dataset"),
    N_("Estimation"),
    N_("Graphs"),
    N_("Prediction"),
    N_("Printing"),
    N_("Programming"),
    N_("Statistics"),
    N_("Tests"),
    N_("Transformations"),
    N_("Utilities")
};

/* ......................................................... */

static int extra_command_number (const char *s)
{
    int i;

    for (i=0; gui_help_items[i].code; i++)
	if (!strcmp(s, gui_help_items[i].string))
	    return gui_help_items[i].code;
    return 0;
}

/* ......................................................... */

static char *help_string_from_cmd (int cmd)
{
    int i;

    for (i=0; gui_help_items[i].code; i++)
	if (cmd == gui_help_items[i].code)
	    return gui_help_items[i].string;
    return NULL;    
}

/* ......................................................... */

static int real_helpfile_init (int cli)
{
    FILE *fp;
    char *helpfile, *headstr;
    struct help_head_t **heads = NULL;
    char testline[MAXLEN], topicword[32];
    int i, g, pos, match, nheads = 0, topic = 0;
    int length = 0, memfail = 0;

    helpfile = (cli)? paths.cmd_helpfile : paths.helpfile;

    /* first pass: find length and number of topics */
    fp = fopen(helpfile, "r");
    if (fp == NULL) {
	fprintf(stderr, _("help file %s is not accessible\n"), helpfile);
	return -1;
    }

    while (!memfail && fgets(testline, MAXLEN-1, fp)) {
	if (*testline == '@') {
	    chopstr(testline);
	    match = 0;
	    for (i=0; i<nheads; i++) {
		if (!strcmp(testline + 1, (heads[i])->name)) {
		    match = 1;
		    (heads[i])->ntopics += 1;
		    break;
		}
	    }
	    if (!match) {
		heads = realloc(heads, (nheads + 2) * sizeof *heads);
		if (heads != NULL) { 
		    heads[nheads] = malloc(sizeof **heads);
		    if (heads[nheads] != NULL) {
			headstr = testline + 1;
			(heads[nheads])->name = malloc(strlen(headstr) + 1);
			if ((heads[nheads])->name != NULL) {
			    strcpy((heads[nheads])->name, headstr);
			    (heads[nheads])->ntopics = 1;
			    nheads++;
			} else memfail = 1;
		    } else memfail = 1;
		} else memfail = 1;
	    }
	} else length++;
    }
    fclose(fp);

    if (memfail) return -1;

    for (i=0; i<nheads; i++) {
	(heads[i])->topics = malloc((heads[i])->ntopics * sizeof(int));
	if ((heads[i])->topics == NULL) memfail = 1;
	(heads[i])->pos = malloc((heads[i])->ntopics * sizeof(int));
	if ((heads[i])->pos == NULL) memfail = 1; 
	(heads[i])->ntopics = 0;
    }
    heads[i] = NULL;

    if (memfail) return -1;

    /* second pass, assemble the topic list */
    fp = fopen(helpfile, "r");
    i = 0;
    pos = 0;
    g = 0;
    while (!memfail && fgets(testline, MAXLEN-1, fp)) {
	if (topic == 1) 
	    sscanf(testline, "%31s", topicword);
	if (*testline == '@') {
	    chopstr(testline);
	    match = -1;
	    for (i=0; i<nheads; i++) {
		if (!strcmp(testline + 1, (heads[i])->name)) {
		    match = i;
		    break;
		}
	    }
	    if (match >= 0) {
		int t, m = (heads[match])->ntopics;

		t = command_number(topicword);
		if (t) (heads[match])->topics[m] = t;
		else (heads[match])->topics[m] = 
			 extra_command_number(topicword);
		(heads[match])->pos[m] = pos - 1;
		(heads[match])->ntopics += 1;
	    }		
	} else pos++;
	if (*testline == '#') topic = 1;
	else topic = 0;
    }
    fclose(fp);

    if (cli) cli_heads = heads;
    else gui_heads = heads;

    return length;
}

/* ......................................................... */

void helpfile_init (void)
{
    gui_help_length = real_helpfile_init(0);
    script_help_length = real_helpfile_init(1);
}

/* ......................................................... */

static char *get_gui_help_string (int pos)
{
    int i, j;

    for (i=0; gui_heads[i] != NULL; i++) 
	for (j=0; j<(gui_heads[i])->ntopics; j++)
	    if (pos == (gui_heads[i])->pos[j])
		return help_string_from_cmd((gui_heads[i])->topics[j]);
    return NULL;
}

/* ........................................................... */

static void add_help_topics (windata_t *hwin, int script)
{
    int i, j;
    GtkItemFactoryEntry helpitem;
    const gchar *mpath = N_("/_Topics");
    struct help_head_t **heads = (script)? cli_heads : gui_heads;

    helpitem.path = NULL;

    /* See if there are any topics to add */
    if (heads == NULL) return;

    /* put the topics under the menu heading */
    for (i=0; heads[i] != NULL; i++) {
	if (helpitem.path == NULL)
	    helpitem.path = mymalloc(80);
	helpitem.accelerator = NULL;
	helpitem.callback_action = 0; 
	helpitem.item_type = "<Branch>";
	sprintf(helpitem.path, "%s/%s", mpath, _((heads[i])->name));
	helpitem.callback = NULL; 
	gtk_item_factory_create_item(hwin->ifac, &helpitem, NULL, 1);
	for (j=0; j<(heads[i])->ntopics; j++) {
	    helpitem.accelerator = NULL;
	    helpitem.callback_action = (heads[i])->pos[j]; 
	    helpitem.item_type = NULL;
	    if ((heads[i])->topics[j] < NC) {
		sprintf(helpitem.path, "%s/%s/%s", 
			mpath, _((heads[i])->name), 
			commands[(heads[i])->topics[j]]);
	    } else {
		sprintf(helpitem.path, "%s/%s/%s", 
			mpath, _((heads[i])->name), 
			get_gui_help_string((heads[i])->pos[j]));
	    }
	    helpitem.callback = (script)? do_script_help : do_gui_help; 
	    gtk_item_factory_create_item(hwin->ifac, &helpitem, NULL, 1);
	}
    }
    free(helpitem.path);
}

/* ........................................................... */

static windata_t *helpwin (int script) 
{
    windata_t *vwin = NULL;

    if (script) {
	vwin = view_file(paths.cmd_helpfile, 0, 0, 77, 400, 
			 CLI_HELP, help_items);
	add_help_topics(vwin, 1);
    } else {
	vwin = view_file(paths.helpfile, 0, 0, 77, 400, 
			 HELP, help_items);
	add_help_topics(vwin, 0);
    }
    return vwin;
}

/* ........................................................... */

void context_help (GtkWidget *widget, gpointer data)
{
    int i, j, help_code = GPOINTER_TO_INT(data);
    int pos = 0;

    for (i=0; gui_heads[i] != NULL; i++) {
	for (j=0; j<(gui_heads[i])->ntopics; j++)
	    if (help_code == (gui_heads[i])->topics[j])
		pos = (gui_heads[i])->pos[j];
    }
    /* fallback */
    if (!pos) {
	char *helpstr = help_string_from_cmd(help_code);
	int altcode;

	if (helpstr != NULL) {
	    altcode = extra_command_number(helpstr);
	    for (i=0; gui_heads[i] != NULL; i++)
		for (j=0; j<(gui_heads[i])->ntopics; j++)
		    if (altcode == (gui_heads[i])->topics[j])
			pos = (gui_heads[i])->pos[j];
	}
    }
    do_gui_help(NULL, pos, NULL);
}

/* ........................................................... */

static void real_do_help (guint pos, int cli)
{
    double frac;
    gfloat adj;
    static GtkWidget *gui_help_view;
    static GtkWidget *script_help_view;
    GtkWidget *w = (cli)? script_help_view : gui_help_view;

    if (w == NULL) {
	windata_t *hwin = helpwin(cli);

	if (hwin != NULL) {
	    if (cli) w = script_help_view = hwin->w;
	    else w = gui_help_view = hwin->w;
	}
	gtk_signal_connect(GTK_OBJECT(w), "destroy",
			   GTK_SIGNAL_FUNC(gtk_widget_destroyed),
			   (cli)? &script_help_view : &gui_help_view);	
    } else {
	gdk_window_show(w->parent->window);
	gdk_window_raise(w->parent->window);
    }
    
    frac = (double) pos * (double) GTK_TEXT(w)->vadj->upper;
    frac /= (double) (cli)? script_help_length : gui_help_length;
    adj = 0.999 * frac;

    gtk_adjustment_set_value(GTK_TEXT(w)->vadj, adj);
}

/* ........................................................... */

void do_gui_help (gpointer data, guint pos, GtkWidget *widget) 
{
    real_do_help(pos, 0);
}

/* ........................................................... */

void do_script_help (gpointer data, guint pos, GtkWidget *widget) 
{
    real_do_help(pos, 1);
}

/* ........................................................... */

static int pos_from_cmd (int cmd)
{
    int i, j;

    for (i=0; cli_heads[i] != NULL; i++)
	for (j=0; j<(cli_heads[i])->ntopics; j++)
	    if (cmd == (cli_heads[i])->topics[j])
		return (cli_heads[i])->pos[j];
    return 0;
}

/* ........................................................... */

void edit_script_help (GtkWidget *widget, GdkEventButton *b,
		       windata_t *vwin)
{
    if (!vwin->help_active) { /* command help not activated */
	return;
    } else {
	gchar *text;
	guint pt = GTK_EDITABLE(vwin->w)->current_pos;
	int len = gtk_text_get_length(GTK_TEXT(vwin->w));
	int pos = 0;

	text = gtk_editable_get_chars(GTK_EDITABLE(vwin->w), 
				      0, (pt + 9 > len)? -1 : pt + 8);

	if (text != NULL && strlen(text) > 0) {
	    char *p, *q;
	    char word[9];

	    p = q = text + pt;
	    if (pt > 0)
		while (p - text && !isspace(*(p-1))) p--;
	    if (pt < strlen(text))
		while (*q && !isspace(*q)) q++;
	    *word = '\0';
	    strncat(word, p, (q - p > 8)? 8 : q - p);
	    pos = pos_from_cmd(command_number(word));
	} 
	
	real_do_help(pos, 1);
	g_free(text);
	gdk_window_set_cursor(GTK_TEXT(vwin->w)->text_area, NULL);
	vwin->help_active = 0;
    }
}

/* ........................................................... */

void menu_find (gpointer data, guint db, GtkWidget *widget)
{
    if (db) 
	find_string_dialog(find_in_clist, data);
    else 
	find_string_dialog(find_in_help, data);
}

/* ........................................................... */

void datafile_find (GtkWidget *widget, gpointer data)
{
    find_string_dialog(find_in_clist, data);
}

/* ........................................................... */

void find_var (gpointer p, guint u, GtkWidget *w)
{
    find_string_dialog(find_in_clist, mdata);
}

/* .................................................................. */

static void find_in_help (GtkWidget *widget, gpointer data)
{
    int found = 0, i, linecount = 0;
    int help_length;
    char *haystack;
    windata_t *vwin = 
	(windata_t *) gtk_object_get_data(GTK_OBJECT(data), "windat");

    haystack = gtk_editable_get_chars(GTK_EDITABLE(vwin->w), 0,
	gtk_text_get_length(GTK_TEXT(vwin->w)));

    if (vwin->role == CLI_HELP) help_length = script_help_length;
    else help_length = gui_help_length;

    if (needle) g_free(needle);

    needle = gtk_editable_get_chars(GTK_EDITABLE (find_entry), 0, -1);
    found = GTK_EDITABLE(vwin->w)->selection_end_pos;

    found = look_for_string(haystack, needle, found);

    if (found >= 0) {
	gtk_text_freeze(GTK_TEXT(vwin->w));
        gtk_text_set_point (GTK_TEXT(vwin->w), found);
        gtk_text_insert (GTK_TEXT(vwin->w), NULL, NULL, NULL, " ", 1);
        gtk_text_backward_delete (GTK_TEXT(vwin->w), 1);
	gtk_text_thaw(GTK_TEXT(vwin->w));
        gtk_editable_select_region (GTK_EDITABLE(vwin->w), 
				    found, found + strlen(needle));
	for (i=0; i<found; i++) 
	    if (haystack[i] == '\n') linecount++;
	gtk_adjustment_set_value(GTK_TEXT(vwin->w)->vadj, 
				 (gfloat) (linecount - 2) *
				 GTK_TEXT(vwin->w)->vadj->upper / help_length);
    } else infobox(_("String was not found."));

    g_free(haystack);
}

/* .................................................................. */

static void close_find_dialog (GtkWidget *widget, gpointer data)
{
    find_window = NULL;
}

/* .................................................................. */

static void find_in_text (GtkWidget *widget, gpointer data)
{
    int found = 0;
    char *haystack;
    windata_t *vwin = 
	(windata_t *) gtk_object_get_data(GTK_OBJECT(data), "windat");

    haystack = gtk_editable_get_chars(GTK_EDITABLE(vwin->w), 0,
	gtk_text_get_length(GTK_TEXT(vwin->w)));

    if (needle) g_free(needle);

    needle = gtk_editable_get_chars(GTK_EDITABLE(find_entry), 0, -1);
    found = GTK_EDITABLE(vwin->w)->selection_end_pos;

    found = look_for_string(haystack, needle, found);

    if (found >= 0) {
	gtk_text_set_point(GTK_TEXT(vwin->w), found);
	gtk_editable_set_position(GTK_EDITABLE(vwin->w), found);
        gtk_editable_select_region(GTK_EDITABLE(vwin->w), 
				   found, found + strlen(needle));
    } else infobox(_("String was not found."));

    g_free(haystack);
}

/* .................................................................. */

static void find_in_clist (GtkWidget *w, gpointer data)
{
    int start, found = 0, n, i;
    gchar *tmp; 
    char haystack[MAXLEN];
    windata_t *dbdat;

    dbdat = (windata_t *) gtk_object_get_data(GTK_OBJECT(data), "windat");

    if (needle) g_free(needle);
    needle = gtk_editable_get_chars(GTK_EDITABLE (find_entry), 0, -1);
    lower(needle);

    start = dbdat->active_var + 1;
    n = GTK_CLIST(dbdat->listbox)->rows;

    for (i=start; i<n; i++) {  
	/* try looking in column 1 first */
	gtk_clist_get_text(GTK_CLIST(dbdat->listbox), i, 1, &tmp);
	strcpy(haystack, tmp);
	lower(haystack);
	found = look_for_string(haystack, needle, 0);
	if (found >= 0) break;
	else { /* try column 0? */
	    gtk_clist_get_text(GTK_CLIST(dbdat->listbox), i, 0, &tmp);
	    strcpy(haystack, tmp);
	    lower(haystack);
	    found = look_for_string(haystack, needle, 0);
	    if (found >= 0) break;
	}
    }

    if (found >= 0) {
	gtk_clist_moveto(GTK_CLIST(dbdat->listbox), i, 0, 0, .1);
	gtk_clist_select_row(GTK_CLIST(dbdat->listbox), i, 0);
	dbdat->active_var = i;
    } else {
	gtk_clist_select_row(GTK_CLIST(dbdat->listbox), 0, 0);
	dbdat->active_var = 0;
	infobox(_("String was not found."));
    }
}

/* .................................................................. */

static int look_for_string (const char *haystack, const char *needle, 
			    int start)
{
    int pos;
    int hlen = strlen(haystack);
    int nlen = strlen(needle);

    for (pos = start; pos < hlen; pos++) {
        if (strncmp(&haystack[pos], needle, nlen) == 0) 
             return pos;
    }
    return -1;
}

/* .................................................................. */
 
static void cancel_find (GtkWidget *widget, gpointer data)
{
    if (find_window != NULL) {
	gtk_widget_destroy(GTK_WIDGET(data));
	find_window = NULL;
    }
}

/* .................................................................. */

static void parent_find (GtkWidget *finder, windata_t *caller)
{
    if (caller != mdata) {
	GtkWidget *w = NULL;

	if (caller->dialog != NULL) {
	    w = caller->dialog;
	} else if (caller->w != NULL) {
	    w = caller->w;
	}

	if (w == NULL) return;

	gtk_window_set_transient_for(GTK_WINDOW(finder),
				     GTK_WINDOW(w));
	gtk_signal_connect(GTK_OBJECT(w), "destroy",
			   GTK_SIGNAL_FUNC(cancel_find),
			   finder);
    }
}

/* .................................................................. */

static void find_string_dialog (void (*findfunc)(), gpointer data)
{
    GtkWidget *label;
    GtkWidget *button;
    GtkWidget *hbox;
    windata_t *mydat = (windata_t *) data;

    if (find_window) {
	gtk_object_set_data(GTK_OBJECT(find_window), "windat", mydat); 
	parent_find(find_window, mydat);
	return;
    }

    find_window = gtk_dialog_new();
    gtk_object_set_data(GTK_OBJECT(find_window), "windat", mydat);
    parent_find(find_window, mydat);

    gtk_signal_connect (GTK_OBJECT (find_window), "destroy",
	                GTK_SIGNAL_FUNC (close_find_dialog),
	                find_window);

    gtk_window_set_title (GTK_WINDOW (find_window), _("gretl: find"));
    gtk_container_border_width (GTK_CONTAINER (find_window), 5);

    hbox = gtk_hbox_new(TRUE, TRUE);
    label = gtk_label_new(_(" Find what:"));
    gtk_widget_show (label);
    find_entry = gtk_entry_new();

    if (needle) {
	gtk_entry_set_text(GTK_ENTRY (find_entry), needle);
	gtk_entry_select_region (GTK_ENTRY (find_entry), 0, 
				 strlen (needle));
    }
    gtk_signal_connect(GTK_OBJECT (find_entry), 
			"activate", 
			GTK_SIGNAL_FUNC(findfunc),
	                find_window);
    gtk_widget_show (find_entry);

    gtk_box_pack_start (GTK_BOX(hbox), label, TRUE, TRUE, 0);
    gtk_box_pack_start (GTK_BOX(hbox), find_entry, TRUE, TRUE, 0);
    gtk_widget_show (hbox);

    gtk_box_pack_start(GTK_BOX (GTK_DIALOG (find_window)->vbox), 
                        hbox, TRUE, TRUE, 0);

    gtk_box_set_spacing(GTK_BOX (GTK_DIALOG (find_window)->action_area), 15);
    gtk_box_set_homogeneous(GTK_BOX 
			     (GTK_DIALOG (find_window)->action_area), TRUE);
    gtk_window_set_position(GTK_WINDOW (find_window), GTK_WIN_POS_MOUSE);

    /* find button -- make this the default */
    button = gtk_button_new_with_label (_("Find next"));
    GTK_WIDGET_SET_FLAGS(button, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX (GTK_DIALOG (find_window)->action_area), 
		       button, TRUE, TRUE, FALSE);
    gtk_signal_connect(GTK_OBJECT (button), "clicked",
		       GTK_SIGNAL_FUNC (findfunc), find_window);
    gtk_widget_grab_default(button);
    gtk_widget_show(button);

    /* cancel button */
    button = gtk_button_new_with_label (_("Cancel"));
    GTK_WIDGET_SET_FLAGS(button, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX (GTK_DIALOG (find_window)->action_area), 
		       button, TRUE, TRUE, FALSE);
    gtk_signal_connect(GTK_OBJECT (button), "clicked",
		       GTK_SIGNAL_FUNC (cancel_find), find_window);
    gtk_widget_show(button);

    gtk_widget_grab_focus(find_entry);
    gtk_widget_show (find_window);
}

/* ........................................................... */

void text_find_callback (GtkWidget *w, gpointer data)
{
    find_string_dialog(find_in_text, data);
}

/* ........................................................... */

void colorize_tooltips (GtkTooltips *tip)
{
    GdkColor t_back;
    GtkStyle *style;

    if (gdk_color_parse("light yellow", &t_back)) {
	gtk_tooltips_force_window(tip);
	if (gdk_color_alloc(gtk_widget_get_colormap(tip->tip_window), &t_back)) {
	    style = gtk_style_copy(gtk_widget_get_style(tip->tip_window));
	    style->bg[GTK_STATE_NORMAL] = t_back;
	    gtk_widget_set_style(tip->tip_window, style);
	} 
    } 
}

static GtkTooltips *gretl_tips;

void gretl_tooltips_init (void)
{
    gretl_tips = gtk_tooltips_new();
    gtk_tooltips_enable(gretl_tips); /* redundant? */
    /* colorize_tooltips(gretl_tips); */
}

void gretl_tooltips_add (GtkWidget *w, const gchar *str)
{
    gtk_tooltips_set_tip(gretl_tips, w, str, NULL);
}
