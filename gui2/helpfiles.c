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
#include "treeutils.h"

#ifdef ENABLE_NLS
static int translated_helpfile = -1;
static char *english_gui_helpfile;
static char *english_script_helpfile;
static int english_gui_help_length;
static int english_script_help_length;
#endif

/* helpfile stuff */
struct help_head_t {
    char *name;
    int *topics;
    int *pos;
    int ntopics;
};

static int gui_help_length, script_help_length;
static struct help_head_t **cli_heads, **gui_heads;

static windata_t *helpwin (int script, int english);

/* searching stuff */
static int look_for_string (const char *haystack, const char *needle, 
			    int start);
static void find_in_text (GtkWidget *widget, gpointer data);
static void find_in_listbox (GtkWidget *widget, gpointer data);
static void find_string_dialog (void (*findfunc)(), gpointer data);

static GtkWidget *find_window = NULL;
static GtkWidget *find_entry;
static char *needle;

GtkItemFactoryEntry help_items[] = {
    { N_("/_Topics"), NULL, NULL, 0, "<Branch>", GNULL },    
    { N_("/_Find"), NULL, menu_find, 0, NULL, GNULL },
    { NULL, NULL, NULL, 0, NULL, GNULL }
};

#ifdef ENABLE_NLS
GtkItemFactoryEntry english_help_items[] = {
    { N_("/_Find"), NULL, menu_find, 0, NULL, GNULL },
    { NULL, NULL, NULL, 0, NULL, GNULL }
};
#endif

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
    { GR_3D,      "3-D" },
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
    { MODELTABLE, "model" },
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

#ifdef ENABLE_NLS
static void set_english_help_file (int script)
{
    char *helpfile, *tmp, *p;
    FILE *fp;

    if (script) helpfile = paths.cmd_helpfile;
    else helpfile = paths.helpfile;

    tmp = malloc(strlen(helpfile) + 1);

    if (tmp != NULL) {
	strcpy(tmp, helpfile);
#ifdef G_OS_WIN32
	p = strrchr(tmp, '_');
	if (p) strcpy(p, ".txt");
#else
	p = strrchr(tmp, '.');
	if (p) *p = 0;
#endif
	if (script) {
	    english_script_helpfile = tmp;
	} else {
	    english_gui_helpfile = tmp;
	}
	fp = fopen(tmp, "r");
	if (fp != NULL) {
	    char testline[MAXLEN];
	    int len;

	    len = 0;
	    while (fgets(testline, MAXLEN-1, fp)) {
		if (*testline != '@') len++;
	    }
	    fclose(fp);
	    if (script) english_script_help_length = len;
	    else english_gui_help_length = len;
	}
    }
}

static void set_translated_helpfile (void)
{
    char *p;

#ifdef G_OS_WIN32
    p = strrchr(paths.helpfile, '_');
    if (p && strncmp(p, "_hlp", 4)) { 
	translated_helpfile = 1;
    } else {
	translated_helpfile = 0;
    }
#else
    p = strrchr(paths.helpfile, '.');
    if (p && strcmp(p, ".hlp")) { 
	translated_helpfile = 1;
    } else {
	translated_helpfile = 0;
    }
#endif

    if (translated_helpfile == 1) {
	set_english_help_file(0);
	set_english_help_file(1);
    }
}
#endif

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

#ifdef ENABLE_NLS
    if (translated_helpfile < 0) 
	set_translated_helpfile();
#endif

    /* first pass: find length and number of topics */
    fp = fopen(helpfile, "r");
    if (fp == NULL) {
	fprintf(stderr, I_("help file %s is not accessible\n"), helpfile);
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
		if (t > 0) (heads[match])->topics[m] = t;
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

#ifdef ENABLE_NLS
static void english_help_callback (gpointer p, int script, 
				   GtkWidget *w)
{
    helpwin(script, 1);
}

static void add_english_help_item (windata_t *hwin, int script)
{
    GtkItemFactoryEntry helpitem;
    gchar mpath[] = "/_English";

    helpitem.accelerator = NULL;
    helpitem.callback_action = script; 
    helpitem.item_type = NULL;
    helpitem.path = mpath;
    helpitem.callback = english_help_callback; 
    gtk_item_factory_create_item(hwin->ifac, &helpitem, NULL, 1);
}
#endif

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
			gretl_commands[(heads[i])->topics[j]]);
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

static windata_t *helpwin (int script, int english) 
{
    windata_t *vwin = NULL;
    GtkItemFactoryEntry *items = help_items;
    char *helpfile;
    int helpcode;

    if (script) {
	helpfile = paths.cmd_helpfile;
	helpcode = CLI_HELP;
    } else {
	helpfile = paths.helpfile;
	helpcode = GUI_HELP;
    }

    if (helpfile == NULL) return NULL;

#ifdef ENABLE_NLS
    if (english) {
	if (script) {
	    helpcode = CLI_HELP_ENGLISH;
	    helpfile = english_script_helpfile;
	} else {
	    helpcode = GUI_HELP_ENGLISH;
	    helpfile = english_gui_helpfile;
	}
	items = english_help_items;
    }
#endif

    vwin = view_file(helpfile, 0, 0, 80, 400, helpcode, items);

    if (!english) 
	add_help_topics(vwin, script);

#ifdef ENABLE_NLS
    if (translated_helpfile && !english) 
	add_english_help_item(vwin, script);
#endif

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
    static GtkWidget *gui_help_view;
    static GtkWidget *script_help_view;
    GtkWidget *w = (cli)? script_help_view : gui_help_view;
    GtkTextBuffer *buf;
    GtkTextIter iter;
    GtkTextMark *vis;

    if (w == NULL) {
	windata_t *hwin = helpwin(cli, 0);

	if (hwin != NULL) {
	    if (cli) w = script_help_view = hwin->w;
	    else w = gui_help_view = hwin->w;
	}
	g_signal_connect(G_OBJECT(w), "destroy",
			 G_CALLBACK(gtk_widget_destroyed),
			 (cli)? &script_help_view : &gui_help_view);	
    } else {
	gdk_window_show(w->parent->window);
	gdk_window_raise(w->parent->window);
    }
    
    buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(w));
    gtk_text_buffer_get_iter_at_line_index(buf, &iter, pos, 0);
    vis = gtk_text_buffer_create_mark(buf, "vis", &iter, FALSE);
    gtk_text_view_scroll_to_mark(GTK_TEXT_VIEW(w), vis, 0.0, TRUE, 0.1, 0.0);    
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

gint edit_script_help (GtkWidget *widget, GdkEventButton *b,
		       windata_t *vwin)
{
    if (!vwin->help_active) { /* command help not activated */
	return FALSE;
    } else {
	gchar *text = NULL;
	GtkTextBuffer *buf;
	GtkTextIter iter;
	int pos = 0;

	buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->w));
	gtk_text_buffer_get_iter_at_mark(buf, &iter,
					 gtk_text_buffer_get_insert(buf));

	if (gtk_text_iter_inside_word(&iter)) {
	    GtkTextIter w_start, w_end;

	    w_start = iter;
	    w_end = iter;

	    if (!gtk_text_iter_starts_word(&iter)) {
		gtk_text_iter_backward_word_start(&w_start);
	    }
	    if (!gtk_text_iter_ends_word(&iter)) {
		gtk_text_iter_forward_word_end(&w_end);
	    }
	    text = gtk_text_buffer_get_text(buf, &w_start, &w_end, FALSE);
	} 

	if (text != NULL && strlen(text) > 0) {
	    char word[9];
	    int cmdnum;

	    *word = 0;
	    strncat(word, text, 8);
	    cmdnum = command_number(word);
	    pos = pos_from_cmd(cmdnum);
	} 
	
	real_do_help(pos, 1);
	g_free(text);
	text_set_cursor(vwin->w, 0);
	vwin->help_active = 0;
    }
    return FALSE;
}

/* ........................................................... */

void menu_find (gpointer data, guint db, GtkWidget *widget)
{
    if (db) 
	find_string_dialog(find_in_listbox, data);
    else 
	find_string_dialog(find_in_text, data);
}

/* ........................................................... */

void datafile_find (GtkWidget *widget, gpointer data)
{
    find_string_dialog(find_in_listbox, data);
}

/* ........................................................... */

void find_var (gpointer p, guint u, GtkWidget *w)
{
    find_string_dialog(find_in_listbox, mdata);
}

/* .................................................................. */

static gint close_find_dialog (GtkWidget *widget, gpointer data)
{
    find_window = NULL;
    return FALSE;
}

/* .................................................................. */

static gboolean real_find_in_text (GtkTextView *view, const gchar* str, 
				   gboolean from_cursor)
{
    GtkTextBuffer *buf;
    GtkTextIter iter;
    gboolean found = FALSE;
    gboolean wrapped = FALSE;
    GtkTextSearchFlags search_flags;

    buf = gtk_text_view_get_buffer (view);

    search_flags = GTK_TEXT_SEARCH_VISIBLE_ONLY | GTK_TEXT_SEARCH_TEXT_ONLY;

 text_search_wrap:
	
    if (from_cursor) {
	GtkTextIter sel_bound;
		
	gtk_text_buffer_get_iter_at_mark (buf,			
					  &iter,
					  gtk_text_buffer_get_mark (buf,
								    "insert"));
	gtk_text_buffer_get_iter_at_mark (buf,			
					  &sel_bound,
					  gtk_text_buffer_get_mark (buf,
								    "selection_bound"));
	gtk_text_iter_order (&sel_bound, &iter);		
    } else {		
	gtk_text_buffer_get_iter_at_offset (buf, &iter, 0);
    }

    if (*str != '\0') {
	GtkTextIter match_start, match_end;
	GtkTextMark *vis;

	found = gtk_text_iter_forward_search (&iter, str, search_flags,
					      &match_start, &match_end,
					      NULL);	
	if (found) {
	    gtk_text_buffer_place_cursor (buf, &match_start);
	    gtk_text_buffer_move_mark_by_name (buf, "selection_bound", &match_end);
	    vis = gtk_text_buffer_create_mark (buf, "vis", &match_end, FALSE);
	    gtk_text_view_scroll_to_mark (view, vis, 0.0, TRUE, 0.1, 0.0);
	} else if (from_cursor && !wrapped) {
	    /* try wrapping */
	    from_cursor = FALSE;
	    wrapped = TRUE;
	    goto text_search_wrap;
	}
    }

    if (found && wrapped) infobox(_("Search wrapped"));

    return found;
}

/* .................................................................. */

static void find_in_text (GtkWidget *widget, gpointer data)
{
    gboolean found;
    windata_t *vwin = 
	(windata_t *) g_object_get_data(G_OBJECT(data), "windat");

    needle = gtk_editable_get_chars(GTK_EDITABLE(find_entry), 0, -1);

    found = real_find_in_text(GTK_TEXT_VIEW(vwin->w), needle, TRUE);

    if (!found) infobox(_("String was not found."));
}

/* .................................................................. */

static void find_in_listbox (GtkWidget *w, gpointer data)
{
    int found = 0, wrapped = 0, minvar = 0;
    gchar *tmp, *pstr; 
    char haystack[MAXLEN];
    windata_t *win;
    GtkTreeModel *model;
    GtkTreeIter iter, iterhere;

    win = (windata_t *) g_object_get_data(G_OBJECT(data), "windat");
    if (win == mdata) {
	/* searching in the main gretl window: start on line 1, not 0 */
	minvar = 1;
    }

#if 0
    fprintf(stderr, "find_in_listbox: win at %p, active_var = %d\n", 
	    (void *) win, win->active_var);
#endif

    if (needle) g_free(needle);
    needle = gtk_editable_get_chars(GTK_EDITABLE(find_entry), 0, -1);
    lower(needle);

    model = gtk_tree_view_get_model (GTK_TREE_VIEW(win->listbox));

    /* try searching downward from the current line plus one */
    pstr = g_strdup_printf("%d", win->active_var);
    gtk_tree_model_get_iter_from_string (model, &iter, pstr);
    g_free(pstr);
    iterhere = iter;
    if (!gtk_tree_model_iter_next(model, &iter)) iter = iterhere;

 search_wrap:

    while (1) {
	/* try looking in column 1 first */
	gtk_tree_model_get(model, &iter, 1, &tmp, -1);
	strcpy(haystack, tmp);
	g_free(tmp);
	lower(haystack);
	found = look_for_string(haystack, needle, 0);
	if (found >= 0) break;
	/* then column 0 */
	gtk_tree_model_get (model, &iter, 0, &tmp, -1);
	strcpy(haystack, tmp);
	g_free(tmp);
	lower(haystack);
	found = look_for_string(haystack, needle, 0);
	if (found >= 0) break;
	if (!gtk_tree_model_iter_next(model, &iter)) break;
    }

    if (found < 0 && win->active_var > minvar && !wrapped) {
	/* try wrapping to start */
	gtk_tree_model_get_iter_first(model, &iter);
	if (minvar > 0 && !gtk_tree_model_iter_next(model, &iter)) {
	    ; /* do nothing: there's only one line in the box */
	} else {
	    wrapped = 1;
	    goto search_wrap;
	}
    }
    
    if (found >= 0) {
	GtkTreePath *path = gtk_tree_model_get_path(model, &iter);

	gtk_tree_view_scroll_to_cell(GTK_TREE_VIEW(win->listbox),
				     path, NULL, FALSE, 0, 0);
	gtk_tree_view_set_cursor(GTK_TREE_VIEW(win->listbox),
				 path, NULL, FALSE);
	win->active_var = tree_path_get_row_number(path);
	gtk_tree_path_free(path);
	if (wrapped) infobox(_("Search wrapped"));
    } else {
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
    GtkWidget *w = NULL;

    if (caller->dialog != NULL) {
	w = caller->dialog;
    } else if (caller->w != NULL) {
	w = caller->w;
    }

    if (w != NULL) {
	gtk_window_set_transient_for(GTK_WINDOW(finder),
				     GTK_WINDOW(w));
	gtk_window_set_destroy_with_parent(GTK_WINDOW(finder), TRUE);
    }
}

/* .................................................................. */

static void find_string_dialog (void (*findfunc)(), gpointer data)
{
    GtkWidget *label;
    GtkWidget *button;
    GtkWidget *hbox;
    windata_t *mydat = (windata_t *) data;

    if (find_window != NULL) {
	g_object_set_data(G_OBJECT(find_window), "windat", mydat);
	parent_find(find_window, mydat);
	return;
    }

    find_window = gtk_dialog_new();
    g_object_set_data(G_OBJECT(find_window), "windat", mydat);
    parent_find(find_window, mydat);

    g_signal_connect (G_OBJECT (find_window), "destroy",
		      G_CALLBACK (close_find_dialog),
		      find_window);

    gtk_window_set_title (GTK_WINDOW (find_window), _("gretl: find"));
    gtk_container_set_border_width (GTK_CONTAINER (find_window), 5);

    hbox = gtk_hbox_new(TRUE, TRUE);
    label = gtk_label_new(_(" Find what:"));
    gtk_widget_show (label);
    find_entry = gtk_entry_new();

    if (needle) {
	gtk_entry_set_text(GTK_ENTRY (find_entry), needle);
	gtk_editable_select_region (GTK_EDITABLE (find_entry), 0, 
				    strlen (needle));
    }
    g_signal_connect(G_OBJECT (find_entry), 
		     "activate", 
		     G_CALLBACK (findfunc),
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
    g_signal_connect(G_OBJECT (button), "clicked",
		     G_CALLBACK (findfunc), find_window);
    gtk_widget_grab_default(button);
    gtk_widget_show(button);

    /* cancel button */
    button = gtk_button_new_with_label (_("Cancel"));
    GTK_WIDGET_SET_FLAGS(button, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX (GTK_DIALOG (find_window)->action_area), 
		       button, TRUE, TRUE, FALSE);
    g_signal_connect(G_OBJECT (button), "clicked",
		     G_CALLBACK (cancel_find), find_window);
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

static GtkTooltips *gretl_tips;

void gretl_tooltips_init (void)
{
    gretl_tips = gtk_tooltips_new();
    gtk_tooltips_enable(gretl_tips); /* redundant? */
}

void gretl_tooltips_add (GtkWidget *w, const gchar *str)
{
    gtk_tooltips_set_tip(gretl_tips, w, str, NULL);
}
