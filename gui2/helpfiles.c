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
#include "textbuf.h"
#include "webget.h"

#ifdef OLD_GTK
# include "dlgutils.h"
#else
# include "treeutils.h"
#endif

#define HDEBUG 0

static int translated_helpfile = -1;
static char *en_gui_helpfile;
static char *en_cli_helpfile;

/* helpfile stuff */
struct help_head_t {
    char *name;
    int *topics;
    char **topicnames;
    int *pos;
    int ntopics;
};

static struct help_head_t **cli_heads, **gui_heads;
static struct help_head_t **en_cli_heads, **en_gui_heads;

static windata_t *helpwin (int cli, int en);
static void real_do_help (int hcode, int pos, int cli, int en);

/* searching stuff */
static void find_in_text (GtkWidget *widget, gpointer data);
static void find_in_listbox (GtkWidget *widget, gpointer data);
static void find_string_dialog (void (*findfunc)(), gpointer data);
#ifdef OLD_GTK
static int clist_start_row;
#endif

static GtkWidget *find_window = NULL;
static GtkWidget *find_entry;
static char *needle;

#ifndef OLD_GTK

GtkItemFactoryEntry help_menu_items[] = {
    { N_("/_Topics"), NULL, NULL, 0, "<Branch>", GNULL },    
    { N_("/_Find"), NULL, NULL, 0, "<Branch>", GNULL },   
    { N_("/Find/_Find in window"), NULL, menu_find, 0, "<StockItem>", GTK_STOCK_FIND },
    { NULL, NULL, NULL, 0, NULL, GNULL }
};

#else

GtkItemFactoryEntry help_menu_items[] = {
    { N_("/_Topics"), NULL, NULL, 0, "<Branch>" },    
    { N_("/_Find"), NULL, menu_find, 0, NULL },
    { NULL, NULL, NULL, 0, NULL}
};

#endif

struct gui_help_item {
    int code;
    char *string;
};

static struct gui_help_item gui_help_items[] = {
    { 0,              "nothing" },
    { GR_PLOT,        "graphing" },
    { GR_XY,          "graphing" },
    { GR_DUMMY,       "factorized" },
    { BXPLOT,         "boxplots" },
    { GR_BOX,         "boxplots" },
    { GR_NBOX,        "boxplots" },
    { GR_3D,          "3-D" },
    { ONLINE,         "online" },
    { MARKERS,        "markers" },
    { EXPORT,         "export" },
    { SMPLBOOL,       "sampling" },
    { SMPLDUM,        "sampling" },
    { COMPACT,        "compact" },
    { EXPAND,         "expand" },
    { VSETMISS,       "missing" },
    { GSETMISS,       "missing" },
    { GUI_HELP,       "dialog" },
    { MODELTABLE,     "modeltab" },
    { GRAPHPAGE,      "graphpag" },
    { SETSEED,        "seed" },
    { KERNEL_DENSITY, "density" },
    { HCCME,          "hccme" },
    { IRF_BOOT,       "irfboot" },
    { HTEST,          "gui-htest" },
    { MODEL_RESTR,    "restrict-model" },
    { SYS_RESTR,      "restrict-system" },
    { VECM_RESTR,     "restrict-vecm" },
    { -1,          NULL },
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

static int extra_command_number (const char *s)
{
    int i;

    for (i=1; gui_help_items[i].code > 0; i++) {
	if (!strcmp(s, gui_help_items[i].string)) {
	    return gui_help_items[i].code;
	}
    }

    return -1;
}

static char *help_string_from_cmd (int cmd)
{
    int i;

    for (i=1; gui_help_items[i].code > 0; i++) {
	if (cmd == gui_help_items[i].code) {
	    return gui_help_items[i].string;
	}
    }

    return NULL;    
}

static void set_en_help_file (int gui)
{
    char *helpfile, *tmp, *p;

    if (gui) {
	helpfile = paths.helpfile;
    } else {
	helpfile = paths.cmd_helpfile;
    }

    tmp = malloc(strlen(helpfile) + 1);

    if (tmp != NULL) {
	strcpy(tmp, helpfile);
#ifdef G_OS_WIN32
	p = strrchr(tmp, '_');
	if (p != NULL) strcpy(p, ".txt");
#else
	p = strrchr(tmp, '.');
	if (p != NULL) *p = 0;
#endif
	if (gui) {
	    en_gui_helpfile = tmp;
	} else {
	    en_cli_helpfile = tmp;
	}
    }
}

static void set_translated_helpfile (void)
{
    char *p;

#ifdef G_OS_WIN32
    p = strrchr(paths.helpfile, '_');
    if (p != NULL && strncmp(p, "_hlp", 4)) { 
	translated_helpfile = 1;
    } else {
	translated_helpfile = 0;
    }
#else
    p = strrchr(paths.helpfile, '.');
    if (p != NULL && strcmp(p, ".hlp")) { 
	translated_helpfile = 1;
    } else {
	translated_helpfile = 0;
    }
#endif

    if (translated_helpfile == 1) {
	set_en_help_file(0);
	set_en_help_file(1);
    }
}

static int match_heading (struct help_head_t **heads, int nh,
			  const char *str)
{
    int i, match = -1;

    if (heads == NULL) return -1;

    for (i=0; i<nh; i++) {
	if (!strcmp(str, heads[i]->name)) {
#if HDEBUG
	    fprintf(stderr, "str='%s', heads[%d].name='%s', matched\n",
		    str, i, heads[i]->name);
#endif
	    match = i;
	    break;
	}
    }

    return match;
}

static int add_help_heading (struct help_head_t ***pheads, 
			     int *pnh, const char *str)
{
    struct help_head_t **heads;
    int nh = *pnh;
    int err = 0;

    heads = realloc(*pheads, (nh + 2) * sizeof *heads);
    if (heads == NULL) return 1;

    heads[nh] = malloc(sizeof **heads);

    if (heads[nh] != NULL) {
	heads[nh]->name = malloc(strlen(str) + 1);
	if (heads[nh]->name != NULL) {
	    strcpy(heads[nh]->name, str);
#if HDEBUG
	    fprintf(stderr, "str='%s', heads[%d].name added new\n", str, nh);
#endif
	    heads[nh]->ntopics = 1;
	    nh++;
	} else err = 1;
    } else err = 1;

    if (!err) {
	*pheads = heads;
	*pnh = nh;
    }

    return err;
}

static int allocate_heads_info (struct help_head_t **heads, int nh, int gui)
{
    int *topics = NULL, *pos = NULL;
    char **topicnames = NULL;
    int i, nt, err = 0;

    for (i=0; i<nh && !err; i++) {
	nt = heads[i]->ntopics;
	if (nt == 0) continue;

	topics = malloc(nt * sizeof *topics);
	if (topics == NULL) err = 1;

	pos = malloc(nt * sizeof *pos);
	if (pos == NULL) err = 1; 

	if (gui) {
	    topicnames = malloc(nt * sizeof *topicnames);
	    if (topicnames == NULL) err = 1;
	} 

	if (!err) {
	    heads[i]->topics = topics;
	    heads[i]->topicnames = topicnames;
	    heads[i]->pos = pos;
	}

	heads[i]->ntopics = 0;
    }

    /* sentinel */
    heads[i] = NULL;

    return err;
}

char *quoted_help_string (const char *s)
{
    const char *p, *q;

    p = strchr(s, '"');
    q = strrchr(s, '"');

    if (p != NULL && q != NULL && q - p > 1) {
	p++;
	return g_strndup(p, q - p);
    }

    return g_strdup("Missing string");
}

static int add_topic_to_heading (struct help_head_t **heads, 
				 int nh, char *s, int pos, 
				 int gui)
{
    char word[32], section[32];
    int n, m, nt;

    if (sscanf(s + 1, "%31s %31s", word, section) != 2) {
	return 1;
    }

    if (!strcmp(section, "Obsolete")) {
	return 0;
    }

#if HDEBUG
    fprintf(stderr, "add_topic_to_heading: '%s' (%s), pos=%d\n", 
	    word, section, pos);
#endif

    m = match_heading(heads, nh, section);
    nt = heads[m]->ntopics;

    n = gretl_command_number(word);
    if (n <= 0) {
	n = extra_command_number(word);
    }

    if (n > 0) {
	heads[m]->topics[nt] = n;
    } else {
	return 1;
    }

    if (gui) {
	heads[m]->topicnames[nt] = quoted_help_string(s);
#if HDEBUG > 1
	fprintf(stderr, "Set heads[%d]->topicnames[%d] = \n"
		"  quoted_help_string(%s) = '%s'\n",
		m, nt, s, heads[m]->topicnames[nt]);
#endif
    }

    heads[m]->pos[nt] = pos;
    heads[m]->ntopics += 1;

    return 0;
}

static int assemble_topic_list (struct help_head_t **heads, int nh,
				int gui, FILE *fp)
{
    char test[256];
    int pos = 0;

#if HDEBUG
    fprintf(stderr, "\n*** starting assemble_topic_list\n");
#endif

    while (fgets(test, sizeof test, fp)) {
	if (*test == '#') {
	    add_topic_to_heading(heads, nh, test, pos, gui);
	} 
	pos += strlen(test);
    }

    return 0;
}

static int 
get_help_topics (struct help_head_t ***pheads, int *pnh, FILE *fp)
{
    char test[256], section[16];    
    int nh = 0;
    int match, err = 0;

    while (!err && fgets(test, sizeof test, fp)) {
	if (*test == '#') {
	    sscanf(test + 1, "%*s %15s", section);
	    if (strcmp(section, "Obsolete")) {
		match = match_heading(*pheads, nh, section);
		if (match >= 0) {
		    (*pheads)[match]->ntopics += 1;
		} else {
		    err = add_help_heading(pheads, &nh, section);
		}
	    }
	} 
    }

    *pnh = nh;

    return err;
}

#if 0

static int get_byte_positions (struct help_head_t **heads, int nh, FILE *fp)
{
    char test[9];
    int n, c, cbak = 0;
    int found;
    int i, j, k, b = 0;

    for (i=0; i<nh; i++) {
	for (j=0; j<heads[i]->ntopics; j++) {
	    found = 0;
	    while (!found) {
		c = getc(fp);
		if (c == EOF) {
		    return 1;
		}
		if (c == '#' && cbak == '\n') {
		    
		    heads[i]->pos[j] = b;
		    found = 1;
		}
		cbak = c;
		b++;
	    }
	}
    }

    return 0;
}

int get_helpfile_structure (int code)
{
    char line[64];
    int nh;
    FILE *fp = NULL;

    fp = fopen(fname, "r");
    if (fp == NULL) {
	return 1;
    }

    if (fgets(line, sizeof line, fp) == NULL) {
	goto bailout;
    }

    if (!sscanf(line, "headings %d", &nh) || nh <= 0) {
	goto bailout;
    }

    heads = malloc(nh * sizeof *heads);
    if (heads == NULL) {
	goto bailout;
    }

    for (i=0; i<nh; i++) {
	if (fgets(line, sizeof line, fp) == NULL) {
	    goto bailout;
	}

    }	


 bailout:

    fclose(fp);

    return 1;
}

#endif

static int real_helpfile_init (int gui, int en)
{
    struct help_head_t **heads = NULL;
    char *helpfile;
    FILE *fp;
    int nh = 0;
    int err = 0;

    if (en) {
	helpfile = (gui)? en_gui_helpfile : en_cli_helpfile;
    } else {
	helpfile = (gui)? paths.helpfile : paths.cmd_helpfile;
    }

#if HDEBUG
    fprintf(stderr, "\n*** real_helpfile_init: gui = %d, en = %d\n"
	    " helpfile = '%s'\n", gui, en, helpfile);
    fflush(stderr);
#endif

    fp = gretl_fopen(helpfile, "r");
    if (fp == NULL) {
	fprintf(stderr, I_("help file %s is not accessible\n"), helpfile);
	return 1;
    }

    /* first pass: find number of topics */
    err = get_help_topics(&heads, &nh, fp);

#if HDEBUG
    fprintf(stderr, "real_helpfile_init: found %d headings\n", nh);
#endif

    if (!err) {
	err = allocate_heads_info(heads, nh, gui);
    }

    if (err) {
	fclose(fp);
	return err;
    }

#if HDEBUG
    fprintf(stderr, "real_helpfile_init: done allocate_heads_info\n");
#endif

    /* second pass, assemble the topic list */
    rewind(fp);
    assemble_topic_list(heads, nh, gui, fp);
    fclose(fp);

    if (gui) {
	if (en) {
	    en_gui_heads = heads;
	} else {
	    gui_heads = heads;
	}
    } else {
	if (en) {
	    en_cli_heads = heads;
	} else {
	    cli_heads = heads;
	}
    }

    return err;
}

void helpfile_init (void)
{
    if (translated_helpfile < 0) { 
	set_translated_helpfile();
    }

    real_helpfile_init(0, 0);
    real_helpfile_init(1, 0);

    if (translated_helpfile == 1) { 
	real_helpfile_init(0, 1);
	real_helpfile_init(1, 1);
    }
}

static void do_gui_help (gpointer p, guint pos, GtkWidget *w) 
{
    int hcode = GPOINTER_TO_INT(p);

    real_do_help(hcode, pos, 0, 0);
}

static void do_en_gui_help (gpointer p, guint pos, GtkWidget *w) 
{
    int hcode = GPOINTER_TO_INT(p);

    real_do_help(hcode, pos, 0, 1);
}

static void do_cli_help (gpointer p, guint pos, GtkWidget *w) 
{
    int hcode = GPOINTER_TO_INT(p);

    real_do_help(hcode, pos, 1, 0);
}

static void do_en_cli_help (gpointer p, guint pos, GtkWidget *w) 
{
    int hcode = GPOINTER_TO_INT(p);

    real_do_help(hcode, pos, 1, 1);
}

static char *get_gui_help_string (int pos)
{
    int i, j;

    for (i=0; gui_heads[i] != NULL; i++) { 
	for (j=0; j<gui_heads[i]->ntopics; j++) {
	    if (pos == gui_heads[i]->pos[j]) {
		return help_string_from_cmd(gui_heads[i]->topics[j]);
	    }
	}
    }

    return NULL;
}

static int gui_pos_from_cmd (int cmd, int en)
{
    struct help_head_t **heads;
    char *helpstr;
    int i, j;

    heads = (en)? en_gui_heads : gui_heads;

    if (heads == NULL) {
	return -1;
    }

    for (i=0; heads[i] != NULL; i++) {
	for (j=0; j<heads[i]->ntopics; j++) {
	    if (cmd == heads[i]->topics[j]) {
		return heads[i]->pos[j];
	    }
	}
    }

    /* special for gui-specific help items */
    helpstr = help_string_from_cmd(cmd);
    if (helpstr != NULL) {
	int altcode = extra_command_number(helpstr);

	for (i=0; heads[i] != NULL; i++) {
	    for (j=0; j<heads[i]->ntopics; j++) {
		if (altcode == heads[i]->topics[j]) {
		    return heads[i]->pos[j];
		}
	    }
	}
    }
    
    return -1;
}

static int cli_pos_from_cmd (int cmd, int en)
{
    struct help_head_t **heads;
    int i, j;

    heads = (en)? en_cli_heads : cli_heads;

    if (heads == NULL) {
	return -1;
    }
    
    for (i=0; heads[i] != NULL; i++) {
	for (j=0; j<heads[i]->ntopics; j++) {
	    if (cmd == heads[i]->topics[j]) {
		return heads[i]->pos[j];
	    }
	}
    }

    return -1;
}

static void en_help_callback (gpointer p, int cli, GtkWidget *w)
{
    windata_t *hwin = (windata_t *) p;
    int pos, hc = hwin->active_var;

    if (cli) {
	pos = cli_pos_from_cmd(hc, 1);
	if (pos < 0) {
	    pos = 0;
	}
    } else {
	pos = gui_pos_from_cmd(hc, 1);
    }

    real_do_help(hc, pos, cli, 1);
}

static void add_en_help_item (windata_t *hwin, int cli)
{
    GtkItemFactoryEntry item;

    item.accelerator = NULL;

    item.item_type = "<Branch>";
    item.callback = NULL;
    item.callback_action = 0; 
    item.path = "/_English";
    gtk_item_factory_create_item(hwin->ifac, &item, NULL, 1);

    item.item_type = NULL;
    item.callback = en_help_callback; 
    item.callback_action = cli; 
    item.path = g_strdup_printf("/English/%s", _("Show English help"));
    gtk_item_factory_create_item(hwin->ifac, &item, hwin, 1);
    g_free(item.path);
}

static void add_help_topics (windata_t *hwin, int cli, int en)
{
    int i, j;
    GtkItemFactoryEntry hitem;
    const gchar *mpath = N_("/_Topics");
    struct help_head_t **hds;

    if (en) {
	hds = (cli)? en_cli_heads : en_gui_heads;
    } else {
	hds = (cli)? cli_heads : gui_heads;
    }

    hitem.accelerator = NULL;

    /* See if there are any topics to add */
    if (hds == NULL) return;

#ifndef OLD_GTK
    if (cli) {
	/* Add general index as "topic" */
	hitem.callback_action = 1; 
	hitem.item_type = NULL;
	hitem.path = g_strdup_printf("%s/%s", mpath, _("Index"));
	hitem.callback = (en)? do_en_cli_help : do_cli_help;
	gtk_item_factory_create_item(hwin->ifac, &hitem, 
				     GINT_TO_POINTER(0), 
				     1);
	g_free(hitem.path);
    }
#endif

    /* put the topics under the menu heading */
    for (i=0; hds[i] != NULL; i++) {

	if (hds[i]->ntopics == 0) continue;

	hitem.callback_action = 0; 
	hitem.item_type = "<Branch>";
	hitem.path = g_strdup_printf("%s/%s", mpath, _(hds[i]->name));
	hitem.callback = NULL; 

	gtk_item_factory_create_item(hwin->ifac, &hitem, NULL, 1);
	g_free(hitem.path);

	for (j=0; j<hds[i]->ntopics; j++) {
	    int tnum = -1;
	    int topic_ok = 1;

	    hitem.callback_action = hds[i]->pos[j]; 
	    hitem.item_type = NULL;
	    hitem.path = NULL;

	    tnum = hds[i]->topics[j];

	    if (hds[i]->topicnames != NULL) {
		hitem.path = 
		    g_strdup_printf("%s/%s/%s", 
				    mpath, _(hds[i]->name), 
				    hds[i]->topicnames[j]);
#if HDEBUG
		fprintf(stderr, "(1) Built help topic path from\n"
			" '%s', '%s' and '%s'\n", mpath, _(hds[i]->name),
			hds[i]->topicnames[j]);
#endif
	    } else {
		if (tnum < NC) {
		    /* a regular gretl command */
		    hitem.path = 
			g_strdup_printf("%s/%s/%s", 
					mpath, _(hds[i]->name), 
					gretl_command_word(tnum));
#if HDEBUG
		    fprintf(stderr, "(2) Built help topic path from\n"
			    " '%s', '%s' and '%s'\n", mpath, _(hds[i]->name),
			    gretl_command_word(tnum));
#endif
		} else if (!cli) {
		    /* a gui special item? */
		    char *gstr = get_gui_help_string(hds[i]->pos[j]);

		    if (gstr != NULL) {
			hitem.path = 
			    g_strdup_printf("%s/%s/%s", 
					    mpath, _(hds[i]->name), gstr);
#if HDEBUG
			fprintf(stderr, "(3) Built help topic path from\n"
				" '%s', '%s' and '%s'\n", mpath, _(hds[i]->name),
				gstr);
#endif
		    } else {
			topic_ok = 0;
		    }
		} else {
		    topic_ok = 0;
		}
	    }

	    if (topic_ok) {
		if (en) {
		    hitem.callback = (cli)? do_en_cli_help : do_en_gui_help; 
		} else {
		    hitem.callback = (cli)? do_cli_help : do_gui_help; 
		}
		gtk_item_factory_create_item(hwin->ifac, &hitem, 
					     GINT_TO_POINTER(tnum), 
					     1);
		g_free(hitem.path);
	    }
	}
    }
}

static windata_t *helpwin (int cli, int en) 
{
    windata_t *vwin = NULL;
    char *helpfile = NULL;
    int role;

    if (cli) {
	helpfile = (en)? en_cli_helpfile : paths.cmd_helpfile;
	role = (en)? CLI_HELP_EN : CLI_HELP;
    } else {
	helpfile = (en)? en_gui_helpfile : paths.helpfile;
	role = (en)? GUI_HELP_EN : GUI_HELP;
    }

    if (helpfile == NULL) return NULL;

    vwin = view_help_file(helpfile, role, help_menu_items);

    add_help_topics(vwin, cli, en);

    if (translated_helpfile && !en) {
	add_en_help_item(vwin, cli);
    }    

    return vwin;
}

void context_help (GtkWidget *widget, gpointer data)
{
    int pos, hcode = GPOINTER_TO_INT(data);
    int en = 0;

    pos = gui_pos_from_cmd(hcode, 0);

    if (pos < 0 && translated_helpfile) {
	en = 1;
	pos = gui_pos_from_cmd(hcode, en);
    }

#if HDEBUG
    fprintf(stderr, "context_help: hcode=%d, pos=%d, en=%d\n", 
	    hcode, pos, en);
#endif

    real_do_help(hcode, pos, 0, en);
}

static gboolean nullify_hwin (GtkWidget *w, windata_t **phwin)
{
    *phwin = NULL;
    return FALSE;
}

static void real_do_help (int hcode, int pos, int cli, int en)
{
    static windata_t *gui_hwin;
    static windata_t *cli_hwin;
    static windata_t *en_gui_hwin;
    static windata_t *en_cli_hwin;

    windata_t *hwin;

    if (pos < 0) {
	dummy_call();
	return;
    }

#if HDEBUG
    fprintf(stderr, "real_do_help: hcode=%d, pos=%d, cli=%d, en=%d\n",
	    hcode, pos, cli, en);
    fprintf(stderr, "gui_hwin = %p\n", (void *) gui_hwin);
    fprintf(stderr, "cli_hwin = %p\n", (void *) cli_hwin);
    fprintf(stderr, "en_gui_hwin = %p\n", (void *) en_gui_hwin);
    fprintf(stderr, "en_cli_hwin = %p\n", (void *) en_cli_hwin);
#endif

    if (en) {
	hwin = (cli)? en_cli_hwin : en_gui_hwin;
    } else {
	hwin = (cli)? cli_hwin : gui_hwin;
    }

    if (hwin == NULL) {
	hwin = helpwin(cli, en);
	if (hwin != NULL) {
	    windata_t **phwin;

	    if (en) {
		if (cli) {
		    en_cli_hwin = hwin;
		    phwin = &en_cli_hwin;
		} else {
		    en_gui_hwin = hwin;
		    phwin = &en_gui_hwin;
		}
	    } else {
		if (cli) {
		    cli_hwin = hwin;
		    phwin = &cli_hwin;
		} else {
		    gui_hwin = hwin;
		    phwin = &gui_hwin;
		}
	    }		
	    g_signal_connect(G_OBJECT(hwin->w), "destroy",
			     G_CALLBACK(nullify_hwin), phwin);
	}
    } else {
	gdk_window_show(hwin->w->parent->window);
	gdk_window_raise(hwin->w->parent->window);
    }

#if HDEBUG
    fprintf(stderr, "real_do_help: doing set_help_topic_buffer:\n"
	    " hwin=%p, hcode=%d, pos=%d, en=%d\n",
	    (void *) hwin, hcode, pos, en);
#endif

    set_help_topic_buffer(hwin, hcode, pos, en);
}

/* called from main menu in gretl.c; also used as callback from
   command ref index in textbuf.c (gtk2 version only)
*/

void plain_text_cmdref (gpointer p, guint cmdnum, GtkWidget *w)
{
    /* pos = 0 gives index of commands */
    int pos = 0;
    int en = 0;
    int cli = 1;

    if (w == NULL && p != NULL) {
	en = GPOINTER_TO_INT(p);
    }

    if (cmdnum > 0) {
	pos = cli_pos_from_cmd(cmdnum, en);
	if (pos < 0 && !en && translated_helpfile == 1) {
	    en = 1;
	    pos = cli_pos_from_cmd(cmdnum, en);
	}
    }

    real_do_help(cmdnum, pos, cli, en);
} 

#ifndef OLD_GTK

gint edit_script_help (GtkWidget *widget, GdkEventButton *b,
		       windata_t *vwin)
{
    if (!window_help_is_active(vwin)) { 
	/* command help not activated */
	return FALSE;
    } else {
	gchar *text = NULL;
	int pos = -1;
	int hcode = 0;
	int en = 0;
	GtkTextBuffer *buf;
	GtkTextIter iter;

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

	if (text != NULL && *text != '\0') {
	    char word[9];

	    *word = '\0';
	    strncat(word, text, 8);
	    hcode = gretl_command_number(word);
	    pos = cli_pos_from_cmd(hcode, en);
	    if (pos < 0 && translated_helpfile) {
		en = 1;
		pos = cli_pos_from_cmd(hcode, en);
	    }
	} 

	g_free(text);
	unset_window_help_active(vwin);
	text_set_cursor(vwin->w, 0);

	real_do_help(hcode, pos, 1, en);
    }

    return FALSE;
}

#else /* now old gtk */

gint edit_script_help (GtkWidget *widget, GdkEventButton *b,
		       windata_t *vwin)
{
    if (!window_help_is_active(vwin)) { 
	/* command help not activated */
	return FALSE;
    } else {
	gchar *text = NULL;
	int pos = -1;
	int hcode = 0;
	int en = 0;
	int pt = GTK_EDITABLE(vwin->w)->current_pos;
	int len = gtk_text_get_length(GTK_TEXT(vwin->w));

	text = gtk_editable_get_chars(GTK_EDITABLE(vwin->w), 
				      0, (pt + 9 > len)? -1 : pt + 8);

	if (text != NULL && *text != '\0') {
	    char *p, *q;
	    char word[9];

	    p = q = text + pt;
	    if (pt > 0) {
		while (p - text && !isspace(*(p-1))) p--;
	    }
	    if (pt < (int) strlen(text)) {
		while (*q && !isspace(*q)) q++;
	    }
	    *word = '\0';
	    strncat(word, p, (q - p > 8)? 8 : q - p);
	    hcode = gretl_command_number(word);
	    pos = cli_pos_from_cmd(hcode, en);
	    if (pos < 0 && translated_helpfile) {
		en = 1;
		pos = cli_pos_from_cmd(hcode, en);
	    }
	} 

	g_free(text);
	unset_window_help_active(vwin);
	gdk_window_set_cursor(GTK_TEXT(vwin->w)->text_area, NULL);

	real_do_help(hcode, pos, 1, en);
    }

    return FALSE;
}

#endif /* gtk versions */

void menu_find (gpointer data, guint db, GtkWidget *widget)
{
    if (db) {
	find_string_dialog(find_in_listbox, data);
    } else {
	find_string_dialog(find_in_text, data);
    }
}

void datafile_find (GtkWidget *widget, gpointer data)
{
    find_string_dialog(find_in_listbox, data);
}

void find_var (gpointer p, guint u, GtkWidget *w)
{
    find_string_dialog(find_in_listbox, mdata);
}

static gint close_find_dialog (GtkWidget *widget, gpointer data)
{
    find_window = NULL;
    return FALSE;
}

static int look_for_string (const char *haystack, const char *needle, 
			    int start)
{
    int hlen = strlen(haystack);
    int nlen = strlen(needle);
    int pos;

    for (pos = start; pos < hlen; pos++) {
        if (strncmp(&haystack[pos], needle, nlen) == 0) { 
             return pos;
	}
    }

    return -1;
}

#ifdef OLD_GTK

static int is_all_lower (const char *s)
{
    int ret = 1;

    while (*s) {
	if (!islower(*s)) {
	    ret = 0;
	    break;
	}
	s++;
    }

    return ret;
}

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

    if (is_all_lower(needle)) {
	lower(haystack);
    }

    found = look_for_string(haystack, needle, found);

    if (found >= 0) {
	gtk_text_set_point(GTK_TEXT(vwin->w), found);
	gtk_editable_set_position(GTK_EDITABLE(vwin->w), found);
        gtk_editable_select_region(GTK_EDITABLE(vwin->w), 
				   found, found + strlen(needle));
    } else {
	infobox(_("String was not found."));
    }

    g_free(haystack);
}

#else /* !OLD_GTK */

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
		
	gtk_text_buffer_get_iter_at_mark(buf,			
					 &iter,
					 gtk_text_buffer_get_mark(buf,
								  "insert"));
	gtk_text_buffer_get_iter_at_mark(buf,			
					 &sel_bound,
					 gtk_text_buffer_get_mark(buf,
								  "selection_bound"));
	gtk_text_iter_order(&sel_bound, &iter);		
    } else {		
	gtk_text_buffer_get_iter_at_offset(buf, &iter, 0);
    }

    if (*str != '\0') {
	GtkTextIter match_start, match_end;
	GtkTextMark *vis;

	found = gtk_text_iter_forward_search(&iter, str, search_flags,
					     &match_start, &match_end,
					     NULL);	
	if (found) {
	    gtk_text_buffer_place_cursor(buf, &match_start);
	    gtk_text_buffer_move_mark_by_name(buf, "selection_bound", &match_end);
	    vis = gtk_text_buffer_create_mark(buf, "vis", &match_end, FALSE);
	    gtk_text_view_scroll_to_mark(view, vis, 0.0, TRUE, 0.1, 0.0);
	} else if (from_cursor && !wrapped) {
	    /* try wrapping */
	    from_cursor = FALSE;
	    wrapped = TRUE;
	    goto text_search_wrap;
	}
    }

    if (found && wrapped) {
	infobox(_("Search wrapped"));
    }

    return found;
}

static void find_in_text (GtkWidget *widget, gpointer data)
{
    gboolean found;
    windata_t *vwin = 
	(windata_t *) g_object_get_data(G_OBJECT(data), "windat");

    needle = gtk_editable_get_chars(GTK_EDITABLE(find_entry), 0, -1);

    found = real_find_in_text(GTK_TEXT_VIEW(vwin->w), needle, TRUE);

    if (!found) {
	infobox(_("String was not found."));
    }
}

#endif /* new vs old gtk */

static void find_in_listbox (GtkWidget *w, gpointer data)
{
    int found = 0, wrapped = 0, minvar = 0;
    gchar *tmp; 
    char haystack[MAXLEN];
    windata_t *win;
#ifndef OLD_GTK
    gchar *pstr;
    GtkTreeModel *model;
    GtkTreeIter iter, iterhere;
#else
    int end, i;
#endif

    win = (windata_t *) g_object_get_data(G_OBJECT(data), "windat");
    if (win == mdata) {
	/* searching in the main gretl window: start on line 1, not 0 */
	minvar = 1;
    }

    if (needle) g_free(needle);
    needle = gtk_editable_get_chars(GTK_EDITABLE(find_entry), 0, -1);
    lower(needle);

#ifndef OLD_GTK
    model = gtk_tree_view_get_model(GTK_TREE_VIEW(win->listbox));

    /* try searching downward from the current line plus one */
    pstr = g_strdup_printf("%d", win->active_var);
    gtk_tree_model_get_iter_from_string (model, &iter, pstr);
    g_free(pstr);
    iterhere = iter;
    if (!gtk_tree_model_iter_next(model, &iter)) {
	iter = iterhere;
    }

 search_wrap:

    while (1) {
	/* try looking in column 1 first */
	gtk_tree_model_get(model, &iter, 1, &tmp, -1);
	strcpy(haystack, tmp);
	g_free(tmp);
	lower(haystack);
	found = look_for_string(haystack, needle, 0);
	if (found >= 0) break;
	if (win == mdata) {
	    /* then column 2 */
	    gtk_tree_model_get(model, &iter, 2, &tmp, -1);
	} else {
	    /* then column 0 */
	    gtk_tree_model_get(model, &iter, 0, &tmp, -1);
	}
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

#else /* GTK version switch */

    end = GTK_CLIST(win->listbox)->rows;

 search_wrap: 

    for (i=clist_start_row; i<end; i++) {  
	/* try looking in column 1 first */
	gtk_clist_get_text(GTK_CLIST(win->listbox), i, 1, &tmp);
	strcpy(haystack, tmp);
	lower(haystack);
	found = look_for_string(haystack, needle, 0);
	if (found >= 0) break;
	if (win == mdata) {
	    /* try column 2? */
	    gtk_clist_get_text(GTK_CLIST(win->listbox), i, 2, &tmp);
	} else {
	    /* try column 0? */
	    gtk_clist_get_text(GTK_CLIST(win->listbox), i, 0, &tmp);
	}
	strcpy(haystack, tmp);
	lower(haystack);
	found = look_for_string(haystack, needle, 0);
	if (found >= 0) break;
    }

    if (found < 0 && win->active_var > minvar && !wrapped) {
	/* try wrapping to start */
	end = win->active_var;
	clist_start_row = minvar;
	wrapped = 1;
	goto search_wrap;
    }    

    if (found >= 0) {
	if (wrapped) infobox(_("Search wrapped"));
	gtk_clist_freeze(GTK_CLIST(win->listbox));
	gtk_clist_moveto(GTK_CLIST(win->listbox), i, 0, 0, .1);
	gtk_clist_unselect_all(GTK_CLIST(win->listbox));
	gtk_clist_select_row(GTK_CLIST(win->listbox), i, 0);
	GTK_CLIST(win->listbox)->focus_row = i;
	gtk_clist_thaw(GTK_CLIST(win->listbox));
	win->active_var = i;
	clist_start_row = i + 1;
    } else {
	infobox(_("String was not found."));
    }
#endif /* OLD_GTK */
}

static void cancel_find (GtkWidget *widget, gpointer data)
{
    if (find_window != NULL) {
	gtk_widget_destroy(GTK_WIDGET(data));
	find_window = NULL;
    }
}

static void parent_find (GtkWidget *finder, windata_t *caller)
{
#ifndef OLD_GTK
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
#else
    if (caller != mdata) {
	GtkWidget *w = NULL;

	if (caller->dialog != NULL) {
	    w = caller->dialog;
	} else if (caller->w != NULL) {
	    w = caller->w;
	}

	if (w != NULL) {
	    gtk_window_set_transient_for(GTK_WINDOW(finder),
					 GTK_WINDOW(w));
	    gtk_signal_connect(GTK_OBJECT(w), "destroy",
			       GTK_SIGNAL_FUNC(cancel_find),
			       finder);
	}
    }
#endif
}

static void find_string_dialog (void (*findfunc)(), gpointer data)
{
    GtkWidget *label;
    GtkWidget *button;
    GtkWidget *hbox;
    windata_t *mydat = (windata_t *) data;

#ifdef OLD_GTK
    clist_start_row = 0;
#endif

    if (find_window != NULL) {
	g_object_set_data(G_OBJECT(find_window), "windat", mydat);
	parent_find(find_window, mydat);
	gdk_window_raise(find_window->window);
	return;
    }

    find_window = gtk_dialog_new();
    g_object_set_data(G_OBJECT(find_window), "windat", mydat);
    parent_find(find_window, mydat);

    g_signal_connect(G_OBJECT(find_window), "destroy",
		     G_CALLBACK(close_find_dialog),
		     find_window);

    gtk_window_set_title(GTK_WINDOW(find_window), _("gretl: find"));
    gtk_container_set_border_width(GTK_CONTAINER(find_window), 5);

    hbox = gtk_hbox_new(TRUE, TRUE);
    label = gtk_label_new(_(" Find what:"));
    gtk_widget_show (label);
    find_entry = gtk_entry_new();

    if (needle) {
	gtk_entry_set_text(GTK_ENTRY(find_entry), needle);
	gtk_editable_select_region(GTK_EDITABLE(find_entry), 0, 
				   strlen(needle));
    }

    g_signal_connect(G_OBJECT(find_entry), "activate", 
		     G_CALLBACK(findfunc), find_window);

    gtk_widget_show(find_entry);

    gtk_box_pack_start(GTK_BOX(hbox), label, TRUE, TRUE, 0);
    gtk_box_pack_start(GTK_BOX(hbox), find_entry, TRUE, TRUE, 0);
    gtk_widget_show(hbox);

    gtk_box_pack_start(GTK_BOX (GTK_DIALOG(find_window)->vbox), 
		       hbox, TRUE, TRUE, 0);

    gtk_box_set_spacing(GTK_BOX(GTK_DIALOG(find_window)->action_area), 15);
    gtk_box_set_homogeneous(GTK_BOX 
			    (GTK_DIALOG(find_window)->action_area), TRUE);
    gtk_window_set_position(GTK_WINDOW(find_window), GTK_WIN_POS_MOUSE);

    /* find button -- make this the default */
    button = standard_button(GTK_STOCK_FIND);
    GTK_WIDGET_SET_FLAGS(button, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(find_window)->action_area), 
		       button, TRUE, TRUE, FALSE);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(findfunc), find_window);
    gtk_widget_grab_default(button);
    gtk_widget_show(button);

    /* cancel button */
    button = standard_button(GTK_STOCK_CANCEL);
    GTK_WIDGET_SET_FLAGS(button, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(find_window)->action_area), 
		       button, TRUE, TRUE, FALSE);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(cancel_find), find_window);
    gtk_widget_show(button);

    gtk_widget_grab_focus(find_entry);
    gtk_widget_show(find_window);
}

void text_find_callback (GtkWidget *w, gpointer data)
{
    find_string_dialog(find_in_text, data);
}

#ifdef OLD_GTK

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

#endif

static GtkTooltips *gretl_tips;

void gretl_tooltips_init (void)
{
    gretl_tips = gtk_tooltips_new();
    gtk_tooltips_enable(gretl_tips); /* redundant? */
#ifdef OLD_GTK
    colorize_tooltips(gretl_tips); 
#endif
}

void gretl_tooltips_add (GtkWidget *w, const gchar *str)
{
    gtk_tooltips_set_tip(gretl_tips, w, str, NULL);
}

#ifdef OSX_PKG
static void osx_help (int uguide)
{
    char *prefix, *syscmd;
   
    prefix = getenv("GTK_EXE_PREFIX");
    if (prefix == NULL) {
        errbox("Couldn't find the manual");
	return;
    }
    
    syscmd = g_strdup_printf("%s/bin/manual.sh %s", prefix, (uguide)? "uguide" : "ref");
    system(syscmd);
    g_free(syscmd);
}
#endif 

enum {
    EN_LETTER,
    EN_A4,
    ITALIAN,
    SPANISH
};

static char *full_doc_path (char *path, const char *fname)
{
    strcpy(path, paths.gretldir);
    strcat(path, "doc");
    strcat(path, SLASHSTR);
    strcat(path, fname);

    return path;
}

static int maybe_grab_pdf (int uguide, int i, char *fullpath)
{
    const char *guide_files[] = {
	"gretl-guide.pdf",
	"gretl-guide-a4.pdf",
	"gretl-guide-it.pdf",
	"gretl-guide-es.pdf"
    };
    const char *ref_files[] = {
	"gretl-ref.pdf",
	"gretl-ref-a4.pdf",
	"gretl-ref-it.pdf",
	"gretl-ref-es.pdf"  
    };
    const char *fname;
    FILE *fp;
    int err = 0;

    if (i < 0 || i > 3) {
	i = 0;
    }

    if (uguide) {
	fname = guide_files[i];
    } else {
	fname = ref_files[i];
    }

    full_doc_path(fullpath, fname);

    /* see if file exists locally */
    fp = fopen(fullpath, "r");
    if (fp != NULL) {
	fclose(fp);
	return 0;
    }

    /* if not, grab from server */
    err = retrieve_manfile(fname, fullpath, errtext);
    if (err) {
	if (*errtext) {
	    errbox(errtext);
	} else {
	    errbox("Failed to download file");
	}
    }

    return err;
}

void display_pdf_help (gpointer p, guint uguide, GtkWidget *w)
{
    char fullpath[FILENAME_MAX];
    int pref, err = 0;

#ifdef OSX_PKG
    osx_help(uguide);
    return;
#endif  

    pref = get_manpref();

    err = maybe_grab_pdf(uguide, pref, fullpath);
    if (err) {
	return;
    }

#ifdef G_OS_WIN32
    if ((int) ShellExecute(NULL, "open", fullpath, NULL, NULL, SW_SHOW) <= 32) {
	DWORD dw = GetLastError();
	win_show_error(dw);
    }
#else
    gretl_fork(viewpdf, fullpath);
#endif
}
