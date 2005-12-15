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

#ifndef G_OS_WIN32
# include <sys/stat.h>
# include <sys/types.h>
# include <fcntl.h>
# include <unistd.h>
# include <dirent.h>
#endif

#define HDEBUG 0

static int translated_helpfile = -1;
static char *en_gui_helpfile;
static char *en_cli_helpfile;

typedef struct help_head_t help_head;

struct help_head_t {
    char *name;        /* name of heading, e.g. "Testing" */
    int *topics;       /* array of topics under heading by command number */
    char **topicnames; /* array of descriptive topic names (GUI help only) */
    int *pos;          /* array of byte offsets into file for topics */
    int ntopics;       /* number of topics under this heading */
};

static help_head **cli_heads, **gui_heads;
static help_head **en_cli_heads, **en_gui_heads;

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

/* extract the word following "# " at the beginning of a new line
   in a gretl help file (up to a max length of len) 
*/

static int get_following_word (char *s, int len, int *b, FILE *fp)
{
    int c, i = 0;

    memset(s, 0, len);

    /* skip one space */
    if ((c = getc(fp)) == EOF) {
	return 1;
    }

    while ((c = getc(fp)) != EOF && i < len) {
	if (c == ' ' || c == '\n') {
	    i++;
	    break;
	} else {
	    s[i++] = c;
	}
    }

    *b += i + 1;

    return (*s == 0);
}

static int help_topic_number (const char *word, int gui)
{
    int hnum = gretl_command_number(word);

    if (hnum == 0 && gui) {
	hnum = extra_command_number(word);
    }

    return hnum;
}

static int 
get_help_word (char *s, help_head **heads, int nh, int gui, int *pi, int *pj)
{
    int hnum;
    int i, j;

    hnum = help_topic_number(s, gui);
    if (hnum == 0) {
	*pi = -1;
	*pj = -1;
	return 0;
    }

    for (i=0; i<nh; i++) {
	for (j=0; j<heads[i]->ntopics; j++) {
	    if (hnum == heads[i]->topics[j]) {
		*pi = i;
		*pj = j;
		return 0;
	    }
	}
    }

    return 1;
}

static int get_byte_positions (help_head **heads, int nh, int gui, FILE *fp)
{
    char word[16];
    int m = 0, n = 0;
    int c, cbak = 0;
    int i, j, b = 0;
    int err = 0;

    rewind(fp);

    for (i=0; i<nh; i++) {
	n += heads[i]->ntopics;
    }

    while ((c = getc(fp)) != EOF) {
	if (c == '#' && cbak == '\n') {
	    int pos = b;

	    err = get_following_word(word, 16, &b, fp);

	    if (!err) {
		/* look up 'test' and set pos */
		err = get_help_word(word, heads, nh, gui, &i, &j);
	    }

	    if (!err && i >= 0 && j >= 0) {
		heads[i]->pos[j] = pos;
		m++;
	    }

#if HDEBUG
	    if (i >= 0 && j >= 0) {
		fprintf(stderr, "%s: setting heads[%d]->pos[%d] = %d\n", word, i, j, pos);
	    } else {
		fprintf(stderr, "%s: pos = %d\n", word, pos);
	    }
#endif

	    if (err || m == n) {
		break;
	    }
	}
	b++;
	cbak = c;
    }

#if HDEBUG
    fprintf(stderr, "get_byte_positions: n = %d, m = %d, returning %d\n", 
	    n, m, err);
#endif

    return err;
}

static void free_help_head (help_head *head)
{
    int i;

    free(head->name);
    free(head->topics);
    free(head->pos);
    
    if (head->topicnames != NULL) {
	for (i=0; i<head->ntopics; i++) {
	    free(head->topicnames[i]);
	}
	free(head->topicnames);
    }

    free(head);
}

static int head_allocate_topicnames (help_head *head)
{
    int i;

    head->topicnames = malloc(head->ntopics * sizeof *head->topicnames);
    
    if (head->topicnames == NULL) {
	return 1;
    }

    for (i=0; i<head->ntopics; i++) {
	head->topicnames[i] = NULL;
    }

    return 0;
}

static help_head *help_head_new (const char *name, int nt, int gui)
{
    help_head *head = malloc(sizeof *head);

    if (head != NULL) {
	head->ntopics = nt;
	head->name = NULL;
	head->topics = NULL;
	head->pos = NULL;
	head->topicnames = NULL;
    }

    head->name = gretl_strdup(name);
    head->topics = malloc(nt * sizeof *head->topics);
    head->pos = malloc(nt * sizeof *head->pos);

    if (head->name == NULL || head->topics == NULL ||
	head->pos == NULL) {
	free_help_head(head);
	head = NULL;
    } else if (gui && head_allocate_topicnames(head)) {
	free_help_head(head);
	head = NULL;
    }

    return head;
}

static help_head **allocate_heads (int nh)
{
    help_head **heads;
    int i;

    heads = malloc(nh * sizeof *heads);

    if (heads != NULL) {
	for (i=0; i<nh; i++) {
	    heads[i] = NULL;
	}
    }

    return heads;
}

static void free_help_heads (help_head **heads, int nh)
{
    int i;

    if (heads == NULL) return;

    for (i=0; i<nh; i++) {
	if (heads[i] != NULL) {
	    free_help_head(heads[i]);
	}
    }

    free(heads);
}

static int add_topic_to_head (help_head *head, int j, const char *word,
			      char *label, int gui)
{
    int hnum, err = 0;

    hnum = help_topic_number(word, gui);

    if (hnum > 0) {
	head->topics[j] = hnum;
	if (label != NULL && head->topicnames != NULL) {
	    head->topicnames[j] = label;
	}
    } else {
	err = 1;
    }

#if HDEBUG
    fprintf(stderr, "add_topic_to_head: topic %d: word = %s: hnum = %d\n", 
	    j+1, word, hnum);
#endif

    return err;
}

static int 
get_helpfile_structure (help_head ***pheads, int gui, const char *fname)
{
    help_head **heads = NULL;
    FILE *fp;
    char line[128];
    char tmp[32];
    int i, j, nh2;
    int nh = 0;
    int err = 0;

    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	fprintf(stderr, I_("help file %s is not accessible\n"), fname);
	return 1;
    }

    if (fgets(line, sizeof line, fp) == NULL) {
	err = 1;
	goto bailout;
    }
    
    if (!sscanf(line, "headings %d", &nh) || nh <= 0) {
	err = 1;
	goto bailout;
    }

    heads = allocate_heads(nh + 1); /* NULL sentinel at end */
    if (heads == NULL) {
	err = 1;
	goto bailout;
    }

#if HDEBUG
    fprintf(stderr, "Found got %d topic headings\n", nh);
#endif

    nh2 = 0;

    /* loop across the headings */

    for (i=0; i<nh; i++) {
	int nt;

	if (fgets(line, sizeof line, fp) == NULL) {
	    err = 1;
	    goto bailout;
	}

	if (string_is_blank(line)) {
	    break;
	}

	/* heading name plus number of following topics */
	if (sscanf(line, "%31s %d", tmp, &nt) != 2) {
	    err = 1;
	    goto bailout;
	}	    

	if (nt == 0 || !strcmp(tmp, "Obsolete")) {
	    continue;
	}

#if HDEBUG
	fprintf(stderr, "Heading %d (%s): got %d topics\n", nh2+1, tmp, nt);
#endif

	/* heading with at least one topic */
	heads[nh2] = help_head_new(tmp, nt, gui);
	if (heads[nh2] == NULL) {
	    err = 1;
	}

	/* get topics under heading */
	for (j=0; j<nt && !err; j++) {
	    char *label = NULL;

	    if (fgets(line, sizeof line, fp) == NULL) {
		err = 1;
	    }

	    if (!err) {
		err = sscanf(line, "%15s", tmp) != 1;
	    }

	    if (gui && !err) {
		label = quoted_help_string(line);
		if (label == NULL) {
		    err = 1;
		}
	    }

	    if (!err) {
		err = add_topic_to_head(heads[nh2], j, tmp, label, gui);
		if (err && label != NULL) {
		    free(label);
		}
	    }
	}

	if (err) {
	    goto bailout;
	}

	nh2++;
    }

    /* shrink array? */
    if (nh2 < nh) {
	heads = realloc(heads, (nh2 + 1) * sizeof *heads);
	nh = nh2;
    }

#if HDEBUG
    fprintf(stderr, "\nNumber of non-empty headings = %d\n\n", nh);
#endif

    /* find bytes offsets for all topics */
    get_byte_positions(heads, nh, gui, fp);

 bailout:

    fclose(fp);

    if (err) {
	fprintf(stderr, "*** get_helpfile_structure: err = %d\n", err);
	free_help_heads(heads, nh);
    } else {
	*pheads = heads;
    }

    return err;
}

static int real_helpfile_init (int gui, int en)
{
    help_head **heads = NULL;
    const char *helpfile;
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

    err = get_helpfile_structure(&heads, gui, helpfile);

    if (err) {
	sprintf(errtext, _("help file %s is not accessible\n"), helpfile);
	errbox(errtext);
    } else {
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
    help_head **heads;
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
    help_head **heads;
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
    help_head **hds;

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

#ifdef G_OS_WIN32

static char *full_doc_path (char *path, const char *fname, int *err)
{
    strcpy(path, paths.gretldir);
    strcat(path, "doc");
    strcat(path, SLASHSTR);
    strcat(path, fname);

    return path;
}

#else

static int helpfile_exists (char *path, const char *fname, int i)
{
    FILE *fp = NULL;
    int ret = 0;

    if (i == 0) {
	/* standard location */
	sprintf(path, "%sdoc/%s", paths.gretldir, fname);
    } else if (i == 1) {
	sprintf(path, "%s../doc/%s", paths.gretldir, fname);
	/* "system" location */
    } else {
	/* "user" location */
	sprintf(path, "%sdoc/%s", paths.userdir, fname);
    }

    fp = fopen(path, "r");

    if (fp != NULL) {
	fclose(fp);
	ret = 1;
    }

    return ret;
}

static int helpfile_is_writable (char *path, const char *fname, int i)
{
    FILE *fp;
    int err = 0;
    int ret = 0;

    if (i == 0) {
	/* standard location */
	sprintf(path, "%sdoc/%s", paths.gretldir, fname);
    } else if (i == 1) {
	sprintf(path, "%s../doc/%s", paths.gretldir, fname);
	/* "system" location */
    } else {
	/* "user" location */
	DIR *dir;

	sprintf(path, "%sdoc", paths.userdir);
	dir = opendir(path);

	if (dir != NULL) {
	    closedir(dir);
	} else if (mkdir(path, 0755)) {
	    err = 1;
	}
	if (!err) {
	    sprintf(path, "%sdoc/%s", paths.userdir, fname);
	}
    }

    if (!err) {
	fp = fopen(path, "w");
	if (fp != NULL) {
	    fclose(fp);
	    remove(path);
	    ret = 1;
	}
    }

    return ret;
}

static char *full_doc_path (char *path, const char *fname)
{
    int i;

    for (i=0; i<3; i++) {
	if (helpfile_exists(path, fname, i)) {
	    return path;
	} else if (helpfile_is_writable(path, fname, i)) {
	    return path;
	}
    }

    return NULL;
}

#endif

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

    if (full_doc_path(fullpath, fname) == NULL) {
	errbox("Failed to download file");
	return 1;
    }

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
