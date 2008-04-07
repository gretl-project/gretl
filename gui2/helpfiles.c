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

/* helpfiles.c for gretl */

#include "gretl.h"
#include "textbuf.h"
#include "gretl_www.h"
#include "treeutils.h"
#include "dlgutils.h"

#ifdef G_OS_WIN32
# include "gretlwin32.h"
#else
# include <sys/stat.h>
# include <sys/types.h>
# include <fcntl.h>
# include <unistd.h>
# include <dirent.h>
#endif

#define HDEBUG 0

enum {
    HELP_CLI =   1 << 0,
    HELP_EN  =   1 << 1,
    HELP_FUNCS = 1 << 2
};

static int translated_helpfile = -1;
static char *en_gui_helpfile;
static char *en_cli_helpfile;

typedef struct help_head_t help_head;

struct help_head_t {
    char *name;        /* name of heading, e.g. "Testing" */
    int *topics;       /* array of topics under heading by command number */
    char **topicnames; /* array of descriptive topic names (GUI help only) */
    int *pos;          /* array of byte offsets into file for topics */
    int ntalloc;       /* number of topics initally allocated */
    int ntopics;       /* actual number of topics under this heading */
};

static help_head **cli_heads, **gui_heads;
static help_head **en_cli_heads, **en_gui_heads;

static windata_t *make_helpwin (int flags);
static void real_do_help (int hcode, int pos, int flags);
static void en_text_cmdref (gpointer p, guint u, GtkWidget *w);

/* searching stuff */
static void find_in_text (GtkWidget *widget, gpointer data);
static void find_in_listbox (GtkWidget *widget, gpointer data);
static void find_string_dialog (void (*findfunc)(), windata_t *vwin);

static GtkWidget *find_window = NULL;
static GtkWidget *find_entry;
static char *needle;

GtkItemFactoryEntry help_menu_items[] = {
    { N_("/_Topics"), NULL, NULL, 0, "<Branch>", GNULL },    
    { N_("/_Find"), NULL, NULL, 0, "<Branch>", GNULL },   
    { N_("/Find/_Find in window"), NULL, menu_find, 0, "<StockItem>", GTK_STOCK_FIND },
    { NULL, NULL, NULL, 0, NULL, GNULL }
};

struct gui_help_item {
    int code;
    char *string;
};

static struct gui_help_item gui_help_items[] = {
    { 0,              "nothing" },
    { GR_PLOT,        "graphing" },
    { GR_XY,          "graphing" },
    { GR_DUMMY,       "factorized" },
    { GR_XYZ,         "controlled" },
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
    { GENR_RANDOM,    "genrand" },
    { SEED_RANDOM,    "genseed" },
    { KERNEL_DENSITY, "density" },
    { HCCME,          "hccme" },
    { IRF_BOOT,       "irfboot" },
    { HTEST,          "gui-htest" },
    { HTESTNP,        "gui-htest-np" },
    { MODEL_RESTR,    "restrict-model" },
    { SYS_RESTR,      "restrict-system" },
    { VECM_RESTR,     "restrict-vecm" },
    { LAGS_DIALOG,    "lags-dialog" },
    { COPY_FORMATS,   "copy-formats" },
    { MINIBUF,        "minibuffer" },
    { VLAGSEL,        "VAR-lagselect" },
    { VAROMIT,        "VAR-omit" },
    { PANEL_MODE,     "panel-mode" },
    { PANEL_WLS,      "panel-wls" },
    { PANEL_B,        "panel-between" },
    { BOOTSTRAP,      "bootstrap" },
    { TRANSPOS,       "transpos" },
    { DATASORT,       "datasort" },
    { WORKDIR,        "working-dir" },
    { DFGLS,          "dfgls" },
    { -1,             NULL },
};

enum {
    COMPAT_CORC = GUI_CMD_MAX + 1,
    COMPAT_FCASTERR,
    COMPAT_HILU,
    COMPAT_PWE
};

static struct gui_help_item compat_help_items[] = {
    { GUI_CMD_MAX,     "nothing" },
    { COMPAT_CORC,     "corc" },
    { COMPAT_FCASTERR, "fcasterr" },
    { COMPAT_HILU,     "hilu" },
    { COMPAT_PWE,      "pwe" },
    { -1,              NULL }
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

static int compat_command_number (const char *s)
{
    int i;

    for (i=1; compat_help_items[i].code > 0; i++) {
	if (!strcmp(s, compat_help_items[i].string)) {
	    return compat_help_items[i].code;
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

    for (i=1; compat_help_items[i].code > 0; i++) {
	if (cmd == compat_help_items[i].code) {
	    return compat_help_items[i].string;
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

    if (hnum <= 0) {
	hnum = compat_command_number(word);
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
    int i, j = 0, b = 0;
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
	for (i=0; i<head->ntalloc; i++) {
	    free(head->topicnames[i]);
	}
	free(head->topicnames);
    }

    free(head);
}

static int head_allocate_topicnames (help_head *head)
{
    int i, n = head->ntopics;

    head->topicnames = malloc(n * sizeof *head->topicnames);
    
    if (head->topicnames == NULL) {
	return 1;
    }

    for (i=0; i<n; i++) {
	head->topicnames[i] = NULL;
    }
    
    head->ntalloc = n;

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

static help_head **allocate_heads (int n)
{
    help_head **heads;
    int i;

    heads = malloc(n * sizeof *heads);

    if (heads != NULL) {
	for (i=0; i<n; i++) {
	    heads[i] = NULL;
	}
    }

    return heads;
}

static void free_help_heads (help_head **heads, int nh)
{
    int i;

    if (heads == NULL) {
	return;
    }

    for (i=0; i<nh; i++) {
	if (heads[i] != NULL) {
	    free_help_head(heads[i]);
	}
    }

    free(heads);
}

static void add_topic_to_head (help_head *head, int *i, const char *word,
			       char *label, int gui)
{
    int hnum = help_topic_number(word, gui);

    if (hnum > 0) {
	head->topics[*i] = hnum;
	if (label != NULL && head->topicnames != NULL) {
	    head->topicnames[*i] = label;
	}
	*i += 1;
#if HDEBUG
	fprintf(stderr, "add_topic_to_head: topic %d: word = %s: hnum = %d\n", 
		*i, word, hnum);
#endif
    } else {
	fprintf(stderr, "helpfile error: word='%s' not recognized\n", word);
	head->ntopics -= 1;
	if (label != NULL) {
	    free(label);
	}
    }
}

static int 
get_helpfile_structure (help_head ***pheads, int gui, const char *fname)
{
    help_head **heads = NULL;
    FILE *fp;
    char line[128];
    char tmp[32];
    int i, j, k, nh2;
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
	fprintf(stderr, "Couldn't reading number of headings in helpfile\n");
	err = 1;
	goto bailout;
    }

    heads = allocate_heads(nh + 1); /* NULL sentinel at end */
    if (heads == NULL) {
	err = 1;
	goto bailout;
    }

#if HDEBUG
    fprintf(stderr, "Found %d topic headings\n", nh);
#endif

    nh2 = 0;

    /* loop across the headings */

    for (i=0; i<nh; i++) {
	int nt;

	if (fgets(line, sizeof line, fp) == NULL) {
	    fprintf(stderr, "Couldn't read line from helpfile: i=%d\n", i);
	    err = 1;
	    goto bailout;
	}

	if (string_is_blank(line)) {
	    break;
	}

	/* heading name plus number of following topics */
	if (sscanf(line, "%31s %d", tmp, &nt) != 2) {
	    fprintf(stderr, "Couldn't read heading and number of topics: i=%d\n", i);
	    err = 1;
	    goto bailout;
	}	    

	if (nt == 0 || !strcmp(tmp, "Obsolete")) {
	    continue;
	}

#if HDEBUG
	fprintf(stderr, "Heading %d (\"%s\"): got %d topics\n", nh2+1, tmp, nt);
#endif

	/* heading with at least one topic */
	heads[nh2] = help_head_new(tmp, nt, gui);
	if (heads[nh2] == NULL) {
	    fprintf(stderr, "help_head_new() failed\n");
	    err = 1;
	}

	k = 0;

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
		add_topic_to_head(heads[nh2], &k, tmp, label, gui);
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
	fprintf(stderr, "*** get_helpfile_structure: err = %d for\n"
		"filename = '%s'\n", err, fname);
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
	warnbox(_("help file %s is not up to date\n"), helpfile);
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

/* apparatus for genr functions help */

typedef struct func_finder_ func_finder;

struct func_finder_ {
    int pos;
    char word[12];
};

static func_finder *ffinder;
static int n_help_funcs;

/* read the functions helpfile and build a mapping from
   function name to offset in file */

static int make_func_help_mapping (void)
{
    int err = 0;

    if (ffinder == NULL) {
	gchar *fname = NULL;
	gchar *s, *buf = NULL;
	char word[12];
	int i = 0, pos = 0;
	int err;

	fname = g_strdup_printf("%s%s", paths.gretldir, _("genrgui.hlp"));

	err = gretl_file_get_contents(fname, &buf);
	if (err) {
	    g_free(fname);
	    return err;
	}

	s = buf;
	while (*s) {
	    if (*s == '\n' && *(s+1) == '#' && *(s+2) != '#') {
		n_help_funcs++;
	    }
	    s++;
	}

	ffinder = mymalloc(n_help_funcs * sizeof *ffinder);

	if (ffinder != NULL) {
	    s = buf;
	    while (*s) {
		if (*s == '\n' && *(s+1) == '#' && 
		    *(s+2) != '#' && *(s+2) != '\0') {
		    if (sscanf(s+2, "%10s", word)) {
			ffinder[i].pos = pos + 1;
			strcpy(ffinder[i].word, word);
			i++;
		    }
		}
		s++;
		pos++;
	    }
	} else {
	    err = 1;
	    n_help_funcs = 0;
	}

	g_free(buf);
	g_free(fname);
    }

    return err;
}

int function_help_index_from_word (const char *s)
{
    int i, idx = 0;

    if (n_help_funcs == 0) {
	make_func_help_mapping();
    }

    for (i=0; i<n_help_funcs; i++) {
	if (!strcmp(ffinder[i].word, s)) {
	    idx = i+1;
	    break;
	}
    }

    return idx;
}

static int function_help_pos_from_word (const char *s)
{
    int i, pos = 0;

    if (n_help_funcs == 0) {
	make_func_help_mapping();
    }

    for (i=0; i<n_help_funcs; i++) {
	if (!strcmp(ffinder[i].word, s)) {
	    pos = ffinder[i].pos;
	}
    }

    return pos;
}

static int help_pos_from_function_index (int fnum)
{
    int pos = 0;

    if (n_help_funcs == 0) {
	make_func_help_mapping();
    }

    if (fnum > 0 && fnum <= n_help_funcs) {
	pos = ffinder[fnum-1].pos;
    } 

    return pos;
}

/* end apparatus for genr functions help */

static void do_gui_help (gpointer p, guint pos, GtkWidget *w) 
{
    int hcode = GPOINTER_TO_INT(p);

    real_do_help(hcode, pos, 0);
}

static void do_en_gui_help (gpointer p, guint pos, GtkWidget *w) 
{
    int hcode = GPOINTER_TO_INT(p);

    real_do_help(hcode, pos, HELP_EN);
}

static void do_cli_help (gpointer p, guint pos, GtkWidget *w) 
{
    int hcode = GPOINTER_TO_INT(p);

    real_do_help(hcode, pos, HELP_CLI);
}

static void do_en_cli_help (gpointer p, guint pos, GtkWidget *w) 
{
    int hcode = GPOINTER_TO_INT(p);

    real_do_help(hcode, pos, HELP_CLI | HELP_EN);
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

static char *get_compat_help_string (int pos)
{
    int i, j;

    for (i=0; cli_heads[i] != NULL; i++) { 
	for (j=0; j<cli_heads[i]->ntopics; j++) {
	    if (pos == cli_heads[i]->pos[j]) {
		return help_string_from_cmd(cli_heads[i]->topics[j]);
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
    int flags = HELP_EN;

    if (cli) {
	flags |= HELP_CLI;
	pos = cli_pos_from_cmd(hc, 1);
	if (pos < 0) {
	    pos = 0;
	}
    } else {
	pos = gui_pos_from_cmd(hc, 1);
    }

    real_do_help(hc, pos, flags);
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

static void add_help_topics (windata_t *hwin, int flags)
{
    GtkItemFactoryEntry hitem;
    const gchar *mpath = N_("/_Topics");
    help_head **hds = NULL;
    int en = (flags & HELP_EN);
    int cli = (flags & HELP_CLI);
    int funcs = (flags & HELP_FUNCS);
    int i, j;

    if (en) {
	hds = (cli)? en_cli_heads : (funcs)? NULL : en_gui_heads;
    } else {
	hds = (cli)? cli_heads : (funcs)? NULL : gui_heads;
    }

    hitem.accelerator = NULL;

    if (cli || funcs) {
	/* Add general index as (first) "topic" */
	hitem.callback_action = 0; 
	hitem.item_type = NULL;
	hitem.path = g_strdup_printf("%s/%s", mpath, (en)? "Index" : _("Index"));
	if (funcs) {
	    hitem.callback = genr_funcs_ref;
	} else {
	    hitem.callback = (en)? en_text_cmdref : plain_text_cmdref;
	}
	gtk_item_factory_create_item(hwin->ifac, &hitem, 
				     GINT_TO_POINTER(0), 
				     1);
	g_free(hitem.path);
    }

    /* Are there any actual topics to add? */
    if (hds == NULL) {
	return;
    }

    /* put the topics under the menu heading */
    for (i=0; hds[i] != NULL; i++) {
	const char *headname = (en)? hds[i]->name : _(hds[i]->name);

	if (hds[i]->ntopics == 0) {
	    continue;
	}

	hitem.callback_action = 0; 
	hitem.item_type = "<Branch>";
	hitem.path = g_strdup_printf("%s/%s", mpath, headname);
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
				    mpath, headname, 
				    hds[i]->topicnames[j]);
#if HDEBUG
		fprintf(stderr, "(1) Built help topic path from\n"
			" '%s', '%s' and '%s'\n", mpath, headname,
			hds[i]->topicnames[j]);
#endif
	    } else {
		if (tnum < NC) {
		    /* a regular gretl command */
		    hitem.path = 
			g_strdup_printf("%s/%s/%s", 
					mpath, headname, 
					gretl_command_word(tnum));
#if HDEBUG
		    fprintf(stderr, "(2) Built help topic path from\n"
			    " '%s', '%s' and '%s'\n", mpath, headname,
			    gretl_command_word(tnum));
#endif
		} else if (!cli) {
		    /* a gui special item? */
		    char *gstr = get_gui_help_string(hds[i]->pos[j]);

		    if (gstr != NULL) {
			hitem.path = 
			    g_strdup_printf("%s/%s/%s", 
					    mpath, headname, gstr);
#if HDEBUG
			fprintf(stderr, "(3) Built help topic path from\n"
				" '%s', '%s' and '%s'\n", mpath, headname,
				gstr);
#endif
		    } else {
			topic_ok = 0;
		    }
		} else {
		    char *cstr = get_compat_help_string(hds[i]->pos[j]);

		    if (cstr != NULL) {
			hitem.path = 
			    g_strdup_printf("%s/%s/%s", 
					    mpath, headname, cstr);
		    } else {		    
			topic_ok = 0;
		    }
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

static char *funcs_helpfile (void)
{
    static char fname[MAXLEN];

    if (*fname == '\0') {
	sprintf(fname, "%s%s", paths.gretldir, _("genrgui.hlp"));
    }

    return fname;
}

static windata_t *make_helpwin (int flags) 
{
    windata_t *vwin = NULL;
    char *helpfile = NULL;
    int en = (flags & HELP_EN);
    int cli = (flags & HELP_CLI);
    int funcs = (flags & HELP_FUNCS);
    int role;

    if (funcs) {
	helpfile = funcs_helpfile();
	role = FUNCS_HELP;
    } else if (cli) {
	helpfile = (en)? en_cli_helpfile : paths.cmd_helpfile;
	role = (en)? CLI_HELP_EN : CLI_HELP;
    } else {
	helpfile = (en)? en_gui_helpfile : paths.helpfile;
	role = (en)? GUI_HELP_EN : GUI_HELP;
    }

    if (helpfile != NULL) {
	vwin = view_help_file(helpfile, role, help_menu_items);
    }

    if (vwin != NULL) {
	add_help_topics(vwin, flags);
	if (!funcs && translated_helpfile && !en) {
	    add_en_help_item(vwin, cli);
	} 
    }   

    return vwin;
}

void context_help (GtkWidget *widget, gpointer data)
{
    int pos, hcode = GPOINTER_TO_INT(data);
    int flags = 0;

    pos = gui_pos_from_cmd(hcode, 0);

    if (pos < 0 && translated_helpfile) {
	flags = HELP_EN;
	pos = gui_pos_from_cmd(hcode, 1);
    }

#if HDEBUG
    fprintf(stderr, "context_help: hcode=%d, pos=%d, en=%d\n", 
	    hcode, pos, (flags & HELP_EN)? 1 : 0);
#endif

    real_do_help(hcode, pos, flags);
}

static gboolean nullify_hwin (GtkWidget *w, windata_t **phwin)
{
    *phwin = NULL;
    return FALSE;
}

static void real_do_help (int hcode, int pos, int flags)
{
    static windata_t *gui_hwin;
    static windata_t *cli_hwin;
    static windata_t *en_gui_hwin;
    static windata_t *en_cli_hwin;
    static windata_t *funcs_hwin;

    windata_t *hwin;
    int cli = (flags & HELP_CLI);
    int en = (flags & HELP_EN);
    int funcs = (flags & HELP_FUNCS);

    if (pos < 0) {
	dummy_call();
	return;
    }

#if HDEBUG
    fprintf(stderr, "real_do_help: hcode=%d, pos=%d, cli=%d, en=%d, funcs=%d\n",
	    hcode, pos, cli, en, funcs);
    fprintf(stderr, "gui_hwin = %p\n", (void *) gui_hwin);
    fprintf(stderr, "cli_hwin = %p\n", (void *) cli_hwin);
    fprintf(stderr, "en_gui_hwin = %p\n", (void *) en_gui_hwin);
    fprintf(stderr, "en_cli_hwin = %p\n", (void *) en_cli_hwin);
    fprintf(stderr, "funcs_hwin = %p\n", (void *) funcs_hwin);
#endif

    if (en) {
	hwin = (cli)? en_cli_hwin : en_gui_hwin;
    } else {
	hwin = (funcs)? funcs_hwin : (cli)? cli_hwin : gui_hwin;
    }

    if (hwin == NULL) {
	hwin = make_helpwin(flags);
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
		if (funcs) {
		    funcs_hwin = hwin;
		    phwin = &funcs_hwin;
		} else if (cli) {
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
	gtk_window_present(GTK_WINDOW(gtk_widget_get_toplevel(hwin->w)));
    }

#if HDEBUG
    fprintf(stderr, "real_do_help: doing set_help_topic_buffer:\n"
	    " hwin=%p, hcode=%d, pos=%d, en=%d\n",
	    (void *) hwin, hcode, pos, en);
#endif

    if (hwin != NULL) {
	set_help_topic_buffer(hwin, hcode, pos, en);
    }
}

/* called from main menu in gretl.c; also used as callback from
   command ref index in textbuf.c
*/

void plain_text_cmdref (gpointer p, guint cmdnum, GtkWidget *w)
{
    /* pos = 0 gives index of commands */
    int pos = 0;
    int flags = HELP_CLI;
    int en = 0;

    if (w == NULL && p != NULL) {
	en = GPOINTER_TO_INT(p);
    }

#if HDEBUG
    fprintf(stderr, "plain_text_cmdref: p=%p, cmdnum=%d, w=%p\n", 
	    (void *) p, cmdnum, (void *) w);
#endif

    if (cmdnum > 0) {
	pos = cli_pos_from_cmd(cmdnum, en);
	if (pos < 0 && !en && translated_helpfile == 1) {
	    /* no translated entry: fall back on English */
	    en = 1;
	    pos = cli_pos_from_cmd(cmdnum, 1);
	}
    }

    if (en) {
	flags |= HELP_EN;
    }

    real_do_help(cmdnum, pos, flags);
} 

/* called from main menu in gretl.c; also used as callback from
   command ref index in textbuf.c
*/

void genr_funcs_ref (gpointer p, guint fnum, GtkWidget *w)
{
    int pos = help_pos_from_function_index(fnum);

    real_do_help(fnum, pos, HELP_FUNCS);    
}

static void en_text_cmdref (gpointer p, guint u, GtkWidget *w)
{
    real_do_help(0, 0, HELP_CLI | HELP_EN);
} 

static int help_pos_from_string (const char *s, int *en, 
				 int *hcode, int *flags)
{
    char word[12];
    int pos;

    *word = '\0';
    strncat(word, s, 11);

    *hcode = gretl_command_number(word);
    pos = cli_pos_from_cmd(*hcode, *en);

    if (pos < 0 && translated_helpfile) {
	pos = cli_pos_from_cmd(*hcode, 1);
	if (pos >= 0) {
	    *en = 1;
	}
    }

    if (pos < 0) {
	/* try function instead of command */
	pos = function_help_pos_from_word(word);
	if (pos > 0) {
	    *flags |= HELP_FUNCS;
	}
    }

    return pos;
}

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
	int flags = HELP_CLI;
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

	    /* dollar accessors! */
	    if (text != NULL) {
		GtkTextIter dstart = w_start;
		gchar *dtest = NULL;

		if (gtk_text_iter_backward_char(&dstart)) {
		    dtest = gtk_text_buffer_get_text(buf, &dstart, 
						     &w_start, FALSE);
		    if (*dtest == '$') {
			gchar *s = g_strdup_printf("$%s", text);
			
			g_free(text);
			text = s;
		    }
		    g_free(dtest);
		}
	    }

	    /* "coint2" command! */
	    if (text != NULL && !strcmp(text, "coint") && 
		gtk_text_iter_forward_char(&w_end)) {
		gchar *s = gtk_text_buffer_get_text(buf, &w_start, &w_end, FALSE);

		if (s != NULL) {
		    if (!strcmp(s, "coint2")) {
			g_free(text);
			text = s;
		    } else {
			g_free(s);
		    }
		}
	    }
	} 

	if (text != NULL && *text != '\0') {
	    pos = help_pos_from_string(text, &en, &hcode, &flags);
	} 

	g_free(text);
	unset_window_help_active(vwin);
	text_set_cursor(vwin->w, 0);

	if (pos <= 0) {
	    warnbox(_("Sorry, help not found"));
	} else {
	    if (en) {
		flags |= HELP_EN;
	    }
	    real_do_help(hcode, pos, flags);
	}
    }

    return FALSE;
}

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

/* case-insensitive search in text buffer */

static gboolean real_find_in_text (GtkTextView *view, const gchar *s, 
				   gboolean from_cursor)
{
    GtkTextBuffer *buf;
    GtkTextIter iter, start, end;
    GtkTextMark *vis;
    int found = 0;
    int wrapped = 0;
    int n = strlen(s);
    gchar *got;

    buf = gtk_text_view_get_buffer(view);

 text_search_wrap:
	
    if (from_cursor) {
	GtkTextIter sel_bound;
		
	gtk_text_buffer_get_iter_at_mark(buf, &iter,
					 gtk_text_buffer_get_mark(buf,
								  "insert"));
	gtk_text_buffer_get_iter_at_mark(buf, &sel_bound,
					 gtk_text_buffer_get_mark(buf,
								  "selection_bound"));
	gtk_text_iter_order(&sel_bound, &iter);		
    } else {		
	gtk_text_buffer_get_iter_at_offset(buf, &iter, 0);
    }

    start = end = iter;

    if (!gtk_text_iter_forward_chars(&end, n)) {
	return 0;
    }

    while (!found) {
	got = gtk_text_buffer_get_text(buf, &start, &end, FALSE);
	if (g_ascii_strcasecmp(got, s) == 0) {
	    found = 1;
	}
	g_free(got);
	if (found || !gtk_text_iter_forward_char(&start) ||
	    !gtk_text_iter_forward_char(&end)) {
	    break;
	}
    }

    if (found) {
	gtk_text_buffer_place_cursor(buf, &start);
	gtk_text_buffer_move_mark_by_name(buf, "selection_bound", &end);
	vis = gtk_text_buffer_create_mark(buf, "vis", &end, FALSE);
	gtk_text_view_scroll_to_mark(view, vis, 0.0, TRUE, 0.1, 0.0);
    } else if (from_cursor && !wrapped) {
	/* try wrapping */
	from_cursor = FALSE;
	wrapped = 1;
	goto text_search_wrap;
    }

    if (found && wrapped) {
	infobox(_("Search wrapped"));
    }

    return found;
}

static void find_in_text (GtkWidget *widget, gpointer data)
{
    windata_t *vwin = g_object_get_data(G_OBJECT(data), "windat");
    gboolean found;

    needle = gtk_editable_get_chars(GTK_EDITABLE(find_entry), 0, -1);
    if (needle == NULL || *needle == '\0') {
	return;
    }

    found = real_find_in_text(GTK_TEXT_VIEW(vwin->w), needle, TRUE);

    if (!found) {
	infobox(_("String was not found."));
    }
}

static void 
get_tree_model_haystack (GtkTreeModel *mod, GtkTreeIter *iter, int col,
			 char *haystack)
{
    gchar *tmp;

    gtk_tree_model_get(mod, iter, col, &tmp, -1);
    if (tmp != NULL) {
	strcpy(haystack, tmp);
	lower(haystack);
	g_free(tmp);
    } else {
	*haystack = '\0';
    }
}

static void find_in_listbox (GtkWidget *w, gpointer p)
{
    windata_t *win = g_object_get_data(G_OBJECT(p), "windat");
    int minvar, wrapped = 0;
    char haystack[MAXLEN];
    char pstr[16];
    GtkTreeModel *model;
    GtkTreeIter iter;
    int found = -1;

    /* if searching in the main gretl window, start on line 1 */
    minvar = (win == mdata)? 1 : 0;

    if (needle != NULL) {
	g_free(needle);
	needle = NULL;
    }

    needle = gtk_editable_get_chars(GTK_EDITABLE(find_entry), 0, -1);
    lower(needle);

    model = gtk_tree_view_get_model(GTK_TREE_VIEW(win->listbox));

    /* if possible, try searching downward from the current line 
       plus one; failing start, start from the top */
    sprintf(pstr, "%d", win->active_var);
    if (!gtk_tree_model_get_iter_from_string(model, &iter, pstr) ||
	!gtk_tree_model_iter_next(model, &iter)) {
	gtk_tree_model_get_iter_first(model, &iter);
    }

 search_wrap:

    while (found < 0) {
	/* try looking in column 1 first */
	get_tree_model_haystack(model, &iter, 1, haystack);

	found = look_for_string(haystack, needle, 0);

	if (found < 0) {
	    if (win == mdata) {
		/* then column 2 */
		get_tree_model_haystack(model, &iter, 2, haystack);
	    } else {
		/* then column 0 */
		get_tree_model_haystack(model, &iter, 0, haystack);
	    }
	    found = look_for_string(haystack, needle, 0);
	}

	if (found >= 0 || !gtk_tree_model_iter_next(model, &iter)) {
	    break;
	}
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
	if (wrapped) {
	    infobox(_("Search wrapped"));
	}
    } else {
	infobox(_("String was not found."));
    }
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

static void find_string_dialog (void (*findfunc)(), windata_t *vwin)
{
    GtkWidget *label;
    GtkWidget *button;
    GtkWidget *hbox;

    if (find_window != NULL) {
	g_object_set_data(G_OBJECT(find_window), "windat", vwin);
	parent_find(find_window, vwin);
	gdk_window_raise(find_window->window);
	return;
    }

    find_window = gretl_dialog_new(_("gretl: find"), NULL, 0);
    g_object_set_data(G_OBJECT(find_window), "windat", vwin);
    parent_find(find_window, vwin);

    g_signal_connect(G_OBJECT(find_window), "destroy",
		     G_CALLBACK(close_find_dialog),
		     find_window);

    hbox = gtk_hbox_new(TRUE, TRUE);
    label = gtk_label_new(_(" Find what:"));
    gtk_widget_show (label);
    find_entry = gtk_entry_new();

    if (needle != NULL) {
	gtk_entry_set_text(GTK_ENTRY(find_entry), needle);
	gtk_editable_select_region(GTK_EDITABLE(find_entry), 0, -1);
    }

    g_signal_connect(G_OBJECT(find_entry), "activate", 
		     G_CALLBACK(findfunc), find_window);

    gtk_widget_show(find_entry);

    gtk_box_pack_start(GTK_BOX(hbox), label, TRUE, TRUE, 0);
    gtk_box_pack_start(GTK_BOX(hbox), find_entry, TRUE, TRUE, 0);
    gtk_widget_show(hbox);

    gtk_box_pack_start(GTK_BOX (GTK_DIALOG(find_window)->vbox), 
		       hbox, TRUE, TRUE, 0);

    hbox = GTK_DIALOG(find_window)->action_area;

    /* cancel button */
    button = cancel_button(hbox);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(cancel_find), find_window);
    gtk_widget_show(button);

    /* find button */
    button = gtk_button_new_from_stock(GTK_STOCK_FIND);
    GTK_WIDGET_SET_FLAGS(button, GTK_CAN_DEFAULT);
    gtk_container_add(GTK_CONTAINER(hbox), button);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(findfunc), find_window);
    gtk_widget_grab_default(button);
    gtk_widget_show(button);

    gtk_widget_grab_focus(find_entry);
    gtk_widget_show(find_window);
}

void text_find_callback (GtkWidget *w, windata_t *vwin)
{
    find_string_dialog(find_in_text, vwin);
}

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

enum {
    EN_LETTER,
    EN_A4,
    ITALIAN,
    SPANISH
};

static int get_writable_path (char *path, const char *fname)
{
    static int sysdoc_writable = -1;
    static int userdoc_writable = -1;

    FILE *fp;
    int err = 0;

    if (sysdoc_writable < 0) {
	sysdoc_writable = 0;
	sprintf(path, "%sdoc", paths.gretldir);
	if (gretl_mkdir(path) == 0) {
	    sprintf(path, "%sdoc%c%s", paths.gretldir, SLASH, fname);
	    fp = gretl_fopen(path, "w");
	    if (fp != NULL) {
		sysdoc_writable = 1;
		fclose(fp);
		remove(path);
	    } 
	} 
    }

    if (!sysdoc_writable && userdoc_writable < 0) {
	userdoc_writable = 0;
	sprintf(path, "%sdoc", paths.dotdir);
	if (gretl_mkdir(path) == 0) {
	    sprintf(path, "%sdoc%c%s", paths.dotdir, SLASH, fname);
	    fp = gretl_fopen(path, "w");
	    if (fp != NULL) {
		userdoc_writable = 1;
		fclose(fp);
		remove(path);
	    } 
	} 
    }

    if (!sysdoc_writable && !userdoc_writable) {
	err = 1;
    }

    return err;
}

static int find_or_download_pdf (int uguide, int i, char *fullpath)
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

    fname = (uguide)? guide_files[i] : ref_files[i];

    /* is the file available in public dir? */
    sprintf(fullpath, "%sdoc%c%s", paths.gretldir, SLASH, fname);
    fp = gretl_fopen(fullpath, "r");
    if (fp != NULL) {
	fclose(fp);
	return 0;
    }

    /* or maybe in user dir? */
    sprintf(fullpath, "%sdoc%c%s", paths.dotdir, SLASH, fname);
    fp = gretl_fopen(fullpath, "r");
    if (fp != NULL) {
	fclose(fp);
	return 0;
    }

    /* check for download location */
    err = get_writable_path(fullpath, fname);

    /* do actual download */
    if (!err) {
	err = retrieve_manfile(fname, fullpath);
    }

    if (err) {
	const char *buf = gretl_errmsg_get();

	if (*buf) {
	    errbox(buf);
	} else {
	    errbox(_("Failed to download file"));
	}
    }

    return err;
}

void display_pdf_help (gpointer p, guint uguide, GtkWidget *w)
{
    char fname[FILENAME_MAX];
    int pref, err = 0;

    pref = get_manpref();

    err = find_or_download_pdf(uguide, pref, fname);
    if (err) {
	return;
    }

#if defined(G_OS_WIN32)
    win32_open_file(fname);
#elif defined(OSX_BUILD)
    osx_open_file(fname);
#else
    gretl_fork("viewpdf", fname);
#endif
}
