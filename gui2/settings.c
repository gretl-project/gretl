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

/* settings.c for gretl */

#include "gretl.h"
#include "webget.h"

#include <unistd.h>

#if defined(USE_GNOME) && !defined(OLD_GTK)
# define GNOME2
#endif

#ifdef GNOME2
# include <gconf/gconf-client.h>
#endif

#ifdef G_OS_WIN32
# include <windows.h>
# include "gretlwin32.h"
#else
# include <dirent.h>
# include <sys/stat.h>
# include <sys/types.h>
# include <fcntl.h>
# include <errno.h>
# ifndef OLD_GTK
#  include "gtkfontselhack.h"
# endif
#endif

#if !defined(G_OS_WIN32) && !defined(USE_GNOME)
char rcfile[MAXLEN];
#endif

extern GtkWidget *toolbar_box;

extern int want_toolbar;
extern char calculator[MAXSTR];
extern char editor[MAXSTR];
extern char Rcommand[MAXSTR];
extern char dbproxy[21];

#ifdef HAVE_TRAMO
extern char tramo[MAXSTR];
extern char tramodir[MAXSTR];
#endif

int use_proxy;

/* filelist stuff */
#define MAXRECENT 4

static void printfilelist (int filetype, FILE *fp);

static char datalist[MAXRECENT][MAXSTR], *datap[MAXRECENT];
static char sessionlist[MAXRECENT][MAXSTR], *sessionp[MAXRECENT];
static char scriptlist[MAXRECENT][MAXSTR], *scriptp[MAXRECENT];

static void make_prefs_tab (GtkWidget *notebook, int tab);
static void apply_changes (GtkWidget *widget, gpointer data);
#ifndef G_OS_WIN32
static void read_rc (void);
#endif

/* font handling */
#ifdef G_OS_WIN32
static char fixedfontname[MAXLEN] = "Courier New 10";
#else
# ifndef OLD_GTK
#  ifdef OSX_BUILD
static char fixedfontname[MAXLEN] = "Monospace 12";
#  else
static char fixedfontname[MAXLEN] = "Monospace 10";
#  endif
# else
static char fixedfontname[MAXLEN] = 
"-b&h-lucidatypewriter-medium-r-normal-sans-12-*-*-*-*-*-*-*";
# endif
#endif

#if defined(G_OS_WIN32)
static char appfontname[MAXLEN] = "tahoma 8";
#elif !defined(USE_GNOME) && !defined(OLD_GTK)
# ifdef OSX_BUILD
static char appfontname[MAXLEN] = "Luxi Sans 12";
# else
static char appfontname[MAXLEN] = "Sans 10";
# endif
#endif

#ifndef OLD_GTK
PangoFontDescription *fixed_font;
#else
GdkFont *fixed_font;
#endif

static int usecwd;
int olddat;
int useqr;
char gpcolors[32];
static char datapage[24];
static char scriptpage[24];

#ifdef G_OS_WIN32
int wimp;
#endif

#ifdef ENABLE_NLS
static int lcnumeric = 1;
#endif

#if defined(HAVE_AUDIO) && !defined(G_OS_WIN32)
char midiplayer[MAXSTR];
#endif

enum {
    ROOTSET = 1 << 0,
    USERSET = 1 << 1,
    BOOLSET = 1 << 2,
    INVISET = 1 << 3,
    FIXSET  = 1 << 4  /* setting fixed by admin (Windows network use) */
};

typedef struct {
    char *key;         /* config file variable name */
    char *description; /* How the field will show up in the options dialog */
    char *link;        /* in case of radio button pair, alternate string */
    void *var;         /* pointer to variable */
    char type;         /* ROOTSET user string
			  USERSET root string
			  BOOLSET boolean (user) 
			  INVISET "invisible" (user) string 
		       */
    int len;           /* storage size for string variable (also see Note) */
    short tab;         /* which tab (if any) does the item fall under? */
    GtkWidget *widget;
} RCVARS;

/* Note: actually "len" above is overloaded: if an rc_var is of type
   BOOLSET and not part of a radio group, then a non-zero value for
   len will link the var's toggle button with the sensitivity of the
   preceding rc_var's entry field.  For example, the "use_proxy" button
   controls the sensitivity of the "dbproxy" entry widget. */

RCVARS rc_vars[] = {
    {"gretldir", N_("Main gretl directory"), NULL, paths.gretldir, 
     ROOTSET, MAXLEN, 1, NULL},
    {"userdir", N_("User's gretl directory"), NULL, paths.userdir, 
     USERSET, MAXLEN, 1, NULL},
    {"expert", N_("Expert mode (no warnings)"), NULL, &expert, 
     BOOLSET, 0, 1, NULL},
    {"updater", N_("Tell me about gretl updates"), NULL, &updater, 
     BOOLSET, 0, 1, NULL},
    {"toolbar", N_("Show gretl toolbar"), NULL, &want_toolbar, 
     BOOLSET, 0, 1, NULL},
#ifdef ENABLE_NLS
    {"lcnumeric", N_("Use locale setting for decimal point"), NULL, &lcnumeric, 
     BOOLSET, 0, 1, NULL},
#endif
#ifdef G_OS_WIN32
    {"wimp", N_("Emulate Windows look"), NULL, &wimp, BOOLSET, 0, 1, NULL},    
#endif
    {"gnuplot", N_("Command to launch gnuplot"), NULL, paths.gnuplot, 
     ROOTSET, MAXLEN, 3, NULL},
    {"Rcommand", N_("Command to launch GNU R"), NULL, Rcommand, 
     ROOTSET, MAXSTR, 3, NULL},
    {"viewdvi", N_("Command to view DVI files"), NULL, viewdvi, 
     ROOTSET, MAXSTR, 3, NULL},
#ifndef G_OS_WIN32
    {"viewps", N_("Command to view postscript files"), NULL, viewps, 
     ROOTSET, MAXSTR, 3, NULL},
#endif
#if defined(HAVE_AUDIO) && !defined(G_OS_WIN32)
    {"midiplayer", N_("Program to play MIDI files"), NULL, midiplayer, 
     USERSET, MAXSTR, 3, NULL},
#endif
    {"calculator", N_("Calculator"), NULL, calculator, USERSET, MAXSTR, 3, NULL},
#ifdef SELECT_EDITOR
    {"editor", N_("Editor"), NULL, editor, USERSET, MAXSTR, 3, NULL},
#endif
#ifdef HAVE_X12A
    {"x12a", N_("path to x12arima"), NULL, paths.x12a, ROOTSET, MAXSTR, 3, NULL},
#endif
#ifdef HAVE_TRAMO
    {"tramo", N_("path to tramo"), NULL, tramo, ROOTSET, MAXSTR, 3, NULL},
#endif
#ifdef G_OS_WIN32
    {"x12adir", N_("X-12-ARIMA working directory"), NULL, paths.x12adir, 
     ROOTSET, MAXSTR, 3, NULL},
#endif
#ifdef G_OS_WIN32
    {"tramodir", N_("TRAMO working directory"), NULL, tramodir, 
     ROOTSET, MAXSTR, 3, NULL},
#endif
    {"binbase", N_("gretl database directory"), NULL, paths.binbase, 
     USERSET, MAXLEN, 2, NULL},
    {"ratsbase", N_("RATS data directory"), NULL, paths.ratsbase, 
     USERSET, MAXLEN, 2, NULL},
    {"dbhost_ip", N_("Database server IP"), NULL, paths.dbhost_ip, 
     USERSET, 16, 2, NULL},
    {"dbproxy", N_("HTTP proxy (ipnumber:port)"), NULL, dbproxy, 
     USERSET, 21, 2, NULL},
    {"useproxy", N_("Use HTTP proxy"), NULL, &use_proxy, BOOLSET, 1, 2, NULL},
    {"usecwd", N_("Use current working directory as default"), 
     N_("Use gretl user directory as default"), &usecwd, BOOLSET, 0, 4, NULL},
    {"olddat", N_("Use \".dat\" as default datafile suffix"), 
     N_("Use \".gdt\" as default suffix"), &olddat, BOOLSET, 0, 5, NULL},
    {"useqr", N_("Use QR decomposition"), 
     N_("Use Cholesky decomposition"), &useqr, BOOLSET, 0, 1, NULL},
    {"Fixed_font", N_("Fixed font"), NULL, fixedfontname, USERSET, MAXLEN, 0, NULL},
#if !defined(USE_GNOME) && !defined(OLD_GTK)
    {"App_font", N_("Menu font"), NULL, appfontname, USERSET, MAXLEN, 0, NULL},
#endif
    {"DataPage", "Default data page", NULL, datapage, INVISET, 
     sizeof datapage, 0, NULL},
    {"ScriptPage", "Default script page", NULL, scriptpage, INVISET, 
     sizeof scriptpage, 0, NULL},    
    {"Png_font", N_("PNG graph font"), NULL, paths.pngfont, INVISET, 16, 0, NULL},
    {"Gp_colors", N_("Gnuplot colors"), NULL, gpcolors, INVISET, 
     sizeof gpcolors, 0, NULL},
    {NULL, NULL, NULL, NULL, 0, 0, 0, NULL}
};

/* ........................................................... */

void set_datapage (const char *str)
{
    strcpy(datapage, str);
}

void set_scriptpage (const char *str)
{
    strcpy(scriptpage, str);
}

const char *get_datapage (void)
{
    return datapage;
}

const char *get_scriptpage (void)
{
    return scriptpage;
}

/* ........................................................... */

void set_fixed_font (void)
{
#ifndef OLD_GTK
    if (fixed_font != NULL) 
	pango_font_description_free(fixed_font);

    fixed_font = pango_font_description_from_string(fixedfontname);
#else
    fixed_font = gdk_font_load(fixedfontname);
#endif
}

/* ........................................................... */

#if !defined(USE_GNOME) && !defined(OLD_GTK)

const char *get_app_fontname (void)
{
    return appfontname;
}

void set_app_font (const char *fontname)
{
    GtkSettings *settings;

    if (fontname != NULL && *fontname == 0) return;

    settings = gtk_settings_get_default();

    if (fontname == NULL) {
	g_object_set(G_OBJECT(settings), "gtk-font-name", appfontname, NULL);
    } else {
	GtkWidget *w;
	PangoFontDescription *pfd;
	PangoContext *pc;
	PangoFont *pfont;

	w = gtk_label_new(NULL);
	pfd = pango_font_description_from_string(fontname);
	pc = gtk_widget_get_pango_context(w);
	pfont = pango_context_load_font(pc, pfd);

	if (pfont != NULL) {
	    strcpy(appfontname, fontname);
	    g_object_set(G_OBJECT(settings), "gtk-font-name", appfontname, NULL);
	}

	gtk_widget_destroy(w);
	pango_font_description_free(pfd);
    }
}

#endif

/* .................................................................. */

static void slash_terminate (char *path)
{
    if (path == NULL || *path == '\0') return;

    if (path[strlen(path) - 1] != SLASH) {
	strcat(path, SLASHSTR);
    }
}

void get_default_dir (char *s)
{
    *s = '\0';

    if (usecwd) {
	char *test = getcwd(s, MAXLEN);

	if (test == NULL) {
	    strcpy(s, paths.userdir);
	} 
    } else {
	strcpy(s, paths.userdir);   
    } 
    
    slash_terminate(s);
}

/* ........................................................... */

#if defined(HAVE_TRAMO) || defined(HAVE_X12A)

# ifdef HAVE_TRAMO
void set_tramo_ok (int set)
{
    static int ok;

    if (set >= 0) ok = set;
    if (mdata != NULL) {
	flip(mdata->ifac, "/Variable/TRAMO analysis", ok);
    }
}
# endif /* HAVE_TRAMO */

# ifdef HAVE_X12A
void set_x12a_ok (int set)
{
    static int ok;

    if (set >= 0) ok = set;
    if (mdata != NULL) {
	flip(mdata->ifac, "/Variable/X-12-ARIMA analysis", ok);
    }
}
# endif /* HAVE_X12A */

# ifdef G_OS_WIN32

static const char *get_reg_base (const char *key)
{
    if (strncmp(key, "x12a", 4) == 0) {
        return "x12arima";
    }
    if (strncmp(key, "tramo", 5) == 0) {
        return "tramo";
    }
    return "gretl";
}

static void set_tramo_x12a_dirs (void)
{
    set_tramo_ok(check_for_prog(tramo));
    set_x12a_ok(check_for_prog(paths.x12a));
}

# else /* not G_OS_WIN32 */

static int my_mkdir (const char *dirname)
{
    int err = 0;
    extern int errno;

    errno = 0;

    if (mkdir(dirname, 0755)) {
	if (errno != EEXIST) { 
	    fprintf(stderr, "%s: %s\n", dirname, strerror(errno));
	    err = 1;
	}
    }
    return err;
}

static void set_tramo_x12a_dirs (void)
{
    char dirname[MAXLEN];
    DIR *test;

#  ifdef HAVE_TRAMO 
    set_tramo_ok(check_for_prog(tramo));
    if (*tramodir == '\0') {
	build_path(paths.userdir, "tramo", tramodir, NULL);
    }
#  endif
#  ifdef HAVE_X12A
    set_x12a_ok(check_for_prog(paths.x12a));
    if (*paths.x12adir == '\0') {
	build_path(paths.userdir, "x12arima", paths.x12adir, NULL);
    }
#  endif

    /* don't make dir structure (yet) if userdir doesn't exist */
    test = opendir(paths.userdir);
    if (test == NULL) {
	return;
    } else {
	closedir(test);
    }
#  ifdef HAVE_X12A
    my_mkdir(paths.x12adir);
#  endif
#  ifdef HAVE_TRAMO
    if (my_mkdir(tramodir)) return;
    sprintf(dirname, "%s/output", tramodir);
    my_mkdir(dirname);
    sprintf(dirname, "%s/graph", tramodir);
    if (my_mkdir(dirname)) return;
    sprintf(dirname, "%s/graph/acf", tramodir);
    my_mkdir(dirname);
    sprintf(dirname, "%s/graph/filters", tramodir);
    my_mkdir(dirname);
    sprintf(dirname, "%s/graph/forecast", tramodir);
    my_mkdir(dirname);
    sprintf(dirname, "%s/graph/series", tramodir);
    my_mkdir(dirname);
    sprintf(dirname, "%s/graph/spectra", tramodir);
    my_mkdir(dirname);
#  endif /* HAVE_TRAMO */
}

# endif /* G_OS_WIN32 */

#endif /* tramo || x12a */

#ifdef G_OS_WIN32

int check_for_prog (const char *prog)
{
    int ret = 1;
    char tmp[MAXLEN];
    WIN32_FIND_DATA find_data;
    HANDLE hfind;

    if (prog == NULL || *prog == 0) return 0;

    hfind = FindFirstFile(prog, &find_data);
    if (hfind == INVALID_HANDLE_VALUE) ret = 0;
    FindClose(hfind);

    if (ret == 0) {
	char *p;

	ret = SearchPath(NULL, prog, NULL, MAXLEN, tmp, &p);
    }

    return ret;
}

#else

int check_for_prog (const char *prog)
{
    char tmp[MAXLEN];
    int ret = 0;

    if (prog == NULL || *prog == 0) return 0;

    if (!strcmp(prog, "latex")) {
	strcpy(tmp, "latex x.tex > /dev/null");
    } else {
	sprintf(tmp, "%s > /dev/null 2>&1", prog);
    }

    ret = gretl_spawn_quiet(tmp) == 0;

    if (!strcmp(prog, "latex")) remove("x.log");

    return ret;
}

#endif

/* ........................................................... */

#if !defined(G_OS_WIN32) && !defined(USE_GNOME)
void set_rcfile (void) 
{
    char *tmp;

    tmp = getenv("HOME");
    strcpy(rcfile, tmp);
# if defined(OSX_PKG)
    strcat(rcfile, "/.gretlosxrc");
# elif defined(OLD_GTK)
    strcat(rcfile, "/.gretlrc");        
# else
    strcat(rcfile, "/.gretl2rc");
# endif
    read_rc(); 
}
#endif

#ifdef USE_GNOME
void set_rcfile (void)
{
    read_rc();
}
#endif

/* .................................................................. */

void options_dialog (gpointer data) 
{
    GtkWidget *tempwid, *dialog, *notebook;

    dialog = gtk_dialog_new();
    gtk_window_set_title(GTK_WINDOW(dialog), _("gretl: options"));
    gtk_container_set_border_width(GTK_CONTAINER(GTK_DIALOG(dialog)->vbox), 10);
    gtk_container_set_border_width 
	(GTK_CONTAINER(GTK_DIALOG(dialog)->action_area), 5);
    gtk_box_set_spacing(GTK_BOX(GTK_DIALOG(dialog)->vbox), 2);
    gtk_box_set_spacing(GTK_BOX(GTK_DIALOG(dialog)->action_area), 15);
    gtk_box_set_homogeneous(GTK_BOX(GTK_DIALOG(dialog)->action_area), TRUE);
    gtk_window_set_position(GTK_WINDOW(dialog), GTK_WIN_POS_MOUSE);

#ifndef OLD_GTK
    g_signal_connect (G_OBJECT(dialog), "delete_event", 
		      G_CALLBACK(delete_widget), 
		      dialog);
#else
    gtk_signal_connect_object(GTK_OBJECT(dialog), "delete_event", 
			      GTK_SIGNAL_FUNC(gtk_widget_destroy), 
			      GTK_OBJECT(dialog));
#endif

    notebook = gtk_notebook_new();
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), notebook, 
		       TRUE, TRUE, 0);
    gtk_widget_show(notebook);

    make_prefs_tab(notebook, 1);
    make_prefs_tab(notebook, 2);
    make_prefs_tab(notebook, 3);
    make_prefs_tab(notebook, 4);
    make_prefs_tab(notebook, 5);
   
    tempwid = standard_button(GTK_STOCK_OK);
    GTK_WIDGET_SET_FLAGS(tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->action_area), 
		       tempwid, TRUE, TRUE, 0);
#ifndef OLD_GTK
    g_signal_connect(G_OBJECT(tempwid), "clicked", 
		     G_CALLBACK(apply_changes), NULL);
    g_signal_connect(G_OBJECT(tempwid), "clicked", 
		     G_CALLBACK(delete_widget), 
		     dialog);
#else
    gtk_signal_connect(GTK_OBJECT(tempwid), "clicked", 
		       GTK_SIGNAL_FUNC(apply_changes), NULL);
    gtk_signal_connect_object(GTK_OBJECT(tempwid), "clicked", 
			      GTK_SIGNAL_FUNC(gtk_widget_destroy), 
			      GTK_OBJECT(dialog));
#endif
    gtk_widget_show(tempwid);

    tempwid = standard_button(GTK_STOCK_CANCEL);
    GTK_WIDGET_SET_FLAGS(tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->action_area), 
		       tempwid, TRUE, TRUE, 0);
#ifndef OLD_GTK
    g_signal_connect(G_OBJECT(tempwid), "clicked", 
		     G_CALLBACK(delete_widget), 
		     dialog);
#else
    gtk_signal_connect_object(GTK_OBJECT(tempwid), "clicked", 
			      GTK_SIGNAL_FUNC(gtk_widget_destroy), 
			      GTK_OBJECT(dialog));
#endif
    gtk_widget_show(tempwid);

    tempwid = standard_button(GTK_STOCK_APPLY);
    GTK_WIDGET_SET_FLAGS(tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->action_area), 
		       tempwid, TRUE, TRUE, 0);
#ifndef OLD_GTK
    g_signal_connect(G_OBJECT(tempwid), "clicked", 
		     G_CALLBACK(apply_changes), NULL);
#else
    gtk_signal_connect(GTK_OBJECT(tempwid), "clicked", 
		       GTK_SIGNAL_FUNC(apply_changes), NULL);
#endif
    gtk_widget_grab_default(tempwid);
    gtk_widget_show(tempwid);

    gtk_widget_show(dialog);
}

/* .................................................................. */

static void flip_sensitive (GtkWidget *w, gpointer data)
{
    GtkWidget *entry = GTK_WIDGET(data);
    
    gtk_widget_set_sensitive(entry, GTK_TOGGLE_BUTTON(w)->active);
}

/* .................................................................. */

void filesel_set_path_callback (const char *setting, char *strvar)
{
    int i = 0;

    while (rc_vars[i].key != NULL) {
	if (rc_vars[i].var == (void *) strvar) {
	    /* FIXME: utf-8 issues here?? */
	    gtk_entry_set_text(GTK_ENTRY(rc_vars[i].widget), 
			       setting);
	    break;
	}
	i++;
    }
}

/* .................................................................. */

static void browse_button_callback (GtkWidget *w, RCVARS *rc)
{
    file_selector(_(rc->description), SET_PATH, rc->var);
}

/* .................................................................. */

static GtkWidget *make_path_browse_button (RCVARS *rc)
{
    GtkWidget *b;

    b = gtk_button_new_with_label(_("Browse..."));
#ifndef OLD_GTK
    g_signal_connect(G_OBJECT(b), "clicked",
		     G_CALLBACK(browse_button_callback), 
		     rc);
#else
    gtk_signal_connect(GTK_OBJECT(b), "clicked",
                       GTK_SIGNAL_FUNC(browse_button_callback), 
		       rc);
#endif
    return b;
}

/* .................................................................. */

static gboolean takes_effect_on_restart (void)
{
    infobox(_("This change will take effect when you restart gretl"));
    return FALSE;
}

/* .................................................................. */

static void make_prefs_tab (GtkWidget *notebook, int tab) 
{
    GtkWidget *box, *b_table, *s_table, *tempwid = NULL;
    int i, s_len, b_len, b_col;
    int s_count = 0, b_count = 0;
    RCVARS *rc = NULL;
   
    box = gtk_vbox_new(FALSE, 0);
    gtk_container_set_border_width(GTK_CONTAINER(box), 10);
    gtk_widget_show(box);

    if (tab == 1)
	tempwid = gtk_label_new(_("General"));
    else if (tab == 2)
	tempwid = gtk_label_new(_("Databases"));
    else if (tab == 3)
	tempwid = gtk_label_new(_("Programs"));
    else if (tab == 4)
	tempwid = gtk_label_new(_("Open/Save path"));
    else if (tab == 5)
	tempwid = gtk_label_new(_("Data files"));
    
    gtk_widget_show(tempwid);
    gtk_notebook_append_page(GTK_NOTEBOOK(notebook), box, tempwid);   

    s_len = 1;
    s_table = gtk_table_new(s_len, 2, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(s_table), 5);
    gtk_table_set_col_spacings(GTK_TABLE(s_table), 5);
    gtk_box_pack_start(GTK_BOX(box), s_table, FALSE, FALSE, 0);
    gtk_widget_show(s_table);

    b_len = b_col = 0;
    b_table = gtk_table_new(1, 2, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(b_table), 5);
    gtk_table_set_col_spacings(GTK_TABLE(b_table), 5);
    gtk_box_pack_start(GTK_BOX(box), b_table, FALSE, FALSE, 10);
    gtk_widget_show(b_table);

    for (i=0; rc_vars[i].key != NULL; i++) {
	rc = &rc_vars[i];
	if (rc->tab == tab) {
	    if ((rc->type & BOOLSET) && rc->link == NULL) { 
		/* simple boolean variable */
		b_count++;

		rc->widget = gtk_check_button_new_with_label(_(rc->description));
		gtk_table_attach_defaults(GTK_TABLE (b_table), rc->widget, 
					  b_col, b_col + 1, b_len, b_len + 1);

		if (*(int *)(rc->var)) {
		    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(rc->widget), TRUE);
		} else {
		    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(rc->widget), FALSE);
		}

		/* special case: warning */
		if (!strcmp(rc->key, "wimp") || !strcmp(rc->key, "lcnumeric")) {
#ifndef OLD_GTK
		    g_signal_connect(G_OBJECT(rc->widget), "toggled",
				     G_CALLBACK(takes_effect_on_restart), 
				     NULL);
#else
		    gtk_signal_connect(GTK_OBJECT(rc->widget), "toggled",
				       GTK_SIGNAL_FUNC(takes_effect_on_restart), 
				       NULL);
#endif
		}

		/* special case: link between toggle and preceding entry */
		if (rc->len && !(rc->type & FIXSET)) {
		    gtk_widget_set_sensitive(rc_vars[i-1].widget,
					     GTK_TOGGLE_BUTTON(rc->widget)->active);
#ifndef OLD_GTK
		    g_signal_connect(G_OBJECT(rc->widget), "clicked",
				     G_CALLBACK(flip_sensitive),
				     rc_vars[i-1].widget);
#else
		    gtk_signal_connect(GTK_OBJECT(rc->widget), "clicked",
				       GTK_SIGNAL_FUNC(flip_sensitive),
				       rc_vars[i-1].widget);
#endif
		} 

		gtk_widget_show(rc->widget);
		b_col++;

		if (b_col == 2) {
		    b_col = 0;
		    b_len++;
		    gtk_table_resize(GTK_TABLE(b_table), b_len + 1, 2);
		}

		if (rc->type & FIXSET) {
		    gtk_widget_set_sensitive(rc->widget, FALSE);
		    if (rc->len) {
			gtk_widget_set_sensitive(rc_vars[i-1].widget, FALSE);
		    }
		}

	    } else if (rc->type & BOOLSET) { 
		/* radio-button dichotomy */
		int val = *(int *)(rc->var);
		GtkWidget *button;
		GSList *group;

		/* do we have some padding to do? */
		if (b_col == 1) {
		    tempwid = gtk_label_new("   ");
		    gtk_table_attach_defaults(GTK_TABLE(b_table), tempwid, 
					      b_col, b_col + 1, b_len, b_len + 1);
		    b_col = 0;
		    b_len++;
		    gtk_table_resize(GTK_TABLE(b_table), b_len + 1, 2);
		}

		b_count++;
		b_len += 3;
		gtk_table_resize(GTK_TABLE(b_table), b_len + 1, 2);

		/* first a separator for the group */
		tempwid = gtk_hseparator_new();
		gtk_table_attach_defaults(GTK_TABLE(b_table), tempwid, 
					  b_col, b_col + 1, b_len - 3, b_len - 2);  
		gtk_widget_show(tempwid);

		/* then a first button */
		button = gtk_radio_button_new_with_label(NULL, _(rc->link));
		gtk_table_attach_defaults 
		    (GTK_TABLE (b_table), button, b_col, b_col + 1, 
		     b_len - 2, b_len - 1);    
		if (!val) {
		    gtk_toggle_button_set_active 
			(GTK_TOGGLE_BUTTON(button), TRUE);
		}
		gtk_widget_show(button);
		group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));

		/* and a second button */
		rc->widget = gtk_radio_button_new_with_label(group, 
							  _(rc->description));
		gtk_table_attach_defaults(GTK_TABLE(b_table), rc->widget, 
					  b_col, b_col + 1, b_len - 1, b_len);  
		if (val) {
		    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(rc->widget), TRUE);
		}
		gtk_widget_show(rc->widget);

		if (rc->type & FIXSET) {
		    gtk_widget_set_sensitive(button, FALSE);
		    gtk_widget_set_sensitive(rc->widget, FALSE);
		}
	    } else if (!(rc->type & INVISET)) { /* string variable */
		s_count++;
		s_len++;
		gtk_table_resize(GTK_TABLE (s_table), s_len, 
				 (tab == 3)? 3 : 2);
		tempwid = gtk_label_new(_(rc->description));
		gtk_misc_set_alignment(GTK_MISC (tempwid), 1, 0.5);
		gtk_table_attach_defaults(GTK_TABLE (s_table), 
					  tempwid, 0, 1, s_len - 1, s_len);
		gtk_widget_show(tempwid);

		rc->widget = gtk_entry_new();
		gtk_table_attach_defaults(GTK_TABLE (s_table), 
					  rc->widget, 1, 2, s_len-1, s_len);
		gtk_entry_set_text(GTK_ENTRY(rc->widget), rc->var);
		gtk_widget_show(rc->widget);

		/* program browse button */
		if (tab == 3 && strstr(rc->description, "directory") == NULL) {
		    tempwid = make_path_browse_button(rc);
		    gtk_table_attach_defaults(GTK_TABLE(s_table), 
					      tempwid, 2, 3, s_len-1, s_len);
		    gtk_widget_show(tempwid);
		}
		if (rc->type & FIXSET) {
		    gtk_widget_set_sensitive(rc->widget, FALSE);
		    gtk_widget_set_sensitive(tempwid, FALSE);
		}
	    }
	} /* end if (rc->tab == tab) */
    } /* end of loop over rc_vars[i].key */

    if (b_count == 0) {
	gtk_widget_destroy(b_table);
    }
    if (s_count == 0) {
	gtk_widget_destroy(s_table);
    }
}

/* .................................................................. */

#ifdef ENABLE_NLS
static void set_lcnumeric (void)
{
    if (lcnumeric) {
#ifdef G_OS_WIN32
	char *lang = getenv("LANG");
	char *set = NULL;

	if (lang != NULL && !strcmp(lang, "es")) {
	    set = setlocale(LC_NUMERIC, "Spanish");
	    if (set == NULL) {
		set = setlocale(LC_NUMERIC, "es");
	    }
	}
	else if (lang != NULL && !strcmp(lang, "fr")) {
	    set = setlocale(LC_NUMERIC, "French");
	    if (set == NULL) {
		set = setlocale(LC_NUMERIC, "fr");
	    }	    
	}
	else if (lang != NULL && !strcmp(lang, "it")) {
	    set = setlocale(LC_NUMERIC, "Italian");
	    if (set == NULL) {
		set = setlocale(LC_NUMERIC, "it");
	    }
	}

	if (set == NULL) {
	    setlocale(LC_NUMERIC, "");
	}
	putenv("LC_NUMERIC=");
#else
	putenv("LC_NUMERIC=");
	setlocale(LC_NUMERIC, "");
#endif
    } else {
	putenv("LC_NUMERIC=C");
	setlocale(LC_NUMERIC, "C");
    }

    reset_local_decpoint();
}
#endif

/* .................................................................. */

static void set_gp_colors (void)
{
    char cstr[4][8];
    int i, nc;

    nc = sscanf(gpcolors, "%7s %7s %7s %7s", 
		cstr[0], cstr[1], cstr[2], cstr[3]);

    for (i=0; i<nc; i++) {
	set_gnuplot_pallette(i, cstr[i]);
    }
}

/* .................................................................. */

static void apply_changes (GtkWidget *widget, gpointer data) 
{
    const gchar *tempstr;
    extern void show_toolbar (void);
    int i = 0;

    for (i=0; rc_vars[i].key != NULL; i++) {
	if (rc_vars[i].widget != NULL) {
	    if (rc_vars[i].type == BOOLSET) {
		if (GTK_TOGGLE_BUTTON(rc_vars[i].widget)->active)
		    *(int *)(rc_vars[i].var) = TRUE;
		else *(int *)(rc_vars[i].var) = FALSE;
	    } 
	    if (rc_vars[i].type == USERSET || rc_vars[i].type == ROOTSET) {
		tempstr = gtk_entry_get_text
		    (GTK_ENTRY(rc_vars[i].widget));
		if (tempstr != NULL && *tempstr != '\0') 
		    strncpy(rc_vars[i].var, tempstr, rc_vars[i].len - 1);
	    }
	}
    }
    
    write_rc();

    if (toolbar_box == NULL && want_toolbar)
	show_toolbar();
    else if (toolbar_box != NULL && !want_toolbar) {
	gtk_widget_destroy(toolbar_box);
	toolbar_box = NULL;
    }

    set_use_qr(useqr);

#if defined(HAVE_TRAMO) || defined(HAVE_X12A)
    set_tramo_x12a_dirs();
#endif

    proxy_init(dbproxy);
}

/* .................................................................. */

#ifndef GNOME2

static void str_to_boolvar (char *s, void *b)
{
    if (strcmp(s, "true") == 0 || strcmp(s, "1") == 0) {
	*(int *)b = TRUE;
    } else {
	*(int *)b = FALSE;
    }	
}

static void boolvar_to_str (void *b, char *s)
{
    if (*(int *)b) strcpy(s, "true");
    else strcpy(s, "false");
}

#endif

/* .................................................................. */

static void initialize_file_lists (void)
{
    int i;

    /* initialize lists of recently opened files */
    for (i=0; i<MAXRECENT; i++) { 
	datalist[i][0] = 0;
	sessionlist[i][0] = 0;
	scriptlist[i][0] = 0;
    }
}

/* .................................................................. */

#ifdef GNOME2

void write_rc (void) 
{
    GConfClient *client;   
    char key[MAXSTR];
    int i;

    client = gconf_client_get_default();

    for (i=0; rc_vars[i].key != NULL; i++) {
	sprintf(key, "/apps/gretl/%s", rc_vars[i].key);
	if (rc_vars[i].type == BOOLSET) {
	    gboolean val = *(gboolean *) rc_vars[i].var;

	    gconf_client_set_bool (client, key, val, NULL);
	} else {
	    gconf_client_set_string (client, key, rc_vars[i].var, NULL);
	}
    }

    g_object_unref(G_OBJECT(client));

    printfilelist(FILE_LIST_DATA, NULL);
    printfilelist(FILE_LIST_SESSION, NULL); 
    printfilelist(FILE_LIST_SCRIPT, NULL);

    set_paths(&paths, 0, 1);
}

static void read_rc (void) 
{
    GConfClient *client;
    GError *error = NULL;
    GSList *flist = NULL;
    gchar *value;
    char key[MAXSTR];
    int i;
    static char *sections[] = {
	"recent_data_files",
	"recent_session_files",
	"recent_script_files"
    };	

    client = gconf_client_get_default();

    for (i=0; rc_vars[i].key != NULL; i++) {
	sprintf(key, "/apps/gretl/%s", rc_vars[i].key);
	if (rc_vars[i].type == BOOLSET) {
	    gboolean val;

	    val = gconf_client_get_bool (client, key, &error);
	    if (error) {
		fprintf(stderr, "Error reading %s\n", rc_vars[i].key);
		g_clear_error(&error);
	    } else {
		*(int *) rc_vars[i].var = val;
	    }
	} else {
	    value = gconf_client_get_string (client, key, &error);

	    if (error) {
		fprintf(stderr, "Error reading %s\n", rc_vars[i].key);
		g_clear_error(&error);
	    } else if (value != NULL) {
		char *s = (char *) rc_vars[i].var;

		*s = 0;
		strncat(s, value, rc_vars[i].len - 1);
		g_free(value);
	    }
	}
    }

    initialize_file_lists();

    for (i=0; i<3; i++) {
	int j;

	sprintf(key, "/apps/gretl/%s", sections[i]);
	flist = gconf_client_get_list (client, key,
				       GCONF_VALUE_STRING, NULL);
	if (flist != NULL) {
	    for (j=0; j<MAXRECENT; j++) {
		if (i == 0) strcpy(datalist[j], flist->data);
		else if (i == 1) strcpy(sessionlist[j], flist->data);
		else if (i == 2) strcpy(scriptlist[j], flist->data);
		flist = flist->next;
	    }
	    g_slist_free(flist);
	    flist = NULL;
	}
    }

    g_object_unref(G_OBJECT(client));

    set_use_qr(useqr);
    set_gp_colors();

    set_paths(&paths, 0, 1); /* 0 = not defaults, 1 = gui */

# if defined(HAVE_TRAMO) || defined(HAVE_X12A)
    set_tramo_x12a_dirs();
# endif

# ifdef ENABLE_NLS
    set_lcnumeric();
# endif
}

#elif defined(USE_GNOME)

void write_rc (void) 
{
    char key[MAXSTR];
    char val[6];
    int i;

    for (i=0; rc_vars[i].key != NULL; i++) {
	sprintf(key, "/gretl/%s/%s", rc_vars[i].description, rc_vars[i].key);
	if (rc_vars[i].type == BOOLSET) {
	    boolvar_to_str(rc_vars[i].var, val);
	    gnome_config_set_string(key, val);
	} else {
	    gnome_config_set_string(key, rc_vars[i].var);
	}
    }

    printfilelist(FILE_LIST_DATA, NULL);
    printfilelist(FILE_LIST_SESSION, NULL);
    printfilelist(FILE_LIST_SCRIPT, NULL);

    gnome_config_sync();

    set_paths(&paths, 0, 1);
}

static void read_rc (void) 
{
    gchar *value = NULL;
    char gpath[MAXSTR];
    int i;

    for (i=0; rc_vars[i].key != NULL; i++) {
	sprintf(gpath, "/gretl/%s/%s", rc_vars[i].description, 
		rc_vars[i].key);
	value = gnome_config_get_string(gpath);
	if (value != NULL) {
	    if (rc_vars[i].type == BOOLSET) {
		str_to_boolvar(value, rc_vars[i].var);
	    } else {
		strncpy(rc_vars[i].var, value, rc_vars[i].len - 1);
	    }
	    g_free(value);
	}
    }

    initialize_file_lists();

    /* get recent file lists */
    for (i=0; i<MAXRECENT; i++) {
	sprintf(gpath, "/gretl/recent data files/%d", i);
	if ((value = gnome_config_get_string(gpath)) != NULL) { 
	    strcpy(datalist[i], value);
	    g_free(value);
	}
	else break;
    }    
    for (i=0; i<MAXRECENT; i++) {
	sprintf(gpath, "/gretl/recent session files/%d", i);
	if ((value = gnome_config_get_string(gpath)) != NULL) { 
	    strcpy(sessionlist[i], value);
	    g_free(value);
	}
	else break;
    } 
    for (i=0; i<MAXRECENT; i++) {
	sprintf(gpath, "/gretl/recent script files/%d", i);
	if ((value = gnome_config_get_string(gpath)) != NULL) { 
	    strcpy(scriptlist[i], value);
	    g_free(value);
	}
	else break;
    }

    set_use_qr(useqr);
    set_gp_colors();

    set_paths(&paths, 0, 1); /* 0 = not defaults, 1 = gui */

# if defined(HAVE_TRAMO) || defined(HAVE_X12A)
    set_tramo_x12a_dirs();
# endif

# ifdef ENABLE_NLS
    set_lcnumeric();
# endif
}

/* end of gnome versions, now win32 */
#elif defined(G_OS_WIN32)

void write_rc (void) 
{
    int i = 0;
    char val[6];

    for (i=0; rc_vars[i].key != NULL; i++) {

	if (rc_vars[i].type & FIXSET) continue;

	if (rc_vars[i].type == BOOLSET) {
	    boolvar_to_str(rc_vars[i].var, val);
	    write_reg_val(HKEY_CURRENT_USER, 
			  "gretl", 
			  rc_vars[i].key, 
			  val);
	} else {
	    if (rc_vars[i].type == ROOTSET) {
		write_reg_val(HKEY_CLASSES_ROOT, 
			      get_reg_base(rc_vars[i].key),
			      rc_vars[i].key, 
			      rc_vars[i].var);
	    } else {
		write_reg_val(HKEY_CURRENT_USER, 
			      get_reg_base(rc_vars[i].key),
			      rc_vars[i].key, 
			      rc_vars[i].var);
	    }
	}
    }

    printfilelist(FILE_LIST_DATA, NULL); 
    printfilelist(FILE_LIST_SESSION, NULL);
    printfilelist(FILE_LIST_SCRIPT, NULL);

    set_paths(&paths, 0, 1);
}

static int get_network_settings (void)
{
    const char *inifile;
    FILE *fp;
    int gotini = 0;

    inifile = get_network_cfg_filename();

    if (*inifile && (fp = fopen(inifile, "r"))) {
	int j, calldrive = tolower(inifile[0]);
	char line[MAXLEN], key[32], linevar[MAXLEN];

	while (fgets(line, MAXLEN, fp)) {
	    int gotvar = 0;
	    char *p = line;

	    while (isspace(*p)) p++;

	    if (*p == '#') continue;
	    if (sscanf(p, "%31s", key) == 1) {
		strcpy(linevar, p + strlen(key) + 3); 
		chopstr(linevar); 
		gotvar = 0;
		for (j=0; rc_vars[j].key != NULL; j++) {
		    if (!strcmp(key, rc_vars[j].key)) {
			if (rc_vars[j].type == BOOLSET) {
			    str_to_boolvar(linevar, rc_vars[j].var);
			} else {
			    if (!strcmp(key, "gretldir") && 
				tolower(linevar[0]) != calldrive) {
				gotini = 0;
				goto network_quit;
			    } else {
				char *var = (char *) rc_vars[j].var;

				*var = '\0';
				strncat(var, linevar, rc_vars[j].len - 1);
			    }
			}
			rc_vars[j].type |= FIXSET;
			gotvar = gotini = 1;
		    }
		    if (gotvar) break;
		}
	    }
	}
    network_quit:
	fclose(fp);
    }

    return gotini;
}

void read_rc (void) 
{
    int i = 0;
    char rpath[MAXSTR], value[MAXSTR];

    if (get_network_settings() && *paths.userdir != '\0') {
	win32_make_user_dirs();
	for (i=0; rc_vars[i].key != NULL; i++) {
	    if (rc_vars[i].var == tramodir ||
		rc_vars[i].var == paths.x12adir) {
		rc_vars[i].type |= FIXSET;
	    }
	}
    } 

    for (i=0; rc_vars[i].key != NULL; i++) {
	int err = 0;

	if (rc_vars[i].type & FIXSET) {
	    continue;
	}

	if (rc_vars[i].type == ROOTSET) {
	    err = read_reg_val (HKEY_CLASSES_ROOT, 
				get_reg_base(rc_vars[i].key),
				rc_vars[i].key, 
				value);
	} else {
	    err = read_reg_val (HKEY_CURRENT_USER, 
				get_reg_base(rc_vars[i].key),
				rc_vars[i].key, 
				value);
	}
	    
	if (!err && *value != '\0') {
	    if (rc_vars[i].type == BOOLSET) {
		str_to_boolvar(value, rc_vars[i].var);
	    } else {
		strncpy((char *) rc_vars[i].var, value, 
			rc_vars[i].len - 1);
	    }
	}
    }

    initialize_file_lists();

    /* get recent file lists */
    for (i=0; i<MAXRECENT; i++) {
	sprintf(rpath, "recent data files\\%d", i);
	if (read_reg_val(HKEY_CURRENT_USER, "gretl", rpath, value) == 0) 
	    strcpy(datalist[i], value);
	else break;
    }    
    for (i=0; i<MAXRECENT; i++) {
	sprintf(rpath, "recent session files\\%d", i);
	if (read_reg_val(HKEY_CURRENT_USER, "gretl", rpath, value) == 0) 
	    strcpy(sessionlist[i], value);
	else break;
    } 
    for (i=0; i<MAXRECENT; i++) {
	sprintf(rpath, "recent script files\\%d", i);
	if (read_reg_val(HKEY_CURRENT_USER, "gretl", rpath, value) == 0) 
	    strcpy(scriptlist[i], value);
	else break;
    }

    set_use_qr(useqr);
    set_gp_colors();

    set_paths(&paths, 0, 1);

# if defined(HAVE_TRAMO) || defined(HAVE_X12A)
    set_tramo_x12a_dirs();
# endif

    set_fixed_font();
    set_app_font(NULL);

# ifdef ENABLE_NLS
    set_lcnumeric();
# endif
}

#else /* end of gnome and win32 versions, now plain GTK */

void write_rc (void) 
{
    FILE *rc;
    char val[6];
    int i;

    rc = fopen(rcfile, "w");
    if (rc == NULL) {
	errbox(_("Couldn't open config file for writing"));
	return;
    }

    fprintf(rc, "# gretl config file (note: not used by gnome version)\n");

    for (i=0; rc_vars[i].var != NULL; i++) {
	fprintf(rc, "# %s\n", rc_vars[i].description);
	if (rc_vars[i].type == BOOLSET) {
	    boolvar_to_str(rc_vars[i].var, val);
	    fprintf(rc, "%s = %s\n", rc_vars[i].key, val);
	} else {
	    fprintf(rc, "%s = %s\n", rc_vars[i].key, (char *) rc_vars[i].var);
	}
    }

    printfilelist(FILE_LIST_DATA, rc);
    printfilelist(FILE_LIST_SESSION, rc);
    printfilelist(FILE_LIST_SCRIPT, rc);

    fclose(rc);
    set_paths(&paths, 0, 1);
}

static void read_rc (void) 
{
    FILE *fp;
    int i, j;
    char line[MAXLEN], key[32], linevar[MAXLEN];
    int gotrecent = 0;

    fp = fopen(rcfile, "r");
    if (fp == NULL) return;

    i = 0;
    while (rc_vars[i].var != NULL) {
	if (fgets(line, MAXLEN, fp) == NULL) break;
	if (line[0] == '#') continue;
	if (!strncmp(line, "recent ", 7)) {
	    gotrecent = 1;
	    break;
	}
	if (sscanf(line, "%s", key) == 1) {
	    strcpy(linevar, line + strlen(key) + 3); 
	    chopstr(linevar); 
	    for (j=0; rc_vars[j].key != NULL; j++) {
		if (!strcmp(key, rc_vars[j].key)) {
		    if (rc_vars[j].type == BOOLSET) {
			str_to_boolvar(linevar, rc_vars[j].var);
		    } else {
			strcpy(rc_vars[j].var, linevar);
		    }
		    break;
		}
	    }
	}
	i++;
    }

    initialize_file_lists();

    if (gotrecent || (fgets(line, MAXLEN, fp) != NULL && 
		      strncmp(line, "recent data files:", 18) == 0)) {
	i = 0;
	while (fgets(line, MAXLEN, fp) && i<MAXRECENT) {
	    if (strncmp(line, "recent session files:", 21) == 0)
		break;
	    chopstr(line);
	    if (*line) strcpy(datalist[i++], line);
	}
    }

    if (strncmp(line, "recent session files:", 21) == 0) {
	i = 0;
	while (fgets(line, MAXLEN, fp) && i<MAXRECENT) {
	    if (strncmp(line, "recent script files:", 20) == 0)
		break;
	    chopstr(line);
	    if (*line) strcpy(sessionlist[i++], line);
	}
    }

    if (strncmp(line, "recent script files:", 20) == 0) {
	i = 0;
	while (fgets(line, MAXLEN, fp) && i<MAXRECENT) {
	    chopstr(line);
	    if (*line) strcpy(scriptlist[i++], line);
	}
    }

    fclose(fp);

    set_use_qr(useqr);
    set_gp_colors();

    set_paths(&paths, 0, 1);

# if defined(HAVE_TRAMO) || defined(HAVE_X12A)
    set_tramo_x12a_dirs();
# endif

# ifdef ENABLE_NLS
    set_lcnumeric();
# endif
}

#endif /* end of "plain gtk" versions of read_rc, write_rc */

/* font selection: non-Windows, gtk-2.0 version first */

#ifndef G_OS_WIN32
# ifndef OLD_GTK

static void font_selection_ok (GtkWidget *w, GtkFontSelectionHackDialog *fs)
{
    guint which = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(fs), "which"));
    gchar *fontname;

    fontname = gtk_font_selection_hack_dialog_get_font_name(fs);

    if (*fontname == 0) {
	g_free(fontname);
	gtk_widget_destroy(GTK_WIDGET(fs));
	return;
    }

    if (which == FIXED_FONT_SELECTION) {
	strcpy(fixedfontname, fontname);
	set_fixed_font();
	write_rc();
    } 

#  ifdef PLOT_FONT_SELECTOR
    else if (which == GRAPH_FONT_SELECTION) {
	GtkWidget *fentry = g_object_get_data(G_OBJECT(fs), "font_entry");

	gtk_entry_set_text(GTK_ENTRY(fentry), fontname);
    }
#  endif

#  ifndef USE_GNOME /* gnome handles the app font */
    else if (which == APP_FONT_SELECTION) {
	set_app_font(fontname);
	write_rc();
    }
#  endif

    g_free(fontname);
    gtk_widget_destroy(GTK_WIDGET(fs));
}

static void fontsel_quit (GtkWidget *w, gpointer p)
{
    gtk_main_quit();
}

void font_selector (gpointer data, guint which, GtkWidget *widget)
{
    static GtkWidget *fontsel = NULL;
    int filter = GTK_FONT_HACK_LATIN;
    char *title = NULL;
    const char *fontname = NULL;

#  ifdef USE_GNOME
    if (which == APP_FONT_SELECTION) return; /* shouldn't happen */
#  endif

    if (fontsel != NULL) {
	if (!GTK_WIDGET_VISIBLE(fontsel)) gtk_widget_show (fontsel);
        gdk_window_raise(fontsel->window);
        return;
    }

    if (which == FIXED_FONT_SELECTION) {
	title = _("Font for gretl output windows");
	filter = GTK_FONT_HACK_LATIN_MONO;
	fontname = fixedfontname;
    }
#  ifndef USE_GNOME
    else if (which == APP_FONT_SELECTION) {
	title = _("Font for menus and labels");
	fontname = appfontname;
    }
#  endif
#  ifdef PLOT_FONT_SELECTOR
    else if (which == GRAPH_FONT_SELECTION) {
	fontname = paths.pngfont;
    }
#  endif

    fontsel = gtk_font_selection_hack_dialog_new(title);
    gtk_font_selection_hack_dialog_set_filter
	(GTK_FONT_SELECTION_HACK_DIALOG (fontsel), filter);
    gtk_font_selection_hack_dialog_set_font_name 
	(GTK_FONT_SELECTION_HACK_DIALOG (fontsel), fontname); 
    g_object_set_data(G_OBJECT(fontsel), "which", GINT_TO_POINTER(which));

#  ifdef PLOT_FONT_SELECTOR
    if (which == GRAPH_FONT_SELECTION) {
	 g_object_set_data(G_OBJECT(fontsel), "font_entry", widget);
    }
#  endif

    gtk_window_set_position (GTK_WINDOW (fontsel), GTK_WIN_POS_MOUSE);

    g_signal_connect (G_OBJECT(fontsel), "destroy",
		      G_CALLBACK(gtk_widget_destroyed),
		      &fontsel);
    g_signal_connect (G_OBJECT(fontsel), "destroy",
		      G_CALLBACK(fontsel_quit),
		      NULL);
    g_signal_connect (G_OBJECT(GTK_FONT_SELECTION_HACK_DIALOG(fontsel)->ok_button),
		      "clicked", 
		      G_CALLBACK(font_selection_ok),
		      fontsel);
    g_signal_connect(G_OBJECT(GTK_FONT_SELECTION_HACK_DIALOG(fontsel)->cancel_button),
		     "clicked", 
		     G_CALLBACK(delete_widget),
		     fontsel);

    gtk_widget_show (fontsel);

    gtk_main();
}

# else /* done gtk 2, now gtk 1.2 */

static void font_selection_ok (GtkWidget *w, GtkFontSelectionDialog *fs)
{
    gchar *fstring = gtk_font_selection_dialog_get_font_name(fs);

    if (strlen(fstring)) {
        strcpy(fixedfontname, fstring);
        gdk_font_unref(fixed_font);
        fixed_font = gdk_font_load(fixedfontname);
        write_rc();
    }
    g_free(fstring);
    gtk_widget_destroy(GTK_WIDGET (fs));
}

/* .................................................................. */

void font_selector (gpointer data, guint u, GtkWidget *w)
{
    static GtkWidget *fontsel = NULL;
    gchar *spacings[] = { "c", "m", NULL };

    if (!fontsel) {
	fontsel = gtk_font_selection_dialog_new 
	    (_("Font for gretl output windows"));

	gtk_window_set_position (GTK_WINDOW (fontsel), GTK_WIN_POS_MOUSE);

	gtk_font_selection_dialog_set_filter(GTK_FONT_SELECTION_DIALOG (fontsel),
					     GTK_FONT_FILTER_BASE, GTK_FONT_ALL,
					     NULL, NULL, NULL, NULL, 
					     spacings, NULL);

	gtk_font_selection_dialog_set_font_name 
	    (GTK_FONT_SELECTION_DIALOG (fontsel), fixedfontname);

	gtk_signal_connect (GTK_OBJECT(fontsel), "destroy",
			    GTK_SIGNAL_FUNC(gtk_widget_destroyed),
			    &fontsel);

	gtk_signal_connect (GTK_OBJECT 
			    (GTK_FONT_SELECTION_DIALOG 
			     (fontsel)->ok_button),
			    "clicked", GTK_SIGNAL_FUNC(font_selection_ok),
			    GTK_FONT_SELECTION_DIALOG (fontsel));

	gtk_signal_connect_object (GTK_OBJECT 
				   (GTK_FONT_SELECTION_DIALOG 
				    (fontsel)->cancel_button),
				   "clicked", 
				   GTK_SIGNAL_FUNC(gtk_widget_destroy),
				   GTK_OBJECT (fontsel));
    }
    if (!GTK_WIDGET_VISIBLE (fontsel)) gtk_widget_show (fontsel);
    else gtk_widget_destroy (fontsel);
}

# endif /* non-Windows, gtk version branches */

#else /* end non-win32 font selection, start win32 */

static const char *font_weight_string (int weight)
{
    if (weight >= FW_THIN && weight <= FW_LIGHT) return " Light";
    if (weight >= FW_NORMAL && weight <= FW_DEMIBOLD) return "";
    if (weight >= FW_BOLD) return " Bold";
    return "";
}

static void fontname_to_win32 (const char *src, int fixed,
			       char *name, int *pts)
{
    int sz;

    if (sscanf(src, "%31[^0123456789]%d", name, &sz) == 2) {
	size_t i, n = strlen(name);

	for (i=n-1; i>0; i--) {
	    if (name[i] == ' ') name[i] = 0;
	    else break;
	}
	*pts = sz * 10; /* measured in tenths of a point */
    } else {
	*name = 0;
	strncat(name, src, 31);
	if (fixed) *pts = 100;
	else *pts = 80;
    }
}

void font_selector (gpointer data, guint which, GtkWidget *widget)
{
    CHOOSEFONT cf;            /* common dialog box structure */
    LOGFONT lf;               /* logical font structure */
    char fontname[48];

    ZeroMemory(&cf, sizeof cf);
    cf.lStructSize = sizeof cf;
    cf.Flags = CF_SCREENFONTS | CF_TTONLY | CF_LIMITSIZE | CF_INITTOLOGFONTSTRUCT;
    cf.nSizeMin = 6;
    cf.nSizeMax = 24;

    ZeroMemory(&lf, sizeof lf);

    cf.Flags |= CF_NOSCRIPTSEL;
    
    lf.lfWeight = FW_REGULAR;
    lf.lfCharSet = DEFAULT_CHARSET;

    if (which == FIXED_FONT_SELECTION) {
	cf.Flags |= CF_FIXEDPITCHONLY;
	fontname_to_win32(fixedfontname, 1, lf.lfFaceName, &(cf.iPointSize));
    } 
    else if (which == APP_FONT_SELECTION) {
	fontname_to_win32(appfontname, 0, lf.lfFaceName, &(cf.iPointSize));
    } 
# ifdef PLOT_FONT_SELECTOR
    else if (which == GRAPH_FONT_SELECTION) {
	fontname_to_win32(paths.pngfont, 0, lf.lfFaceName, &(cf.iPointSize));
    }
# endif

    cf.lpLogFont = &lf;

    if (ChooseFont(&cf) == TRUE && *(cf.lpLogFont->lfFaceName)) {
	sprintf(fontname, "%s%s %d", cf.lpLogFont->lfFaceName, 
		font_weight_string(cf.lpLogFont->lfWeight), 
		cf.iPointSize / 10);
	if (which == FIXED_FONT_SELECTION) {
	    strcpy(fixedfontname, fontname);
	    set_fixed_font();
	    write_rc();
	} 
	else if (which == APP_FONT_SELECTION) {
	    set_app_font(fontname);
	    write_rc();
	}
# ifdef PLOT_FONT_SELECTOR
	else if (which == GRAPH_FONT_SELECTION) {
	    gtk_entry_set_text(GTK_ENTRY(widget), fontname);
	}
# endif
    }
}

#endif /* end win32 font selection */

/* .................................................................. */

void init_fileptrs (void)
{
    int i;
    
    for (i=0; i<MAXRECENT; i++) {
	datap[i] = datalist[i];
	sessionp[i] = sessionlist[i];
	scriptp[i] = scriptlist[i];
    }
}

/* .................................................................. */

static void clear_files_list (int filetype, char **filep)
{
    GtkWidget *w;
    char tmpname[MAXSTR];
    gchar itempath[80];
    int i, pindex = -1;
    const gchar *fpath[] = {
	N_("/File/Open data"), 
	N_("/Session"),
	N_("/File/Open command file")
    };

    if (filetype == FILE_LIST_DATA) 
	pindex = 0;
    else if (filetype == FILE_LIST_SESSION)
	pindex = 1;
    else if (filetype == FILE_LIST_SCRIPT)
	pindex = 2;

    if (pindex == -1) return;

    for (i=0; i<MAXRECENT; i++) {
#ifndef OLD_GTK
	sprintf(itempath, "%s/%d. %s", fpath[pindex],
		i+1, endbit(tmpname, filep[i], 0)); 
#else
	sprintf(itempath, "%s/%d. %s", fpath[pindex],
		i+1, endbit(tmpname, filep[i], -1)); 
#endif
	w = gtk_item_factory_get_widget(mdata->ifac, itempath);
	if (w != NULL) {
	    gtk_item_factory_delete_item(mdata->ifac, itempath);
	}
    }
}

/* .................................................................. */

static char *cut_multiple_slashes (char *fname)
{
    int i, n = strlen(fname);
#ifdef G_OS_WIN32
    /* may be ok for a filename to start with a double backslash */
    int start = 1;
#else
    int start = 0;
#endif

    for (i=start; i<n-1; i++) {
	if (fname[i] == SLASH && fname[i+1] == SLASH) {
	    memmove(fname + i, fname + i + 1, strlen(fname + i + 1) + 1);
	    i--;
	    n--;
	}
    }

    return fname;
}

void mkfilelist (int filetype, char *fname)
{
    char *tmp[MAXRECENT-1];
    char **filep;
    int i, match = -1;

    cut_multiple_slashes(fname);

    if (filetype == FILE_LIST_DATA) filep = datap;
    else if (filetype == FILE_LIST_SESSION) filep = sessionp;
    else if (filetype == FILE_LIST_SCRIPT) filep = scriptp;
    else return;

    /* see if this file is already on the list */
    for (i=0; i<MAXRECENT; i++) {
        if (strcmp(filep[i], fname) == 0) {
            match = i;
            break;
        }
    }
    if (match == 0) return; /* file is on top: no change in list */

    /* clear menu files list before rebuilding */
    clear_files_list(filetype, filep);
    
    /* save pointers to current order */
    for (i=0; i<MAXRECENT-1; i++) tmp[i] = filep[i];

    /* copy fname into array, if not already present */
    if (match == -1) {
        for (i=1; i<MAXRECENT; i++) {
            if (filep[i][0] == '\0') {
                strcpy(filep[i], fname);
                match = i;
                break;
	    }
	    if (match == -1) {
		match = MAXRECENT - 1;
		strcpy(filep[match], fname);
	    }
	}
    } 

    /* set first pointer to new file */
    filep[0] = filep[match];

    /* rearrange other pointers */
    for (i=1; i<=match; i++) {
	filep[i] = tmp[i-1];
    }

    add_files_to_menu(filetype);
}

/* .................................................................. */

void delete_from_filelist (int filetype, const char *fname)
{
    char *tmp[MAXRECENT];
    char **filep;
    int i, match = -1;

    if (filetype == FILE_LIST_DATA) filep = datap;
    else if (filetype == FILE_LIST_SESSION) filep = sessionp;
    else if (filetype == FILE_LIST_SCRIPT) filep = scriptp;
    else return;

    /* save pointers to current order */
    for (i=0; i<MAXRECENT; i++) {
	tmp[i] = filep[i];
	if (!strcmp(filep[i], fname)) match = i;
    }

    if (match == -1) return;

    /* clear menu files list before rebuilding */
    clear_files_list(filetype, filep);

    for (i=match; i<MAXRECENT-1; i++) {
	filep[i] = tmp[i+1];
    }

    filep[MAXRECENT-1] = tmp[match];
    filep[MAXRECENT-1][0] = '\0';

    add_files_to_menu(filetype);
    /* need to save to file at this point? */
}

/* .................................................................. */

char *endbit (char *dest, char *src, int addscore)
{
    /* take last part of src filename */
    if (strrchr(src, SLASH))
	strcpy(dest, strrchr(src, SLASH) + 1);
    else
	strcpy(dest, src);

    if (addscore != 0) {
	/* then either double (1) or delete (-1) any underscores */
	char mod[MAXSTR];
	size_t i, j, n;

	n = strlen(dest);
	j = 0;
	for (i=0; i<=n; i++) {
	    if (dest[i] != '_')
		mod[j++] = dest[i];
	    else {
		if (addscore == 1) {
		    mod[j++] = '_';
		    mod[j++] = dest[i];
		} 
	    }
	}
	strcpy(dest, mod);
    }
    return dest;
}

/* .................................................................. */

#ifdef GNOME2

static void printfilelist (int filetype, FILE *fp)
     /* param fp is ignored */
{
    GConfClient *client;
    GSList *flist = NULL;
    gchar *key;
    int i;
    char **filep;
    static char *sections[] = {
	"recent_data_files",
	"recent_session_files",
	"recent_script_files"
    };

    switch (filetype) {
    case FILE_LIST_DATA: filep = datap; break;
    case FILE_LIST_SESSION: filep = sessionp; break;
    case FILE_LIST_SCRIPT: filep = scriptp; break;
    default: return;
    }

    client = gconf_client_get_default();

    for (i=0; i<MAXRECENT; i++) {
	flist = g_slist_append (flist, g_strdup(filep[i]));
    }

    key = g_strdup_printf("/apps/gretl/%s", sections[filetype - 1]);

    gconf_client_set_list (client, key, GCONF_VALUE_STRING, 
			   flist, NULL);

    g_free(key);
    g_slist_free(flist);
    g_object_unref(G_OBJECT(client));
}

#elif defined(USE_GNOME)

static void printfilelist (int filetype, FILE *fp)
     /* fp is ignored */
{
    int i;
    char **filep;
    char gpath[MAXLEN];
    static char *section[] = {
	"recent data files",
	"recent session files",
	"recent script files"
    };

    switch (filetype) {
    case FILE_LIST_DATA: filep = datap; break;
    case FILE_LIST_SESSION: filep = sessionp; break;
    case FILE_LIST_SCRIPT: filep = scriptp; break;
    default: return;
    }

    for (i=0; i<MAXRECENT; i++) {
	sprintf(gpath, "/gretl/%s/%d", section[filetype - 1], i);
	gnome_config_set_string(gpath, filep[i]);
    }
}

#elif defined(G_OS_WIN32)

static void printfilelist (int filetype, FILE *fp)
     /* param fp is ignored */
{
    int i;
    char **filep;
    char rpath[MAXLEN];
    static char *sections[] = {
	"recent data files",
	"recent session files",
	"recent script files"
    };

    switch (filetype) {
    case FILE_LIST_DATA: filep = datap; break;
    case FILE_LIST_SESSION: filep = sessionp; break;
    case FILE_LIST_SCRIPT: filep = scriptp; break;
    default: return;
    }

    for (i=0; i<MAXRECENT; i++) {
	if (filep[i] == NULL) continue;
	sprintf(rpath, "%s\\%d", sections[filetype - 1], i);
	write_reg_val(HKEY_CURRENT_USER, "gretl", rpath, filep[i]);
    }
}

#else /* "plain" version follows */

static void printfilelist (int filetype, FILE *fp)
{
    int i;
    char **filep;

    if (filetype == FILE_LIST_DATA) {
	fprintf(fp, "recent data files:\n");
	filep = datap;
    } else if (filetype == FILE_LIST_SESSION) {
	fprintf(fp, "recent session files:\n");
	filep = sessionp;
    } else if (filetype == FILE_LIST_SCRIPT) {
	fprintf(fp, "recent script files:\n");
	filep = scriptp;
    } else 
	return;

    for (i=0; i<MAXRECENT; i++) {
	if (filep[i][0]) 
	    fprintf(fp, "%s\n", filep[i]);
	else break;
    }
}

#endif 

/* ........................................................... */

static void set_data_from_filelist (gpointer data, guint i, 
				    GtkWidget *widget)
{
    strcpy(trydatfile, datap[i]);
    if (strstr(trydatfile, ".csv")) delimiter_dialog();
    verify_open_data(NULL, 0);
}

/* ........................................................... */

static void set_session_from_filelist (gpointer data, guint i, 
				       GtkWidget *widget)
{
    strcpy(tryscript, sessionp[i]);
    verify_open_session(NULL);
}

/* ........................................................... */

static void set_script_from_filelist (gpointer data, guint i, 
				      GtkWidget *widget)
{
    strcpy(tryscript, scriptp[i]);
    do_open_script();
}

/* .................................................................. */

void add_files_to_menu (int filetype)
{
    int i;
    char **filep, tmp[MAXSTR];
    void (*callfunc)();
    GtkItemFactoryEntry fileitem;
    GtkWidget *w;
    const gchar *msep[] = {
	N_("/File/Open data/sep"),
	N_("/Session/sep"),
	N_("/File/Open command file/sep")
    };
    const gchar *mpath[] = {
	N_("/File/Open data"),
	N_("/Session"),
	N_("/File/Open command file")
    };

    fileitem.path = NULL;

    if (filetype == FILE_LIST_DATA) {
	callfunc = set_data_from_filelist;
	filep = datap;
    } else if (filetype == FILE_LIST_SESSION) {
	callfunc = set_session_from_filelist;
	filep = sessionp;
    } else if (filetype == FILE_LIST_SCRIPT) {
	callfunc = set_script_from_filelist;
	filep = scriptp;
    }
    else return;

    /* See if there are any files to add */
    if (filep[0][0] == '\0') {
	return;
    } else {
	gchar *itemtype = "<Separator>";
	GtkWidget *w;

	/* is a separator already in place? */
	w = gtk_item_factory_get_widget(mdata->ifac, msep[filetype - 1]);
	if (w == NULL) {
	    fileitem.path = g_strdup(msep[filetype - 1]);
	    fileitem.accelerator = NULL;
	    fileitem.callback = NULL;
	    fileitem.callback_action = 0;
	    fileitem.item_type = itemtype;
	    gtk_item_factory_create_item(mdata->ifac, &fileitem, NULL, 1);
	    g_free(fileitem.path);
	}
    }

    /* put the files under the menu separator */
    for (i=0; i<MAXRECENT; i++) {
	if (filep[i][0]) {
	    fileitem.accelerator = NULL;
	    fileitem.callback_action = i; 
	    fileitem.item_type = NULL;
	    fileitem.path = g_strdup_printf("%s/%d. %s", mpath[filetype - 1],
					    i+1, endbit(tmp, filep[i], 1));
	    fileitem.callback = callfunc; 
	    gtk_item_factory_create_item(mdata->ifac, &fileitem, NULL, 1);
	    g_free(fileitem.path);
	    w = gtk_item_factory_get_widget_by_action(mdata->ifac, i);
	    if (w != NULL) {
		gretl_tooltips_add(w, filep[i]);
	    } 
	} else break;
    }
}

/* graph color selection apparatus */

static double scale_round (double val)
{
#ifndef OLD_GTK
    return val * 255.0 / 65535.0;
#else
    return val * 255.0;
#endif
}

#define XPMROWS 19
#define XPMCOLS 17

static GtkWidget *get_image_for_color (const char *colstr)
{
#ifndef OLD_GTK
    GdkPixbuf *icon;
#else
    GdkPixmap *pixmap;
    GdkBitmap *mask;
    GtkStyle *style;
#endif
    GtkWidget *image;
    static char **xpm = NULL;
    int i;

    if (xpm == NULL) {
	xpm = malloc(XPMROWS * sizeof *xpm);
	if (xpm != NULL) {
	    for (i=0; i<XPMROWS; i++) {
		xpm[i] = malloc(XPMCOLS * sizeof **xpm);
		if (xpm[i] == NULL) {
		    int j;

		    for (j=0; j<i; j++) free(xpm[j]);
		    free(xpm);
		    xpm = NULL;
		}
		if (i == 0) {
		    strcpy(xpm[i], "16 16 2 1");
		} else if (i == 1) {
		    strcpy(xpm[i], "X      c #000000");
		} else if (i == 2) {
		    strcpy(xpm[i], ".      c #000000");
		} else if (i == 3 || i == XPMROWS - 1) {
		    strcpy(xpm[i], "................");
		} else {
		    strcpy(xpm[i], ".XXXXXXXXXXXXXX.");
		}
	    }
	}
    }

    if (xpm == NULL) return NULL;

    for (i=0; i<6; i++) {
	xpm[1][10+i] = colstr[1+i];
    }    

#ifndef OLD_GTK
    icon = gdk_pixbuf_new_from_xpm_data((const char **) xpm);
    image = gtk_image_new_from_pixbuf(icon);
#else
    style = gtk_widget_get_style(mdata->w);
    pixmap = gdk_pixmap_create_from_xpm_d(mdata->w->window,
					  &mask, 
					  &style->bg[GTK_STATE_NORMAL], 
					  xpm);
    image = gtk_pixmap_new(pixmap, mask);
#endif
    
    return image;
}

static void color_select_callback (GtkWidget *button, GtkWidget *w)
{
    GtkWidget *csel;
    GtkWidget *color_button, *image;
#ifndef OLD_GTK
    GdkColor color;
#else
    gdouble color[4];
#endif
    char color_string[12];
    gint i;

#ifndef OLD_GTK
    color_button = g_object_get_data(G_OBJECT(w), "color_button");
#else
    color_button = gtk_object_get_data(GTK_OBJECT(w), "color_button");
#endif

    csel = GTK_COLOR_SELECTION_DIALOG(w)->colorsel;

#ifndef OLD_GTK
    gtk_color_selection_get_current_color(GTK_COLOR_SELECTION(csel), &color);
    sprintf(color_string, "x%02x%02x%02x",
	    (guint) (scale_round (color.red)),
	    (guint) (scale_round (color.green)),
	    (guint) (scale_round (color.blue)));

    i = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "colnum"));
#else
    gtk_color_selection_get_color(GTK_COLOR_SELECTION(csel), color);
    sprintf(color_string, "x%02x%02x%02x",
	    (guint) (scale_round (color[0])),
	    (guint) (scale_round (color[1])),
	    (guint) (scale_round (color[2])));

    i = GPOINTER_TO_INT(gtk_object_get_data(GTK_OBJECT(w), "colnum"));
#endif

    set_gnuplot_pallette(i, color_string);

    sprintf(gpcolors, "%s %s %s %s", 
	    get_gnuplot_pallette(0, 0),
	    get_gnuplot_pallette(1, 0),
	    get_gnuplot_pallette(2, 0),
	    get_gnuplot_pallette(3, 0));

#ifndef OLD_GTK
    image = g_object_get_data(G_OBJECT(color_button), "image");
#else
    image = gtk_object_get_data(GTK_OBJECT(color_button), "image");
#endif
    gtk_widget_destroy(image);
    image = get_image_for_color(color_string);
    gtk_widget_show(image);
    gtk_container_add(GTK_CONTAINER(color_button), image);
  
    gtk_widget_destroy(w);
}

static void color_cancel (GtkWidget *button, GtkWidget *w)
{
    gtk_widget_destroy(w);
}

GtkWidget *color_patch_button (int colnum)
{
    GtkWidget *image, *button;
    const char *colstr;

    if (colnum == COLOR_MAX) {
	colstr = get_gnuplot_pallette(0, PLOT_FREQ_SIMPLE);
    } else {
	colstr = get_gnuplot_pallette(colnum, 0);
    }

    image = get_image_for_color(colstr);

    if (image == NULL) {
	button = gtk_button_new_with_label(_("Select color"));
    } else {
	button = gtk_button_new();
	gtk_container_add(GTK_CONTAINER(button), image);
#ifndef OLD_GTK
	g_object_set_data(G_OBJECT(button), "image", image);
#else
	gtk_object_set_data(GTK_OBJECT(button), "image", image);
#endif
    }	

    return button;
}

#ifdef OLD_GTK

static int colstr_to_color (const char *colstr, gdouble *color)
{
    char s[3];
    int ci, i;

    for (i=0; i<3; i++) {
	*s = '\0';
	strncat(s, colstr + 2*i + 1, 2);
	if (sscanf(s, "%x", &ci) == 1) {
	    color[i] = ci / 255.0;
	} else {
	    color[i] = 0.0;
	}
    }

    return 0;
}

void gnuplot_color_selector (GtkWidget *w, gpointer p)
{
    GtkWidget *cdlg;
    GtkWidget *button;
    gint i = GPOINTER_TO_INT(p);
    gdouble color[4];
    const gchar *colstr;

    if (i == COLOR_MAX) {
	colstr = get_gnuplot_pallette(0, PLOT_FREQ_SIMPLE);
    } else {
	colstr = get_gnuplot_pallette(i, 0); 
    }

    colstr_to_color(colstr, color);

    cdlg = gtk_color_selection_dialog_new("gretl color selection");

    gtk_object_set_data(GTK_OBJECT(cdlg), "colnum", GINT_TO_POINTER(i));
    gtk_object_set_data(GTK_OBJECT(cdlg), "color_button", w);

    gtk_color_selection_set_color(GTK_COLOR_SELECTION
				  (GTK_COLOR_SELECTION_DIALOG(cdlg)->colorsel),
				  color);

    button = GTK_COLOR_SELECTION_DIALOG(cdlg)->ok_button;
    gtk_signal_connect(GTK_OBJECT(button), "clicked", 
		       GTK_SIGNAL_FUNC(color_select_callback), cdlg);

    button = GTK_COLOR_SELECTION_DIALOG(cdlg)->cancel_button;
    gtk_signal_connect(GTK_OBJECT(button), "clicked", 
		       GTK_SIGNAL_FUNC(color_cancel), cdlg);
    
    gtk_widget_show(cdlg);
}

#else  /* !OLD_GTK */

void gnuplot_color_selector (GtkWidget *w, gpointer p)
{
    GtkWidget *cdlg;
    GtkWidget *button;
    gint i = GPOINTER_TO_INT(p);
    const char *colstr;
    char my_colstr[8];
    GdkColor color;

    if (i == COLOR_MAX) {
	colstr = get_gnuplot_pallette(0, PLOT_FREQ_SIMPLE);
    } else {
	colstr = get_gnuplot_pallette(i, 0); 
    } 

    strcpy(my_colstr, colstr);
    *my_colstr = '#';
    gdk_color_parse(my_colstr, &color);

    cdlg = gtk_color_selection_dialog_new("gretl color selection");

    g_object_set_data(G_OBJECT(cdlg), "colnum", GINT_TO_POINTER(i));
    g_object_set_data(G_OBJECT(cdlg), "color_button", w);

    gtk_color_selection_set_current_color(GTK_COLOR_SELECTION
					  (GTK_COLOR_SELECTION_DIALOG(cdlg)->colorsel),
					  &color);					  

    button = GTK_COLOR_SELECTION_DIALOG(cdlg)->ok_button;
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(color_select_callback), cdlg);

    button = GTK_COLOR_SELECTION_DIALOG(cdlg)->cancel_button;
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(color_cancel), cdlg);
    
    gtk_widget_show(cdlg);
}

#endif /* gtk versions */

/* end graph color selection apparatus */

#ifndef G_OS_WIN32

static int validate_dir (const char *dirname)
{
    DIR *test;
    int err = 0;

    test = opendir(dirname);

    if (test != NULL) {
	fprintf(stderr, "Working dir already exists, OK\n");
	closedir(test);
    } else {
	err = mkdir(dirname, 0755);
	if (err) {
	    sprintf(errtext, _("Couldn't create directory '%s'"), dirname);
	    errbox(errtext);
	} else {
	    infobox(_("Working directory created OK"));
	}
    }

    return err;
}

static void real_set_userdir (GtkWidget *widget, dialog_t *ddata)
{
    const gchar *dirname;

    dirname = gtk_entry_get_text (GTK_ENTRY (ddata->edit));

    if (validate_dir(dirname)) {
	return;
    } else {
	size_t n = strlen(dirname);

	strcpy(paths.userdir, dirname);
	if (dirname[n - 1] != '/') {
	    strcat(paths.userdir, "/");
	}
#if defined(HAVE_TRAMO) || defined(HAVE_X12A)
	set_tramo_x12a_dirs();
#endif
	close_dialog(ddata);
    }
}

void first_time_set_user_dir (void)
{
    DIR *test;

    /* see if the already-specified userdir exists */
    if (*paths.userdir != '\0') {
	test = opendir(paths.userdir);
	if (test != NULL) {
	    closedir(test);
	    return;
	}
    }
	
    /* user dir is not specified, or doesn't exist */
    edit_dialog (_("gretl: working directory"), 
                 _("You seem to be using gretl for the first time.\n"
		   "Please enter a directory for gretl user files."),
                 paths.userdir, 
                 real_set_userdir, NULL, 
                 CREATE_USERDIR, 0);
}

#endif /* G_OS_WIN32 */
