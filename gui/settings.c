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
#include <unistd.h>

#ifndef USE_GNOME
char rcfile[MAXLEN];
#endif

extern GtkTooltips *gretl_tips;
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

#ifdef HAVE_X12A
extern char x12a[MAXSTR];
extern char x12adir[MAXSTR];
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
static void read_rc (void);

/* font handling */
static char fontspec[MAXLEN] = 
"-b&h-lucidatypewriter-medium-r-normal-sans-12-*-*-*-*-*-*-*";
GdkFont *fixed_font;

static int usecwd;
int olddat;
int jwdata;
#ifdef ENABLE_NLS
static int lcnumeric = 1;
#endif

typedef struct {
    char *key;         /* config file variable name */
    char *description; /* How the field will show up in the options dialog */
    char *link;        /* in case of radio button pair, alternate string */
    void *var;         /* pointer to variable */
    char type;         /* 'U' (user) or 'R' (root) for string, 'B' for boolean */
    int len;           /* storage size for string variable (also see Note) */
    short tab;         /* which tab (if any) does the item fall under? */
    GtkWidget *widget;
} RCVARS;

/* Note: actually "len" above is overloaded: if an rc_var is of type 'B'
   (boolean) and not part of a radio group, then a non-zero value for
   len will link the var's toggle button with the sensitivity of the
   preceding rc_var's entry field.  For example, the "use_proxy" button
   controls the sensitivity of the "dbproxy" entry widget. */

RCVARS rc_vars[] = {
    {"gretldir", N_("Main gretl directory"), NULL, paths.gretldir, 
     'R', MAXLEN, 1, NULL},
    {"userdir", N_("User's gretl directory"), NULL, paths.userdir, 
     'U', MAXLEN, 1, NULL},
    {"expert", N_("Expert mode (no warnings)"), NULL, &expert, 
     'B', 0, 1, NULL},
    {"updater", N_("Tell me about gretl updates"), NULL, &updater, 
     'B', 0, 1, NULL},
    {"toolbar", N_("Show gretl toolbar"), NULL, &want_toolbar, 
     'B', 0, 1, NULL},
#ifdef ENABLE_NLS
    {"lcnumeric", N_("Use locale setting for decimal point"), NULL, &lcnumeric, 
     'B', 0, 1, NULL},
#endif
    {"gnuplot", N_("Command to launch gnuplot"), NULL, paths.gnuplot, 
     'R', MAXLEN, 3, NULL},
    {"Rcommand", N_("Command to launch GNU R"), NULL, Rcommand, 
     'R', MAXSTR, 3, NULL},
    {"viewdvi", N_("Command to view DVI files"), NULL, viewdvi, 
     'R', MAXSTR, 3, NULL},
    {"calculator", N_("Calculator"), NULL, calculator, 
     'U', MAXSTR, 3, NULL},
    {"editor", N_("Editor"), NULL, editor, 
     'U', MAXSTR, 3, NULL},
#ifdef HAVE_X12A
    {"x12a", N_("path to x12arima"), NULL, x12a, 
     'R', MAXSTR, 3, NULL},
#endif
#ifdef HAVE_TRAMO
    {"tramo", N_("path to tramo"), NULL, tramo, 
     'R', MAXSTR, 3, NULL},
#endif
#ifdef HAVE_X12A
    {"x12adir", N_("X-12-ARIMA working directory"), NULL, x12adir, 
     'R', MAXSTR, 3, NULL},
#endif
#ifdef HAVE_TRAMO
    {"tramodir", N_("TRAMO working directory"), NULL, tramodir, 
     'R', MAXSTR, 3, NULL},
#endif
    {"binbase", N_("gretl database directory"), NULL, paths.binbase, 
     'U', MAXLEN, 2, NULL},
    {"ratsbase", N_("RATS data directory"), NULL, paths.ratsbase, 
     'U', MAXLEN, 2, NULL},
    {"dbhost_ip", N_("Database server IP"), NULL, paths.dbhost_ip, 
     'U', 16, 2, NULL},
    {"dbproxy", N_("HTTP proxy (ipnumber:port)"), NULL, dbproxy, 
     'U', 21, 2, NULL},
    {"useproxy", N_("Use HTTP proxy"), NULL, &use_proxy, 
     'B', 1, 2, NULL},
    {"usecwd", N_("Use current working directory as default"), 
     N_("Use gretl user directory as default"), &usecwd, 'B', 0, 4, NULL},
    {"olddat", N_("Use \".dat\" as default datafile suffix"), 
     N_("Use \".gdt\" as default suffix"), &olddat, 'B', 0, 5, NULL},
    {"jwdata", N_("Toolbar folder icon opens Wooldridge data"), 
     N_("Toolbar folder icon opens Ramanathan data"), 
     &jwdata, 'B', 0, 5, NULL},
    {"fontspec", N_("Fixed font"), NULL, fontspec, 'U', MAXLEN, 0, NULL},
    {NULL, NULL, NULL, NULL, 0, 0, 0, NULL}   
};

/* ........................................................... */

void load_fixed_font (void)
{
    /* get a monospaced font for various windows */
    fixed_font = gdk_font_load(fontspec);
}

/* .................................................................. */

void get_default_dir (char *s)
{
    char *test = NULL;

    if (usecwd) {
	test = getcwd(s, MAXLEN);
	if (test == NULL) {
	    strcpy(s, paths.userdir);
	} else {
	    strcat(s, SLASHSTR);
	}
    } else {
	strcpy(s, paths.userdir); 
    }   
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
# endif

# ifdef HAVE_X12A
void set_x12a_ok (int set)
{
    static int ok;

    if (set >= 0) ok = set;
    if (mdata != NULL) {
	flip(mdata->ifac, "/Variable/X-12-ARIMA analysis", ok);
    }
}
# endif

static int check_for_prog (const char *prog)
{
    char tmp[MAXLEN];
    int ret;

    if (prog == NULL || *prog == 0) return 0;

    sprintf(tmp, "%s > /dev/null 2>&1", prog);
    ret = system(tmp);
    return (ret == 0);
}

static void set_tramo_x12a_dirs (void)
{
    char cmd[MAXLEN];

# ifdef HAVE_TRAMO 
    set_tramo_ok(check_for_prog(tramo));
    if (*tramodir == 0) {
	build_path(paths.userdir, "tramo", tramodir, NULL);
    }
# endif
# ifdef HAVE_X12A
    set_x12a_ok(check_for_prog(x12a));
    if (*x12adir == 0) {
	build_path(paths.userdir, "x12a", x12adir, NULL);
    }
# endif

    /* make directory structure */
# ifdef HAVE_X12A
    sprintf(cmd, "mkdir -p %s", x12adir);
    system(cmd);
# endif
# ifdef HAVE_TRAMO
    sprintf(cmd, "mkdir -p %s/output", tramodir);
    system(cmd);
    sprintf(cmd, "mkdir -p %s/graph/acf", tramodir);
    system(cmd);
    sprintf(cmd, "mkdir -p %s/graph/filters", tramodir);
    system(cmd);
    sprintf(cmd, "mkdir -p %s/graph/forecast", tramodir);
    system(cmd);
    sprintf(cmd, "mkdir -p %s/graph/series", tramodir);
    system(cmd);
    sprintf(cmd, "mkdir -p %s/graph/spectra", tramodir);
    system(cmd);
# endif /* HAVE_TRAMO */
}
#endif /* tramo or x12a */

/* ........................................................... */

#ifdef USE_GNOME
void set_rcfile (void)
{
    read_rc();
}
#else
void set_rcfile (void) 
{
    char *tmp;

    tmp = getenv("HOME");
    strcpy(rcfile, tmp);
    strcat(rcfile, "/.gretlrc");
    read_rc(); 
}
#endif

/* .................................................................. */

void options_dialog (gpointer data) 
{
    GtkWidget *tempwid, *dialog, *notebook;

    dialog = gtk_dialog_new ();
    gtk_window_set_title (GTK_WINDOW (dialog), _("gretl: options"));
    gtk_container_border_width 
	(GTK_CONTAINER (GTK_DIALOG (dialog)->vbox), 10);
    gtk_container_border_width 
	(GTK_CONTAINER (GTK_DIALOG (dialog)->action_area), 5);
    gtk_box_set_spacing (GTK_BOX (GTK_DIALOG (dialog)->vbox), 2);
    gtk_box_set_spacing (GTK_BOX (GTK_DIALOG (dialog)->action_area), 15);
    gtk_box_set_homogeneous (GTK_BOX (GTK_DIALOG (dialog)->action_area), TRUE);
    gtk_window_set_position (GTK_WINDOW (dialog), GTK_WIN_POS_MOUSE);
    gtk_signal_connect_object 
	(GTK_OBJECT (dialog), "delete_event", GTK_SIGNAL_FUNC 
	 (gtk_widget_destroy), GTK_OBJECT (dialog));
   
    notebook = gtk_notebook_new ();
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (dialog)->vbox), 
			notebook, TRUE, TRUE, 0);
    gtk_widget_show (notebook);

    make_prefs_tab (notebook, 1);
    make_prefs_tab (notebook, 2);
    make_prefs_tab (notebook, 3);
    make_prefs_tab (notebook, 4);
    make_prefs_tab (notebook, 5);
   
    tempwid = gtk_button_new_with_label ("OK");
    GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG 
				 (dialog)->action_area), 
			tempwid, TRUE, TRUE, 0);
    gtk_signal_connect (GTK_OBJECT (tempwid), "clicked", 
			GTK_SIGNAL_FUNC (apply_changes), NULL);
    gtk_signal_connect_object (GTK_OBJECT (tempwid), "clicked", 
			       GTK_SIGNAL_FUNC (gtk_widget_destroy), 
			       GTK_OBJECT (dialog));
    gtk_widget_show (tempwid);

    tempwid = gtk_button_new_with_label (_("  Cancel  "));
    GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG 
				 (dialog)->action_area), 
			tempwid, TRUE, TRUE, 0);
    gtk_signal_connect_object (GTK_OBJECT (tempwid), "clicked", 
			       GTK_SIGNAL_FUNC (gtk_widget_destroy), 
			       GTK_OBJECT (dialog));
    gtk_widget_show (tempwid);

    tempwid = gtk_button_new_with_label (_("Apply"));
    GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG 
				 (dialog)->action_area), 
			tempwid, TRUE, TRUE, 0);
    gtk_signal_connect (GTK_OBJECT (tempwid), "clicked", 
			GTK_SIGNAL_FUNC (apply_changes), NULL);
    gtk_widget_grab_default (tempwid);
    gtk_widget_show (tempwid);

    gtk_widget_show (dialog);
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

    strcpy(strvar, setting);

    while (rc_vars[i].key != NULL) {
	if (rc_vars[i].var == (void *) strvar) {
	    gtk_entry_set_text(GTK_ENTRY(rc_vars[i].widget), 
			       strvar);
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
    gtk_signal_connect(GTK_OBJECT(b), "clicked",
                       GTK_SIGNAL_FUNC(browse_button_callback), 
		       rc);
    return b;
}

/* .................................................................. */

static void make_prefs_tab (GtkWidget *notebook, int tab) 
{
    GtkWidget *box, *b_table, *s_table, *tempwid = NULL;
    int i, s_len, b_len, b_col;
    int s_count = 0, b_count = 0;
    RCVARS *rc = NULL;
   
    box = gtk_vbox_new (FALSE, 0);
    gtk_container_border_width (GTK_CONTAINER (box), 10);
    gtk_widget_show (box);

    if (tab == 1)
	tempwid = gtk_label_new (_("General"));
    else if (tab == 2)
	tempwid = gtk_label_new (_("Databases"));
    else if (tab == 3)
	tempwid = gtk_label_new (_("Programs"));
    else if (tab == 4)
	tempwid = gtk_label_new (_("Open/Save path"));
    else if (tab == 5)
	tempwid = gtk_label_new (_("Data files"));
    
    gtk_widget_show (tempwid);
    gtk_notebook_append_page (GTK_NOTEBOOK (notebook), box, tempwid);   

    s_len = 1;
    s_table = gtk_table_new (s_len, 2, FALSE);
    gtk_table_set_row_spacings (GTK_TABLE (s_table), 5);
    gtk_table_set_col_spacings (GTK_TABLE (s_table), 5);
    gtk_box_pack_start (GTK_BOX (box), s_table, FALSE, FALSE, 0);
    gtk_widget_show (s_table);
   
    b_len = b_col = 0;
    b_table = gtk_table_new (1, 2, FALSE);
    gtk_table_set_row_spacings (GTK_TABLE (b_table), 2);
    gtk_table_set_col_spacings (GTK_TABLE (b_table), 5);
    gtk_box_pack_start (GTK_BOX (box), b_table, FALSE, FALSE, 10);
    gtk_widget_show (b_table);

    i = 0;
    while (rc_vars[i].key != NULL) {
	rc = &rc_vars[i];
	if (rc->tab == tab) {
	    if (rc->type == 'B' && rc->link == NULL) { 
		/* simple boolean variable */
		b_count++;
		tempwid = gtk_check_button_new_with_label 
		    (_(rc->description));
		gtk_table_attach_defaults 
		    (GTK_TABLE (b_table), tempwid, b_col, b_col + 1, 
		     b_len, b_len + 1);
		if (*(int *)(rc->var))
		    gtk_toggle_button_set_active 
			(GTK_TOGGLE_BUTTON (tempwid), TRUE);
		else
		    gtk_toggle_button_set_active 
			(GTK_TOGGLE_BUTTON (tempwid), FALSE);
		/* special case: link between toggle and preceding entry */
		if (rc->len) {
		    gtk_widget_set_sensitive(rc_vars[i-1].widget,
					     GTK_TOGGLE_BUTTON(tempwid)->active);
		    gtk_signal_connect(GTK_OBJECT(tempwid), "clicked",
				       GTK_SIGNAL_FUNC(flip_sensitive),
				       rc_vars[i-1].widget);
		}
		/* end link to entry */
		gtk_widget_show (tempwid);
		rc->widget = tempwid;
		b_col++;
		if (b_col == 2) {
		    b_col = 0;
		    b_len++;
		    gtk_table_resize (GTK_TABLE (b_table), b_len + 1, 2);
		}
	    } else if (rc->type == 'B') { /* radio-button dichotomy */
		int val = *(int *)(rc->var);
		GSList *group;

		b_count++;
		b_len += 3;
		gtk_table_resize (GTK_TABLE(b_table), b_len + 1, 2);

		tempwid = gtk_radio_button_new_with_label(NULL, 
							  _(rc->description));
		gtk_table_attach_defaults 
		    (GTK_TABLE (b_table), tempwid, b_col, b_col + 1, 
		     b_len - 3, b_len - 2);    
		if (val) 
		    gtk_toggle_button_set_active 
			(GTK_TOGGLE_BUTTON(tempwid), TRUE);
		gtk_widget_show (tempwid);
		rc->widget = tempwid;
		group = gtk_radio_button_group(GTK_RADIO_BUTTON(tempwid));
		tempwid = gtk_radio_button_new_with_label(group, _(rc->link));
		gtk_table_attach_defaults 
		    (GTK_TABLE (b_table), tempwid, b_col, b_col + 1, 
		     b_len - 2, b_len - 1);  
		if (!val) {
		    gtk_toggle_button_set_active
			(GTK_TOGGLE_BUTTON(tempwid), TRUE);
		}
		gtk_widget_show (tempwid);
                tempwid = gtk_hseparator_new ();
                gtk_table_attach_defaults 
                    (GTK_TABLE (b_table), tempwid, b_col, b_col + 1, 
                     b_len - 1, b_len);  
                gtk_widget_show (tempwid);

	    } else { /* string variable */
		s_count++;
		s_len++;
		gtk_table_resize (GTK_TABLE (s_table), s_len, 
				  (tab == 3)? 3 : 2);
		tempwid = gtk_label_new (_(rc->description));
		gtk_misc_set_alignment (GTK_MISC (tempwid), 1, 0.5);
		gtk_table_attach_defaults (GTK_TABLE (s_table), 
					   tempwid, 0, 1, s_len-1, s_len);
		gtk_widget_show(tempwid);

		tempwid = gtk_entry_new();
		gtk_table_attach_defaults (GTK_TABLE (s_table), 
					   tempwid, 1, 2, s_len-1, s_len);
		gtk_entry_set_text(GTK_ENTRY(tempwid), rc->var);
		gtk_widget_show(tempwid);
		rc->widget = tempwid;
		/* program browse button */
		if (tab == 3 && strstr(rc->description, "directory") == NULL) {
		    tempwid = make_path_browse_button(rc);
		    gtk_table_attach_defaults(GTK_TABLE(s_table), 
					      tempwid, 2, 3, s_len-1, s_len);
		    gtk_widget_show(tempwid);
		}
	    } 
	}
	i++;
    }

    if (b_count == 0) gtk_widget_destroy(b_table);
    if (s_count == 0) gtk_widget_destroy(s_table);
}

/* .................................................................. */

#ifdef ENABLE_NLS
static void set_lcnumeric (void)
{
    if (lcnumeric) {
	putenv("LC_NUMERIC=");
	setlocale(LC_NUMERIC, "");
    } else {
	putenv("LC_NUMERIC=C");
	setlocale(LC_NUMERIC, "C");
    }
    reset_local_decpoint();
}
#endif

/* .................................................................. */

static void apply_changes (GtkWidget *widget, gpointer data) 
{
    const gchar *tempstr;
    extern void show_toolbar (void);
    int i = 0;
#ifdef ENABLE_NLS
    int lcnum_bak = lcnumeric;
#endif

    while (rc_vars[i].key != NULL) {
	if (rc_vars[i].widget != NULL) {
	    if (rc_vars[i].type == 'B') {
		if (GTK_TOGGLE_BUTTON(rc_vars[i].widget)->active)
		    *(int *)(rc_vars[i].var) = TRUE;
		else *(int *)(rc_vars[i].var) = FALSE;
	    } 
	    if (rc_vars[i].type == 'U' || rc_vars[i].type == 'R') {
		tempstr = gtk_entry_get_text
		    (GTK_ENTRY(rc_vars[i].widget));
		if (tempstr != NULL && strlen(tempstr)) 
		    strncpy(rc_vars[i].var, tempstr, rc_vars[i].len - 1);
	    }
	}
	i++;
    }
    write_rc();
    if (toolbar_box == NULL && want_toolbar)
	show_toolbar();
    else if (toolbar_box != NULL && !want_toolbar) {
	gtk_widget_destroy(toolbar_box);
	toolbar_box = NULL;
    }

#ifdef ENABLE_NLS
    set_lcnumeric();
    if (lcnumeric != lcnum_bak) 
	infobox(_("Please restart gretl to ensure consistent results"));
#endif

#if defined(HAVE_TRAMO) || defined (HAVE_X12A)
    set_tramo_x12a_dirs();
#endif

    proxy_init(dbproxy);
}

/* .................................................................. */

static void str_to_boolvar (char *s, void *b)
{
    if (strcmp(s, "true") == 0 || strcmp(s, "1") == 0)
	*(int *)b = TRUE;
    else
	*(int *)b = FALSE;	
}

/* .................................................................. */

static void boolvar_to_str (void *b, char *s)
{
    if (*(int *)b) strcpy(s, "true");
    else strcpy(s, "false");
}

/* .................................................................. */

#ifdef USE_GNOME

void write_rc (void) 
{
    char gpath[MAXSTR];
    char val[6];
    int i = 0;

    while (rc_vars[i].key != NULL) {
	sprintf(gpath, "/gretl/%s/%s", rc_vars[i].description, rc_vars[i].key);
	if (rc_vars[i].type == 'B') {
	    boolvar_to_str(rc_vars[i].var, val);
	    gnome_config_set_string(gpath, val);
	} else
	    gnome_config_set_string(gpath, rc_vars[i].var);
	i++;
    }
    printfilelist(1, NULL); /* data files */
    printfilelist(2, NULL); /* session files */
    printfilelist(3, NULL); /* script files */    
    gnome_config_sync();
    set_paths(&paths, 0, 1);
}

static void read_rc (void) 
{
    int i = 0;
    gchar *value = NULL;
    char gpath[MAXSTR];

    while (rc_vars[i].key != NULL) {
	sprintf(gpath, "/gretl/%s/%s", 
		rc_vars[i].description, 
		rc_vars[i].key);
	if ((value = gnome_config_get_string(gpath)) != NULL) {
	    if (rc_vars[i].type == 'B')
		str_to_boolvar(value, rc_vars[i].var);
	    else
		strncpy(rc_vars[i].var, value, rc_vars[i].len - 1);
	    g_free(value);
	}
	i++;
    }

    /* initialize lists of recently opened files */
    for (i=0; i<MAXRECENT; i++) { 
	datalist[i][0] = 0;
	sessionlist[i][0] = 0;
	scriptlist[i][0] = 0;
    }
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

    set_paths(&paths, 0, 1); /* 0 = not defaults, 1 = gui */

#if defined(HAVE_TRAMO) || defined(HAVE_X12A)
    set_tramo_x12a_dirs();
#endif

#ifdef ENABLE_NLS
    set_lcnumeric();
#endif
}

#else /* plain GTK */

void write_rc (void) 
{
    FILE *rc;
    int i;
    char val[6];

    rc = fopen(rcfile, "w");
    if (rc == NULL) {
	errbox(_("Couldn't open config file for writing"));
	return;
    }
    fprintf(rc, "# config file written by gretl: do not edit\n");
    i = 0;
    while (rc_vars[i].var != NULL) {
	fprintf(rc, "# %s\n", rc_vars[i].description);
	if (rc_vars[i].type == 'B') {
	    boolvar_to_str(rc_vars[i].var, val);
	    fprintf(rc, "%s = %s\n", rc_vars[i].key, val);
	} else
	    fprintf(rc, "%s = %s\n", rc_vars[i].key, (char *) rc_vars[i].var);
	i++;
    }
    printfilelist(1, rc); /* data files */
    printfilelist(2, rc); /* session files */
    printfilelist(3, rc); /* script files */
    fclose(rc);
    set_paths(&paths, 0, 1);
}

static void read_rc (void) 
{
    FILE *rc;
    int i, j;
    char line[MAXLEN], key[32], linevar[MAXLEN];
    int gotrecent = 0;

    if ((rc = fopen(rcfile, "r")) == NULL) return;

    i = 0;
    while (rc_vars[i].var != NULL) {
	if (fgets(line, MAXLEN, rc) == NULL) 
	    break;
	if (line[0] == '#') 
	    continue;
	if (!strncmp(line, "recent ", 7)) {
	    gotrecent = 1;
	    break;
	}
	if (sscanf(line, "%s", key) == 1) {
	    strcpy(linevar, line + strlen(key) + 3); 
	    chopstr(linevar); 
	    for (j=0; rc_vars[j].key != NULL; j++) {
		if (!strcmp(key, rc_vars[j].key)) {
		    if (rc_vars[j].type == 'B')
			str_to_boolvar(linevar, rc_vars[j].var);
		    else
			strcpy(rc_vars[j].var, linevar);
		}
	    }
	}
	i++;
    }

    /* get lists of recently opened files */
    for (i=0; i<MAXRECENT; i++) { 
	datalist[i][0] = 0;
	sessionlist[i][0] = 0;
	scriptlist[i][0] = 0;
    }
    if (gotrecent || (fgets(line, MAXLEN, rc) != NULL && 
	strncmp(line, "recent data files:", 18) == 0)) {
	i = 0;
	while (fgets(line, MAXLEN, rc) && i<MAXRECENT) {
	    if (strncmp(line, "recent session files:", 21) == 0)
		break;
	    chopstr(line);
	    if (strlen(line)) 
		strcpy(datalist[i++], line);
	}
    }
    if (strncmp(line, "recent session files:", 21) == 0) {
	i = 0;
	while (fgets(line, MAXLEN, rc) && i<MAXRECENT) {
	    if (strncmp(line, "recent script files:", 20) == 0)
		break;
	    chopstr(line);
	    if (strlen(line)) 
		strcpy(sessionlist[i++], line);
	}
    }
    if (strncmp(line, "recent script files:", 20) == 0) {
	i = 0;
	while (fgets(line, MAXLEN, rc) && i<MAXRECENT) {
	    chopstr(line);
	    if (strlen(line)) 
		strcpy(scriptlist[i++], line);
	}
    }
    fclose(rc);
    set_paths(&paths, 0, 1);

#if defined(HAVE_TRAMO) || defined(HAVE_X12A)
    set_tramo_x12a_dirs();
#endif

#ifdef ENABLE_NLS
    set_lcnumeric();
#endif
}

#endif /* end of "plain gtk" versions of read_rc, write_rc */

/* .................................................................. */

static void font_selection_ok (GtkWidget *w, GtkFontSelectionDialog *fs)
{
    gchar *fstring = gtk_font_selection_dialog_get_font_name(fs);

    if (strlen(fstring)) {
        strcpy(fontspec, fstring);
        gdk_font_unref(fixed_font);
        fixed_font = gdk_font_load(fontspec);
        write_rc();
    }
    g_free(fstring);
    gtk_widget_destroy(GTK_WIDGET (fs));
}

/* .................................................................. */

void font_selector (void)
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
	    (GTK_FONT_SELECTION_DIALOG (fontsel), fontspec);

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
/* .................................................................. */

void allocate_fileptrs (void)
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
    int i;
    const gchar *pathstart[] = {
	N_("/File/Open data"), 
	N_("/Session/Open"),
	N_("/File/Open command file")
    };

    for (i=0; i<MAXRECENT; i++) {
	sprintf(itempath, "%s/%d. %s", pathstart[filetype - 1],
		i+1, endbit(tmpname, filep[i], -1));
	w = gtk_item_factory_get_widget(mdata->ifac, itempath);
	if (w != NULL) 
	    gtk_item_factory_delete_item(mdata->ifac, itempath);
    }
}

/* .................................................................. */

void mkfilelist (int filetype, const char *fname)
{
    char *tmp[MAXRECENT-1];
    char **filep;
    int i, match = -1;

    if (filetype == 1) filep = datap;
    else if (filetype == 2) filep = sessionp;
    else if (filetype == 3) filep = scriptp;
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
    for (i=1; i<=match; i++) filep[i] = tmp[i-1];

    add_files_to_menu(filetype);
}

/* .................................................................. */

void delete_from_filelist (int filetype, const char *fname)
{
    char *tmp[MAXRECENT];
    char **filep;
    int i, match = -1;

    if (filetype == 1) filep = datap;
    else if (filetype == 2) filep = sessionp;
    else if (filetype == 3) filep = scriptp;
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

#if defined(USE_GNOME)

static void printfilelist (int filetype, FILE *fp)
     /* fp is ignored */
{
    int i;
    char **filep;
    char gpath[MAXLEN];
    static char *section[] = {"recent data files",
			      "recent session files",
			      "recent script files"};

    switch (filetype) {
    case 1: filep = datap; break;
    case 2: filep = sessionp; break;
    case 3: filep = scriptp; break;
    default: return;
    }

    for (i=0; i<MAXRECENT; i++) {
	sprintf(gpath, "/gretl/%s/%d", section[filetype - 1], i);
	gnome_config_set_string(gpath, filep[i]);
    }
}

#else /* "plain" version follows */

static void printfilelist (int filetype, FILE *fp)
{
    int i;
    char **filep;

    if (filetype == 1) {
	fprintf(fp, "recent data files:\n");
	filep = datap;
    } else if (filetype == 2) {
	fprintf(fp, "recent session files:\n");
	filep = sessionp;
    } else if (filetype == 3) {
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
    do_open_script(NULL, NULL);
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

    if (filetype == 1) {
	callfunc = set_data_from_filelist;
	filep = datap;
    } else if (filetype == 2) {
	callfunc = set_session_from_filelist;
	filep = sessionp;
    } else if (filetype == 3) {
	callfunc = set_script_from_filelist;
	filep = scriptp;
    }
    else return;

    /* See if there are any files to add */
    if (filep[0][0] == '\0') return;
    else {
	gchar *itemtype = "<Separator>";
	GtkWidget *w;

	/* is a separator already in place? */
	w = gtk_item_factory_get_widget(mdata->ifac, msep[filetype - 1]);
	if (w == NULL) {
	    fileitem.path = mymalloc(80);
	    strcpy(fileitem.path, msep[filetype - 1]);
	    fileitem.accelerator = NULL;
	    fileitem.callback = NULL;
	    fileitem.callback_action = 0;
	    fileitem.item_type = itemtype;
	    gtk_item_factory_create_item(mdata->ifac, &fileitem, NULL, 1);
	}
    }

    /* put the files under the menu separator */
    for (i=0; i<MAXRECENT; i++) {
	if (filep[i][0]) {
	    if (fileitem.path == NULL) fileitem.path = mymalloc(80);
	    fileitem.accelerator = NULL;
	    fileitem.callback_action = i; 
	    fileitem.item_type = NULL;
	    sprintf(fileitem.path, "%s/%d. %s", mpath[filetype - 1],
		    i+1, endbit(tmp, filep[i], 1));
	    fileitem.callback = callfunc; 
	    gtk_item_factory_create_item(mdata->ifac, &fileitem, NULL, 1);
	    w = gtk_item_factory_get_widget_by_action(mdata->ifac, i);
	    if (w != NULL)
		gtk_tooltips_set_tip(gretl_tips, w, filep[i], NULL);
	} else break;
    }
    free(fileitem.path);
}
