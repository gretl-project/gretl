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
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

/* fileselect.c for gretl -- use the gtkextra file selector under X11,
   the native MS file selector under MS Windows */

#include "gretl.h"

extern GtkItemFactoryEntry script_items[];
extern GtkItemFactoryEntry sample_script_items[];
extern void do_save_graph (const char *fname, const char *savestr);
extern int ps_print_plots (const char *fname, int flag, gpointer data);
extern int plot_to_xpm (const char *fname, gpointer data);
#ifdef GNUPLOT_PNG
extern void save_this_graph (GPT_SPEC *plot, const char *fname);
#endif

extern int olddat; /* gui_utils.c */

static char remember_dir[MAXLEN];

struct extmap {
    int action;
    char *ext;
};

static struct extmap action_map[] = {
    {SAVE_DATA, ".gdt"},
    {SAVE_GZDATA, ".gdt"},
    {SAVE_BIN1, ".gdt"},
    {SAVE_BIN2, ".gdt"},
    {SAVE_CMDS, ".inp"},
    {SAVE_SCRIPT, ".inp"},
    {SAVE_CONSOLE, ".inp"},
    {SAVE_MODEL, ".txt"},
    {SAVE_SESSION, ".gretl"},
    {SAVE_GP_CMDS, ".gp"},
    {SAVE_BOXPLOT_EPS, ".eps"},
    {SAVE_BOXPLOT_PS, ".ps"},
    {SAVE_BOXPLOT_XPM, ".xpm"},
    {EXPORT_CSV, ".csv"},
    {EXPORT_R, ".R"},
    {EXPORT_R_ALT, ".R"},
    {EXPORT_OCTAVE, ".m"},
    {SAVE_OUTPUT, ".txt"},
    {SAVE_TEX_TAB, ".tex"},
    {SAVE_TEX_EQ, ".tex"},
    {OPEN_DATA, ".gdt"},
    {OPEN_SCRIPT, ".inp"},
    {OPEN_SESSION, ".gretl"},
    {OPEN_CSV,  ".csv"},
    {APPEND_CSV,  ".csv"},
    {OPEN_BOX, ".box"},
    {OPEN_GNUMERIC, ".gnumeric"},
    {APPEND_GNUMERIC, ".gnumeric"},
    {OPEN_EXCEL, ".xls"},
    {APPEND_EXCEL, ".xls"},
    {OP_MAX, NULL}
};

/* ........................................................... */

static int action_to_flag (const int action)
{
    switch (action) {
    case SAVE_GZDATA: return OPT_Z;
    case SAVE_BIN1: return OPT_S;
    case SAVE_BIN2: return OPT_O;
    case EXPORT_OCTAVE: return OPT_M;
    case EXPORT_R: return OPT_R;
    case EXPORT_R_ALT: return OPT_R_ALT;
    case EXPORT_CSV: return OPT_C;
    default: return 0;
    }
}

/* ........................................................... */

static const char *get_gp_ext (const char *termtype)
{
    if (!strcmp(termtype, "postscript")) return ".eps";
    else if (!strcmp(termtype, "fig")) return ".fig";
    else if (!strcmp(termtype, "latex")) return ".tex";
    else if (!strcmp(termtype, "png")) return ".png";
    else if (!strcmp(termtype, "plot commands")) return ".gp";
    else return "*";
}

/* ........................................................... */

static int is_data_action (int i)
{
    if (i == SAVE_DATA || i == SAVE_GZDATA || i == SAVE_BIN1 || 
	i == SAVE_BIN2 || i == OPEN_DATA)
	return 1;
    else
	return 0;
}

/* ........................................................... */

static int dat_ext (char *str, int err)
{
    char *suff;

    if (str == NULL) return 0;
    suff = strrchr(str, '.');
    if (suff != NULL && 
	(!strcmp(suff, ".dat") || !strcmp(suff, ".gdt"))) {
	if (err) 
	    errbox(_("The suffix you selected should be used\n"
		   "only for gretl datafiles"));
	return 1;
    }
    return 0;
}

/* ........................................................... */

static const char *get_ext (int action, gpointer data)
{
    const char *s = NULL;

    if (olddat && is_data_action(action)) 
	return ".dat";

    if (action == SAVE_GNUPLOT || action == SAVE_THIS_GRAPH) {
	GPT_SPEC *plot = (GPT_SPEC *) data;
	s = get_gp_ext(plot->termtype);
    }
    else if (action == SAVE_LAST_GRAPH) 
	s = get_gp_ext(data);
    else {
	int i;

	for (i=0; i < sizeof action_map / sizeof *action_map; i++) {
	    if (action == action_map[i].action) {
		s = action_map[i].ext;
		break;
	    }
	}
    }
    return s;
}

/* ........................................................... */

static void maybe_add_ext (char *fname, int action, gpointer data)
{
    FILE *fp;
    const char *ext = NULL;

    /* don't mess with the name of a previously existing file */
    if ((fp = fopen(fname, "r")) && fgetc(fp) != EOF) {
	fclose(fp);
	return;
    }    

    /* don't mess with a filename that already has an extension */
    if (dotpos(fname) != strlen(fname)) return;
    
    /* otherwise add an appropriate extension */
    ext = get_ext(action, data);
    if (ext != NULL && strlen(ext) > 1) 
	strcat(fname, ext);
}

/* ........................................................... */

static void script_set_title (windata_t *vwin, const char *fname)
{
    gchar *title;
    const char *p = strrchr(fname, SLASH);

    title = g_strdup_printf("gretl: %s", p? p + 1 : fname);
    gtk_window_set_title(GTK_WINDOW(vwin->dialog), title);
    strcpy(vwin->fname, fname);
    g_free(title);
}

/* ........................................................... */

static void save_editable_content (int action, const char *fname,
				   windata_t *vwin)
{
    FILE *fp;
    gchar *buf;

    if ((fp = fopen(fname, "w")) == NULL) {
	errbox(_("Couldn't open file for writing"));
	return;
    }

    buf = gtk_editable_get_chars(GTK_EDITABLE(vwin->w), 0, -1);
    if (buf == NULL) {
	errbox("Couldn't retrieve buffer");
	fclose(fp);
	return;
    }

    fprintf(fp, "%s", buf);
    g_free(buf);
    infobox(_("File saved OK"));
    fclose(fp);
    
    if (action == SAVE_SCRIPT) {
	strcpy(scriptfile, fname);
	mkfilelist(3, scriptfile);
	vwin->active_var = 0; /* zero out "changed" flag */
	script_set_title(vwin, fname);
    }
}

/* ........................................................... */

          /* MS Windows version of file selection code */

/* ........................................................... */

#ifdef G_OS_WIN32 

#include <windows.h>

struct winfilter {
    const char *descrip;
    const char *pat;
} winfilter;

struct win32_filtermap {
    int action;
    struct winfilter filter;
};

static struct winfilter get_gp_filter (const char *termtype)
{
    static struct winfilter gpfilters[] = {
	{ N_("postscript files (*.eps)"), "*.eps" },
	{ N_("xfig files (*.fig)"), "*.fig" },
	{ N_("LaTeX files (*.tex)"), "*.tex" },
	{ N_("PNG files (*.png)"), "*.png" },
	{ N_("gnuplot files (*.gp)"), "*.gp" },
	{ N_("all files (*.*)"), "*.*" }
    };

    if (!strcmp(termtype, "postscript")) 
	return gpfilters[0];
    else if (!strcmp(termtype, "fig")) 
	return gpfilters[1];
    else if (!strcmp(termtype, "latex")) 
	return gpfilters[2];
    else if (!strcmp(termtype, "png")) 
	return gpfilters[3];
    else if (!strcmp(termtype, "plot commands")) 
	return gpfilters[4];
    else 
	return gpfilters[5];
}

static struct winfilter get_filter (int action, gpointer data)
{
    int i;
    struct winfilter filter;
    static struct win32_filtermap map[] = {
	{SAVE_DATA, { N_("gretl data files (*.gdt)"), "*.gdt" }},
	{SAVE_GZDATA, { N_("gretl data files (*.gdt)"), "*.gdt" }},
	{SAVE_BIN1, { N_("gretl data files (*.gdt)"), "*.gdt" }},
	{SAVE_BIN2, { N_("gretl data files (*.gdt)"), "*.gdt" }},
	{SAVE_CMDS, { N_("gretl command files (*.inp)"), "*.inp" }},
	{SAVE_SCRIPT, { N_("gretl script files (*.inp)"), "*.inp" }},
	{SAVE_CONSOLE, { N_("gretl command files (*.inp)"), "*.inp" }},
	{SAVE_MODEL, { N_("text files (*.txt)"), "*.txt" }},
	{SAVE_SESSION, { N_("session files (*.gretl)"), "*.gretl" }},
	{SAVE_BOXPLOT_EPS, { N_("postscript files (*.eps)"), "*.eps" }},
	{SAVE_BOXPLOT_PS, { N_("postscript files (*.ps)"), "*.ps" }},
	{SAVE_GP_CMDS, { N_("gnuplot files (*.gp)"), "*.gp" }},
	{EXPORT_CSV, { N_("CSV files (*.csv)"), "*.csv" }},
	{EXPORT_R, { N_("GNU R files (*.R)"), "*.R" }},
	{EXPORT_R_ALT, { N_("GNU R files (*.R)"), "*.R" }},
	{EXPORT_OCTAVE, { N_("GNU Octave files (*.m)"), "*.m" }},
	{SAVE_OUTPUT, { N_("text files (*.txt)"), "*.txt" }},
	{SAVE_TEX_TAB, { N_("TeX files (*.tex)"), "*.tex" }},
	{SAVE_TEX_EQ, { N_("TeX files (*.tex)"), "*.tex" }},
	{OPEN_DATA, { N_("gretl data files (*.gdt)"), "*.gdt" }},
	{OPEN_SCRIPT, { N_("gretl script files (*.inp)"), "*.inp" }},
	{OPEN_SESSION, { N_("session files (*.gretl)"), "*.gretl" }},
	{OPEN_CSV, { N_("CSV files (*.csv)"), "*.csv" }},
	{APPEND_CSV, { N_("CSV files (*.csv)"), "*.csv" }},
	{OPEN_BOX, { N_("BOX data files (*.box)"), "*.box" }},
	{OPEN_GNUMERIC, { N_("Gnumeric files (*.gnumeric)"), "*.gnumeric" }},
	{APPEND_GNUMERIC, { N_("Gnumeric files (*.gnumeric)"), "*.gnumeric" }},
	{OPEN_EXCEL, { N_("Excel files (*.xls)"), "*.xls" }},
	{APPEND_EXCEL, { N_("Excel files (*.xls)"), "*.xls" }},
    };
    static struct winfilter olddat_filter = {
	N_("gretl data files (*.dat)"), "*.dat"
    };

    if (olddat && is_data_action(action)) 
	return olddat_filter;

    if (action == SAVE_GNUPLOT || action == SAVE_THIS_GRAPH) {
	GPT_SPEC *plot = (GPT_SPEC *) data;
	return get_gp_filter(plot->termtype);
    }
    if (action == SAVE_LAST_GRAPH) {
	return get_gp_filter(data);
    }
    for (i=0; i< sizeof map/sizeof *map; i++) {
	if (action == map[i].action) {
	    filter = map[i].filter;
	    break;
	}
    }
    return filter;
}

static char *make_winfilter (int action, gpointer data)
{
    char *p = mymalloc(128);
    char *start = p;
    gchar *trf = NULL;
    struct winfilter filter;

    if (p == NULL) return NULL;

    filter = get_filter(action, data);

    if (nls_on) {
	gint wrote;

	trf = g_locale_from_utf8 (_(filter.descrip), -1, NULL, &wrote, NULL);
	strcpy(p, trf);
    } else strcpy(p, _(filter.descrip));

    p += strlen(p) + 1;
    strcpy(p, filter.pat);
    p += strlen(p) + 1;
    strcpy(p, _("all files (*.*)")); /* ENCODING? */
    p += strlen(p) + 1;
    strcpy(p, "*.*");
    p += strlen(p) + 1;
    *p = '\0';

    if (nls_on) g_free(trf);

    return start;
}

/* ........................................................... */

void file_selector (char *msg, int action, gpointer data) 
{
    OPENFILENAME of;
    int retval;
    char fname[MAXLEN], endname[64], startd[MAXLEN];
    char *filter;
    gchar *trmsg;

    fname[0] = '\0';
    endname[0] = '\0';
    if (remember_dir[0] != '\0')
	strcpy(startd, remember_dir);
    else
	get_default_dir(startd);

    /* special case: default save of data */
    if ((action == SAVE_DATA || action == SAVE_GZDATA) && paths.datfile[0]
	&& !strcmp(paths.datfile + strlen(paths.datfile) - 4, 
		(olddat)? ".dat" : ".gdt")) {
	strcpy(fname, paths.datfile + slashpos(paths.datfile) + 1);
	get_base(startd, paths.datfile, SLASH);
    }

    if (nls_on) {
	gint wrote;

	trmsg = g_locale_from_utf8 (msg, -1, NULL, &wrote, NULL);
    } else trmsg = msg;

    /* initialize file dialog info struct */
    memset(&of, 0, sizeof of);
#ifdef OPENFILENAME_SIZE_VERSION_400
    of.lStructSize = OPENFILENAME_SIZE_VERSION_400;
#else
    of.lStructSize = sizeof of;
#endif
    of.hwndOwner = NULL;
    filter = make_winfilter(action, data);
    of.lpstrFilter = filter;
    of.lpstrCustomFilter = NULL;
    of.nFilterIndex = 1;
    of.lpstrFile = fname;
    of.nMaxFile = sizeof fname;
    of.lpstrFileTitle = endname;
    of.nMaxFileTitle = sizeof endname;
    of.lpstrInitialDir = startd;
    of.lpstrTitle = trmsg;
    of.lpstrDefExt = NULL;
    of.Flags = OFN_HIDEREADONLY;

    if (action < END_OPEN)
	retval = GetOpenFileName(&of);
    else  /* a file save action */
	retval = GetSaveFileName(&of);

    free(filter);
    if (nls_on) g_free(trmsg);

    if (!retval) {
	if (CommDlgExtendedError())
	    errbox(_("File dialog box error"));
	return;
    }
	
    strncpy(remember_dir, fname, slashpos(fname));

    if (action == OPEN_DATA || action == OPEN_CSV || 
	action == OPEN_BOX || action == OPEN_GNUMERIC || action == OPEN_EXCEL) {
	strcpy(trydatfile, fname);
	verify_open_data(NULL, action);
    }
    else if (action == APPEND_CSV || action == APPEND_GNUMERIC || 
	     action == APPEND_EXCEL) {
	strcpy(trydatfile, fname);
	do_open_data(NULL, NULL, action);
    }
    else if (action == OPEN_SCRIPT) {
	int spos;

	strcpy(tryscript, fname);

	if (view_file(tryscript, 1, 0, 78, 370, EDIT_SCRIPT, 
		      script_items) != NULL) {
	    strcpy(scriptfile, tryscript);
	    mkfilelist(3, scriptfile);
	    spos = slashpos(scriptfile);
	    if (spos) strncpy(paths.currdir, scriptfile, spos + 1);
	}
    }
    else if (action == OPEN_SESSION) {
	int pub = !strncmp(tryscript, paths.scriptdir, strlen(paths.scriptdir));

	strcpy(tryscript, fname);

	if (saved_objects(tryscript)) {
	    verify_open_session(NULL);
	    return;
	} 
	if (view_file(tryscript, 1, 0, 78, 370, 
		      pub ? VIEW_SCRIPT : EDIT_SCRIPT, 
		      pub ? sample_script_items : script_items))
		strcpy(scriptfile, tryscript);
    }

    if (action < END_OPEN) return;

    /* now for the save options */

    if (action > SAVE_BIN2 && dat_ext(fname, 1)) return;

    maybe_add_ext(fname, action, data);

    if (action >= SAVE_DATA && action < END_SAVE_DATA) {
	int overwrite = 0;

	if (!strcmp(fname, paths.datfile)) overwrite = 1;
	do_store(fname, action_to_flag(action), overwrite);
    }
    else if (action == SAVE_GNUPLOT) {
	int err = 0;
	GPT_SPEC *plot = (GPT_SPEC *) data;

	err = go_gnuplot(plot, fname, &paths);
	if (err == 0) infobox(_("graph saved"));
	else if (err == 1) errbox(_("gnuplot command failed"));
	else if (err == 2) infobox(_("There were missing observations"));
    }
#ifdef GNUPLOT_PNG
    else if (action == SAVE_THIS_GRAPH) {
	GPT_SPEC *plot = (GPT_SPEC *) data;

	save_this_graph(plot, fname);
    }
#endif
    else if (action == SAVE_BOXPLOT_EPS || action == SAVE_BOXPLOT_PS) {
	int err;

	err = ps_print_plots(fname, action, data);
	if (!err) infobox(_("boxplots saved"));
	else errbox(_("boxplot save failed"));
    }
    else if (action == SAVE_LAST_GRAPH) {
	char *savestr = (char *) data;
	
	do_save_graph(fname, savestr);
    }    
    else if (action == SAVE_SESSION) {
	save_session(fname);
    }
    else if (action == SAVE_TEX_TAB || action == SAVE_TEX_EQ) {
	MODEL *pmod = (MODEL *) data;
	do_save_tex(fname, action, pmod); 
    }
    else {
	windata_t *vwin = (windata_t *) data;

	save_editable_content(action, fname, vwin);

    }
}

#else /* End of MS Windows file selection code */

/* ........................................................... */

static void filesel_callback (GtkWidget *w, gpointer data) 
{
    GtkIconFileSel *fs = GTK_ICON_FILESEL(data);
    gint action = GPOINTER_TO_INT(gtk_object_get_data
				  (GTK_OBJECT(data), "action"));
    char fname[MAXLEN];
    char *test, *path;
    FILE *fp = NULL;
    gpointer extdata = NULL;

    test = gtk_entry_get_text(GTK_ENTRY(fs->file_entry));
    if (test == NULL || *test == '\0') 
	return;    
    path = gtk_file_list_get_path(GTK_FILE_LIST(fs->file_list));
    sprintf(fname, "%s%s", path, test);

    /* do some elementary checking */
    if (action < END_OPEN) {
	if ((fp = fopen(fname, "r")) == NULL) {
	    errbox(_("Couldn't open the specified file"));
	    return;
	} else fclose(fp);
    } 

    strcpy(remember_dir, path);

    if (action == OPEN_DATA || action == OPEN_CSV || 
	action == OPEN_BOX || action == OPEN_GNUMERIC || action == OPEN_EXCEL) {
	strcpy(trydatfile, fname);
	gtk_widget_destroy(GTK_WIDGET(fs));  
	verify_open_data(NULL, action);
	return;
    }
    else if (action == APPEND_CSV || action == APPEND_GNUMERIC || 
	     action == APPEND_EXCEL) {
	strcpy(trydatfile, fname);
	gtk_widget_destroy(GTK_WIDGET(fs)); 
	do_open_data(NULL, NULL, action);
	return;
    }
    else if (action == OPEN_SCRIPT) {
	int spos;

	strcpy(tryscript, fname);

	if (view_file(tryscript, 1, 0, 78, 370, 
		      EDIT_SCRIPT, script_items) != NULL) {
	    strcpy(scriptfile, tryscript);
	    mkfilelist(3, scriptfile);
	    spos = slashpos(scriptfile);
	    if (spos) strncpy(paths.currdir, scriptfile, spos + 1);
	}

    }
    else if (action == OPEN_SESSION) {
	int pub = !strncmp(tryscript, paths.scriptdir, strlen(paths.scriptdir));

	strcpy(tryscript, fname);

	if (saved_objects(tryscript)) {
	    verify_open_session(NULL);
	    gtk_widget_destroy(GTK_WIDGET(fs));    
	    return;
	} 
	if (view_file(tryscript, 1, 0, 78, 370, 
		      pub ? VIEW_SCRIPT : EDIT_SCRIPT, 
		      pub ? sample_script_items : script_items))
		strcpy(scriptfile, tryscript);
    }

    if (action < END_OPEN) {
	gtk_widget_destroy(GTK_WIDGET(fs));    
	return;
    }

    /* now for the save options */

    if (action == SAVE_GNUPLOT || action == SAVE_LAST_GRAPH || 
	action == SAVE_THIS_GRAPH) 
	extdata = gtk_object_get_data(GTK_OBJECT(fs), "graph");

    maybe_add_ext(fname, action, extdata); 

    if (action > SAVE_BIN2 && dat_ext(fname, 1)) {
	gtk_widget_destroy(GTK_WIDGET(fs));
	return;
    }

    if (action >= SAVE_DATA && action < END_SAVE_DATA) {
	int overwrite = 0;

	if (!strcmp(fname, paths.datfile)) overwrite = 1;
	do_store(fname, action_to_flag(action), overwrite);
    }
    else if (action == SAVE_GNUPLOT) {
	int err = 0;
	GPT_SPEC *plot = gtk_object_get_data(GTK_OBJECT(fs), "graph");

	err = go_gnuplot(plot, fname, &paths);
	if (err == 0) infobox(_("graph saved"));
	else if (err == 1) errbox(_("gnuplot command failed"));
	else if (err == 2) infobox(_("There were missing observations"));
    }
#ifdef GNUPLOT_PNG
    else if (action == SAVE_THIS_GRAPH) {
	GPT_SPEC *plot = gtk_object_get_data(GTK_OBJECT(fs), "graph");

	save_this_graph(plot, fname);
    }
#endif
    else if (action == SAVE_BOXPLOT_EPS || action == SAVE_BOXPLOT_PS) {
	int err;

	err = ps_print_plots(fname, action,
			     gtk_object_get_data(GTK_OBJECT(fs), "graph"));
	if (!err) infobox(_("boxplots saved"));
	else errbox(_("boxplot save failed"));
    }
    else if (action == SAVE_BOXPLOT_XPM) {
	int err;

	err = plot_to_xpm(fname, gtk_object_get_data(GTK_OBJECT(fs), "graph"));
	if (!err) infobox(_("boxplots saved"));
	else errbox(_("boxplot save failed"));
    }
    else if (action == SAVE_LAST_GRAPH) {
	char *savestr = gtk_object_get_data(GTK_OBJECT(fs), "graph");
	
	do_save_graph(fname, savestr);
    }    
    else if (action == SAVE_SESSION) {
	save_session(fname);
    }
    else if (action == SAVE_TEX_TAB || action == SAVE_TEX_EQ) {
	MODEL *pmod;
	pmod = (MODEL *) gtk_object_get_data(GTK_OBJECT(fs), "model");
	do_save_tex(fname, action, pmod); 
    }
    else {
	windata_t *vwin = gtk_object_get_data(GTK_OBJECT(fs), "text");

	save_editable_content(action, fname, vwin);

    }

    gtk_widget_destroy(GTK_WIDGET(fs));    
}

/* ........................................................... */

static void extra_get_filter (int action, gpointer data, char *suffix)
{
    
    const char *ext = get_ext(action, data);

    if (ext == NULL) 
	strcpy(suffix, "*");
    else
	sprintf(suffix, "*%s", ext);
}

/* ........................................................... */

void file_selector (char *msg, int action, gpointer data) 
{
    GtkWidget *filesel;
    int gotdir = 0;
    char suffix[8], startdir[MAXLEN];

    if (remember_dir[0] != '\0')
	strcpy(startdir, remember_dir);
    else
	get_default_dir(startdir);

    filesel = gtk_icon_file_selection_new(msg);

    if (strstr(startdir, "/."))
	gtk_icon_file_selection_show_hidden(GTK_ICON_FILESEL(filesel), TRUE);
    else 
	gtk_icon_file_selection_show_hidden(GTK_ICON_FILESEL(filesel), FALSE);

    gtk_object_set_data(GTK_OBJECT(filesel), "action", GINT_TO_POINTER(action));

    extra_get_filter(action, data, suffix);
    gtk_icon_file_selection_set_filter(GTK_ICON_FILESEL(filesel), suffix);

    gtk_signal_connect(GTK_OBJECT(GTK_ICON_FILESEL(filesel)->ok_button),
		       "clicked", 
		       GTK_SIGNAL_FUNC(filesel_callback), filesel);

    if (action > END_OPEN) /* a file save action */
	gtk_object_set_data(GTK_OBJECT(filesel), "text", data);

    /* special cases */

    if (action == SAVE_GNUPLOT || action == SAVE_THIS_GRAPH  
	|| action == SAVE_LAST_GRAPH ||
	action == SAVE_BOXPLOT_EPS || action == SAVE_BOXPLOT_PS ||
	action == SAVE_BOXPLOT_XPM) 
	gtk_object_set_data(GTK_OBJECT(filesel), "graph", data);

    else if (action == SAVE_TEX_TAB || action == SAVE_TEX_EQ) 
	gtk_object_set_data(GTK_OBJECT(filesel), "model", data);

    else if ((action == SAVE_DATA || action == SAVE_GZDATA) 
	     && paths.datfile[0] && dat_ext(paths.datfile, 0)) {
	char *fname = paths.datfile + slashpos(paths.datfile) + 1;
	char startd[MAXLEN];

	gtk_entry_set_text(GTK_ENTRY(GTK_ICON_FILESEL(filesel)->file_entry),
			   fname);
	if (get_base(startd, paths.datfile, SLASH) == 1) {
	    gtk_icon_file_selection_open_dir(GTK_ICON_FILESEL(filesel), startd);
	    gotdir = 1;
	}
    }

    if (!gotdir)
	gtk_icon_file_selection_open_dir(GTK_ICON_FILESEL(filesel), startdir);

    gtk_signal_connect(GTK_OBJECT(GTK_ICON_FILESEL(filesel)), "destroy",
		       gtk_main_quit, NULL);
    gtk_signal_connect_object(GTK_OBJECT(GTK_ICON_FILESEL
					 (filesel)->cancel_button),
			      "clicked", (GtkSignalFunc) gtk_widget_destroy,
			      GTK_OBJECT (filesel));

    gtk_widget_show(filesel);
    gtk_main(); /* make file selector modal */
}

#endif /* end of non-MS Windows code */

