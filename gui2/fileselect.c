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
#include "fileselect.h"
#include "boxplots.h"
#include "gpt_control.h"

extern GtkItemFactoryEntry script_items[];
extern GtkItemFactoryEntry sample_script_items[];

extern int olddat; /* settings.c */

static char remember_dir[MAXLEN];

struct extmap {
    int action;
    char *ext;
};

static struct extmap action_map[] = {
    {SAVE_DATA, ".gdt"},
    {SAVE_DATA_AS, ".gdt"},
    {SAVE_GZDATA, ".gdt"},
    {SAVE_BIN1, ".gdt"},
    {SAVE_BIN2, ".gdt"},
    {SAVE_CMDS, ".inp"},
    {SAVE_SCRIPT, ".inp"},
    {SAVE_CONSOLE, ".inp"},
    {SAVE_MODEL, ".txt"},
    {SAVE_SESSION, ".gretl"},
    {SAVE_GP_CMDS, ".plt"},
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
    {SAVE_TEX_TAB_FRAG, ".tex"},
    {SAVE_TEX_EQ_FRAG, ".tex"},
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
    {OPEN_DES, ".des"},
    {FILE_OP_MAX, NULL}
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
    else if (!strcmp(termtype, "plot commands")) return ".plt";
    else return "*";
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

    if (olddat && IS_DAT_ACTION(action)) { 
	return ".dat";
    }

    if (action == SAVE_GNUPLOT || action == SAVE_THIS_GRAPH) {
	GPT_SPEC *plot = (GPT_SPEC *) data;
	s = get_gp_ext(plot->termtype);
    }
    else if (action == SAVE_LAST_GRAPH) {
	s = get_gp_ext(data);
    } else {
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

static int check_maybe_add_ext (char *fname, int action, gpointer data)
{
    FILE *fp;
    const char *ext = NULL;

    if (fname == NULL) return 1;

    /* don't mess if the fname is really a dir */
    if (isdir(fname)) return 1;

    /* don't mess with the name of a previously existing file */
    if ((fp = fopen(fname, "r")) && fgetc(fp) != EOF) {
	fclose(fp);
	return 0;
    }    

    /* don't mess with a filename that already has an extension */
    if (dotpos(fname) != strlen(fname)) return 0;
    
    /* otherwise add an appropriate extension */
    ext = get_ext(action, data);
    if (ext != NULL && strlen(ext) > 1) {
	strcat(fname, ext);
    }

    return 0;
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
#ifdef ENABLE_NLS
    gsize bytes;
    gchar *trbuf;
#endif

    buf = textview_get_text(GTK_TEXT_VIEW(vwin->w));

    if (buf == NULL) {
	errbox("Couldn't retrieve buffer");
	return;
    }

    if ((fp = fopen(fname, "w")) == NULL) {
	errbox(_("Couldn't open file for writing"));
	g_free(buf);
	return;
    }

#ifdef ENABLE_NLS
    trbuf = g_locale_from_utf8(buf, -1, NULL, &bytes, NULL);
    fprintf(fp, "%s", trbuf);
    g_free(trbuf);
#else
    fprintf(fp, "%s", buf);
#endif

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

static void set_startdir (int action, char *startdir)
{
    if (*remember_dir != '\0')
	strcpy(startdir, remember_dir);
    else
	get_default_dir(startdir);

#ifndef G_OS_WIN32
    if (startdir[strlen(startdir)-1] != '/') strcat(startdir, "/");
#endif
}

/* ........................................................... */

static void filesel_open_script (const char *fname)
{
    int spos;

    strcpy(tryscript, fname);

    if (view_file(tryscript, 1, 0, 78, 370, EDIT_SCRIPT, 
		  script_items) != NULL) {
	strcpy(scriptfile, tryscript);
	mkfilelist(3, scriptfile);
	spos = slashpos(scriptfile);
	if (spos) {
	    paths.currdir[0] = 0;
	    strncat(paths.currdir, scriptfile, spos + 1);
	}
    }
}

/* ........................................................... */

static void filesel_open_session (const char *fname)
{
    int pub = !strncmp(tryscript, paths.scriptdir, strlen(paths.scriptdir));

    strcpy(tryscript, fname);

    if (saved_objects(tryscript)) {
	verify_open_session(NULL);
    } else if (view_file(tryscript, 1, 0, 78, 370, 
			 pub ? VIEW_SCRIPT : EDIT_SCRIPT, 
			 pub ? sample_script_items : script_items)) {
	strcpy(scriptfile, tryscript);
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
	{ N_("gnuplot files (*.plt)"), "*.plt" },
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
	{SAVE_DATA_AS, { N_("gretl data files (*.gdt)"), "*.gdt" }},
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
	{SAVE_GP_CMDS, { N_("gnuplot files (*.plt)"), "*.plt" }},
	{EXPORT_CSV, { N_("CSV files (*.csv)"), "*.csv" }},
	{EXPORT_R, { N_("GNU R files (*.R)"), "*.R" }},
	{EXPORT_R_ALT, { N_("GNU R files (*.R)"), "*.R" }},
	{EXPORT_OCTAVE, { N_("GNU Octave files (*.m)"), "*.m" }},
	{SAVE_OUTPUT, { N_("text files (*.txt)"), "*.txt" }},
	{SAVE_TEX_TAB, { N_("TeX files (*.tex)"), "*.tex" }},
	{SAVE_TEX_EQ, { N_("TeX files (*.tex)"), "*.tex" }},
	{SAVE_TEX_TAB_FRAG, { N_("TeX files (*.tex)"), "*.tex" }},
	{SAVE_TEX_EQ_FRAG, { N_("TeX files (*.tex)"), "*.tex" }},
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
	{OPEN_DES, { N_("DES files (*.des)"), "*.des" }}
    };
    static struct winfilter olddat_filter = {
	N_("gretl data files (*.dat)"), "*.dat"
    };

    if (olddat && IS_DAT_ACTION(action)) {
	return olddat_filter;
    }

    if (action == SAVE_GNUPLOT || action == SAVE_THIS_GRAPH) {
	GPT_SPEC *plot = (GPT_SPEC *) data;

	return get_gp_filter(plot->termtype);
    }
    if (action == SAVE_LAST_GRAPH) {
	return get_gp_filter(data);
    }
    for (i=0; i< sizeof map / sizeof *map; i++) {
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
    struct winfilter filter;

    if (p == NULL) return NULL;

    filter = get_filter(action, data);

    strcpy(p, I_(filter.descrip));

    p += strlen(p) + 1;
    strcpy(p, filter.pat);
    p += strlen(p) + 1;
    strcpy(p, I_("all files (*.*)"));
    p += strlen(p) + 1;
    strcpy(p, "*.*");
    p += strlen(p) + 1;
    *p = '\0';

    return start;
}

/* ........................................................... */

void file_selector (char *msg, int action, gpointer data) 
{
    OPENFILENAME of;
    int retval;
    char fname[MAXLEN], endname[64], startdir[MAXLEN];
    char *filter;
    gchar *trmsg;

    fname[0] = '\0';
    endname[0] = '\0';

    set_startdir(action, startdir);

    /* special case: default save of data */
    if ((action == SAVE_DATA || action == SAVE_GZDATA) && paths.datfile[0]
	&& !strcmp(paths.datfile + strlen(paths.datfile) - 4, 
		(olddat)? ".dat" : ".gdt")) {
	strcpy(fname, paths.datfile + slashpos(paths.datfile) + 1);
	get_base(startdir, paths.datfile, SLASH);
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
    of.lpstrInitialDir = startdir;
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

    if (OPEN_DATA_ACTION(action)) {
	strcpy(trydatfile, fname);
	verify_open_data(NULL, action);
    }
    else if (APPEND_DATA_ACTION(action)) {
	strcpy(trydatfile, fname);
	do_open_data(NULL, NULL, action);
    }
    else if (action == OPEN_SCRIPT) {
	filesel_open_script(fname);
    }
    else if (action == OPEN_SESSION) {
	filesel_open_session(fname);
    }

    if (action < END_OPEN) return;

    /* now for the save options */

    if (action > SAVE_BIN2 && dat_ext(fname, 1)) return;

    if (check_maybe_add_ext(fname, action, data))
	return;

    if (SAVE_DATA_ACTION(action)) {
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
    else if (SAVE_TEX_ACTION(action)) {
	MODEL *pmod = (MODEL *) data;

	do_save_tex(fname, action, pmod); 
    }
    else {
	windata_t *vwin = (windata_t *) data;

	save_editable_content(action, fname, vwin);

    }
}

#else /* End of MS Windows file selection code, start GTK */

/* ........................................................... */

static void filesel_callback (GtkWidget *w, gpointer data) 
{
    GtkWidget *fs = GTK_WIDGET(data);
    gint action = GPOINTER_TO_INT(g_object_get_data
				  (G_OBJECT(data), "action"));
    char fname[MAXLEN];
    const gchar *path;
    FILE *fp = NULL;
    gpointer extdata = NULL;

    path = gtk_file_selection_get_filename(GTK_FILE_SELECTION(fs));
    if (path == NULL || *path == '\0') return;
    strcpy(fname, path);

    /* do some elementary checking */
    if (action < END_OPEN) {
	if ((fp = fopen(fname, "r")) == NULL) {
	    errbox(_("Couldn't open the specified file"));
	    return;
	} else fclose(fp);
    } 

    strcpy(remember_dir, path);

    if (OPEN_DATA_ACTION(action)) {
	strcpy(trydatfile, fname);
	gtk_widget_destroy(GTK_WIDGET(fs));  
	verify_open_data(NULL, action);
	return;
    }
    else if (APPEND_DATA_ACTION(action)) {
	strcpy(trydatfile, fname);
	gtk_widget_destroy(GTK_WIDGET(fs)); 
	do_open_data(NULL, NULL, action);
	return;
    }
    else if (action == OPEN_SCRIPT) {
	filesel_open_script(fname);
    }
    else if (action == OPEN_SESSION) {
	filesel_open_session(fname);
    }

    if (action < END_OPEN) {
	gtk_widget_destroy(GTK_WIDGET(fs));    
	return;
    }

    /* now for the save options */

    if (action == SAVE_GNUPLOT || action == SAVE_LAST_GRAPH || 
	action == SAVE_THIS_GRAPH) 
	extdata = g_object_get_data(G_OBJECT(fs), "graph");

    if (check_maybe_add_ext(fname, action, extdata))
	return;

    if (action > SAVE_BIN2 && dat_ext(fname, 1)) {
	gtk_widget_destroy(GTK_WIDGET(fs));
	return;
    }

    if (SAVE_DATA_ACTION(action)) {
	int overwrite = 0;

	if (!strcmp(fname, paths.datfile)) overwrite = 1;
	do_store(fname, action_to_flag(action), overwrite);
    }
    else if (action == SAVE_GNUPLOT) {
	int err = 0;
	GPT_SPEC *plot = g_object_get_data(G_OBJECT(fs), "graph");

	err = go_gnuplot(plot, fname, &paths);
	if (err == 0) infobox(_("graph saved"));
	else if (err == 1) errbox(_("gnuplot command failed"));
	else if (err == 2) infobox(_("There were missing observations"));
    }
#ifdef GNUPLOT_PNG
    else if (action == SAVE_THIS_GRAPH) {
	GPT_SPEC *plot = g_object_get_data(G_OBJECT(fs), "graph");

	save_this_graph(plot, fname);
    }
#endif
    else if (action == SAVE_BOXPLOT_EPS || action == SAVE_BOXPLOT_PS) {
	int err;

	err = ps_print_plots(fname, action,
			     g_object_get_data(G_OBJECT(fs), "graph"));
	if (!err) infobox(_("boxplots saved"));
	else errbox(_("boxplot save failed"));
    }
    else if (action == SAVE_BOXPLOT_XPM) {
	int err;

	err = plot_to_xpm(fname, g_object_get_data(G_OBJECT(fs), "graph"));
	if (!err) infobox(_("boxplots saved"));
	else errbox(_("boxplot save failed"));
    }
    else if (action == SAVE_LAST_GRAPH) {
	char *savestr = g_object_get_data(G_OBJECT(fs), "graph");
	
	do_save_graph(fname, savestr);
    }    
    else if (action == SAVE_SESSION) {
	save_session(fname);
    }
    else if (SAVE_TEX_ACTION(action)) {
	MODEL *pmod;
	pmod = (MODEL *) g_object_get_data(G_OBJECT(fs), "model");
	do_save_tex(fname, action, pmod); 
    }
    else {
	windata_t *vwin = g_object_get_data(G_OBJECT(fs), "text");

	save_editable_content(action, fname, vwin);

    }

    gtk_widget_destroy(GTK_WIDGET(fs));    
}

/* ........................................................... */

static void extra_get_filter (int action, gpointer data, char *suffix)
{
    
    const char *ext = get_ext(action, data);

    if (ext == NULL) { 
	strcpy(suffix, "*");
    } else {
	sprintf(suffix, "*%s", ext);
    }
}

/* ........................................................... */

#include <glob.h> /* POSIX */

static void
gtk_file_selection_glob_populate (GtkFileSelection *fs, 
				  const gchar *dirname, const gchar *suffix)
{
    glob_t globbuf;
    GtkTreeIter iter;
    GtkListStore *file_model;
    gchar *pattern;
    size_t dirlen;
    gint i;

    g_return_if_fail (GTK_IS_FILE_SELECTION (fs));

    if (dirname == NULL || *dirname == 0 || suffix == NULL || *suffix == 0)
	return;

    pattern = g_build_path(G_DIR_SEPARATOR_S, dirname, suffix, NULL);
    dirlen = strlen(pattern) - strlen(suffix);

    glob(pattern, 0, NULL, &globbuf);

    file_model = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(fs->file_list)));

    gtk_list_store_clear (file_model);

    for (i=0; i<globbuf.gl_pathc; i++) {
	gtk_list_store_append (file_model, &iter);
	gtk_list_store_set (file_model, &iter, 0, globbuf.gl_pathv[i] + dirlen, -1);
    }

    globfree(&globbuf);
    g_free(pattern);
}

void file_selector (char *msg, int action, gpointer data) 
{
    GtkWidget *filesel;
    char suffix[16], startdir[MAXLEN];
    int do_glob = 1;

    set_startdir(action, startdir);

    filesel = gtk_file_selection_new(msg);

    g_object_set_data(G_OBJECT(filesel), "action", GINT_TO_POINTER(action));

    g_signal_connect(G_OBJECT(GTK_FILE_SELECTION(filesel)->ok_button),
		     "clicked", 
		     G_CALLBACK(filesel_callback), filesel);

    if (action > END_OPEN) { /* a file save action */
	g_object_set_data(G_OBJECT(filesel), "text", data);
    }

    gtk_file_selection_set_filename(GTK_FILE_SELECTION(filesel), startdir);

    /* special cases */

    if (SAVE_GRAPH_ACTION(action)) {
	g_object_set_data(G_OBJECT(filesel), "graph", data);
    }

    else if (SAVE_TEX_ACTION(action)) {
	g_object_set_data(G_OBJECT(filesel), "model", data);
    }

    else if ((action == SAVE_DATA || action == SAVE_GZDATA) 
	     && paths.datfile[0]) {
	gtk_file_selection_set_filename(GTK_FILE_SELECTION(filesel), 
					paths.datfile);
	if (!strstr(paths.datfile, startdir)) do_glob = 0;
    }

    else if ((action == SAVE_SESSION) 
	     && scriptfile[0]) {
	gtk_file_selection_set_filename(GTK_FILE_SELECTION(filesel), 
					scriptfile);
	if (!strstr(scriptfile, startdir)) do_glob = 0;
    }

    /* end special cases */

    g_signal_connect(G_OBJECT(filesel), "destroy",
		     gtk_main_quit, NULL);
    g_signal_connect(G_OBJECT(GTK_FILE_SELECTION(filesel)->cancel_button),
		     "clicked", G_CALLBACK(delete_widget), filesel);

    gtk_widget_show(filesel);

    if (do_glob) {
	extra_get_filter(action, data, suffix);
	gtk_file_selection_glob_populate (GTK_FILE_SELECTION(filesel), 
					  startdir, suffix);
    }

    gtk_main(); /* make file selector modal */
}

#endif /* end of non-MS Windows code */

