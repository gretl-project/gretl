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

extern int olddat; /* gui_utils.c */

static char remember_dir[MAXLEN];

struct extmap {
    int action;
    char *ext;
};

static struct extmap action_map[] = {
    {SAVE_DATA, ".gdt"},
    {SAVE_GZDATA, ".gz"},
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
    {OPEN_BOX, ".box"},
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

static char *get_gp_ext (const char *termtype)
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
    if (i == SAVE_DATA || i == SAVE_BIN1 || i == SAVE_BIN2 ||
	i == OPEN_DATA)
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
	    errbox("The suffix you selected should be used\n"
		   "only for gretl datafiles");
	return 1;
    }
    return 0;
}

/* ........................................................... */

static char *get_ext (int action, gpointer data)
{
    char *s = NULL;

    if (olddat && is_data_action(action)) 
	return ".dat";

    if (action == SAVE_GNUPLOT) {
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
    char *ext = NULL;

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

          /* MS Windows version of file selection code */

/* ........................................................... */

#ifdef G_OS_WIN32 

#include <windows.h>

static char *get_gp_filter (const char *termtype)
{
    if (!strcmp(termtype, "postscript")) 
	return "postscript files\0*.eps\0";
    else if (!strcmp(termtype, "fig")) 
	return "xfig files\0*.fig\0";
    else if (!strcmp(termtype, "latex")) 
	return "LaTeX files\0*.tex\0";
    else if (!strcmp(termtype, "png")) 
	return "PNG files\0*.png\0";
    else if (!strcmp(termtype, "plot commands")) 
	return "gnuplot files\0*.gp\0";
    else return "all files\0*\0";
}

struct win32_filtermap {
    int action;
    char *filter;
};

static char *get_filter (int action, gpointer data)
{
    int i;
    char *filter;
    static struct win32_filtermap map[] = {
	{SAVE_DATA, "gretl data files (*.gdt)\0*.gdt\0all files\0*\0"},
	{SAVE_GZDATA, "compressed data files (*.gz)\0*.gz\0all files\0*\0"},
	{SAVE_BIN1, "gretl data files (*.gdt)\0*.gdt\0all files\0*\0"},
	{SAVE_BIN2, "gretl data files (*.gdt)\0*.gdt\0all files\0*\0"},
	{SAVE_CMDS, "gretl command files (*.inp)\0*.inp\0all files\0*\0"},
	{SAVE_SCRIPT, "gretl script files (*.inp)\0*.inp\0all files\0*\0"},
	{SAVE_CONSOLE, "gretl command files (*.inp)\0*.inp\0all files\0*\0"},
	{SAVE_MODEL, "text files (*.txt)\0*.txt\0all files\0*\0"},
	{SAVE_SESSION, "session files (*.gretl)\0*.gretl\0all files\0*\0"},
	{SAVE_BOXPLOT_EPS, "postscript files (*.eps)\0*.eps\0all files\0*\0"},
	{SAVE_BOXPLOT_PS, "postscript files (*.ps)\0*.ps\0all files\0*\0"},
	{SAVE_LAST_GRAPH, "all files\0*\0"},
	{SAVE_GP_CMDS, "gnuplot files (*.gp)\0*.gp\0all files\0*\0"},
	{EXPORT_CSV, "CSV files (*.csv)\0*.csv\0all files\0*\0"},
	{EXPORT_R, "GNU R files (*.R)\0*.R\0all files\0*\0"},
	{EXPORT_R_ALT, "GNU R files (*.R)\0*.R\0all files\0*\0"},
	{EXPORT_OCTAVE, "GNU Octave files (*.m)\0*.m\0all files\0*\0"},
	{SAVE_OUTPUT, "text files (*.txt)\0*.txt\0all files\0*\0"},
	{SAVE_TEX_TAB, "TeX files (*.tex)\0*.tex\0all files\0*\0"},
	{SAVE_TEX_EQ, "TeX files (*.tex)\0*.tex\0all files\0*\0"},
	{OPEN_DATA, "gretl data files (*.gdt)\0*.gdt*\0all files\0*\0"},
	{OPEN_SCRIPT, "gretl script files (*.inp)\0*.inp\0all files\0*\0"},
	{OPEN_SESSION, "session files (*.gretl)\0*.gretl\0all files\0*\0"},
	{OPEN_CSV,  "CSV files (*.csv)\0*.csv\0all files\0*\0"},
	{OPEN_BOX, "BOX data files (*.box)\0*.box\0all files\0*\0"}};

    if (olddat && is_data_action(action)) 
	return "gretl data files (*.dat)\0*.dat\0all files\0*\0";

    if (action == SAVE_GNUPLOT) {
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

/* ........................................................... */

void file_selector (char *msg, int action, gpointer data) 
{
    OPENFILENAME of;
    int retval, gotdir = 0;
    char fname[MAXLEN], endname[64], startd[MAXLEN];
    char title[48];

    fname[0] = '\0';
    endname[0] = '\0';
    if (remember_dir[0] != '\0')
	strcpy(startd, remember_dir);
    else
	get_default_dir(startd);

    /* special case: default save of data */
    if (action == SAVE_DATA && paths.datfile[0] &&
	!strcmp(paths.datfile + strlen(paths.datfile) - 4, 
		(olddat)? ".dat" : ".gdt")) {
	strcpy(fname, paths.datfile + slashpos(paths.datfile) + 1);
	get_base(startd, paths.datfile, SLASH);
    }

    /* initialize file dialog info struct */
    memset(&of, 0, sizeof of);
#ifdef OPENFILENAME_SIZE_VERSION_400
    of.lStructSize = OPENFILENAME_SIZE_VERSION_400;
#else
    of.lStructSize = sizeof of;
#endif
    of.hwndOwner = NULL;
    of.lpstrFilter = get_filter(action, data);
    of.lpstrCustomFilter = NULL;
    of.nFilterIndex = 1;
    of.lpstrFile = fname;
    of.nMaxFile = sizeof fname;
    of.lpstrFileTitle = endname;
    of.nMaxFileTitle = sizeof endname;
    of.lpstrInitialDir = startd;
    of.lpstrTitle = msg;
    of.lpstrDefExt = NULL;
    of.Flags = OFN_HIDEREADONLY;

    if (action < END_OPEN)
	retval = GetOpenFileName(&of);
    else  /* a file save action */
	retval = GetSaveFileName(&of);

    if (!retval) {
	if (CommDlgExtendedError())
	    errbox("File dialog box error");
	return;
    }
	
    strcpy(title, "gretl: ");
    strncat(title, of.lpstrFileTitle, 40);

    strncpy(remember_dir, fname, slashpos(fname));

    if (action == OPEN_DATA || action == OPEN_CSV || action == OPEN_BOX) {
	strcpy(trydatfile, fname);
	verify_open_data(NULL);
    }
    else if (action == OPEN_SCRIPT) {
	int spos;

	strcpy(tryscript, fname);

	if (view_file(tryscript, 1, 0, 78, 370, title, script_items) != NULL) {
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
	if (view_file(tryscript, 1, 0, 78, 370, title, pub ? 
		      sample_script_items : script_items))
		strcpy(scriptfile, tryscript);
    }

    if (action < END_OPEN) return;

    /* now for the save options */

    if (action > SAVE_BIN2 && dat_ext(fname, 1)) return;

    maybe_add_ext(fname, action, data);

    if (action >= SAVE_DATA && action < END_SAVE_DATA) {
	do_store(fname, action_to_flag(action));
    }
    else if (action == SAVE_GNUPLOT) {
	int err = 0;
	GPT_SPEC *plot = (GPT_SPEC *) data;

	err = go_gnuplot(plot, fname, &paths);
	if (err == 0) infobox("graph saved");
	else if (err == 1) errbox("gnuplot command failed");
	else if (err == 2) infobox("There were missing observations");
    }
    else if (action == SAVE_BOXPLOT_EPS || action == SAVE_BOXPLOT_PS) {
	int err;

	err = ps_print_plots(fname, action, data);
	if (!err) infobox("boxplots saved");
	else errbox("boxplot save failed");
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
    else { /* save contents of an editable text window */
	GtkWidget *editwin;
	FILE *fp;

	editwin = (GtkWidget *) data;
	if ((fp = fopen(fname, "w")) == NULL) {
	    errbox("Couldn't open file for writing");
	    return;
	}
	fprintf(fp, "%s", 
		gtk_editable_get_chars(GTK_EDITABLE(editwin), 0, -1));
	fclose(fp);
	infobox("File saved OK");
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
    char *test, *path, title[48];
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
	    errbox("Couldn't open specified file");
	    return;
	} else fclose(fp);
    } 

    strcpy(title, "gretl: ");
    strncat(title, test, 40);
    strcpy(remember_dir, path);

    if (action == OPEN_DATA || action == OPEN_CSV || action == OPEN_BOX) {
	strcpy(trydatfile, fname);
	verify_open_data(NULL);
    }
    else if (action == OPEN_SCRIPT) {
	int spos;

	strcpy(tryscript, fname);

	if (view_file(tryscript, 1, 0, 78, 370, title, script_items) != NULL) {
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
	if (view_file(tryscript, 1, 0, 78, 370, title, pub ? 
		      sample_script_items : script_items))
		strcpy(scriptfile, tryscript);
    }

    if (action < END_OPEN) {
	gtk_widget_destroy(GTK_WIDGET(fs));    
	return;
    }

    /* now for the save options */
    if (action == SAVE_GNUPLOT || action == SAVE_LAST_GRAPH) 
	extdata = gtk_object_get_data(GTK_OBJECT(fs), "graph");

    maybe_add_ext(fname, action, extdata); 

    if (action > SAVE_BIN2 && dat_ext(fname, 1)) {
	gtk_widget_destroy(GTK_WIDGET(fs));
	return;
    }

    if (action >= SAVE_DATA && action < END_SAVE_DATA) {
	do_store(fname, action_to_flag(action));
    }
    else if (action == SAVE_GNUPLOT) {
	int err = 0;
	GPT_SPEC *plot = gtk_object_get_data(GTK_OBJECT(fs), "graph");

	err = go_gnuplot(plot, fname, &paths);
	if (err == 0) infobox("graph saved");
	else if (err == 1) errbox("gnuplot command failed");
	else if (err == 2) infobox("There were missing observations");
    }
    else if (action == SAVE_BOXPLOT_EPS || action == SAVE_BOXPLOT_PS) {
	int err;

	err = ps_print_plots(fname, action,
			     gtk_object_get_data(GTK_OBJECT(fs), "graph"));
	if (!err) infobox("boxplots saved");
	else errbox("boxplot save failed");
    }
    else if (action == SAVE_BOXPLOT_XPM) {
	int err;

	err = plot_to_xpm(fname, gtk_object_get_data(GTK_OBJECT(fs), "graph"));
	if (!err) infobox("boxplots saved");
	else errbox("boxplot save failed");
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
    else { /* save contents of an editable text window */
	GtkWidget *editwin;
	FILE *fp;

	editwin = gtk_object_get_data(GTK_OBJECT(fs), "text");
	if ((fp = fopen(fname, "w")) == NULL) {
	    errbox("Couldn't open file for writing");
	    gtk_widget_destroy(GTK_WIDGET(fs));
	    return;
	}
	fprintf(fp, "%s", 
		gtk_editable_get_chars(GTK_EDITABLE(editwin), 0, -1));
	fclose(fp);
	infobox("File saved OK");
    }
    gtk_widget_destroy(GTK_WIDGET(fs));    
}

/* ........................................................... */

static void extra_get_filter (int action, gpointer data, char *suffix)
{
    
    char *ext = get_ext(action, data);

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

    if (action == SAVE_GNUPLOT || action == SAVE_LAST_GRAPH
	|| action == SAVE_BOXPLOT_EPS || action == SAVE_BOXPLOT_PS ||
	action == SAVE_BOXPLOT_XPM) 
	gtk_object_set_data(GTK_OBJECT(filesel), "graph", data);

    else if (action == SAVE_TEX_TAB || action == SAVE_TEX_EQ) 
	gtk_object_set_data(GTK_OBJECT(filesel), "model", data);

    else if (action == SAVE_DATA && paths.datfile[0] &&
	     dat_ext(paths.datfile, 0)) {
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

    gtk_signal_connect_object(GTK_OBJECT(GTK_ICON_FILESEL
					 (filesel)->cancel_button),
			      "clicked", (GtkSignalFunc) gtk_widget_destroy,
			      GTK_OBJECT (filesel));

    gtk_widget_show(filesel);
}

#endif /* end of non-MS Windows code */

