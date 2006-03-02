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
#include "boxplots.h"
#include "gpt_control.h"
#include "session.h"
#include "textbuf.h"
#include "textutil.h"
#include "filelists.h"
#include "fileselect.h"

#if (GTK_MAJOR_VERSION >= 2) && (GTK_MINOR_VERSION >= 4)
# ifndef G_OS_WIN32
#  define USE_GTK_CHOOSER
# endif
#endif

#ifdef OLD_GTK
# include "menustate.h"
#endif

#define IS_DAT_ACTION(i) (i == SAVE_DATA || \
                          i == SAVE_DATA_AS || \
                          i == SAVE_GZDATA || \
                          i == SAVE_BIN1 || \
                          i == SAVE_BIN2 || \
                          i == OPEN_DATA)

#define OPEN_DATA_ACTION(i)  (i == OPEN_DATA || \
                              i == OPEN_CSV || \
                              i == OPEN_ASCII || \
                              i == OPEN_OCTAVE || \
	                      i == OPEN_BOX || \
                              i == OPEN_GNUMERIC || \
	                      i == OPEN_EXCEL || \
                              i == OPEN_WF1 || \
                              i == OPEN_DTA)

#define APPEND_DATA_ACTION(i) (i == APPEND_DATA || \
                               i == APPEND_CSV || \
                               i == APPEND_OCTAVE || \
                               i == APPEND_GNUMERIC || \
                               i == APPEND_EXCEL || \
                               i == APPEND_ASCII || \
                               i == APPEND_WF1 || \
                               i == APPEND_DTA)

#define SAVE_GRAPH_ACTION(i) (i == SAVE_GNUPLOT || \
                              i == SAVE_THIS_GRAPH || \
                              i == SAVE_LAST_GRAPH || \
                              i == SAVE_BOXPLOT_EPS || \
                              i == SAVE_BOXPLOT_PS || \
                              i == SAVE_BOXPLOT_XPM)

#define EXPORT_ACTION(a,s) ((a == EXPORT_OCTAVE || \
                             a == EXPORT_R || \
                             a == EXPORT_CSV || \
                             a == EXPORT_DAT) && s != FSEL_DATA_PRN)

#ifdef REMEMBER_DIR
static char remember_dir[MAXLEN];
#endif

struct extmap {
    int action;
    char *ext;
};

static struct extmap action_map[] = {
    { SAVE_DATA,         ".gdt" },
    { SAVE_DATA_AS,      ".gdt" },
    { SAVE_GZDATA,       ".gdt" },
    { SAVE_BIN1,         ".gdt" },
    { SAVE_BIN2,         ".gdt" },
    { SAVE_DBDATA,       ".bin" },
    { SAVE_SCRIPT,       ".inp" },
    { SAVE_CONSOLE,      ".inp" },
    { SAVE_SESSION,      ".gretl" },
    { SAVE_GP_CMDS,      ".plt" },
    { SAVE_BOXPLOT_EPS,  ".eps" },
    { SAVE_BOXPLOT_PS,   ".ps" },
    { SAVE_BOXPLOT_XPM,  ".xpm" },
    { EXPORT_CSV,        ".csv" },
    { EXPORT_R,          ".R" },
    { OPEN_OCTAVE,       ".m" },
    { APPEND_OCTAVE,     ".m" },
    { EXPORT_OCTAVE,     ".m" },
    { EXPORT_DAT,        ".dat" },
    { SAVE_OUTPUT,       ".txt" },
    { SAVE_TEX,          ".tex" },
    { SAVE_RTF,          ".rtf" },
    { OPEN_DATA,         ".gdt" },
    { APPEND_DATA,       ".gdt" },    
    { OPEN_SCRIPT,       ".inp" },
    { OPEN_SESSION,      ".gretl" },
    { OPEN_CSV,          ".csv" },
    { APPEND_CSV,        ".csv" },
    { OPEN_ASCII,        ".txt" },
    { APPEND_ASCII,      ".txt" },
    { OPEN_BOX,          ".box" },
    { OPEN_GNUMERIC,     ".gnumeric" },
    { APPEND_GNUMERIC,   ".gnumeric" },
    { OPEN_EXCEL,        ".xls" },
    { APPEND_EXCEL,      ".xls" },
    { OPEN_WF1,          ".wf1" },
    { APPEND_WF1,        ".wf1" },
    { OPEN_DTA,          ".dta" },
    { APPEND_DTA,        ".dta" },
    { FILE_OP_MAX,       NULL }
};

static gretlopt save_action_to_opt (int action, gpointer p)
{
    gretlopt opt = OPT_NONE;

    switch (action) {
    case SAVE_GZDATA:   opt = OPT_Z; break;
    case SAVE_BIN1:     opt = OPT_S; break;
    case SAVE_BIN2:     opt = OPT_O; break;
    case SAVE_DBDATA:   opt = OPT_D; break;
    case EXPORT_OCTAVE: opt = OPT_M; break;
    case EXPORT_R:      opt = OPT_R; break;
    case EXPORT_CSV:    opt = OPT_C; break;
    case EXPORT_DAT:    opt = OPT_G; break; /* PcGive */
    default: break;
    }

    if (p != NULL) {
	opt |= GPOINTER_TO_INT(p);
    }

    return opt;
}

static const char *get_gp_ext (const char *termtype)
{
    if (!strncmp(termtype, "postscript", 10))    return ".eps";
    else if (!strcmp(termtype, "PDF"))           return ".pdf";
    else if (!strcmp(termtype, "fig"))           return ".fig";
    else if (!strcmp(termtype, "latex"))         return ".tex";
    else if (!strncmp(termtype, "png", 3))       return ".png";
    else if (!strncmp(termtype, "emf", 3))       return ".emf";
    else if (!strcmp(termtype, "plot commands")) return ".plt";
    else return "*";
}

static int dat_ext (const char *str, int showerr)
{
    const char *suff;
    int err = 0;

    if (str == NULL) {
	return 0;
    }

    suff = strrchr(str, '.');

    if (suff != NULL && (!strcmp(suff, ".dat") || !strcmp(suff, ".gdt"))) {
	if (showerr) {
	    errbox(_("The suffix you selected should be used\n"
		   "only for gretl datafiles"));
	}
	err = 1;
    }

    return err;
}

static const char *get_ext (int action, gpointer data)
{
    const char *s = NULL;

    if (using_olddat() && IS_DAT_ACTION(action)) { 
	return ".dat";
    }

    if (action == SAVE_GNUPLOT || action == SAVE_THIS_GRAPH) {
	GPT_SPEC *plot = (GPT_SPEC *) data;
	s = get_gp_ext(plot->termtype);
    } else if (action == SAVE_LAST_GRAPH) {
	s = get_gp_ext(data);
    } else {
	size_t i;

	for (i=0; i < sizeof action_map / sizeof *action_map; i++) {
	    if (action == action_map[i].action) {
		s = action_map[i].ext;
		break;
	    }
	}
    }

    return s;
}

static int check_maybe_add_ext (char *fname, int action, gpointer data)
{
    FILE *fp;
    const char *ext = NULL;

    if (fname == NULL) {
	return 1;
    }

    /* don't mess if the fname is really a dir */
    if (isdir(fname)) {
	return 1;
    }

    /* don't mess with the name of a previously existing file */
    fp = gretl_fopen(fname, "r");
    if (fp != NULL && fgetc(fp) != EOF) {
	fclose(fp);
	return 0;
    }    

    /* don't mess with a filename that already has an extension */
    if (dotpos(fname) != strlen(fname)) {
	return 0;
    }
    
    /* otherwise add an appropriate extension */
    ext = get_ext(action, data);
    if (ext != NULL && strlen(ext) > 1) {
	strcat(fname, ext);
    }

    return 0;
}

static void script_window_update (windata_t *vwin, const char *fname)
{
    gchar *title;
    const char *p = strrchr(fname, SLASH);

    /* update the window title */
    title = g_strdup_printf("gretl: %s", p? p + 1 : fname);
    gtk_window_set_title(GTK_WINDOW(vwin->dialog), title);
    strcpy(vwin->fname, fname);
    g_free(title);

    if (vwin->role == VIEW_LOG || vwin->role == VIEW_SCRIPT) {
	vwin->role = EDIT_SCRIPT;
    }

    /* make the window editable */
#ifndef OLD_GTK
    if (!gtk_text_view_get_editable(GTK_TEXT_VIEW(vwin->w))) {
	file_view_set_editable(vwin);
    }
#else
    if (vwin->role == EDIT_SCRIPT) {
	file_view_set_editable(vwin);
    } 
#endif
}

static void save_editable_content (int action, const char *fname,
				   windata_t *vwin)
{
    FILE *fp;
    gchar *buf;
#if defined(ENABLE_NLS) && !defined(OLD_GTK)
    gchar *trbuf;
#endif

#ifndef OLD_GTK
    buf = textview_get_text(GTK_TEXT_VIEW(vwin->w));
#else
    buf = gtk_editable_get_chars(GTK_EDITABLE(vwin->w), 0, -1);
#endif

    if (buf == NULL) {
	errbox("Couldn't retrieve buffer");
	return;
    }

    fp = gretl_fopen(fname, "w");
    if (fp == NULL) {
	errbox(_("Couldn't open file for writing"));
	g_free(buf);
	return;
    }

#if defined(ENABLE_NLS) && !defined(OLD_GTK)
    trbuf = my_locale_from_utf8(buf);
    if (trbuf != NULL) {
	system_print_buf(trbuf, fp);
	g_free(trbuf);
    }
#else
    system_print_buf(buf, fp);
#endif

    g_free(buf);
    fclose(fp);
    
    if (action == SAVE_SCRIPT) {
	strcpy(scriptfile, fname);
	mkfilelist(FILE_LIST_SCRIPT, scriptfile);
	script_window_update(vwin, fname);
    }
}

static void set_startdir (char *startdir, int action)
{
#ifdef REMEMBER_DIR
    if (*remember_dir != '\0') {
	strcpy(startdir, remember_dir);
    } else {
	get_default_dir(startdir, action);
    }
#else
    get_default_dir(startdir, action);
#endif

#ifndef G_OS_WIN32
    if (startdir[strlen(startdir) - 1] != '/') {
	strcat(startdir, "/");
    }
#endif
}

static void filesel_save_prn_buffer (PRN *prn, const char *fname)
{
    FILE *fp = gretl_fopen(fname, "w");

    if (fp == NULL) {
	errbox(_("Couldn't write to %s"), fname);
    } else {
	const char *buf = gretl_print_get_buffer(prn);

	fputs(buf, fp);
	fclose(fp);
    }
}

static void filesel_open_script (const char *fname)
{
    int spos;

    strcpy(tryscript, fname);

    if (view_file(tryscript, 1, 0, 78, 370, EDIT_SCRIPT) != NULL) {
	strcpy(scriptfile, tryscript);
	mkfilelist(FILE_LIST_SCRIPT, scriptfile);
	spos = slashpos(scriptfile);
	if (spos) {
	    paths.currdir[0] = 0;
	    strncat(paths.currdir, scriptfile, spos + 1);
	}
    }
}

static void filesel_open_session (const char *fname)
{
    int pub = !strncmp(tryscript, paths.scriptdir, strlen(paths.scriptdir));

    strcpy(tryscript, fname);

    if (saved_objects(tryscript)) {
	verify_open_session(NULL);
    } else if (view_file(tryscript, 1, 0, 78, 370, 
			 pub ? VIEW_SCRIPT : EDIT_SCRIPT)) {
	strcpy(scriptfile, tryscript);
    }
}

static char *suggested_savename (const char *fname)
{
    const char *ss = strrchr(fname, SLASH);
    char *s, *sfx;

    if (ss == NULL) {
	s = g_strdup(fname);
    } else {
	s = g_strdup(ss + 1);
    }

    sfx = strrchr(s, '.');

    if (sfx != NULL && (strlen(sfx) == 4 || !strcmp(sfx, ".gnumeric"))) {
	const char *test = (using_olddat())? ".dat" : ".gdt";

	if (strcmp(test, sfx)) {
	    strcpy(sfx, test);
	}
    }

    return s;
}

static char *suggested_exportname (const char *fname, int action)
{
    const char *ss = strrchr(fname, SLASH);
    char *s, *sfx;

    if (ss == NULL) {
	s = g_strdup(fname);
    } else {
	s = g_strdup(ss + 1);
    }

    sfx = strrchr(s, '.');

    if (sfx != NULL && (strlen(sfx) == 4 || !strcmp(sfx, ".gnumeric"))) {
	const char *test;

	switch(action) {
	case EXPORT_OCTAVE:
	    test = ".m";
	    break;
	case EXPORT_R:
	    test = ".R";
	    break;
	case EXPORT_CSV:
	    test = ".csv";
	    break;
	case EXPORT_DAT:
	    test = ".dat";
	    break;
	default:
	    test = NULL;
	    break;
	}

	if (test != NULL && strcmp(test, sfx)) {
	    strcpy(sfx, test);
	}
    }

    return s;
}

#ifdef REMEMBER_DIR
static void remember_this_dir (const char *fname)
{
    int spos = slashpos(fname);
    
    if (spos > 0 && spos < sizeof remember_dir) {
	*remember_dir = '\0';
	strncat(remember_dir, fname, spos);
    }
}
#endif

static void
file_selector_process_result (const char *in_fname, int action, FselDataSrc src,
			      gpointer data)
{
    char fname[FILENAME_MAX];

    *fname = 0;
    strncat(fname, in_fname, FILENAME_MAX - 1);

    if (action < END_OPEN) {
	FILE *fp = gretl_fopen(fname, "r");

	if (fp == NULL) {
	    errbox(_("Couldn't open %s"), fname);
	    return;
	} else {
	    fclose(fp);
	}
    } 

#ifdef REMEMBER_DIR
    if (action != SET_PATH) {
	remember_this_dir(fname);
    }
#endif

    if (OPEN_DATA_ACTION(action)) {
	strcpy(trydatfile, fname);
	verify_open_data(NULL, action);
    } else if (APPEND_DATA_ACTION(action)) {
	strcpy(trydatfile, fname);
	do_open_data(NULL, NULL, action);
    } else if (action == OPEN_SCRIPT) {
	filesel_open_script(fname);
    } else if (action == OPEN_SESSION) {
	filesel_open_session(fname);
    } else if (action == OPEN_MARKERS) {
	do_add_markers(fname);
    }

    if (action < END_OPEN) {
	return;
    }

    /* now for the save/export options */

    if (action > SAVE_BIN2 && action != EXPORT_DAT && dat_ext(fname, 1)) { 
	return;
    }

    if (check_maybe_add_ext(fname, action, data)) {
	return;
    }

    if (src == FSEL_DATA_PRN) {
	if (action == SAVE_TEX) {
	    save_latex((PRN *) data, fname);
	} else {
	    filesel_save_prn_buffer((PRN *) data, fname);
	}
    } else if (SAVE_DATA_ACTION(action)) {
	int overwrite = 0;

	if (!strcmp(fname, paths.datfile) || action == SAVE_DBDATA) {
	    overwrite = 1;
	}
	do_store(fname, save_action_to_opt(action, data), overwrite);
    } else if (action == SAVE_GNUPLOT) {
	int err = 0;
	GPT_SPEC *plot = (GPT_SPEC *) data;

	err = go_gnuplot(plot, fname);
	if (err == 1) {
	    errbox(_("gnuplot command failed"));
	} else if (err == 2) {
	    infobox(_("There were missing observations"));
	}
    } else if (action == SAVE_THIS_GRAPH) {
	GPT_SPEC *plot = (GPT_SPEC *) data;

	save_this_graph(plot, fname);
    } else if (action == SAVE_BOXPLOT_EPS || action == SAVE_BOXPLOT_PS) {
	int err;

	err = ps_print_plots(fname, action, data);
	if (err) {
	    errbox(_("boxplot save failed"));
	}
    } else if (action == SAVE_SESSION) {
	save_session(fname);
    } else if (action == SET_PATH) {
	char *strvar = (char *) data;

	filesel_set_path_callback(fname, strvar);
    } else {
	windata_t *vwin = (windata_t *) data;

	save_editable_content(action, fname, vwin);
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
	{ N_("Windows metafiles (*.emf)"), "*.emf" },
	{ N_("gnuplot files (*.plt)"), "*.plt" },
	{ N_("all files (*.*)"), "*.*" }
    };

    if (!strncmp(termtype, "postscript", 10)) 
	return gpfilters[0];
    else if (!strcmp(termtype, "fig")) 
	return gpfilters[1];
    else if (!strcmp(termtype, "latex")) 
	return gpfilters[2];
    else if (!strncmp(termtype, "png", 3)) 
	return gpfilters[3];
    else if (!strncmp(termtype, "emf", 3)) 
	return gpfilters[4];
    else if (!strcmp(termtype, "plot commands")) 
	return gpfilters[5];
    else 
	return gpfilters[6];
}

static struct winfilter get_filter (int action, gpointer data)
{
    int i;
    struct winfilter filter;
    static struct win32_filtermap map[] = {
	{SAVE_DATA,    { N_("gretl data files (*.gdt)"), "*.gdt" }},
	{SAVE_DATA_AS, { N_("gretl data files (*.gdt)"), "*.gdt" }},
	{SAVE_GZDATA,  { N_("gretl data files (*.gdt)"), "*.gdt" }},
	{SAVE_BIN1,    { N_("gretl data files (*.gdt)"), "*.gdt" }},
	{SAVE_BIN2,    { N_("gretl data files (*.gdt)"), "*.gdt" }},
	{SAVE_DBDATA,  { N_("gretl database files (*.bin)"), "*.bin" }},
	{SAVE_SCRIPT,  { N_("gretl script files (*.inp)"), "*.inp" }},
	{SAVE_CONSOLE, { N_("gretl command files (*.inp)"), "*.inp" }},
	{SAVE_SESSION, { N_("session files (*.gretl)"), "*.gretl" }},
	{SAVE_BOXPLOT_EPS, { N_("postscript files (*.eps)"), "*.eps" }},
	{SAVE_BOXPLOT_PS,  { N_("postscript files (*.ps)"), "*.ps" }},
	{SAVE_GP_CMDS, { N_("gnuplot files (*.plt)"), "*.plt" }},
	{EXPORT_CSV,   { N_("CSV files (*.csv)"), "*.csv" }},
	{EXPORT_R,     { N_("GNU R files (*.R)"), "*.R" }},
	{EXPORT_OCTAVE, { N_("GNU Octave files (*.m)"), "*.m" }},
	{OPEN_OCTAVE,  { N_("GNU Octave files (*.m)"), "*.m" }},
	{APPEND_OCTAVE, { N_("GNU Octave files (*.m)"), "*.m" }},
	{EXPORT_DAT,   { N_("PcGive files (*.dat)"), "*.dat" }},
	{SAVE_OUTPUT,  { N_("text files (*.txt)"), "*.txt" }},
	{SAVE_TEX,     { N_("TeX files (*.tex)"), "*.tex" }},
	{SAVE_RTF,     { N_("RTF files (*.rtf)"), "*.rtf" }},
	{OPEN_DATA,    { N_("gretl data files (*.gdt)"), "*.gdt" }},
	{APPEND_DATA,  { N_("gretl data files (*.gdt)"), "*.gdt" }},
	{OPEN_SCRIPT,  { N_("gretl script files (*.inp)"), "*.inp" }},
	{OPEN_SESSION, { N_("session files (*.gretl)"), "*.gretl" }},
	{OPEN_CSV,     { N_("CSV files (*.csv)"), "*.csv" }},
	{APPEND_CSV,   { N_("CSV files (*.csv)"), "*.csv" }},
	{OPEN_ASCII,   { N_("ASCII files (*.txt)"), "*.txt" }},
	{APPEND_ASCII, { N_("ASCII files (*.txt)"), "*.txt" }},
	{OPEN_BOX,     { N_("BOX data files (*.box)"), "*.box" }},
	{OPEN_GNUMERIC,   { N_("Gnumeric files (*.gnumeric)"), "*.gnumeric" }},
	{APPEND_GNUMERIC, { N_("Gnumeric files (*.gnumeric)"), "*.gnumeric" }},
	{OPEN_EXCEL,   { N_("Excel files (*.xls)"), "*.xls" }},
	{APPEND_EXCEL, { N_("Excel files (*.xls)"), "*.xls" }},
	{OPEN_WF1,     { N_("Eviews workfiles (*.wf1)"), "*.wf1" }},
	{APPEND_WF1,   { N_("Eviews workfiles (*.wf1)"), "*.wf1" }},
	{OPEN_DTA,     { N_("Stata files (*.dta)"), "*.dta" }},
	{APPEND_DTA,   { N_("Stata files (*.dta)"), "*.dta" }},
	{SET_PATH,     { N_("program files (*.exe)"), "*.exe" }}
    };

    static struct winfilter olddat_filter = {
	N_("gretl data files (*.dat)"), "*.dat"
    };    

    static struct winfilter default_filter = {
	N_("all files (*.*)"), "*.*" 
    };

    if (using_olddat() && IS_DAT_ACTION(action)) {
	return olddat_filter;
    }

    if (action == SAVE_GNUPLOT || action == SAVE_THIS_GRAPH) {
	GPT_SPEC *plot = (GPT_SPEC *) data;

	return get_gp_filter(plot->termtype);
    }

    if (action == SAVE_LAST_GRAPH) {
	return get_gp_filter(data);
    }

    filter = default_filter;

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

    if (p == NULL) {
	return NULL;
    }

    filter = get_filter(action, data);

    strcpy(p, I_(filter.descrip));
    p += strlen(p) + 1;
    strcpy(p, filter.pat);

    if (strncmp(filter.descrip, "all", 3)) {
	p += strlen(p) + 1;
	strcpy(p, I_("all files (*.*)"));
	p += strlen(p) + 1;
	strcpy(p, "*.*");
    }

    p += strlen(p) + 1;
    *p = '\0';

    return start;
}

void file_selector (const char *msg, int action, FselDataSrc src, gpointer data) 
{
    OPENFILENAME of;
    int retval;
    char fname[MAXLEN], endname[64], startdir[MAXLEN];
    char *filter = NULL;
    gchar *trmsg = NULL;

    *fname = '\0';
    *endname = '\0';

    set_startdir(startdir, action);

    /* special cases */
    if ((action == SAVE_DATA || action == SAVE_GZDATA) && paths.datfile[0]) {
	char *savename = suggested_savename(paths.datfile);

	strcpy(fname, savename);
	g_free(savename);
	if (!(data_status & BOOK_DATA)) {
	    get_base(startdir, paths.datfile, SLASH);
	}
    } else if (EXPORT_ACTION(action, src) && paths.datfile[0]) {
	char *savename = suggested_exportname(paths.datfile, action);

	strcpy(fname, savename);
	g_free(savename);
	get_base(startdir, paths.datfile, SLASH);
    } else if (action == SET_PATH) {
	char *strvar = (char *) data;

	if (strvar != NULL && *strvar != '\0') {
	    if (get_base(startdir, strvar, SLASH)) {
		strcpy(fname, strvar + slashpos(strvar) + 1);
	    } 
	}
    }	

    if (doing_nls()) {
	trmsg = my_locale_from_utf8(msg);
    } else {
	trmsg = g_strdup(msg);
    }

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

    if (action < END_OPEN) {
	retval = GetOpenFileName(&of);
    } else {
	/* a file save action */
	retval = GetSaveFileName(&of);
    }

    free(filter);
    g_free(trmsg);

    if (!retval) {
	if (CommDlgExtendedError()) {
	    errbox(_("File dialog box error"));
	}
	return;
    }

    my_filename_to_utf8(fname);

    file_selector_process_result(fname, action, src, data);
}

#else /* End of MS Windows file selection code, start GTK */

/* ........................................................... */

# ifndef USE_GTK_CHOOSER

struct fsinfo_t {
    GtkWidget *w;
    char fname[FILENAME_MAX];
};

static void filesel_callback (GtkWidget *w, struct fsinfo_t *fsinfo) 
{
# ifdef OLD_GTK
    GtkIconFileSel *fsel = GTK_ICON_FILESEL(fsinfo->w);
    char *test;
# endif
    const gchar *path;

    *fsinfo->fname = '\0';

# ifndef OLD_GTK
    path = gtk_file_selection_get_filename(GTK_FILE_SELECTION(fsinfo->w));
    if (path == NULL || *path == '\0' || isdir(path)) {
	return;
    }
    strcpy(fsinfo->fname, path);
# else
    test = gtk_entry_get_text(GTK_ENTRY(fsel->file_entry));
    if (test == NULL || *test == '\0') {
	return;
    }
    path = gtk_file_list_get_path(GTK_FILE_LIST(fsel->file_list));
    sprintf(fsinfo->fname, "%s%s", path, test);
# endif

    gtk_widget_destroy(GTK_WIDGET(fsinfo->w));
}

# endif /* !USE_GTK_CHOOSER */

static char *get_filter_suffix (int action, gpointer data, char *suffix)
{
    
    const char *ext = get_ext(action, data);

    if (ext == NULL) { 
	strcpy(suffix, "*");
    } else {
	sprintf(suffix, "*%s", ext);
    }

    return suffix;
}

# ifndef OLD_GTK

# ifdef USE_GTK_CHOOSER

static GtkFileFilter *get_file_filter (int action, gpointer data)
{
    GtkFileFilter *filter;
    char suffix[16];
    
    filter = gtk_file_filter_new();
    get_filter_suffix(action, data, suffix);
    gtk_file_filter_add_pattern(filter, suffix);

    return filter;
}

void file_selector (const char *msg, int action, FselDataSrc src, gpointer data) 
{
    GtkWidget *filesel;
    char startdir[MAXLEN];
    GtkFileFilter *filter;

    set_startdir(startdir, action);

    /* FIXME: parent window below should not always be mdata->w,
       in particular when action == SAVE_THIS_GRAPH
    */

    if (action > END_OPEN && action != SET_PATH) {
	filesel = gtk_file_chooser_dialog_new(msg, NULL, /* GTK_WINDOW(mdata->w), */
					      GTK_FILE_CHOOSER_ACTION_SAVE,
					      GTK_STOCK_SAVE, GTK_RESPONSE_ACCEPT,
					      GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
					      NULL);
    } else {
	filesel = gtk_file_chooser_dialog_new(msg, GTK_WINDOW(mdata->w), 
					      GTK_FILE_CHOOSER_ACTION_OPEN,
					      GTK_STOCK_OPEN, GTK_RESPONSE_ACCEPT,
					      GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
					      NULL);
    }

    gtk_dialog_set_default_response(GTK_DIALOG(filesel), GTK_RESPONSE_ACCEPT);

    filter = get_file_filter(action, data);
    gtk_file_chooser_set_filter(GTK_FILE_CHOOSER(filesel), filter);
    gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(filesel), startdir);

    /* special cases */

    if ((action == SAVE_DATA || action == SAVE_GZDATA) && paths.datfile[0]) {
	char *savename = suggested_savename(paths.datfile);

	gtk_file_chooser_set_current_name(GTK_FILE_CHOOSER(filesel), 
					  savename);
	g_free(savename);
    } else if (EXPORT_ACTION(action, src) && paths.datfile[0]) {
	char *savename = suggested_exportname(paths.datfile, action);

	gtk_file_chooser_set_current_name(GTK_FILE_CHOOSER(filesel), 
					  savename);
	g_free(savename);
    } else if ((action == SAVE_SESSION) && *scriptfile != '\0') {
	gtk_file_chooser_set_filename(GTK_FILE_CHOOSER(filesel), 
				      scriptfile);
    } else if (action == SET_PATH) {
	char *strvar = (char *) data;

	if (strvar != NULL && slashpos(strvar) > 0) {
	    gtk_file_chooser_set_filename(GTK_FILE_CHOOSER(filesel), 
					  strvar);
	} else {
	    gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(filesel), 
						"/usr/bin");
	}	    
    } 

    if (gtk_dialog_run(GTK_DIALOG(filesel)) == GTK_RESPONSE_ACCEPT) {
	char *fname;

	fname = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(filesel));
	gtk_widget_destroy(filesel);
	file_selector_process_result(fname, action, src, data);
	g_free(fname);
    } else {
	gtk_widget_destroy(filesel);
    }
}

# else

#include <glob.h> /* POSIX */

static void
gtk_file_selection_glob_populate (GtkFileSelection *fs, 
				  const gchar *dirname, const gchar *suffix)
{
    glob_t globbuf;
    GtkTreeIter iter;
    GtkListStore *file_model;
    gchar *pattern;
    size_t i, dirlen;

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

void file_selector (const char *msg, int action, FselDataSrc src, gpointer data) 
{
    struct fsinfo_t fsinfo;
    GtkWidget *filesel;
    char suffix[16], startdir[MAXLEN];
    int do_glob = 1;

    *fsinfo.fname = '\0';

    set_startdir(startdir, action);

    filesel = gtk_file_selection_new(msg);
    fsinfo.w = filesel;
    g_signal_connect(G_OBJECT(GTK_FILE_SELECTION(filesel)->ok_button),
		     "clicked", 
		     G_CALLBACK(filesel_callback), &fsinfo);

    gtk_file_selection_set_filename(GTK_FILE_SELECTION(filesel), startdir);

    /* special cases */

    if ((action == SAVE_DATA || action == SAVE_GZDATA) && paths.datfile[0]) {
	char *savename = suggested_savename(paths.datfile);

	gtk_file_selection_set_filename(GTK_FILE_SELECTION(filesel), 
					savename);
	if (!strstr(savename, startdir)) do_glob = 0;
	g_free(savename);
    } else if (EXPORT_ACTION(action, src) && paths.datfile[0]) {
	char *savename = suggested_exportname(paths.datfile, action);

	gtk_file_selection_set_filename(GTK_FILE_SELECTION(filesel), 
					savename);
	if (!strstr(savename, startdir)) do_glob = 0;
	g_free(savename);
    } else if ((action == SAVE_SESSION) && *scriptfile != '\0') {
	gtk_file_selection_set_filename(GTK_FILE_SELECTION(filesel), 
					scriptfile);
	if (!strstr(scriptfile, startdir)) do_glob = 0;
    } else if (action == SET_PATH) {
	char *strvar = (char *) data;

	if (strvar != NULL && slashpos(strvar) > 0) {
	    gtk_file_selection_set_filename(GTK_FILE_SELECTION(filesel), 
					    strvar);
	} else {
	    gtk_file_selection_set_filename(GTK_FILE_SELECTION(filesel), 
					    "/usr/bin/");
	}	    
	do_glob = 0;
    }	

    /* end special cases */

    g_signal_connect(G_OBJECT(filesel), "destroy",
		     gtk_main_quit, NULL);
    g_signal_connect(G_OBJECT(GTK_FILE_SELECTION(filesel)->cancel_button),
		     "clicked", G_CALLBACK(delete_widget), filesel);

    gtk_widget_show(filesel);

    if (do_glob) {
	get_filter_suffix(action, data, suffix);
	gtk_file_selection_glob_populate (GTK_FILE_SELECTION(filesel), 
					  startdir, suffix);
    }

    gtk_window_set_modal(GTK_WINDOW(filesel), TRUE);
    gtk_main(); 

    if (*fsinfo.fname != '\0') {
	file_selector_process_result(fsinfo.fname, action, src, data);
    } 
}

# endif

# else /* gtk version diffs continue */

void file_selector (const char *msg, int action, FselDataSrc src, gpointer data) 
{
    struct fsinfo_t fsinfo;
    GtkWidget *filesel;
    int gotdir = 0;
    char suffix[8], startdir[MAXLEN];

    *fsinfo.fname = '\0';

    set_startdir(startdir, action);

    filesel = gtk_icon_file_selection_new(msg);

    if (strstr(startdir, "/.")) {
	gtk_icon_file_selection_show_hidden(GTK_ICON_FILESEL(filesel), TRUE);
    } else {
	gtk_icon_file_selection_show_hidden(GTK_ICON_FILESEL(filesel), FALSE);
    }

    get_filter_suffix(action, data, suffix);
    gtk_icon_file_selection_set_filter(GTK_ICON_FILESEL(filesel), suffix);

    fsinfo.w = filesel;
    gtk_signal_connect(GTK_OBJECT(GTK_ICON_FILESEL(filesel)->ok_button),
		       "clicked", 
		       GTK_SIGNAL_FUNC(filesel_callback), &fsinfo);

    /* special cases */

    if ((action == SAVE_DATA || action == SAVE_GZDATA) && paths.datfile[0]) {
	char *savename = suggested_savename(paths.datfile);
	char startd[MAXLEN];

	gtk_entry_set_text(GTK_ENTRY(GTK_ICON_FILESEL(filesel)->file_entry),
			   savename);
	if (!(data_status & BOOK_DATA) && get_base(startd, paths.datfile, SLASH)) {
	    gtk_icon_file_selection_open_dir(GTK_ICON_FILESEL(filesel), startd);
	    gotdir = 1;
	}
	g_free(savename);
    } else if (EXPORT_ACTION(action, src) && paths.datfile[0]) {
	char *savename = suggested_exportname(paths.datfile, action);
	char startd[MAXLEN];

	gtk_entry_set_text(GTK_ENTRY(GTK_ICON_FILESEL(filesel)->file_entry),
			   savename);
	if (get_base(startd, paths.datfile, SLASH)) {
	    gtk_icon_file_selection_open_dir(GTK_ICON_FILESEL(filesel), startd);
	    gotdir = 1;
	}
	g_free(savename);
    } else if (action == SET_PATH) {
	char *strvar = (char *) data;
	char startd[MAXLEN];

	if (get_base(startd, strvar, SLASH) == 1) {
	    gtk_icon_file_selection_open_dir(GTK_ICON_FILESEL(filesel), startd);
	    gtk_entry_set_text(GTK_ENTRY(GTK_ICON_FILESEL(filesel)->file_entry),
			       strvar + slashpos(strvar) + 1);
	} else {
	    gtk_icon_file_selection_open_dir(GTK_ICON_FILESEL(filesel), 
					     "/usr/bin/");
	}
	gotdir = 1;
    }

    if (!gotdir) {
	gtk_icon_file_selection_open_dir(GTK_ICON_FILESEL(filesel), startdir);
    }

    gtk_signal_connect(GTK_OBJECT(GTK_ICON_FILESEL(filesel)), "destroy",
                       gtk_main_quit, NULL);
    gtk_signal_connect_object(GTK_OBJECT(GTK_ICON_FILESEL(filesel)->cancel_button),
			      "clicked", (GtkSignalFunc) gtk_widget_destroy,
			      GTK_OBJECT(filesel));

    gtk_widget_show(filesel);
    gretl_set_window_modal(filesel);
    gtk_main(); 

    if (*fsinfo.fname != '\0') {
	file_selector_process_result(fsinfo.fname, action, src, data);
    } 
}

# endif /* old gtk */

#endif /* end of non-MS Windows code */



