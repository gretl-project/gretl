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

#include "gretl.h"
#include "textbuf.h"
#include "textutil.h"
#include "fileselect.h"
#include "winstack.h"
#include "tabwin.h"
#include "base_utils.h"

#ifndef GRETL_EDIT
#include "gui_utils.h"
#include "filelists.h"
#include "gpt_control.h"
#include "gpt_dialog.h"
#include "session.h"
#include "fnsave.h"
#include "database.h"
#include "datafiles.h"
#include "guiprint.h"
#include "graphics.h"
#include "bootstrap.h"
#endif

#include <sys/stat.h>
#include <unistd.h>

#ifdef G_OS_WIN32
# include "gretlwin32.h"
#endif

#ifdef __APPLE__
# include "osx_open.h"
#endif

#define EXPORT_OTHER(a) (a == EXPORT_OCTAVE ||	\
			 a == EXPORT_R ||	\
			 a == EXPORT_CSV ||	\
			 a == EXPORT_DAT ||	\
			 a == EXPORT_DTA ||	\
			 a == EXPORT_JM)

#define SET_DIR_ACTION(i) (i == SET_DIR || i == SET_WDIR || \
                           i == SET_FDIR || i == SET_DBDIR)

struct extmap {
    int action;
    char *ext;
};

static struct extmap action_map[] = {
    { EXPORT_DB,         ".bin" },
    { SAVE_SCRIPT,       ".inp" },
    { SAVE_FUNCTIONS_AS, ".inp" },
    { SAVE_CMD_LOG,      ".inp" },
    { SAVE_CONSOLE,      ".txt" },
    { SAVE_SESSION,      ".gretl" },
    { SAVE_GP_CMDS,      ".plt" },
    { SAVE_R_CMDS,       ".R" },
    { SAVE_OX_CMDS,      ".ox" },
    { SAVE_OCTAVE_CMDS,  ".m" },
    { SAVE_PYTHON_CMDS,  ".py" },
    { SAVE_STATA_CMDS,   ".do" },
    { SAVE_JULIA_CODE,   ".jl" },
    { SAVE_DYNARE_CODE,  ".mod" },
    { SAVE_LPSOLVE_CODE, ".lp" },
    { SAVE_SPEC_FILE,    ".spec" },
    { SAVE_X13_SPC,      ".spc" },
    { SAVE_FUNCTIONS,    ".gfn" },
    { SAVE_MARKERS,      ".txt" },
    { SAVE_LABELS,       ".txt" },
    { SAVE_GFN_SPEC,     ".spec" },
    { SAVE_GFN_ZIP,      ".zip" },
    { SAVE_MAP,          ".geojson" },
    { WRITE_MAP,         ".geojson" },
    { EXPORT_CSV,        ".csv" },
    { EXPORT_R,          ".txt" },
    { EXPORT_OCTAVE,     ".m" },
    { EXPORT_DAT,        ".dat" },
    { EXPORT_DTA,        ".dta" },
    { SAVE_OUTPUT,       ".txt" },
    { SAVE_TEX,          ".tex" },
    { SAVE_RTF,          ".rtf" },
    { SAVE_XML,          ".xml" },
    { SAVE_TEXT,         ".txt" },
    { SAVE_HELP_TEXT,    ".txt" },
    { OPEN_SCRIPT,       ".inp" },
    { OPEN_SESSION,      ".gretl" },
    { OPEN_RATS_DB,      ".rat" },
    { OPEN_PCGIVE_DB,    ".bn7" },
    { OPEN_GFN,          ".gfn" },
    { OPEN_SPEC,         ".spec" },
    { OPEN_BARS,         ".txt" },
    { SELECT_PDF,        ".pdf" },
    { FILE_OP_MAX,       NULL }
};

static int gdtb_save;
static char alt_gfndir[FILENAME_MAX];

#if GRETL_EDIT

static const char *get_extension_for_action (int action, gpointer data)
{
    int i;

    for (i=0; i < G_N_ELEMENTS(action_map); i++) {
	if (action == action_map[i].action) {
	    return action_map[i].ext;
	}
    }

    return NULL;
}

#else

static const char *get_extension_for_action (int action, gpointer data)
{
    const char *s = NULL;

    if (gdtb_save) {
	return ".gdtb";
    }

    if (action == SAVE_DATA || action == EXPORT_GDT ||
	action == SAVE_DATA_AS || action == SAVE_BOOT_DATA) {
	return ".gdt";
    } else if (action == EXPORT_GDTB) {
	return ".gdtb";
    } else if (action == SAVE_GNUPLOT || action == SAVE_GRAPHIC) {
	int ttype = gp_term_code(data, action);

	if (ttype == GP_TERM_EPS)      s = ".eps";
	else if (ttype == GP_TERM_PDF) s = ".pdf";
	else if (ttype == GP_TERM_FIG) s = ".fig";
	else if (ttype == GP_TERM_TEX) s = ".tex";
	else if (ttype == GP_TERM_PNG) s = ".png";
	else if (ttype == GP_TERM_EMF) s = ".emf";
	else if (ttype == GP_TERM_SVG) s = ".svg";
	else if (ttype == GP_TERM_HTM) s = ".html";
	else if (ttype == GP_TERM_PLT) s = ".plt";
    } else {
	int i;

	for (i=0; i < G_N_ELEMENTS(action_map); i++) {
	    if (action == action_map[i].action) {
		s = action_map[i].ext;
		break;
	    }
	}
    }

    return s;
}

#endif /* GRETL_EDIT or not */

static int fname_has_extension (const char *fname)
{
    const char *p;
    int ret = 0;

#ifdef G_OS_WIN32
    p = strrchr(fname, '\\');
    if (p == NULL) {
	p = strrchr(fname, '/');
    }
#else
    p = strrchr(fname, '/');
#endif

    if (p == NULL) {
	p = fname;
    }

    p = strrchr(p, '.');

    if (p != NULL) {
	/* got a '.' in the basename: we'll count fname
	   as having an extension unless the trailing
	   portion is an integer
	*/
	if (!integer_string(p+1)) {
	    ret = 1;
	}
    }

    return ret;
}

static int post_process_savename (char *fname, int action, gpointer data)
{
    int err = 0;

    if (fname == NULL) {
	err = 1;
    } else if (gretl_isdir(fname)) {
	/* fname is actually a directory */
	err = !SET_DIR_ACTION(action);
    } else {
	/* get the recommended extension for this action */
	const char *ext = get_extension_for_action(action, data);

	if (ext != NULL && !has_suffix(fname, ext)) {
	    /* We found a recommended extension, and it's not
	       already present on @fname; so we should probably
	       append the extension -- unless perhaps we're
	       looking at an action where there isn't exactly
	       a canonical extension and @fname already has some
	       sort of extension in place.
	    */
	    if (strcmp(ext, ".txt") && strcmp(ext, ".plt")) {
		/* append canonical extension */
		strcat(fname, ext);
	    } else if (!fname_has_extension(fname)) {
		strcat(fname, ext);
	    }
	}
    }

    return err;
}

static void script_window_update (windata_t *vwin,
				  const char *fname)
{
    gchar *basename;

    /* update internal filename record */
    strcpy(vwin->fname, fname);

    /* get basename for display */
    basename = g_path_get_basename(fname);

    if (window_is_tab(vwin)) {
	/* update the tab label */
	tabwin_tab_set_title(vwin, basename);
    } else {
	/* update the window title */
	const gchar *t0 = gtk_window_get_title(GTK_WINDOW(vwin->main));
	gchar *t1 = g_strdup_printf("gretl: %s", basename);

	if (strcmp(t0, t1)) {
	    gtk_window_set_title(GTK_WINDOW(vwin->main), t1);
	    window_list_revise_label(vwin->main, basename);
	}
	g_free(t1);
    }

    g_free(basename);

    if (vwin->role == EDIT_GP) {
	/* plot file no longer under session control */
	vwin->flags &= ~VWIN_SESSION_GRAPH;
	vwin->data = NULL;
    }

    mark_vwin_content_saved(vwin);
}

#ifndef GRETL_EDIT

static void handle_geoplot_save (const char *buf,
				 const char *fname,
				 windata_t *vwin)
{
    gchar *grfpath = session_graph_get_filename(vwin->data);
    gchar *tmpname, *datname;
    GError *gerr = NULL;
    gboolean ok;
    int err = 0;

    datname = g_strdup_printf("%s.dat", grfpath);
    tmpname = gretl_make_dotpath("geoplot_save.tmp");

    /* write @buf into a temporary file */
    ok = g_file_set_contents(tmpname, buf, -1, &gerr);

    if (!ok) {
	if (gerr != NULL) {
	    errbox(gerr->message);
	    g_error_free(gerr);
	} else {
	    gui_errmsg(E_FOPEN);
	}
    } else {
	err = transcribe_geoplot_file(tmpname, fname, datname);
	if (err) {
	    gui_errmsg(err);
	}
	gretl_remove(tmpname);
    }

    g_free(grfpath);
    g_free(tmpname);
    g_free(datname);
}

#endif /* not GRETL_EDIT */

static void
save_editable_content (int action, const char *fname, windata_t *vwin)
{
    gchar *buf;
    FILE *fp;

    if (editing_hansl(vwin->role)) {
        buf = textview_get_hansl(vwin->text, 1);
    } else {
        buf = textview_get_text(vwin->text);
    }

    if (buf == NULL) {
	errbox("Couldn't retrieve buffer");
	return;
    }

#ifndef GRETL_EDIT
    if (vwin->flags & VWIN_SESSION_GRAPH && vwin->data != NULL) {
	handle_geoplot_save(buf, fname, vwin);
	g_free(buf);
	return;
    }
#endif

    fp = gretl_fopen(fname, "w");
    if (fp == NULL) {
	file_write_errbox(fname);
	g_free(buf);
	return;
    }

    system_print_buf(buf, fp);
    g_free(buf);
    fclose(fp);

    if (action == SAVE_SCRIPT) {
	strcpy(scriptfile, fname);
#ifndef GRETL_EDIT
	mkfilelist(FILE_LIST_SCRIPT, scriptfile, 0);
#endif
	script_window_update(vwin, fname);
    } else if (action == SAVE_GP_CMDS ||
	       action == SAVE_R_CMDS ||
	       action == SAVE_OX_CMDS ||
	       action == SAVE_OCTAVE_CMDS ||
	       action == SAVE_PYTHON_CMDS ||
	       action == SAVE_STATA_CMDS ||
	       action == SAVE_JULIA_CODE ||
	       action == SAVE_DYNARE_CODE ||
	       action == SAVE_LPSOLVE_CODE ||
	       action == SAVE_X13_SPC) {
	script_window_update(vwin, fname);
    }
}

static void filesel_save_prn_buffer (PRN *prn, const char *fname)
{
    int fmt = prn_format(prn);
    FILE *fp;

    /* if @prn is in RTF, use binary mode to avoid messing
       with pre-set line endings
    */

    if ((fmt & GRETL_FORMAT_RTF) || fmt == GRETL_FORMAT_RTF_TXT) {
	fp = gretl_fopen(fname, "wb");
    } else {
	fp = gretl_fopen(fname, "w");
    }

    if (fp == NULL) {
	file_write_errbox(fname);
    } else {
	const char *buf = gretl_print_get_buffer(prn);

	fputs(buf, fp);
	fclose(fp);
    }
}

static void script_open_choice (const char *fname, windata_t *vwin,
				int role, int foreign)
{
    gchar *opts[] = {
	NULL, NULL,
    };
    char origname[MAXLEN];
    char newname[MAXLEN];
    int save = 0;
    int *cptr = NULL;
    gchar *chktext = NULL;
    int resp;

    gretl_basename(newname, fname, 0);
    if (strstr(vwin->fname, gretl_dotdir()) &&
	strstr(vwin->fname, "script_tmp")) {
	sprintf(origname, "[%s]", _("untitled"));
    } else {
	gretl_basename(origname, vwin->fname, 0);
    }
    opts[0] = g_strdup_printf(_("Open %s in new window"), newname);
    opts[1] = g_strdup_printf(_("Replace %s"), origname);

    if (vwin_content_changed(vwin)) {
	save = 1;
	cptr = &save;
	chktext = g_strdup_printf(_("Save changes to %s"), origname);
    }

    resp = radio_dialog_with_check(NULL, _("Open script"),
				   (const char **) opts, 2, 0, 0,
				   cptr, _(chktext),
				   vwin->main);

    if (resp == 0) {
	/* new window */
	if (foreign) {
	    view_script(fname, 1, role);
	} else {
	    if (view_script(fname, 1, EDIT_HANSL) != NULL) {
		strcpy(scriptfile, fname);
#ifndef GRETL_EDIT
		mkfilelist(FILE_LIST_SCRIPT, scriptfile, 0);
#endif
	    }
	}
    } else if (resp == 1) {
	/* replace */
	if (save) {
	    vwin_save_callback(NULL, vwin);
	}
	sourceview_insert_file(vwin, fname);
	vwin->role = role;
	mark_vwin_content_saved(vwin);
	vwin_set_filename(vwin, fname);
    }

    g_free(opts[0]);
    g_free(opts[1]);
    g_free(chktext);
}

static void filesel_open_script (const char *fname, windata_t *vwin)
{
    int role = EDIT_HANSL;
    int foreign = 0;

    if (vwin != NULL && vwin->role == EDIT_PKG_SAMPLE) {
	sourceview_insert_file(vwin, fname);
	return;
    }

    if (has_suffix(fname, ".R")) {
	role = EDIT_R;
    } else if (has_suffix(fname, ".plt")) {
	role = EDIT_GP;
    } else if (has_suffix(fname, ".ox")) {
	role = EDIT_OX;
    } else if (has_suffix(fname, ".m")) {
	role = EDIT_OCTAVE;
    } else if (has_suffix(fname, ".py")) {
	role = EDIT_PYTHON;
    } else if (has_suffix(fname, ".do") || has_suffix(fname, ".ado")) {
	role = EDIT_STATA;
    } else if (has_suffix(fname, ".jl")) {
	role = EDIT_JULIA;
    } else if (has_suffix(fname, ".mod")) {
	role = EDIT_DYNARE;
    } else if (has_suffix(fname, ".lp")) {
	role = EDIT_LPSOLVE;
    }

    if (role >= EDIT_GP && role < EDIT_MAX) {
	foreign = 1;
    }

    if (vwin != NULL && !window_is_tab(vwin)) {
	/* We're called from an existing (and untabbed) script
	   editor window: give choice of opening a new window
	   or replacing the current file.
	*/
	script_open_choice(fname, vwin, role, foreign);
    } else if (foreign) {
	view_script(fname, 1, role);
    } else if (view_script(fname, 1, EDIT_HANSL) != NULL) {
	strcpy(scriptfile, fname);
#ifndef GRETL_EDIT
	mkfilelist(FILE_LIST_SCRIPT, scriptfile, 0);
#endif
    }
}

#ifndef GRETL_EDIT

static void filesel_open_session (const char *fname)
{
    if (gretl_is_pkzip_file(fname)) {
	set_tryfile(fname);
	verify_open_session();
    } else {
	/* old script-style session file? */
	windata_t *vwin;

	if (has_system_prefix(fname, SCRIPT_SEARCH)) {
	    vwin = view_script(fname, 0, VIEW_SCRIPT);
	} else {
	    vwin = view_script(fname, 1, EDIT_HANSL);
	}

	if (vwin != NULL) {
	    strcpy(scriptfile, fname);
	}
    }
}

/* suggested_savename: pertains to the case where some data has
   been imported from a file with name @fname and the user has
   now chosen "Save data" (implicitly in native format). We'd
   like to offer a suggestion for the name of the file to save.
   A simple case would be, e.g., "foo" from imported "foo.xls".
   We'll not add a suffix because this will depend on the user's
   selection of gretl format, gdt vs gdtb.
*/

static char *suggested_savename (const char *fname)
{
    const char *ss = path_last_slash_const(fname);
    char *s, *sfx;

    if (ss == NULL) {
	s = g_strdup(fname);
    } else {
	s = g_strdup(ss + 1);
    }

    sfx = strrchr(s, '.');

    if (sfx != NULL) {
	/* if what follows '.' really looks like a suffix,
	   trim it off
	*/
	if (strlen(sfx) == 4 ||
	    !strcmp(sfx, ".xlsx") ||
	    !strcmp(sfx, ".gnumeric")) {
	    *sfx = '\0';
	}
    }

    return s;
}

static gchar *suggested_exportname (const char *fname, int action)
{
    const char *ss = path_last_slash_const(fname);
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
	    test = ".txt";
	    break;
	case EXPORT_CSV:
	    test = ".csv";
	    break;
	case EXPORT_DAT:
	case EXPORT_JM:
	    test = ".dat";
	    break;
	case EXPORT_DTA:
	    test = ".dta";
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

static gchar *suggested_mapname (const char *fname)
{
    if (has_suffix(fname, ".shp")) {
	gchar *targ = g_malloc0(strlen(fname) + 5);

	switch_ext(targ, fname, "geojson");
	return targ;
    } else {
	return NULL;
    }
}

static void bootstrap_save_callback (const char *fname)
{
    int err = bootstrap_save_data(fname);

    if (err) {
	gui_errmsg(err);
    }
}

static void os_open_other (const char *fname)
{
#if defined(G_OS_WIN32)
    win32_open_file(fname);
#elif defined(__APPLE__)
    osx_open_file(fname);
#else
    gretl_fork("xdg-open", fname, NULL);
#endif
}

#endif /* not GRETL_EDIT */

static void maybe_set_fsel_status (int action, FselDataSrc src,
				   gpointer p, int val)
{
    if (src == FSEL_DATA_STATUS && p != NULL) {
	* (int *) p = val;
    }
#ifndef GRETL_EDIT
    else if (action == OPEN_BARS && val) {
	/* error or cancel */
	set_plotbars_filename(NULL, p);
    }
#endif
}

static void
file_selector_process_result (const char *in_fname, int action,
			      FselDataSrc src, gpointer data)
{
    char fname[FILENAME_MAX];
    int err = 0;

    *fname = '\0';
    strncat(fname, in_fname, FILENAME_MAX - 1);

    if (action < END_OPEN) {
	/* opening a file: check that the file is accessible */
	err = gretl_test_fopen(fname, "r");
	if (err) {
	    file_read_errbox(fname);
	    maybe_set_fsel_status(action, src, data, E_FOPEN);
	    return;
	}
    }

#ifndef GRETL_EDIT
    if (action == OPEN_ANY) {
	/* designed for browsing an addon's "examples" directory */
	if (has_suffix(fname, ".inp")) {
	    action = OPEN_SCRIPT;
	} else if (has_suffix(fname, ".gdt") ||
		   has_suffix(fname, ".gdtb") ||
		   has_suffix(fname, ".csv") ||
		   has_suffix(fname, ".csv.gz") ||
		   has_suffix(fname, ".json") ||
		   has_suffix(fname, ".geojson") ||
		   has_suffix(fname, ".shp")) {
	    action = OPEN_DATA;
	} else if (strstr(fname, "README") || has_suffix(fname, ".txt")) {
	    view_file(fname, 0, 0, 78, 370, VIEW_FILE);
	    return;
	} else {
	    os_open_other(fname);
	    return;
	}
    }
#endif

    if (action == OPEN_SCRIPT) {
	if (src == FSEL_DATA_VWIN) {
	    filesel_open_script(fname, data);
	} else {
	    filesel_open_script(fname, NULL);
	}
    }
#ifndef GRETL_EDIT
    else if (action == OPEN_DATA) {
	set_tryfile(fname);
	verify_open_data(NULL, action, FALSE);
    } else if (action == APPEND_DATA) {
	set_tryfile(fname);
	do_open_data(NULL, action);
    } else if (action == OPEN_SCRIPT) {
	if (src == FSEL_DATA_VWIN) {
	    filesel_open_script(fname, data);
	} else {
	    filesel_open_script(fname, NULL);
	}
    } else if (action == OPEN_SESSION) {
	filesel_open_session(fname);
    } else if (action == OPEN_MARKERS) {
	do_add_markers(fname);
    } else if (action == OPEN_LABELS) {
	do_add_labels(fname);
    } else if (action == OPEN_BARS) {
	set_plotbars_filename(fname, data);
    } else if (action == OPEN_GFN) {
	edit_function_package(fname);
    } else if (action == OPEN_SPEC) {
	view_script(fname, 1, EDIT_SPEC);
    } else if (action == OPEN_RATS_DB) {
	open_rats_window(fname);
    } else if (action == OPEN_PCGIVE_DB) {
	open_bn7_window(fname);
    } else if (action == UPLOAD_PKG) {
	upload_specified_package(fname);
    } else if (action == INSTALL_PKG) {
	do_local_pkg_install(fname);
    }
#endif /* not GRETL_EDIT */

    if (action < END_OPEN) {
	return;
    }

    /* now for the save/export options */

    if (post_process_savename(fname, action, data)) {
	return;
    }

#ifdef GRETL_EDIT
    if (src == FSEL_DATA_PRN) {
	filesel_save_prn_buffer((PRN *) data, fname);
    } else if (action == SET_PROG ||
	       action == SET_DIR ||
	       action == SET_OTHER) {
	set_path_callback(data, fname);
    } else {
	windata_t *vwin = (windata_t *) data;

	save_editable_content(action, fname, vwin);
    }
#else
    if (action == SAVE_TEX) {
	if (src == FSEL_DATA_PRN) {
	    save_latex((PRN *) data, fname);
	} else {
	    save_latex(NULL, fname);
	}
    } else if (src == FSEL_DATA_PRN) {
	filesel_save_prn_buffer((PRN *) data, fname);
    } else if (SAVE_DATA_ACTION(action)) {
	err = do_store(fname, action, data);
    } else if (action == SAVE_GNUPLOT) {
	save_graph_to_file(data, fname);
    } else if (action == SAVE_GRAPHIC) {
	save_graphic_to_file(data, fname);
    } else if (action == SAVE_SESSION) {
	err = save_session(fname);
    } else if (action == SAVE_CMD_LOG) {
	err = save_session_commands(fname);
    } else if (action == SAVE_FUNCTIONS) {
	err = save_function_package(fname, data);
    } else if (action == SAVE_FUNCTIONS_AS) {
	err = save_function_package_as_script(fname, data);
    } else if (action == SAVE_GFN_SPEC) {
	err = save_function_package_spec(fname, data);
    } else if (action == SAVE_GFN_ZIP) {
	err = save_function_package_zipfile(fname, data);
    } else if (action == SAVE_MAP || action == WRITE_MAP) {
	err = do_store(fname, action, data);
    } else if (action == SELECT_PDF) {
	err = set_package_pdfname(fname, data);
    } else if (action == SAVE_BOOT_DATA) {
	bootstrap_save_callback(fname);
    } else if (action == SAVE_MARKERS) {
	err = do_save_markers(fname);
    } else if (action == SAVE_LABELS) {
	err = do_save_labels(fname);
    } else if (action == SET_PROG ||
	       action == SET_DIR ||
	       action == SET_OTHER) {
	set_path_callback(data, fname);
    } else if (action == SET_WDIR) {
	set_working_dir_callback(data, fname);
    } else if (action == SET_FDIR) {
	set_alternate_gfn_dir(data, fname);
    } else if (action == SET_DBDIR) {
	set_db_dir_callback(data, fname);
    } else {
	windata_t *vwin = (windata_t *) data;

	save_editable_content(action, fname, vwin);
    }
#endif /* GRETL_EDIT or not */

    maybe_set_fsel_status(action, src, data, err);
}

static char *get_filter_suffix (int action, gpointer data, char *suffix)
{

    const char *ext = get_extension_for_action(action, data);

    if (ext == NULL) {
	strcpy(suffix, "*");
    } else {
	sprintf(suffix, "*%s", ext);
    }

    return suffix;
}

static void maybe_upcase_filter_pattern (GtkFileFilter *filter,
					 const char *s)
{
    char *p, tmp[16];
    int changed = 0;

    strcpy(tmp, s);
    p = tmp + 1;

    while (*p) {
	if (islower((unsigned char) *p)) {
	    *p = toupper(*p);
	    changed = 1;
	}
	p++;
    }

    if (changed) {
	gtk_file_filter_add_pattern(filter, tmp);
    }
}

static GtkFileFilter *get_file_filter (int action, gpointer data)
{
    GtkFileFilter *filter;
    char suffix[16];

    filter = gtk_file_filter_new();
    get_filter_suffix(action, data, suffix);
    gtk_file_filter_add_pattern(filter, suffix);

    /* support Windows stupidity */
    maybe_upcase_filter_pattern(filter, suffix);

    return filter;
}

static void set_default_progs_path (GtkFileChooser *fsel)
{
#ifdef G_OS_WIN32
    char *progs = program_files_path();

    if (progs != NULL) {
	gtk_file_chooser_set_current_folder(fsel, progs);
	free(progs);
    }
#else
    gtk_file_chooser_set_current_folder(fsel, "/usr/bin");
#endif
}

static void set_default_other_path (GtkFileChooser *fsel)
{
    const gchar *home = g_get_home_dir();

    if (home != NULL) {
	gtk_file_chooser_set_current_folder(fsel, home);
    }
}

static void filesel_maybe_set_current_name (GtkFileChooser *filesel,
					    int action, FselDataSrc src,
					    gpointer data)
{
    gchar *currname;

#ifdef GRETL_EDIT
    if (action == SET_PROG || action == SET_DIR || action == SET_OTHER) {
	currname = (gchar *) data;
	if (currname != NULL && g_path_is_absolute(currname)) {
	    gtk_file_chooser_set_filename(filesel, currname);
	} else if (action == SET_PROG) {
	    set_default_progs_path(filesel);
	} else if (action == SET_OTHER) {
	    set_default_other_path(filesel);
	}
    }
#else
    if (action == SAVE_DATA && *datafile != '\0') {
	currname = suggested_savename(datafile);
	gtk_file_chooser_set_current_name(filesel, currname);
	g_free(currname);
    } else if (EXPORT_OTHER(action) && src != FSEL_DATA_PRN &&
	       *datafile != '\0') {
	/* formulate export name based on current datafile */
	currname = suggested_exportname(datafile, action);
	gtk_file_chooser_set_current_name(filesel, currname);
	g_free(currname);
    } else if (action == SAVE_MAP && *datafile != '\0') {
	currname = suggested_mapname(datafile);
	if (currname != NULL) {
	    gtk_file_chooser_set_current_name(filesel, currname);
	    g_free(currname);
	}
    } else if (action == SAVE_CMD_LOG) {
	gtk_file_chooser_set_current_name(filesel, "session.inp");
    } else if (action == SET_PROG || action == SET_DIR || action == SET_OTHER) {
	currname = (gchar *) data;
	if (currname != NULL && g_path_is_absolute(currname)) {
	    gtk_file_chooser_set_filename(filesel, currname);
	} else if (action == SET_PROG) {
	    set_default_progs_path(filesel);
	} else if (action == SET_OTHER) {
	    set_default_other_path(filesel);
	}
    } else if (action == SAVE_FUNCTIONS ||
	       action == SAVE_FUNCTIONS_AS ||
	       action == SAVE_GFN_SPEC ||
	       action == SAVE_GFN_ZIP) {
	char fname[MAXLEN];

	*fname = '\0';
	get_default_package_name(fname, data, action);
	if (*fname != '\0') {
	    gtk_file_chooser_set_current_name(filesel, fname);
	}
    } else if (action == SELECT_PDF) {
	char fname[MAXLEN];

	*fname = '\0';
	get_default_package_name(fname, data, action);
	if (*fname != '\0') {
	    gtk_file_chooser_set_filename(filesel, fname);
	}
    } else if (action == SAVE_MARKERS) {
	gtk_file_chooser_set_current_name(filesel, "markers.txt");
    } else if (action == SAVE_LABELS) {
	gtk_file_chooser_set_current_name(filesel, "labels.txt");
    }
#endif /* GRETL_EDIT or not */
}

/* local extension of enum GretlFileType in dataio.h */
enum {
    GRETL_TXT = GRETL_UNRECOGNIZED + 1,
    GRETL_SHP,
    GRETL_GEOJSON,
    GRETL_ALL
};

struct filter_info {
    const char *desc;
    const char *pat;
    int id;
};

static struct filter_info data_filters[] = {
    { N_("Gretl datafiles (*.gdt, *.gdtb)"), "*.gdt", GRETL_XML_DATA },
    { N_("CSV files (*.csv)"), "*.csv", GRETL_CSV },
    { N_("ASCII files (*.txt)"), "*.txt",  GRETL_TXT },
    { N_("Gnumeric files (*.gnumeric)"), "*.gnumeric", GRETL_GNUMERIC },
    { N_("Open Document files (*.ods)"), "*.ods", GRETL_ODS },
    { N_("Excel files (*.xls)"), "*.xls", GRETL_XLS },
    { N_("Excel files (*.xlsx)"), "*.xlsx", GRETL_XLSX },
    { N_("Stata files (*.dta)"), "*.dta", GRETL_DTA },
    { N_("Eviews files (*.wf1)"), "*.wf1", GRETL_WF1 },
    { N_("SPSS files (*.sav)"), "*.sav", GRETL_SAV },
    { N_("SAS xport files (*.xpt)"), "*.xpt", GRETL_SAS },
    { N_("JMulTi files (*.dat)"), "*.dat", GRETL_JMULTI },
    { N_("Shapefiles (*.shp)"), "*.shp", GRETL_SHP },
    { N_("GeoJSON files (*.geojson, *.json)"), "*.geojson", GRETL_GEOJSON },
    { N_("all files (*.*)"), "*" , GRETL_ALL }
};

static struct filter_info script_filters[] = {
    { N_("gretl script files (*.inp)"), "*.inp", EDIT_HANSL },
    { N_("GNU R files (*.R)"), "*.R", EDIT_R },
    { N_("gnuplot files (*.plt)"), "*.plt", EDIT_GP },
    { N_("GNU Octave files (*.m)"), "*.m", EDIT_OCTAVE },
    { N_("Python files (*.py)"), "*.py", EDIT_PYTHON },
    { N_("Julia files (*.jl)"), "*.jl", EDIT_JULIA },
    { N_("Ox files (*.ox)"), "*.ox", EDIT_OX },
    { N_("Stata files (*.do, *.ado)"), "*.do", EDIT_STATA },
    { N_("Dynare files (*.mod)"), "*.mod", EDIT_DYNARE },
    { N_("lpsolve files (*.lp)"), "*.lp", EDIT_LPSOLVE }
};

static int n_data_filters = G_N_ELEMENTS(data_filters);
static int n_script_filters = G_N_ELEMENTS(script_filters);
static int open_filter_index;

static GtkFileFilter *add_filter_by_index (GtkWidget *filesel,
					   int action, int i)
{
    GtkFileFilter *filt = gtk_file_filter_new();
    struct filter_info *fi;

    if (action == OPEN_SCRIPT) {
	fi = &script_filters[i];
    } else {
	fi = &data_filters[i];
    }

    gtk_file_filter_set_name(filt, _(fi->desc));
    gtk_file_filter_add_pattern(filt, fi->pat);
    if (action == OPEN_DATA) {
	if (fi->id == GRETL_XML_DATA) {
	    gtk_file_filter_add_pattern(filt, "*.gdtb");
	} else if (fi->id == GRETL_GEOJSON) {
	    gtk_file_filter_add_pattern(filt, "*.json");
	} else if (fi->id == GRETL_CSV) {
	    gtk_file_filter_add_pattern(filt, "*.csv.gz");
	} else {
	    maybe_upcase_filter_pattern(filt, fi->pat);
	}
	g_object_set_data(G_OBJECT(filt), "filter-index",
			  GINT_TO_POINTER(i));
    } else if (action == OPEN_SCRIPT) {
	if (fi->id == EDIT_STATA) {
	    gtk_file_filter_add_pattern(filt, "*.ado");
	}
    }

    gtk_file_chooser_add_filter(GTK_FILE_CHOOSER(filesel), filt);

    return filt;
}

static void filesel_add_data_filters (GtkWidget *filesel)
{
    GtkFileFilter *flt;
    int i;

    for (i=0; i<n_data_filters; i++) {
	flt = add_filter_by_index(filesel, OPEN_DATA, i);
	if (i == open_filter_index) {
	    gtk_file_chooser_set_filter(GTK_FILE_CHOOSER(filesel), flt);
	}
    }
}

static void filesel_add_script_filters (GtkWidget *filesel,
					int role)
{
    GtkFileFilter *flt;
    int i;

    for (i=0; i<n_script_filters; i++) {
	flt = add_filter_by_index(filesel, OPEN_SCRIPT, i);
	if (script_filters[i].id == role) {
	    gtk_file_chooser_set_filter(GTK_FILE_CHOOSER(filesel), flt);
	}
    }
}

static GtkFileFilter *filesel_add_data_filter (GtkWidget *filesel,
					       int id)
{
    int i;

    for (i=0; i<n_data_filters; i++) {
	if (data_filters[i].id == id) {
	    return add_filter_by_index(filesel, OPEN_DATA, i);
	}
    }
    return NULL;
}

static GtkFileFilter *filesel_add_script_filter (GtkWidget *filesel,
						 int id)
{
    int i;

    for (i=0; i<n_script_filters; i++) {
	if (script_filters[i].id == id) {
	    return add_filter_by_index(filesel, OPEN_DATA, i);
	}
    }
    return NULL;
}

static GtkFileFilter *filesel_add_filter (GtkWidget *filesel,
					  const char *desc,
					  const char *pat)
{
    GtkFileFilter *filt = gtk_file_filter_new();

    gtk_file_filter_set_name(filt, _(desc));
    gtk_file_filter_add_pattern(filt, pat);
    maybe_upcase_filter_pattern(filt, pat);
    gtk_file_chooser_add_filter(GTK_FILE_CHOOSER(filesel), filt);

    return filt;
}

static void filesel_set_filter_patterns (GtkWidget *filesel,
					 const char *pat1,
					 const char *pat2)
{
    GtkFileFilter *filt = gtk_file_filter_new();

    gtk_file_filter_add_pattern(filt, pat1);
    gtk_file_filter_add_pattern(filt, pat2);
    gtk_file_chooser_set_filter(GTK_FILE_CHOOSER(filesel), filt);
}

static void save_filter_changed (GtkComboBox *cb, GtkWidget *w)
{
    GtkFileChooser *fc = GTK_FILE_CHOOSER(w);
    GtkFileFilter *ff = gtk_file_chooser_get_filter(fc);

    if (ff != NULL) {
	GtkWidget *ew = gtk_file_chooser_get_extra_widget(fc);
	const char *s = gtk_file_filter_get_name(ff);

	if (ew != NULL && s != NULL) {
	    gtk_widget_set_sensitive(ew, strstr(s, ".gdtb") == NULL);
	}
    }
}

void find_filter_combo (GtkWidget *w, gpointer p)
{
    GtkWidget **pw = p;

    if (GTK_IS_COMBO_BOX(w)) {
	*pw = w;
    } else if (*pw == NULL && GTK_IS_CONTAINER(w)) {
	gtk_container_foreach(GTK_CONTAINER(w), find_filter_combo, pw);
    }
}

static void conditionalize_compression (GtkWidget *filesel)
{
    GtkWidget *ca = gtk_dialog_get_content_area(GTK_DIALOG(filesel));
    GList *L = gtk_container_get_children(GTK_CONTAINER(ca));

    if (GTK_IS_BOX(L->data)) {
	GtkWidget *combo = NULL;

	gtk_container_foreach(GTK_CONTAINER(L->data),
			      find_filter_combo, &combo);
	if (combo != NULL) {
	    g_signal_connect(G_OBJECT(combo), "changed",
			     G_CALLBACK(save_filter_changed),
			     filesel);
	}
    }

    g_list_free(L);
}

/* return non-zero if we add more than one selectable filter */

static int filesel_set_filters (GtkWidget *filesel, int action,
				FselDataSrc src, gpointer data)
{
    int multi = 1, role = 0;

    if (src == FSEL_DATA_VWIN) {
	windata_t *vwin = (windata_t *) data;

	role = vwin->role;
    }

    if (role == EDIT_PKG_SAMPLE) {
	GtkFileFilter *filter = get_file_filter(action, data);

	gtk_file_chooser_set_filter(GTK_FILE_CHOOSER(filesel), filter);
	multi = 0;
    } else if (action == OPEN_DATA || action == APPEND_DATA) {
	filesel_add_data_filters(filesel);
    } else if (action == OPEN_SCRIPT) {
	filesel_add_script_filters(filesel, role);
    } else if (action == OPEN_LABELS || action == OPEN_BARS) {
	filesel_add_data_filter(filesel, GRETL_TXT);
	filesel_add_data_filter(filesel, GRETL_ALL);
    } else if (action == SAVE_DATA || action == SAVE_DATA_AS ||
	       action == SAVE_BOOT_DATA) {
	filesel_add_filter(filesel, N_("Gretl datafiles (*.gdt)"), "*.gdt");
	filesel_add_filter(filesel, N_("Gretl binary datafiles (*.gdtb)"), "*.gdtb");
    } else if (action == OPEN_ANY) {
	filesel_add_data_filter(filesel, GRETL_ALL);
	filesel_add_script_filter(filesel, EDIT_HANSL);
	filesel_add_data_filter(filesel, GRETL_XML_DATA);
	filesel_add_data_filter(filesel, GRETL_SHP);
	filesel_add_data_filter(filesel, GRETL_GEOJSON);
    } else if (action == UPLOAD_PKG || action == INSTALL_PKG) {
	/* could add OPEN_GFN here?? */
	filesel_set_filter_patterns(filesel, "*.gfn", "*.zip");
    } else {
	GtkFileFilter *filter = get_file_filter(action, data);

	gtk_file_chooser_set_filter(GTK_FILE_CHOOSER(filesel), filter);
	multi = 0;
    }

    return multi;
}

static void remember_folder (GtkFileChooser *chooser, char *savedir)
{
    gchar *folder = gtk_file_chooser_get_current_folder(chooser);

    if (folder != NULL) {
	strcpy(savedir, folder);
	g_free(folder);
    }
}

#ifndef GRETL_EDIT

static void record_compress_level (GtkWidget *b, gpointer p)
{
    int level;

    level = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(b));
    set_optval_int(STORE, OPT_Z, level); /* data compression */
}

static void add_compression_level_option (GtkWidget *filesel)
{
    GtkWidget *hbox, *label, *spin;
    guint64 dsetsize;
    int deflt = 1;

    dsetsize = (dataset->v - 1) * dataset->n * 8;
    if (dsetsize < 30 * 1024) {
	/* default to no compression for small data */
	deflt = 0;
    }

    hbox = gtk_hbox_new(FALSE, 5);
    label = gtk_label_new(_("Compression level (0 = none)"));
    spin = gtk_spin_button_new_with_range(0, 9, 1);

    gtk_spin_button_set_value(GTK_SPIN_BUTTON(spin), deflt);
    gtk_spin_button_set_digits(GTK_SPIN_BUTTON(spin), 0);
    g_signal_connect(G_OBJECT(spin), "destroy",
		     G_CALLBACK(record_compress_level), NULL);

    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), spin, FALSE, FALSE, 5);
    gtk_widget_show_all(hbox);
    gtk_file_chooser_set_extra_widget(GTK_FILE_CHOOSER(filesel),
				      hbox);
}

#endif /* not GRETL_EDIT */

static void check_native_data_save_filter (GtkWidget *filesel)
{
    GtkFileFilter *filter;
    const gchar *filter_name;

    /* we need to know if .gdt or .gdtb was chosen */

    filter = gtk_file_chooser_get_filter(GTK_FILE_CHOOSER(filesel));
    filter_name = gtk_file_filter_get_name(filter);
    if (strstr(filter_name, "gdtb") != NULL) {
	gdtb_save = 1;
    }
}

static void peek_data_open_filter (GtkWidget *filesel)
{
    GtkFileFilter *filter;

    filter = gtk_file_chooser_get_filter(GTK_FILE_CHOOSER(filesel));
    open_filter_index = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(filter),
							  "filter-index"));
}

static void set_up_filesel_filters (GtkWidget *filesel,
				    int action,
				    FselDataSrc src,
				    gpointer data,
				    gboolean plain_open)
{
    int multi = filesel_set_filters(filesel, action, src, data);

    if (plain_open && !multi) {
	/* For an "Open" action with a single filename filter:
	   make the dialog's title bar more informative by
	   showing the required extension. (When multiple
	   filters are offered the user will already have
	   plenty of information available.)
	*/
	const char *s = get_extension_for_action(action, data);

	if (s != NULL) {
	    gchar *title = g_strdup_printf(_("Open %s file"), s);

	    gtk_window_set_title(GTK_WINDOW(filesel), title);
	    g_free(title);
	}
    }
}

static void gtk_file_selector (int action, FselDataSrc src,
			       gpointer data, GtkWidget *parent,
			       const char *dirname)
{
    static char savedir[MAXLEN];
    char startdir[MAXLEN];
    GtkWidget *filesel;
    GtkFileChooserAction fsel_action;
    gchar *title = NULL;
    const gchar *okstr;
    int remember = get_keep_folder();
    gboolean plain_open = 0;
    gint response;

    gdtb_save = 0;

    if (SET_DIR_ACTION(action)) {
	fsel_action = GTK_FILE_CHOOSER_ACTION_SELECT_FOLDER;
	okstr = GTK_STOCK_OK;
	title = g_strdup_printf("gretl: %s", _("select directory"));
	remember = 0;
    } else if (action == SET_PROG ||
	       action == SET_OTHER ||
	       action == UPLOAD_PKG ||
	       action == SELECT_PDF) {
	fsel_action = GTK_FILE_CHOOSER_ACTION_OPEN;
	okstr = GTK_STOCK_OK;
	if (action == SELECT_PDF) {
	    title = g_strdup_printf("gretl: %s", _("select PDF file"));
	} else {
	    title = g_strdup_printf("gretl: %s", _("select file"));
	}
	remember = 0;
    } else if (action == INSTALL_PKG) {
	fsel_action = GTK_FILE_CHOOSER_ACTION_OPEN;
	okstr = _("Select package"); /* ?? */
	title = g_strdup_printf("gretl: %s", _("install function package"));
    } else if (action < END_OPEN) {
	fsel_action = GTK_FILE_CHOOSER_ACTION_OPEN;
	okstr = GTK_STOCK_OPEN;
	title = g_strdup_printf("gretl: %s", _("open file"));
	plain_open = 1;
    } else {
	/* it's a save action of some sort */
	fsel_action = GTK_FILE_CHOOSER_ACTION_SAVE;
	okstr = GTK_STOCK_SAVE;
	if (action == SAVE_FUNCTIONS ||
	    action == SAVE_DATA_PKG ||
	    action == SAVE_GFN_SPEC ||
	    action == SAVE_GFN_ZIP) {
	    remember = 0;
	}
	if (action == SAVE_MAP) {
	    title = g_strdup_printf("gretl: %s", _("save map as geojson"));
	} else if (action == WRITE_MAP) {
	    title = g_strdup_printf("gretl: %s", _("write map as geojson"));
	} else {
	    title = g_strdup_printf("gretl: %s", _("save file"));
	}
    }

#ifdef GRETL_EDIT
    if (dirname != NULL) {
	strcpy(startdir, dirname);
    } else if (remember && *savedir != '\0') {
	strcpy(startdir, savedir);
    } else {
	get_default_dir_for_action(startdir, action);
    }
#else
    if (dirname != NULL) {
	strcpy(startdir, dirname);
    } else if (remember && *savedir != '\0') {
	strcpy(startdir, savedir);
    } else if ((action == UPLOAD_PKG || action == SET_FDIR) &&
	       *alt_gfndir != '\0') {
	strcpy(startdir, alt_gfndir);
    } else if (action == SAVE_FUNCTIONS) {
	/* we come here only if the user has not chosen
	   to "install" a newly saved package */
	get_default_dir_for_action(startdir, 0);
    } else if (action == SAVE_GFN_SPEC ||
	       action == SAVE_GFN_ZIP) {
	get_gfn_dir(startdir, data);
    } else if (action == SELECT_PDF) {
	get_gfn_pdf_dir(startdir, data);
    } else {
	get_default_dir_for_action(startdir, action);
    }
#endif

    if (parent == NULL) {
	/* by default the file selector is parented by the
	   gretl main window */
	parent = mdata->main;
    }

    filesel = gtk_file_chooser_dialog_new(title,
					  GTK_WINDOW(parent),
					  fsel_action,
					  GTK_STOCK_CANCEL,
					  GTK_RESPONSE_CANCEL,
					  okstr,
					  GTK_RESPONSE_ACCEPT,
					  NULL);
#ifndef GRETL_EDIT
    if (action == SAVE_DATA ||
	action == SAVE_DATA_AS ||
	action == SAVE_BOOT_DATA ||
	action == EXPORT_GDT) {
	add_compression_level_option(filesel);
    }
#endif

    if (fsel_action == GTK_FILE_CHOOSER_ACTION_SAVE) {
	gtk_file_chooser_set_do_overwrite_confirmation(GTK_FILE_CHOOSER(filesel),
						       TRUE);
    }

    gtk_dialog_set_default_response(GTK_DIALOG(filesel), GTK_RESPONSE_ACCEPT);

    if (SET_DIR_ACTION(action)) {
	gtk_file_chooser_set_filename(GTK_FILE_CHOOSER(filesel), startdir);
    } else {
	if (action != SET_PROG && action != SET_OTHER) {
	    set_up_filesel_filters(filesel, action, src, data, plain_open);
	}
	gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(filesel), startdir);
	filesel_maybe_set_current_name(GTK_FILE_CHOOSER(filesel), action,
				       src, data);
    }

    if (action == SAVE_DATA || action == SAVE_DATA_AS || action == SAVE_BOOT_DATA) {
	conditionalize_compression(filesel);
    }

    response = gtk_dialog_run(GTK_DIALOG(filesel));

    g_free(title);

    if (response == GTK_RESPONSE_ACCEPT) {
	gchar *fname;

	if (action == SAVE_DATA || action == SAVE_DATA_AS) {
	    check_native_data_save_filter(filesel);
	} else if (action == OPEN_DATA) {
	    /* we'll want to record this choice */
	    peek_data_open_filter(filesel);
	}

	fname = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(filesel));

	if (remember) {
	    remember_folder(GTK_FILE_CHOOSER(filesel), savedir);
	} else if (action == UPLOAD_PKG || action == SET_FDIR) {
	    strcpy(alt_gfndir, fname);
	}

	gtk_widget_destroy(filesel);
	filesel = NULL;

	if (fname != NULL) {
	    file_selector_process_result(fname, action, src, data);
	    g_free(fname);
	}
    } else {
	maybe_set_fsel_status(action, src, data, GRETL_CANCEL);
    }

    if (filesel != NULL) {
	gtk_widget_destroy(filesel);
    }
}

void file_selector (int action,
		    FselDataSrc src,
		    gpointer data)
{
    GtkWidget *w = NULL;

    if (src == FSEL_DATA_VWIN) {
	windata_t *vwin = (windata_t *) data;

	w = vwin_toplevel(vwin);
    }

    gtk_file_selector(action, src, data, w, NULL);
}

void file_selector_with_parent (int action,
				FselDataSrc src,
				gpointer data,
				GtkWidget *parent)
{
    gtk_file_selector(action, src, data, parent, NULL);
}

void file_selector_with_startdir (int action,
				  const char *startdir,
				  GtkWidget *parent)
{
    gtk_file_selector(action, 0, NULL, parent, startdir);
}
