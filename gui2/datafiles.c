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

/* datafiles.c : for gretl */

#define COLL_DEBUG 0

#include "gretl.h"
#include "datafiles.h"
#include "database.h"
#include "filelists.h"
#include "gretl_www.h"
#include "menustate.h"
#include "fnsave.h"
#include "fncall.h"
#include "treeutils.h"
#include "selector.h"
#include "toolbar.h"
#include "winstack.h"

#include "gretl_xml.h"
#include "gretl_func.h"

#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include <errno.h>
#include <time.h>

#define N_BROWSER_TYPES 5
#define BROWSER_BUSY    1
#define BROWSER_OK      0

static GtkWidget *browsers[N_BROWSER_TYPES];

static GtkWidget *files_window (windata_t *fdata);
static GtkWidget *files_notebook (windata_t *fdata, int code);
static int populate_notebook_filelists (windata_t *fdata, 
					GtkWidget *notebook,
					int code);

typedef struct _file_collection file_collection;

struct _file_collection {
    char *path;
    char *descfile;
    char *title;
    int which;
    GtkWidget *page;
};

enum {
    STACK_PUSH,
    STACK_POP_DATA,
    STACK_POP_PS,
    STACK_RESET_DATA,
    STACK_RESET_PS,
    STACK_SORT_DATA,
    STACK_SORT_PS,
    STACK_DESTROY
};

enum {
    RECOG_OK,
    RECOG_ERROR,
    RECOG_NOT
};

enum {
    COLL_NONE,
    COLL_DATA,
    COLL_PS
};

struct fpkg_response {
    int col1_width;
    int try_server;
};

#define REMOTE_ACTION(c) (c == REMOTE_DB || c == REMOTE_FUNC_FILES)

static void 
fpkg_response_init (struct fpkg_response *f, gpointer p)
{
    f->col1_width = 0;

    if (p == NULL || p == mdata) {
	/* called from main menu */
	f->try_server = 0;
    } else {
	f->try_server = -1;
    }
}

static char *full_path (char *s1, const char *s2)
{
    static char fpath[FILENAME_MAX];
    int n = strlen(s1);

    if (s1[n-1] == '.') {
	s1[n-1] = '\0';
	n--;
    }
    
    if (s1[n-1] == SLASH) {
	sprintf(fpath, "%s%s", s1, s2);
    } else {
	sprintf(fpath, "%s%c%s", s1, SLASH, s2);
    }

#if COLL_DEBUG
    fprintf(stderr, "full_path: got '%s' from '%s' + '%s'\n",
	    fpath, s1, s2);
#endif

    return fpath;
}

static int recognized_collection (file_collection *coll)
{
    const file_collection recognized_data[] = {
	{ "wooldridge", "jw_descriptions", "Wooldridge", COLL_DATA, NULL },
	{ "gujarati", "dg_descriptions", "Gujarati", COLL_DATA, NULL },
	{ "pwt56", "descriptions", "PWT 56", COLL_DATA, NULL },
	{ NULL, NULL, NULL, COLL_DATA, NULL }
    }; 
    const file_collection recognized_ps[] = {
	{ "pwt56", "ps_descriptions", "PWT 56", COLL_PS, NULL },
	{ NULL, NULL, NULL, COLL_PS, NULL }
    }; 
    int i;

    for (i=0; recognized_data[i].path != NULL; i++) {
	if (strstr(coll->path, recognized_data[i].path) &&
	    !strcmp(coll->descfile, recognized_data[i].descfile)) {
	    coll->title = malloc(strlen(recognized_data[i].title) + 1);
	    if (coll->title == NULL) {
		return RECOG_ERROR;
	    }
	    strcpy(coll->title, recognized_data[i].title);
	    coll->which = COLL_DATA;
	    return RECOG_OK;
	}
    }

    for (i=0; recognized_ps[i].path != NULL; i++) {
	if (strstr(coll->path, recognized_ps[i].path) &&
	    !strcmp(coll->descfile, recognized_ps[i].descfile)) {
	    coll->title = malloc(strlen(recognized_ps[i].title) + 1);
	    if (coll->title == NULL) {
		return RECOG_ERROR;
	    }
	    strcpy(coll->title, recognized_ps[i].title);
	    coll->which = COLL_PS;
	    return RECOG_OK;
	}
    }

    return RECOG_NOT;
} 

static int get_title_from_descfile (file_collection *coll)
{
    FILE *fp;
    char *test;
    char line[64], title[24];
    int err = 0;

    test = full_path(coll->path, coll->descfile);

    fp = gretl_fopen(test, "r");

    if (fp == NULL) {
	err = 1;
    } else if (fgets(line, sizeof line, fp) == NULL) {
	err = 1;
    } else {
	chopstr(line);

	if (sscanf(line, "# %23[^:]", title) != 1) {
	    err = 1;
	} else {
	    coll->title = g_strdup(title);
	    if (coll->title == NULL) {
		err = 1;
	    }
	}
    } 

    if (fp != NULL) {
	fclose(fp);
    }

    return err;
}

static void free_file_collection (file_collection *coll)
{
    free(coll->path);
    free(coll->descfile);
    if (coll->title != NULL) {
	free(coll->title);
    }
    free(coll);
}  

static file_collection *file_collection_new (const char *path,
					     const char *descfile,
					     int *err)
{
    file_collection *coll;
    int chk;

    *err = 0;

    coll = malloc(sizeof *coll);
    if (coll == NULL) {
	*err = 1;
	return NULL;
    }

    coll->title = NULL;
    
    coll->path = gretl_strdup(path);
    if (coll->path == NULL) {
	free(coll);
	*err = 1;
	return NULL;
    }

    coll->descfile = gretl_strdup(descfile);
    if (coll->descfile == NULL) {
	free(coll->path);
	free(coll);
	*err = 1;
	return NULL;
    }  

    if (strstr(coll->descfile, "ps_")) {
	coll->which = COLL_PS;
    } else {
	coll->which = COLL_DATA;
    }

    chk = recognized_collection(coll);

    if (chk == RECOG_ERROR) {
	free_file_collection(coll);
	*err = 1;
	coll = NULL;
    } else if (chk == RECOG_NOT) {
	chk = get_title_from_descfile(coll);
	if (chk) {
	    free_file_collection(coll);
	    coll = NULL;
	}
    }
    
    return coll;
}

static int compare_colls (const void *a, const void *b)
{
    const file_collection *ca = *(const file_collection **) a;
    const file_collection *cb = *(const file_collection **) b;
     
    return strcmp(ca->title, cb->title);
}

static void collection_stack_sort (file_collection **colls, int n)
{
    file_collection *tmp;
    int i;

    for (i=0; i<n; i++) {
	if (!strcmp(colls[i]->title, "Gretl")) {
	    if (i > 0) {
		tmp = colls[0];
		colls[0] = colls[i];
		colls[i] = tmp;
	    }
	    break;
	}
    }

    qsort(colls + 1, n - 1, sizeof *colls, compare_colls);
}

static file_collection *collection_stack (file_collection *coll, int op)
{
    static file_collection **datacoll;
    static file_collection **pscoll;
    static int n_data;
    static int n_data_popped;
    static int n_ps;
    static int n_ps_popped;
    file_collection *ret = NULL;
    int j;

    if (op == STACK_PUSH && coll != NULL) {
	if (coll->which == COLL_DATA) {
	    datacoll = realloc(datacoll, (n_data + 1) * sizeof *datacoll);
	    if (datacoll != NULL) {
		datacoll[n_data++] = coll;
		ret = coll;
	    }
	} else if (coll->which == COLL_PS) {
	    pscoll = realloc(pscoll, (n_ps + 1) * sizeof *pscoll);
	    if (pscoll != NULL) {
		pscoll[n_ps++] = coll;
		ret = coll;
	    }
	}
    } else if (op == STACK_POP_DATA && n_data_popped < n_data) {
	ret = datacoll[n_data_popped++];
    } else if (op == STACK_POP_PS && n_ps_popped < n_ps) {
	ret = pscoll[n_ps_popped++];
    } else if (op == STACK_RESET_DATA) {
	n_data_popped = 0;
    } else if (op == STACK_RESET_PS) {
	n_ps_popped = 0;
    } else if (op == STACK_SORT_DATA) {
	collection_stack_sort(datacoll, n_data);
    } else if (op == STACK_SORT_PS) {
	collection_stack_sort(pscoll, n_ps);
    } else if (op == STACK_DESTROY) {
        for (j=0; j<n_data; j++) {
	    free_file_collection(datacoll[j]);
	}
        free(datacoll);
        datacoll = NULL;
        n_data = 0;
        n_data_popped = 0;

        for (j=0; j<n_ps; j++) {
	    free_file_collection(pscoll[j]);
	}
        free(pscoll);
        pscoll = NULL;
        n_ps = 0;
        n_ps_popped = 0;
    } 

    return ret;
}

static int push_collection (file_collection *collection)
{
    return (collection_stack(collection, STACK_PUSH) == NULL);
}

static file_collection *pop_data_collection (void)
{
    return collection_stack(NULL, STACK_POP_DATA);
}

static file_collection *pop_ps_collection (void)
{
    return collection_stack(NULL, STACK_POP_PS);
}

void destroy_file_collections (void)
{
    collection_stack(NULL, STACK_DESTROY);
}

static void reset_data_stack (void)
{
    collection_stack(NULL, STACK_RESET_DATA);
}

static void reset_ps_stack (void)
{
    collection_stack(NULL, STACK_RESET_PS);
}

static void sort_data_stack (void)
{
    collection_stack(NULL, STACK_SORT_DATA);
}

static void sort_ps_stack (void)
{
    collection_stack(NULL, STACK_SORT_PS);
}

static int get_file_collections_from_dir (const char *dname, DIR *dir)
{
    file_collection *coll;
    struct dirent *dirent;
    size_t n;
    int err = 0;

    while (!err && (dirent = readdir(dir))) { 
	if (strstr(dirent->d_name, "descriptions")) {
#if COLL_DEBUG
	    fprintf(stderr, "   %s: looking at '%s'\n", dname, dirent->d_name);
#endif
	    n = strlen(dirent->d_name);
	    if (!strcmp(dirent->d_name + n - 12, "descriptions")) {
		coll = file_collection_new(dname, dirent->d_name, &err);
		if (coll != NULL && !err) {
		    err = push_collection(coll);
		}
	    }
	}
    }

    return err;
}

static int dont_go_there (const char *s)
{
    int ret = 0;

    if (!strcmp(s, "..") || strstr(s, ".inp") || strstr(s, ".gdt") || 
	strstr(s, ".gretl") || strstr(s, ".hdr")) {
	ret = 1;
    }

    return ret;
}

static int seek_file_collections (int location)
{
    char *tmp = NULL;
    DIR *dir, *try;  
    struct dirent *dirent;
    char *subdir;
    int i = 0, err = 0;

    if (location == DATA_SEARCH) {
	tmp = g_strdup_printf("%sdata", paths.gretldir);
    } else if (location == SCRIPT_SEARCH) {
	tmp = g_strdup_printf("%sscripts", paths.gretldir);
    } else if (location == USER_SEARCH) {
	tmp = g_strdup(paths.workdir);
	trim_slash(tmp);
    } else {
	return 1;
    }

 user_search_2:

    dir = opendir(tmp);
    if (dir == NULL) {
	return 1;
    }

#if COLL_DEBUG
    fprintf(stderr, "seeking file collections in '%s'\n", tmp);
#endif

    while (!err && (dirent = readdir(dir))) {
	if (!dont_go_there(dirent->d_name)) {
	    if (strcmp(dirent->d_name, ".")) {
		subdir = full_path(tmp, dirent->d_name);
	    } else {
		subdir = tmp;
	    }
	    try = opendir(subdir);
	    if (try != NULL) {
#if COLL_DEBUG
		fprintf(stderr, " trying in subdir '%s'\n", subdir);
#endif
		err = get_file_collections_from_dir(subdir, try);
#if COLL_DEBUG
		fprintf(stderr, " result: err = %d\n", err);
#endif
		closedir(try);
	    }
	}
    }

    closedir(dir);
    g_free(tmp);

    if (location == USER_SEARCH && i++ == 0) {
	tmp = gretl_default_workdir(&paths);
	if (tmp != NULL) {
	    goto user_search_2;
	}
    }

    return err;
}

#if COLL_DEBUG
static void print_collection (const file_collection *coll)
{
    printf("path = '%s'\n", coll->path);
    printf("descfile = '%s'\n", coll->descfile);
    if (coll->title != NULL && *coll->title != '\0') {
	printf("title = '%s'\n", coll->title);
    }
}

static void print_data_collections (void)
{
    file_collection *coll;

    printf("\n*** Data collections:\n");
    while ((coll = pop_data_collection())) {
	print_collection(coll);
    }
    reset_data_stack();
}

static void print_script_collections (void)
{
    file_collection *coll;

    printf("\n*** Script collections:\n");

    while ((coll = pop_ps_collection())) {
	print_collection(coll);
    }
    reset_ps_stack();
}
#endif

static int build_file_collections (void)
{
    static int built;
    static int err;

    if (!built) {
	err = seek_file_collections(DATA_SEARCH);
	if (!err) {
	    err = seek_file_collections(SCRIPT_SEARCH);
	}
	if (!err) {
	    err = seek_file_collections(USER_SEARCH);
	}
	if (!err) {
	    sort_data_stack();
	    sort_ps_stack();
	}
	built = 1;
    }

#if COLL_DEBUG
    print_data_collections();
    print_script_collections();
#endif

    return err;
}

char *strip_extension (char *s)
{
    char *p = strrchr(s, '.');
    
    if (p != NULL && 
	(!strcmp(p, ".gdt") || !strcmp(p, ".inp") ||
	 !strcmp(p, ".bin") || !strcmp(p, ".gfn") ||
	 !strcmp(p, ".bn7"))) {
	*p = '\0';
    }

    return s;
}

static int read_file_descriptions (windata_t *win, gpointer p)
{
    FILE *fp;
    GtkListStore *store;
    GtkTreeSelection *selection;
    GtkTreeIter iter;
    char line[MAXLEN];
    char *index;
    file_collection *coll = (file_collection *) p;

    index = full_path(coll->path, coll->descfile);

    fp = gretl_fopen(index, "r");
    if (fp == NULL) {
	return 1;
    }

    store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(win->listbox)));
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);
    
    while (fgets(line, sizeof line, fp)) {
	char fname[24], descrip[80], data[64];

	if (*line == '#') continue;

	if (win->role == TEXTBOOK_DATA) {
	    if (sscanf(line, " \"%23[^\"]\",\"%79[^\"]\"", 
		       fname, descrip) == 2) {
		gtk_list_store_append(store, &iter);
		gtk_list_store_set(store, &iter, 
				   0, strip_extension(fname), 
				   1, descrip, -1);
	    }
	} else { /* script files */
	    if (sscanf(line, " \"%23[^\"]\",\"%79[^\"]\",\"%63[^\"]\"", 
		       fname, descrip, data) == 3) {
		gtk_list_store_append(store, &iter);
		gtk_list_store_set(store, &iter, 
				   0, strip_extension(fname), 
				   1, descrip, 
				   2, data, -1);
	    }
	}
    }

    fclose(fp);

    /* select the first row */
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);
    selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(win->listbox));
    gtk_tree_selection_select_iter(selection, &iter);
    
    return 0;
}

static void show_datafile_info (GtkWidget *w, gpointer data)
{
    char hdrname[MAXLEN];
    windata_t *vwin = (windata_t *) data;
    char *descrip;
    PRN *prn;
    gchar *fname;
    file_collection *coll;

    tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), vwin->active_var,
			 0, &fname);
    coll = g_object_get_data(G_OBJECT(vwin->listbox), "coll");

    build_path(hdrname, coll->path, fname, ".gdt");
    g_free(fname);

    descrip = gretl_get_gdt_description(hdrname);

    if (descrip != NULL) {
	prn = gretl_print_new_with_buffer(descrip);
	view_buffer(prn, 80, 320, _("gretl: data header"), INFO, NULL);
    } else {
	errbox(_("Failed to retrieve description of data"));
	fprintf(stderr, I_("didn't get description from %s\n"), hdrname);
    }
}

void browser_open_data (GtkWidget *w, gpointer data)
{
    windata_t *vwin = (windata_t *) data;
    gchar *datname;
    file_collection *coll;

    tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), vwin->active_var, 
			 0, &datname);
    coll = g_object_get_data(G_OBJECT(vwin->listbox), "coll");

    build_path(tryfile, coll->path, datname, ".gdt");
    g_free(datname);

    set_datapage(coll->title);

    verify_open_data(vwin, OPEN_DATA);
} 

void browser_open_ps (GtkWidget *w, gpointer data)
{
    windata_t *vwin = (windata_t *) data;
    gchar *fname;
    file_collection *coll;

    tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), vwin->active_var, 
			 0, &fname);
    coll = g_object_get_data(G_OBJECT(vwin->listbox), "coll");

    build_path(scriptfile, coll->path, fname, ".inp");
    g_free(fname);

    /* close the calling window */
    gtk_widget_destroy(GTK_WIDGET(vwin->main));

    set_scriptpage(coll->title);

    view_file(scriptfile, 0, 0, 78, 370, VIEW_SCRIPT);
} 

enum {
    VIEW_FN_PKG_INFO,
    VIEW_FN_PKG_CODE,
    LOAD_FN_PKG,
    EDIT_FN_PKG,
    DELETE_FN_PKG,
    CALL_FN_PKG
};

static int gui_load_user_functions (const char *fname)
{
    int err;

    err = load_user_function_file(fname);
    if (err) {
	gui_errmsg(err);
    } 

    return err;
}

static int gui_delete_fn_pkg (const char *fname, windata_t *vwin)
{
    char *msg = g_strdup_printf(_("Really delete %s?"), fname);
    int err;
    
    if (yes_no_dialog(_("gretl: delete"), msg, 0) != GRETL_YES) {
        return 0;
    }

    err = remove(fname);
    if (err) {
	file_write_errbox(fname);
    } else {
	GtkTreeModel *mod;
	GtkTreeIter iter;
	int i = 0;

	mod = gtk_tree_view_get_model(GTK_TREE_VIEW(vwin->listbox));
	gtk_tree_model_get_iter_first(mod, &iter);
	while (i < vwin->active_var) {
	    gtk_tree_model_iter_next(mod, &iter);
	    i++;
	}
	gtk_list_store_remove(GTK_LIST_STORE(mod), &iter);
    }

    return err;
}

windata_t *gui_show_function_info (const char *fname, int role)
{
    char *pkgname = NULL;
    windata_t *vwin = NULL;
    PRN *prn;
    int err;

    if (bufopen(&prn)) {
	return NULL;
    }

    if (role == VIEW_FUNC_INFO) {
	err = get_function_file_info(fname, prn, &pkgname);
    } else {
	err = get_function_file_code(fname, prn, &pkgname);
    }
	
    if (err) {
	gretl_print_destroy(prn);
	gui_errmsg(err);
    } else {
	gchar *title;

	title = g_strdup_printf("gretl: %s", (pkgname)? 
				pkgname : _("function code"));
	vwin = view_buffer(prn, 78, 350, title, role, NULL);
	strcpy(vwin->fname, fname);
	if (strstr(fname, "dltmp")) {
	    set_window_delete_filename(vwin);
	}
	g_free(title);
    }

    free(pkgname);

    return vwin;
}

static void browser_functions_handler (windata_t *vwin, int task)
{
    char fnfile[FILENAME_MAX];
    gchar *fname;
    gchar *dir;
    int dircol = 0;
    int err = 0;

    if (vwin->role == FUNC_FILES) {
	dircol = 4;
    } else if (vwin->role != REMOTE_FUNC_FILES) {
	dummy_call();
	return;
    }

    tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), vwin->active_var, 
			 0, &fname);

    if (dircol != 0) {
	tree_view_get_string(GTK_TREE_VIEW(vwin->listbox), vwin->active_var, 
			     dircol, &dir);
	build_path(fnfile, dir, fname, ".gfn");
	g_free(dir);
    } else {
	strcpy(fnfile, fname);
    }

    g_free(fname);

    if (task == LOAD_FN_PKG) {
	err = gui_load_user_functions(fnfile);
    } else if (task == DELETE_FN_PKG) {
	err = gui_delete_fn_pkg(fnfile, vwin);
    } else if (task == VIEW_FN_PKG_INFO) {
	gui_show_function_info(fnfile, VIEW_FUNC_INFO);
    } else if (task == VIEW_FN_PKG_CODE) {
	gui_show_function_info(fnfile, VIEW_FUNC_CODE);
    } else if (task == EDIT_FN_PKG) {
	edit_function_package(fnfile, &err);
    } else if (task == CALL_FN_PKG) {
	call_function_package(fnfile, vwin->main, &err);
    }

    if (!err && (task == EDIT_FN_PKG || task == CALL_FN_PKG)) {
	maybe_update_func_files_window(0);
    }
}

void browser_load_func (GtkWidget *w, gpointer data)
{
    windata_t *vwin = (windata_t *) data;

    browser_functions_handler(vwin, LOAD_FN_PKG);
} 

void browser_edit_func (GtkWidget *w, gpointer data)
{
    windata_t *vwin = (windata_t *) data;

    browser_functions_handler(vwin, EDIT_FN_PKG);
}

void browser_call_func (GtkWidget *w, gpointer data)
{
    windata_t *vwin = (windata_t *) data;

    browser_functions_handler(vwin, CALL_FN_PKG);
}

static void show_function_info (GtkWidget *w, gpointer data)
{
    windata_t *vwin = (windata_t *) data;

    browser_functions_handler(vwin, VIEW_FN_PKG_INFO);
}

static void show_function_code (GtkWidget *w, gpointer data)
{
    windata_t *vwin = (windata_t *) data;

    browser_functions_handler(vwin, VIEW_FN_PKG_CODE);
}

static void browser_del_func (GtkWidget *w, gpointer data)
{
    windata_t *vwin = (windata_t *) data;

    browser_functions_handler(vwin, DELETE_FN_PKG);
} 

static void set_browser_status (windata_t *vwin, int status)
{
    int i = vwin->role - TEXTBOOK_DATA;

    if (status == BROWSER_BUSY) {
	browsers[i] = vwin->main;
    } else {
	browsers[i] = NULL;
    }
} 

static void free_browser (GtkWidget *w, gpointer data)
{
    windata_t *vwin = (windata_t *) data;

    set_browser_status(vwin, BROWSER_OK);
    free_windata(NULL, data);
}

static gpointer get_browser_ptr (int role)
{
    int i = role - TEXTBOOK_DATA;

    return (gpointer) &(browsers[i]);
}

static gpointer get_browser (int role)
{
    int i = role - TEXTBOOK_DATA;

    return browsers[i];
}

static int browser_busy (guint code)
{
    int ret = 0;

    if (code >= TEXTBOOK_DATA && code <= REMOTE_DB) {
	GtkWidget *w;

	w = browsers[code - TEXTBOOK_DATA];
	if (w != NULL) {
	    gdk_window_raise(w->window);
	    ret = 1;
	}
    }

    return ret;
}

windata_t *get_local_viewer (int remote_role)
{
    windata_t *vwin = NULL;
    GtkWidget *w = NULL;

    if (remote_role == REMOTE_DB) {
	w = get_browser(NATIVE_DB);
    } else if (remote_role == REMOTE_FUNC_FILES) {
	w = get_browser(FUNC_FILES);
    }

    if (w != NULL) {
	vwin = g_object_get_data(G_OBJECT(w), "vwin");
    }

    return vwin;
}

static void new_package_callback (GtkWidget *w , gpointer p)
{
    data_save_selection_wrapper(SAVE_FUNCTIONS, NULL);
}

static void close_files_viewer (GtkWidget *w, windata_t *vwin)
{
    gtk_widget_destroy(vwin->main);
} 

static void build_datafiles_popup (windata_t *vwin)
{
    if (vwin->popup == NULL) {
	vwin->popup = gtk_menu_new();

	add_popup_item(_("Info"), vwin->popup, 
		       G_CALLBACK(show_datafile_info), 
		       vwin);
	add_popup_item(_("Open"), vwin->popup, 
		       G_CALLBACK(browser_open_data), 
		       vwin);
    }
}   

static void build_funcfiles_popup (windata_t *vwin)
{
    if (vwin->popup == NULL) {
	vwin->popup = gtk_menu_new();

	if (vwin->role == FUNC_FILES) {
	    add_popup_item(_("Edit"), vwin->popup, 
			   G_CALLBACK(browser_edit_func), 
			   vwin);
	    add_popup_item(_("Info"), vwin->popup, 
			   G_CALLBACK(show_function_info), 
			   vwin);
	    add_popup_item(_("View code"), vwin->popup, 
			   G_CALLBACK(show_function_code), 
			   vwin);
	    add_popup_item(_("Execute"), vwin->popup, 
			   G_CALLBACK(browser_call_func), 
			   vwin);
	    add_popup_item(_("Delete"), vwin->popup, 
			   G_CALLBACK(browser_del_func), 
			   vwin);
	    add_popup_item(_("New"), vwin->popup, 
			   G_CALLBACK(new_package_callback), 
			   vwin);
	} else {
	    add_popup_item(_("Info"), vwin->popup, 
			   G_CALLBACK(file_info_from_server), 
			   vwin);
	    add_popup_item(_("Install"), vwin->popup, 
			   G_CALLBACK(install_file_from_server), 
			   vwin);
	    add_popup_item(_("Execute"), vwin->popup, 
			   G_CALLBACK(browser_call_func), 
			   vwin);
	}
    }
}

static void build_db_popup (windata_t *vwin)
{
    if (vwin->popup == NULL) {
	vwin->popup = gtk_menu_new();

	if (vwin->role == NATIVE_DB) {
	    add_popup_item(_("List series"), vwin->popup, 
			   G_CALLBACK(open_db_index), 
			   vwin);
	    add_popup_item(_("Find..."), vwin->popup, 
			   G_CALLBACK(listbox_find), 
			   vwin);
	} else {
	    add_popup_item(_("List series"), vwin->popup, 
			   G_CALLBACK(open_remote_db_index), 
			   vwin);
	    add_popup_item(_("Install"), vwin->popup, 
			   G_CALLBACK(install_file_from_server), 
			   vwin);
	    add_popup_item(_("Find..."), vwin->popup, 
			   G_CALLBACK(listbox_find), 
			   vwin);
	}
    }
}

static void show_server_dbs (GtkWidget *w, gpointer p)
{
    display_files(REMOTE_DB, p);
}

static void show_local_dbs (GtkWidget *w, gpointer p)
{
    display_files(NATIVE_DB, p);
}

static void show_server_funcs (GtkWidget *w, gpointer p)
{
    display_files(REMOTE_FUNC_FILES, p);
}

static void show_local_funcs (GtkWidget *w, gpointer p)
{
    display_files(FUNC_FILES, p);
}

enum {
    BTN_EDIT = 1,
    BTN_INFO,
    BTN_CODE,
    BTN_INDX,
    BTN_INST,
    BTN_EXEC,
    BTN_DEL,
    BTN_WWW,
    BTN_HOME,
    BTN_NEW,
    BTN_FIND,
    BTN_OPEN,
    BTN_CLOSE
};

static GretlToolItem files_items[] = {
    { N_("Open"),           GTK_STOCK_OK,         NULL,                          BTN_OPEN },
    { N_("Edit"),           GTK_STOCK_EDIT,       G_CALLBACK(browser_edit_func), BTN_EDIT },
    { N_("Info"),           GTK_STOCK_INFO,       NULL,                          BTN_INFO },
    { N_("View code"),      GTK_STOCK_PROPERTIES, G_CALLBACK(show_function_code), BTN_CODE },
    { N_("List series"),    GTK_STOCK_INDEX,      NULL,                          BTN_INDX },
    { N_("Install"),        GTK_STOCK_SAVE,       G_CALLBACK(install_file_from_server), BTN_INST },
    { N_("Execute"),        GTK_STOCK_EXECUTE,    G_CALLBACK(browser_call_func), BTN_EXEC },
    { N_("Delete"),         GTK_STOCK_DELETE,     G_CALLBACK(browser_del_func),  BTN_DEL },
    { N_("Look on server"), GTK_STOCK_NETWORK,    NULL,                          BTN_WWW },
    { N_("Local machine"),  GTK_STOCK_HOME,       NULL,                          BTN_HOME },
    { N_("New"),            GTK_STOCK_NEW,        G_CALLBACK(new_package_callback), BTN_NEW },
    { N_("Find..."),        GTK_STOCK_FIND,       G_CALLBACK(listbox_find),       BTN_FIND },
    { N_("Close"),          GTK_STOCK_CLOSE,      G_CALLBACK(close_files_viewer), BTN_CLOSE }
};

static int n_files_items = G_N_ELEMENTS(files_items);

#define common_item(f) (f == BTN_FIND || f == BTN_CLOSE)

#define local_funcs_item(f) (f == BTN_EDIT || f == BTN_NEW || \
			     f == BTN_DEL || f == BTN_CODE)

static int files_item_get_callback (GretlToolItem *item, int role)
{
    if (common_item(item->flag)) {
	return 1;
    } else if (local_funcs_item(item->flag)) {
	return (role == FUNC_FILES);
    } else if (item->flag == BTN_INST) {
	return (role == REMOTE_DB || role == REMOTE_FUNC_FILES);
    } else if (item->flag == BTN_EXEC) {
	return (role == FUNC_FILES || role == REMOTE_FUNC_FILES);
    }

    item->func = NULL;

    if (item->flag == BTN_OPEN) {
	/* open: only data files and scripts */
	if (role == TEXTBOOK_DATA) {
	    item->func = G_CALLBACK(browser_open_data);
	} else if (role == PS_FILES) {
	    item->func = G_CALLBACK(browser_open_ps);
	} 
    } else if (item->flag == BTN_INFO) {
	if (role == TEXTBOOK_DATA) {
	    item->func = G_CALLBACK(show_datafile_info);
	} else if (role == FUNC_FILES) {
	    item->func = G_CALLBACK(show_function_info);
	} else if (role == REMOTE_FUNC_FILES) {
	    item->func = G_CALLBACK(file_info_from_server);
	} 
    } else if (item->flag == BTN_INDX) {
	/* index: databases only */
	if (role == NATIVE_DB) {
	    item->func = G_CALLBACK(open_db_index);
	} else if (role == REMOTE_DB) {
	    item->func = G_CALLBACK(open_remote_db_index);
	} 
    } else if (item->flag == BTN_WWW) {
	if (role == FUNC_FILES) {
	    item->func = G_CALLBACK(show_server_funcs);
	} else if (role == NATIVE_DB) {
	    item->func = G_CALLBACK(show_server_dbs);
	} 
    } else if (item->flag == BTN_HOME) {
	/* home: show only for on-server items */
	if (role == REMOTE_FUNC_FILES) {
	    item->func = G_CALLBACK(show_local_funcs);	
	} else if (role == REMOTE_DB) {
	    item->func = G_CALLBACK(show_local_dbs);
	} 
    } 

    return (item->func != NULL);
}

static void make_filesbar (windata_t *vwin)
{
    GtkWidget *hbox;
    GretlToolItem *item;
    int i;

    hbox = gtk_hbox_new(FALSE, 0);
    gtk_box_pack_start(GTK_BOX(vwin->vbox), hbox, FALSE, FALSE, 0);

    vwin->mbar = gretl_toolbar_new();

    for (i=0; i<n_files_items; i++) {
	item = &files_items[i];
	if (files_item_get_callback(item, vwin->role)) {
	    gretl_toolbar_insert(vwin->mbar, item, item->func, vwin, -1);
	}
    }

    gtk_box_pack_start(GTK_BOX(hbox), vwin->mbar, FALSE, FALSE, 0);
    gtk_widget_show_all(hbox);
}

static gchar *files_title (int code)
{
    static char hname[48];
    gchar *ret = NULL;

    if (*hname == '\0') {
#if GTK_MAJOR_VERSION == 2 && GTK_MINOR_VERSION < 8
	/* no g_get_host_name() */
	strcpy(hname, _("local machine"));
#else
	const gchar *s = g_get_host_name();

	if (s != NULL && strlen(s) < 48) {
	    strcpy(hname, s);
	} else {
	    strcpy(hname, _("local machine"));
	}
#endif
    }

    if (code == NATIVE_DB) {
	ret = g_strdup_printf(_("gretl: databases on %s"), hname);
    } else {
	ret = g_strdup_printf(_("gretl: function packages on %s"), hname);
    }

    return ret;
}

static void  
db_window_handle_drag  (GtkWidget *widget,
			GdkDragContext *context,
			gint x,
			gint y,
			GtkSelectionData *data,
			guint info,
			guint time,
			gpointer p)
{
    /* handle drag of pointer from remote database window */
    if (info == GRETL_REMOTE_DB_PTR && data != NULL && 
	data->type == GDK_SELECTION_TYPE_INTEGER) {
	install_file_from_server(NULL, *(void **) data->data);
    }
}

static void  
pkg_window_handle_drag  (GtkWidget *widget,
			 GdkDragContext *context,
			 gint x,
			 gint y,
			 GtkSelectionData *data,
			 guint info,
			 guint time,
			 gpointer p)
{
    /* handle drag of pointer from remote function package window */
    if (info == GRETL_REMOTE_FNPKG_PTR && data != NULL && 
	data->type == GDK_SELECTION_TYPE_INTEGER) {
	install_file_from_server(NULL, *(void **) data->data);
    }
}

static void set_up_viewer_drag_target (windata_t *vwin)
{
    GCallback callback;
    int i;

    if (vwin->role == NATIVE_DB) {
	i = GRETL_REMOTE_DB_PTR;
	callback = G_CALLBACK(db_window_handle_drag);
    } else if (vwin->role == FUNC_FILES) {
	i = GRETL_REMOTE_FNPKG_PTR;
	callback = G_CALLBACK(pkg_window_handle_drag);
    } else {
	return;
    }

    gtk_drag_dest_set(vwin->listbox,
		      GTK_DEST_DEFAULT_ALL,
		      &gretl_drag_targets[i], 1,
		      GDK_ACTION_COPY);

    g_signal_connect(G_OBJECT(vwin->listbox), "drag_data_received",
		     callback, NULL);
}

void display_files (int code, gpointer p)
{
    GtkWidget *filebox;
    windata_t *vwin;
    struct fpkg_response fresp;
    gchar *title = NULL;
    int err = 0;

    if (browser_busy(code)) {
	/* an appropriate window is already open */
	return;
    }

    if (code == FUNC_FILES || code == NATIVE_DB) {
	title = files_title(code);
    } else if (code == PS_FILES) {
	title = g_strdup(_("gretl: practice files"));
    } else if (code == TEXTBOOK_DATA) {
	title = g_strdup(_("gretl: data files"));
    } else if (code == REMOTE_DB) {
	title = g_strdup(_("gretl: databases on server"));
    } else if (code == REMOTE_FUNC_FILES) {
	title = g_strdup(_("gretl: function packages on server"));
    }

    vwin = gretl_browser_new(code, title, 0);
    g_free(title);

    fpkg_response_init(&fresp, p);
    g_signal_connect(G_OBJECT(vwin->main), "destroy",
		     G_CALLBACK(free_browser),
		     vwin);

    set_browser_status(vwin, BROWSER_BUSY);

    if (code == REMOTE_DB) {
	gtk_window_set_default_size(GTK_WINDOW(vwin->main), 640, 480);
    }

    /* set up grids */
    vwin->vbox = gtk_vbox_new(FALSE, 1);
    gtk_box_set_spacing(GTK_BOX(vwin->vbox), 4);
    gtk_container_set_border_width(GTK_CONTAINER(vwin->vbox), 4);
    gtk_container_add(GTK_CONTAINER(vwin->main), vwin->vbox);

    make_filesbar(vwin);

    if (code == TEXTBOOK_DATA || code == PS_FILES) {
	/* we'll need more than one tab */
	filebox = files_notebook(vwin, code);
    } else {
	filebox = files_window(vwin);
    }

    gtk_box_pack_start(GTK_BOX(vwin->vbox), filebox, TRUE, TRUE, 0);

    g_object_set_data(G_OBJECT(vwin->main), "vwin", vwin);

    if (code == TEXTBOOK_DATA) { 
	file_collection *coll;

	build_datafiles_popup(vwin);
	while ((coll = pop_data_collection())) {
	    g_signal_connect(G_OBJECT(coll->page), "button_press_event",
			     G_CALLBACK(popup_menu_handler), 
			     vwin->popup);
	}
	reset_data_stack();
    } else if (code == FUNC_FILES || code == REMOTE_FUNC_FILES) {
	build_funcfiles_popup(vwin);
	g_signal_connect(G_OBJECT(vwin->listbox), "button_press_event",
			 G_CALLBACK(popup_menu_handler), 
			 vwin->popup);
    } else if (code == NATIVE_DB || code == REMOTE_DB) {
	build_db_popup(vwin);
	g_signal_connect(G_OBJECT(vwin->listbox), "button_press_event",
			 G_CALLBACK(popup_menu_handler), 
			 vwin->popup);
    }	

    if (REMOTE_ACTION(code)) {
	GtkWidget *hbox;

	hbox = gtk_hbox_new(FALSE, 0);
	gtk_box_pack_start(GTK_BOX(vwin->vbox), hbox, FALSE, FALSE, 0);
	vwin->status = gtk_label_new(_("Network status: OK"));
	gtk_label_set_justify(GTK_LABEL(vwin->status), GTK_JUSTIFY_LEFT);
	gtk_box_pack_start(GTK_BOX(hbox), vwin->status, FALSE, FALSE, 0);
    } 

    /* put stuff into list box(es) */
    if (code == TEXTBOOK_DATA || code == PS_FILES) {
	err = populate_notebook_filelists(vwin, filebox, code);
    } else if (code == FUNC_FILES) {
	err = populate_filelist(vwin, &fresp);
	if (fresp.col1_width > 15) {
	    /* widen the file box if need be */
	    gint w, h;

	    gtk_widget_get_size_request(filebox, &w, &h);
	    w += 100;
	    gtk_widget_set_size_request(filebox, w, h);
	}
    } else {
	err = populate_filelist(vwin, NULL);
    }

    if (err) {
	gtk_widget_destroy(vwin->main);
    } else {
	gtk_widget_show_all(vwin->main); 
	gtk_widget_grab_focus(vwin->listbox);
	if (code == NATIVE_DB || code == FUNC_FILES) {
	    set_up_viewer_drag_target(vwin);
	}
    }

    if (err && code == FUNC_FILES && fresp.try_server == 1) {
	/* no function packages on local machine */
	display_files(REMOTE_FUNC_FILES, p);
    }
}

static int display_files_code (const gchar *s)
{
    if (!strcmp(s, "DisplayDataFiles"))
	return TEXTBOOK_DATA;
    if (!strcmp(s, "DisplayScripts"))
	return PS_FILES;
    if (!strcmp(s, "NativeDB"))
	return NATIVE_DB;
    if (!strcmp(s, "RemoteDB"))
	return REMOTE_DB;
    if (!strcmp(s, "LocalGfn"))
	return FUNC_FILES;
    if (!strcmp(s, "RemoteGfn"))
	return REMOTE_FUNC_FILES;
    return 0;
}

/* make a browser window to display a set of files: textbook
   data files, practice scripts, databases...  
*/

void show_files (GtkAction *action, gpointer p)
{
    int code = display_files_code(gtk_action_get_name(action));

    display_files(code, p);
}

static int get_func_info (const char *fname, const char *fndir,
			  char **pdesc, char **pver)
{
    char fullname[FILENAME_MAX];
    int err = 0;

    build_path(fullname, fndir, fname, NULL);
    *pdesc = get_function_file_header(fullname, pver, &err);
    if (err) {
	gui_errmsg(err);
    }

    return err;
}

static int
read_fn_files_in_dir (int role, DIR *dir, const char *fndir, 
		      GtkListStore *store, GtkTreeIter *iter,
		      int *maxlen)
{
    struct dirent *dirent;
    char fullname[MAXLEN];
    gboolean loaded;
    char *fname;
    char *descrip;
    char *version;
    int n, nfn = 0;
    int err;

    while ((dirent = readdir(dir)) != NULL) {
	if (!strcmp(dirent->d_name, ".") ||
	    !strcmp(dirent->d_name, "..")) {
	    continue;
	}
	fname = g_strdup(dirent->d_name);
	if (fname == NULL) {
	    break;
	}
	n = strlen(fname);
	if (!g_ascii_strcasecmp(fname + n - 4, ".gfn")) {
	    err = get_func_info(fname, fndir, &descrip, &version);
	    if (!err) {
		gtk_list_store_append(store, iter);
		build_path(fullname, fndir, fname, NULL);
		n -= 4;
		fname[n] = '\0';
		if (n > *maxlen) {
		    *maxlen = n;
		}
		loaded = function_package_is_loaded(fullname);
		gtk_list_store_set(store, iter, 
				   0, fname, 
				   1, version,
				   2, descrip, 
				   3, loaded, 
				   4, fndir, -1);
		free(descrip);
		free(version);
		nfn++;
	    }
	}
	g_free(fname);
    }

    return nfn;
}

static gint 
compare_pkgnames (GtkTreeModel *model, GtkTreeIter *a, GtkTreeIter *b,
		  gpointer p)
{
    gchar *t1, *t2;
    gint ret;

    gtk_tree_model_get(model, a, 0, &t1, -1);
    gtk_tree_model_get(model, b, 0, &t2, -1);

    lower(t1);
    lower(t2);

    ret = strcmp(t1, t2);

    g_free(t1);
    g_free(t2);

    return ret;    
}

static void sort_pkglist (windata_t *vwin)
{
    GtkTreeModel *model;

    model = gtk_tree_view_get_model(GTK_TREE_VIEW(vwin->listbox));

    gtk_tree_sortable_set_sort_column_id(GTK_TREE_SORTABLE(model), 
					 0, GTK_SORT_ASCENDING);
    gtk_tree_sortable_set_sort_func(GTK_TREE_SORTABLE(model), 0,
				    compare_pkgnames, NULL, NULL);
}

gint populate_func_list (windata_t *vwin, struct fpkg_response *fresp)
{
    GtkListStore *store;
    GtkTreeIter iter;
    char fndir[FILENAME_MAX];
    char *tmp;
    DIR *dir;
    int maxlen = 0;
    int nfn = 0;

    store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(vwin->listbox)));
    gtk_list_store_clear(store);
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    /* pick up any function files in system dir */
    build_path(fndir, paths.gretldir, "functions", NULL);
    dir = opendir(fndir);
    if (dir != NULL) {
	nfn += read_fn_files_in_dir(vwin->role, dir, fndir, store, &iter,
				    &maxlen);
	closedir(dir);
    }

    /* plus any function files in the user's working dir */
    build_path(fndir, paths.workdir, "functions", NULL);
    dir = opendir(fndir);
    if (dir != NULL) {
	nfn += read_fn_files_in_dir(vwin->role, dir, fndir, store, &iter,
				    &maxlen);
	closedir(dir);
    }

    /* plus any function files in the user's personal data dir */
    build_path(fndir, paths.dotdir, "functions", NULL);
    dir = opendir(fndir);
    if (dir != NULL) {
	nfn += read_fn_files_in_dir(vwin->role, dir, fndir, store, &iter,
				    &maxlen);
	closedir(dir);
    } 

    /* plus any in the default working dir, if not already searched */
    tmp = gretl_default_workdir(&paths);
    if (tmp != NULL) {
	build_path(fndir, tmp, "functions", NULL);
	dir = opendir(fndir);
	if (dir != NULL) {
	    nfn += read_fn_files_in_dir(vwin->role, dir, fndir, store, &iter,
				    &maxlen);
	    closedir(dir);
	}
	free(tmp);
    }

    if (nfn == 0) {
	if (fresp->try_server == 0) {
	    int resp;

	    resp = yes_no_dialog(_("gretl: function packages"),
				 _("No gretl function packages were found on this computer.\n"
				   "Do you want to take a look on the gretl server?"),
				 0);
	    if (resp == GRETL_YES) {
		fresp->try_server = 1;
	    }
	} else {
	    warnbox(_("No gretl function packages were found on this computer."));
	}
	return 1;
    } else {
	sort_pkglist(vwin);
    }

    return 0;
}

/* update loaded status after run, edit calls */

void maybe_update_func_files_window (int editing)
{
    GtkWidget *w = get_browser(FUNC_FILES);
    windata_t *vwin = NULL;
    GtkWidget *lbox = NULL;

    if (w != NULL) {
	vwin = g_object_get_data(G_OBJECT(w), "vwin");
    }

    if (vwin != NULL) {
	lbox = vwin->listbox;
    }

    if (editing && lbox != NULL) {
	GtkListStore *store;

	store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(lbox)));	
	gtk_list_store_clear(store);
	populate_func_list(vwin, NULL);
	return;
    }

    if (lbox != NULL) {
	GtkListStore *store = 
	    GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(lbox)));
	GtkTreeIter iter;
	char fullname[MAXLEN];
	gchar *fname;
	gchar *fdir;
	gboolean loaded;

	if (!gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter)) {
	    return;
	}

	while (1) {
	    gtk_tree_model_get(GTK_TREE_MODEL(store), &iter, 0, &fname, 
			       4, &fdir, -1);
	    sprintf(fullname, "%s%c%s.gfn", fdir, SLASH, fname);
	    g_free(fname);
	    g_free(fdir);
	    loaded = function_package_is_loaded(fullname);
	    gtk_list_store_set(store, &iter, 3, loaded, -1);
	    if (!gtk_tree_model_iter_next(GTK_TREE_MODEL(store), &iter)) {
		break;
	    }
	}
    }
}

gint populate_filelist (windata_t *vwin, gpointer p)
{
    if (vwin->role == NATIVE_DB) {
	return populate_dbfilelist(vwin);
    } else if (vwin->role == REMOTE_DB) {
	return populate_remote_db_list(vwin);
    } else if (vwin->role == REMOTE_FUNC_FILES) {
	return populate_remote_func_list(vwin);
    } else if (vwin->role == FUNC_FILES) {
	return populate_func_list(vwin, p);
    } else {
	return read_file_descriptions(vwin, p);
    }
}

static GtkWidget *files_window (windata_t *vwin) 
{
    const char *data_titles[] = {
	N_("File"), 
	N_("Summary")
    };
    const char *ps_titles[] = {
	N_("Script"), 
	N_("Topic"), 
	N_("Data")
    };
    const char *db_titles[] = {
	N_("Database"), 
	N_("Source")
    };
    const char *remote_db_titles[] = {
	N_("Database"), 
	N_("Source"), 
	N_("Local status")
    };
    const char *func_titles[] = {
	N_("Package"), 
	N_("Version"),
	N_("Summary"), 
	N_("Loaded?")
    };
    const char *remote_func_titles[] = {
	N_("Package"), 
	N_("Version"),
	N_("Summary"), 
	N_("Local status")
    };

    GType types_2[] = {
	G_TYPE_STRING,
	G_TYPE_STRING
    };
    GType types_3[] = {
	G_TYPE_STRING,
	G_TYPE_STRING,
	G_TYPE_STRING
    };
    GType types_4[] = {
	G_TYPE_STRING,
	G_TYPE_STRING,
	G_TYPE_STRING,
	G_TYPE_STRING
    };
    GType types_5[] = {
	G_TYPE_STRING,
	G_TYPE_STRING,
	G_TYPE_STRING,
	G_TYPE_BOOLEAN,
	G_TYPE_STRING
    };

    const char **titles = data_titles;
    GType *types = types_2;

    int full_width = 500, file_height = 260;
    int hidden_col = 0;
    int use_tree = 0;

    GtkWidget *box;
    int cols = 2;

    switch (vwin->role) {
    case NATIVE_DB:
	titles = db_titles;
	cols = 3;
	hidden_col = TRUE;
	break;
    case REMOTE_DB:
	titles = remote_db_titles;
	cols = 3;
	full_width = 580;
	use_tree = 1;
	break;
    case PS_FILES:
	titles = ps_titles;
	cols = 3;
	full_width = 560;
	file_height = 300;
	break;
    case FUNC_FILES:
	titles = func_titles;
	cols = 5;
	types = types_5;
	hidden_col = TRUE;
	full_width = 560;
	file_height = 320;
	break;
    case REMOTE_FUNC_FILES:
	titles = remote_func_titles;
	cols = 4;
	types = types_4;
	full_width = 580;
	break;
    default:
	break;
    }

    if (cols == 3) {
	types = types_3;
    }

    full_width *= gui_scale;
    file_height *= gui_scale;

    box = gtk_vbox_new(FALSE, 0);
    gtk_widget_set_size_request(box, full_width, file_height);
    vwin_add_list_box(vwin, GTK_BOX(box), cols, hidden_col, 
		      types, titles, use_tree);
    gtk_widget_show(box);

    return box;
}

static void 
switch_file_page_callback (GtkNotebook *notebook, GtkNotebookPage *page,
			   guint page_num, windata_t *vwin)
{
    gpointer p = g_object_get_data(G_OBJECT(notebook), "browse_ptr");

    if (p != NULL) {
	GtkWidget *w = *(GtkWidget **) p;

	if (w != NULL) {
	    char wnum[4];

	    sprintf(wnum, "%d", (int) page_num);
	    vwin->listbox = g_object_get_data(G_OBJECT(notebook), wnum);
	}
    }
}

/* below: Construct a set of notebook pages for either data files
   (Ramanathan, Wooldridge, etc.) or practice scripts.  Creates the
   pages but does not yet fill them out.
*/

static GtkWidget *files_notebook (windata_t *vwin, int code)
{
    GtkWidget *notebook;
    GtkWidget *listpage;
    GtkWidget *label;
    char wnum[4];
    file_collection *coll;
    int j;

    if (code != TEXTBOOK_DATA && code != PS_FILES) {
	return NULL;
    }

    build_file_collections(); /* FIXME check for errors */

    notebook = gtk_notebook_new();
    gtk_notebook_set_scrollable(GTK_NOTEBOOK(notebook), TRUE);

    j = 0;

    while (1) {
	if (code == TEXTBOOK_DATA) {
	    coll = pop_data_collection();
	} else {
	    coll = pop_ps_collection();
	}

	if (coll == NULL) break;

	listpage = files_window(vwin);
	label = gtk_label_new(coll->title);
	gtk_widget_show(label);
	gtk_notebook_append_page(GTK_NOTEBOOK(notebook), listpage, label);
	coll->page = vwin->listbox;
	sprintf(wnum, "%d", j);
	g_object_set_data(G_OBJECT(notebook), wnum, coll->page);
	g_object_set_data(G_OBJECT(coll->page), "coll", coll);
	j++;
    }

    if (code == TEXTBOOK_DATA) {
	reset_data_stack();
    } else {
	reset_ps_stack();
    }

    g_object_set_data(G_OBJECT(GTK_NOTEBOOK(notebook)), "browse_ptr",
		      get_browser_ptr(vwin->role));
    g_signal_connect(G_OBJECT(GTK_NOTEBOOK(notebook)), "switch-page",
		     G_CALLBACK(switch_file_page_callback),
		     vwin);

    if (gtk_notebook_get_n_pages(GTK_NOTEBOOK(notebook)) > 5) {
	gtk_notebook_popup_enable(GTK_NOTEBOOK(notebook));
    }

    gtk_widget_show(notebook);

    return notebook;
}

/* below: fill out a set of notebook pages (for data files
   or scripts files), entering the details into the page.
*/

static int populate_notebook_filelists (windata_t *vwin, 
					GtkWidget *notebook,
					int code)
{
    file_collection *coll;
    const char *title;
    int gotpref = 0;
    int gotcol = 0;
    int j;

    if (vwin == NULL) {
	return 1;
    }

    while (1) {
	if (code == TEXTBOOK_DATA) {
	    coll = pop_data_collection();
	} else {
	    coll = pop_ps_collection();
	}

	if (coll != NULL) {
	    gotcol = 1;
	} else {
	    break;
	}

	vwin->listbox = coll->page;
	populate_filelist(vwin, coll);
    }

    if (!gotcol) {
	return 1;
    }

    j = 0;

    if (code == TEXTBOOK_DATA) {
	reset_data_stack();
	title = get_datapage();
	if (*title != '\0') {
	    while ((coll = pop_data_collection())) {
		if (!strcmp(coll->title, title)) {
		    gotpref = 1;
		    break;
		}
		j++;
	    }
	    if (!gotpref) {
		reset_data_stack();
		coll = pop_data_collection();
		j = 0;
	    }
	} else {
	    coll = pop_data_collection();
	}
    } else {
	reset_ps_stack();
	title = get_scriptpage();
	if (*title != '\0') {
	    while ((coll = pop_ps_collection())) {
		if (!strcmp(coll->title, title)) {
		    gotpref = 1;
		    break;
		}
		j++;
	    }
	    if (!gotpref) {
		reset_ps_stack();
		coll = pop_ps_collection();
		j = 0;
	    }	    
	} else {
	    coll = pop_ps_collection();
	}
    }

    vwin->listbox = coll->page;

    if (code == TEXTBOOK_DATA) {
	reset_data_stack();
    } else {
	reset_ps_stack();
    }

    gtk_notebook_set_current_page(GTK_NOTEBOOK(notebook), j);
    gtk_widget_grab_focus(vwin->listbox);

    return 0;
} 
