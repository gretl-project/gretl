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

#define FULL_XML_HEADERS

#include "libgretl.h"
#include "version.h"
#include "build.h"
#include "gretl_xml.h"

/* Note: here's the canonical listing of gretl addons */

static const char *addon_names[] = {
    "SVAR", "gig", "HIP", "ivpanel",
    "dbnomics", "extra", "geoplot",
    "regls", "logging", "KFgui", NULL
};

typedef struct addon_basics_ {
    char path[MAXLEN];
    char *version;
    char *date;
    char *descrip;
    int ok;
} addon_basics;

static void clear_addon_basics (addon_basics *ab)
{
    ab->path[0] = '\0';
    free(ab->version);
    free(ab->date);
    free(ab->descrip);
    ab->version = ab->date = ab->descrip = NULL;
    ab->ok = 0;
}

static int get_n_addons (void)
{
    static int n;

    if (n == 0) {
	int i;

	for (i=0; addon_names[i] != NULL; i++) {
	    n++;
	}
    }

    return n;
}

/* Determine if @pkgname is the name of an addon: @pkgname may be
   given with or without the ".gfn" suffix.
*/

int is_gretl_addon (const char *pkgname)
{
    int i, n = get_n_addons();

    if (has_suffix(pkgname, ".gfn")) {
	int n = strlen(pkgname) - 4;

	for (i=0; i<n; i++) {
	    if (!strncmp(pkgname, addon_names[i], n)) {
		return 1;
	    }
	}
    } else {
	for (i=0; i<n; i++) {
	    if (!strcmp(pkgname, addon_names[i])) {
		return 1;
	    }
	}
    }

    return 0;
}

/* Return a NULL-terminated array of addons names; optionally, if @n
   is non-NULL, supply the number of addons.
*/

const char **get_addon_names (int *n)
{
    if (n != NULL) {
	*n = get_n_addons();
    }
    return addon_names;
}

/* Given the full path to an addon gfn file, supply its version, date
   and description. The strings are newly allocated and belong to the
   caller. Return non-zero on error.
*/

static int addon_basics_from_gfn (addon_basics *ab)
{
    xmlDocPtr doc = NULL;
    xmlNodePtr node, n1, n2;
    int targ = 3;
    int got = 0;

    if (gretl_stat(ab->path, NULL) != 0) {
	return E_FOPEN;
    }

    gretl_xml_open_doc_root(ab->path, "gretl-functions", &doc, &node);
    if (doc == NULL || node == NULL) {
	return E_FOPEN;
    }

    n1 = node->xmlChildrenNode;

    while (n1 != NULL && got < targ) {
	if (!xmlStrcmp(n1->name, (XUC) "gretl-function-package")) {
	    n2 = n1->xmlChildrenNode;
	    while (n2 != NULL && got < targ) {
		if (!xmlStrcmp(n2->name, (XUC) "version")) {
		    gretl_xml_node_get_trimmed_string(n2, doc, &ab->version);
		    got++;
		} else if (!xmlStrcmp(n2->name, (XUC) "date")) {
		    gretl_xml_node_get_trimmed_string(n2, doc, &ab->date);
		    got++;
		} else if (!xmlStrcmp(n2->name, (XUC) "description")) {
		    gretl_xml_node_get_trimmed_string(n2, doc, &ab->descrip);
		    got++;
		}
		n2 = n2->next;
	    }
	}
	n1 = n1->next;
    }

    if (doc != NULL) {
	xmlFreeDoc(doc);
    }

    return (got < targ)? E_DATA: 0;
}

static char *get_user_path (char *targ, const char *pgkname,
			    const char *gfnname)
{
#ifdef OS_OSX
    return gretl_build_path(targ, gretl_app_support_dir(),
			    "functions", pgkname, gfnname, NULL);
#else
    return gretl_build_path(targ, gretl_dotdir(), "functions",
			    pgkname, gfnname, NULL);
#endif
}

static int select_addon_file (addon_basics *sys_ab,
			      addon_basics *usr_ab)
{
    int err = 0;

    if (sys_ab->ok) {
	sys_ab->ok = !strcmp(sys_ab->version, GRETL_VERSION);
    }
    if (usr_ab->ok) {
	usr_ab->ok = !strcmp(usr_ab->version, GRETL_VERSION);
    }

    if (sys_ab->ok && usr_ab->ok) {
	/* both match gretl version, compare dates */
	guint32 sed = get_epoch_day(sys_ab->date);
	guint32 ued = get_epoch_day(usr_ab->date);

	if (sed >= ued) {
	    usr_ab->ok = 0;
	} else {
	    sys_ab->ok = 0;
	}
    } else if (sys_ab->ok == 0 && usr_ab->ok == 0) {
	err = E_DATA;
    }

    return err;
}

static void report_addon_result (addon_basics *ab,
				 int err, PRN *prn)
{
    pprintf(prn, " %s '%s'\n", _("try"), ab->path);
    if (!err) {
	pprintf(prn, "  %s %s (%s)\n", _("found version"),
		ab->version, ab->date);
    } else {
	pprintf(prn, "  %s\n", _("not found"));
    }
}

/* Build a plain text index of the installed addons, holding name,
   version, date and full path, one addon per line. We allow for the
   fact that a given addon might exist both in the "system" location
   and in the user's personal filespace.  (This may happen if a user
   lacking write-permission for the system location installs or
   updates an addon.) In the case of such duplicates we determine
   which version is newer and enter its details in the index file.

   This update routine is called automatically when (a) a user
   installs an addon (FIXME: is this true for all ways an addon can be
   installed?) or (b) gretl figures out that it has been updated (new
   release or snapshot). It can also be called explicitly by the user,
   via the command "pkg index addons".

   if @prn is non-NULL verbose output will be printed, otherwise the
   function operates silently.
*/

int update_addons_index (PRN *prn)
{
    gchar *idxname = gretl_make_dotpath("addons.idx");
    addon_basics sys_ab = {0};
    addon_basics usr_ab = {0};
    char gfnname[64];
    int verbose = (prn != NULL);
    int n_addons = get_n_addons();
    FILE *fp;
    int i;

    fp = gretl_fopen(idxname, "wb");
    if (fp == NULL) {
	g_free(idxname);
	return E_FOPEN;
    }

    for (i=0; i<n_addons; i++) {
	int err = 0;

	/* construct the gfn name */
	sprintf(gfnname, "%s.gfn", addon_names[i]);
	if (verbose) {
	    pprintf(prn, _("check for %s\n"), addon_names[i]);
	}

	/* (1) build and check the system path */
	gretl_build_path(sys_ab.path, gretl_home(), "functions",
			 addon_names[i], gfnname, NULL);
	err = addon_basics_from_gfn(&sys_ab);
	if (!err) {
	    sys_ab.ok = 1;
	}
	if (verbose) {
	    report_addon_result(&sys_ab, err, prn);
	}

	/* (2) build and check the userspace path */
	get_user_path(usr_ab.path, addon_names[i], gfnname);
	err = addon_basics_from_gfn(&usr_ab);
	if (!err) {
	    usr_ab.ok = 1;
	}
	if (verbose) {
	    report_addon_result(&usr_ab, err, prn);
	}

	if (sys_ab.ok || usr_ab.ok) {
	    /* carry out further checks */
	    err = select_addon_file(&sys_ab, &usr_ab);
	} else {
	    err = E_DATA;
	}

	if (!err) {
	    /* write line to addons.idx */
	    addon_basics *ab = sys_ab.ok ? &sys_ab : &usr_ab;

	    fprintf(fp, "%s %s %s \"%s\" \"%s\"\n", addon_names[i],
		    ab->version, ab->date, ab->descrip, ab->path);
	}
	if (verbose) {
	    if (err) {
		pprintf(prn, " %s\n", _("no valid version found"));
	    } else {
		pprintf(prn, " %s %s (%s)\n", _("indexed version"),
			sys_ab.ok ? sys_ab.version : usr_ab.version,
			sys_ab.ok ? sys_ab.date : usr_ab.date);
	    }
	}

	/* clear for next addon */
	clear_addon_basics(&sys_ab);
	clear_addon_basics(&usr_ab);
    }

    fclose(fp);
    g_free(idxname);

    return 0;
}

/* Determine whether gretl has been updated since it was last run. We
   do this by comparing the build_date string in the program itself
   with that previously saved in the gretl config file.
*/

int gretl_is_updated (const char *prev_build)
{
    int b_curr, b_prev;
    int y, m, d;

    sscanf(BUILD_DATE, "%d-%d-%d", &y, &m, &d);
    b_curr = 10000*y + 100*m + d;

    sscanf(prev_build, "%d-%d-%d", &y, &m, &d);
    b_prev = 10000*y + 100*m + d;

    return b_curr > b_prev;
}

/* Get the full path to the .gfn file for a given addon. We first try
   for this via the simple plain text index file addons.idx. If that's
   not found we construct the index from scratch.

   The path that's returned is newly allocated and should be freed by
   the caller.
*/

char *gretl_addon_get_path (const char *addon)
{
    gchar *idxname = gretl_make_dotpath("addons.idx");
    FILE *fp = gretl_fopen(idxname, "rb");
    char *ret = NULL;

    if (fp == NULL) {
	update_addons_index(NULL);
	fp = gretl_fopen(idxname, "rb");
    }

    if (fp != NULL) {
	char *s, line[512];
	int nsp, n = strlen(addon);
	int got = 0;

	while (fgets(line, sizeof line, fp) && !got) {
	    if (!strncmp(addon, line, n)) {
		s = line;
		nsp = 0;
		while (*s && !got) {
		    if (*s == ' ') nsp++;
		    if (nsp == 2) {
			ret = gretl_strdup(s + 1);
			g_strchomp(ret);
			got = 1;
		    }
		    s++;
		}
	    }
	}
	fclose(fp);
    } else {
	fprintf(stderr, "failed to read addons.idx\n");
    }

    g_free(idxname);

    return ret;
}

/* Retrieve the path to an addon's "examples" sub-dir.  Note that it's
   not required that every addon has such.
*/

char *get_addon_examples_dir (const char *addon)
{
    char epath[MAXLEN];
    char *s, *path = gretl_addon_get_path(addon);
    char *ret = NULL;

    if (path != NULL) {
	s = strrslash(path);
	if (s != NULL) {
	    *s = '\0';
	}
	gretl_build_path(epath, path, "examples", NULL);
	if (g_file_test(epath, G_FILE_TEST_IS_DIR)) {
	    ret = gretl_strdup(epath);
	}
	free(path);
    }

    return ret;
}

/* Retrieve the path to an addon's PDF documentation, which is
   required of every addon. We accept the @addon argument with or
   without the ".pdf" suffix.
*/

char *get_addon_pdf_path (const char *addon)
{
    char *s, *path = NULL;
    char *ret = NULL;

    if (has_suffix(addon, ".pdf")) {
	/* strip off the suffix */
	gchar *tmp = g_strndup(addon, strlen(addon) - 4);

	path = gretl_addon_get_path(tmp);
	g_free(tmp);
    } else {
	path = gretl_addon_get_path(addon);
    }

    if (path != NULL && has_suffix(path, ".gfn")) {
	/* should be the case */
	s = strrchr(path, '.');
	*s = '\0';
	strcpy(s, ".pdf");
    }

    if (path != NULL && gretl_stat(path, NULL) == 0) {
	/* success */
	ret = path;
	path = NULL;
    }

    free(path);

    return ret;
}

int get_addon_basic_info (const char *addon,
			  char **version,
			  char **date,
			  char **descrip)
{
    gchar *idxname = gretl_make_dotpath("addons.idx");
    FILE *fp = gretl_fopen(idxname, "rb");
    int err = 0;

    if (fp == NULL) {
	err = update_addons_index(NULL);
	if (err) {
	    return err;
	} else {
	    fp = gretl_fopen(idxname, "rb");
	    if (fp == NULL) {
		err = E_FOPEN;
	    }
	}
    }

    if (fp != NULL) {
	char line[1024];
	char name[16];
	char verstr[16];
	char datestr[16];
	char descstr[32];
	int got = 0;

	while (fgets(line, sizeof line, fp) && !got) {
	    if (sscanf(line, "%s %s %s \"%31[^\"]", name, verstr,
		       datestr, descstr) == 3) {
		if (!strcmp(name, addon)) {
		    got = 1;
		    *version = gretl_strdup(verstr);
		    *date = gretl_strdup(datestr);
		    *descrip = gretl_strdup(descstr);
		}
	    }
	}
	fclose(fp);
	if (!got) {
	    fprintf(stderr, "addons.idx: couldn't find '%s'\n", addon);
	    err = E_DATA;
	}
    } else {
	fprintf(stderr, "failed to read addons.idx\n");
    }

    g_free(idxname);

    return err;
}
