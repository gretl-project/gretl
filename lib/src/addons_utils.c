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
    "dbnomics", "extra", "geoplot", NULL
};

/* Determine if @pkgname is the name of an addon:
   @pkgname may be given with or without the ".gfn"
   suffix.
*/

int is_gretl_addon (const char *pkgname)
{
    int i;

    if (has_suffix(pkgname, ".gfn")) {
	int n = strlen(pkgname) - 4;

	for (i=0; addon_names[i] != NULL; i++) {
	    if (!strncmp(pkgname, addon_names[i], n)) {
		return 1;
	    }
	}
    } else {
	for (i=0; addon_names[i] != NULL; i++) {
	    if (!strcmp(pkgname, addon_names[i])) {
		return 1;
	    }
	}
    }

    return 0;
}

/* Return a NULL-terminated array of addons names. */

const char **get_addon_names (int *n)
{
    if (n != NULL) {
	int i;

	*n = 0;
	for (i=0; addon_names[i] != NULL; i++) {
	    *n = *n + 1;
	}
    }
    return addon_names;
}

/* Given a full path to an addon package, return its
   version string. If @data is non-NULL, also return its
   release date via this pointer. In both cases the
   strings are allocated, and belong to the caller.
*/

char *get_addon_version (const char *fname, char **date)
{
    char *version = NULL;
    xmlDocPtr doc = NULL;
    xmlNodePtr node, n1, n2;
    int targ, got = 0;

    if (gretl_stat(fname, NULL) != 0) {
	/* not found */
	return NULL;
    }

    gretl_xml_open_doc_root(fname, "gretl-functions", &doc, &node);
    if (doc == NULL || node == NULL) {
	return NULL;
    }

    targ = (date != NULL)? 2 : 1;
    n1 = node->xmlChildrenNode;

    while (n1 != NULL && got < targ) {
	if (!xmlStrcmp(n1->name, (XUC) "gretl-function-package")) {
	    n2 = n1->xmlChildrenNode;
	    while (n2 != NULL && got < targ) {
		if (!xmlStrcmp(n2->name, (XUC) "version")) {
		    gretl_xml_node_get_trimmed_string(n2, doc, &version);
		    got++;
		} else if (date != NULL && !xmlStrcmp(n2->name, (XUC) "date")) {
		    gretl_xml_node_get_trimmed_string(n2, doc, date);
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

    return version;
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

int update_addons_index (void)
{
    gchar *idxname = gretl_make_dotpath("addons.idx");
    char *pkgver1 = NULL;
    char *pkgver2 = NULL;
    char syspath[MAXLEN];
    char usrpath[MAXLEN];
    char gfnname[64];
    double v1, v2;
    FILE *fp;
    int i;

    fp = gretl_fopen(idxname, "wb");
    if (fp == NULL) {
	g_free(idxname);
	return E_FOPEN;
    }

    for (i=0; addon_names[i] != NULL; i++) {
	sprintf(gfnname, "%s.gfn", addon_names[i]);
	gretl_build_path(syspath, gretl_home(), "functions",
			 addon_names[i], gfnname, NULL);
	pkgver1 = get_addon_version(syspath, NULL);
	v1 = (pkgver1 != NULL)? dot_atof(pkgver1) : 0;
	get_user_path(usrpath, addon_names[i], gfnname);
	pkgver2 = get_addon_version(usrpath, NULL);
	v2 = (pkgver2 != NULL)? dot_atof(pkgver2) : 0;
	if (v1 > v2) {
	    fprintf(fp, "%s %s %s\n", addon_names[i], pkgver1, syspath);
	} else if (v2 > 0) {
	    fprintf(fp, "%s %s %s\n", addon_names[i], pkgver2, usrpath);
	}
	free(pkgver1);
	free(pkgver2);
    }

    fclose(fp);
    g_free(idxname);

    return 0;
}

static int gretl_is_updated (const char *prev_build)
{
    int b_curr, b_prev;
    int y, m, d;

    sscanf(BUILD_DATE, "%d-%d-%d", &y, &m, &d);
    b_curr = 10000*y + 100*m + d;

    sscanf(prev_build, "%d-%d-%d", &y, &m, &d);
    b_prev = 10000*y + 100*m + d;

    return b_curr > b_prev;
}

/* On starting gretl, see if we're updated and if so
   try rebuilding the addons index.
*/

int maybe_update_addons_index (const char *prev_build)
{
    if (gretl_is_updated(prev_build)) {
	return update_addons_index();
    }
    return 0;
}

char *gretl_addon_get_path (const char *addon)
{
    gchar *idxname = gretl_make_dotpath("addons.idx");
    FILE *fp = gretl_fopen(idxname, "rb");
    char *ret = NULL;

    if (fp == NULL) {
	update_addons_index();
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
    }

    g_free(idxname);

    return ret;
}

char *get_addon_examples_dir (const char *addon)
{
    char epath[MAXLEN];
    char *s, *path = gretl_addon_get_path(addon);
    char *ret = NULL;

    if (path != NULL) {
	s = strrchr(path, SLASH);
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
