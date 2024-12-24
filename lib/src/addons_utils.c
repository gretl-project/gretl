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

#define AU_DEBUG 0

/* Note: here's the canonical listing of gretl addons. If a new addon
   is added this must be updated
*/

static const char *addon_names[] = {
    "SVAR", "gig", "HIP", "ivpanel",
    "dbnomics", "extra", "geoplot",
    "regls", "logging", "KFgui", NULL
};

static int n_addons = (G_N_ELEMENTS(addon_names)) - 1;

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

/* Determine if @pkgname is the name of an addon: @pkgname may be
   given with or without the ".gfn" suffix.
*/

int is_gretl_addon (const char *pkgname)
{
    int i, n = strlen(pkgname);

    if (has_suffix(pkgname, ".gfn")) {
        n -= 4;
    }

    for (i=0; i<n_addons; i++) {
        if (!strncmp(pkgname, addon_names[i], n)) {
            return 1;
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
        *n = n_addons;
    }
    return addon_names;
}

/* Given the full path to an addon gfn file, supply its version, date
   and description. The strings are newly allocated and belong to the
   caller. Set the "ok" member of @ab to 1 if successful.
*/

static void get_addon_basics_from_gfn (addon_basics *ab)
{
    xmlDocPtr doc = NULL;
    xmlNodePtr node, n1, n2;
    int targ = 3;
    int got = 0;

    if (gretl_stat(ab->path, NULL) != 0) {
        return;
    }

    gretl_xml_open_doc_root(ab->path, "gretl-functions", &doc, &node);
    if (doc == NULL || node == NULL) {
        return;
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

    ab->ok = (got == targ);
}

static char *get_user_path (char *targ, const char *pgkname,
                            const char *gfnname)
{
#ifdef __APPLE__
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

    if (sys_ab->ok == 0 && usr_ab->ok == 0) {
        err = E_DATA;
    } else if (sys_ab->ok && usr_ab->ok) {
        /* both match gretl version, compare dates */
        guint32 sed = get_epoch_day(sys_ab->date);
        guint32 ued = get_epoch_day(usr_ab->date);

        if (sed >= ued) {
            /* sys more recent: invalidate usr */
            usr_ab->ok = 0;
        } else {
            /* usr more recent: invalidate sys */
            sys_ab->ok = 0;
        }
    }

    return err;
}

static void report_addon_result (addon_basics *ab, PRN *prn)
{
    pprintf(prn, " %s '%s'\n", _("try"), ab->path);
    if (ab->ok) {
        pprintf(prn, "  %s %s (%s)\n", _("found version"),
                ab->version, ab->date);
    } else {
        pprintf(prn, "  %s\n", _("not found"));
    }
}

/* Build a plain text index of the installed addons, holding name,
   version, date and full path, one addon per line. We allow for the
   fact that a given addon might exist both in the "system" location
   and in the user's personal filespace.  In case of such duplicates
   we determine which version is newer and enter its details in the
   index file.

   This update routine is called automatically when (a) a user
   installs an addons package (addons.tar.gz) or (b) gretl figures out
   that it has been updated (new release or snapshot). It can also be
   called explicitly by the user, via the command "pkg index addons".

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
    FILE *fp;
    int i;

#if AU_DEBUG
    fprintf(stderr, "*** update_addons_index called ***\n");
#endif

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
        get_addon_basics_from_gfn(&sys_ab);
        if (verbose) {
            report_addon_result(&sys_ab, prn);
        }

        /* (2) build and check the userspace path */
        get_user_path(usr_ab.path, addon_names[i], gfnname);
        get_addon_basics_from_gfn(&usr_ab);
        if (verbose) {
            report_addon_result(&usr_ab, prn);
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
    int err = 0;
    PRN *auxprn;
    
#if AU_DEBUG
    fprintf(stderr, "gretl_addon_get_path: '%s'\n", addon);
    auxprn = gretl_print_new(GRETL_PRINT_STDERR, &err);
#else
    auxprn = NULL;
#endif

    if (fp == NULL) {
        err = update_addons_index(auxprn);
        if (!err) {
            fp = gretl_fopen(idxname, "rb");
        }
    }

    gretl_print_destroy(auxprn);
    
    if (fp == NULL) {
        fprintf(stderr, "failed to read addons.idx\n");
    } else {
        char *s, line[1024];
        int n = strlen(addon);
        int nq = 0;

        while (fgets(line, sizeof line, fp)) {
            if (!strncmp(addon, line, n)) {
                char *p = NULL;
                char *q = NULL;

                s = line + n;
                while (*s) {
                    if (*s == '"') {
                        nq++;
                        if (nq == 3) {
                            p = s + 1;
                        } else if (nq == 4) {
                            q = s;
                        }
                    }
                    s++;
                }
                if (p != NULL && q != NULL) {
                    ret = gretl_strndup(p, q - p);
                }
            }
        }
        fclose(fp);
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

/* called from gui/database.c only */

int get_addon_basic_info (const char *addon,
                          char **date,
                          char **descrip,
                          char **fname)
{
    gchar *idxname = gretl_make_dotpath("addons.idx");
    FILE *fp = gretl_fopen(idxname, "rb");
    int err = 0;

    if (fp == NULL) {
        err = update_addons_index(NULL);
        if (!err) {
            fp = gretl_fopen(idxname, "rb");
        }
    }

    if (fp == NULL) {
        err = E_FOPEN;
    } else {
        char line[1024];
        char name[16];
        char dstr[16];
        char desc[64];
        char path[MAXLEN];
        int got = 0;

        while (fgets(line, sizeof line, fp) && !got) {
            if (sscanf(line, "%s %*s %s \"%63[^\"]\" \"%511[^\"]",
                       name, dstr, desc, path) == 4) {
                if (!strcmp(name, addon) &&
                    gretl_test_fopen(path, "r") == 0) {
                    got = 1;
                    *date = gretl_strdup(dstr);
                    *descrip = gretl_strdup(desc);
                    *fname = gretl_strdup(path);
                }
            }
        }
        fclose(fp);
        if (!got) {
            fprintf(stderr, "addons.idx: couldn't find '%s'\n", addon);
            err = E_DATA;
        }
    }

    g_free(idxname);

    return err;
}

