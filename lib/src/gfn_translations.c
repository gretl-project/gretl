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
#include "gretl_xml.h"
#include "gfn_translations.h"
#include "gretl_string_table.h"

static int verbose = 0;

typedef struct msg_ {
    char *id;     /* English text */
    char **strs;  /* One or more translations */
} msg;

struct Translations_ {
    int n_msgs;   /* number of messages handled */
    int n_langs;  /* number of languages (other than English) */
    msg *msgs;    /* array of structs */
    char **langs; /* array of labuage IDs */
};

void destroy_translations (Translations *T)
{
    int i, j;

    for (i=0; i<T->n_langs; i++) {
        free(T->langs[i]);
    }
    free(T->langs);

    for (i=0; i<T->n_msgs; i++) {
        free(T->msgs[i].id);
        for (j=0; j<T->n_langs; j++) {
            free(T->msgs[i].strs[j]);
        }
        free(T->msgs[i].strs);
    }
    free(T->msgs);

    free(T);
}

static Translations *allocate_translations (int n_msgs, int n_langs)
{
    Translations *T = malloc(sizeof *T);

    if (T != NULL) {
        int i;

        T->n_msgs = n_msgs;
        T->n_langs = n_langs;
        T->msgs = malloc(n_msgs * sizeof *T->msgs);
        for (i=0; i<n_msgs; i++) {
            T->msgs[i].strs = NULL;
        }
        T->langs = NULL;
    }

    return T;
}

const char *get_gfn_translation (Translations *T,
                                 const char *id)
{
    char *lang = get_built_in_string_by_name("lang");
    int i, L = -1;

    for (i=0; i<T->n_langs && L < 0; i++) {
        /* FIXME? */
        if (!strcmp(lang, T->langs[i]) ||
            !strncmp(lang, T->langs[i], 2)) {
            L = i;
        }
    }

    if (L >= 0) {
        for (i=0; i<T->n_msgs; i++) {
            if (!strcmp(id, T->msgs[i].id)) {
                return T->msgs[i].strs[L];
            }
        }
    }

    return id;
}

/* This is first called with @n_langs non-NULL but @langs and @strs
   NULL, just to get a count of the elements to be stored.

   It is then called with @langs and @strs non-NULL to record the
   language IDs and translated strings.
*/

static int msg_get_data (xmlNodePtr np, xmlDocPtr doc, int i,
                         int *n_langs, char **langs, char **strs)
{
    xmlNodePtr cur = np->xmlChildrenNode;
    char *lang;
    char *str;
    int j = 0;
    int err = 0;

    while (cur != NULL && !err) {
        if (!xmlStrcmp(cur->name, (XUC) "s")) {
            lang = str = NULL;
            if (!gretl_xml_get_prop_as_string(cur, "lang", &lang)) {
                err = E_DATA;
            } else if (!gretl_xml_node_get_string(cur, doc, &str)) {
                err = E_DATA;
            } else if (langs != NULL && strs != NULL) {
                if (i == 0) {
                    /* attach lang string */
                    langs[j] = lang;
                    if (verbose) {
                        printf("attach: lang[%d] = %s\n", j, langs[j]);
                    }
                } else {
                    /* just checking langs */
                    if (strcmp(lang, langs[j])) {
                        printf(" language ids are not consistent\n");
                    }
                    free(lang);
                }
                strs[j] = str;
            } else {
                /* just counting */
                if (verbose) {
                    printf(" msg: lang '%s', str '%s'\n", lang, str);
                }
                *n_langs += 1;
                free(lang);
                free(str);
            }
            j++;
        }
        cur = cur->next;
    }

    return err;
}

Translations *read_translations_element (xmlNodePtr root,
                                         xmlDocPtr doc)
{
    Translations *T = NULL;
    xmlNodePtr cur;
    int n_msgs = 0;
    int n_langs = 0;
    char *id = NULL;
    int i, nli;
    int err = 0;

    /* first pass through the XML */

    cur = root->xmlChildrenNode;
    i = 0;
    while (cur != NULL && !err) {
        if (!xmlStrcmp(cur->name, (XUC) "msg")) {
            nli = 0;
            if (gretl_xml_get_prop_as_string(cur, "id", &id)) {
                ++n_msgs;
                if (verbose) {
                    printf("id %d: '%s'\n", n_msgs, id);
                }
                err = msg_get_data(cur, doc, i, &nli, NULL, NULL);
                if (n_langs == 0) {
                    n_langs = nli;
                } else if (nli != n_langs) {
                    printf(" number of_langs is not consistent\n");
                    err = E_DATA;
                }
                free(id);
            } else {
                err = E_DATA;
            }
            i++;
        }
        cur = cur->next;
    }

    if (!err && (n_msgs == 0 || n_langs == 0)) {
        err = E_DATA;
    }

    /* second pass */

    if (!err) {
        char **langs = calloc(n_langs, sizeof *langs);
        int i;

        if (verbose) {
            printf("Got %d msg elements, each with %d langs\n", n_msgs, n_langs);
        }
        T = allocate_translations(n_msgs, n_langs);

        cur = root->xmlChildrenNode;
        for (i=0; i<n_msgs && !err; i++) {
            if (!xmlStrcmp(cur->name, (XUC) "msg")) {
                char **strs = calloc(n_langs, sizeof *strs);

                gretl_xml_get_prop_as_string(cur, "id", &id);
                err = msg_get_data(cur, doc, i, NULL, langs, strs);
                if (!err) {
                    T->msgs[i].id = id;
                    T->msgs[i].strs = strs;
                }
            }
            cur = cur->next;
        }
        if (!err) {
            T->langs = langs;
        }
    }

    return T;
}

Translations *read_translations_file (const char *fname, int *err)
{
    Translations *T = NULL;
    xmlDocPtr doc = NULL;
    xmlNodePtr root;

    *err = gretl_xml_open_doc_root(fname, "translations", &doc, &root);
    if (!*err) {
        T = read_translations_element(root, doc);
    }

    if (doc != NULL) {
	xmlFreeDoc(doc);
    }

    return T;
}

void write_translations (Translations *T, PRN *prn)
{
    if (T != NULL) {
        msg *m;
        int i, j;

        pputs(prn, "<translations>\n");
        for (i=0; i<T->n_msgs; i++) {
            m = &T->msgs[i];
            pprintf(prn, " <msg id=\"%s\">\n", m->id);
            for (j=0; j<T->n_langs; j++) {
                pprintf(prn, "  <s lang=\"%s\">%s</s>\n", T->langs[j],
                        m->strs[j]);
            }
            pputs(prn, " </msg>\n");
        }
        pputs(prn, "</translations>\n");
    }
}
