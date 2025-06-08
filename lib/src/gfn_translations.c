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

// static int verbose = 1;

struct Translation_ {
    char *lang;       /* language */
    GHashTable *msgs; /* table */
};

#define SUPPORT_MULTIPLE 0 /* not yet */

#if SUPPORT_MULTIPLE

/* Apparatus for dealing with multiple translations: not used yet */

struct Translations_ {
    int n_langs;     /* number of languages (other than English) */
    Translation **T; /* individual language translations */
    int active;      /* ID number of active translation */
};

typedef struct Translations_ Translations;

void destroy_translations (Translations *TT)
{
    int i;

    for (i=0; i<TT->n_langs; i++) {
        destroy_translation(TT->T[i]);
    }
    free(TT);
}

static Translations *allocate_translations (int n_langs)
{
    Translations *TT = malloc(sizeof *TT);

    if (TT != NULL) {
        TT->n_langs = n_langs;
        TT->T = malloc(n_langs * sizeof *TT->T);
        TT->active = -1;
    }

    return TT;
}

const char *get_gfn_translation (Translations *TT,
                                 const char *id)
{
    int i;

    if (TT->active < 0) {
        char *lang = get_built_in_string_by_name("lang");

        for (i=0; i<TT->n_langs; i++) {
            /* FIXME comparison? */
            if (!strcmp(lang, TT->T[i]->lang) ||
                !strncmp(lang, TT->T[i]->lang, 2)) {
                TT->active = i;
                break;
            }
        }
    }

    if (TT->active >= 0) {
        Translation *T = TT->T[TT->active];
        const char *tr = g_hash_table_lookup(T->msgs, id);

        if (tr != NULL) {
            return tr;
        }
    }

    return id;
}

Translations *read_translations_element (xmlNodePtr root,
                                         xmlDocPtr doc)
{
    return NULL; /* FIXME */
}

void write_translations (Translations *TT, PRN *prn)
{
    ; /* FIXME */
}

#endif /* SUPPORT_MULTIPLE */

void destroy_translation (Translation *T)
{
    free(T->lang);
    g_hash_table_destroy(T->msgs);
    free(T);
}

static Translation *allocate_translation (char *lang)
{
    Translation *T = malloc(sizeof *T);

    if (T != NULL) {
        T->lang = lang;
        T->msgs = g_hash_table_new_full(g_str_hash, g_str_equal,
                                        free, free);
    }

    return T;
}

const char *get_gfn_translation (Translation *T,
                                 const char *id)
{
    static char *lang;
    const char *ret = NULL;

    if (lang == NULL) {
        lang = get_built_in_string_by_name("lang");
    }

    if (!strcmp(lang, T->lang) || !strncmp(lang, T->lang, 2)) {
        ret = g_hash_table_lookup(T->msgs, id);
    }

    return ret != NULL ? ret : id;
}

Translation *read_translation_element (xmlNodePtr root,
                                       xmlDocPtr doc)
{
    Translation *T = NULL;
    xmlNodePtr cur;
    char *lang = NULL;
    int i;

    if (!gretl_xml_get_prop_as_string(root, "lang", &lang)) {
        fprintf(stderr, "translation: @lang property is missing\n");
        return NULL;
    }

    T = allocate_translation(lang);
    if (T == NULL) {
        fprintf(stderr, "translation: allocation failed\n");
        return NULL;
    }

    i = 0;
    cur = root->xmlChildrenNode;
    while (cur != NULL) {
        if (!xmlStrcmp(cur->name, (XUC) "msg")) {
            char *en = NULL;
            char *tr = NULL;

            if (gretl_xml_get_prop_as_string(cur, "en", &en)) {
                gretl_xml_node_get_string(cur, doc, &tr);
                g_hash_table_insert(T->msgs, en, tr);
                i++;
            }
        }
        cur = cur->next;
    }

    if (i == 0) {
        destroy_translation(T);
        T = NULL;
    }

    return T;
}

Translation *read_translations_file (const char *fname, int *err)
{
    Translation *T = NULL;
    xmlDocPtr doc = NULL;
    xmlNodePtr root;

    *err = gretl_xml_open_doc_root(fname, "translation", &doc, &root);
    if (!*err) {
        T = read_translation_element(root, doc);
    }

    if (doc != NULL) {
	xmlFreeDoc(doc);
    }

    return T;
}

static void merge_tr_msg (gpointer key,
                          gpointer value,
                          gpointer data)
{
    Translation *T0 = data;
    char *newkey, *newval;

    if (!g_hash_table_contains(T0->msgs, key)) {
        newkey = gretl_strdup((const char *) key);
        newval = gretl_strdup((const char *) value);
        g_hash_table_insert(T0->msgs, newkey, newval);
    } else {
        char *val0 = g_hash_table_lookup(T0->msgs, key);

        if (strcmp(val0, (const char *) value)) {
            newkey = gretl_strdup((const char *) key);
            newval = gretl_strdup((const char *) value);
            g_hash_table_replace(T0->msgs, newkey, newval);
        }
    }
}

Translation *update_translation (Translation *T0, const char *trbuf)
{
    Translation *ret = NULL;
    xmlDocPtr doc = NULL;
    xmlNodePtr root;
    int err;

    err = gretl_xml_read_buffer(trbuf, "translation", &doc, &root);
    if (err) {
        fprintf(stderr, "update_translation: error on reading buffer\n");
        return NULL;
    }

    if (T0 == NULL) {
        ret = read_translation_element(root, doc);
    } else {
        /* merge content of @trbuf into @T0 */
        Translation *T1 = read_translation_element(root, doc);

        g_hash_table_foreach(T1->msgs, merge_tr_msg, T0);
        destroy_translation(T1);
        ret = T0;
    }

    if (doc != NULL) {
	xmlFreeDoc(doc);
    }

    return ret;
}

static void write_tr_msg (gpointer key,
                          gpointer value,
                          gpointer data)
{
    pprintf((PRN *) data, " <msg en=\"%s\">%s</msg>\n",
            (const char *) key, (const char *) value);
}

void write_translation (Translation *T, PRN *prn)
{
    if (T != NULL) {
        pprintf(prn, "<translation lang=\"%s\">\n", T->lang);
        g_hash_table_foreach(T->msgs, write_tr_msg, prn);
        pputs(prn, "</translation>\n");
    }
}

/* Here we're looking at T_(...) and we want the "..." bit, unquoted */

gchar *get_translatable_content (const char **ps)
{
    gchar *ret = NULL;
    const char *s = *ps;
    const char *p;
    int quoted = 0;

    s += 3; /* skip "T_(" */
    if (*s == '"') {
        quoted = 1;
        s++;
    }
    p = s;
    while (*p) {
        if (quoted) {
            if (*p == '"') {
                break;
            }
        } else if (*p == ')') {
            break;
        }
        p++;
    }

    ret = g_strndup(s, p - s);

    if (*p == '"') {
        p++;
    }
    if (*p == ')') {
        p++;
    }
    *ps = p;

    return ret;
}
