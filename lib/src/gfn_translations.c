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

typedef struct msg_ {
    char *en;  /* English text */
    char *tr;  /* translation */
} msg;

struct Translation_ {
    char *lang;  /* language */
    int n_msgs;  /* number of messages handled */
    msg *msgs;   /* array of structs */
};

#if 0

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

        for (i=0; i<T->n_msgs; i++) {
            if (!strcmp(id, T->msgs[i].en)) {
                return T->msgs[i].tr;
            }
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

#endif

void destroy_translation (Translation *T)
{
    int i;

    free(T->lang);

    for (i=0; i<T->n_msgs; i++) {
        free(T->msgs[i].en);
        free(T->msgs[i].tr);
    }

    free(T->msgs);
    free(T);
}

static Translation *allocate_translation (char *lang, int n_msgs)
{
    Translation *T = malloc(sizeof *T);

    if (T != NULL) {
        int i;

        T->lang = lang;
        T->n_msgs = n_msgs;
        T->msgs = malloc(n_msgs * sizeof *T->msgs);
        for (i=0; i<n_msgs; i++) {
            T->msgs[i].en = NULL;
            T->msgs[i].tr = NULL;
        }
    }

    return T;
}

const char *get_gfn_translation (Translation *T,
                                 const char *id)
{
    static char *lang;

    if (lang == NULL) {
        lang = get_built_in_string_by_name("lang");
    }

    if (!strcmp(lang, T->lang) || !strncmp(lang, T->lang, 2)) {
        int i;

        for (i=0; i<T->n_msgs; i++) {
            if (!strcmp(id, T->msgs[i].en)) {
                return T->msgs[i].tr;
            }
        }
    }

    return id;
}

Translation *read_translation_element (xmlNodePtr root,
                                       xmlDocPtr doc)
{
    Translation *T = NULL;
    xmlNodePtr cur;
    int n_msgs = 0;
    char *lang = NULL;
    int i;

    if (!gretl_xml_get_prop_as_string(root, "lang", &lang)) {
        fprintf(stderr, "translation: @lang property is missing\n");
        return NULL;
    }

    /* first pass: count the @msg elements */

    cur = root->xmlChildrenNode;
    while (cur != NULL) {
        if (!xmlStrcmp(cur->name, (XUC) "msg")) {
            n_msgs++;
        }
        cur = cur->next;
    }

    if (n_msgs > 0) {
        T = allocate_translation(lang, n_msgs);
    }
    if (T == NULL) {
        return NULL;
    }

    /* second pass: populate the Translations struct */

    i = 0;
    cur = root->xmlChildrenNode;
    while (cur != NULL) {
        if (!xmlStrcmp(cur->name, (XUC) "msg")) {
            char *en = NULL;
            char *tr = NULL;

            if (gretl_xml_get_prop_as_string(cur, "en", &en) &&
                gretl_xml_node_get_string(cur, doc, &tr)) {
                if (*en != '\0' && *tr != '\0') {
                    T->msgs[i].en = en;
                    T->msgs[i].tr = tr;
                    i++;
                }
            }
        }
        cur = cur->next;
    }

    if (i == 0) {
        destroy_translation(T);
        T = NULL;
    } else if (i < n_msgs) {
        fprintf(stderr, "n_msgs = %d but only %d were usable\n", n_msgs, i);
        T->n_msgs = i;
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

void write_translation (Translation *T, PRN *prn)
{
    if (T != NULL) {
        msg *m;
        int i;

        pprintf(prn, "<translation lang=\"%s\">\n", T->lang);
        for (i=0; i<T->n_msgs; i++) {
            m = &T->msgs[i];
            pprintf(prn, " <msg en=\"%s\">%s</msg>\n", m->en, m->tr);
        }
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
