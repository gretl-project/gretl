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

#include "libgretl.h"
#include "gretl_typemap.h"

struct type_map {
    GretlType std;     /* a "standard" type: non-array, non-reference */
    GretlType stdref;  /* the "reference" form of the standard type, if any */
    GretlType plural;  /* the associated array type, if any */
    GretlType plref;   /* reference form of array type, if any */
};

static struct type_map gretl_type_map[] = {
    { GRETL_TYPE_MATRIX,   GRETL_TYPE_MATRIX_REF,
      GRETL_TYPE_MATRICES, GRETL_TYPE_MATRICES_REF},
    { GRETL_TYPE_BUNDLE,   GRETL_TYPE_BUNDLE_REF,
      GRETL_TYPE_BUNDLES,  GRETL_TYPE_BUNDLES_REF},
    { GRETL_TYPE_STRING,   GRETL_TYPE_STRING_REF,
      GRETL_TYPE_STRINGS,  GRETL_TYPE_STRINGS_REF},
    { GRETL_TYPE_ARRAY,    GRETL_TYPE_ARRAY_REF,
      GRETL_TYPE_ARRAYS,   GRETL_TYPE_ARRAYS_REF},
    { GRETL_TYPE_LIST,     GRETL_TYPE_LIST_REF,
      GRETL_TYPE_LISTS,    GRETL_TYPE_LISTS_REF},
    { GRETL_TYPE_SERIES,   GRETL_TYPE_SERIES_REF, 0, 0},
    { GRETL_TYPE_DOUBLE,   GRETL_TYPE_SCALAR_REF, 0, 0}
};

static int N_TYPES = G_N_ELEMENTS(gretl_type_map);

GretlType gretl_type_get_plural (GretlType type)
{
    int i;

    if (type == 0) return 0;

    for (i=0; i<N_TYPES; i++) {
	if (type == gretl_type_map[i].std) {
	    return gretl_type_map[i].plural;
	}
    }

    return 0;
}

GretlType gretl_type_get_singular (GretlType type)
{
    int i;

    if (type == 0) return 0;

    for (i=0; i<N_TYPES; i++) {
	if (type == gretl_type_map[i].plural) {
	    return gretl_type_map[i].std;
	}
    }

    return 0;
}

GretlType gretl_type_get_ref_type (GretlType type)
{
    int i;

    if (type == 0) return 0;

    for (i=0; i<N_TYPES; i++) {
	if (type == gretl_type_map[i].std) {
	    return gretl_type_map[i].stdref;
	}
	if (type == gretl_type_map[i].plural) {
	    return gretl_type_map[i].plref;
	}
	if (type == gretl_type_map[i].stdref ||
	    type == gretl_type_map[i].plref) {
	    /* allow pass-through */
	    return type;
	}
    }

    return 0;
}

GretlType gretl_type_get_plain_type (GretlType type)
{
    int i;

    if (type == 0) return 0;

    for (i=0; i<N_TYPES; i++) {
	if (type == gretl_type_map[i].stdref) {
	    return gretl_type_map[i].std;
	}
	if (type == gretl_type_map[i].plref) {
	    return gretl_type_map[i].plural;
	}
	if (type == gretl_type_map[i].std ||
	    type == gretl_type_map[i].plural) {
	    /* allow pass-through */
	    return type;
	}
    }

    return 0;
}

int gretl_is_array_type (GretlType type)
{
    int i;

    if (type == 0) return 0;

    for (i=0; i<N_TYPES; i++) {
	if (type == gretl_type_map[i].plural) {
	    return 1;
	}
    }

    return 0;
}

int gretl_is_array_ref_type (GretlType type)
{
    int i;

    if (type == 0) return 0;

    for (i=0; i<N_TYPES; i++) {
	if (type == gretl_type_map[i].plref) {
	    return 1;
	}
    }

    return 0;
}

int gretl_is_arrayable_type (GretlType type)
{
    return gretl_type_get_plural(type) > 0;
}

/**
 * gretl_type_get_name:
 * @type: a gretl type.
 *
 * Returns: the name of @type, or "invalid" on failure.
 */

const char *gretl_type_get_name (GretlType type)
{
    switch (type) {
    case GRETL_TYPE_BOOL:       return "bool";
    case GRETL_TYPE_INT:        return "int";
    case GRETL_TYPE_UINT32:     return "unsigned32";
    case GRETL_TYPE_UINT64:     return "unsigned64";
    case GRETL_TYPE_OBS:        return "obs";
    case GRETL_TYPE_DOUBLE:     return "scalar";
    case GRETL_TYPE_SERIES:     return "series";
    case GRETL_TYPE_USERIES:    return "series";
    case GRETL_TYPE_MATRIX:     return "matrix";
    case GRETL_TYPE_LIST:       return "list";
    case GRETL_TYPE_BUNDLE:     return "bundle";
    case GRETL_TYPE_ARRAY:      return "array";
    case GRETL_TYPE_STRING:     return "string";
    case GRETL_TYPE_SCALAR_REF: return "scalar *";
    case GRETL_TYPE_SERIES_REF: return "series *";
    case GRETL_TYPE_MATRIX_REF: return "matrix *";
    case GRETL_TYPE_BUNDLE_REF: return "bundle *";
    case GRETL_TYPE_ARRAY_REF:  return "array *";
    case GRETL_TYPE_STRING_REF: return "string *";
    case GRETL_TYPE_LIST_REF:   return "list *";

    case GRETL_TYPE_STRINGS:      return "strings";
    case GRETL_TYPE_MATRICES:     return "matrices";
    case GRETL_TYPE_BUNDLES:      return "bundles";
    case GRETL_TYPE_LISTS:        return "lists";
    case GRETL_TYPE_ARRAYS:       return "arrays";

    case GRETL_TYPE_STRINGS_REF:  return "strings *";
    case GRETL_TYPE_MATRICES_REF: return "matrices *";
    case GRETL_TYPE_BUNDLES_REF:  return "bundles *";
    case GRETL_TYPE_LISTS_REF:    return "lists *";
    case GRETL_TYPE_ARRAYS_REF:   return "arrays *";

    case GRETL_TYPE_DATE:         return "date"; /* ODBC special */

    case GRETL_TYPE_VOID:       return "void";
    case GRETL_TYPE_NUMERIC:    return "numeric";
    case GRETL_TYPE_NONE:       return "null";
    case GRETL_TYPE_ANY:        return "any";
    default:
	return "invalid";
    }
}

GretlType gretl_type_from_string (const char *s)
{
    const char *p;

    if (!strncmp(s, "matrix", 6)) {
	p = s + 6;
	if (*p == '\0') {
	    return GRETL_TYPE_MATRIX;
	} else if (!strcmp(p, " *") || !strcmp(p, "ref")) {
	    return GRETL_TYPE_MATRIX_REF;
	}
    } else if (!strncmp(s, "series", 6)) {
	p = s + 6;
	if (*p == '\0') {
	    return GRETL_TYPE_SERIES;
	} else if (!strcmp(p, " *") || !strcmp(p, "ref")) {
	    return GRETL_TYPE_SERIES_REF;
	}
    } else if (!strncmp(s, "scalar", 6)) {
	p = s + 6;
	if (*p == '\0') {
	    return GRETL_TYPE_DOUBLE;
	} else if (!strcmp(p, " *") || !strcmp(p, "ref")) {
	    return GRETL_TYPE_SCALAR_REF;
	}
    } else if (!strncmp(s, "bundle", 6)) {
	p = s + 6;
	if (*p == 's') {
	    p++;
	    if (*p == '\0') {
		return GRETL_TYPE_BUNDLES;
	    } else if (!strcmp(p, " *") || !strcmp(p, "ref")) {
		return GRETL_TYPE_BUNDLES_REF;
	    }
	} else if (*p == '\0') {
	    return GRETL_TYPE_BUNDLE;
	} else if (!strcmp(p, " *") || !strcmp(p, "ref")) {
	    return GRETL_TYPE_BUNDLE_REF;
	}
    } else if (!strncmp(s, "array", 5)) {
	p = s + 5;
	if (*p == 's') {
	    p++;
	    if (*p == '\0') {
		return GRETL_TYPE_ARRAYS;
	    } else if (!strcmp(p, " *") || !strcmp(p, "ref")) {
		return GRETL_TYPE_ARRAYS_REF;
	    }
	} else if (*p == '\0') {
	    return GRETL_TYPE_ARRAY;
	} else if (!strcmp(p, " *") || !strcmp(p, "ref")) {
	    return GRETL_TYPE_ARRAY_REF;
	}
    } else if (!strncmp(s, "string", 6)) {
	p = s + 6;
	if (*p == 's') {
	    p++;
	    if (*p == '\0') {
		return GRETL_TYPE_STRINGS;
	    } else if (!strcmp(p, " *") || !strcmp(p, "ref")) {
		return GRETL_TYPE_STRINGS_REF;
	    }
	} else if (*p == '\0') {
	    return GRETL_TYPE_STRING;
	} else if (!strcmp(p, " *") || !strcmp(p, "ref")) {
	    return GRETL_TYPE_STRING_REF;
	}
    } else if (!strncmp(s, "matrices", 8)) {
	p = s + 8;
	if (*p == '\0') {
	    return GRETL_TYPE_MATRICES;
	} else if (!strcmp(p, " *") || !strcmp(p, "ref")) {
	    return GRETL_TYPE_MATRICES_REF;
	}
    } else if (!strncmp(s, "list", 4)) {
	p = s + 4;
	if (*p == 's') {
	    p++;
	    if (*p == '\0') {
		return GRETL_TYPE_LISTS;
	    } else if (!strcmp(p, " *") || !strcmp(p, "ref")) {
		return GRETL_TYPE_LISTS_REF;
	    }
	} else if (*p == '\0') {
	    return GRETL_TYPE_LIST;
	} else if (!strcmp(p, " *") || !strcmp(p, "ref")) {
	    return GRETL_TYPE_LIST_REF;
	}
    } else if (!strncmp(s, "numeric", 7)) {
	return GRETL_TYPE_NUMERIC;
    } else {
	/* aliases */
	if (!strcmp(s, "bool"))     return GRETL_TYPE_BOOL;
	if (!strcmp(s, "boolean"))  return GRETL_TYPE_BOOL;
	if (!strcmp(s, "int"))      return GRETL_TYPE_INT;
	if (!strcmp(s, "unsigned")) return GRETL_TYPE_UINT32;
	if (!strcmp(s, "obs"))      return GRETL_TYPE_OBS;
    }

    return GRETL_TYPE_NONE;
}

struct lookup {
    const char *word;
    GretlType type;
};

static const struct lookup gentypes[] = {
    { "series",   GRETL_TYPE_SERIES },
    { "scalar",   GRETL_TYPE_DOUBLE },
    { "matrix",   GRETL_TYPE_MATRIX },
    { "string",   GRETL_TYPE_STRING },
    { "list",     GRETL_TYPE_LIST },
    { "bundle",   GRETL_TYPE_BUNDLE },
    { "strings",  GRETL_TYPE_STRINGS },
    { "matrices", GRETL_TYPE_MATRICES },
    { "bundles",  GRETL_TYPE_BUNDLES },
    { "lists",    GRETL_TYPE_LISTS },
    { "arrays",   GRETL_TYPE_ARRAYS },
};

GretlType gretl_get_gen_type (const char *s)
{
    static GHashTable *ht;
    gpointer p;
    GretlType t = 0;

    if (s == NULL) {
	/* the clean-up signal */
	if (ht != NULL) {
	    g_hash_table_destroy(ht);
	    ht = NULL;
	}
	return 0;
    }

    if (ht == NULL) {
	int i, n = G_N_ELEMENTS(gentypes);

	ht = g_hash_table_new(g_str_hash, g_str_equal);
	for (i=0; i<n; i++) {
	    g_hash_table_insert(ht, (gpointer) gentypes[i].word,
				GINT_TO_POINTER(gentypes[i].type));
	}
    }

    p = g_hash_table_lookup(ht, s);
    if (p != NULL) {
	t = GPOINTER_TO_INT(p);
    }

    return t;
}

/* Note: this must agree with the doc for the typeof() function */

int gretl_type_get_order (GretlType type)
{
    if (gretl_scalar_type(type)) {
        return 1;
    } else if (type == GRETL_TYPE_SERIES) {
        return 2;
    } else if (type == GRETL_TYPE_MATRIX) {
        return 3;
    } else if (type == GRETL_TYPE_STRING) {
        return 4;
    } else if (type == GRETL_TYPE_BUNDLE) {
        return 5;
    } else if (type == GRETL_TYPE_ARRAY) {
        return 6;
    } else if (type == GRETL_TYPE_LIST) {
        return 7;
    } else {
        return 0;
    }
}

int gretl_is_scalar_type (GretlType type)
{
    return type == GRETL_TYPE_BOOL ||
	type == GRETL_TYPE_INT ||
	type == GRETL_TYPE_UINT32 ||
        type == GRETL_TYPE_UINT64 ||
	type == GRETL_TYPE_DOUBLE;
}

int gretl_is_series_type (GretlType type)
{
    return type == GRETL_TYPE_SERIES ||
	type == GRETL_TYPE_USERIES;
}

void gretl_typemap_cleanup (void)
{
    gretl_get_gen_type(NULL);
}

gchar *name_conflict_message (const char *name, GretlType type)
{
    gchar *s = NULL;

    if (type == GRETL_TYPE_SERIES) {
	s = g_strdup_printf(_("A series named %s already exists"), name);
    } else if (type == GRETL_TYPE_MATRIX) {
	s = g_strdup_printf(_("A matrix named %s already exists"), name);
    } else if (type == GRETL_TYPE_DOUBLE) {
	s = g_strdup_printf(_("A scalar named %s already exists"), name);
    } else if (type == GRETL_TYPE_LIST) {
	s = g_strdup_printf(_("A list named %s already exists"), name);
    } else if (type == GRETL_TYPE_STRING) {
	s = g_strdup_printf(_("A string named %s already exists"), name);
    } else if (type == GRETL_TYPE_BUNDLE) {
	s = g_strdup_printf(_("A bundle named %s already exists"), name);
    } else if (type == GRETL_TYPE_ARRAY) {
	s = g_strdup_printf(_("An array named %s already exists"), name);
    }

    return s;
}

/* Used in gui/fncall.c, where (until it becomes time to initiate an
   actual function call) we don't make a hard distinction between
   plain and pointer types.
*/

int gretl_type_mismatch (GretlType t1, GretlType t2)
{
    if (t1 == t2) {
	return 0;
    } else {
	t1 = gretl_type_get_plain_type(t1);
	t2 = gretl_type_get_plain_type(t2);
	if (t1 == t2) {
	    return 0;
	}
    }

    return 1;
}
