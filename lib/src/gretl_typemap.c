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
    { GRETL_TYPE_STRING,   0,
      GRETL_TYPE_STRINGS,  GRETL_TYPE_STRINGS_REF},
    { GRETL_TYPE_LIST,     0, 
      GRETL_TYPE_LISTS,    GRETL_TYPE_LISTS_REF},
    { GRETL_TYPE_SERIES,   GRETL_TYPE_SERIES_REF, 0, 0},
    { GRETL_TYPE_DOUBLE,   GRETL_TYPE_SCALAR_REF, 0, 0},
};

GretlType gretl_type_get_plural (GretlType type)
{
    int i, n = G_N_ELEMENTS(gretl_type_map);

    if (type == 0) return 0;

    for (i=0; i<n; i++) {
	if (type == gretl_type_map[i].std) {
	    return gretl_type_map[i].plural;
	}
    }

    return 0;
}

GretlType gretl_type_get_singular (GretlType type)
{
    int i, n = G_N_ELEMENTS(gretl_type_map);

    if (type == 0) return 0;

    for (i=0; i<n; i++) {
	if (type == gretl_type_map[i].plural) {
	    return gretl_type_map[i].std;
	}
    }

    return 0;
}

GretlType gretl_type_get_ref_type (GretlType type)
{
    int i, n = G_N_ELEMENTS(gretl_type_map);

    if (type == 0) return 0;

    for (i=0; i<n; i++) {
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
    int i, n = G_N_ELEMENTS(gretl_type_map);

    if (type == 0) return 0;

    for (i=0; i<n; i++) {
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
    case GRETL_TYPE_OBS:        return "obs";
    case GRETL_TYPE_DOUBLE:     return "scalar";
    case GRETL_TYPE_SERIES:     return "series";
    case GRETL_TYPE_USERIES:    return "series";
    case GRETL_TYPE_MATRIX:     return "matrix";	
    case GRETL_TYPE_LIST:       return "list";
    case GRETL_TYPE_BUNDLE:     return "bundle";
    case GRETL_TYPE_ARRAY:      return "array";
    case GRETL_TYPE_SCALAR_REF: return "scalar *";
    case GRETL_TYPE_SERIES_REF: return "series *";
    case GRETL_TYPE_MATRIX_REF: return "matrix *";
    case GRETL_TYPE_BUNDLE_REF: return "bundle *";
    case GRETL_TYPE_STRING:     return "string";

    case GRETL_TYPE_STRINGS:      return "strings";
    case GRETL_TYPE_MATRICES:     return "matrices";
    case GRETL_TYPE_BUNDLES:      return "bundles";
    case GRETL_TYPE_LISTS:        return "lists";
	
    case GRETL_TYPE_STRINGS_REF:  return "strings *";
    case GRETL_TYPE_MATRICES_REF: return "matrices *";
    case GRETL_TYPE_BUNDLES_REF:  return "bundles *";
    case GRETL_TYPE_LISTS_REF:    return "lists *";

    case GRETL_TYPE_VOID:       return "void";
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
	}
    } else {
	/* aliases */
	if (!strcmp(s, "bool"))     return GRETL_TYPE_BOOL;
	if (!strcmp(s, "boolean"))  return GRETL_TYPE_BOOL;
	if (!strcmp(s, "int"))      return GRETL_TYPE_INT;
	if (!strcmp(s, "obs"))      return GRETL_TYPE_OBS;
	if (!strcmp(s, "array"))    return GRETL_TYPE_ARRAY;
    }

    return GRETL_TYPE_NONE;
}

