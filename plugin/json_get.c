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

/* parsing of JSON buffer using the json-glib library */

#include "libgretl.h"
#include "version.h"
#include "gretl_typemap.h"

#include <glib/gprintf.h>
#include <glib-object.h>
#include <json-glib/json-glib.h>

#define handled_type(t) (t == G_TYPE_STRING || \
			 t == G_TYPE_DOUBLE || \
	                 t == G_TYPE_BOOLEAN || \
			 t == G_TYPE_INT64)

#define non_empty_array(a) (a != NULL && json_array_get_length(a) > 0)
#define null_node(n) (n == NULL || json_node_is_null(n))

/* We don't want to leak memory allocated for @node, but it
   should not be freed directly if it's the root node of
   @parser.
*/

static void json_deallocate (JsonNode *node, JsonParser *parser)
{
    if (node != NULL) {
	if (parser != NULL && node == json_parser_get_root(parser)) {
	    ; /* don't touch it! */
	} else {
	    json_node_free(node);
	}
    }
    if (parser != NULL) {
	g_object_unref(parser);
    }
}

/* common starting point for our json_get* functions */

static JsonNode *get_root_for_path (JsonNode *node,
				    const char *pathstr,
				    int allow_empty,
				    int *err)
{
    JsonPath *path = json_path_new();
    GError *gerr = NULL;
    JsonNode *match = NULL;

    if (json_path_compile(path, pathstr, &gerr)) {
	/* compilation went OK */
	match = json_path_match(path, node);
	if (null_node(match)) {
	    if (match != NULL) {
		json_node_free(match);
		match = NULL;
	    }
	    if (!allow_empty) {
		*err = E_DATA;
	    }
	}
    } else {
	/* compilation failed */
	if (gerr != NULL) {
	    gretl_errmsg_sprintf("jsonget: failed to compile JsonPath: %s",
				 gerr->message);
	    g_error_free(gerr);
	} else {
	    gretl_errmsg_set("jsonget: failed to compile JsonPath");
	}
	*err = E_DATA;
    }

    g_object_unref(path);

    return match;
}

static JsonNode *get_root_for_data (const char *data,
				    const char *path,
				    JsonParser **pjp,
				    int allow_empty,
				    int *err)
{
    GError *gerr = NULL;
    JsonParser *parser;
    JsonNode *root = NULL;
    JsonNode *ret = NULL;

    parser = json_parser_new();
    if (parser == NULL) {
	gretl_errmsg_set("json_get_bundle: couldn't allocate parser");
	*err = E_ALLOC;
	return NULL;
    }

    json_parser_load_from_data(parser, data, -1, &gerr);

    if (gerr != NULL) {
	gretl_errmsg_sprintf("Couldn't parse JSON input: %s",
			     gerr->message);
	g_error_free(gerr);
	g_object_unref(parser);
	*err = E_DATA;
    } else {
	/* check status of @root */
	root = json_parser_get_root(parser);
	if (root == NULL || json_node_is_null(root)) {
	    gretl_errmsg_set("jsonget: got null root node");
	    g_object_unref(parser);
	    *err = E_DATA;
	}
    }

    if (!*err) {
	if (path != NULL) {
	    ret = get_root_for_path(root, path, allow_empty, err);
	    /* clean-up on fail? */
	} else {
	    /* note: this node belong to @parser */
	    ret = root;
	}
    }

    if (pjp != NULL) {
	/* pass back a pointer to enable correct clean-up */
	*pjp = parser;
    }

    return ret;
}

static int output_json_node_value (JsonNode *node,
				   PRN *prn)
{
    GType type = 0;
    int err = 0;

    if (null_node(node)) {
	gretl_errmsg_set("jsonget: got a null node");
	return E_DATA;
    }

    type = json_node_get_value_type(node);

#if 0
    fprintf(stderr, "jsonget: node type %s\n", g_type_name(type));
#endif

    if (!handled_type(type)) {
	const char *s = g_type_name(type);

	gretl_errmsg_sprintf("jsonget: unhandled object type '%s'", s);
	err = E_DATA;
    } else if (type == G_TYPE_STRING) {
	const gchar *s = json_node_get_string(node);

	if (s != NULL) {
	    pputs(prn, s);
	} else {
	    err = E_DATA;
	}
    } else if (type == G_TYPE_DOUBLE) {
	double x = json_node_get_double(node);

	pprintf(prn, "%.15g", x);
    } else if (type == G_TYPE_INT64) {
	gint64 k = json_node_get_int(node);
	double x = (double) k;

	pprintf(prn, "%.15g", x);
    } else {
	int k = json_node_get_boolean(node);
	double x = (double) k;

	pprintf(prn, "%g", x);
    }

    return err;
}

/* for passing deep into callbacks */
struct jsdata {
    int *n_objects;
    int *err;
    PRN *prn;
};

static void show_obj_value (gpointer data, gpointer p)
{
    JsonNode *node = data;
    struct jsdata *jsd = p;

    if (JSON_NODE_HOLDS_ARRAY(node)) {
	fprintf(stderr, " show_obj_value: got array!\n");
    }

    if (node != NULL && !*jsd->err) {
	*jsd->err = output_json_node_value(node, jsd->prn);
	if (!*jsd->err) {
	    *jsd->n_objects += 1;
	    pputc(jsd->prn, '\n');
	}
    }
}

static int excavate_json_object (JsonNode *node,
				 int *n_objects,
				 PRN *prn)
{
    JsonObject *obj = json_node_get_object(node);
    int err = 0;

    if (obj == NULL) {
	return E_DATA;
    }

    if (json_object_get_size(obj) > 0) {
	GList *list = json_object_get_values(obj);
	struct jsdata jsd = {
	    n_objects,
	    &err,
	    prn
	};

	list = json_object_get_values(obj);
	g_list_foreach(list, show_obj_value, &jsd);
	g_list_free(list);
    }

    return err;
}

static int real_json_get (JsonNode *match,
			  int *n_objects,
			  int allow_empty,
			  PRN *prn)
{
    JsonNode *node;
    GType ntype;
    int err = 0;

    *n_objects = 0;

    /* in case we get floating-point output */
    gretl_push_c_numeric_locale();

    if (JSON_NODE_HOLDS_ARRAY(match)) {
	JsonArray *array = json_node_get_array(match);
	int len = 0, index = 0;

	if (non_empty_array(array)) {
	    len = json_array_get_length(array);
	    node = json_array_get_element(array, index);
	} else {
	    node = NULL;
	}

    repeat:

	if (null_node(node)) {
	    gretl_errmsg_set("jsonget: failed to match JsonPath");
	    ntype = 0;
	    err = allow_empty ? 0 : E_DATA;
	    goto bailout;
	} else {
	    ntype = json_node_get_value_type(node);
	}

	if (node != NULL && !handled_type(ntype)) {
	    if (JSON_NODE_HOLDS_ARRAY(node)) {
		/* recurse on array type */
		array = json_node_get_array(node);
		if (non_empty_array(array)) {
		    node = json_array_get_element(array, 0);
		    goto repeat;
		}
	    } else if (json_node_get_node_type(node) == JSON_NODE_OBJECT) {
		err = excavate_json_object(node, n_objects, prn);
		if (!err) {
		    if (index < len - 1) {
			node = json_array_get_element(array, ++index);
			goto repeat;
		    }
		}
	    } else {
		gretl_errmsg_sprintf("jsonget: unhandled array type '%s'",
				     g_type_name(ntype));
		err = E_DATA;
	    }
	} else if (array != NULL) {
	    int i, n = json_array_get_length(array);

	    for (i=0; i<n && !err; i++) {
		node = json_array_get_element(array, i);
		err = output_json_node_value(node, prn);
		if (!err) {
		    *n_objects += 1;
		    if (n > 1) {
			pputc(prn, '\n');
		    }
		}
	    }
	}
    } else {
	/* not an array-holding node */
	err = output_json_node_value(match, prn);
	if (!err) {
	    *n_objects += 1;
	}
    }

 bailout:

    gretl_pop_c_numeric_locale();

    return err;
}

/*
  @data: JSON buffer.
  @path: the JsonPath to the target info.
  @n_objects: location to receive the number of pieces
  of information retrieved, or NULL.
  @err: location to receive error code.

  On success, returns an allocated string. If the "target"
  is an array, the members are printed one per line. This
  function handles target types of double, int or string;
  in the case of doubles or ints, their string representation
  is returned (using the C locale for doubles).
*/

char *json_get_string (const char *data, const char *path,
		       int *n_objects, int *err)
{
    JsonNode *root;
    JsonParser *parser = NULL;
    int allow_empty;
    char *ret = NULL;
    int n = 0;

    if (data == NULL || path == NULL) {
	if (n_objects != NULL) {
	    *n_objects = 0;
	}
	*err = E_DATA;
	return NULL;
    }

    allow_empty = n_objects != NULL;
    root = get_root_for_data(data, path, &parser, allow_empty, err);

    if (!*err) {
	PRN *prn = gretl_print_new(GRETL_PRINT_BUFFER, err);

	if (!*err) {
	    *err = real_json_get(root, &n, allow_empty, prn);
	    if (!*err) {
		if (n == 0 && allow_empty) {
		    ret = gretl_strdup("");
		} else {
		    ret = gretl_print_steal_buffer(prn);
		}
	    }
	    gretl_print_destroy(prn);
	}
    }

    json_deallocate(root, parser);

    if (*err) {
	fprintf(stderr, "json_get: err = %d\n", *err);
    }

    if (n_objects != NULL) {
	*n_objects = n;
    }

    return ret;
}

/* start code subserving json_get_bundle() */

#define JB_DEBUG 0

struct jbundle_ {
    gretl_bundle *b0;
    gretl_bundle *bcurr;
    gchar ***a;
    int nlev;
    int level;
};

typedef struct jbundle_ jbundle;

static void free_pathbits (gchar ***a, int nlev)
{
    int i;

    for (i=0; i<nlev; i++) {
	if (a[i] != NULL) {
	    g_strfreev(a[i]);
	}
    }
    g_free(a);
}

#if JB_DEBUG

static void print_pathbits (gchar ***a, int nlev)
{
    int i, j, n;

    fprintf(stderr, "parsed pathbits:\n");
    for (i=0; i<nlev; i++) {
	n = g_strv_length(a[i]);
	fprintf(stderr, "i=%d: ", i);
	for (j=0; j<n; j++) {
	    fprintf(stderr, "'%s'", a[i][j]);
	}
	fputc('\n', stderr);
    }
    fputc('\n', stderr);
}

#endif

/* Given a user-suppled @path string, parse it into chunks
   specific to levels in the JSON tree; each of these should
   be a single name or {name1,name2,...} or "*".
*/

static int jb_make_pathbits (jbundle *jb, const char *path)
{
    gchar **S = g_strsplit(path, "/", -1);
    gchar ***a = NULL;
    int i, nlev, n;
    int err = 0;

    nlev = g_strv_length(S);
    if (nlev == 0) {
	return 0;
    }

    a = g_malloc0(nlev * sizeof *a);

    for (i=0; i<nlev && !err; i++) {
	g_strstrip(S[i]);
	if (*S[i] == '{') {
	    n = strlen(S[i]);
	    if (S[i][n-1] == '}') {
		S[i][0] = S[i][n-1] = ' ';
		g_strstrip(S[i]);
		a[i] = g_strsplit(S[i], ",", -1);
	    } else {
		err = E_PARSE;
	    }
	} else {
	    a[i] = g_malloc(2 * sizeof **a);
	    a[i][0] = g_strdup(S[i]);
	    a[i][1] = NULL;
	}
    }

    g_strfreev(S);

    if (err) {
	free_pathbits(a, nlev);
    } else {
#if JB_DEBUG
	print_pathbits(a, nlev);
#endif
	jb->a = a;
	jb->nlev = nlev;
    }

    return err;
}

static int jb_do_object (JsonReader *reader, jbundle *jb,
			 gretl_array *a);
static int jb_do_array (JsonReader *reader, jbundle *jb,
			gretl_array *a);
static int jb_do_value (JsonReader *reader, jbundle *jb,
			gretl_array *a, int i);

/* Check whether a given JSON element is wanted in the
   context of json_get_bundle(). We treat every element
   as wanted by default, but if we have a path spec for
   the current structural level (jb->level) and the
   current node has a name, we need to check for a match.
*/

static int is_wanted (jbundle *jb, JsonReader *reader)
{
    int i = jb->level - 1;
    int ret = 1;

    if (jb->a != NULL && i < jb->nlev) {
	/* there's a relevant path spec */
	const gchar *name = json_reader_get_member_name(reader);

	if (name != NULL) {
	    int j, n = g_strv_length(jb->a[i]);

#if JB_DEBUG
	    fprintf(stderr, "test for inclusion of %s at level %d: ",
		    name, jb->level);
#endif
	    ret = 0;
	    if (strlen(jb->a[i][0]) == 0 || !strcmp(jb->a[i][0], "*")) {
		/* everything "matches" */
		ret = 1;
	    } else {
		/* look for a specific match */
		for (j=0; j<n; j++) {
		    if (!strcmp(name, jb->a[i][j])) {
			ret = 1;
			break;
		    }
		}
	    }
#if JB_DEBUG
	    fprintf(stderr, "%s\n", ret ? "yes" : "no");
#endif
	}
    }

    return ret;
}

/* Add a new bundle to the tree -- either as a named
   member of the bundle jb->bcurr, or, if @a is non-NULL,
   as an anonymous element in an array of bundles.
*/

static int jb_add_bundle (jbundle *jb, const char *name,
			  gretl_array *a, int i)
{
    gretl_bundle *b = gretl_bundle_new();
    int err;

    if (b == NULL) {
	err = E_ALLOC;
    } else if (a != NULL) {
	err = gretl_array_set_bundle(a, i, b, 0);
    } else {
	if (name == NULL || *name == '\0') {
	    gretl_errmsg_set("JSON object member name is missing");
	    err = E_DATA;
	} else {
	    err = gretl_bundle_donate_data(jb->bcurr, name, b,
					   GRETL_TYPE_BUNDLE, 0);
	}
    }

    if (err) {
	gretl_bundle_destroy(b);
    } else {
	jb->bcurr = b;
    }

    return err;
}

/* Process a JSON object node: it becomes a gretl bundle */

static int jb_do_object (JsonReader *reader, jbundle *jb,
			 gretl_array *a)
{
    gchar **S = NULL;
    int i, n, err = 0;

    n = json_reader_count_members(reader);
    S = json_reader_list_members(reader);

#if JB_DEBUG
    if (a == NULL) {
	const gchar *name = json_reader_get_member_name(reader);
	fprintf(stderr, "level %d: got object, name '%s', %d member(s)\n",
		jb->level, name == NULL ? "NULL" : name, n);
    } else {
	fprintf(stderr, "level %d: got object (array element), %d member(s)\n",
		jb->level, n);
    }
# if 0
    for (i=0; i<n; i++) {
	fprintf(stderr, "  %s\n", S[i]);
    }
# endif
#endif

    for (i=0; i<n && !err; i++) {
	json_reader_read_member(reader, S[i]);
	if (json_reader_is_object(reader)) {
	    int lsave = jb->level;

	    jb->level += 1;
	    if (is_wanted(jb, reader)) {
		gretl_bundle *bsave = jb->bcurr;

		err = jb_add_bundle(jb, S[i], NULL, 0);
		if (!err) {
		    err = jb_do_object(reader, jb, NULL);
		}
		jb->bcurr = bsave;
	    }
	    jb->level = lsave;
	} else if (json_reader_is_array(reader)) {
	    int lsave = jb->level;

	    jb->level += 1;
	    if (is_wanted(jb, reader)) {
		err = jb_do_array(reader, jb, NULL);
	    }
	    jb->level = lsave;
	} else if (json_reader_is_value(reader)) {
	    int lsave = jb->level;

	    jb->level += 1;
	    if (is_wanted(jb, reader)) {
		err = jb_do_value(reader, jb, NULL, 0);
	    }
	    jb->level = lsave;
	}
	json_reader_end_member(reader);
    }

    g_strfreev(S);

    return err;
}

/* Switch the type of the array @a from its original value
   to @targ, but only if we haven't already added elements
   that make the type of @a immutable.
*/

static int jb_transmute_array (gretl_array *a,
			       GretlType targ,
			       GretlType *pt)
{
    int err;

    err = gretl_array_set_type(a, targ);
    if (err) {
	gretl_errmsg_set("JSON array: can't mix types");
	fprintf(stderr, "jb_transmute_array: array type was %s, trying to change to %s\n",
		gretl_type_get_name(*pt), gretl_type_get_name(targ));
    } else {
	*pt = targ;
    }

    return err;
}

/* Process a JSON array node: we'll construct an array of
   strings, bundles, or arrays.
*/

static int jb_do_array (JsonReader *reader, jbundle *jb,
			gretl_array *a0)
{
    GretlType atype;
    const gchar *name;
    gretl_array *a;
    gboolean ok;
    int i, n, err = 0;

    n = json_reader_count_elements(reader);
    name = json_reader_get_member_name(reader);

#if JB_DEBUG
    fprintf(stderr, "level %d: got array, name '%s', %d element(s)\n",
	    jb->level, name == NULL ? "NULL" : name, n);
#endif

    /* Packing an array into a bundle requires a key, so
       it's a problem if an array node doesn't have a name.
    */
    if (name == NULL || name[0] == '\0') {
	name = "anon";
    }

    /* assume an array of strings by default */
    atype = GRETL_TYPE_STRINGS;
    a = gretl_array_new(atype, n, &err);

    for (i=0; i<n && !err; i++) {
	ok = json_reader_read_element(reader, i);
	if (!ok) {
	    gretl_errmsg_set("JSON array: couldn't read element");
	    err = E_DATA;
	    break;
	}
	if (json_reader_is_value(reader)) {
	    if (atype != GRETL_TYPE_STRINGS) {
		gretl_errmsg_set("JSON array: can't mix types");
		fprintf(stderr, "JSON array type %s, got 'value' member\n",
			gretl_type_get_name(atype));
		err = E_DATA;
	    } else {
		err = jb_do_value(reader, jb, a, i);
	    }
	} else if (json_reader_is_object(reader)) {
	    if (atype != GRETL_TYPE_BUNDLES) {
		/* try switching to bundles */
		err = jb_transmute_array(a, GRETL_TYPE_BUNDLES, &atype);
	    }
	    if (!err) {
		gretl_bundle *bsave = jb->bcurr;

		err = jb_add_bundle(jb, NULL, a, i);
		if (!err) {
		    int lsave = jb->level;

		    jb->level += 1;
		    err = jb_do_object(reader, jb, a);
		    jb->level = lsave;
		}
		jb->bcurr = bsave;
	    }
	} else if (json_reader_is_array(reader)) {
	    if (atype != GRETL_TYPE_ARRAYS) {
		/* try switching to arrays */
		err = jb_transmute_array(a, GRETL_TYPE_ARRAYS, &atype);
		if (err && atype == GRETL_TYPE_STRINGS) {
		    /* 2020-04-21: we got an array @a whose first element
		       was a string, followed by an array. This seems
		       malformed and we'll try just skipping @a.
		    */
		    fprintf(stderr, "skipping malformed array\n");
		    gretl_array_destroy(a);
		    a = NULL;
		    err = 0;
		    gretl_error_clear();
		    goto end_element;
		}
	    }
	    if (!err) {
		int lsave = jb->level;

		jb->level += 1;
		err = jb_do_array(reader, jb, a);
		jb->level = lsave;
	    }
	} else {
	    gretl_errmsg_set("JSON array: unrecognized type");
	    err = E_DATA;
	}

    end_element:
	json_reader_end_element(reader);
    }

    if (err) {
	gretl_array_destroy(a);
    } else if (a0 != NULL && a != NULL) {
	/* nesting arrays */
	int idx = gretl_array_get_next_index(a0);

	if (idx < 0) {
	    gretl_array_destroy(a);
	    err = E_DATA;
	} else {
#if JB_DEBUG
	    fprintf(stderr, "nesting array, level %d, name %s\n",
		    jb->level, name);
#endif
	    err = gretl_array_set_array(a0, idx, a, 0);
	}
    } else if (a != NULL) {
	/* sticking array into bundle */
	err = gretl_bundle_donate_data(jb->bcurr, name, a,
				       GRETL_TYPE_ARRAY, 0);
    }

    return err;
}

/* Process a JSON value node: we convert all array values
   to strings */

static int jb_do_value (JsonReader *reader, jbundle *jb,
			gretl_array *a, int i)
{
    JsonNode *node = json_reader_get_value(reader);
    const gchar *name, *typename;
    char tmp[32];
    GType type;
    int err = 0;

    name = json_reader_get_member_name(reader);
    type = json_node_get_value_type(node);
    typename = g_type_name(type);

#if JB_DEBUG
    fprintf(stderr, "  got value: name='%s', type %s (%ld)\n",
	    name, typename, type);
#endif

    if (a == NULL && type == 0) {
	/* bundle member: skip null nodes */
	return 0;
    }

    if (a == NULL && (name == NULL || name[0] == '\0')) {
	name = "anon";
    }

    if (type == G_TYPE_INT64) {
	int k = (int) json_reader_get_int_value(reader);

	if (a != NULL) {
	    sprintf(tmp, "%d", k);
	    gretl_array_set_string(a, i, tmp, 1);
	} else {
	    gretl_bundle_set_int(jb->bcurr, name, k);
	}
    } else if (type == G_TYPE_DOUBLE) {
	gdouble x = json_reader_get_double_value(reader);

	if (a != NULL) {
	    sprintf(tmp, "%.15g", x);
	    gretl_array_set_string(a, i, tmp, 1);
	} else {
	    gretl_bundle_set_scalar(jb->bcurr, name, x);
	}
    } else if (type == G_TYPE_STRING) {
	const gchar *s = json_reader_get_string_value(reader);

	if (a != NULL) {
	    gretl_array_set_string(a, i, (char *) s, 1);
	} else {
	    gretl_bundle_set_string(jb->bcurr, name, s);
	}
    } else if (type == G_TYPE_BOOLEAN) {
	int k = (int) json_reader_get_boolean_value(reader);

	if (a != NULL) {
	    sprintf(tmp, "%d", k);
	    gretl_array_set_string(a, i, tmp, 1);
	} else {
	    gretl_bundle_set_int(jb->bcurr, name, k);
	}
    } else if (type == 0) {
	/* in array context: null object -> empty string */
	gretl_array_set_string(a, i, "", 1);
    } else {
	gretl_errmsg_sprintf("Unhandled JSON value of type %s\n",
			     typename);
	err = E_DATA;
    }

    return err;
}

/* end code subserving json_get_bundle() */

/*
  @data: JSON buffer.
  @path: array of strings identifying JSON objects to include,
  or NULL to retrieve all.
  @err: location to receive error code.

  On success, returns an allocated gretl_bundle whose
  structure mirrors that of the root JSON object.
*/

gretl_bundle *json_get_bundle (const char *data,
			       const char *path,
			       int *err)
{
    gretl_bundle *ret = NULL;
    jbundle jb = {0};
    JsonNode *root;
    JsonParser *parser = NULL;
    JsonReader *reader;

    if (data == NULL) {
	gretl_errmsg_set("json_get_bundle: no data supplied");
	*err = E_DATA;
	return NULL;
    }

    root = get_root_for_data(data, NULL, &parser, 1, err);
    if (*err) {
	return NULL;
    }

    if (path != NULL) {
	if (*path == '/') {
	    path++;
	}
	*err = jb_make_pathbits(&jb, path);
	if (*err) {
	    return NULL;
	}
    }

    jb.bcurr = jb.b0 = gretl_bundle_new();

    reader = json_reader_new(root);
    gretl_push_c_numeric_locale();

    if (json_reader_is_object(reader)) {
	*err = jb_do_object(reader, &jb, NULL);
    } else if (json_reader_is_array(reader)) {
	*err = jb_do_array(reader, &jb, NULL);
    } else if (json_reader_is_value(reader)) {
	*err = jb_do_value(reader, &jb, NULL, 0);
    }

    gretl_pop_c_numeric_locale();
    g_object_unref(reader);

    json_deallocate(root, parser);
    if (jb.a != NULL) {
	free_pathbits(jb.a, jb.nlev);
    }

    if (*err) {
	/* trash bundle on failure */
	gretl_bundle_destroy(jb.b0);
    } else {
	ret = jb.b0;
    }

    return ret;
}

static char *path_first (const char *path)
{
    const char *p = strchr(path, '/');
    char *ret;

    if (p == NULL) {
	ret = g_strdup(path);
    } else {
	ret = g_strndup(path, p - path);
    }

    return ret;
}

/* public function, referenced in geneval.c (I think) */

gretl_array *json_get_array (const char *data,
			     const char *path,
			     int *err)
{
    gretl_array *ret = NULL;
    gretl_bundle *b;

    b = json_get_bundle(data, path, err);

    if (!*err) {
	/* FIXME handling of path -> key */
	gchar *key = path_first(path);
	GretlType type = 0;
	void *ptr;

	ptr = gretl_bundle_steal_data(b, key, &type, NULL, err);
	if (ptr == NULL) {
	    gretl_errmsg_sprintf("json_get_array: got NULL for key '%s'",
				 key);
	    *err = E_DATA;
	} else if (type != GRETL_TYPE_ARRAY) {
	    gretl_errmsg_sprintf("json_get_array: got wrong type %s for key '%s'",
				 gretl_type_get_name(type), key);
	    *err = E_DATA;
	} else {
	    ret = ptr;
	}
	gretl_bundle_destroy(b);
	g_free(key);
    }

    return ret;
}

static int filter_bundle_tree (gretl_bundle *b, gretl_array *A)
{
    gretl_array *K, *ai;
    gretl_bundle *bj;
    void *child;
    GretlType type;
    char **keys;
    int addit = 1;
    int i, j, nb, n = 0;
    int err = 0;

    K = gretl_bundle_get_keys(b, NULL);
    keys = gretl_array_get_strings(K, &n);

    for (i=0; i<n; i++) {
	if (!strcmp(keys[i], "children") ||
	    !strcmp(keys[i], "category_tree")) {
	    /* bundle has children: not terminal */
	    addit = 0;
	    break;
	}
    }

    if (addit) {
	/* push copy of bundle onto array */
	err = gretl_array_append_bundle(A, b, 1);
    }

    for (i=0; i<n && !err; i++) {
	child = gretl_bundle_get_data(b, keys[i], &type, NULL, NULL);
	if (type == GRETL_TYPE_BUNDLE) {
	    filter_bundle_tree((gretl_bundle *) child, A);
	} else if (type == GRETL_TYPE_ARRAY) {
	    ai = (gretl_array *) child;
	    type = gretl_array_get_content_type(ai);
	    if (type == GRETL_TYPE_BUNDLE) {
		nb = gretl_array_get_length(ai);
		for (j=0; j<nb; j++) {
		    bj = gretl_array_get_bundle(ai, j);
		    filter_bundle_tree(bj, A);
		}
	    }
	}
    }

    gretl_array_destroy(K);

    return err;
}

/* Given a bundle @b produced by json_get_bundle(), this function
   constructs an array of bundles holding all and only the
   "terminal" bundles within @b. A "terminal" bundle is one
   which contains no member named "children" (or "category tree").

   We want this for handling the results from the category_tree
   of a dbnomics provider: if we're looking to list datasets we
   don't want to include any intermediate nodes that are not
   themselves datasets but rather groupings of datasets. The
   latter can be distinguished by the presence of a "children"
   node; an actual dataset node never has "children".
*/

gretl_array *json_bundle_get_terminals (gretl_bundle *b, int *err)
{
    gretl_array *a;

    a = gretl_array_new(GRETL_TYPE_BUNDLES, 0, err);

    if (!*err) {
	*err = filter_bundle_tree(b, a);
    }

    if (*err) {
	gretl_array_destroy(a);
	a = NULL;
    }

    return a;
}
