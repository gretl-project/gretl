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

#include <glib-object.h>
#include <json-glib/json-glib.h>

#define handled_type(t) (t == G_TYPE_STRING || \
			 t == G_TYPE_DOUBLE || \
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
    } else {
	gint64 k = json_node_get_int(node);
	double x = (double) k;

	pprintf(prn, "%.15g", x);
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
    gretl_bundle *curr;
    char **excludes;
    int n_exclude;
};

typedef struct jbundle_ jbundle;

static int jb_do_object (JsonReader *reader, jbundle *jb);
static int jb_do_array (JsonReader *reader, jbundle *jb);
static int jb_do_value (JsonReader *reader, jbundle *jb,
			gretl_array *a, int i);

static int is_excluded (jbundle *jb, JsonReader *reader)
{
    if (jb->n_exclude > 0) {
	const gchar *name = json_reader_get_member_name(reader);
	int i;

	if (name == NULL) {
	    return 0;
	}
	for (i=0; i<jb->n_exclude; i++) {
	    if (!strcmp(name, jb->excludes[i])) {
		return 1;
	    }
	}
    }

    return 0;
}

/* Add a new bundle to the tree -- either as a named
   member of the bundle jb->curr, or, if @a is non-NULL,
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
	err = gretl_bundle_donate_data(jb->curr, name, b,
				       GRETL_TYPE_BUNDLE, 0);
    }

    if (!err) {
	jb->curr = b;
    }

    return err;
}

/* Process a JSON object node: it becomes a gretl bundle */

static int jb_do_object (JsonReader *reader, jbundle *jb)
{
    const gchar *name;
    gchar **S = NULL;
    int i, n, err = 0;

    n = json_reader_count_members(reader);
    S = json_reader_list_members(reader);

#if JB_DEBUG
    name = json_reader_get_member_name(reader);
    fprintf(stderr, "got object, %d members, name '%s'\n", n,
	    name == NULL ? "NULL" : name);
    for (i=0; i<n; i++) {
	fprintf(stderr, "  %s\n", S[i]);
    }
#endif

    for (i=0; i<n && !err; i++) {
	json_reader_read_member(reader, S[i]);
	if (json_reader_is_object(reader)) {
	    if (!is_excluded (jb, reader)) {
		gretl_bundle *bsave = jb->curr;

		err = jb_add_bundle(jb, S[i], NULL, 0);
		if (!err) {
		    err = jb_do_object(reader, jb);
		}
		jb->curr = bsave;
	    }
	} else if (json_reader_is_array(reader)) {
	    err = jb_do_array(reader, jb);
	} else if (json_reader_is_value(reader)) {
	    err = jb_do_value(reader, jb, NULL, 0);
	}
	json_reader_end_member(reader);
    }

    g_strfreev(S);

    return err;
}

static int jb_transmute_array (gretl_array **pa,
			       GretlType *pt,
			       int n, int ns)
{
    int err = 0;

    if (ns > 0) {
	gretl_errmsg_set("JSON array: can't mix types");
	err = E_DATA;
    } else {
	gretl_array_destroy(*pa);
	*pt = GRETL_TYPE_BUNDLES;
	*pa = gretl_array_new(*pt, n, &err);
    }

    return err;
}

/* Process a JSON array node: we'll construct either an
   array of strings or an array of bundles. Since gretl
   arrays cannot be nested we're somewhat more restrictive
   here than in the JSON spec.
*/

static int jb_do_array (JsonReader *reader, jbundle *jb)
{
    GretlType atype;
    const gchar *name;
    gretl_array *a;
    gboolean ok;
    int ns = 0;
    int i, n, err = 0;

    n = json_reader_count_elements(reader);
    name = json_reader_get_member_name(reader);

#if JB_DEBUG
    fprintf(stderr, "got array, %d element(s), name %s\n", n,
	    name == NULL ? "NULL" : name);
#endif

    /* Arrays can be packed only into bundles, and that
       requires a key, so it's a problem if an array
       node doesn't have a name. Hopefully this should
       occur only if the "root" node is itself an array.
    */
    if (name == NULL || name[0] == '\0') {
	name = "anon";
    }

    /* we'll assume an array of strings by default */
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
	    if (atype == GRETL_TYPE_BUNDLES) {
		/* we already switched to bundles! */
		gretl_errmsg_set("JSON array: can't mix types");
		err = E_DATA;
	    } else {
		err = jb_do_value(reader, jb, a, i);
		ns++;
	    }
	} else if (json_reader_is_object(reader)) {
	    if (atype != GRETL_TYPE_BUNDLES) {
		/* try switching to bundles */
		err = jb_transmute_array(&a, &atype, n, ns);
	    }
	    if (!err) {
		gretl_bundle *bsave = jb->curr;

		err = jb_add_bundle(jb, NULL, a, i);
		if (!err) {
		    err = jb_do_object(reader, jb);
		}
		jb->curr = bsave;
	    }
	} else if (json_reader_is_array(reader)) {
	    /* the gretl_array type cannot be nested */
	    gretl_errmsg_set("JSON array: arrays cannot be nested");
	    err = E_DATA;
	} else {
	    gretl_errmsg_set("JSON array: unrecognized type");
	    err = E_DATA;
	}
	json_reader_end_element(reader);
    }

    if (!err) {
	err = gretl_bundle_donate_data(jb->curr, name, a,
				       GRETL_TYPE_ARRAY, 0);
    } else if (a != NULL) {
	gretl_array_destroy(a);
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
    fprintf(stderr, "  got value: name='%s', type %s\n", name, typename);
#endif

    if (a == NULL && (name == NULL || name[0] == '\0')) {
	name = "anon";
    }

    if (type == G_TYPE_INT64) {
	int k = (int) json_reader_get_int_value(reader);

	if (a != NULL) {
	    sprintf(tmp, "%d", k);
	    gretl_array_set_string(a, i, tmp, 1);
	} else {
	    gretl_bundle_set_int(jb->curr, name, k);
	}
    } else if (type == G_TYPE_DOUBLE) {
	gdouble x = json_reader_get_double_value(reader);

	if (a != NULL) {
	    sprintf(tmp, "%.15g", x);
	    gretl_array_set_string(a, i, tmp, 1);
	} else {
	    gretl_bundle_set_scalar(jb->curr, name, x);
	}
    } else if (type == G_TYPE_STRING) {
	const gchar *s = json_reader_get_string_value(reader);

	if (a != NULL) {
	    gretl_array_set_string(a, i, (char *) s, 1);
	} else {
	    gretl_bundle_set_string(jb->curr, name, s);
	}
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
  @excludes: array of strings identifying JSON objects to exclude,
  or NULL to retrieve all.
  @err: location to receive error code.

  On success, returns an allocated gretl_bundle whose
  structure mirrors that of the root JSON object.
*/

gretl_bundle *json_get_bundle (const char *data,
			       gretl_array *excludes,
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

    if (excludes != NULL) {
	jb.excludes = gretl_array_get_strings(excludes, &jb.n_exclude);
    }

    jb.b0 = gretl_bundle_new();
    jb.curr = jb.b0;

    reader = json_reader_new(root);
    gretl_push_c_numeric_locale();

    if (json_reader_is_object(reader)) {
	*err = jb_do_object(reader, &jb);
    } else if (json_reader_is_array(reader)) {
	*err = jb_do_array(reader, &jb);
    } else if (json_reader_is_value(reader)) {
	*err = jb_do_value(reader, &jb, NULL, 0);
    }

    gretl_pop_c_numeric_locale();
    g_object_unref(reader);

    json_deallocate(root, parser);

    if (*err) {
	/* trash bundle on failure */
	gretl_bundle_destroy(jb.b0);
    } else {
	ret = jb.b0;
    }

    return ret;
}
