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

static int non_empty_array (JsonArray *array)
{
    return array != NULL && json_array_get_length(array) > 0;
}

static int null_node (JsonNode *node)
{
    return node == NULL || json_node_is_null(node);
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

static int real_json_get (JsonParser *parser, const char *pathstr,
			  int *n_objects, int allow_fail,
			  PRN *prn)
{
    GError *gerr = NULL;
    JsonNode *match, *node;
    JsonPath *path;
    GType ntype;
    int err = 0;

    *n_objects = 0;

    node = json_parser_get_root(parser);

    if (node == NULL || json_node_is_null(node)) {
	gretl_errmsg_set("jsonget: got null root node");
	return E_DATA;
    }

    path = json_path_new();

    if (!json_path_compile(path, pathstr, &gerr)) {
	if (gerr != NULL) {
	    gretl_errmsg_sprintf("jsonget: failed to compile JsonPath: %s",
				 gerr->message);
	    g_error_free(gerr);
	} else {
	    gretl_errmsg_set("jsonget: failed to compile JsonPath");
	}
	g_object_unref(path);
	return E_DATA;
    }

    match = json_path_match(path, node);

    if (null_node(match)) {
	if (match != NULL) {
	    json_node_free(match);
	}
	g_object_unref(path);
	return allow_fail ? 0 : E_DATA;
    }

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
	    err = allow_fail ? 0 : E_DATA;
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

    json_node_free(match);
    g_object_unref(path);

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
    GError *gerr = NULL;
    JsonParser *parser;
    char *ret = NULL;
    int n = 0;

    if (data == NULL || path == NULL) {
	if (n_objects != NULL) {
	    *n_objects = 0;
	}
	*err = E_DATA;
	return NULL;
    }

    parser = json_parser_new();
    if (parser == NULL) {
	gretl_errmsg_set("json_parser_new returned NULL!\n");
	*err = 1;
	return NULL;
    }

    json_parser_load_from_data(parser, data, -1, &gerr);

    if (gerr != NULL) {
	gretl_errmsg_sprintf("Couldn't parse JSON input: %s",
			     gerr->message);
	g_error_free(gerr);
	*err = E_DATA;
    } else {
	PRN *prn = gretl_print_new(GRETL_PRINT_BUFFER, err);

	if (!*err) {
	    int allow_fail = n_objects != NULL;

	    *err = real_json_get(parser, path, &n, allow_fail, prn);
	    if (!*err) {
		if (n == 0 && allow_fail) {
		    ret = gretl_strdup("");
		} else {
		    ret = gretl_print_steal_buffer(prn);
		}
	    }
	    gretl_print_destroy(prn);
	}
    }

    if (*err) {
	fprintf(stderr, "json_get: err = %d\n", *err);
    }

    if (n_objects != NULL) {
	*n_objects = n;
    }

    g_object_unref(parser);

    return ret;
}

/* start code subserving json_get_bundle() */

#define JB_DEBUG 0

struct jbundle_ {
    gretl_bundle *b0;
    gretl_bundle *curr;
    char **targets;
    int n_targets;
};

/* The idea with "targets" is that we could ignore some stuff
   if it's not wanted, and just add the target elements.
   But this is not implemented at all right now.
*/

typedef struct jbundle_ jbundle;

static int jb_do_object (JsonReader *reader, jbundle *jb);
static int jb_do_array (JsonReader *reader, jbundle *jb);
static int jb_do_value (JsonReader *reader, jbundle *jb,
			gretl_array *a, int i);

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

static int jb_do_object (JsonReader *reader, jbundle *jb)
{
    gretl_bundle *btop = jb->curr;
    gchar **S = NULL;
    int i, n, err = 0;

    n = json_reader_count_members(reader);

#if JB_DEBUG
    const gchar *name = json_reader_get_member_name(reader);
    fprintf(stderr, "got object, %d members, name '%s'\n", n,
	    name == NULL ? "NULL" : name);
#endif

    S = json_reader_list_members(reader);

    for (i=0; i<n && !err; i++) {
	json_reader_read_member(reader, S[i]);
	if (json_reader_is_object(reader)) {
	    err = jb_add_bundle(jb, S[i], NULL, 0);
	    if (!err) {
		err = jb_do_object(reader, jb);
	    }
	    jb->curr = btop;
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

static int jb_transmute_array (gretl_array **pa, int n)
{
    int err;

    gretl_array_destroy(*pa);
    *pa = gretl_array_new(GRETL_TYPE_BUNDLES, n, &err);
    return err;
}

/* process a JSON array node: we'll construct either an
   array of strings or an array of bundles */

static int jb_do_array (JsonReader *reader, jbundle *jb)
{
    gretl_bundle *btop = jb->curr;
    const gchar *name;
    gretl_array *a;
    gboolean ok;
    int ns = 0, nb = 0;
    int i, n, err = 0;

    n = json_reader_count_elements(reader);
    name = json_reader_get_member_name(reader);

#if JB_DEBUG
    fprintf(stderr, "got array, %d elements, name %s\n", n,
	   name == NULL ? "NULL" : name);
#endif

    /* we'll assume strings by default */
    a = gretl_array_new(GRETL_TYPE_STRINGS, n, &err);

    for (i=0; i<n && !err; i++) {
	ok = json_reader_read_element(reader, i);
	if (ok && json_reader_is_value(reader)) {
	    if (nb > 0) {
		/* we already switched to bundles! */
		fprintf(stderr, "element %d: can't mix types in array!\n", i);
		err = E_DATA;
	    } else {
		err = jb_do_value(reader, jb, a, i);
		ns++;
	    }
	} else if (json_reader_is_object(reader)) {
	    if (ns > 0) {
		fprintf(stderr, "element %d: can't mix types in array!\n", i);
		err = E_DATA;
	    } else {
		/* switch to bundles if we haven't already got
		   any string elements
		*/
		err = jb_transmute_array(&a, n);
		jb_add_bundle(jb, NULL, a, i);
		err = jb_do_object(reader, jb);
		jb->curr = btop;
		nb++;
	    }
	} else if (json_reader_is_array(reader)) {
	    /* the gretl_array type cannot be nested */
	    fprintf(stderr, "element %d: nesting of arrays is not handled!\n", i);
	    err = E_DATA;
	} else {
	    /* ?? */
	    fprintf(stderr, "array: element %d not handled!\n", i);
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

/* process a JSON value node: we convert all values to strings */

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
    fprintf(stderr, "  got value (%s), type %s\n", name, typename);
#endif

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
  @targets: array of strings identifying members to retrieve,
  or NULL to retrieve all.
  @err: location to receive error code.

  On success, returns an allocated gretl_bundle whose
  structure mirrors that of the root JSON object.
*/

gretl_bundle *json_get_bundle (const char *data,
			       gretl_array *targets,
			       int *err)
{
    jbundle jb = {0};
    GError *gerr = NULL;
    JsonParser *parser;
    JsonReader *reader;

    if (data == NULL) {
	gretl_errmsg_set("json_get_bundle: no data supplied");
	*err = E_DATA;
	return NULL;
    }

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
	*err = E_DATA;
	return NULL;
    }

    if (targets != NULL) {
	/* doesn't do anything yet! */
	jb.targets = gretl_array_get_strings(targets, &jb.n_targets);
	fprintf(stderr, "json_get_bundle: found %d targets\n",
		jb.n_targets);
    }

    jb.b0 = gretl_bundle_new();
    jb.curr = jb.b0;

    gretl_push_c_numeric_locale();

    reader = json_reader_new(json_parser_get_root(parser));
    if (json_reader_is_object(reader)) {
	*err = jb_do_object(reader, &jb);
    }

    gretl_pop_c_numeric_locale();

    g_object_unref(reader);
    g_object_unref(parser);

    return jb.b0;
}
