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
	gretl_errmsg_sprintf("jsonget: unhandled object type '%s'",
			     g_type_name(type));
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
	    json_node_unref(match);
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

char *json_get (const char *data, const char *path, int *n_objects,
		int *err)
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
