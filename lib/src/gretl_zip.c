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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <glib.h>

#include "config.h"
#include "gretl_zip.h"

#ifdef USE_GSF

#include <gsf/gsf-utils.h>

#include <gsf/gsf-input-stdio.h>
#include <gsf/gsf-infile.h>
#include <gsf/gsf-infile-stdio.h>
#include <gsf/gsf-infile-zip.h>

#include <gsf/gsf-output-stdio.h>
#include <gsf/gsf-outfile-stdio.h>
#include <gsf/gsf-outfile.h>
#include <gsf/gsf-outfile-zip.h>

#define ZDEBUG 1

#define CHUNK 8192

#define gsf_is_dir(i) (GSF_IS_INFILE(i) && gsf_infile_num_children(GSF_INFILE(i)) >= 0)

static void ensure_gsf_init (void)
{
    static int initted;

    if (!initted) {
	gsf_init();
	initted = 1;
    }
}

static int transcribe_gsf_data (GsfInput *input, GsfOutput *output)
{
    guint8 const *data;
    size_t len;

    while ((len = gsf_input_remaining(input)) > 0) {
	if (len > CHUNK) {
	    len = CHUNK;
	}
	if ((data = gsf_input_read(input, len, NULL)) == NULL) {
	    return 1;
	}
	if (!gsf_output_write(output, len, data)) {
	    return 1;
	}
    }

    return 0;
}

static int clone_to_zip (GsfInput *input, GsfOutput *output)
{
    int err = 0;

    if (gsf_input_size(input) > 0) {
	err = transcribe_gsf_data(input, output);
    } 

    if (GSF_IS_INFILE(input) && gsf_infile_num_children(GSF_INFILE(input)) > 0) {
	GsfInfile *in = GSF_INFILE(input);
	GsfOutfile *out = GSF_OUTFILE(output);
	int i, n = gsf_infile_num_children(in);
	char const *name;
	GsfInput *src;
	GsfOutput *dest;
	gboolean is_dir;
	GDateTime *modtime;

	for (i=0; i<n && !err; i++) {
	    src = gsf_infile_child_by_index(in, i);
	    name = gsf_infile_name_by_index(in, i);
#if ZDEBUG
	    fprintf(stderr, "clone_to_zip: i = %d (%s)\n", i, name); 
#endif
	    is_dir = gsf_is_dir(src);
	    modtime = gsf_input_get_modtime(src);
	    dest = gsf_outfile_new_child_full(out,
					      gsf_infile_name_by_index(in, i),
					      is_dir,
					      "modtime", modtime,
					      NULL);
	    err = clone_to_zip(src, dest);
	}
    }

    gsf_output_close(output);
    g_object_unref(G_OBJECT(output));
    g_object_unref(G_OBJECT(input));

    return err;
}

int gretl_make_zipfile (const char *fname, const char *path,
			GError **gerr)
{
    GsfInfile *infile;
    GsfOutput *output, *ziproot;
    GsfOutfile *outfile = NULL;
    int err = 0;

    ensure_gsf_init();

#if ZDEBUG
    fprintf(stderr, "gretl_make_zipfile (gsf); fname='%s', path='%s'\n", fname, path);
#endif

    infile = gsf_infile_stdio_new(path, gerr);

    if (infile == NULL) {
	err = 1;
    } else {
	output = gsf_output_stdio_new(fname, gerr);
	if (output == NULL) {
	    err = 1;
	}
    }  
    
    if (!err) {
	outfile = gsf_outfile_zip_new(output, gerr);
	g_object_unref(G_OBJECT(output));
	if (outfile == NULL) {
	    err = 1;
	}
    }  

    if (!err) {
	ziproot = gsf_outfile_new_child(outfile, path, 1);
	if (ziproot == NULL) {
	    fprintf(stderr, "failed to create ziproot for '%s'\n", path);
	    err = 1;
	} else {
	    err = clone_to_zip(GSF_INPUT(infile), ziproot);
	}
    }

    if (outfile != NULL) {
	gsf_output_close(GSF_OUTPUT(outfile));
	g_object_unref(outfile);
    }

#if ZDEBUG
    fprintf(stderr, "*** gretl_make_zipfile: returning %d\n", err);
#endif

    return err;
}

static int clone_to_stdio (GsfInfile *in, GsfOutput *output)
{
    GsfInput *input = GSF_INPUT(in);
    int err = 0;

#if ZDEBUG
    fprintf(stderr, "clone_to_stdio: in = %p (%s), output = %p\n", (void *) in, 
	    gsf_input_name(input), (void *) output);
#endif

    if (gsf_input_size(input) > 0) {
	transcribe_gsf_data(input, output);
    } else {
	GsfOutfile *out = output ? GSF_OUTFILE(output) : NULL;
	int i, n = gsf_infile_num_children(in);
	char const *name;
	gboolean is_dir;
	GDateTime *modtime;

	for (i=0; i<n && !err; i++) {
	    GError *gerr = NULL;

	    input = gsf_infile_child_by_index(in, i);
	    if (input == NULL) {
		err = 1;
		break;
	    }

	    name = gsf_infile_name_by_index(in, i);
	    is_dir = gsf_is_dir(input);

	    if (out == NULL) {
		if (is_dir) {
		    out = gsf_outfile_stdio_new(name, &gerr);
		    output = GSF_OUTPUT(out);
		    err = clone_to_stdio(GSF_INFILE(input), output);
		    g_object_unref(input);
		    gsf_output_close(output);
		    g_object_unref(output);
		    output = NULL;
		    out = NULL;
		} else {
		    output = gsf_output_stdio_new(name, &gerr);
		    err = clone_to_stdio(GSF_INFILE(input), output);
		    g_object_unref(input);
		    gsf_output_close(output);
		    g_object_unref(output);
		    output = NULL;
		    out = NULL;
		}
	    } else {
		modtime = gsf_input_get_modtime(input);
		output = gsf_outfile_new_child_full(out, name, is_dir, 
						    "modtime", modtime,
						    NULL);
		err = clone_to_stdio(GSF_INFILE(input), output);
		g_object_unref(input);
		gsf_output_close(output);
		g_object_unref(output);
		output = NULL;
	    }

	    if (gerr != NULL) {
		fprintf(stderr, "error: %s\n", gerr->message);
		g_error_free(gerr);
	    }
	}
    }

    return err;
}

int gretl_unzip_file (const char *fname, GError **gerr)
{
    GsfInput *input;
    GsfInfile *infile;
    int err = 0;

    ensure_gsf_init();
    input = gsf_input_stdio_new(fname, gerr);

#if ZDEBUG
    fprintf(stderr, "*** gretl_unzip_file (gsf): fname='%s', input=%p\n", 
	    fname, (void *) input);
#endif

    if (input == NULL) {
	err = 1;
    } else {
	infile = gsf_infile_zip_new(input, gerr);
	g_object_unref(G_OBJECT(input));
	if (infile == NULL) {
	    err = 1;
	} else {
	    err = clone_to_stdio(infile, NULL);
	    g_object_unref(infile);
	}
    }

#if ZDEBUG
    fprintf(stderr, "*** gretl_unzip_file: returning %d\n", err);
#endif

    return err;
}

int gretl_unzip_session_file (const char *fname, gchar **zdirname, GError **gerr)
{
    GsfInput *input;
    GsfInfile *infile;
    int err = 0;

    ensure_gsf_init();
    input = gsf_input_stdio_new(fname, gerr);

    if (input == NULL) {
	err = 1;
    } else {
	infile = gsf_infile_zip_new(input, gerr);
	g_object_unref(G_OBJECT(input));
	if (infile == NULL) {
	    err = 1;
	} else {
	    *zdirname = g_strdup(gsf_infile_name_by_index(infile, 0));
	    err = clone_to_stdio(infile, NULL);
	    g_object_unref(infile);
	}
    }

#if ZDEBUG
    fprintf(stderr, "*** gretl_unzip_session_file: zdirname = '%s', err = %d\n", 
	    *zdirname, err);
#endif

    return err;
}

#else /* native, using gretlzip plugin, not libgsf */

#include "plugins.h"

int gretl_unzip_file (const char *fname, GError **gerr)
{
    int (*zfunc) (const char *, GError **);
    void *handle = NULL;
    int err = 0;

    zfunc = get_plugin_function("gretl_native_unzip_file", &handle);
    if (zfunc == NULL) {
	/* error message handled by get_plugin_function() */
        return 1;
    }

    err = (*zfunc)(fname, gerr);

    close_plugin(handle);    

    return err;
}

/*
 * @fname: full path to zipfile to be created.
 * @path: path relative to userdir for files to be picked up
 * and zipped.
 */

int gretl_make_zipfile (const char *fname, const char *path,
			GError **gerr)
{
    int (*zfunc) (const char *, const char *, GError **);
    void *handle = NULL;
    int err = 0;

    zfunc = get_plugin_function("gretl_native_make_zipfile", &handle);
    if (zfunc == NULL) {
        return 1;
    }

    err = (*zfunc)(fname, path, gerr);

    close_plugin(handle);    

    return err;
}

int gretl_unzip_session_file (const char *fname, gchar **zdirname, GError **gerr)
{
    int (*zfunc) (const char *, gchar **, GError **);
    void *handle = NULL;
    int err = 0;

    zfunc = get_plugin_function("gretl_native_unzip_session_file", &handle);
    if (zfunc == NULL) {
        return 1;
    }

    err = (*zfunc)(fname, zdirname, gerr);

    close_plugin(handle);    

    return err;
}

#endif


