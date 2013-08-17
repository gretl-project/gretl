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
	    g_warning("error reading ?");
	    return 1;
	}
	if (!gsf_output_write(output, len, data)) {
	    g_warning("error writing ?");
	    return 1;
	}
    }

    return 0;
}

static void clone_recursive (GsfInput *input, GsfOutput *output)
{
    if (gsf_input_size(input) > 0) {
	transcribe_gsf_data(input, output);
    } 

    if (GSF_IS_INFILE(input) && gsf_infile_num_children(GSF_INFILE(input)) > 0) {
	GsfInfile *in = GSF_INFILE(input);
	GsfOutfile *out = GSF_OUTFILE(output);
	GsfInput *src;
	GsfOutput *dest;
	gboolean is_dir;
	GDateTime *modtime;
	int i;

	for (i=0; i<gsf_infile_num_children(in); i++) {
#if ZDEBUG
	    fprintf(stderr, "gsf_clone_recursive: i = %d (%s)\n", i, 
		    gsf_infile_name_by_index(in, i));
#endif
	    src = gsf_infile_child_by_index(in, i);
	    is_dir = gsf_is_dir(src);
	    modtime = gsf_input_get_modtime(src);
	    dest = gsf_outfile_new_child_full(out,
					      gsf_infile_name_by_index(in, i),
					      is_dir,
					      "modtime", modtime,
					      NULL);
	    clone_recursive(src, dest);
	}
    }

    gsf_output_close(output);
    g_object_unref(G_OBJECT(output));
    g_object_unref(G_OBJECT(input));
}

int gretl_make_zipfile (const char *fname, const char *path,
			GError **gerr)
{
    GsfInfile *infile;
    GsfOutput *output, *ziproot;
    GsfOutfile *outfile = NULL;
    int ok = 1;

    ensure_gsf_init();

#if ZDEBUG
    fprintf(stderr, "gretl_make_zipfile (gsf); fname='%s', path='%s'\n", fname, path);
#endif

    infile = gsf_infile_stdio_new(path, gerr);

    if (infile == NULL) {
	ok = 0;
    }

    if (ok) {
	output = gsf_output_stdio_new(fname, gerr);
	if (output == NULL) {
	    ok = 0;
	}
    }  
    
    if (ok) {
	outfile = gsf_outfile_zip_new(output, gerr);
	g_object_unref(G_OBJECT(output));
	if (outfile == NULL) {
	    ok = 0;
	}
    }  

    if (ok) {
	ziproot = gsf_outfile_new_child(outfile, path, 1);
	if (ziproot == NULL) {
	    g_warning("failed to create ziproot for '%s'\n", path);
	    ok = 0;
	} else {
	    clone_recursive(GSF_INPUT(infile), ziproot);
	}
    }

    if (outfile != NULL) {
	gsf_output_close(GSF_OUTPUT(outfile));
	g_object_unref(outfile);
    }

    fprintf(stderr, "*** gretl_make_zipfile: returning %d\n", !ok);

    return !ok;
}

int gretl_unzip_file (const char *fname, GError **gerr)
{
    GsfInput *input;
    GsfInfile *infile;
    GsfOutfile *outfile = NULL;
    int ok = 1;

    ensure_gsf_init();

    input = gsf_input_stdio_new(fname, gerr);

    fprintf(stderr, "*** gretl_unzip_file (gsf): fname='%s', input=%p\n", 
	    fname, (void *) input);

    if (input == NULL) {
	ok = 0;
    } else {
	infile = gsf_infile_zip_new(input, gerr);
	g_object_unref(G_OBJECT(input));
	if (infile != NULL) {
	    /* FIXME!! */
	    outfile = gsf_outfile_stdio_new("foo", gerr);
	    if (outfile == NULL) {
		ok = 0;
	    } else {
		clone_recursive(GSF_INPUT(infile), GSF_OUTPUT(outfile));
	    }
	} else {
	    ok = 0;
	}
    }

#if 0
    if (outfile != NULL) { /* ?? */
	gsf_output_close(GSF_OUTPUT(outfile));
	g_object_unref(outfile);
    }
#endif

    fprintf(stderr, "*** gretl_unzip_file: returning %d\n", !ok);

    return !ok;
}

int gretl_unzip_session_file (const char *fname, gchar **zdirname, GError **gerr)
{
    GsfInput *input;
    GsfInfile *infile;
    GsfOutfile *outfile = NULL;
    int ok = 1;

    check_gsf_init();

    input = gsf_input_stdio_new(fname, gerr);
    fprintf(stderr, "*** gsf: input = %p\n", (void *) input);

    if (input == NULL) {
	ok = 0;
    } else {
	infile = gsf_infile_zip_new(input, gerr);
	g_object_unref(G_OBJECT(input));

	if (infile == NULL) {
	    ok = 0;
	} else {
	    *zdirname = g_strdup(gsf_infile_name_by_index(infile, 0));
	    /* FIXME!! */
	    clone_recursive(GSF_INPUT(infile), NULL);
	    g_object_unref(infile);
	}
    }

    fprintf(stderr, "*** gretl_unzip_session_file: zdirname = '%s', ok = %d\n", 
	    *zdirname, ok);

    return !ok;
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


