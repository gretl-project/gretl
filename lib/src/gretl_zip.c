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

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#ifdef WIN32
# include "winconfig.h"
#endif

#include "gretl_zip.h"
#include "strutils.h"

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

#define CHUNK 32768

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

static int clone_gsf_tree (GsfInput *input, GsfOutput *output)
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
	    if (src == NULL) {
		err = 1;
		break;
	    }

	    name = gsf_infile_name_by_index(in, i);
	    is_dir = gsf_is_dir(src);
#if ZDEBUG
	    fprintf(stderr, "clone_gsf_tree: i = %d (%s, %s)\n", i, name, 
		    is_dir ? "directory" : "file"); 
#endif
	    modtime = gsf_input_get_modtime(src);
	    dest = gsf_outfile_new_child_full(out, name, is_dir,
					      "modtime", modtime,
					      NULL);
	    if (dest == NULL) {
		err = 1;
	    } else {
		err = clone_gsf_tree(src, dest);
	    }
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
    fprintf(stderr, "gretl_make_zipfile (gsf):\n fname='%s'\n path='%s'\n", 
	    fname, path);
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
	    err = clone_gsf_tree(GSF_INPUT(infile), ziproot);
	}
    }

    if (outfile != NULL) {
	gsf_output_close(GSF_OUTPUT(outfile));
	g_object_unref(G_OBJECT(outfile));
    }

#if ZDEBUG
    fprintf(stderr, "*** gretl_make_zipfile: returning %d\n", err);
#endif

    return err;
}

#if 0 /* libgsf not ready yet: doesn't support setting of
	 zlib compression level as of version 1.14.30
      */

int gretl_zip_datafile (const char *fname, const char *path,
			int level, GError **gerr)
{
    GsfOutput *output;
    GsfOutfile *outfile = NULL;
    int err = 0;

    ensure_gsf_init();

#if ZDEBUG
    fprintf(stderr, "gretl_zip_datafile (gsf):\n fname='%s'\n", fname);
#endif

    output = gsf_output_stdio_new(fname, gerr);
    if (output == NULL) {
	err = 1;
    }
     
    if (!err) {
	outfile = gsf_outfile_zip_new(output, gerr);
	g_object_unref(G_OBJECT(output));
	if (outfile == NULL) {
	    err = 1;
	}
    } 

    if (!err) {
	const char *names[] = { "data.xml", "data.bin" };
	char fullname[FILENAME_MAX];
	GsfInput *zinp;
	GsfOutput *zout;
	int i;

	for (i=0; i<2 && !err; i++) {
	    build_path(fullname, path, names[i], NULL);
	    zinp = gsf_input_stdio_new(fullname, gerr);
	    if (zinp == NULL) {
		err = 1;
	    } else {
		/* property not present in libgsf <= 1.14.30 */
		zout = gsf_outfile_new_child_full(outfile, names[i], 0,
						  "deflate-level", level,
						  NULL);
		err = transcribe_gsf_data(zinp, zout);
		gsf_output_close(zout);
		g_object_unref(G_OBJECT(zout));
		g_object_unref(G_OBJECT(zinp));
	    }
	}
    }

    if (outfile != NULL) {
	gsf_output_close(GSF_OUTPUT(outfile));
	g_object_unref(G_OBJECT(outfile));
    }

#if ZDEBUG
    fprintf(stderr, "*** gretl_zip_datafile: returning %d\n", err);
#endif

    return err;
}

#else

/* use native gretlzip plugin code for now */

int gretl_zip_datafile (const char *fname, const char *path,
			int level, GError **gerr)
{
    int (*zfunc) (const char *, const char *, int, GError **);
    void *handle = NULL;
    int err = 0;

    zfunc = get_plugin_function("gretl_native_zip_datafile", &handle);
    if (zfunc == NULL) {
        return 1;
    }

    err = (*zfunc)(fname, path, level, gerr);

    close_plugin(handle);    

    return err;
}

#endif /* work-around for missing functionality in libgsf */

static int real_unzip_file (const char *fname, const char *path,
			    gchar **zdirname, GError **gerr)
{
    GsfInput *input;
    GsfInfile *infile;
    GsfOutfile *outfile;
    int err = 0;

#if ZDEBUG
    fprintf(stderr, "*** real_unzip_file (gsf):\n fname='%s'\n", fname);
#endif

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
	    if (zdirname != NULL) {
		*zdirname = g_strdup(gsf_infile_name_by_index(infile, 0));
	    }
	    if (path != NULL) {
		outfile = gsf_outfile_stdio_new(path, gerr);
	    } else {
		outfile = gsf_outfile_stdio_new(".", gerr);
	    }
	    if (outfile == NULL) {
		err = 1;
	    } else {
		err = clone_gsf_tree(GSF_INPUT(infile), GSF_OUTPUT(outfile));
	    }
	}
    }

#if ZDEBUG
    fprintf(stderr, "*** real_unzip_file: returning %d\n", err);
#endif

    return err;
}

int gretl_unzip_file (const char *fname, GError **gerr)
{
    return real_unzip_file(fname, NULL, NULL, gerr);
}

int gretl_unzip_session_file (const char *fname, gchar **zdirname, GError **gerr)
{
    return real_unzip_file(fname, NULL, zdirname, gerr);
}

int gretl_unzip_datafile (const char *fname, const char *path, GError **gerr)
{
    return real_unzip_file(fname, path, NULL, gerr);
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

int gretl_zip_datafile (const char *fname, const char *path,
			int level, GError **gerr)
{
    int (*zfunc) (const char *, const char *, int, GError **);
    void *handle = NULL;
    int err = 0;

    zfunc = get_plugin_function("gretl_native_zip_datafile", &handle);
    if (zfunc == NULL) {
        return 1;
    }

    err = (*zfunc)(fname, path, level, gerr);

    close_plugin(handle);    

    return err;
}

int gretl_unzip_datafile (const char *fname, const char *path,
			  GError **gerr)
{
    int (*zfunc) (const char *, const char *, GError **);
    void *handle = NULL;
    int err = 0;

    zfunc = get_plugin_function("gretl_native_unzip_datafile", &handle);
    if (zfunc == NULL) {
        return 1;
    }

    err = (*zfunc)(fname, path, gerr);

    close_plugin(handle);    

    return err;
}

#endif /* zip/unzip variants */


