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

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#ifdef WIN32
# include "winconfig.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <glib.h>

#include "gretl_zip.h"
#include "strutils.h"

static int handle_zip_error (const char *fname,
			     GError *gerr, int err,
			     const char *action)
{
    if (gerr != NULL) {
	gretl_errmsg_sprintf("%s: %s", fname, gerr->message);
	g_error_free(gerr);
    } else if (err) {
	gretl_errmsg_sprintf("%s: error %s", fname, action);
    }

    return err;
}

#ifdef USE_GSF /* libgsf-1 >= 1.14.31 */

#include <gsf/gsf.h>

#define ZDEBUG 0
#define CHUNK 32768

#define gsf_is_dir(i) (GSF_IS_INFILE(i) && \
		       gsf_infile_num_children(GSF_INFILE(i)) >= 0)

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

static int gretl_gsf_make_zipfile (const char *fname,
				   const char *path)
{
    GsfInfile *infile;
    GsfOutput *output, *ziproot;
    GsfOutfile *outfile = NULL;
    GError *gerr = NULL;
    int err = 0;

    ensure_gsf_init();

#if ZDEBUG
    fprintf(stderr, "gretl_make_zipfile (gsf):\n fname='%s'\n path='%s'\n", 
	    fname, path);
#endif

    infile = gsf_infile_stdio_new(path, &gerr);

    if (infile == NULL) {
	err = 1;
    } else {
	output = gsf_output_stdio_new(fname, &gerr);
	if (output == NULL) {
	    err = 1;
	}
    }  
    
    if (!err) {
	outfile = gsf_outfile_zip_new(output, &gerr);
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

    return handle_zip_error(fname, gerr, err, "zipping");
}

static int gretl_gsf_zip_datafile (const char *fname,
				   const char *path,
				   int level)
{
    GsfOutput *output;
    GsfOutfile *outfile = NULL;
    GError *gerr = NULL;
    int err = 0;

    ensure_gsf_init();

#if ZDEBUG
    fprintf(stderr, "gretl_zip_datafile (gsf):\n fname='%s'\n", fname);
#endif

    output = gsf_output_stdio_new(fname, &gerr);
    if (output == NULL) {
	err = 1;
    }
     
    if (!err) {
	outfile = gsf_outfile_zip_new(output, &gerr);
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
	    zinp = gsf_input_stdio_new(fullname, &gerr);
	    if (zinp == NULL) {
		err = 1;
	    } else {
		/* note: property not present in libgsf <= 1.14.30 */
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
    fprintf(stderr, "*** gretl_zip_datafile (gsf): returning %d\n", err);
#endif

    return handle_zip_error(fname, gerr, err, "zipping");
}

/* @fname: the name of the file to be unzipped.
   @path: if non-NULL, the location into which the unzipping
   should be done (if NULL, use pwd).
   @zdirname: if non-NULL, set its content to the name of
   the top-level directory within the zipfile (via g_strdup).
*/

static int gretl_gsf_unzip (const char *fname,
			    const char *path,
			    gchar **zdirname)
{
    GsfInput *input;
    GsfInfile *infile;
    GsfOutfile *outfile;
    GError *gerr;
    int err = 0;

#if ZDEBUG
    fprintf(stderr, "*** gretl_gsf_unzip_file:\n fname='%s'\n", fname);
#endif

    ensure_gsf_init();
    input = gsf_input_stdio_new(fname, &gerr);

    if (input == NULL) {
	err = 1;
    } else {
	infile = gsf_infile_zip_new(input, &gerr);
	g_object_unref(G_OBJECT(input));
	if (infile == NULL) {
	    err = 1;
	} else {
	    if (zdirname != NULL) {
		*zdirname = g_strdup(gsf_infile_name_by_index(infile, 0));
	    }
	    if (path != NULL) {
		outfile = gsf_outfile_stdio_new(path, &gerr);
	    } else {
		outfile = gsf_outfile_stdio_new(".", &gerr);
	    }
	    if (outfile == NULL) {
		err = 1;
	    } else {
		err = clone_gsf_tree(GSF_INPUT(infile), GSF_OUTPUT(outfile));
	    }
	}
    }

#if ZDEBUG
    fprintf(stderr, "*** gretl_gsf_unzip_file: returning %d\n", err);
#endif

    return handle_zip_error(fname, gerr, err, "unzipping");
}

#else /* native, using gretlzip plugin, not libgsf */

#include "plugins.h"

static int gretl_plugin_unzip (const char *fname)
{
    int (*zfunc) (const char *, GError **);
    GError *gerr = NULL;
    int err = 0;

    zfunc = get_plugin_function("gretl_native_unzip");
    if (zfunc == NULL) {
	/* error message handled by get_plugin_function() */
        return 1;
    }

    err = (*zfunc)(fname, &gerr);

    return handle_zip_error(fname, gerr, err, "unzipping");
}

/*
 * @fname: full path to zipfile to be created.
 * @path: path relative to userdir for files to be picked up
 * and zipped.
 */

static int gretl_plugin_make_zipfile (const char *fname,
				      const char *path)
{
    int (*zfunc) (const char *, const char *, GError **);
    GError *gerr = NULL;
    int err = 0;

    zfunc = get_plugin_function("gretl_native_make_zipfile");
    if (zfunc == NULL) {
        return 1;
    }

    err = (*zfunc)(fname, path, &gerr);

    return handle_zip_error(fname, gerr, err, "zipping");
}

static int gretl_plugin_unzip_session_file (const char *fname,
					    gchar **zdirname)
{
    int (*zfunc) (const char *, gchar **, GError **);
    GError *gerr = NULL;
    int err = 0;

    zfunc = get_plugin_function("gretl_native_unzip_session_file");
    if (zfunc == NULL) {
        return 1;
    }

    err = (*zfunc)(fname, zdirname, &gerr);

    return handle_zip_error(fname, gerr, err, "unzipping");
}

static int gretl_plugin_zip_datafile (const char *fname,
				      const char *path,
				      int level)
{
    int (*zfunc) (const char *, const char *, int, GError **);
    GError *gerr = NULL;
    int err = 0;

    zfunc = get_plugin_function("gretl_native_zip_datafile");
    if (zfunc == NULL) {
        return 1;
    }

    err = (*zfunc)(fname, path, level, &gerr);

    return handle_zip_error(fname, gerr, err, "zipping");
}

static int gretl_plugin_unzip_into (const char *fname,
				    const char *path)
{
    int (*zfunc) (const char *, const char *, GError **);
    GError *gerr = NULL;
    int err = 0;

    zfunc = get_plugin_function("gretl_native_unzip_into");
    if (zfunc == NULL) {
        return 1;
    }

    err = (*zfunc)(fname, path, &gerr);

    return handle_zip_error(fname, gerr, err, "unzipping");
}

#endif /* zip/unzip variants */

/**
 * gretl_unzip:
 * @fname: name of the file to unzip.
 *
 * Unzips @fname in the current directory, preserving any
 * internal directory structure.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_unzip (const char *fname)
{
#if USE_GSF
    return gretl_gsf_unzip(fname, NULL, NULL);
#else
    return gretl_plugin_unzip(fname);
#endif    
}

/**
 * gretl_unzip_into:
 * @fname: name of the file to unzip.
 * @dirname: the name of the directory in which unzipping
 * should take place.
 *
 * Unzips @fname in the specified directory, preserving any
 * internal directory structure.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_unzip_into (const char *fname, const char *dirname)
{
#if USE_GSF
    return gretl_gsf_unzip(fname, dirname, NULL);
#else
    return gretl_plugin_unzip_into(fname, dirname);
#endif 
}

/**
 * gretl_make zipfile:
 * @fname: name of the zip file to create.
 * @path: the path to the content which should be zipped.
 *
 * Creates a zip file of the specified name, whose content
 * is determined (recursively) by @path.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_make_zipfile (const char *fname, const char *path)
{
#if USE_GSF
    return gretl_gsf_make_zipfile(fname, path);    
#else
    return gretl_plugin_make_zipfile(fname, path);       
#endif 
}

/**
 * gretl_unzip_session_file:
 * @fname: name of the file to unzip.
 * @zdirname: location to receive the name of the top-level
 * directory within the zip archive.
 *
 * Specialized (slightly) unzipper for gretl session files.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_unzip_session_file (const char *fname, gchar **zdirname)
{
#if USE_GSF
    return gretl_gsf_unzip(fname, NULL, zdirname);
#else
    return gretl_plugin_unzip_session_file(fname, zdirname);    
#endif 
}

/**
 * gretl_zip_datafile:
 * @fname: name of the zip file to create.
 * @path: the path to the content which should be zipped.
 * @level: the zlib compression level to apply.
 *
 * Creates a zip file of the specified name, with content
 * given by @path, using the specified compression level.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_zip_datafile (const char *fname, const char *path,
			int level)
{
#if USE_GSF
    return gretl_gsf_zip_datafile(fname, path, level);
#else
    return gretl_plugin_zip_datafile(fname, path, level);
#endif 
}
