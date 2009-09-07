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

#undef USE_LIBGSF

#ifdef USE_LIBGSF

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

#define CHUNK 4096

static int gsf_transcribe_data (GsfInput *input, GsfOutput *output)
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
	    g_warning ("error writing ?");
	    return 1;
	}
    }

    return 0;
}

static void clone_to_stdio (GsfInfile *in, GsfOutput *output)
{
    GsfInput *input = GSF_INPUT(in);

    if (gsf_input_size(input) > 0) {
	gsf_transcribe_data(input, output);
    } else {
	GsfOutfile *out = GSF_OUTFILE(output);
	int i, n = gsf_infile_num_children(in);
	char const *name;
	char *display_name;
	GError *err = NULL;
	time_t mtime;

#if ZDEBUG
	g_print("*** gsf_infile_num_children = %d\n", n);
#endif

	for (i=0 ; i<n; i++) {
	    gboolean is_dir;
	    
	    name = gsf_infile_name_by_index(in, i);
	    display_name = (name != NULL)? g_filename_display_name(name) : NULL;

	    input = gsf_infile_child_by_index(in, i);
	    if (input == NULL) {
		g_print("Error opening '%s', index = %d\n",
			display_name ? display_name : "?", i);
		g_free(display_name);
		continue;
	    } else {
		g_print("Opened '%s', index = %d\n",
			display_name ? display_name : "?", i);
	    }		

	    is_dir = gsf_infile_num_children(GSF_INFILE(input)) >= 0;

	    g_print("%s: size=%ld, %s\n",
		    display_name ? display_name : "?",
		    (long) gsf_input_size (input),
		    is_dir ? "directory" : "file");
	    g_free(display_name);

	    if (output == NULL) {
		out = gsf_outfile_stdio_new(name, &err);
		output = GSF_OUTPUT(out);
	    } else {    
		g_object_get(G_OBJECT(input), "modification-time", &mtime, NULL);
		output = gsf_outfile_new_child_full(out, name, is_dir, 
						    "modification-time", mtime,
						    NULL);
	    }

	    clone_to_stdio(GSF_INFILE(input), output);

	    g_object_unref(input);
	    gsf_output_close(output);
	    g_object_unref(output);
	}
    }
}

#define gsf_is_dir(i) (GSF_IS_INFILE(i) && gsf_infile_num_children(GSF_INFILE(i)) >= 0)

static void clone_to_zip (GsfInput *input, GsfOutput *output)
{
    if (gsf_input_size(input) > 0) {
	gsf_transcribe_data(input, output);
    } 

    if (GSF_IS_INFILE(input) && gsf_infile_num_children(GSF_INFILE(input)) > 0) {
	GsfInfile *in = GSF_INFILE(input);
	GsfOutfile *out = GSF_OUTFILE(output);
	GsfInput *src;
	GsfOutput *dst;
	gboolean is_dir;
	time_t mtime;
	int i;

	for (i=0; i<gsf_infile_num_children(in); i++) {
	    src = gsf_infile_child_by_index(in, i);
	    is_dir = gsf_is_dir(src);
	    g_object_get(G_OBJECT(src), "modification-time", &mtime, NULL);
	    dst = gsf_outfile_new_child_full(out,
					     gsf_infile_name_by_index(in, i),
					     is_dir,
					     "modification-time", mtime,
					     NULL);
	    clone_to_zip(src, dst);
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

    gsf_init();

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
	ziproot = gsf_outfile_new_child(outfile, fname, 1);
	if (ziproot == NULL) {
	    g_warning("failed to create ziproot for '%s'\n", fname);
	    ok = 0;
	} else {
	    clone_to_zip(GSF_INPUT(infile), ziproot);
	}
    }

    if (outfile != NULL) {
	gsf_output_close(GSF_OUTPUT(outfile));
	g_object_unref(outfile);
    }

    gsf_shutdown();

    return !ok;
}

int gretl_unzip_file (const char *fname, GError **gerr)
{
    GsfInput *input;
    GsfInfile *infile;
    int ok = 1;

    fprintf(stderr, "*** gretl_unzip_file: calling gsf_init\n");

    gsf_init();

    input = gsf_input_stdio_new(fname, gerr);

    fprintf(stderr, "*** gretl_unzip_file: input = %p\n", (void *) input);

    if (input == NULL) {
	ok = 0;
    } else {
	infile = gsf_infile_zip_new(input, gerr);
	g_object_unref(G_OBJECT(input));
	if (infile == NULL) {
	    clone_to_stdio(infile, NULL);
	    g_object_unref(infile);
	} else {
	    ok = 0;
	}
    }

    gsf_shutdown();

    fprintf(stderr, "*** gretl_unzip_file: returning %d\n", !ok);

    return !ok;
}

gchar *gretl_zipfile_get_topdir (const char *fname)
{
    GsfInput *input;
    GsfInfile *infile;
    GError *err = NULL;
    gchar *topdir = NULL;

    gsf_init();

    input = gsf_input_stdio_new(fname, &err);

    if (input == NULL) {
	g_return_val_if_fail(err != NULL, 0);
	g_warning("'%s' error: %s", fname, err->message);
	g_error_free(err);
    } else {
	infile = gsf_infile_zip_new(input, &err);
	g_object_unref(G_OBJECT(input));

	if (infile == NULL) {
	    g_error_free(err);
	} else {
	    topdir = g_strdup(gsf_infile_name_by_index(infile, 0));
	    g_object_unref(infile);
	}
    }

    gsf_shutdown();

    fprintf(stderr, "*** gretl_zipfile_get_topdir: returning '%s'\n", topdir);

    return topdir;
}

#else

#include "zipunzip.h"

int gretl_unzip_file (const char *fname, GError **gerr)
{
    /* for verbose operation, make 3rd arg ZIP_VERBOSE or ZIP_TRACE */
    return zipfile_extract_files(fname, NULL, 0, gerr);
}

/*
 * @fname: full path to zipfile to be created.
 * @path: path relative to userdir for files to be picked up
 * and zipped.
 */

int gretl_make_zipfile (const char *fname, const char *path,
			GError **gerr)
{
    const char *array[2] = { path, NULL };
    int err;

    err = zipfile_archive_files(fname, array, 9, 
				ZIP_RECURSE_DIRS,
				gerr);
    return err;
}

gchar *gretl_zipfile_get_topdir (const char *fname)
{
    zipinfo *zinfo;
    gchar *topdir = NULL;

    zinfo = zipfile_get_info(fname, 0, NULL);

    if (zinfo != NULL) {
	int i, n, gotit = 0;
	const gchar *s;

	for (i=0; i<zinfo->nfiles && !gotit; i++) {
	    s = zinfo->fnames[i];
	    if (s != NULL) {
		n = strlen(s);
		if (n > 13 && !strcmp(s + n - 11, "session.xml")) {
		    topdir = g_strndup(s, n - 11);
		    if (topdir != NULL) {
			n = strlen(topdir);
			if (topdir[n-1] == '/' || topdir[n-1] == '\\') {
			    topdir[n-1] = '\0';
			}
		    }
		}
	    }
	}
	zipinfo_destroy(zinfo);
    }

    return topdir;
}

#endif


