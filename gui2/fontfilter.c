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

#include "gretl.h"
#include "fontfilter.h"

#define FDEBUG 0

static GtkWidget *font_test_widget;
static PangoContext *font_test_context;
static PangoLanguage *font_test_lang;

int gretl_font_filter_init (void)
{
    font_test_widget = gtk_label_new(NULL);  
    font_test_context = gtk_widget_get_pango_context(font_test_widget); 

    if (font_test_context == NULL) {
	gtk_widget_destroy(font_test_widget);
	font_test_widget = NULL;
	return 1;
    }    

    font_test_lang = pango_language_from_string("eng");
    if (font_test_lang == NULL) {
	gretl_font_filter_cleanup();
	return 1;
    }

    return 0;
}

void gretl_font_filter_cleanup (void)
{
    if (font_test_context != NULL) {
	g_object_unref(G_OBJECT(font_test_context));
	font_test_context = NULL;
    }
    if (font_test_widget != NULL) {
	gtk_widget_destroy(font_test_widget);
	font_test_widget = NULL;
    }
}

enum {
    HACK_LATIN_FONT = 1 << 0,
    HACK_MONO_FONT  = 1 << 1,
    HACK_WEIRD_FONT = 1 << 2
};

/* We test a font for "latin text" compatibility, via the heuristic
   of seeing if it contains the letters 'i' and 'W' in English.  Given the
   latin text characteristic, we then see if the font is monospaced
   using pango_font_family_is_monospace().
*/

static int get_font_characteristics (PangoFontFamily *family,
				     PangoFontDescription *desc)
{
    PangoFont *font;
    int ret = 0;

    font = pango_context_load_font(font_test_context, desc);

    if (font != NULL) {
	PangoCoverage *coverage;

	coverage = pango_font_get_coverage(font, font_test_lang);
	if (coverage != NULL) {
	    if (pango_coverage_get(coverage, 'i') == PANGO_COVERAGE_EXACT &&
		pango_coverage_get(coverage, 'W') == PANGO_COVERAGE_EXACT) {
		ret = HACK_LATIN_FONT;
		if (pango_font_family_is_monospace(family)) {
		    ret |= HACK_MONO_FONT;
		}
	    }
	    pango_coverage_unref(coverage);
	}
	g_object_unref(font);
    }

#if FDEBUG
    fprintf(stderr, " latin %s, monospaced %s\n",
	    (ret & HACK_LATIN_FONT)? "yes" : "no",
	    (ret & HACK_MONO_FONT)? "yes" : "no");
#endif

    return ret;
}

static void font_progress_bar (int i, int nf)
{
    static int (*show_progress) (gint64, gint64, int) = NULL;
    static void *handle = NULL;
    static int show = 1;

    if (show && handle == NULL) {
	show_progress = gui_get_plugin_function("show_progress", 
						&handle);
	if (show_progress == NULL) {
	    show = 0;
	} else {
	    (*show_progress)(0L, nf, SP_FONT_INIT);
	}
    }

    if (show) {
	if (i == nf - 1) {
	    (*show_progress)(0, nf, SP_FINISH);
	    show = 0;
	    close_plugin(handle);
	    handle = NULL;
	} else if (i > 0) {
	    (*show_progress)(1, nf, SP_NONE);
	}
    }
}

#define weird_font(s) (strstr(s, "arioso") || \
                       strstr(s, "dings") || \
                       strstr(s, "ancery"))

/* If we're on the first run, check the font's characteristics and
   cache the info.  Otherwise just read off the information we cached
   previously.  In either case, return non-zero iff the font validates
   in respect of the criterion in @filter.
*/

int validate_font_family (PangoFontFamily *family,
			  const gchar *famname, 
			  gint i, gint nf,
			  gint filter, gint *err)
{
    static char *fcache;
    static int build_cache;
    int ret;

    if (i >= nf) {
	fprintf(stderr, "validate_font_family: got excess font!\n");
	return 0;
    }

    if (fcache == NULL) {
	fcache = calloc(nf, 1);
	if (fcache == NULL) {
	    *err = 1;
	    return 0;
	} 
	build_cache = 1;
	gretl_font_filter_init();
    }

    if (build_cache) {
#if FDEBUG
	fprintf(stderr, "Checking font family %d, '%s'\n", i, famname);
#endif
	font_progress_bar(i, nf);
	if (weird_font(famname)) {
	    fcache[i] = HACK_WEIRD_FONT;
	} else {
	    gchar *fontname = g_strdup_printf("%s 10", famname);
	    PangoFontDescription *desc = 
		pango_font_description_from_string(fontname);

	    if (desc != NULL) {
		fcache[i] = get_font_characteristics(family, desc);
		pango_font_description_free(desc);
	    }
	    g_free(fontname);
	}
    } 

    if (build_cache && i == nf - 1) {
	gretl_font_filter_cleanup();
	build_cache = 0;
    }

    if (filter == FONT_FILTER_LATIN) {
	ret = fcache[i] & HACK_LATIN_FONT;
    } else if (filter == FONT_FILTER_LATIN_MONO) {
	ret = fcache[i] & HACK_MONO_FONT;
    } else {
	ret = !(fcache[i] & HACK_WEIRD_FONT);
    }

    return ret;
}

int validate_single_font (const PangoFontFamily *family,
			  gint filter)
{
    const gchar *famname = 
	pango_font_family_get_name((PangoFontFamily *) family);
    int ret, fc = 0;

#if FDEBUG
    fprintf(stderr, "Checking font '%s'\n", famname);
#endif

    if (weird_font(famname)) {
	fc = HACK_WEIRD_FONT;
    } else {
	gchar *font = g_strdup_printf("%s 10", famname);
	PangoFontDescription *desc = pango_font_description_from_string(font);

	if (desc != NULL) {
	    fc = get_font_characteristics((PangoFontFamily *) family, desc);
	    pango_font_description_free(desc);
	}
	g_free(font);
    }


    if (filter == FONT_FILTER_LATIN) {
	ret = fc & HACK_LATIN_FONT;
    } else if (filter == FONT_FILTER_LATIN_MONO) {
	ret = fc & HACK_MONO_FONT;
    } else {
	ret = !(fc & HACK_WEIRD_FONT);
    }

    return ret;
}
