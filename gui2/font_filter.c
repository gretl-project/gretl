static GtkWidget *font_test_widget;
static PangoContext *font_test_context;
static PangoLayout *font_test_layout;
static PangoLanguage *font_test_lang;

static int create_font_test_rig (void)
{
    font_test_widget = gtk_label_new(NULL);  
    font_test_context = gtk_widget_get_pango_context(font_test_widget); 

    if (font_test_context == NULL) {
	gtk_widget_destroy(font_test_widget);
	return 1;
    }    

    font_test_layout = pango_layout_new(font_test_context); 
    if (font_test_layout == NULL) {
	g_object_unref(G_OBJECT(font_test_context));
	gtk_widget_destroy(font_test_widget);
	return 1;
    }	

    font_test_lang = pango_language_from_string("eng");

    return 0;
}

static void destroy_font_test_rig (void)
{
    g_object_unref(G_OBJECT(font_test_layout));
    g_object_unref(G_OBJECT(font_test_context));
    gtk_widget_destroy(font_test_widget);    
}

enum {
    LATIN_FONT = 1 << 0,
    MONO_FONT  = 1 << 1,
    WEIRD_FONT = 1 << 2
};

/* We test a font for "latin text" compatibility, via the heuristic
   of seeing if it contains the letters 'i' and 'W' in English.  Given the
   latin text characteristic, we then see if the font is monospaced
   by checking whether or not the letters 'i' and 'W' come out the
   same width. */

static int get_font_characteristics (PangoFontDescription *desc)
{
    PangoCoverage *coverage = NULL;
    PangoFont *pfont = NULL;
    int ret = 0;

    pfont = pango_context_load_font(font_test_context, desc);
    if (pfont == NULL) {
	return 0;
    }

    coverage = pango_font_get_coverage(pfont, font_test_lang);
    if (coverage == NULL) {
	return 0;
    }

    if (pango_coverage_get(coverage, 'i') == PANGO_COVERAGE_EXACT &&
	pango_coverage_get(coverage, 'W') == PANGO_COVERAGE_EXACT) {
	int x1, x2;

	ret = LATIN_FONT;
	
	if (font_test_context != NULL && font_test_layout != NULL) {
	    pango_context_set_font_description(font_test_context, desc);
	    pango_layout_set_text(font_test_layout, "i", 1);
	    pango_layout_get_pixel_size(font_test_layout, &x1, NULL);
	    pango_layout_set_text(font_test_layout, "W", 1);
	    pango_layout_get_pixel_size(font_test_layout, &x2, NULL);

	    if (x1 == x2) {
		ret |= MONO_FONT;
	    }
	}
    }

    pango_coverage_unref(coverage);

#if FDEBUG
    fprintf(stderr, " latin %s, monospaced %s\n",
	    (ret & LATIN_FONT)? "yes" : "no",
	    (ret & MONO_FONT)? "yes" : "no");
#endif

    return ret;
}

int font_has_minus (PangoFontDescription *desc)
{
    GtkWidget *widget;
    PangoContext *context = NULL;
    PangoLayout *layout = NULL;
    PangoLanguage *lang = NULL;
    PangoCoverage *coverage = NULL;
    PangoFont *pfont = NULL;
    int ret = 0;

    widget = gtk_label_new(NULL);  
    context = gtk_widget_get_pango_context(widget); 

    if (context == NULL) {
	gtk_widget_destroy(widget);
	return 0;
    }    

    layout = pango_layout_new(context); 
    lang = pango_language_from_string("eng");

    if (layout != NULL && lang != NULL) {
	pfont = pango_context_load_font(context, desc);
	if (pfont != NULL) {
	    coverage = pango_font_get_coverage(pfont, lang);
	    if (coverage != NULL) {
		/* U+2212 = minus sign */
		ret = (pango_coverage_get(coverage, 0x2212) == PANGO_COVERAGE_EXACT);
	    }
	}
    } 

    pango_coverage_unref(coverage);
    g_object_unref(G_OBJECT(layout));
    g_object_unref(G_OBJECT(context));
    gtk_widget_destroy(widget);    

    return ret;
}

static void font_progress_bar (int i, int nf)
{
    static int (*show_progress) (long, long, int) = NULL;
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
	    (*show_progress)(0L, nf, SP_FINISH);
	    show = 0;
	    close_plugin(handle);
	    handle = NULL;
	} else if (i > 0) {
	    (*show_progress)(1L, nf, SP_NONE);
	}
    }
}

#define weird_font(s) (strstr(s, "arioso") || \
                       strstr(s, "dings") || \
                       strstr(s, "ancery"))

/* If we're on the first run, check the font's characteristics and
   cache the info.  Otherwise just read off the information we cached
   previously.  In either case, return non-zero iff the font validates
   in respect of the criterion in filter.
*/

static int validate_font_family (const gchar *famname, 
				 gint i, gint nf,
				 gint filter, gint *err)
{
    static char *fcache;
    static int build_cache;
    int ret;

    if (fcache == NULL) {
	fcache = calloc(nf, 1);
	if (fcache == NULL) {
	    *err = 1;
	    return 1;
	} 
	build_cache = 1;
	create_font_test_rig();
    }

    if (build_cache) {
#if FDEBUG
	fprintf(stderr, "Checking font family %d, '%s'\n", i, famname);
#endif

	font_progress_bar(i, nf);

	if (weird_font(famname)) {
	    fcache[i] = WEIRD_FONT;
	} else {
	    gchar *font = g_strdup_printf("%s 10", famname);
	    PangoFontDescription *desc = pango_font_description_from_string(font);

	    if (desc != NULL) {
		fcache[i] = get_font_characteristics(desc);
		pango_font_description_free(desc);
	    }
	    g_free(font);
	}
    } 

    if (build_cache && i == nf - 1) {
	destroy_font_test_rig();
	build_cache = 0;
    }

    ret = (filter == GTK_FONT_HACK_LATIN)? (fcache[i] & LATIN_FONT) :
	(filter == GTK_FONT_HACK_LATIN_MONO)? (fcache[i] & MONO_FONT) :
	!(fcache[i] & WEIRD_FONT);

    return ret;
}

