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

static gboolean font_is_latin_monospaced (PangoFontDescription *desc)
{
    int x1, x2;

    pango_context_set_font_description(font_test_context, desc);

    pango_layout_set_text(font_test_layout, "i", 1);
    pango_layout_get_pixel_size(font_test_layout, &x1, NULL);
    pango_layout_set_text(font_test_layout, "W", 1);
    pango_layout_get_pixel_size(font_test_layout, &x2, NULL);

    return (x1 == x2);
}

static gboolean font_is_latin_text_font (PangoFontDescription *desc)
{
    PangoCoverage *coverage = NULL;
    PangoFont *pfont = NULL;
    gboolean ok = FALSE;

    pfont = pango_context_load_font(font_test_context, desc);
    if (pfont == NULL) return FALSE;

    coverage = pango_font_get_coverage(pfont, font_test_lang);
    if (coverage == NULL) return FALSE;

    if (pango_coverage_get(coverage, 'i') == PANGO_COVERAGE_EXACT &&
	pango_coverage_get(coverage, 'W') == PANGO_COVERAGE_EXACT) {
	ok = TRUE;
    } 

    return ok;
}

static gboolean weird_font (const gchar *name)
{
    if (strstr(name, "arioso") || strstr(name, "dings") || strstr(name, "chancery"))
	return TRUE;
    else
	return FALSE;
}

struct ok_font {
    gchar *name;
    int mono;
};

static struct ok_font *ok_fonts;
static int n_ok;

static int push_ok_font (const char *name, int mono)
{
    struct ok_font *tmp;

    tmp = realloc(ok_fonts, (n_ok + 1) * sizeof *ok_fonts);
    if (tmp == NULL) {
	return 1;
    }

    ok_fonts = tmp;
    ok_fonts[n_ok].name = g_strdup(name);
    ok_fonts[n_ok].mono = mono;

    n_ok++;
    
    return 0;
}

static int font_is_ok (const char *name, int mono)
{
    int i;

    for (i=0; i<n_ok; i++) {
	if (ok_fonts[i].mono == mono &&
	    !strcmp(name, ok_fonts[i].name)) {
	    return 1;
	}
    }

    return 0;
}

/* We can test for a font for "latin text" compatibility, via the heuristic
   of seeing if it contains the letters 'i' and 'W' in English.  Given the
   latin text characteristic, we can then see if the font is monospaced
   by checking whether or not the letters 'i' and 'W' come out the
   same width. */

#define SHOW_PROGRESS 1

static gboolean validate_font_family (const gchar *familyname, 
				      gint filter,
				      gint n_families,
				      gboolean cache_built,
				      int *err)
{
    static void *handle = NULL;
    static int show = 1, n_done = 0;
    static int (*show_progress) (long, long, int) = NULL;

    int latin = 0;
    int mono = 0;

    if (!cache_built) {
	gchar *fontname = NULL;
	PangoFontDescription *desc = NULL;

#if FDEBUG
	fprintf(dbg, "Checking font family '%s'\n", familyname);
	fflush(dbg);
#endif

	if (show && handle == NULL) {
	    show_progress = gui_get_plugin_function("show_progress", 
						    &handle);
	    if (show_progress == NULL) {
		show = 0;
	    } else {
		(*show_progress)(0L, n_families, SP_FONT_INIT);
	    }
	}
	if (show && n_done > 0 && n_done < n_families) {
	    (*show_progress)(1L, n_families, SP_NONE);
	}
	n_done++;

	if (show && (n_done == n_families - 1)) {
	    (*show_progress)(0L, n_families, SP_FINISH);
	    show = 0;
	    close_plugin(handle);
	    handle = NULL;
	}

	if (weird_font(familyname)) {
	    return FALSE;
	}

	fontname = g_strdup_printf("%s 10", familyname);
	desc = pango_font_description_from_string(fontname);

	if (desc != NULL) {
#if FDEBUG
	    fprintf(dbg, "Testing %s... ", fontname);
	    fflush(dbg);
#endif

	    latin = font_is_latin_text_font(desc);
	    if (latin) {
		mono = font_is_latin_monospaced(desc);
		*err = push_ok_font(familyname, mono);
	    }

#if FDEBUG
	    fprintf(dbg, "Latin %s, monospaced %s\n",
		    latin? "yes" : "no", mono? "yes" : "no");
	    fflush(dbg);
#endif

	    pango_font_description_free(desc);
	}

	g_free(fontname);

    } else if (filter != GTK_FONT_HACK_NONE) { 
#if FDEBUG
	fprintf(dbg, "Looking up cached info on '%s'\n", familyname);
	fflush(dbg);
#endif
	return font_is_ok(familyname, filter == GTK_FONT_HACK_LATIN_MONO);
    }

    if (filter == GTK_FONT_HACK_LATIN) {
	return latin;
    } else if (filter == GTK_FONT_HACK_LATIN_MONO) {
	return mono;
    } else if (filter == GTK_FONT_HACK_NONE) {
	return TRUE;
    }

    return FALSE; /* not reachable */
}

