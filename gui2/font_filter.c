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

/* We can test for a font for "latin text" compatibility, via the heuristic
   of seeing if it contains the letters 'i' and 'W' in English.  Given the
   latin text characteristic, we can then see if the font is monospaced
   by checking whether or not the letters 'i' and 'W' come out the
   same width. */

#define SHOW_PROGRESS 1

static gboolean validate_font_family (const gchar *familyname, 
				      gint filter,
				      gint n_families,
				      gboolean cache_built)
{
    static gchar **latin_families;
    static gboolean *monospaced;
    static gint n_latin;
    gboolean is_latin = FALSE, is_mono = FALSE;
#ifdef SHOW_PROGRESS
    static void *handle = NULL;
    static int show = 1, n_done = 0;
    static int (*show_progress) (long, long, int) = NULL;
#endif

    if (!cache_built) { /* we haven't set up the cache yet */
	gchar *fontname = NULL;
	PangoFontDescription *desc = NULL;

#ifdef FONT_FILTER_DEBUG
	fprintf(dbg, "Checking font family '%s'\n", familyname);
	fflush(dbg);
#endif

#ifdef SHOW_PROGRESS
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
# ifdef FONT_FILTER_DEBUG
	fprintf(stderr, "n_families=%d, n_done=%d\n", n_families, n_done);
# endif
	if (show && (n_done == n_families - 1)) {
# ifdef FONT_FILTER_DEBUG
	    fprintf(stderr, "doing SP_FINISH\n");
# endif
	    (*show_progress)(0L, n_families, SP_FINISH);
	    show = 0;
	    close_plugin(handle);
	    handle = NULL;
	}
#endif

	if (weird_font(familyname)) return FALSE;

	fontname = g_strdup_printf("%s 10", familyname);
	desc = pango_font_description_from_string(fontname);

	if (desc != NULL) {
	    int memerr = 0;

# ifdef FONT_FILTER_DEBUG
	    fprintf(dbg, "Got pango_font_description for '%s'\n", fontname);
	    fprintf(dbg, "Doing font_is_latin_text_font() test\n");
	    fflush(dbg);
# endif
	    if (font_is_latin_text_font(desc)) {
# ifdef FONT_FILTER_DEBUG
		fprintf(dbg, "%s passed font_is_latin_text_font() test\n", fontname);
		fflush(dbg);
# endif
		/* extend the cache */
		if (latin_families == NULL) {
		    latin_families = malloc(sizeof *latin_families);
		} else {
		    latin_families = realloc(latin_families, (n_latin + 1) * 
					     sizeof *latin_families);
		}
		if (latin_families == NULL) memerr = 1;

		if (monospaced == NULL) {
		    monospaced = malloc(sizeof *monospaced);
		} else {
		    monospaced = realloc(monospaced, (n_latin + 1) * sizeof *monospaced);
		}
		if (monospaced == NULL) memerr = 1;

		if (!memerr) {
		    latin_families[n_latin] = g_strdup(familyname);
		    is_latin = TRUE;
#ifdef FONT_FILTER_DEBUG
		    fprintf(dbg, "Doing monospace test for '%s'\n", fontname);
		    fflush(dbg);
#endif
		    if (font_is_latin_monospaced(desc)) {
			is_mono = monospaced[n_latin] = TRUE;
		    } else {
			is_mono = monospaced[n_latin] = FALSE;
		    }
		    n_latin++;
		}
	    }
	    pango_font_description_free(desc);
	}

	g_free(fontname);

    } else if (filter != GTK_FONT_HACK_NONE) { /* refer to the cached information */
	gint i;

	for (i=0; i<n_latin; i++) {
	    if (strcmp(familyname, latin_families[i]) == 0) {
		is_latin = TRUE;
		is_mono = monospaced[i];
		break;
	    }
	    if (latin_families[i][0] > familyname[0]) break;
	}
    }

    if (filter == GTK_FONT_HACK_LATIN) return is_latin;
    else if (filter == GTK_FONT_HACK_LATIN_MONO) return (is_latin && is_mono);
    else if (filter == GTK_FONT_HACK_NONE) return TRUE;
    return FALSE; /* not reachable */
}

