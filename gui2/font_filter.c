static gboolean font_is_latin_monospaced (PangoFontDescription *desc)
{
    PangoLayout *layout;
    PangoContext *context;
    GtkWidget *wdg;
    int x1, x2;

    wdg = gtk_label_new(NULL);

    context = gtk_widget_get_pango_context(wdg);
    if (context == NULL) {
	gtk_widget_destroy(wdg);
	return FALSE;
    }

    pango_context_set_font_description(context, desc);

    layout = pango_layout_new(context);
    pango_layout_set_text(layout, "i", 1);
    pango_layout_get_pixel_size(layout, &x1, NULL);
    pango_layout_set_text(layout, "W", 1);
    pango_layout_get_pixel_size(layout, &x2, NULL);

    g_object_unref(G_OBJECT(layout));
    g_object_unref(G_OBJECT(context));
    gtk_widget_destroy(wdg);

    return (x1 == x2);
}

static gboolean font_is_latin_text_font (PangoFontDescription *desc,
					 PangoContext *context) 
{
    PangoCoverage *coverage;
    PangoLanguage *lang;
    PangoFont *pfont;
    gboolean ok = FALSE;

    pfont = pango_context_load_font(context, desc);
    if (pfont == NULL) return FALSE;

    lang = pango_language_from_string("en_US");
    if (lang == NULL) return FALSE;

    coverage = pango_font_get_coverage(pfont, lang);
    if (coverage == NULL) return FALSE;

    if (pango_coverage_get(coverage, 'A') == PANGO_COVERAGE_EXACT) {
	ok = TRUE;
    } 

    /* what needs to be unref'd now? */
    /* g_object_unref(G_OBJECT(lang)); */

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
   of seeing if it contains the letter 'A' in US English.  Given the
   latin text characteristic, we can then see if the font is monospaced
   by checking whether or not the letters 'i' and 'W' come out the
   same width. */

#define SHOW_PROGRESS 1

static gboolean validate_font_family (const gchar *familyname, 
				      PangoContext *context,
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
    static int (*show_progress) (long, long, int) = NULL;
    static int show = 1, n_done = 0;
#endif

    if (!cache_built) { /* we haven't set up the cache yet */
	gchar *fontname;
	PangoFontDescription *desc;

#ifdef FONT_FILTER_DEBUG
	fprintf(dbg, "Checking font family '%s'\n", familyname);
	fflush(dbg);
#endif

#ifdef SHOW_PROGRESS
	if (show && handle == NULL) {
	    int err = gui_open_plugin("progress_bar-2", &handle);

	    if (err) show = 0;
	    else {
		show_progress = get_plugin_function("show_progress", handle);
		if (show_progress == NULL) show = 0;
		else (*show_progress)(0L, n_families, SP_FONT_INIT);
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
#endif

	if (weird_font(familyname)) return FALSE;

	fontname = g_strdup_printf("%s 10", familyname);
	desc = pango_font_description_from_string(fontname);
	g_free(fontname);

	if (desc != NULL) {
	    if (font_is_latin_text_font(desc, context)) {
		/* extend the cache */
		latin_families = g_realloc(latin_families, (n_latin + 1) * sizeof *latin_families);
		monospaced = g_realloc(monospaced, (n_latin + 1) * sizeof *monospaced);

		latin_families[n_latin] = g_strdup(familyname);
		is_latin = TRUE;
		if (font_is_latin_monospaced(desc)) {
		    is_mono = monospaced[n_latin] = TRUE;
		} else {
		    is_mono = monospaced[n_latin] = FALSE;
		}
		n_latin++;
	    }
	    pango_font_description_free(desc);
	}

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

