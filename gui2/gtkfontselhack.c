/* GTK - The GIMP Toolkit
 * Copyright (C) 1995-1997 Peter Mattis, Spencer Kimball and Josh MacDonald
 *
 * Massively updated for Pango by Owen Taylor, May 2000
 * GtkFontSelection widget for Gtk+, by Damon Chaplin, May 1998.
 * Based on the GnomeFontSelector widget, by Elliot Lee, but major changes.
 * The GnomeFontSelector was derived from app/text_tool.c in the GIMP.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

/*
 * Modified by the GTK+ Team and others 1997-2000.  See the AUTHORS
 * file for a list of people on the GTK+ Team.  See the ChangeLog
 * files for a list of changes.  These files are distributed with
 * GTK+ at ftp://ftp.gtk.org/pub/gtk/. 
 */

/*
 * Modified to allow filtering of fonts for gretl, Allin Cottrell, 2002
 */

#include "gretl.h" /* this must provide standard and GTK headers and
		      deal with gettext */
#include "dlgutils.h"
#include "gtkfontselhack.h"
#include "fontfilter.h"

#define GTK_TYPE_FNTHACK              (gtk_fontsel_hack_get_type ())
#define GTK_FNTHACK(obj)              (G_TYPE_CHECK_INSTANCE_CAST ((obj), GTK_TYPE_FNTHACK, GtkFontselHack))
#define GTK_FNTHACK_CLASS(klass)      (G_TYPE_CHECK_CLASS_CAST ((klass), GTK_TYPE_FNTHACK, GtkFontselHackClass))
#define GTK_IS_FNTHACK(obj)           (G_TYPE_CHECK_INSTANCE_TYPE ((obj), GTK_TYPE_FNTHACK))
#define GTK_IS_FNTHACK_CLASS(klass)   (G_TYPE_CHECK_CLASS_TYPE ((klass), GTK_TYPE_FNTHACK))
#define GTK_FNTHACK_GET_CLASS(obj)    (G_TYPE_INSTANCE_GET_CLASS ((obj), GTK_TYPE_FNTHACK, GtkFontselHackClass))


#define GTK_TYPE_FNTHACK_DIALOG              (gtk_fontsel_hack_dialog_get_type ())
#define GTK_FNTHACK_DIALOG(obj)              (G_TYPE_CHECK_INSTANCE_CAST ((obj), GTK_TYPE_FNTHACK_DIALOG, GtkFontselHackDialog))
#define GTK_FNTHACK_DIALOG_CLASS(klass)      (G_TYPE_CHECK_CLASS_CAST ((klass), GTK_TYPE_FNTHACK_DIALOG, GtkFontselHackDialogClass))
#define GTK_IS_FNTHACK_DIALOG(obj)           (G_TYPE_CHECK_INSTANCE_TYPE ((obj), GTK_TYPE_FNTHACK_DIALOG))
#define GTK_IS_FNTHACK_DIALOG_CLASS(klass)   (G_TYPE_CHECK_CLASS_TYPE ((klass), GTK_TYPE_FNTHACK_DIALOG))
#define GTK_FNTHACK_DIALOG_GET_CLASS(obj)    (G_TYPE_INSTANCE_GET_CLASS ((obj), GTK_TYPE_FNTHACK_DIALOG, GtkFontselHackDialogClass))

#define GTK_TYPE_FONT_FILTER (gtk_font_filter_get_type())

struct _GtkFontselHack
{
    GtkVBox parent_instance;
  
    GtkWidget *font_entry;
    GtkWidget *family_list;
    GtkWidget *font_style_entry;
    GtkWidget *face_list;
    GtkWidget *size_entry;
    GtkWidget *size_list;
    GtkWidget *pixels_button;
    GtkWidget *points_button;
    GtkWidget *filter_button;
    GtkWidget *preview_entry;

    PangoFontFamily *family;	/* Current family */
    PangoFontFace *face;        /* Current face */
  
    gint size;
    FontFilterType filter;
};

struct _GtkFontselHackClass
{
    GtkVBoxClass parent_class;

    /* Padding for future expansion */
    void (*_gtk_reserved1) (void);
    void (*_gtk_reserved2) (void);
    void (*_gtk_reserved3) (void);
    void (*_gtk_reserved4) (void);
};

struct _GtkFontselHackDialog
{
    GtkDialog parent_instance;
    GtkWidget *fontsel;
    GtkWidget *main_vbox;
    GtkWidget *action_area;
    GtkWidget *ok_button;
    /* The 'Apply' button is not shown by default but you can show/hide it. */
    GtkWidget *apply_button;
    GtkWidget *cancel_button;
    /* If the user changes the width of the dialog, we turn auto-shrink off. */
    gint dialog_width;
    gboolean auto_resize;
};

struct _GtkFontselHackDialogClass
{
    GtkDialogClass parent_class;

    /* Padding for future expansion */
    void (*_gtk_reserved1) (void);
    void (*_gtk_reserved2) (void);
    void (*_gtk_reserved3) (void);
    void (*_gtk_reserved4) (void);
};

/* This is the default text shown in the preview entry, though the user
   can set it. Remember that some fonts only have capital letters. */
#define PREVIEW_TEXT N_("abcdefghijk ABCDEFGHIJK")

/* This is the initial and maximum height of the preview entry (it expands
   when large font sizes are selected). Initial height is also the minimum. */
#define INITIAL_PREVIEW_HEIGHT 44
#define MAX_PREVIEW_HEIGHT 300

/* These are the sizes of the font, style & size lists. */
#define FONT_LIST_HEIGHT	136
#define FONT_LIST_WIDTH		190
#define FONT_STYLE_LIST_WIDTH	170
#define FONT_SIZE_LIST_WIDTH	60

/* These are what we use as the standard font sizes, for the size list.
 */
static const guint16 font_sizes[] = {
    6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 20, 22, 24, 26, 28,
    32, 36, 40, 48, 56, 64, 72
};

enum {
    PROP_0,
    PROP_FONT_NAME,
    PROP_PREVIEW_TEXT,
    PROP_FILTER
};

enum {
    FAMILY_COLUMN,
    FAMILY_NAME_COLUMN
};

enum {
    FACE_COLUMN,
    FACE_NAME_COLUMN
};

enum {
    SIZE_COLUMN
};

static void    gtk_fontsel_hack_class_init	   (GtkFontselHackClass *klass);
static void    gtk_fontsel_hack_set_property       (GObject         *object,
						    guint            prop_id,
						    const GValue    *value,
						    GParamSpec      *pspec);
static void    gtk_fontsel_hack_get_property       (GObject         *object,
						    guint            prop_id,
						    GValue          *value,
						    GParamSpec      *pspec);
static void    gtk_fontsel_hack_init		     (GtkFontselHack      *fontsel);
static void    gtk_fontsel_hack_finalize	     (GObject               *object);

/* These are the callbacks & related functions. */
static void     gtk_fontsel_hack_select_font           (GtkTreeSelection *selection,
							gpointer          data);
static void     gtk_fontsel_hack_show_available_fonts  (GtkFontselHack *fs);

static void     gtk_fontsel_hack_show_available_styles (GtkFontselHack *fs);
static void     gtk_fontsel_hack_select_best_style     (GtkFontselHack *fs,
							gboolean          use_first);
static void     gtk_fontsel_hack_select_style          (GtkTreeSelection *selection,
							gpointer          data);

static void     gtk_fontsel_hack_select_best_size      (GtkFontselHack *fs);
static void     gtk_fontsel_hack_show_available_sizes  (GtkFontselHack *fs,
							gboolean          first_time);
static void     gtk_fontsel_hack_size_activate         (GtkWidget        *w,
							gpointer          data);
static gboolean gtk_fontsel_hack_size_focus_out        (GtkWidget        *w,
							GdkEventFocus    *event,
							gpointer          data);
static void     gtk_fontsel_hack_select_size           (GtkTreeSelection *selection,
							gpointer          data);

static void     gtk_fontsel_hack_scroll_on_map         (GtkWidget        *w,
							gpointer          data);

static void     gtk_fontsel_hack_screen_changed        (GtkWidget *widget,
							GdkScreen *previous_screen);

static void     gtk_fontsel_hack_preview_changed       (GtkWidget        *entry,
							GtkFontselHack *fontsel);

/* Misc. utility functions. */
static void gtk_fontsel_hack_load_font          (GtkFontselHack *fs);
static void gtk_fontsel_hack_update_preview     (GtkFontselHack *fs);

/* FontselHackDialog */
static void gtk_fontsel_hack_dialog_class_init  (GtkFontselHackDialogClass *klass);
static void gtk_fontsel_hack_dialog_init        (GtkFontselHackDialog *fontseldiag);

static gint gtk_fontsel_hack_dialog_on_configure (GtkWidget      *widget,
						  GdkEventConfigure *event,
						  GtkFontselHackDialog *fsd);

GType gtk_font_filter_get_type (void)
{
    static GType etype = 0;

    if (etype == 0) {
	static const GEnumValue values[] = {
	    { FONT_HACK_NONE, "FONT_HACK_NONE", "use no font filter" },
	    { FONT_HACK_LATIN, "FONT_HACK_LATIN", "latin text fonts" },
	    { FONT_HACK_LATIN_MONO, "FONT_HACK_LATIN_MONO", "monospaced latin text fonts" },
	    { 0, NULL, NULL }
	};

	etype = g_enum_register_static("GtkFontFilter", values);
    }

    return etype;
}

G_DEFINE_TYPE (GtkFontselHack, gtk_fontsel_hack, GTK_TYPE_VBOX)

static void
gtk_fontsel_hack_class_init (GtkFontselHackClass *klass)
{
    GObjectClass *gobject_class = G_OBJECT_CLASS(klass);
    GtkWidgetClass *widget_class = GTK_WIDGET_CLASS(klass);
  
    gobject_class->set_property = gtk_fontsel_hack_set_property;
    gobject_class->get_property = gtk_fontsel_hack_get_property;

    widget_class->screen_changed = gtk_fontsel_hack_screen_changed;
   
    g_object_class_install_property(gobject_class,
				    PROP_FONT_NAME,
				    g_param_spec_string("font-name",
							_("Font name"),
							_("The X string that represents this font"),
							NULL,
							G_PARAM_READWRITE));
    g_object_class_install_property(gobject_class,
				    PROP_PREVIEW_TEXT,
				    g_param_spec_string("preview-text",
							_("Preview text"),
							_("The text to display in order to demonstrate the selected font"),
							PREVIEW_TEXT,
							G_PARAM_READWRITE));
    g_object_class_install_property(gobject_class,
				    PROP_FILTER,
				    g_param_spec_enum("filter",
						      _("Filter"),
						      _("The filter for acceptable fonts."),
						      GTK_TYPE_FONT_FILTER,
						      FONT_HACK_NONE,
						      G_PARAM_READWRITE));

    gobject_class->finalize = gtk_fontsel_hack_finalize;
}

static void 
gtk_fontsel_hack_set_property (GObject *object,
			       guint prop_id,
			       const GValue *value,
			       GParamSpec *pspec)
{
    GtkFontselHack *fontsel;

    fontsel = GTK_FNTHACK(object);

    switch (prop_id) {
    case PROP_FONT_NAME:
	gtk_fontsel_hack_set_font_name(fontsel, g_value_get_string(value));
	break;
    case PROP_PREVIEW_TEXT:
	gtk_fontsel_hack_set_preview_text(fontsel, g_value_get_string(value));
	break;
    default:
	G_OBJECT_WARN_INVALID_PROPERTY_ID(object, prop_id, pspec);
	break;
    }
}

static void gtk_fontsel_hack_get_property (GObject *object,
					   guint prop_id,
					   GValue *value,
					   GParamSpec *pspec)
{
    GtkFontselHack *fontsel;

    fontsel = GTK_FNTHACK(object);

    switch (prop_id) {
    case PROP_FONT_NAME:
	g_value_set_string(value, gtk_fontsel_hack_get_font_name(fontsel));
	break;
    case PROP_PREVIEW_TEXT:
	g_value_set_string(value, gtk_fontsel_hack_get_preview_text(fontsel));
	break;
    case PROP_FILTER:
	g_value_set_int(value, gtk_fontsel_hack_get_filter(fontsel));
	break;
    default:
	G_OBJECT_WARN_INVALID_PROPERTY_ID(object, prop_id, pspec);
	break;
    }
}

static void
gtk_fontsel_hack_init (GtkFontselHack *fontsel)
{
    GtkWidget *scrolled_win;
    GtkWidget *text_frame;
    GtkWidget *text_box;
    GtkWidget *table, *label;
    GtkWidget *font_label, *style_label;
    GtkListStore *model;
    GtkTreeViewColumn *column;
    GList *focus_chain = NULL;

    gtk_widget_push_composite_child();

    fontsel->size = 12 * PANGO_SCALE;
    fontsel->filter = FONT_HACK_NONE;
  
    /* Create the table of font, style & size */
    table = gtk_table_new(3, 3, FALSE);
    gtk_widget_show(table);
    gtk_table_set_col_spacings(GTK_TABLE (table), 8);
    gtk_box_pack_start(GTK_BOX(fontsel), table, TRUE, TRUE, 0);

    fontsel->size_entry = gtk_entry_new();
    gtk_widget_set_size_request(fontsel->size_entry, 20, -1);
    gtk_widget_show(fontsel->size_entry);
    gtk_table_attach(GTK_TABLE(table), fontsel->size_entry, 2, 3, 1, 2,
		     GTK_FILL, 0, 0, 0);
    g_signal_connect(G_OBJECT(fontsel->size_entry), "activate",
		     (GCallback) gtk_fontsel_hack_size_activate,
		     fontsel);
    g_signal_connect_after(G_OBJECT(fontsel->size_entry), "focus-out-event",
			   (GCallback) gtk_fontsel_hack_size_focus_out,
			   fontsel);
  
    font_label = gtk_label_new_with_mnemonic(_("_Family:"));
    gtk_misc_set_alignment(GTK_MISC(font_label), 0.0, 0.5);
    gtk_widget_show(font_label);
    gtk_table_attach(GTK_TABLE(table), font_label, 0, 1, 0, 1,
		     GTK_FILL, 0, 0, 0);  

    style_label = gtk_label_new_with_mnemonic(_("_Style:"));
    gtk_misc_set_alignment(GTK_MISC(style_label), 0.0, 0.5);
    gtk_widget_show(style_label);
    gtk_table_attach(GTK_TABLE(table), style_label, 1, 2, 0, 1,
		     GTK_FILL, 0, 0, 0);
  
    label = gtk_label_new_with_mnemonic(_("Si_ze:"));
    gtk_label_set_mnemonic_widget(GTK_LABEL(label),
				  fontsel->size_entry);
    gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.5);
    gtk_widget_show(label);
    gtk_table_attach(GTK_TABLE(table), label, 2, 3, 0, 1,
		     GTK_FILL, 0, 0, 0);
  
  
    /* Create the lists  */

    model = gtk_list_store_new(2,
			       G_TYPE_OBJECT,  /* FAMILY_COLUMN */
			       G_TYPE_STRING); /* FAMILY_NAME_COLUMN */
    fontsel->family_list = gtk_tree_view_new_with_model(GTK_TREE_MODEL(model));
    g_object_unref(model);

    column = gtk_tree_view_column_new_with_attributes("Family",
						      gtk_cell_renderer_text_new(),
						      "text", FAMILY_NAME_COLUMN,
						      NULL);
    gtk_tree_view_column_set_sizing(column, GTK_TREE_VIEW_COLUMN_AUTOSIZE);
    gtk_tree_view_append_column(GTK_TREE_VIEW(fontsel->family_list), column);

    gtk_tree_view_set_headers_visible(GTK_TREE_VIEW(fontsel->family_list), FALSE);
    gtk_tree_selection_set_mode(gtk_tree_view_get_selection(GTK_TREE_VIEW(fontsel->family_list)),
				GTK_SELECTION_BROWSE);
  
    gtk_label_set_mnemonic_widget(GTK_LABEL(font_label), fontsel->family_list);

    scrolled_win = gtk_scrolled_window_new(NULL, NULL);
    gtk_scrolled_window_set_shadow_type(GTK_SCROLLED_WINDOW(scrolled_win), GTK_SHADOW_IN);
    gtk_widget_set_size_request(scrolled_win, FONT_LIST_WIDTH, FONT_LIST_HEIGHT);
    gtk_container_add(GTK_CONTAINER(scrolled_win), fontsel->family_list);
    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scrolled_win),
				   GTK_POLICY_AUTOMATIC, GTK_POLICY_ALWAYS);
    gtk_widget_show(fontsel->family_list);
    gtk_widget_show(scrolled_win);

    gtk_table_attach(GTK_TABLE(table), scrolled_win, 0, 1, 1, 3,
		     GTK_EXPAND | GTK_FILL,
		     GTK_EXPAND | GTK_FILL, 0, 0);
    focus_chain = g_list_append(focus_chain, scrolled_win);
  
    model = gtk_list_store_new(2,
			       G_TYPE_OBJECT,  /* FACE_COLUMN */
			       G_TYPE_STRING); /* FACE_NAME_COLUMN */
    fontsel->face_list = gtk_tree_view_new_with_model(GTK_TREE_MODEL(model));
    g_object_unref(model);

    gtk_label_set_mnemonic_widget(GTK_LABEL(style_label), fontsel->face_list);

    column = gtk_tree_view_column_new_with_attributes("Face",
						      gtk_cell_renderer_text_new(),
						      "text", FACE_NAME_COLUMN,
						      NULL);
    gtk_tree_view_column_set_sizing(column, GTK_TREE_VIEW_COLUMN_AUTOSIZE);
    gtk_tree_view_append_column(GTK_TREE_VIEW(fontsel->face_list), column);

    gtk_tree_view_set_headers_visible(GTK_TREE_VIEW(fontsel->face_list), FALSE);
    gtk_tree_selection_set_mode(gtk_tree_view_get_selection(GTK_TREE_VIEW(fontsel->face_list)),
				GTK_SELECTION_BROWSE);
  
    scrolled_win = gtk_scrolled_window_new(NULL, NULL);
    gtk_scrolled_window_set_shadow_type(GTK_SCROLLED_WINDOW(scrolled_win), GTK_SHADOW_IN);
    gtk_widget_set_size_request(scrolled_win, FONT_STYLE_LIST_WIDTH, FONT_LIST_HEIGHT);
    gtk_container_add(GTK_CONTAINER(scrolled_win), fontsel->face_list);
    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scrolled_win),
				   GTK_POLICY_AUTOMATIC, GTK_POLICY_ALWAYS);
    gtk_widget_show(fontsel->face_list);
    gtk_widget_show(scrolled_win);
    gtk_table_attach(GTK_TABLE(table), scrolled_win, 1, 2, 1, 3,
		     GTK_EXPAND | GTK_FILL,
		     GTK_EXPAND | GTK_FILL, 0, 0);
    focus_chain = g_list_append(focus_chain, scrolled_win);
  
    focus_chain = g_list_append(focus_chain, fontsel->size_entry);

    model = gtk_list_store_new(1, G_TYPE_INT);
    fontsel->size_list = gtk_tree_view_new_with_model(GTK_TREE_MODEL(model));
    g_object_unref(model);

    column = gtk_tree_view_column_new_with_attributes("Size",
						      gtk_cell_renderer_text_new(),
						      "text", SIZE_COLUMN,
						      NULL);
    gtk_tree_view_column_set_sizing(column, GTK_TREE_VIEW_COLUMN_AUTOSIZE);
    gtk_tree_view_append_column(GTK_TREE_VIEW(fontsel->size_list), column);

    gtk_tree_view_set_headers_visible(GTK_TREE_VIEW(fontsel->size_list), FALSE);
    gtk_tree_selection_set_mode(gtk_tree_view_get_selection(GTK_TREE_VIEW(fontsel->size_list)),
				GTK_SELECTION_BROWSE);
  
    scrolled_win = gtk_scrolled_window_new(NULL, NULL);
    gtk_scrolled_window_set_shadow_type(GTK_SCROLLED_WINDOW(scrolled_win), GTK_SHADOW_IN);
    gtk_container_add(GTK_CONTAINER(scrolled_win), fontsel->size_list);
    gtk_widget_set_size_request(scrolled_win, -1, FONT_LIST_HEIGHT);
    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scrolled_win),
				   GTK_POLICY_NEVER, GTK_POLICY_ALWAYS);
    gtk_widget_show(fontsel->size_list);
    gtk_widget_show(scrolled_win);
    gtk_table_attach(GTK_TABLE(table), scrolled_win, 2, 3, 2, 3,
		     GTK_FILL, GTK_EXPAND | GTK_FILL, 0, 0);
    focus_chain = g_list_append(focus_chain, scrolled_win);

    gtk_container_set_focus_chain(GTK_CONTAINER(table), focus_chain);
    g_list_free(focus_chain);
  
    /* Insert the fonts. */
    gtk_fontsel_hack_show_available_fonts(fontsel);
  
    g_signal_connect(gtk_tree_view_get_selection(GTK_TREE_VIEW(fontsel->family_list)), "changed",
		     G_CALLBACK(gtk_fontsel_hack_select_font), fontsel);

    g_signal_connect_after(G_OBJECT(fontsel->family_list), "map",
			   G_CALLBACK(gtk_fontsel_hack_scroll_on_map),
			   fontsel);
  
    gtk_fontsel_hack_show_available_styles(fontsel);
  
    g_signal_connect(gtk_tree_view_get_selection(GTK_TREE_VIEW(fontsel->face_list)), "changed",
		     G_CALLBACK(gtk_fontsel_hack_select_style), fontsel);

    gtk_fontsel_hack_show_available_sizes(fontsel, TRUE);
  
    g_signal_connect(gtk_tree_view_get_selection(GTK_TREE_VIEW(fontsel->size_list)), "changed",
		     G_CALLBACK(gtk_fontsel_hack_select_size), fontsel);

    /* create the text entry widget */
    label = gtk_label_new_with_mnemonic(_("_Preview:"));
    gtk_widget_show(label);
  
    text_frame = gtk_frame_new(NULL);
    gtk_frame_set_label_widget(GTK_FRAME(text_frame), label);
  
    gtk_widget_show(text_frame);
    gtk_frame_set_shadow_type(GTK_FRAME(text_frame), GTK_SHADOW_ETCHED_IN);
    gtk_box_pack_start(GTK_BOX(fontsel), text_frame,
		       FALSE, TRUE, 0);
  
    /* This is just used to get a 4-pixel space around the preview entry. */
    text_box = gtk_hbox_new(FALSE, 0);
    gtk_widget_show(text_box);
    gtk_container_add(GTK_CONTAINER(text_frame), text_box);
    gtk_container_set_border_width(GTK_CONTAINER(text_box), 4);
  
    fontsel->preview_entry = gtk_entry_new();
    gtk_label_set_mnemonic_widget(GTK_LABEL(label), fontsel->preview_entry);
  
    gtk_widget_show(fontsel->preview_entry);
    g_signal_connect(G_OBJECT(fontsel->preview_entry), "changed",
		     (GCallback) gtk_fontsel_hack_preview_changed,
		     fontsel);
    gtk_widget_set_size_request(fontsel->preview_entry, -1, INITIAL_PREVIEW_HEIGHT);
    gtk_box_pack_start(GTK_BOX(text_box), fontsel->preview_entry,
		       TRUE, TRUE, 0);

    gtk_fontsel_hack_update_preview(fontsel);

    gtk_widget_pop_composite_child();
}

static void
gtk_fontsel_hack_display_fonts (GtkFontselHack *fontsel)
{
    /* Insert the (possibly revised) list of fonts. */
    gtk_fontsel_hack_show_available_fonts(fontsel);
    gtk_fontsel_hack_show_available_styles(fontsel);
    gtk_fontsel_hack_show_available_sizes(fontsel, TRUE);
}

GtkWidget *
gtk_fontsel_hack_new (void)
{
    GtkFontselHack *fontsel;
  
    fontsel = g_object_new(GTK_TYPE_FNTHACK, NULL);
  
    return GTK_WIDGET(fontsel);
}

static void
gtk_fontsel_hack_finalize (GObject *object)
{
    g_return_if_fail(GTK_IS_FNTHACK(object));
  
    (* G_OBJECT_CLASS (gtk_fontsel_hack_parent_class)->finalize) (object);
}

static void
gtk_fontsel_hack_screen_changed (GtkWidget *widget,
				 GdkScreen *previous_screen)
{
    GtkFontselHack *fontsel = GTK_FNTHACK(widget);

    if (gtk_widget_has_screen(GTK_WIDGET(fontsel))) {
	gtk_fontsel_hack_show_available_fonts(fontsel);
	gtk_fontsel_hack_show_available_sizes(fontsel, TRUE);
	gtk_fontsel_hack_show_available_styles(fontsel);
    }
}

static void
gtk_fontsel_hack_preview_changed (GtkWidget *entry,
				  GtkFontselHack *fontsel)
{
    g_object_notify(G_OBJECT(fontsel), "preview_text");
}

static void
scroll_to_selection (GtkTreeView *tree_view)
{
    GtkTreeSelection *selection = gtk_tree_view_get_selection(tree_view);
    GtkTreeModel *model;
    GtkTreeIter iter;

    if (gtk_tree_selection_get_selected(selection, &model, &iter)) {
	GtkTreePath *path = gtk_tree_model_get_path(model, &iter);

	gtk_tree_view_scroll_to_cell(tree_view, path, NULL, TRUE, 0.5, 0.5);
	gtk_tree_path_free(path);
    }
}

static void
set_cursor_to_iter (GtkTreeView *view,
		    GtkTreeIter *iter)
{
    GtkTreeModel *model = gtk_tree_view_get_model(view);
    GtkTreePath *path = gtk_tree_model_get_path(model, iter);

    if (path != NULL) {
	gtk_tree_view_set_cursor(view, path, 0, FALSE);
	gtk_tree_path_free(path);
    }
}

/* This is called when the list is mapped. Here we scroll to the current
   font if necessary. */
static void
gtk_fontsel_hack_scroll_on_map (GtkWidget *widget,
				gpointer data)
{
    GtkFontselHack *fontsel;
  
    fontsel = GTK_FNTHACK (data);
  
    /* Try to scroll the font family list to the selected item */
    scroll_to_selection(GTK_TREE_VIEW(fontsel->family_list));
      
    /* Try to scroll the font family list to the selected item */
    scroll_to_selection(GTK_TREE_VIEW(fontsel->face_list));
      
    /* Try to scroll the font family list to the selected item */
    scroll_to_selection(GTK_TREE_VIEW(fontsel->size_list));
}

/* This is called when a family is selected in the list. */
static void
gtk_fontsel_hack_select_font (GtkTreeSelection *selection,
			      gpointer data)
{
    GtkFontselHack *fontsel;
    GtkTreeModel *model;
    GtkTreeIter iter;
  
    fontsel = GTK_FNTHACK(data);

    if (gtk_tree_selection_get_selected(selection, &model, &iter)) {
	PangoFontFamily *family;
      
	gtk_tree_model_get(model, &iter, FAMILY_COLUMN, &family, -1);
	if (fontsel->family != family) {
	    fontsel->family = family;
	    gtk_fontsel_hack_show_available_styles(fontsel);
	    gtk_fontsel_hack_select_best_style(fontsel, TRUE);
	}

	g_object_unref(family);
    }
}

static int
cmp_families (const void *a, const void *b)
{
    const char *a_name = pango_font_family_get_name (*(PangoFontFamily **)a);
    const char *b_name = pango_font_family_get_name (*(PangoFontFamily **)b);
  
    return g_utf8_collate(a_name, b_name);
}

static void
gtk_fontsel_hack_show_available_fonts (GtkFontselHack *fontsel)
{
    GtkTreeModel *model;
    GtkListStore *store;
    PangoContext *context; 
    PangoFontFamily **families;
    const char *famname;
    gint nf, i, got_ok;
    GtkTreeIter iter, match_iter;
    gint err = 0;

    model = gtk_tree_view_get_model(GTK_TREE_VIEW(fontsel->family_list));
    store = GTK_LIST_STORE(model);
    gtk_list_store_clear(store);

    context = gtk_widget_get_pango_context(GTK_WIDGET(fontsel));
    pango_context_list_families(context, &families, &nf);
    qsort(families, nf, sizeof *families, cmp_families);

#if FDEBUG
    fprintf(stderr, "Font selector: got %d families\n", nf);
#endif

    fontsel->family = NULL;
    got_ok = 0;

    gtk_tree_model_get_iter_first(model, &iter);

    for (i=0; i<nf && !err; i++) {
	famname = pango_font_family_get_name(families[i]);
	if (famname == NULL) {
	    /* bullet-proofing */
	    continue;
	}

#if FDEBUG
	fprintf(stderr, "Examining font family '%s'\n", famname);
#endif
	if (!validate_font_family(families[i], famname, i, nf, fontsel->filter, &err)) { 
	    continue;
	}

	gtk_list_store_append(store, &iter);
	gtk_list_store_set(store, &iter,
			   FAMILY_COLUMN, families[i],
			   FAMILY_NAME_COLUMN, famname,
			   -1);

	if (!got_ok) {
	    if (fontsel->filter != FONT_HACK_NONE ||
		i == 0 || !g_ascii_strcasecmp(famname, "sans")) {
		got_ok = 1;
		fontsel->family = families[i];
		match_iter = iter;
	    }
	}
    }

    if (fontsel->family != NULL) {
	set_cursor_to_iter(GTK_TREE_VIEW(fontsel->family_list), &match_iter);
    }

    g_free(families);
}

static int
compare_font_descriptions (const PangoFontDescription *a, 
			   const PangoFontDescription *b)
{
    int val = strcmp(pango_font_description_get_family(a), 
		     pango_font_description_get_family(b));

    if (val != 0) {
	return val;
    }

    if (pango_font_description_get_weight(a) != 
	pango_font_description_get_weight(b))
	return pango_font_description_get_weight(a) - 
	    pango_font_description_get_weight(b);

    if (pango_font_description_get_style(a) != 
	pango_font_description_get_style(b))
	return pango_font_description_get_style(a) - 
	    pango_font_description_get_style(b);
  
    if (pango_font_description_get_stretch(a) != 
	pango_font_description_get_stretch(b))
	return pango_font_description_get_stretch(a) - 
	    pango_font_description_get_stretch(b);

    if (pango_font_description_get_variant(a) != 
	pango_font_description_get_variant(b))
	return pango_font_description_get_variant(a) - 
	    pango_font_description_get_variant(b);

    return 0;
}

static int
faces_sort_func (const void *a, const void *b)
{
    PangoFontDescription *desc_a = pango_font_face_describe(*(PangoFontFace **)a);
    PangoFontDescription *desc_b = pango_font_face_describe(*(PangoFontFace **)b);
  
    int ord = compare_font_descriptions(desc_a, desc_b);

    pango_font_description_free(desc_a);
    pango_font_description_free(desc_b);

    return ord;
}

static gboolean
font_description_style_equal (const PangoFontDescription *a,
			      const PangoFontDescription *b)
{
    return (pango_font_description_get_weight(a) == 
	    pango_font_description_get_weight(b) &&
	    pango_font_description_get_style(a) == 
	    pango_font_description_get_style(b) &&
	    pango_font_description_get_stretch(a) == 
	    pango_font_description_get_stretch(b) &&
	    pango_font_description_get_variant(a) == 
	    pango_font_description_get_variant(b));
}

/* This fills the font style list with all the possible style combinations
   for the current font family. */

static void
gtk_fontsel_hack_show_available_styles (GtkFontselHack *fontsel)
{
    gint n_faces, i;
    PangoFontFace **faces;
    PangoFontDescription *old_desc;
    GtkListStore *model;
    GtkTreeIter iter, match_row;
    PangoFontFace *match_face = NULL;
  
    model = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(fontsel->face_list)));
  
    if (fontsel->face) {
	old_desc = pango_font_face_describe(fontsel->face);
    } else {
	old_desc= NULL;
    }

    pango_font_family_list_faces(fontsel->family, &faces, &n_faces);
    qsort(faces, n_faces, sizeof *faces, faces_sort_func);

    gtk_list_store_clear(model);
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(model), &iter);

    for (i=0; i<n_faces; i++) {
	const gchar *str = pango_font_face_get_face_name(faces[i]);

	gtk_list_store_append(model, &iter);
	gtk_list_store_set(model, &iter,
			   FACE_COLUMN, faces[i],
			   FACE_NAME_COLUMN, str,
			   -1);

	if (i == 0) {
	    match_row = iter;
	    match_face = faces[i];
	} else if (old_desc) {
	    PangoFontDescription *tmp_desc = pango_font_face_describe(faces[i]);
	  
	    if (font_description_style_equal(tmp_desc, old_desc)) {
		match_row = iter;
		match_face = faces[i];
	    }
	    pango_font_description_free(tmp_desc);
	}
    }

    if (old_desc != NULL) {
	pango_font_description_free(old_desc);
    }

    fontsel->face = match_face;
    if (match_face) {
	set_cursor_to_iter(GTK_TREE_VIEW(fontsel->face_list), &match_row);
    }

    g_free(faces);
}

/* This selects a style when the user selects a font. It just uses the first
   available style at present. I was thinking of trying to maintain the
   selected style, e.g. bold italic, when the user selects different fonts.
   However, the interface is so easy to use now I'm not sure it's worth it.
   Note: This will load a font. */

static void
gtk_fontsel_hack_select_best_style (GtkFontselHack *fontsel,
				    gboolean use_first)
{
    GtkTreeIter iter;
    GtkTreeModel *model;

    model = gtk_tree_view_get_model (GTK_TREE_VIEW (fontsel->face_list));

    if (gtk_tree_model_get_iter_first (model, &iter)) {
	set_cursor_to_iter(GTK_TREE_VIEW(fontsel->face_list), &iter);
	scroll_to_selection(GTK_TREE_VIEW(fontsel->face_list));
    }

    gtk_fontsel_hack_show_available_sizes(fontsel, FALSE);
    gtk_fontsel_hack_select_best_size(fontsel);
}


/* This is called when a style is selected in the list. */

static void
gtk_fontsel_hack_select_style (GtkTreeSelection *selection,
			       gpointer          data)
{
    GtkFontselHack *fontsel = GTK_FNTHACK(data);
    GtkTreeModel *model;
    GtkTreeIter iter;
  
    if (gtk_tree_selection_get_selected(selection, &model, &iter)) {
	PangoFontFace *face;
      
	gtk_tree_model_get(model, &iter, FACE_COLUMN, &face, -1);
	fontsel->face = face;

	g_object_unref(face);
    }

    gtk_fontsel_hack_show_available_sizes(fontsel, FALSE);
    gtk_fontsel_hack_select_best_size(fontsel);
}

static void
gtk_fontsel_hack_show_available_sizes (GtkFontselHack *fontsel,
				       gboolean first_time)
{
    size_t i;
    GtkListStore *model;
    gchar buffer[128];
    gchar *p;
      
    model = GTK_LIST_STORE(gtk_tree_view_get_model (GTK_TREE_VIEW (fontsel->size_list)));

    /* Insert the standard font sizes */
    if (first_time) {
	GtkTreeIter iter;

	gtk_list_store_clear (model);
	gtk_tree_model_get_iter_first (GTK_TREE_MODEL (model), &iter);

	for (i = 0; i < G_N_ELEMENTS (font_sizes); i++) {
	    gtk_list_store_append (model, &iter);
	    gtk_list_store_set (model, &iter, SIZE_COLUMN, font_sizes[i], -1);
	  
	    if (font_sizes[i] * PANGO_SCALE == fontsel->size)
		set_cursor_to_iter (GTK_TREE_VIEW (fontsel->size_list), &iter);
	}
    } else {
	GtkTreeIter iter;
	gboolean found = FALSE;
      
	gtk_tree_model_get_iter_first (GTK_TREE_MODEL (model), &iter);
	for (i = 0; i < G_N_ELEMENTS (font_sizes) && !found; i++) {
	    if (font_sizes[i] * PANGO_SCALE == fontsel->size) {
		set_cursor_to_iter (GTK_TREE_VIEW (fontsel->size_list), &iter);
		found = TRUE;
	    }

	    gtk_tree_model_iter_next (GTK_TREE_MODEL (model), &iter);
	}

	if (!found) {
	    GtkTreeSelection *selection = gtk_tree_view_get_selection (GTK_TREE_VIEW (fontsel->size_list));

	    gtk_tree_selection_unselect_all (selection);
	}
    }

    /* Set the entry to the new size, rounding to 1 digit,
     * trimming of trailing 0's and a trailing period
     */
    sprintf (buffer, "%.1f", fontsel->size / (1.0 * PANGO_SCALE));
    if (strchr (buffer, '.')) {
	p = buffer + strlen (buffer) - 1;
	while (*p == '0')
	    p--;
	if (*p == '.')
	    p--;
	p[1] = '\0';
    }

    /* Compare, to avoid moving the cursor unecessarily */
    if (strcmp (gtk_entry_get_text (GTK_ENTRY (fontsel->size_entry)), buffer) != 0)
	gtk_entry_set_text (GTK_ENTRY (fontsel->size_entry), buffer);
}

static void
gtk_fontsel_hack_select_best_size (GtkFontselHack *fontsel)
{
    gtk_fontsel_hack_load_font(fontsel);  
}

static void
gtk_fontsel_hack_set_size (GtkFontselHack *fontsel,
			   gint new_size)
{
    if (fontsel->size != new_size) {
	fontsel->size = new_size;

	gtk_fontsel_hack_show_available_sizes (fontsel, FALSE);      
	gtk_fontsel_hack_load_font (fontsel);
    }
}

/* If the user hits return in the font size entry, we change to the new font
   size. */

static void
gtk_fontsel_hack_size_activate (GtkWidget *w,
				gpointer data)
{
    GtkFontselHack *fontsel;
    gint new_size;
    const gchar *text;
  
    fontsel = GTK_FNTHACK(data);

    text = gtk_entry_get_text(GTK_ENTRY(fontsel->size_entry));
    new_size = MAX(0.1, atof(text) * PANGO_SCALE + 0.5);

    gtk_fontsel_hack_set_size(fontsel, new_size);
}

static gboolean
gtk_fontsel_hack_size_focus_out (GtkWidget*w,
				 GdkEventFocus *event,
				 gpointer data)
{
    gtk_fontsel_hack_size_activate(w, data);
  
    return TRUE;
}

/* This is called when a size is selected in the list. */

static void
gtk_fontsel_hack_select_size (GtkTreeSelection *selection,
			      gpointer data)
{
    GtkFontselHack *fontsel;
    GtkTreeModel *model;
    GtkTreeIter iter;
    gint new_size;
  
    fontsel = GTK_FNTHACK(data);
  
    if (gtk_tree_selection_get_selected(selection, &model, &iter)) {
	gtk_tree_model_get(model, &iter, SIZE_COLUMN, &new_size, -1);
	gtk_fontsel_hack_set_size(fontsel, new_size * PANGO_SCALE);
    }
}

static void
gtk_fontsel_hack_load_font (GtkFontselHack *fontsel)
{
    gtk_fontsel_hack_update_preview(fontsel);
}

static PangoFontDescription *
gtk_fontsel_hack_get_font_description (GtkFontselHack *fontsel)
{
    PangoFontDescription *font_desc = pango_font_face_describe(fontsel->face);

    pango_font_description_set_size(font_desc, fontsel->size);

    return font_desc;
}

/* This sets the font in the preview entry to the selected font, and
   tries to make sure that the preview entry is a reasonable size,
   i.e. so that the text can be seen with a bit of space to spare. But
   it tries to avoid resizing the entry every time the font changes.
   This is also used to shrink the preview if the font size was
   decreased, but that made it awkward if the users wanted to resize
   the window themselves. */

static void
gtk_fontsel_hack_update_preview (GtkFontselHack *fontsel)
{
#if GTK_MAJOR_VERSION < 3
    GtkRcStyle *rc_style;
#endif
    gint new_height;
    GtkRequisition req, old_requisition;
    GtkWidget *preview_entry = fontsel->preview_entry;
    const gchar *text;

    gtk_widget_get_child_requisition(preview_entry, &old_requisition);

#if GTK_MAJOR_VERSION >= 3
    gtk_widget_modify_font(preview_entry, gtk_fontsel_hack_get_font_description(fontsel));
#else  
    rc_style = gtk_rc_style_new();
    rc_style->font_desc = gtk_fontsel_hack_get_font_description(fontsel);
    gtk_widget_modify_style(preview_entry, rc_style);
    g_object_unref(rc_style);
#endif

    gtk_widget_size_request(preview_entry, NULL);

#if (GTK_MAJOR_VERSION == 2 && GTK_MINOR_VERSION < 20)
    req = preview_entry->requisition;
#else
    gtk_widget_get_requisition(preview_entry, &req);
#endif
  
    /* We don't ever want to be over MAX_PREVIEW_HEIGHT pixels high. */
    new_height = CLAMP(req.height, 
		       INITIAL_PREVIEW_HEIGHT, 
		       MAX_PREVIEW_HEIGHT);

    if (new_height > old_requisition.height || new_height < old_requisition.height - 30)
	gtk_widget_set_size_request(preview_entry, -1, new_height);
  
    /* This sets the preview text, if it hasn't been set already. */
    text = gtk_entry_get_text(GTK_ENTRY (preview_entry));
    if (*text == '\0') {
	gtk_entry_set_text(GTK_ENTRY(preview_entry), _(PREVIEW_TEXT));
    }
    gtk_editable_set_position(GTK_EDITABLE(preview_entry), 0);
}

/*****************************************************************************
 * These functions are the main public interface for getting/setting the font.
 *****************************************************************************/

gchar *
gtk_fontsel_hack_get_font_name (GtkFontselHack *fontsel)
{
    PangoFontDescription *font_desc = 
	gtk_fontsel_hack_get_font_description(fontsel);
    gchar *result;

    result = pango_font_description_to_string(font_desc);
    pango_font_description_free(font_desc);

    return result;
}

/* This sets the current font, selecting the appropriate list rows.
   First we check the fontname is valid and try to find the font
   family - i.e. the name in the main list. If we can't find that,
   then just return.  Next we try to set each of the properties
   according to the fontname.  Finally we select the font family &
   style in the lists. */

gboolean
gtk_fontsel_hack_set_font_name (GtkFontselHack *fontsel,
				const gchar *fontname)
{
    PangoFontFamily *new_family = NULL;
    PangoFontFace *new_face = NULL;
    PangoFontFace *fallback_face = NULL;
    PangoFontDescription *new_desc;
    GtkTreeModel *model;
    GtkTreeIter iter;
    GtkTreeIter match_iter;
    gboolean valid;
  
    g_return_val_if_fail(GTK_IS_FNTHACK(fontsel), FALSE);
  
    new_desc = pango_font_description_from_string(fontname);

    /* Check to make sure that this is in the list of allowed fonts */

    model = gtk_tree_view_get_model(GTK_TREE_VIEW(fontsel->family_list));
    for (valid = gtk_tree_model_get_iter_first (model, &iter);
	 valid;
	 valid = gtk_tree_model_iter_next(model, &iter)) {
	PangoFontFamily *family;
      
	gtk_tree_model_get(model, &iter, FAMILY_COLUMN, &family, -1);
      
	if (g_ascii_strcasecmp(pango_font_family_get_name(family),
			       pango_font_description_get_family(new_desc)) == 0)
	    new_family = family;
      
	g_object_unref(family);
      
	if (new_family)
	    break;
    }

    if (!new_family)
	return FALSE;

    fontsel->family = new_family;
    set_cursor_to_iter(GTK_TREE_VIEW(fontsel->family_list), &iter);
    gtk_fontsel_hack_show_available_styles(fontsel);

    model = gtk_tree_view_get_model (GTK_TREE_VIEW (fontsel->face_list));
    for (valid = gtk_tree_model_get_iter_first(model, &iter);
	 valid;
	 valid = gtk_tree_model_iter_next(model, &iter)) {
	PangoFontFace *face;
	PangoFontDescription *tmp_desc;
      
	gtk_tree_model_get(model, &iter, FACE_COLUMN, &face, -1);
	tmp_desc = pango_font_face_describe(face);
      
	if (font_description_style_equal(tmp_desc, new_desc))
	    new_face = face;
      
	if (!fallback_face) {
	    fallback_face = face;
	    match_iter = iter;
	}
      
	pango_font_description_free(tmp_desc);
	g_object_unref(face);
      
	if (new_face) {
	    match_iter = iter;
	    break;
	}
    }

    if (!new_face)
	new_face = fallback_face;

    fontsel->face = new_face;
    set_cursor_to_iter(GTK_TREE_VIEW(fontsel->face_list), &match_iter);  

    gtk_fontsel_hack_set_size(fontsel, pango_font_description_get_size(new_desc));
  
    g_object_freeze_notify(G_OBJECT(fontsel));
    g_object_notify(G_OBJECT(fontsel), "font_name");
    g_object_thaw_notify(G_OBJECT(fontsel));

    pango_font_description_free(new_desc);

    return TRUE;
}

gint
gtk_fontsel_hack_get_filter (GtkFontselHack *fontsel)
{
    return fontsel->filter;
}

void
gtk_fontsel_hack_set_filter (GtkFontselHack *fontsel, 
			     FontFilterType filter)
{
    fontsel->filter = filter;
    gtk_fontsel_hack_display_fonts(fontsel);
}

/* This returns the text in the preview entry. You should copy the returned
   text if you need it. */

G_CONST_RETURN gchar*
gtk_fontsel_hack_get_preview_text (GtkFontselHack *fontsel)
{
    return gtk_entry_get_text(GTK_ENTRY(fontsel->preview_entry));
}


/* This sets the text in the preview entry. */

void
gtk_fontsel_hack_set_preview_text (GtkFontselHack *fontsel,
				   const gchar *text)
{
    gtk_entry_set_text(GTK_ENTRY(fontsel->preview_entry), text);
}

/*****************************************************************************
 * GtkFontselHackDialog
 *****************************************************************************/

static void gtk_fontsel_hack_dialog_buildable_interface_init     (GtkBuildableIface *iface);
static GObject * gtk_fontsel_hack_dialog_buildable_get_internal_child (GtkBuildable *buildable,
								       GtkBuilder   *builder,
								       const gchar  *childname);

G_DEFINE_TYPE_WITH_CODE (GtkFontselHackDialog, gtk_fontsel_hack_dialog,
			 GTK_TYPE_DIALOG,
			 G_IMPLEMENT_INTERFACE(GTK_TYPE_BUILDABLE,
					       gtk_fontsel_hack_dialog_buildable_interface_init))

static GtkBuildableIface *parent_buildable_iface;

static void
gtk_fontsel_hack_dialog_buildable_interface_init (GtkBuildableIface *iface)
{
    parent_buildable_iface = g_type_interface_peek_parent(iface);
    iface->get_internal_child = gtk_fontsel_hack_dialog_buildable_get_internal_child;
}

static GObject *
gtk_fontsel_hack_dialog_buildable_get_internal_child (GtkBuildable *buildable,
						      GtkBuilder   *builder,
						      const gchar  *childname)
{
    if (strcmp(childname, "ok_button") == 0)
	return G_OBJECT(GTK_FNTHACK_DIALOG(buildable)->ok_button);
    else if (strcmp(childname, "cancel_button") == 0)
	return G_OBJECT(GTK_FNTHACK_DIALOG(buildable)->cancel_button);
    else if (strcmp(childname, "apply_button") == 0)
	return G_OBJECT(GTK_FNTHACK_DIALOG(buildable)->apply_button);
    else if (strcmp(childname, "font_selection") == 0)
	return G_OBJECT(GTK_FNTHACK_DIALOG(buildable)->fontsel);

    return parent_buildable_iface->get_internal_child(buildable, builder, childname);
}

static void
gtk_fontsel_hack_dialog_class_init (GtkFontselHackDialogClass *klass)
{
}

static void
gtk_fontsel_hack_dialog_init (GtkFontselHackDialog *fontseldiag)
{
    GtkDialog *dialog;

    gtk_widget_push_composite_child();

    dialog = GTK_DIALOG(fontseldiag);
  
    fontseldiag->dialog_width = -1;
    fontseldiag->auto_resize = TRUE;
  
    gtk_widget_set_events(GTK_WIDGET(fontseldiag), GDK_STRUCTURE_MASK);
    g_signal_connect(G_OBJECT(fontseldiag), "configure-event",
		     (GCallback) gtk_fontsel_hack_dialog_on_configure,
		     fontseldiag);
  
    gtk_container_set_border_width(GTK_CONTAINER(fontseldiag), 4);
    gtk_window_set_resizable(GTK_WINDOW(fontseldiag), FALSE);
  
    fontseldiag->main_vbox = gtk_dialog_get_content_area(dialog);
  
    fontseldiag->fontsel = gtk_fontsel_hack_new();
    gtk_container_set_border_width(GTK_CONTAINER(fontseldiag->fontsel), 4);
    gtk_widget_show(fontseldiag->fontsel);
    gtk_box_pack_start(GTK_BOX(fontseldiag->main_vbox),
		       fontseldiag->fontsel, TRUE, TRUE, 0);
  
    /* Create the action area */
    fontseldiag->action_area = gtk_dialog_get_action_area(dialog);

    fontseldiag->cancel_button = gtk_dialog_add_button(dialog,
						       GTK_STOCK_CANCEL,
						       GTK_RESPONSE_CANCEL);

    fontseldiag->apply_button = gtk_dialog_add_button(dialog,
						      GTK_STOCK_APPLY,
						      GTK_RESPONSE_APPLY);
    gtk_widget_hide(fontseldiag->apply_button);

    fontseldiag->ok_button = gtk_dialog_add_button(dialog,
						   GTK_STOCK_OK,
						   GTK_RESPONSE_OK);
    gtk_widget_grab_default(fontseldiag->ok_button);
  
    gtk_window_set_title(GTK_WINDOW(fontseldiag), _("Font Selection"));
#if GTK_MAJOR_VERSION < 3
    gtk_dialog_set_has_separator(GTK_DIALOG (dialog), FALSE);
#endif  
    gtk_widget_pop_composite_child();
}

GtkWidget*
gtk_fontsel_hack_dialog_new (const gchar *title)
{
    GtkFontselHackDialog *fontseldiag;

    fontseldiag = g_object_new(GTK_TYPE_FNTHACK_DIALOG, NULL);

    if (title)
	gtk_window_set_title(GTK_WINDOW(fontseldiag), title);

    return GTK_WIDGET(fontseldiag);
}

gchar*
gtk_fontsel_hack_dialog_get_font_name (GtkFontselHackDialog *fsd)
{
    return gtk_fontsel_hack_get_font_name(GTK_FNTHACK(fsd->fontsel));
}

gboolean
gtk_fontsel_hack_dialog_set_font_name (GtkFontselHackDialog *fsd,
				       const gchar *fontname)
{
    return gtk_fontsel_hack_set_font_name(GTK_FNTHACK(fsd->fontsel), fontname);
}

gint
gtk_fontsel_hack_dialog_get_filter (GtkFontselHackDialog *fsd)
{
    return gtk_fontsel_hack_get_filter (GTK_FNTHACK (fsd->fontsel));
}

void
gtk_fontsel_hack_dialog_set_filter (GtkFontselHackDialog *fsd,
				    FontFilterType	  filter)
{
    gtk_fontsel_hack_set_filter(GTK_FNTHACK(fsd->fontsel), filter);
}

G_CONST_RETURN gchar*
gtk_fontsel_hack_dialog_get_preview_text (GtkFontselHackDialog *fsd)
{
    return gtk_fontsel_hack_get_preview_text(GTK_FNTHACK(fsd->fontsel));
}

void
gtk_fontsel_hack_dialog_set_preview_text (GtkFontselHackDialog *fsd,
					  const gchar *text)
{
    gtk_fontsel_hack_set_preview_text(GTK_FNTHACK(fsd->fontsel), text);
}

GtkWidget *gtk_fontsel_hack_dialog_ok_button (GtkWidget *fsd)
{
    return GTK_FNTHACK_DIALOG(fsd)->ok_button;
}

GtkWidget *gtk_fontsel_hack_dialog_cancel_button (GtkWidget *fsd)
{
    return GTK_FNTHACK_DIALOG(fsd)->cancel_button;
}

/* This turns auto-shrink off if the user resizes the width of the dialog.
   It also turns it back on again if the user resizes it back to its normal
   width. */

static gint
gtk_fontsel_hack_dialog_on_configure (GtkWidget *widget,
				      GdkEventConfigure *event,
				      GtkFontselHackDialog *fsd)
{
    /* This sets the initial width. */
    if (fsd->dialog_width == -1) {
	fsd->dialog_width = event->width;
    } else if (fsd->auto_resize && fsd->dialog_width != event->width) {
	fsd->auto_resize = FALSE;
	gtk_window_set_resizable(GTK_WINDOW(fsd), FALSE);
    } else if (!fsd->auto_resize && fsd->dialog_width == event->width) {
	fsd->auto_resize = TRUE;
	gtk_window_set_resizable(GTK_WINDOW(fsd), FALSE);
    }
  
    return FALSE;
}
