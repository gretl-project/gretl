/* -*- Mode: C; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8; coding: utf-8 -*- 
 *  gtksourceview.c
 *
 *  Copyright (C) 2001 - Mikael Hermansson <tyan@linux.se> and
 *  Chris Phelps <chicane@reninet.com>
 *
 *  Copyright (C) 2002 - Jeroen Zwartepoorte
 *
 *  Copyright (C) 2003 - Gustavo Giráldez and Paolo Maggi 
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Library General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU Library General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h> /* For strlen */

#include <gtk/gtk.h>
#include <gdk/gdkkeysyms.h>
#include <pango/pango-tabs.h>

#include "gtksourceview-marshal.h"
#include "gtksourceview.h"

/*
#define ENABLE_DEBUG
*/
#undef ENABLE_DEBUG

/*
#define ENABLE_PROFILE
*/
#undef ENABLE_PROFILE

#ifdef ENABLE_DEBUG
#define DEBUG(x) (x)
#else
#define DEBUG(x)
#endif

#ifdef ENABLE_PROFILE
#define PROFILE(x) (x)
#else
#define PROFILE(x)
#endif

#define COMPOSITE_ALPHA                 225
#define GUTTER_PIXMAP 			16
#define DEFAULT_TAB_WIDTH 		8
#define MIN_NUMBER_WINDOW_WIDTH		20
#define MAX_TAB_WIDTH			32

/* Signals */
enum {
	UNDO,
	REDO,
	LAST_SIGNAL
};

/* Properties */
enum {
	PROP_0,
	PROP_SHOW_LINE_NUMBERS,
	PROP_TABS_WIDTH,
	PROP_AUTO_INDENT,
	PROP_INSERT_SPACES,
	PROP_SMART_HOME_END,
	PROP_INDENT_ON_TAB
};

struct _GtkSourceViewPrivate
{
	guint		 tabs_width;
	gboolean 	 show_line_numbers;
	gboolean	 auto_indent;
	gboolean	 insert_spaces;
	gboolean	 indent_on_tab;
	gboolean	 smart_home_end;
	
	GtkSourceBuffer *source_buffer;
	gint		 old_lines;
};

/* Implement DnD for application/x-color drops */
typedef enum {
	TARGET_COLOR = 200
} GtkSourceViewDropTypes;

static const GtkTargetEntry drop_types[] = {
	{"application/x-color", 0, TARGET_COLOR}
};

static guint signals[LAST_SIGNAL] = { 0 };

static GObjectClass *parent_class = NULL;

/* Prototypes. */
static void	gtk_source_view_class_init 		(GtkSourceViewClass *klass);
static void	gtk_source_view_init 			(GtkSourceView      *view);
static void 	gtk_source_view_finalize 		(GObject            *object);

static void	gtk_source_view_undo 			(GtkSourceView      *view);
static void	gtk_source_view_redo 			(GtkSourceView      *view);

static void 	set_source_buffer 			(GtkSourceView      *view,
			       				 GtkTextBuffer      *buffer);

static void	gtk_source_view_populate_popup 		(GtkTextView        *view,
					    		 GtkMenu            *menu);

static void	gtk_source_view_move_cursor		(GtkTextView        *text_view,
							 GtkMovementStep     step,
							 gint                count,
							 gboolean            extend_selection);

static void 	menu_item_activate_cb 			(GtkWidget          *menu_item,
				  			 GtkTextView        *text_view);

static void 	gtk_source_view_get_lines 		(GtkTextView       *text_view,
				       			 gint               first_y,
				       			 gint               last_y,
				       			 GArray            *buffer_coords,
				       			 GArray            *numbers,
				       			 gint              *countp);
static gint     gtk_source_view_expose 			(GtkWidget         *widget,
							 GdkEventExpose    *event);
static gboolean	gtk_source_view_key_press_event		(GtkWidget         *widget,
							 GdkEventKey       *event);
static gboolean	gtk_source_view_button_press_event	(GtkWidget         *widget,
							 GdkEventButton    *event);
static void 	view_dnd_drop 				(GtkTextView       *view, 
							 GdkDragContext    *context,
							 gint               x,
							 gint               y,
							 GtkSelectionData  *selection_data,
							 guint              info,
							 guint              time,
							 gpointer           data);

static gint	calculate_real_tab_width 		(GtkSourceView     *view, 
							 guint              tab_size,
							 gchar              c);

static void	gtk_source_view_set_property 		(GObject           *object,
							 guint              prop_id,
							 const GValue      *value,
							 GParamSpec        *pspec);
static void	gtk_source_view_get_property		(GObject           *object,
							 guint              prop_id,
							 GValue            *value,
							 GParamSpec        *pspec);

static void     gtk_source_view_style_set               (GtkWidget         *widget,
							 GtkStyle          *previous_style);


/* Private functions. */
static void
gtk_source_view_class_init (GtkSourceViewClass *klass)
{
	GObjectClass	 *object_class;
	GtkTextViewClass *textview_class;
	GtkBindingSet    *binding_set;
	GtkWidgetClass   *widget_class;
	
	object_class 	= G_OBJECT_CLASS (klass);
	textview_class 	= GTK_TEXT_VIEW_CLASS (klass);
	parent_class 	= g_type_class_peek_parent (klass);
	widget_class 	= GTK_WIDGET_CLASS (klass);
	
	object_class->finalize = gtk_source_view_finalize;
	object_class->get_property = gtk_source_view_get_property;
	object_class->set_property = gtk_source_view_set_property;

	widget_class->key_press_event = gtk_source_view_key_press_event;
	widget_class->button_press_event = gtk_source_view_button_press_event;
	widget_class->expose_event = gtk_source_view_expose;
	widget_class->style_set = gtk_source_view_style_set;

	textview_class->populate_popup = gtk_source_view_populate_popup;
	textview_class->move_cursor = gtk_source_view_move_cursor;
	
	klass->undo = gtk_source_view_undo;
	klass->redo = gtk_source_view_redo;

	g_object_class_install_property (object_class,
					 PROP_SHOW_LINE_NUMBERS,
					 g_param_spec_boolean ("show_line_numbers",
							       "Show Line Numbers",
							       "Whether to display line numbers",
							       FALSE,
							       G_PARAM_READWRITE));
	
	g_object_class_install_property (object_class,
					 PROP_TABS_WIDTH,
					 g_param_spec_uint ("tabs_width",
							    "Tabs Width",
							    "Tabs Width",
							    1,
							    MAX_TAB_WIDTH,
							    DEFAULT_TAB_WIDTH,
							    G_PARAM_READWRITE));
	
	g_object_class_install_property (object_class,
					 PROP_AUTO_INDENT,
					 g_param_spec_boolean ("auto_indent",
							       "Auto Indentation",
							       "Whether to enable auto indentation",
							       FALSE,
							       G_PARAM_READWRITE));
	g_object_class_install_property (object_class,
					 PROP_INSERT_SPACES,
					 g_param_spec_boolean ("insert_spaces_instead_of_tabs",
							       "Insert Spaces Instead of Tabs",
							       "Whether to insert spaces instead of tabs",
							       FALSE,
							       G_PARAM_READWRITE));

	g_object_class_install_property (object_class,
					 PROP_SMART_HOME_END,
					 g_param_spec_boolean ("smart_home_end",
							       "Use smart home/end",
							       "HOME and END keys move to first/last "
								 "non whitespace characters on line before going "
								 "to the start/end of the line",
							       TRUE,
							       G_PARAM_READWRITE));

	g_object_class_install_property (object_class,
					 PROP_INDENT_ON_TAB,
					 g_param_spec_boolean ("indent_on_tab",
							       "Indent on tab",
							       "Whether to indent the selected text when the tab key is pressed",
							       FALSE,
							       G_PARAM_READWRITE));


	signals [UNDO] =
		g_signal_new ("undo",
			      G_TYPE_FROM_CLASS (klass),
			      G_SIGNAL_RUN_LAST | G_SIGNAL_ACTION,
			      G_STRUCT_OFFSET (GtkSourceViewClass, undo),
			      NULL,
			      NULL,
			      gtksourceview_marshal_VOID__VOID,
			      G_TYPE_NONE,
			      0);
	signals [REDO] =
		g_signal_new ("redo",
			      G_TYPE_FROM_CLASS (klass),
			      G_SIGNAL_RUN_LAST | G_SIGNAL_ACTION,
			      G_STRUCT_OFFSET (GtkSourceViewClass, redo),
			      NULL,
			      NULL,
			      gtksourceview_marshal_VOID__VOID,
			      G_TYPE_NONE,
			      0);

	binding_set = gtk_binding_set_by_class (klass);

	gtk_binding_entry_add_signal (binding_set,
				      GDK_z,
				      GDK_CONTROL_MASK,
				      "undo", 0);
	gtk_binding_entry_add_signal (binding_set,
				      GDK_z,
				      GDK_CONTROL_MASK | GDK_SHIFT_MASK,
				      "redo", 0);
	gtk_binding_entry_add_signal (binding_set,
				      GDK_F14,
				      0,
				      "undo", 0);
}

static void 
gtk_source_view_set_property (GObject      *object,
			      guint         prop_id,
			      const GValue *value,
			      GParamSpec   *pspec)
{
	GtkSourceView *view;
	
	g_return_if_fail (GTK_IS_SOURCE_VIEW (object));

	view = GTK_SOURCE_VIEW (object);
    
	switch (prop_id)
	{
		case PROP_SHOW_LINE_NUMBERS:
			gtk_source_view_set_show_line_numbers (view,
							       g_value_get_boolean (value));
			break;
			
		case PROP_TABS_WIDTH:
			gtk_source_view_set_tabs_width (view, 
							g_value_get_uint (value));
			break;
			
		case PROP_AUTO_INDENT:
			gtk_source_view_set_auto_indent (view,
							 g_value_get_boolean (value));
			break;
			
		case PROP_INSERT_SPACES:
			gtk_source_view_set_insert_spaces_instead_of_tabs (
							view,
							g_value_get_boolean (value));
			break;
			
		case PROP_SMART_HOME_END:
			gtk_source_view_set_smart_home_end (view,
							    g_value_get_boolean (value));
			break;

		case PROP_INDENT_ON_TAB:
			gtk_source_view_set_indent_on_tab (view,
							   g_value_get_boolean (value));
			break;
		
		default:
			G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
			break;
	}
}

static void 
gtk_source_view_get_property (GObject    *object,
			      guint       prop_id,
			      GValue     *value,
			      GParamSpec *pspec)
{
	GtkSourceView *view;
	
	g_return_if_fail (GTK_IS_SOURCE_VIEW (object));

	view = GTK_SOURCE_VIEW (object);
    
	switch (prop_id)
	{
		case PROP_SHOW_LINE_NUMBERS:
			g_value_set_boolean (value,
					     gtk_source_view_get_show_line_numbers (view));
					     
			break;
			
		case PROP_TABS_WIDTH:
			g_value_set_uint (value,
					  gtk_source_view_get_tabs_width (view));
			break;
			
		case PROP_AUTO_INDENT:
			g_value_set_boolean (value,
					     gtk_source_view_get_auto_indent (view));

			break;
			
		case PROP_INSERT_SPACES:
			g_value_set_boolean (value,
					     gtk_source_view_get_insert_spaces_instead_of_tabs (view));
	
			break;

		case PROP_SMART_HOME_END:
			g_value_set_boolean (value,
					     gtk_source_view_get_smart_home_end (view));
			break;
			
		case PROP_INDENT_ON_TAB:
			g_value_set_boolean (value,
					     gtk_source_view_get_indent_on_tab (view));
			break;

		default:
			G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
			break;
	}
}

static void
gtk_source_view_init (GtkSourceView *view)
{
	GtkTargetList *tl;

	view->priv = g_new0 (GtkSourceViewPrivate, 1);

	view->priv->tabs_width = DEFAULT_TAB_WIDTH;
	view->priv->smart_home_end = TRUE;
	
	gtk_text_view_set_left_margin (GTK_TEXT_VIEW (view), 2);
	gtk_text_view_set_right_margin (GTK_TEXT_VIEW (view), 2);

	tl = gtk_drag_dest_get_target_list (GTK_WIDGET (view));
	g_return_if_fail (tl != NULL);

	gtk_target_list_add_table (tl, drop_types, G_N_ELEMENTS (drop_types));

	g_signal_connect (G_OBJECT (view), 
			  "drag_data_received", 
			  G_CALLBACK (view_dnd_drop), 
			  NULL);
}

static void
gtk_source_view_finalize (GObject *object)
{
	GtkSourceView *view;

	g_return_if_fail (object != NULL);
	g_return_if_fail (GTK_IS_SOURCE_VIEW (object));

	view = GTK_SOURCE_VIEW (object);

	set_source_buffer (view, NULL);

	g_free (view->priv);

	G_OBJECT_CLASS (parent_class)->finalize (object);
}

static void
set_source_buffer (GtkSourceView *view, GtkTextBuffer *buffer)
{
	/* keep our pointer to the source buffer in sync with
	 * textview's, though it would be a lot nicer if GtkTextView
	 * had a "set_buffer" signal */
	/* FIXME: in gtk 2.3 we have a buffer property so we can
	 * connect to the notify signal.  Unfortunately we can't
	 * depend on gtk 2.3 yet (see bug #108353) */
	if (view->priv->source_buffer) 
	{
		g_object_remove_weak_pointer (G_OBJECT (view->priv->source_buffer),
					      (gpointer *) &view->priv->source_buffer);
	}

	if (buffer && GTK_IS_SOURCE_BUFFER (buffer)) 
	{
		view->priv->source_buffer = GTK_SOURCE_BUFFER (buffer);
		g_object_add_weak_pointer (G_OBJECT (buffer),
					   (gpointer *) &view->priv->source_buffer);
	}
	else 
	{
		view->priv->source_buffer = NULL;
	}
}

static void
gtk_source_view_undo (GtkSourceView *view)
{
	GtkTextBuffer *buffer;

	g_return_if_fail (GTK_IS_SOURCE_VIEW (view));

	buffer = gtk_text_view_get_buffer (GTK_TEXT_VIEW (view));

	if (gtk_source_buffer_can_undo (GTK_SOURCE_BUFFER (buffer)))
	{
		gtk_source_buffer_undo (GTK_SOURCE_BUFFER (buffer));
		gtk_text_view_scroll_mark_onscreen (GTK_TEXT_VIEW (view),
						    gtk_text_buffer_get_insert (buffer));
	}
}

static void
gtk_source_view_redo (GtkSourceView *view)
{
	GtkTextBuffer *buffer;

	g_return_if_fail (GTK_IS_SOURCE_VIEW (view));

	buffer = gtk_text_view_get_buffer (GTK_TEXT_VIEW (view));

	if (gtk_source_buffer_can_redo (GTK_SOURCE_BUFFER (buffer)))
	{
		gtk_source_buffer_redo (GTK_SOURCE_BUFFER (buffer));
		gtk_text_view_scroll_mark_onscreen (GTK_TEXT_VIEW (view),
						    gtk_text_buffer_get_insert (buffer));
	}
}

static void
gtk_source_view_populate_popup (GtkTextView *text_view,
				GtkMenu     *menu)
{
	GtkTextBuffer *buffer;
	GtkWidget *menu_item;

	buffer = gtk_text_view_get_buffer (text_view);
	if (!buffer && !GTK_IS_SOURCE_BUFFER (buffer))
		return;

	/* separator */
	menu_item = gtk_menu_item_new ();
	gtk_menu_shell_prepend (GTK_MENU_SHELL (menu), menu_item);
	gtk_widget_show (menu_item);

	/* create redo menu_item. */
	menu_item = gtk_image_menu_item_new_from_stock ("gtk-redo", NULL);
	g_object_set_data (G_OBJECT (menu_item), "gtk-signal", "redo");
	g_signal_connect (G_OBJECT (menu_item), "activate",
			  G_CALLBACK (menu_item_activate_cb), text_view);
	gtk_menu_shell_prepend (GTK_MENU_SHELL (menu), menu_item);
	gtk_widget_set_sensitive (menu_item,
				  gtk_source_buffer_can_redo (GTK_SOURCE_BUFFER (buffer)));
	gtk_widget_show (menu_item);

	/* create undo menu_item. */
	menu_item = gtk_image_menu_item_new_from_stock ("gtk-undo", NULL);
	g_object_set_data (G_OBJECT (menu_item), "gtk-signal", "undo");
	g_signal_connect (G_OBJECT (menu_item), "activate",
			  G_CALLBACK (menu_item_activate_cb), text_view);
	gtk_menu_shell_prepend (GTK_MENU_SHELL (menu), menu_item);
	gtk_widget_set_sensitive (menu_item, 
				  gtk_source_buffer_can_undo (GTK_SOURCE_BUFFER (buffer)));
	gtk_widget_show (menu_item);

}

static void
move_cursor (GtkTextView       *text_view,
	     const GtkTextIter *new_location,
	     gboolean           extend_selection)
{
	GtkTextBuffer *buffer = text_view->buffer;

	if (extend_selection)
		gtk_text_buffer_move_mark_by_name (buffer, "insert",
						   new_location);
	else
		gtk_text_buffer_place_cursor (buffer, new_location);

	gtk_text_view_scroll_mark_onscreen (text_view,
					    gtk_text_buffer_get_insert (buffer));
}

static void
gtk_source_view_move_cursor (GtkTextView    *text_view,
			     GtkMovementStep step,
			     gint            count,
			     gboolean        extend_selection)
{
	GtkSourceView *source_view = GTK_SOURCE_VIEW (text_view);
	GtkTextBuffer *buffer = text_view->buffer;
	GtkTextMark *mark;
	GtkTextIter cur, iter;

	mark = gtk_text_buffer_get_insert (buffer);
	gtk_text_buffer_get_iter_at_mark (buffer, &cur, mark);
	iter = cur;

	if (step == GTK_MOVEMENT_DISPLAY_LINE_ENDS &&
	    source_view->priv->smart_home_end && count == -1)
	{
		/* Find the iter of the first character on the line. */
		gtk_text_iter_set_line_offset (&cur, 0);
		while (!gtk_text_iter_ends_line (&cur))
		{
			gunichar c = gtk_text_iter_get_char (&cur);
			if (g_unichar_isspace (c))
				gtk_text_iter_forward_char (&cur);
			else
				break;
		}

		if (gtk_text_iter_starts_line (&iter) ||
		    !gtk_text_iter_equal (&cur, &iter))
		{
			move_cursor (text_view, &cur, extend_selection);
		}
		else
		{
			gtk_text_iter_set_line_offset (&cur, 0);
			move_cursor (text_view, &cur, extend_selection);
		}
	}
	else if (step == GTK_MOVEMENT_DISPLAY_LINE_ENDS &&
		 source_view->priv->smart_home_end && count == 1)
	{
		/* Find the iter of the last character on the line. */
		if (!gtk_text_iter_ends_line (&cur))
			gtk_text_iter_forward_to_line_end (&cur);
		while (!gtk_text_iter_starts_line (&cur))
		{
			gunichar c;
			gtk_text_iter_backward_char (&cur);
			c = gtk_text_iter_get_char (&cur);
			if (!g_unichar_isspace (c))
			{
				/* We've gone one character too far. */
				gtk_text_iter_forward_char (&cur);
				break;
			}
		}

		if (gtk_text_iter_ends_line (&iter) ||
		    !gtk_text_iter_equal (&cur, &iter))
		{
			move_cursor (text_view, &cur, extend_selection);
		}
		else
		{
			gtk_text_iter_forward_to_line_end (&cur);
			move_cursor (text_view, &cur, extend_selection);
		}
	}
	else
	{
		GTK_TEXT_VIEW_CLASS (parent_class)->move_cursor (text_view,
								 step, count,
								 extend_selection);
	}
}

static void
menu_item_activate_cb (GtkWidget   *menu_item,
		       GtkTextView *text_view)
{
	const gchar *signal;

	signal = g_object_get_data (G_OBJECT (menu_item), "gtk-signal");
	g_signal_emit_by_name (G_OBJECT (text_view), signal);
}

/* This function is taken from gtk+/tests/testtext.c */
static void
gtk_source_view_get_lines (GtkTextView  *text_view,
			   gint          first_y,
			   gint          last_y,
			   GArray       *buffer_coords,
			   GArray       *numbers,
			   gint         *countp)
{
	GtkTextIter iter;
	gint count;
	gint size;
      	gint last_line_num = -1;	

	g_array_set_size (buffer_coords, 0);
	g_array_set_size (numbers, 0);
  
	/* Get iter at first y */
	gtk_text_view_get_line_at_y (text_view, &iter, first_y, NULL);

	/* For each iter, get its location and add it to the arrays.
	 * Stop when we pass last_y
	*/
	count = 0;
  	size = 0;

  	while (!gtk_text_iter_is_end (&iter))
    	{
		gint y, height;
      
		gtk_text_view_get_line_yrange (text_view, &iter, &y, &height);

		g_array_append_val (buffer_coords, y);
		last_line_num = gtk_text_iter_get_line (&iter);
		g_array_append_val (numbers, last_line_num);
      	
		++count;

		if ((y + height) >= last_y)
			break;
      
		gtk_text_iter_forward_line (&iter);
	}

	if (gtk_text_iter_is_end (&iter))
    	{
		gint y, height;
		gint line_num;
      
		gtk_text_view_get_line_yrange (text_view, &iter, &y, &height);

		line_num = gtk_text_iter_get_line (&iter);

		if (line_num != last_line_num)
		{
			g_array_append_val (buffer_coords, y);
			g_array_append_val (numbers, line_num);
			++count;
		}
	}

	*countp = count;
}

static void
gtk_source_view_paint_margin (GtkSourceView *view,
			      GdkEventExpose *event)
{
	GtkTextView *text_view;
	GdkWindow *win;
	PangoLayout *layout;
	GArray *numbers;
	GArray *pixels;
	gchar str [8];  /* we don't expect more than ten million lines ;-) */
	gint y1, y2;
	gint count;
	gint margin_width;
	gint text_width;
	gint i;
	GtkTextIter cur;
	gint cur_line;

	text_view = GTK_TEXT_VIEW (view);

	if (!view->priv->show_line_numbers)
	{
		gtk_text_view_set_border_window_size (GTK_TEXT_VIEW (text_view),
						      GTK_TEXT_WINDOW_LEFT,
						      0);

		return;
	}

	win = gtk_text_view_get_window (text_view,
					GTK_TEXT_WINDOW_LEFT);

	y1 = event->area.y;
	y2 = y1 + event->area.height;

	/* get the extents of the line printing */
	gtk_text_view_window_to_buffer_coords (text_view,
					       GTK_TEXT_WINDOW_LEFT,
					       0,
					       y1,
					       NULL,
					       &y1);

	gtk_text_view_window_to_buffer_coords (text_view,
					       GTK_TEXT_WINDOW_LEFT,
					       0,
					       y2,
					       NULL,
					       &y2);

	numbers = g_array_new (FALSE, FALSE, sizeof (gint));
	pixels = g_array_new (FALSE, FALSE, sizeof (gint));

	/* get the line numbers and y coordinates. */
	gtk_source_view_get_lines (text_view,
				   y1,
				   y2,
				   pixels,
				   numbers,
				   &count);

	/* A zero-lined document should display a "1"; we don't need to worry about
	scrolling effects of the text widget in this special case */
	
	if (count == 0)
	{
		gint y = 0;
		gint n = 0;
		count = 1;
		g_array_append_val (pixels, y);
		g_array_append_val (numbers, n);
	}

	DEBUG ({
		g_message ("Painting line numbers %d - %d",
			   g_array_index (numbers, gint, 0),
			   g_array_index (numbers, gint, count - 1));
	});
	
	/* set size. */
	g_snprintf (str, sizeof (str),
		    "%d", MAX (99, gtk_text_buffer_get_line_count (text_view->buffer)));
	layout = gtk_widget_create_pango_layout (GTK_WIDGET (view), str);

	pango_layout_get_pixel_size (layout, &text_width, NULL);
	
	pango_layout_set_width (layout, text_width);
	pango_layout_set_alignment (layout, PANGO_ALIGN_RIGHT);

	/* determine the width of the left margin. */
	if (view->priv->show_line_numbers)
		margin_width = text_width + 4;
	else
		margin_width = 0;

	g_return_if_fail (margin_width != 0);
	
	gtk_text_view_set_border_window_size (GTK_TEXT_VIEW (text_view),
					      GTK_TEXT_WINDOW_LEFT,
					      margin_width);

	i = 0;
	
	gtk_text_buffer_get_iter_at_mark (text_view->buffer, 
					  &cur, 
					  gtk_text_buffer_get_insert (text_view->buffer));

	cur_line = gtk_text_iter_get_line (&cur) + 1;
	
	while (i < count) 
	{
		gint pos;
		
		gtk_text_view_buffer_to_window_coords (text_view,
						       GTK_TEXT_WINDOW_LEFT,
						       0,
						       g_array_index (pixels, gint, i),
						       NULL,
						       &pos);

		if (view->priv->show_line_numbers) 
		{
			gint line_to_paint = g_array_index (numbers, gint, i) + 1;
			
			if (line_to_paint == cur_line)
			{
				gchar *markup;
				markup = g_strdup_printf ("<b>%d</b>", line_to_paint);
				
				pango_layout_set_markup (layout, markup, -1);

				g_free (markup);
			}
			else
			{
				g_snprintf (str, sizeof (str),
					    "%d", line_to_paint);

				pango_layout_set_markup (layout, str, -1);
			}

			gtk_paint_layout (GTK_WIDGET (view)->style,
					  win,
					  GTK_WIDGET_STATE (view),
					  FALSE,
					  NULL,
					  GTK_WIDGET (view),
					  NULL,
					  text_width + 2, 
					  pos,
					  layout);
		}

		++i;
	}

	g_array_free (pixels, TRUE);
	g_array_free (numbers, TRUE);

	g_object_unref (G_OBJECT (layout));
}

static gint
gtk_source_view_expose (GtkWidget      *widget,
			GdkEventExpose *event)
{
	GtkSourceView *view;
	GtkTextView *text_view;
	gboolean event_handled;
	
	view = GTK_SOURCE_VIEW (widget);
	text_view = GTK_TEXT_VIEW (widget);

	event_handled = FALSE;
	
	/* maintain the our source_buffer pointer synchronized */
	if (text_view->buffer != GTK_TEXT_BUFFER (view->priv->source_buffer) &&
	    GTK_IS_SOURCE_BUFFER (text_view->buffer)) 
	{
		set_source_buffer (view, text_view->buffer);
	}
	
	/* check if the expose event is for the text window first, and
	 * make sure the visible region is highlighted */
	if (event->window == gtk_text_view_get_window (text_view, GTK_TEXT_WINDOW_TEXT) &&
	    view->priv->source_buffer != NULL) 
	{
		GdkRectangle visible_rect;
		GtkTextIter iter1, iter2;
		
		gtk_text_view_get_visible_rect (text_view, &visible_rect);
		gtk_text_view_get_line_at_y (text_view, &iter1,
					     visible_rect.y, NULL);
		gtk_text_iter_backward_line (&iter1);
		gtk_text_view_get_line_at_y (text_view, &iter2,
					     visible_rect.y
					     + visible_rect.height, NULL);
		gtk_text_iter_forward_line (&iter2);

		_gtk_source_buffer_highlight_region (view->priv->source_buffer,
						     &iter1, &iter2, FALSE);
	}

	/* now check for the left window, which contains the margin */
	if (event->window == gtk_text_view_get_window (text_view,
						       GTK_TEXT_WINDOW_LEFT)) 
	{
		gtk_source_view_paint_margin (view, event);
		event_handled = TRUE;
	} 
	else 
	{
		gint lines;

		/* FIXME: could it be a performances problem? - Paolo */
		lines = gtk_text_buffer_get_line_count (text_view->buffer);

		if (view->priv->old_lines != lines)
		{
			GdkWindow *w;
			view->priv->old_lines = lines;

			w = gtk_text_view_get_window (text_view, GTK_TEXT_WINDOW_LEFT);

			if (w != NULL)
				gdk_window_invalidate_rect (w, NULL, FALSE);
		}

		/* Have GtkTextView draw the text first. */
		if (GTK_WIDGET_CLASS (parent_class)->expose_event)
			event_handled = 
				(* GTK_WIDGET_CLASS (parent_class)->expose_event)
				(widget, event);
	}
	
	return event_handled;	
}

/*
 *This is a pretty important function...we call it when the tab_stop is changed,
 *And when the font is changed.
 *NOTE: You must change this with the default font for now...
 *It may be a good idea to set the tab_width for each GtkTextTag as well
 *based on the font that we set at creation time
 *something like style_cache_set_tabs_from_font or the like.
 *Now, this *may* not be necessary because most tabs wont be inside of another highlight,
 *except for things like multi-line comments (which will usually have an italic font which
 *may or may not be a different size than the standard one), or if some random language
 *definition decides that it would be spiffy to have a bg color for "start of line" whitespace
 *"^\(\t\| \)+" would probably do the trick for that.
 */
static gint
calculate_real_tab_width (GtkSourceView *view, guint tab_size, gchar c)
{
	PangoLayout *layout;
	gchar *tab_string;
	gint tab_width = 0;

	if (tab_size == 0)
		return -1;

	tab_string = g_strnfill (tab_size, c);
	layout = gtk_widget_create_pango_layout (GTK_WIDGET (view), tab_string);
	g_free (tab_string);

	if (layout != NULL) {
		pango_layout_get_pixel_size (layout, &tab_width, NULL);
		g_object_unref (G_OBJECT (layout));
	} else
		tab_width = -1;

	return tab_width;
}


/* ----------------------------------------------------------------------
 * Public interface 
 * ---------------------------------------------------------------------- */

/**
 * gtk_source_view_new:
 *
 * Creates a new #GtkSourceView. An empty default buffer will be
 * created for you. If you want to specify your own buffer, consider
 * gtk_source_view_new_with_buffer().
 *
 * Return value: a new #GtkSourceView
 **/
GtkWidget *
gtk_source_view_new ()
{
	GtkWidget *widget;
	GtkSourceBuffer *buffer;

	buffer = gtk_source_buffer_new (NULL);
	widget = gtk_source_view_new_with_buffer (buffer);
	g_object_unref (buffer);
	return widget;
}

/**
 * gtk_source_view_new_with_buffer:
 * @buffer: a #GtkSourceBuffer.
 *
 * Creates a new #GtkSourceView widget displaying the buffer
 * @buffer. One buffer can be shared among many widgets.
 *
 * Return value: a new #GtkTextView.
 **/
GtkWidget *
gtk_source_view_new_with_buffer (GtkSourceBuffer *buffer)
{
	GtkWidget *view;

	g_return_val_if_fail (buffer != NULL && GTK_IS_SOURCE_BUFFER (buffer), NULL);
	
	view = g_object_new (GTK_TYPE_SOURCE_VIEW, NULL);
	gtk_text_view_set_buffer (GTK_TEXT_VIEW (view), GTK_TEXT_BUFFER (buffer));

	return view;
}

GType
gtk_source_view_get_type (void)
{
	static GType our_type = 0;

	if (our_type == 0) {
		static const GTypeInfo our_info = {
			sizeof (GtkSourceViewClass),
			(GBaseInitFunc) NULL,
			(GBaseFinalizeFunc) NULL,
			(GClassInitFunc) gtk_source_view_class_init,
			NULL,	/* class_finalize */
			NULL,	/* class_data */
			sizeof (GtkSourceView),
			0,	/* n_preallocs */
			(GInstanceInitFunc) gtk_source_view_init
		};

		our_type = g_type_register_static (GTK_TYPE_TEXT_VIEW,
						   "GtkSourceView",
						   &our_info, 0);
	}

	return our_type;
}

/**
 * gtk_source_view_get_show_line_numbers:
 * @view: a #GtkSourceView.
 *
 * Returns whether line numbers are displayed beside the text.
 *
 * Return value: %TRUE if the line numbers are displayed.
 **/
gboolean
gtk_source_view_get_show_line_numbers (GtkSourceView *view)
{
	g_return_val_if_fail (view != NULL, FALSE);
	g_return_val_if_fail (GTK_IS_SOURCE_VIEW (view), FALSE);

	return view->priv->show_line_numbers;
}

/**
 * gtk_source_view_set_show_line_numbers:
 * @view: a #GtkSourceView.
 * @show: whether line numbers should be displayed.
 *
 * If %TRUE line numbers will be displayed beside the text.
 *
 **/
void
gtk_source_view_set_show_line_numbers (GtkSourceView *view,
				       gboolean       show)
{
	g_return_if_fail (view != NULL);
	g_return_if_fail (GTK_IS_SOURCE_VIEW (view));

	show = (show != FALSE);

	if (show) 
	{
		if (!view->priv->show_line_numbers) 
		{
			/* Set left margin to minimum width if no margin is 
			   visible yet. Otherwise, just queue a redraw, so the
			   expose handler will automatically adjust the margin. */
			gtk_text_view_set_border_window_size (GTK_TEXT_VIEW (view),
							      GTK_TEXT_WINDOW_LEFT,
							      MIN_NUMBER_WINDOW_WIDTH);

			view->priv->show_line_numbers = show;

			g_object_notify (G_OBJECT (view), "show_line_numbers");
		}
	}
	else 
	{
		if (view->priv->show_line_numbers) 
		{
			view->priv->show_line_numbers = show;

			/* force expose event, which will adjust margin. */
			gtk_widget_queue_draw (GTK_WIDGET (view));

			g_object_notify (G_OBJECT (view), "show_line_numbers");
		}
	}
}

/**
 * gtk_source_view_get_tabs_width:
 * @view: a #GtkSourceView.
 *
 * Returns the width of tabulation in characters.
 *
 * Return value: width of tab.
 **/
guint
gtk_source_view_get_tabs_width (GtkSourceView *view)
{
	g_return_val_if_fail (view != NULL, FALSE);
	g_return_val_if_fail (GTK_IS_SOURCE_VIEW (view), FALSE);

	return view->priv->tabs_width;
}

static gboolean
set_tab_stops_internal (GtkSourceView *view)
{
	PangoTabArray *tab_array;
	gint real_tab_width;

	real_tab_width = calculate_real_tab_width (view, view->priv->tabs_width, ' ');

	if (real_tab_width < 0)
		return FALSE;
	
	tab_array = pango_tab_array_new (1, TRUE);
	pango_tab_array_set_tab (tab_array, 0, PANGO_TAB_LEFT, real_tab_width);

	gtk_text_view_set_tabs (GTK_TEXT_VIEW (view), 
				tab_array);

	pango_tab_array_free (tab_array);

	return TRUE;
}

/**
 * gtk_source_view_set_tabs_width:
 * @view: a #GtkSourceView.
 * @width: width of tab in characters.
 *
 * Sets the width of tabulation in characters.
 *
 **/
void
gtk_source_view_set_tabs_width (GtkSourceView *view,
				guint          width)
{
	guint save_width;
	
	g_return_if_fail (GTK_SOURCE_VIEW (view));
	g_return_if_fail (width <= MAX_TAB_WIDTH);
	g_return_if_fail (width > 0);

	if (view->priv->tabs_width == width)
		return;
	
	gtk_widget_ensure_style (GTK_WIDGET (view));
	
	save_width = view->priv->tabs_width;
	view->priv->tabs_width = width;
	if (set_tab_stops_internal (view))
	{
		g_object_notify (G_OBJECT (view), "tabs_width");
	}
	else
	{
		g_warning ("Impossible to set tabs width.");
		view->priv->tabs_width = save_width;
	}
}

static gchar *
compute_indentation (GtkSourceView *view, 
		     GtkTextIter   *cur)
{
	GtkTextIter start;
	GtkTextIter end;

	gunichar ch;
	gint line;

	line = gtk_text_iter_get_line (cur);

	gtk_text_buffer_get_iter_at_line (gtk_text_view_get_buffer (GTK_TEXT_VIEW (view)), 
					  &start, 
					  line);

	end = start;

	ch = gtk_text_iter_get_char (&end);

	while (g_unichar_isspace (ch) && 
	       (ch != '\n') && 
	       (ch != '\r') && 
	       (gtk_text_iter_compare (&end, cur) < 0))
	{
		if (!gtk_text_iter_forward_char (&end))
			break;

		ch = gtk_text_iter_get_char (&end);
	}

	if (gtk_text_iter_equal (&start, &end))
		return NULL;

	return gtk_text_iter_get_slice (&start, &end);
}

static void
indent_lines (GtkSourceView *view, GtkTextIter *start, GtkTextIter *end)
{
	GtkTextBuffer *buf;
	gint start_line, end_line;
	gint i;
	gchar *tab_buffer = NULL;

	buf = gtk_text_view_get_buffer (GTK_TEXT_VIEW (view));

	start_line = gtk_text_iter_get_line (start);
	end_line = gtk_text_iter_get_line (end);

	if ((gtk_text_iter_get_visible_line_offset (end) == 0) &&
	    (end_line > start_line))
	{
		end_line--;
	}

	if (gtk_source_view_get_insert_spaces_instead_of_tabs (view))
	{
		gint tabs_size;

		tabs_size = view->priv->tabs_width;
		tab_buffer = g_strnfill (tabs_size, ' ');
	}
	else
	{
		tab_buffer = g_strdup ("\t");
	}

	gtk_text_buffer_begin_user_action (buf);

	for (i = start_line; i <= end_line; i++)
	{
		GtkTextIter iter;

		gtk_text_buffer_get_iter_at_line (buf, &iter, i);

		/* don't add indentation on empty lines */
		if (gtk_text_iter_ends_line (&iter))
			continue;

		gtk_text_buffer_insert (buf, &iter, tab_buffer, -1);
	}

	gtk_text_buffer_end_user_action (buf);

	g_free (tab_buffer);

	gtk_text_view_scroll_mark_onscreen (GTK_TEXT_VIEW (view),
					    gtk_text_buffer_get_insert (buf));
}

static void
unindent_lines (GtkSourceView *view, GtkTextIter *start, GtkTextIter *end)
{
	GtkTextBuffer *buf;
	gint start_line, end_line;
	gint i;

	buf = gtk_text_view_get_buffer (GTK_TEXT_VIEW (view));

	start_line = gtk_text_iter_get_line (start);
	end_line = gtk_text_iter_get_line (end);

	if ((gtk_text_iter_get_visible_line_offset (end) == 0) &&
	    (end_line > start_line))
	{
		end_line--;
	}

	gtk_text_buffer_begin_user_action (buf);

	for (i = start_line; i <= end_line; i++)
	{
		GtkTextIter iter, iter2;

		gtk_text_buffer_get_iter_at_line (buf, &iter, i);

		if (gtk_text_iter_get_char (&iter) == '\t')
		{
			iter2 = iter;
			gtk_text_iter_forward_char (&iter2);
			gtk_text_buffer_delete (buf, &iter, &iter2);
		}
		else if (gtk_text_iter_get_char (&iter) == ' ')
		{
			gint spaces = 0;

			iter2 = iter;

			while (!gtk_text_iter_ends_line (&iter2))
			{
				if (gtk_text_iter_get_char (&iter2) == ' ')
					spaces++;
				else
					break;

				gtk_text_iter_forward_char (&iter2);
			}

			if (spaces > 0)
			{
				gint tabs = 0;
				gint tabs_size = view->priv->tabs_width;

				tabs = spaces / tabs_size;
				spaces = spaces - (tabs * tabs_size);

				if (spaces == 0)
					spaces = tabs_size;

				iter2 = iter;

				gtk_text_iter_forward_chars (&iter2, spaces);
				gtk_text_buffer_delete (buf, &iter, &iter2);
			}
		}
	}

	gtk_text_buffer_end_user_action (buf);

	gtk_text_view_scroll_mark_onscreen (GTK_TEXT_VIEW (view),
					    gtk_text_buffer_get_insert (buf));
}

static void
insert_tab_or_spaces (GtkSourceView *view, GtkTextIter *start, GtkTextIter *end)
{
	GtkTextBuffer *buf;
	gchar *tab_buf;

	if (view->priv->insert_spaces)
	{
		gint tabs_size;
		GtkTextIter iter;
		gint cur_pos;
		gint tab_pos;
		gint num_of_equivalent_spaces;

		tabs_size = view->priv->tabs_width; 

		iter = *start;

		cur_pos = gtk_text_iter_get_line_offset (start);
		tab_pos = cur_pos;

		/* If there are tabs in the line take them into account */
		while (tab_pos > 0)
		{
			gunichar c;

			gtk_text_iter_backward_char (&iter);
			c = gtk_text_iter_get_char (&iter);
			if (c == '\t')
				break;
			tab_pos--;
		}

		num_of_equivalent_spaces = tabs_size - (cur_pos - tab_pos) % tabs_size; 

		tab_buf = g_strnfill (num_of_equivalent_spaces, ' ');
	}
	else
	{
		tab_buf = g_strdup ("\t");
	}

	buf = gtk_text_view_get_buffer (GTK_TEXT_VIEW (view));

	gtk_text_buffer_begin_user_action (buf);
	gtk_text_buffer_delete (buf, start, end);
	gtk_text_buffer_insert (buf, start, tab_buf, -1);
	gtk_text_buffer_end_user_action (buf);

	g_free (tab_buf);
}

static gboolean
gtk_source_view_key_press_event (GtkWidget *widget, GdkEventKey *event)
{
	GtkSourceView *view;
	GtkTextBuffer *buf;
	GtkTextIter cur;
	GtkTextMark *mark;
	guint modifiers;
	gint key;

	view = GTK_SOURCE_VIEW (widget);
	buf = gtk_text_view_get_buffer (GTK_TEXT_VIEW (widget));

	/* Be careful when testing for modifier state equality:
	 * caps lock, num lock,etc need to be taken into account */
	modifiers = gtk_accelerator_get_default_mod_mask ();

	key = event->keyval;

	mark = gtk_text_buffer_get_insert (buf);
	gtk_text_buffer_get_iter_at_mark (buf, &cur, mark);

	if ((key == GDK_Return || key == GDK_KP_Enter) &&
	    !(event->state & GDK_SHIFT_MASK) &&
	    view->priv->auto_indent)
	{
		/* Auto-indent means that when you press ENTER at the end of a
		 * line, the new line is automatically indented at the same
		 * level as the previous line.
		 * SHIFT+ENTER allows to avoid autoindentation.
		 */
		gchar *indent = NULL;

		/* Calculate line indentation and create indent string. */
		indent = compute_indentation (view, &cur);

		if (indent != NULL)
		{
			/* Allow input methods to internally handle a key press event. 
			 * If this function returns TRUE, then no further processing should be done 
			 * for this keystroke. */
			if (gtk_im_context_filter_keypress (GTK_TEXT_VIEW(view)->im_context, event))
				return TRUE;

			/* If an input method has inserted some test while handling the key press event,
			 * the cur iterm may be invalid, so get the iter again */
			gtk_text_buffer_get_iter_at_mark (buf, &cur, mark);
	
			/* Insert new line and auto-indent. */
			gtk_text_buffer_begin_user_action (buf);
			gtk_text_buffer_insert (buf, &cur, "\n", 1);
			gtk_text_buffer_insert (buf, &cur, indent, strlen (indent));
			g_free (indent);
			gtk_text_buffer_end_user_action (buf);
			gtk_text_view_scroll_mark_onscreen (GTK_TEXT_VIEW (widget),
							    mark);
			return TRUE;
		}
	}

	/* if tab or shift+tab:
	 * with shift+tab key is GDK_ISO_Left_Tab (depends on X?)
	 */
	if ((key == GDK_Tab || key == GDK_KP_Tab || key == GDK_ISO_Left_Tab) &&
	    ((event->state & modifiers) == 0 ||
	     (event->state & modifiers) == GDK_SHIFT_MASK))
	{
		GtkTextIter s, e;
		gboolean has_selection;

		has_selection = gtk_text_buffer_get_selection_bounds (buf, &s, &e);

		if (view->priv->indent_on_tab)
		{
			/* shift+tab: always unindent */
			if (event->state & GDK_SHIFT_MASK)
			{
				unindent_lines (view, &s, &e);
				return TRUE;
			}

			/* tab: if we have a selection which spans one whole line
			 * or more, we mass indent, if the selection spans less then
			 * the full line just replace the text with \t
			 */
			if (has_selection &&
			    ((gtk_text_iter_starts_line (&s) && gtk_text_iter_ends_line (&e)) ||
			     (gtk_text_iter_get_line (&s) != gtk_text_iter_get_line (&e))))
			{
				indent_lines (view, &s, &e);
				return TRUE;
			}
		}
 
		insert_tab_or_spaces (view, &s, &e);
 		return TRUE;
	}

	return (* GTK_WIDGET_CLASS (parent_class)->key_press_event) (widget, event);
}

static void
extend_selection_to_line (GtkTextBuffer *buf, GtkTextIter *line_start)
{
	GtkTextIter start;
	GtkTextIter end;
	GtkTextIter line_end;

	gtk_text_buffer_get_selection_bounds (buf, &start, &end);

	line_end = *line_start;
	gtk_text_iter_forward_to_line_end (&line_end);

	if (gtk_text_iter_compare (&start, line_start) < 0)
	{
		gtk_text_buffer_select_range (buf, &start, &line_end);
	}
	else if (gtk_text_iter_compare (&end, &line_end) < 0)
	{
		/* if the selection is in this line, extend
		 * the selection to the whole line */
		gtk_text_buffer_select_range (buf, &line_end, line_start);
	}
	else
	{
		gtk_text_buffer_select_range (buf, &end, line_start);
	}
}

static void
select_line (GtkTextBuffer *buf, GtkTextIter *line_start)
{
	GtkTextIter iter;

	iter = *line_start;

	if (!gtk_text_iter_ends_line (&iter))
		gtk_text_iter_forward_to_line_end (&iter);

	/* Select the line, put the cursor at the end of the line */
	gtk_text_buffer_select_range (buf, &iter, line_start);
}

static gboolean
gtk_source_view_button_press_event (GtkWidget *widget, GdkEventButton *event)
{
	GtkSourceView *view;
	GtkTextBuffer *buf;
	int y_buf;
	GtkTextIter line_start;

	view = GTK_SOURCE_VIEW (widget);
	buf = gtk_text_view_get_buffer (GTK_TEXT_VIEW (widget));

	if (view->priv->show_line_numbers &&
	    (event->window == gtk_text_view_get_window (GTK_TEXT_VIEW (view),
						       GTK_TEXT_WINDOW_LEFT)))
	{
		gtk_text_view_window_to_buffer_coords (GTK_TEXT_VIEW (view),
						       GTK_TEXT_WINDOW_LEFT,
						       event->x, event->y,
						       NULL, &y_buf);

		gtk_text_view_get_line_at_y (GTK_TEXT_VIEW (view),
					     &line_start,
					     y_buf,
					     NULL);
	
		if (event->type == GDK_BUTTON_PRESS && (event->button == 1))
		{
			if ((event->state & GDK_CONTROL_MASK) != 0)
			{
				/* Single click + Ctrl -> select the line */
				select_line (buf, &line_start);
			}
			else if ((event->state & GDK_SHIFT_MASK) != 0)
			{
				/* Single click + Shift -> extended current
				   selection to include the clicked line */
				extend_selection_to_line (buf, &line_start);
			}
			else
			{
				gtk_text_buffer_place_cursor (buf, &line_start);
			}
		}
		else if (event->type == GDK_2BUTTON_PRESS && (event->button == 1))
		{
			select_line (buf, &line_start);
		}

		/* consume the event also on right click etc */
		return TRUE;
	}

	return GTK_WIDGET_CLASS (parent_class)->button_press_event (widget, event);
}

/**
 * gtk_source_view_get_auto_indent:
 * @view: a #GtkSourceView.
 *
 * Returns whether auto indentation of text is enabled.
 *
 * Return value: %TRUE if auto indentation is enabled.
 **/
gboolean
gtk_source_view_get_auto_indent (GtkSourceView *view)
{
	g_return_val_if_fail (GTK_IS_SOURCE_VIEW (view), FALSE);

	return view->priv->auto_indent;
}

/**
 * gtk_source_view_set_auto_indent:
 * @view: a #GtkSourceView.
 * @enable: whether to enable auto indentation.
 *
 * If %TRUE auto indentation of text is enabled.
 *
 **/
void
gtk_source_view_set_auto_indent (GtkSourceView *view, gboolean enable)
{
	g_return_if_fail (GTK_IS_SOURCE_VIEW (view));

	enable = (enable != FALSE);

	if (view->priv->auto_indent == enable)
		return;

	view->priv->auto_indent = enable;

	g_object_notify (G_OBJECT (view), "auto_indent");
}

/**
 * gtk_source_view_get_insert_spaces_instead_of_tabs:
 * @view: a #GtkSourceView.
 *
 * Returns whether when inserting a tabulator character it should
 * be replaced by a group of space characters.
 *
 * Return value: %TRUE if spaces are inserted instead of tabs.
 **/
gboolean
gtk_source_view_get_insert_spaces_instead_of_tabs (GtkSourceView *view)
{
	g_return_val_if_fail (GTK_IS_SOURCE_VIEW (view), FALSE);

	return view->priv->insert_spaces;
}

/**
 * gtk_source_view_set_insert_spaces_instead_of_tabs:
 * @view: a #GtkSourceView.
 * @enable: whether to insert spaces instead of tabs.
 *
 * If %TRUE any tabulator character inserted is replaced by a group
 * of space characters.
 *
 **/
void
gtk_source_view_set_insert_spaces_instead_of_tabs (GtkSourceView *view, gboolean enable)
{
	g_return_if_fail (GTK_IS_SOURCE_VIEW (view));

	enable = (enable != FALSE);
	
	if (view->priv->insert_spaces == enable)
		return;
	
	view->priv->insert_spaces = enable;

	g_object_notify (G_OBJECT (view), "insert_spaces_instead_of_tabs");
}

/**
 * gtk_source_view_get_indent_on_tab:
 * @view: a #GtkSourceView.
 *
 * Returns whether when the tab key is pressed the current selection
 * should get indented instead of replaced with the \t character.
 *
 * Return value: %TRUE if the selection is indented when tab is pressed.
 *
 * Since: 1.8
 **/
gboolean
gtk_source_view_get_indent_on_tab (GtkSourceView *view)
{
	g_return_val_if_fail (GTK_IS_SOURCE_VIEW (view), FALSE);

	return view->priv->indent_on_tab;
}

/**
 * gtk_source_view_set_indent_on_tab:
 * @view: a #GtkSourceView.
 * @enable: whether to indent a block when tab is pressed.
 *
 * If %TRUE, when the tab key is pressed and there is a selection, the
 * selected text is indented of one level instead of being replaced with
 * the \t characters. Shift+Tab unindents the selection.
 *
 * Since: 1.8
 **/
void
gtk_source_view_set_indent_on_tab (GtkSourceView *view, gboolean enable)
{
	g_return_if_fail (GTK_IS_SOURCE_VIEW (view));

	enable = (enable != FALSE);

	if (view->priv->indent_on_tab == enable)
		return;

	view->priv->indent_on_tab = enable;

	g_object_notify (G_OBJECT (view), "indent_on_tab");
}

static void 
view_dnd_drop (GtkTextView *view, 
	       GdkDragContext *context,
	       gint x,
	       gint y,
	       GtkSelectionData *selection_data,
	       guint info,
	       guint time,
	       gpointer data)
{

	GtkTextIter iter;

	if (info == TARGET_COLOR) 
	{
		guint16 *vals;
		gchar string[] = "#000000";
		gint buffer_x;
		gint buffer_y;
		
		if (selection_data->length < 0)
			return;

		if ((selection_data->format != 16) || (selection_data->length != 8)) 
		{
			g_warning ("Received invalid color data\n");
			return;
		}

		vals = (guint16 *) selection_data->data;

		vals[0] /= 256;
	        vals[1] /= 256;
		vals[2] /= 256;
		
		g_snprintf (string, sizeof (string), "#%02X%02X%02X", vals[0], vals[1], vals[2]);
		
		gtk_text_view_window_to_buffer_coords (view, 
						       GTK_TEXT_WINDOW_TEXT, 
						       x, 
						       y, 
						       &buffer_x, 
						       &buffer_y);
		gtk_text_view_get_iter_at_location (view, &iter, buffer_x, buffer_y);

		if (gtk_text_view_get_editable (view))
		{
			gtk_text_buffer_insert (gtk_text_view_get_buffer (view), 
						&iter, 
						string, 
						strlen (string));
			gtk_text_buffer_place_cursor (gtk_text_view_get_buffer (view),
						&iter);
		}

		/*
		 * FIXME: Check if the iter is inside a selection
		 * If it is, remove the selection and then insert at
		 * the cursor position - Paolo
		 */

		return;
	}
}

/**
 * gtk_source_view_set_smart_home_end:
 * @view: a #GtkSourceView.
 * @enable: whether to enable smart behavior for HOME and END keys.
 *
 * If %TRUE HOME and END keys will move to the first/last non-space
 * character of the line before moving to the start/end.
 *
 **/
void
gtk_source_view_set_smart_home_end (GtkSourceView *view, gboolean enable)
{
	g_return_if_fail (GTK_IS_SOURCE_VIEW (view));

	enable = (enable != FALSE);

	if (view->priv->smart_home_end == enable)
		return;

	view->priv->smart_home_end = enable;

	g_object_notify (G_OBJECT (view), "smart_home_end");
}

/**
 * gtk_source_view_get_smart_home_end:
 * @view: a #GtkSourceView.
 *
 * Returns whether HOME and END keys will move to the first/last non-space
 * character of the line before moving to the start/end.
 *
 * Return value: %TRUE if smart behavior for HOME and END keys is enabled.
 **/
gboolean
gtk_source_view_get_smart_home_end (GtkSourceView *view)
{
	g_return_val_if_fail (GTK_IS_SOURCE_VIEW (view), FALSE);

	return view->priv->smart_home_end;
}

/**
 * gtk_source_view_style_set:
 * @widget: a #GtkSourceView.
 * @previous_style:
 *
 *
 **/
static void 
gtk_source_view_style_set (GtkWidget *widget, GtkStyle *previous_style)
{
	GtkSourceView *view;
	
	g_return_if_fail (GTK_IS_SOURCE_VIEW (widget));
	
	/* call default handler first */
	if (GTK_WIDGET_CLASS (parent_class)->style_set)
		(* GTK_WIDGET_CLASS (parent_class)->style_set) (widget, previous_style);

	view = GTK_SOURCE_VIEW (widget);
	if (previous_style)
	{
		/* If previous_style is NULL this is the initial
		 * emission and we can't set the tab array since the
		 * text view doesn't have a default style yet */
		
		/* re-set tab stops */
		set_tab_stops_internal (view);
	}
}

