/* GTK - The GIMP Toolkit
 * Copyright (C) 1995-1997 Peter Mattis, Spencer Kimball and Josh MacDonald
 *
 * GtkSpinButton widget for GTK+
 * Copyright (C) 1998 Lars Hamann and Stefan Jeske
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

/*
 * Modified by the GTK+ Team and others 1997-1999.  See the AUTHORS
 * file for a list of people on the GTK+ Team.  See the ChangeLog
 * files for a list of changes.  These files are distributed with
 * GTK+ at ftp://ftp.gtk.org/pub/gtk/. 
 */

#include "gretl.h"

#include "obsbutton.h"
#include <gtk/gtkmain.h>

#define MIN_OBS_BUTTON_WIDTH              30
#define ARROW_SIZE                         11
#define OBS_BUTTON_INITIAL_TIMER_DELAY    200
#define OBS_BUTTON_TIMER_DELAY            20
#define MAX_TEXT_LENGTH                    256
#define MAX_TIMER_CALLS                    5
#define EPSILON                            1e-5

enum {
  ARG_0,
  ARG_SHADOW_TYPE,
  ARG_VALUE
};


static void obs_button_class_init     (ObsButtonClass *klass);
static void obs_button_init           (ObsButton      *obs_button);
static void obs_button_finalize       (GtkObject          *object);
static void obs_button_set_arg        (GtkObject          *object,
				       GtkArg             *arg,
				       guint               arg_id);
static void obs_button_get_arg        (GtkObject          *object,
				       GtkArg             *arg,
				       guint               arg_id);
static void obs_button_map            (GtkWidget          *widget);
static void obs_button_unmap          (GtkWidget          *widget);
static void obs_button_realize        (GtkWidget          *widget);
static void obs_button_unrealize      (GtkWidget          *widget);
static void obs_button_size_request   (GtkWidget          *widget,
				       GtkRequisition     *requisition);
static void obs_button_size_allocate  (GtkWidget          *widget,
				       GtkAllocation      *allocation);
static void obs_button_paint          (GtkWidget          *widget,
				       GdkRectangle       *area);
static void obs_button_draw           (GtkWidget          *widget,
				       GdkRectangle       *area);
static gint obs_button_expose         (GtkWidget          *widget,
				       GdkEventExpose     *event);
static gint obs_button_button_press   (GtkWidget          *widget,
				       GdkEventButton     *event);
static gint obs_button_button_release (GtkWidget          *widget,
				       GdkEventButton     *event);
static gint obs_button_motion_notify  (GtkWidget          *widget,
				       GdkEventMotion     *event);
static gint obs_button_enter_notify   (GtkWidget          *widget,
				       GdkEventCrossing   *event);
static gint obs_button_leave_notify   (GtkWidget          *widget,
				       GdkEventCrossing   *event);
static gint obs_button_focus_out      (GtkWidget          *widget,
				       GdkEventFocus      *event);
static void obs_button_draw_arrow     (ObsButton      *obs_button, 
				       guint               arrow);
static gint obs_button_timer          (ObsButton      *obs_button);
static void obs_button_value_changed  (GtkAdjustment      *adjustment,
				       ObsButton      *obs_button); 
static gint obs_button_key_press      (GtkWidget          *widget,
				       GdkEventKey        *event);
static gint obs_button_key_release    (GtkWidget          *widget,
				       GdkEventKey        *event);
static void obs_button_activate       (GtkEditable        *editable);
static void obs_button_insert_text    (GtkEditable        *editable,
				       const gchar        *new_text,
				       gint                new_text_length,
				       gint               *position);
static void obs_button_real_spin      (ObsButton      *obs_button,
				       gfloat              step);
static void obs_button_set_shadow_type (ObsButton *obs_button,
					GtkShadowType  shadow_type);
static void obs_button_set_adjustment  (ObsButton *obs_button,
					GtkAdjustment *adjustment);


static GtkEntryClass *parent_class = NULL;


GtkType
obs_button_get_type (void)
{
  static guint obs_button_type = 0;

  if (!obs_button_type)
    {
      static const GtkTypeInfo obs_button_info =
      {
	"ObsButton",
	sizeof (ObsButton),
	sizeof (ObsButtonClass),
	(GtkClassInitFunc) obs_button_class_init,
	(GtkObjectInitFunc) obs_button_init,
	/* reserved_1 */ NULL,
        /* reserved_2 */ NULL,
        (GtkClassInitFunc) NULL,
      };

      obs_button_type = gtk_type_unique (GTK_TYPE_ENTRY, &obs_button_info);
    }
  return obs_button_type;
}

static void
obs_button_class_init (ObsButtonClass *class)
{
  GtkObjectClass   *object_class;
  GtkWidgetClass   *widget_class;
  GtkEditableClass *editable_class;

  object_class   = (GtkObjectClass*)   class;
  widget_class   = (GtkWidgetClass*)   class;
  editable_class = (GtkEditableClass*) class; 

  parent_class = gtk_type_class (GTK_TYPE_ENTRY);

  gtk_object_add_arg_type ("ObsButton::shadow_type",
			   GTK_TYPE_SHADOW_TYPE,
			   GTK_ARG_READWRITE,
			   ARG_SHADOW_TYPE);
  gtk_object_add_arg_type ("ObsButton::value",
			   GTK_TYPE_FLOAT,
			   GTK_ARG_READWRITE,
			   ARG_VALUE);
  
  object_class->set_arg = obs_button_set_arg;
  object_class->get_arg = obs_button_get_arg;
  object_class->finalize = obs_button_finalize;

  widget_class->map = obs_button_map;
  widget_class->unmap = obs_button_unmap;
  widget_class->realize = obs_button_realize;
  widget_class->unrealize = obs_button_unrealize;
  widget_class->size_request = obs_button_size_request;
  widget_class->size_allocate = obs_button_size_allocate;
  widget_class->draw = obs_button_draw;
  widget_class->expose_event = obs_button_expose;
  widget_class->button_press_event = obs_button_button_press;
  widget_class->button_release_event = obs_button_button_release;
  widget_class->motion_notify_event = obs_button_motion_notify;
  widget_class->key_press_event = obs_button_key_press;
  widget_class->key_release_event = obs_button_key_release;
  widget_class->enter_notify_event = obs_button_enter_notify;
  widget_class->leave_notify_event = obs_button_leave_notify;
  widget_class->focus_out_event = obs_button_focus_out;

  editable_class->insert_text = obs_button_insert_text;
  editable_class->activate = obs_button_activate;
}

static void
obs_button_set_arg (GtkObject        *object,
		    GtkArg           *arg,
		    guint             arg_id)
{
    ObsButton *obs_button;

    obs_button = OBS_BUTTON (object);
  
    switch (arg_id)
	{
	case ARG_SHADOW_TYPE:
	    obs_button_set_shadow_type (obs_button, GTK_VALUE_ENUM (*arg));
	    break;
	case ARG_VALUE:
	    obs_button_set_value (obs_button, GTK_VALUE_FLOAT (*arg));
	    break;
	default:
	    break;
	}
}

static void
obs_button_get_arg (GtkObject        *object,
		    GtkArg           *arg,
		    guint             arg_id)
{
    ObsButton *obs_button;

    obs_button = OBS_BUTTON (object);
  
    switch (arg_id)
	{
	case ARG_SHADOW_TYPE:
	    GTK_VALUE_ENUM (*arg) = obs_button->shadow_type;
	    break;
	case ARG_VALUE:
	    GTK_VALUE_FLOAT (*arg) = obs_button->adjustment->value;
	    break;
	default:
	    arg->type = GTK_TYPE_INVALID;
	    break;
	}
}

static void
obs_button_init (ObsButton *obs_button)
{
  obs_button->adjustment = NULL;
  obs_button->panel = NULL;
  obs_button->shadow_type = GTK_SHADOW_NONE;
  obs_button->timer = 0;
  obs_button->ev_time = 0;
  obs_button->timer_step = 0.0;
  obs_button->in_child = 2;
  obs_button->click_child = 2;
  obs_button->button = 0;
  obs_button->need_timer = FALSE;
  obs_button->timer_calls = 0;
  obs_button->adjustment = NULL;
}

static void
obs_button_finalize (GtkObject *object)
{
  g_return_if_fail (object != NULL);
  g_return_if_fail (GTK_IS_OBS_BUTTON (object));

  gtk_object_unref (GTK_OBJECT (OBS_BUTTON (object)->adjustment));
  
  GTK_OBJECT_CLASS (parent_class)->finalize (object);
}

static void
obs_button_map (GtkWidget *widget)
{
  g_return_if_fail (widget != NULL);
  g_return_if_fail (GTK_IS_OBS_BUTTON (widget));

  if (GTK_WIDGET_REALIZED (widget) && !GTK_WIDGET_MAPPED (widget))
    {
      GTK_WIDGET_CLASS (parent_class)->map (widget);
      gdk_window_show (OBS_BUTTON (widget)->panel);
    }
}

static void
obs_button_unmap (GtkWidget *widget)
{
  g_return_if_fail (widget != NULL);
  g_return_if_fail (GTK_IS_OBS_BUTTON (widget));

  if (GTK_WIDGET_MAPPED (widget))
    {
      gdk_window_hide (OBS_BUTTON (widget)->panel);
      GTK_WIDGET_CLASS (parent_class)->unmap (widget);
    }
}

static void
obs_button_realize (GtkWidget *widget)
{
  ObsButton *spin;
  GdkWindowAttr attributes;
  gint attributes_mask;
  guint real_width;

  g_return_if_fail (widget != NULL);
  g_return_if_fail (GTK_IS_OBS_BUTTON (widget));
  
  spin = OBS_BUTTON (widget);

  real_width = widget->allocation.width;
  widget->allocation.width -= ARROW_SIZE + 2 * widget->style->klass->xthickness;
  gtk_widget_set_events (widget, gtk_widget_get_events (widget) |
			 GDK_KEY_RELEASE_MASK);
  GTK_WIDGET_CLASS (parent_class)->realize (widget);

  widget->allocation.width = real_width;
  
  attributes.window_type = GDK_WINDOW_CHILD;
  attributes.wclass = GDK_INPUT_OUTPUT;
  attributes.visual = gtk_widget_get_visual (widget);
  attributes.colormap = gtk_widget_get_colormap (widget);
  attributes.event_mask = gtk_widget_get_events (widget);
  attributes.event_mask |= GDK_EXPOSURE_MASK | GDK_BUTTON_PRESS_MASK 
    | GDK_BUTTON_RELEASE_MASK | GDK_LEAVE_NOTIFY_MASK | GDK_ENTER_NOTIFY_MASK 
    | GDK_POINTER_MOTION_MASK | GDK_POINTER_MOTION_HINT_MASK;

  attributes_mask = GDK_WA_X | GDK_WA_Y | GDK_WA_VISUAL | GDK_WA_COLORMAP;

  attributes.x = (widget->allocation.x + widget->allocation.width - ARROW_SIZE -
		  2 * widget->style->klass->xthickness);
  attributes.y = widget->allocation.y + (widget->allocation.height -
					 widget->requisition.height) / 2;
  attributes.width = ARROW_SIZE + 2 * widget->style->klass->xthickness;
  attributes.height = widget->requisition.height;
  
  spin->panel = gdk_window_new (gtk_widget_get_parent_window (widget), 
				&attributes, attributes_mask);
  gdk_window_set_user_data (spin->panel, widget);

  gtk_style_set_background (widget->style, spin->panel, GTK_STATE_NORMAL);
}

static void
obs_button_unrealize (GtkWidget *widget)
{
  ObsButton *spin;

  g_return_if_fail (widget != NULL);
  g_return_if_fail (GTK_IS_OBS_BUTTON (widget));

  spin = OBS_BUTTON (widget);

  GTK_WIDGET_CLASS (parent_class)->unrealize (widget);

  if (spin->panel)
    {
      gdk_window_set_user_data (spin->panel, NULL);
      gdk_window_destroy (spin->panel);
      spin->panel = NULL;
    }
}

static void
obs_button_size_request (GtkWidget      *widget,
			 GtkRequisition *requisition)
{
    int cw, w, width, max_string_len, string_len;

    g_return_if_fail (widget != NULL);
    g_return_if_fail (requisition != NULL);
    g_return_if_fail (GTK_IS_OBS_BUTTON (widget));

    GTK_WIDGET_CLASS (parent_class)->size_request (widget, requisition);

    cw = gdk_char_width(fixed_font, 'x');

    width = MIN_OBS_BUTTON_WIDTH;
    max_string_len = OBSLEN;
    string_len = strlen(datainfo->endobs) + 1;
    w = MIN (string_len, max_string_len) * cw;
    width = MAX (width, w);  
  
    requisition->width = width + ARROW_SIZE 
	+ 2 * widget->style->klass->xthickness;
}

static void
obs_button_size_allocate (GtkWidget     *widget,
			       GtkAllocation *allocation)
{
  GtkAllocation child_allocation;

  g_return_if_fail (widget != NULL);
  g_return_if_fail (GTK_IS_OBS_BUTTON (widget));
  g_return_if_fail (allocation != NULL);

  child_allocation = *allocation;
  if (child_allocation.width > ARROW_SIZE + 2 * widget->style->klass->xthickness)
    child_allocation.width -= ARROW_SIZE + 2 * widget->style->klass->xthickness;

  GTK_WIDGET_CLASS (parent_class)->size_allocate (widget, &child_allocation);

  widget->allocation = *allocation;

  if (GTK_WIDGET_REALIZED (widget))
    {
      child_allocation.width = ARROW_SIZE + 2 * widget->style->klass->xthickness;
      child_allocation.height = widget->requisition.height;  
      child_allocation.x = (allocation->x + allocation->width - ARROW_SIZE - 
			    2 * widget->style->klass->xthickness);
      child_allocation.y = allocation->y + (allocation->height - widget->requisition.height) / 2;

      gdk_window_move_resize (OBS_BUTTON (widget)->panel, 
			      child_allocation.x,
			      child_allocation.y,
			      child_allocation.width,
			      child_allocation.height); 
    }
}

static GtkShadowType
obs_button_get_shadow_type (ObsButton *obs_button)
{
  GtkWidget *widget = GTK_WIDGET (obs_button);
  
  GtkShadowType shadow_type =
    gtk_style_get_prop_experimental (widget->style,
				     "ObsButton::shadow_type", -1);

  if (shadow_type != (GtkShadowType)-1)
    return shadow_type;
  else
    return obs_button->shadow_type;
}

static void
obs_button_paint (GtkWidget    *widget,
		       GdkRectangle *area)
{
  ObsButton *spin;
  GtkShadowType shadow_type;

  g_return_if_fail (widget != NULL);
  g_return_if_fail (GTK_IS_OBS_BUTTON (widget));

  spin = OBS_BUTTON (widget);
  shadow_type = obs_button_get_shadow_type (spin);

  if (GTK_WIDGET_DRAWABLE (widget))
    {
      if (shadow_type != GTK_SHADOW_NONE)
	gtk_paint_box (widget->style, spin->panel,
		       GTK_STATE_NORMAL, shadow_type,
		       area, widget, "spinbutton",
		       0, 0, 
		       ARROW_SIZE + 2 * widget->style->klass->xthickness,
		       widget->requisition.height); 
      else
	{
	  gdk_window_set_back_pixmap (spin->panel, NULL, TRUE);
	  gdk_window_clear_area (spin->panel, area->x, area->y, area->width, area->height);
	}
      obs_button_draw_arrow (spin, GTK_ARROW_UP);
      obs_button_draw_arrow (spin, GTK_ARROW_DOWN);
      
      GTK_WIDGET_CLASS (parent_class)->draw (widget, area);
    }
}

static void
obs_button_draw (GtkWidget    *widget,
		      GdkRectangle *area)
{
  g_return_if_fail (widget != NULL);
  g_return_if_fail (GTK_IS_OBS_BUTTON (widget));
  g_return_if_fail (area != NULL);

  if (GTK_WIDGET_DRAWABLE (widget))
    obs_button_paint (widget, area);
}

static gint
obs_button_expose (GtkWidget      *widget,
			GdkEventExpose *event)
{
  g_return_val_if_fail (widget != NULL, FALSE);
  g_return_val_if_fail (GTK_IS_OBS_BUTTON (widget), FALSE);
  g_return_val_if_fail (event != NULL, FALSE);

  if (GTK_WIDGET_DRAWABLE (widget))
    obs_button_paint (widget, &event->area);

  return FALSE;
}

static void
obs_button_draw_arrow (ObsButton *obs_button, 
			    guint          arrow)
{
  GtkStateType state_type;
  GtkShadowType shadow_type;
  GtkShadowType spin_shadow_type;
  GtkWidget *widget;
  gint x;
  gint y;

  g_return_if_fail (obs_button != NULL);
  g_return_if_fail (GTK_IS_OBS_BUTTON (obs_button));
  
  widget = GTK_WIDGET (obs_button);

  spin_shadow_type = obs_button_get_shadow_type (obs_button);

  if (GTK_WIDGET_DRAWABLE (obs_button))
    {
      if (((arrow == GTK_ARROW_UP &&
	  (obs_button->adjustment->upper - obs_button->adjustment->value
	   <= EPSILON))) ||
	  ((arrow == GTK_ARROW_DOWN &&
	  (obs_button->adjustment->value - obs_button->adjustment->lower
	   <= EPSILON))))
	{
	  shadow_type = GTK_SHADOW_ETCHED_IN;
	  state_type = GTK_STATE_NORMAL;
	}
      else
	{
	  if (obs_button->in_child == arrow)
	    {
	      if (obs_button->click_child == arrow)
		state_type = GTK_STATE_ACTIVE;
	      else
		state_type = GTK_STATE_PRELIGHT;
	    }
	  else
	    state_type = GTK_STATE_NORMAL;
	  
	  if (obs_button->click_child == arrow)
	    shadow_type = GTK_SHADOW_IN;
	  else
	    shadow_type = GTK_SHADOW_OUT;
	}
      if (arrow == GTK_ARROW_UP)
	{
	  if (spin_shadow_type != GTK_SHADOW_NONE)
	    {
	      x = widget->style->klass->xthickness;
	      y = widget->style->klass->ythickness;
	    }
	  else
	    {
	      x = widget->style->klass->xthickness - 1;
	      y = widget->style->klass->ythickness - 1;
	    }
	  gtk_paint_arrow (widget->style, obs_button->panel,
			   state_type, shadow_type, 
			   NULL, widget, "spinbutton",
			   arrow, TRUE, 
			   x, y, ARROW_SIZE, widget->requisition.height / 2 
			   - widget->style->klass->ythickness);
	}
      else
	{
	  if (spin_shadow_type != GTK_SHADOW_NONE)
	    {
	      x = widget->style->klass->xthickness;
	      y = widget->requisition.height / 2;
	    }
	  else
	    {
	      x = widget->style->klass->xthickness - 1;
	      y = widget->requisition.height / 2 + 1;
	    }
	  gtk_paint_arrow (widget->style, obs_button->panel,
			   state_type, shadow_type, 
			   NULL, widget, "spinbutton",
			   arrow, TRUE, 
			   x, y, ARROW_SIZE, widget->requisition.height / 2 
			   - widget->style->klass->ythickness);
	}
    }
}

static gint
obs_button_enter_notify (GtkWidget        *widget,
			 GdkEventCrossing *event)
{
  ObsButton *spin;

  g_return_val_if_fail (widget != NULL, FALSE);
  g_return_val_if_fail (GTK_IS_OBS_BUTTON (widget), FALSE);
  g_return_val_if_fail (event != NULL, FALSE);

  spin = OBS_BUTTON (widget);

  if (event->window == spin->panel)
    {
      gint x;
      gint y;

      gdk_window_get_pointer (spin->panel, &x, &y, NULL);

      if (y <= widget->requisition.height / 2)
	{
	  spin->in_child = GTK_ARROW_UP;
	  if (spin->click_child == 2) 
	    obs_button_draw_arrow (spin, GTK_ARROW_UP);
	}
      else
	{
	  spin->in_child = GTK_ARROW_DOWN;
	  if (spin->click_child == 2) 
	    obs_button_draw_arrow (spin, GTK_ARROW_DOWN);
	}
    }
  return FALSE;
}

static gint
obs_button_leave_notify (GtkWidget        *widget,
			      GdkEventCrossing *event)
{
  ObsButton *spin;

  g_return_val_if_fail (widget != NULL, FALSE);
  g_return_val_if_fail (GTK_IS_OBS_BUTTON (widget), FALSE);
  g_return_val_if_fail (event != NULL, FALSE);

  spin = OBS_BUTTON (widget);

  if (event->window == spin->panel && spin->click_child == 2)
    {
      if (spin->in_child == GTK_ARROW_UP) 
	{
	  spin->in_child = 2;
	  obs_button_draw_arrow (spin, GTK_ARROW_UP);
	}
      else
	{
	  spin->in_child = 2;
	  obs_button_draw_arrow (spin, GTK_ARROW_DOWN);
	}
    }
  return FALSE;
}

static gint
obs_button_focus_out (GtkWidget     *widget,
			   GdkEventFocus *event)
{
  g_return_val_if_fail (widget != NULL, FALSE);
  g_return_val_if_fail (GTK_IS_OBS_BUTTON (widget), FALSE);
  g_return_val_if_fail (event != NULL, FALSE);

  if (GTK_EDITABLE (widget)->editable)
    obs_button_update (OBS_BUTTON (widget));

  return GTK_WIDGET_CLASS (parent_class)->focus_out_event (widget, event);
}

static gint
obs_button_button_press (GtkWidget      *widget,
			      GdkEventButton *event)
{
  ObsButton *spin;

  g_return_val_if_fail (widget != NULL, FALSE);
  g_return_val_if_fail (GTK_IS_OBS_BUTTON (widget), FALSE);
  g_return_val_if_fail (event != NULL, FALSE);

  spin = OBS_BUTTON (widget);

  if (!spin->button)
    {
      if (event->window == spin->panel)
	{
	  if (!GTK_WIDGET_HAS_FOCUS (widget))
	    gtk_widget_grab_focus (widget);
	  gtk_grab_add (widget);
	  spin->button = event->button;
	  
	  if (GTK_EDITABLE (widget)->editable)
	    obs_button_update (spin);
	  
	  if (event->y <= widget->requisition.height / 2)
	    {
	      spin->click_child = GTK_ARROW_UP;
	      if (event->button == 1)
		{
		 obs_button_real_spin (spin, 
					    spin->adjustment->step_increment);
		  if (!spin->timer)
		    {
		      spin->timer_step = spin->adjustment->step_increment;
		      spin->need_timer = TRUE;
		      spin->timer = gtk_timeout_add 
			(OBS_BUTTON_INITIAL_TIMER_DELAY, 
			 (GtkFunction) obs_button_timer, (gpointer) spin);
		    }
		}
	      else if (event->button == 2)
		{
		 obs_button_real_spin (spin, 
					    spin->adjustment->page_increment);
		  if (!spin->timer) 
		    {
		      spin->timer_step = spin->adjustment->page_increment;
		      spin->need_timer = TRUE;
		      spin->timer = gtk_timeout_add 
			(OBS_BUTTON_INITIAL_TIMER_DELAY, 
			 (GtkFunction) obs_button_timer, (gpointer) spin);
		    }
		}
	      obs_button_draw_arrow (spin, GTK_ARROW_UP);
	    }
	  else 
	    {
	      spin->click_child = GTK_ARROW_DOWN;
	      if (event->button == 1)
		{
		  obs_button_real_spin (spin,
					     -spin->adjustment->step_increment);
		  if (!spin->timer)
		    {
		      spin->timer_step = spin->adjustment->step_increment;
		      spin->need_timer = TRUE;
		      spin->timer = gtk_timeout_add 
			(OBS_BUTTON_INITIAL_TIMER_DELAY, 
			 (GtkFunction) obs_button_timer, (gpointer) spin);
		    }
		}      
	      else if (event->button == 2)
		{
		  obs_button_real_spin (spin,
					     -spin->adjustment->page_increment);
		  if (!spin->timer) 
		    {
		      spin->timer_step = spin->adjustment->page_increment;
		      spin->need_timer = TRUE;
		      spin->timer = gtk_timeout_add 
			(OBS_BUTTON_INITIAL_TIMER_DELAY, 
			 (GtkFunction) obs_button_timer, (gpointer) spin);
		    }
		}
	      obs_button_draw_arrow (spin, GTK_ARROW_DOWN);
	    }
	}
      else
	GTK_WIDGET_CLASS (parent_class)->button_press_event (widget, event);
    }
  return FALSE;
}

static gint
obs_button_button_release (GtkWidget      *widget,
				GdkEventButton *event)
{
  ObsButton *spin;

  g_return_val_if_fail (widget != NULL, FALSE);
  g_return_val_if_fail (GTK_IS_OBS_BUTTON (widget), FALSE);
  g_return_val_if_fail (event != NULL, FALSE);

  spin = OBS_BUTTON (widget);

  if (event->button == spin->button)
    {
      guint click_child;

      if (spin->timer)
	{
	  gtk_timeout_remove (spin->timer);
	  spin->timer = 0;
	  spin->timer_calls = 0;
	  spin->need_timer = FALSE;
	}

      if (event->button == 3)
	{
	  if (event->y >= 0 && event->x >= 0 && 
	      event->y <= widget->requisition.height &&
	      event->x <= ARROW_SIZE + 2 * widget->style->klass->xthickness)
	    {
	      if (spin->click_child == GTK_ARROW_UP &&
		  event->y <= widget->requisition.height / 2)
		{
		  gfloat diff;

		  diff = spin->adjustment->upper - spin->adjustment->value;
		  if (diff > EPSILON)
		    obs_button_real_spin (spin, diff);
		}
	      else if (spin->click_child == GTK_ARROW_DOWN &&
		       event->y > widget->requisition.height / 2)
		{
		  gfloat diff;

		  diff = spin->adjustment->value - spin->adjustment->lower;
		  if (diff > EPSILON)
		    obs_button_real_spin (spin, -diff);
		}
	    }
	}		  
      gtk_grab_remove (widget);
      click_child = spin->click_child;
      spin->click_child = 2;
      spin->button = 0;
      obs_button_draw_arrow (spin, click_child);
    }
  else
    GTK_WIDGET_CLASS (parent_class)->button_release_event (widget, event);

  return FALSE;
}

static gint
obs_button_motion_notify (GtkWidget      *widget,
			       GdkEventMotion *event)
{
  ObsButton *spin;

  g_return_val_if_fail (widget != NULL, FALSE);
  g_return_val_if_fail (GTK_IS_OBS_BUTTON (widget), FALSE);
  g_return_val_if_fail (event != NULL, FALSE);

  spin = OBS_BUTTON (widget);
  
  if (spin->button)
    return FALSE;

  if (event->window == spin->panel)
    {
      gint y;

      y = event->y;
      if (event->is_hint)
	gdk_window_get_pointer (spin->panel, NULL, &y, NULL);

      if (y <= widget->requisition.height / 2 && 
	  spin->in_child == GTK_ARROW_DOWN)
	{
	  spin->in_child = GTK_ARROW_UP;
	  obs_button_draw_arrow (spin, GTK_ARROW_UP);
	  obs_button_draw_arrow (spin, GTK_ARROW_DOWN);
	}
      else if (y > widget->requisition.height / 2 && 
	  spin->in_child == GTK_ARROW_UP)
	{
	  spin->in_child = GTK_ARROW_DOWN;
	  obs_button_draw_arrow (spin, GTK_ARROW_UP);
	  obs_button_draw_arrow (spin, GTK_ARROW_DOWN);
	}
      return FALSE;
    }
	  
  return GTK_WIDGET_CLASS (parent_class)->motion_notify_event (widget, event);
}

static gint
obs_button_timer (ObsButton *obs_button)
{
  gboolean retval = FALSE;
  
  GDK_THREADS_ENTER ();

  if (obs_button->timer)
    {
      if (obs_button->click_child == GTK_ARROW_UP)
	obs_button_real_spin (obs_button,	obs_button->timer_step);
      else
	obs_button_real_spin (obs_button,	-obs_button->timer_step);

      if (obs_button->need_timer)
	{
	  obs_button->need_timer = FALSE;
	  obs_button->timer = gtk_timeout_add 
	    (OBS_BUTTON_TIMER_DELAY, (GtkFunction) obs_button_timer, 
	     (gpointer) obs_button);
	}
      else 
	{
	  retval = TRUE;
	}
    }

  GDK_THREADS_LEAVE ();

  return retval;
}

extern gboolean update_obs_label (GtkEditable *entry, gpointer data);

static void default_output (ObsButton *obs_button)
{
    char buf[OBSLEN];
    gpointer data;

    ntodate(buf, (int) obs_button->adjustment->value, datainfo);

    if (strcmp (buf, gtk_entry_get_text(GTK_ENTRY(obs_button)))) {
	gtk_entry_set_text(GTK_ENTRY(obs_button), buf);
	data = gtk_object_get_data(GTK_OBJECT(obs_button), "rset");
	if (data != NULL) {
	    update_obs_label(NULL, data);
	}
    }
}

static void
obs_button_value_changed (GtkAdjustment *adjustment,
			  ObsButton *obs_button)
{
  g_return_if_fail (adjustment != NULL);
  g_return_if_fail (GTK_IS_ADJUSTMENT (adjustment));

  default_output(obs_button);

  obs_button_draw_arrow (obs_button, GTK_ARROW_UP);
  obs_button_draw_arrow (obs_button, GTK_ARROW_DOWN);
}

static gint
obs_button_key_press (GtkWidget     *widget,
		      GdkEventKey   *event)
{
  ObsButton *spin;
  gint key;
  gboolean key_repeat = FALSE;

  g_return_val_if_fail (widget != NULL, FALSE);
  g_return_val_if_fail (GTK_IS_OBS_BUTTON (widget), FALSE);
  g_return_val_if_fail (event != NULL, FALSE);
  
  spin = OBS_BUTTON (widget);
  key = event->keyval;

  key_repeat = (event->time == spin->ev_time);

  if (GTK_EDITABLE (widget)->editable &&
      (key == GDK_Up || key == GDK_Down || 
       key == GDK_Page_Up || key == GDK_Page_Down))
    obs_button_update (spin);

  switch (key)
    {
    case GDK_Up:

      if (GTK_WIDGET_HAS_FOCUS (widget))
	{
	  gtk_signal_emit_stop_by_name (GTK_OBJECT (widget), 
					"key_press_event");
	  if (!key_repeat)
	    spin->timer_step = spin->adjustment->step_increment;

	 obs_button_real_spin (spin, spin->timer_step);

	  return TRUE;
	}
      return FALSE;

    case GDK_Down:

      if (GTK_WIDGET_HAS_FOCUS (widget))
	{
	  gtk_signal_emit_stop_by_name (GTK_OBJECT (widget), 
					"key_press_event");
	  if (!key_repeat)
	    spin->timer_step = spin->adjustment->step_increment;

	 obs_button_real_spin (spin, -spin->timer_step);

	  return TRUE;
	}
      return FALSE;

    case GDK_Page_Up:

      if (event->state & GDK_CONTROL_MASK)
	{
	  gfloat diff = spin->adjustment->upper - spin->adjustment->value;
	  if (diff > EPSILON)
	    obs_button_real_spin (spin, diff);
	}
      else
	obs_button_real_spin (spin, spin->adjustment->page_increment);
      return TRUE;

    case GDK_Page_Down:

      if (event->state & GDK_CONTROL_MASK)
	{
	  gfloat diff = spin->adjustment->value - spin->adjustment->lower;
	  if (diff > EPSILON)
	    obs_button_real_spin (spin, -diff);
	}
      else
	obs_button_real_spin (spin, -spin->adjustment->page_increment);
      return TRUE;

    default:
      break;
    }

  return GTK_WIDGET_CLASS (parent_class)->key_press_event (widget, event);
}

static gint
obs_button_key_release (GtkWidget   *widget,
			     GdkEventKey *event)
{
  ObsButton *spin;

  g_return_val_if_fail (widget != NULL, FALSE);
  g_return_val_if_fail (GTK_IS_OBS_BUTTON (widget), FALSE);
  
  spin = OBS_BUTTON (widget);
  
  spin->ev_time = event->time;
  return TRUE;
}

void 
obs_button_update (ObsButton *obs_button)
{
    ObsButton *other; 
    gfloat val;

    g_return_if_fail (obs_button != NULL);
    g_return_if_fail (GTK_IS_OBS_BUTTON (obs_button));

    other = gtk_object_get_data(GTK_OBJECT(obs_button), "startspin");
    if (other != NULL) {
	obs_button->adjustment->lower = obs_button_get_value(other);
    } else {
	other = gtk_object_get_data(GTK_OBJECT(obs_button), "endspin");
	if (other != NULL) {
	    obs_button->adjustment->upper = obs_button_get_value(other);
	}
    }

    val = dateton(gtk_entry_get_text(GTK_ENTRY (obs_button)), datainfo);

    if (val < obs_button->adjustment->lower)
	val = obs_button->adjustment->lower;
    else if (val > obs_button->adjustment->upper)
	val = obs_button->adjustment->upper;

    if (fabs (val - obs_button->adjustment->value) > EPSILON)
	gtk_adjustment_set_value (obs_button->adjustment, val);
    else {
	default_output(obs_button);
    }
}

static void
obs_button_activate (GtkEditable *editable)
{
  g_return_if_fail (editable != NULL);
  g_return_if_fail (GTK_IS_OBS_BUTTON (editable));

  if (editable->editable)
    obs_button_update (OBS_BUTTON (editable));
}

static void
obs_button_insert_text (GtkEditable *editable,
			const gchar *new_text,
			gint         new_text_length,
			gint        *position)
{
    g_return_if_fail (editable != NULL);
    g_return_if_fail (GTK_IS_OBS_BUTTON (editable));

    GTK_EDITABLE_CLASS (parent_class)->insert_text (editable, new_text,
						    new_text_length, position);
}

static void
obs_button_real_spin (ObsButton *obs_button,
		      gfloat         increment)
{
  GtkAdjustment *adj;
  gfloat new_value = 0.0;

  g_return_if_fail (obs_button != NULL);
  g_return_if_fail (GTK_IS_OBS_BUTTON (obs_button));
  
  adj = obs_button->adjustment;

  new_value = adj->value + increment;

  if (increment > 0)
    {
	new_value = MIN (new_value, adj->upper);
    }
  else if (increment < 0) 
    {
	new_value = MAX (new_value, adj->lower);
    }

  if (fabs (new_value - adj->value) > EPSILON)
    gtk_adjustment_set_value (adj, new_value);
}


/***********************************************************
 ***********************************************************
 ***                  Public interface                   ***
 ***********************************************************
 ***********************************************************/

GtkWidget *
obs_button_new (GtkAdjustment *adjustment)
{
    ObsButton *spin;

    if (adjustment)
	g_return_val_if_fail (GTK_IS_ADJUSTMENT (adjustment), NULL);

    spin = gtk_type_new (GTK_TYPE_OBS_BUTTON);

    obs_button_set_adjustment(spin, adjustment);
    gtk_adjustment_value_changed (adjustment);

    return GTK_WIDGET (spin);
}

/* Callback used when the spin button's adjustment changes.  We need to redraw
 * the arrows when the adjustment's range changes.
 */
static void
adjustment_changed_cb (GtkAdjustment *adjustment, gpointer data)
{
  ObsButton *obs_button;

  obs_button = OBS_BUTTON (data);

  obs_button_draw_arrow (obs_button, GTK_ARROW_UP);
  obs_button_draw_arrow (obs_button, GTK_ARROW_DOWN);
}

static void
obs_button_set_adjustment (ObsButton *obs_button,
			   GtkAdjustment *adjustment)
{
  g_return_if_fail (obs_button != NULL);
  g_return_if_fail (GTK_IS_OBS_BUTTON (obs_button));

  if (obs_button->adjustment != adjustment)
    {
      if (obs_button->adjustment)
        {
          gtk_signal_disconnect_by_data (GTK_OBJECT (obs_button->adjustment),
                                         (gpointer) obs_button);
          gtk_object_unref (GTK_OBJECT (obs_button->adjustment));
        }
      obs_button->adjustment = adjustment;
      if (adjustment)
        {
          gtk_object_ref (GTK_OBJECT (adjustment));
	  gtk_object_sink (GTK_OBJECT (adjustment));
          gtk_signal_connect (GTK_OBJECT (adjustment), "value_changed",
			      (GtkSignalFunc) obs_button_value_changed,
			      (gpointer) obs_button);
	  gtk_signal_connect (GTK_OBJECT (adjustment), "changed",
			      (GtkSignalFunc) adjustment_changed_cb,
			      (gpointer) obs_button);
        }
    }
}

GtkAdjustment *
obs_button_get_adjustment (ObsButton *obs_button)
{
  g_return_val_if_fail (obs_button != NULL, NULL);
  g_return_val_if_fail (GTK_IS_OBS_BUTTON (obs_button), NULL);

  return obs_button->adjustment;
}

void 
obs_button_set_value (ObsButton *obs_button, 
		      gfloat         value)
{
    g_return_if_fail (obs_button != NULL);
    g_return_if_fail (GTK_IS_OBS_BUTTON (obs_button));

    if (fabs (value - obs_button->adjustment->value) > EPSILON)
	gtk_adjustment_set_value (obs_button->adjustment, value);
    else
	default_output(obs_button);
}

gfloat
obs_button_get_value (ObsButton *obs_button)
{
  g_return_val_if_fail (obs_button != NULL, 0.0);
  g_return_val_if_fail (GTK_IS_OBS_BUTTON (obs_button), 0.0);

  return obs_button->adjustment->value;
}

static void
obs_button_set_shadow_type (ObsButton *obs_button,
			    GtkShadowType  shadow_type)
{
    g_return_if_fail (obs_button != NULL);
    g_return_if_fail (GTK_IS_OBS_BUTTON (obs_button));

    if (shadow_type != obs_button->shadow_type)
	{
	    obs_button->shadow_type = shadow_type;
	    if (GTK_WIDGET_DRAWABLE (obs_button))
		gtk_widget_queue_draw (GTK_WIDGET (obs_button));
	}
}

void
obs_button_spin (ObsButton *obs_button,
		 GtkSpinType    direction,
		 gfloat         increment)
{
  GtkAdjustment *adj;
  gfloat diff;

  g_return_if_fail (obs_button != NULL);
  g_return_if_fail (GTK_IS_OBS_BUTTON (obs_button));
  
  adj = obs_button->adjustment;

  /* for compatibility with the 1.0.x version of this function */
  if (increment != 0 && increment != adj->step_increment &&
      (direction == GTK_SPIN_STEP_FORWARD ||
       direction == GTK_SPIN_STEP_BACKWARD))
    {
      if (direction == GTK_SPIN_STEP_BACKWARD && increment > 0)
	increment = -increment;
      direction = GTK_SPIN_USER_DEFINED;
    }

  switch (direction)
    {
    case GTK_SPIN_STEP_FORWARD:

      obs_button_real_spin (obs_button, adj->step_increment);
      break;

    case GTK_SPIN_STEP_BACKWARD:

      obs_button_real_spin (obs_button, -adj->step_increment);
      break;

    case GTK_SPIN_PAGE_FORWARD:

      obs_button_real_spin (obs_button, adj->page_increment);
      break;

    case GTK_SPIN_PAGE_BACKWARD:

      obs_button_real_spin (obs_button, -adj->page_increment);
      break;

    case GTK_SPIN_HOME:

      diff = adj->value - adj->lower;
      if (diff > EPSILON)
	obs_button_real_spin (obs_button, -diff);
      break;

    case GTK_SPIN_END:

      diff = adj->upper - adj->value;
      if (diff > EPSILON)
	obs_button_real_spin (obs_button, diff);
      break;

    case GTK_SPIN_USER_DEFINED:

      if (increment != 0)
	obs_button_real_spin (obs_button, increment);
      break;

    default:
      break;
    }
}
