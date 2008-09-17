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

/*
 * A custon spin-button for gretl, based on gtkspinbutton.
 * Allin Cottrell, November 2003.
 */

#include "gretl.h"

#include <gdk/gdk.h>
#include <gtk/gtkentry.h>
#include <gtk/gtkadjustment.h>
#include <gtk/gtkmain.h>

#include "obsbutton.h"

#define MIN_OBS_BUTTON_WIDTH              30
#define OBS_BUTTON_INITIAL_TIMER_DELAY    200
#define OBS_BUTTON_TIMER_DELAY            20
#define MAX_TIMER_CALLS                    5
#define EPSILON                            1e-10
#define MIN_ARROW_WIDTH			   6

enum {
  PROP_0,
  PROP_VALUE
};

struct _ObsButton
{
    GtkEntry entry;
    GtkAdjustment *adjustment;
    const DATAINFO *pdinfo;
    GdkWindow *panel;
    guint32 timer;
    gdouble timer_step;
    guint in_child : 2;
    guint click_child : 2; /* valid: GTK_ARROW_UP=0, GTK_ARROW_DOWN=1 or 2=NONE/BOTH */
    guint button : 2;
    guint need_timer : 1;
    guint timer_calls : 3;
};

struct _ObsButtonClass
{
    GtkEntryClass parent_class;

    gint (*input)  (ObsButton *obs_button,
		    gdouble       *new_value);
    gint (*output) (ObsButton *obs_button);
    void (*value_changed) (ObsButton *obs_button);

    /* Action signals for keybindings, do not connect to these */
    void (*change_value) (ObsButton *obs_button,
			  GtkScrollType scroll);
};

static void obs_button_class_init     (ObsButtonClass *klass);
static void obs_button_editable_init  (GtkEditableClass   *iface);
static void obs_button_init           (ObsButton      *obs_button);
static void obs_button_finalize       (GObject            *object);
static void obs_button_destroy        (GtkObject          *object);
static void obs_button_set_property   (GObject         *object,
				       guint            prop_id,
				       const GValue    *value,
				       GParamSpec      *pspec);
static void obs_button_get_property   (GObject         *object,
				       guint            prop_id,
				       GValue          *value,
				       GParamSpec      *pspec);
static void obs_button_map            (GtkWidget          *widget);
static void obs_button_unmap          (GtkWidget          *widget);
static void obs_button_realize        (GtkWidget          *widget);
static void obs_button_unrealize      (GtkWidget          *widget);
static void obs_button_size_request   (GtkWidget          *widget,
				       GtkRequisition     *requisition);
static void obs_button_size_allocate  (GtkWidget          *widget,
				       GtkAllocation      *allocation);
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
static void obs_button_grab_notify    (GtkWidget          *widget,
				       gboolean            was_grabbed);
static void obs_button_state_changed  (GtkWidget          *widget,
				       GtkStateType        previous_state);
static void obs_button_draw_arrow     (ObsButton      *obs_button, 
				       GtkArrowType        arrow_type);
static gint obs_button_timer          (ObsButton      *obs_button);
static void obs_button_stop_spinning  (ObsButton      *spin);
static void obs_button_value_changed  (GtkAdjustment      *adjustment,
				       ObsButton      *obs_button); 
static gint obs_button_key_release    (GtkWidget          *widget,
				       GdkEventKey        *event);
static gint obs_button_scroll         (GtkWidget          *widget,
				       GdkEventScroll     *event);
static void obs_button_activate       (GtkEntry           *entry);
static void obs_button_insert_text    (GtkEditable        *editable,
				       const gchar        *new_text,
				       gint                new_text_length,
				       gint               *position);
static void obs_button_real_spin      (ObsButton      *obs_button,
				       gdouble             step);
static void obs_button_real_change_value (ObsButton   *spin,
					  GtkScrollType    scroll);

static gint obs_button_default_input  (ObsButton      *obs_button,
				       gdouble            *new_val);
static gint obs_button_default_output (ObsButton      *obs_button);

static gint obs_button_get_arrow_size     (ObsButton      *obs_button);
static gint obs_button_get_shadow_type    (ObsButton      *spin_button);
static void obs_button_redraw             (ObsButton      *obs_button);
static void obs_button_set_adjustment     (ObsButton      *obs_button,
					   GtkAdjustment *adjustment);


static GtkEntryClass *parent_class = NULL;

#define NO_ARROW 2

GType obs_button_get_type (void)
{
    static GType obs_button_type = 0;

    if (!obs_button_type) {
	static const GTypeInfo obs_button_info = {
	    sizeof (ObsButtonClass),
	    NULL,		/* base_init */
	    NULL,		/* base_finalize */
	    (GClassInitFunc) obs_button_class_init,
	    NULL,		/* class_finalize */
	    NULL,		/* class_data */
	    sizeof (ObsButton),
	    0,		/* n_preallocs */
	    (GInstanceInitFunc) obs_button_init,
	    NULL
	};

	static const GInterfaceInfo editable_info = {
	    (GInterfaceInitFunc) obs_button_editable_init, /* interface_init */
	    NULL, /* interface_finalize */
	    NULL  /* interface_data */
	};

	obs_button_type =
	    g_type_register_static(GTK_TYPE_ENTRY, "ObsButton",
				   &obs_button_info, 0);

	g_type_add_interface_static(obs_button_type,
				    GTK_TYPE_EDITABLE,
				    &editable_info);
    }
    return obs_button_type;
}

static void
obs_button_class_init (ObsButtonClass *class)
{
  GObjectClass     *gobject_class = G_OBJECT_CLASS(class);
  GtkObjectClass   *object_class;
  GtkWidgetClass   *widget_class;
  GtkEntryClass    *entry_class;

  object_class = (GtkObjectClass*)   class;
  widget_class = (GtkWidgetClass*)   class;
  entry_class  = (GtkEntryClass*)    class;

  parent_class = g_type_class_peek_parent (class);

  gobject_class->finalize = obs_button_finalize;

  gobject_class->set_property = obs_button_set_property;
  gobject_class->get_property = obs_button_get_property;

  object_class->destroy = obs_button_destroy;

  widget_class->map = obs_button_map;
  widget_class->unmap = obs_button_unmap;
  widget_class->realize = obs_button_realize;
  widget_class->unrealize = obs_button_unrealize;
  widget_class->size_request = obs_button_size_request;
  widget_class->size_allocate = obs_button_size_allocate;
  widget_class->expose_event = obs_button_expose;
  widget_class->scroll_event = obs_button_scroll;
  widget_class->button_press_event = obs_button_button_press;
  widget_class->button_release_event = obs_button_button_release;
  widget_class->motion_notify_event = obs_button_motion_notify;
  widget_class->key_release_event = obs_button_key_release;
  widget_class->enter_notify_event = obs_button_enter_notify;
  widget_class->leave_notify_event = obs_button_leave_notify;
  widget_class->focus_out_event = obs_button_focus_out;
  widget_class->grab_notify = obs_button_grab_notify;
  widget_class->state_changed = obs_button_state_changed;

  entry_class->activate = obs_button_activate;

  class->input = NULL;
  class->output = NULL;
  class->change_value = obs_button_real_change_value;

  g_object_class_install_property (gobject_class,
                                   PROP_VALUE,
                                   g_param_spec_double ("value",
							_("Value"),
							_("Reads the current value, or sets a new value"),
							-G_MAXDOUBLE,
							G_MAXDOUBLE,
							0.0,
							G_PARAM_READWRITE));  
  
  gtk_widget_class_install_style_property_parser (widget_class,
						  g_param_spec_enum ("shadow_type", "Shadow Type", NULL,
								     GTK_TYPE_SHADOW_TYPE,
								     GTK_SHADOW_IN,
								     G_PARAM_READABLE),
						  gtk_rc_property_parse_enum);
}

static void
obs_button_editable_init (GtkEditableClass *iface)
{
    iface->insert_text = obs_button_insert_text;
}

static void
obs_button_set_property (GObject      *object,
			 guint         prop_id,
			 const GValue *value,
			 GParamSpec   *pspec)
{
    ObsButton *obs_button = OBS_BUTTON(object);

    if (prop_id == PROP_VALUE) {
	obs_button_set_value(obs_button, g_value_get_double(value));
    }
}

static void
obs_button_get_property (GObject      *object,
			 guint         prop_id,
			 GValue       *value,
			 GParamSpec   *pspec)
{
    ObsButton *obs_button = OBS_BUTTON (object);

    if (prop_id == PROP_VALUE) {
	g_value_set_double (value, obs_button->adjustment->value);
    } else {
	G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
    }
}

static void
obs_button_init (ObsButton *obs_button)
{
    obs_button->adjustment = NULL;
    obs_button->panel = NULL;
    obs_button->timer = 0;
    obs_button->timer_step = 0.0;
    obs_button->in_child = NO_ARROW;
    obs_button->click_child = NO_ARROW;
    obs_button->button = 0;
    obs_button->need_timer = FALSE;
    obs_button->timer_calls = 0;
    obs_button->adjustment = NULL;
}

static void
obs_button_finalize (GObject *object)
{
    obs_button_set_adjustment(OBS_BUTTON(object), NULL);
  
    G_OBJECT_CLASS(parent_class)->finalize(object);
}

static void
obs_button_destroy (GtkObject *object)
{
    obs_button_stop_spinning(OBS_BUTTON(object));
  
    GTK_OBJECT_CLASS(parent_class)->destroy(object);
}

static void
obs_button_map (GtkWidget *widget)
{
    if (GTK_WIDGET_REALIZED(widget) && !GTK_WIDGET_MAPPED(widget)) {
	GTK_WIDGET_CLASS(parent_class)->map(widget);
	gdk_window_show(OBS_BUTTON(widget)->panel);
    }
}

static void
obs_button_unmap (GtkWidget *widget)
{
    if (GTK_WIDGET_MAPPED(widget)) {
	gdk_window_hide(OBS_BUTTON(widget)->panel);
	GTK_WIDGET_CLASS(parent_class)->unmap(widget);
    }
}

static void
obs_button_realize (GtkWidget *widget)
{
    ObsButton *obs_button;
    GdkWindowAttr attributes;
    gint attributes_mask;
    guint real_width;
    gint arrow_size;

    obs_button = OBS_BUTTON(widget);
    arrow_size = obs_button_get_arrow_size(obs_button);

    real_width = widget->allocation.width;
    widget->allocation.width -= arrow_size + 2 * widget->style->xthickness;
    gtk_widget_set_events(widget, gtk_widget_get_events(widget) |
			  GDK_KEY_RELEASE_MASK);
    GTK_WIDGET_CLASS(parent_class)->realize(widget);

    widget->allocation.width = real_width;
  
    attributes.window_type = GDK_WINDOW_CHILD;
    attributes.wclass = GDK_INPUT_OUTPUT;
    attributes.visual = gtk_widget_get_visual(widget);
    attributes.colormap = gtk_widget_get_colormap(widget);
    attributes.event_mask = gtk_widget_get_events(widget);
    attributes.event_mask |= GDK_EXPOSURE_MASK | GDK_BUTTON_PRESS_MASK 
	| GDK_BUTTON_RELEASE_MASK | GDK_LEAVE_NOTIFY_MASK | GDK_ENTER_NOTIFY_MASK 
	| GDK_POINTER_MOTION_MASK | GDK_POINTER_MOTION_HINT_MASK;

    attributes_mask = GDK_WA_X | GDK_WA_Y | GDK_WA_VISUAL | GDK_WA_COLORMAP;

    attributes.x = (widget->allocation.x +
		    widget->allocation.width - arrow_size -
		    2 * widget->style->xthickness);
    attributes.y = widget->allocation.y + (widget->allocation.height -
					   widget->requisition.height) / 2;
    attributes.width = arrow_size + 2 * widget->style->xthickness;
    attributes.height = widget->requisition.height;
  
    obs_button->panel = gdk_window_new(gtk_widget_get_parent_window(widget), 
				       &attributes, attributes_mask);
    gdk_window_set_user_data(obs_button->panel, widget);

    gtk_style_set_background(widget->style, obs_button->panel, GTK_STATE_NORMAL);

    obs_button_default_output(obs_button);

    gtk_widget_queue_resize(GTK_WIDGET(obs_button));
}

static void
obs_button_unrealize (GtkWidget *widget)
{
    ObsButton *spin = OBS_BUTTON(widget);

    GTK_WIDGET_CLASS(parent_class)->unrealize(widget);

    if (spin->panel) {
	gdk_window_set_user_data(spin->panel, NULL);
	gdk_window_destroy(spin->panel);
	spin->panel = NULL;
    }
}

static void
obs_button_size_request (GtkWidget      *widget,
			 GtkRequisition *requisition)
{
    GtkEntry *entry;
    ObsButton *obs_button;
    gint arrow_size;

    entry = GTK_ENTRY (widget);
    obs_button = OBS_BUTTON (widget);
    arrow_size = obs_button_get_arrow_size (obs_button);
  
    GTK_WIDGET_CLASS (parent_class)->size_request (widget, requisition);

    if (entry->width_chars < 0) {
	PangoContext *context;
	PangoFontMetrics *metrics;
	gint width;
	gint w;
	gint string_len;
	gint max_string_len;
	gint digit_width;
	gboolean interior_focus;
	gint focus_width;

	gtk_widget_style_get(widget,
			     "interior-focus", &interior_focus,
			     "focus-line-width", &focus_width,
			     NULL);

	context = gtk_widget_get_pango_context(widget);
	metrics = pango_context_get_metrics(context,
					    widget->style->font_desc,
					    pango_context_get_language(context));

	digit_width = pango_font_metrics_get_approximate_digit_width(metrics);
	digit_width = PANGO_SCALE *
	    ((digit_width + PANGO_SCALE - 1) / PANGO_SCALE);

	pango_font_metrics_unref(metrics);
      
	/* Get max of MIN_OBS_BUTTON_WIDTH, size of upper, size of lower */
      
	width = MIN_OBS_BUTTON_WIDTH;
	max_string_len = OBSLEN;
	string_len = strlen(obs_button->pdinfo->endobs) + 1;
	w = PANGO_PIXELS(MIN(string_len, max_string_len) * digit_width);
	width = MAX(width, w);
      
	requisition->width = width + /* INNER_BORDER */ 2 * 2;
	if (!interior_focus)
	    requisition->width += 2 * focus_width;
    }

    requisition->width += arrow_size + 2 * widget->style->xthickness;
}

static void
obs_button_size_allocate (GtkWidget     *widget,
			  GtkAllocation *allocation)
{
    ObsButton *spin;
    GtkAllocation entry_allocation;
    GtkAllocation panel_allocation;
    gint arrow_size;
    gint panel_width;

    g_return_if_fail(GTK_IS_OBS_BUTTON(widget));
    g_return_if_fail(allocation != NULL);

    spin = OBS_BUTTON(widget);
    arrow_size = obs_button_get_arrow_size(spin);
    panel_width = arrow_size + 2 * widget->style->xthickness;
  
    widget->allocation = *allocation;
  
    entry_allocation = *allocation;
    entry_allocation.width -= panel_width;

    if (gtk_widget_get_direction(widget) == GTK_TEXT_DIR_RTL) {
	entry_allocation.x += panel_width;
	panel_allocation.x = allocation->x;
    } else {
	panel_allocation.x = allocation->x + allocation->width - panel_width;
    }

    panel_allocation.width = panel_width;
    panel_allocation.height = MIN(widget->requisition.height, allocation->height);

    panel_allocation.y = allocation->y + (allocation->height -
					  panel_allocation.height) / 2;

    GTK_WIDGET_CLASS(parent_class)->size_allocate(widget, &entry_allocation);

    if (GTK_WIDGET_REALIZED(widget)) {
	gdk_window_move_resize(OBS_BUTTON(widget)->panel, 
			       panel_allocation.x,
			       panel_allocation.y,
			       panel_allocation.width,
			       panel_allocation.height); 
    }

    obs_button_redraw(spin);
}

static gint
obs_button_expose (GtkWidget      *widget,
		   GdkEventExpose *event)
{
    ObsButton *spin;

    g_return_val_if_fail(GTK_IS_OBS_BUTTON(widget), FALSE);
    g_return_val_if_fail(event != NULL, FALSE);

    spin = OBS_BUTTON(widget);

    if (GTK_WIDGET_DRAWABLE(widget)) {
	GtkShadowType shadow_type;
	GdkRectangle rect;

	if (event->window != spin->panel)
	    GTK_WIDGET_CLASS(parent_class)->expose_event(widget, event);

	/* we redraw the panel even if it wasn't exposed. This is
	 * because spin->panel is not a child window of widget->window,
	 * so it will not be invalidated by eg. gtk_widget_queue_draw().
	 */
	rect.x = 0;
	rect.y = 0;

	gdk_drawable_get_size(spin->panel, &rect.width, &rect.height);

	shadow_type = obs_button_get_shadow_type(spin);
      
	gdk_window_begin_paint_rect(spin->panel, &rect);      

	if (shadow_type != GTK_SHADOW_NONE) {
	    gtk_paint_box(widget->style, spin->panel,
			  GTK_STATE_NORMAL, shadow_type,
			  NULL, widget, "spinbutton",
			  rect.x, rect.y, rect.width, rect.height);
	}

	obs_button_draw_arrow(spin, GTK_ARROW_UP);
	obs_button_draw_arrow(spin, GTK_ARROW_DOWN);

	gdk_window_end_paint(spin->panel);
    }
  
    return FALSE;
}

static gboolean
obs_button_at_limit (ObsButton *obs_button,
                     GtkArrowType   arrow)
{
    ObsButton *other;
    gpointer p;
    int minsep = 0;
    gdouble val;

    p = g_object_get_data(G_OBJECT(obs_button), "minsep");
    if (p != NULL) {
	minsep = GPOINTER_TO_INT(p);
    }

    other = g_object_get_data(G_OBJECT(obs_button), "startspin");
    if (other != NULL) {
	/* we're looking at the end spinner */
	val = obs_button_get_value(other);
	obs_button->adjustment->lower = val + minsep;
    } else {
	other = g_object_get_data(G_OBJECT(obs_button), "endspin");
	if (other != NULL) {
	    /* we're looking at the start spinner */
	    val = obs_button_get_value(other);
	    obs_button->adjustment->upper = val - minsep;
	}
    }

    if (arrow == GTK_ARROW_UP &&
	(obs_button->adjustment->upper - obs_button->adjustment->value <= EPSILON))
	return TRUE;
  
    if (arrow == GTK_ARROW_DOWN &&
	(obs_button->adjustment->value - obs_button->adjustment->lower <= EPSILON))
	return TRUE;
  
    return FALSE;
}

static void
obs_button_draw_arrow (ObsButton *obs_button, 
		       GtkArrowType   arrow_type)
{
    GtkStateType state_type;
    GtkShadowType shadow_type;
    GtkWidget *widget;
    gint x;
    gint y;
    gint height;
    gint width;
    gint h, w;

    g_return_if_fail (GTK_IS_OBS_BUTTON (obs_button));
    g_return_if_fail (arrow_type == GTK_ARROW_UP || arrow_type == GTK_ARROW_DOWN);
  
    widget = GTK_WIDGET (obs_button);

    if (GTK_WIDGET_DRAWABLE (widget)) {
	width = obs_button_get_arrow_size (obs_button) + 2 * widget->style->xthickness;

	if (arrow_type == GTK_ARROW_UP) {
	    x = 0;
	    y = 0;

	    height = widget->requisition.height / 2;
	} else {
	    x = 0;
	    y = widget->requisition.height / 2;

	    height = (widget->requisition.height + 1) / 2;
	}

	if (obs_button_at_limit (obs_button, arrow_type)) {
	    shadow_type = GTK_SHADOW_OUT;
	    state_type = GTK_STATE_INSENSITIVE;
	} else {
	    if (obs_button->click_child == arrow_type) {
		state_type = GTK_STATE_ACTIVE;
		shadow_type = GTK_SHADOW_IN;
	    } else {
		if (obs_button->in_child == arrow_type &&
		    obs_button->click_child == NO_ARROW) {
		    state_type = GTK_STATE_PRELIGHT;
		} else {
		    state_type = GTK_WIDGET_STATE (widget);
		}
	      
		shadow_type = GTK_SHADOW_OUT;
	    }
	}
      
	gtk_paint_box (widget->style, obs_button->panel,
		       state_type, shadow_type,
		       NULL, widget,
		       (arrow_type == GTK_ARROW_UP)? "spinbutton_up" : "spinbutton_down",
		       x, y, width, height);

	height = widget->requisition.height;

	if (arrow_type == GTK_ARROW_DOWN) {
	    y = height / 2;
	    height = height - y - 2;
	} else {
	    y = 2;
	    height = height / 2 - 2;
	}

	width -= 3;

	if (widget && gtk_widget_get_direction (widget) == GTK_TEXT_DIR_RTL)
	    x = 2;
	else
	    x = 1;

	w = width / 2;
	w -= w % 2 - 1; /* force odd */
	h = (w + 1) / 2;
      
	x += (width - w) / 2;
	y += (height - h) / 2;
      
	height = h;
	width = w;

	gtk_paint_arrow (widget->style, obs_button->panel,
			 state_type, shadow_type, 
			 NULL, widget, "spinbutton",
			 arrow_type, TRUE, 
			 x, y, width, height);
    }
}

static gint
obs_button_enter_notify (GtkWidget        *widget,
			 GdkEventCrossing *event)
{
    ObsButton *spin = OBS_BUTTON(widget);

    if (event->window == spin->panel) {
	gint x;
	gint y;

	gdk_window_get_pointer(spin->panel, &x, &y, NULL);

	if (y <= widget->requisition.height / 2)
	    spin->in_child = GTK_ARROW_UP;
	else
	    spin->in_child = GTK_ARROW_DOWN;

	obs_button_redraw (spin);
    }
  
    return FALSE;
}

static gint
obs_button_leave_notify (GtkWidget        *widget,
			 GdkEventCrossing *event)
{
    ObsButton *spin = OBS_BUTTON(widget);

    spin->in_child = NO_ARROW;
    obs_button_redraw (spin);
  
    return FALSE;
}

static gint
obs_button_focus_out (GtkWidget     *widget,
		      GdkEventFocus *event)
{
    if (GTK_ENTRY(widget)->editable)
	obs_button_update(OBS_BUTTON(widget));

    return GTK_WIDGET_CLASS(parent_class)->focus_out_event(widget, event);
}

static void
obs_button_grab_notify (GtkWidget *widget,
			gboolean   was_grabbed)
{
    ObsButton *spin = OBS_BUTTON(widget);

    if (!was_grabbed) {
	obs_button_stop_spinning(spin);
	obs_button_redraw(spin);
    }
}

static void
obs_button_state_changed (GtkWidget    *widget,
			  GtkStateType  previous_state)
{
    ObsButton *spin = OBS_BUTTON(widget);

    if (!GTK_WIDGET_IS_SENSITIVE(widget)) {
	obs_button_stop_spinning(spin);    
	obs_button_redraw(spin);
    }
}

static gint
obs_button_scroll (GtkWidget      *widget,
		   GdkEventScroll *event)
{
    ObsButton *spin = OBS_BUTTON (widget);

    if (event->direction == GDK_SCROLL_UP) {
	if (!GTK_WIDGET_HAS_FOCUS (widget))
	    gtk_widget_grab_focus (widget);
	obs_button_real_spin (spin, spin->adjustment->step_increment);
    }
    else if (event->direction == GDK_SCROLL_DOWN) {
	if (!GTK_WIDGET_HAS_FOCUS (widget))
	    gtk_widget_grab_focus (widget);
	obs_button_real_spin (spin, -spin->adjustment->step_increment); 
    }
    else
	return FALSE;

    return TRUE;
}

static void
obs_button_stop_spinning (ObsButton *spin)
{
    if (spin->timer) {
	g_source_remove (spin->timer);
	spin->timer = 0;
	spin->timer_calls = 0;
	spin->need_timer = FALSE;
    }

    spin->button = 0;
    spin->timer = 0;
    spin->timer_step = spin->adjustment->step_increment;
    spin->timer_calls = 0;

    spin->click_child = NO_ARROW;
    spin->button = 0;
}

static void
start_spinning (ObsButton *spin,
		GtkArrowType   click_child,
		gdouble        step)
{
    g_return_if_fail (click_child == GTK_ARROW_UP || click_child == GTK_ARROW_DOWN);
  
    spin->click_child = click_child;
    obs_button_real_spin (spin, click_child == GTK_ARROW_UP ? step : -step);
  
    if (!spin->timer) {
	spin->timer_step = step;
	spin->need_timer = TRUE;
	spin->timer = g_timeout_add (OBS_BUTTON_INITIAL_TIMER_DELAY, 
				     (GSourceFunc) obs_button_timer, (gpointer) spin);
    }

    obs_button_redraw (spin);
}

static gint
obs_button_button_press (GtkWidget      *widget,
			 GdkEventButton *event)
{
    ObsButton *spin = OBS_BUTTON (widget);

    if (!spin->button) {
	if (event->window == spin->panel) {
	    if (!GTK_WIDGET_HAS_FOCUS (widget))
		gtk_widget_grab_focus (widget);
	    spin->button = event->button;
	  
	    if (GTK_ENTRY (widget)->editable)
		obs_button_update (spin);
	  
	    if (event->y <= widget->requisition.height / 2) {
		if (event->button == 1)
		    start_spinning (spin, GTK_ARROW_UP, spin->adjustment->step_increment);
		else if (event->button == 2)
		    start_spinning (spin, GTK_ARROW_UP, spin->adjustment->page_increment);
		else
		    spin->click_child = GTK_ARROW_UP;
	    }
	    else {
		if (event->button == 1)
		    start_spinning (spin, GTK_ARROW_DOWN, spin->adjustment->step_increment);
		else if (event->button == 2)
		    start_spinning (spin, GTK_ARROW_DOWN, spin->adjustment->page_increment);
		else
		    spin->click_child = GTK_ARROW_DOWN;
	    }
	    return TRUE;
	}
	else
	    return GTK_WIDGET_CLASS (parent_class)->button_press_event (widget, event);
    }
    return FALSE;
}

static gint
obs_button_button_release (GtkWidget      *widget,
			   GdkEventButton *event)
{
    ObsButton *spin = OBS_BUTTON (widget);
    gint arrow_size;

    arrow_size = obs_button_get_arrow_size (spin);

    if (event->button == spin->button) {
	int click_child = spin->click_child;

	obs_button_stop_spinning (spin);

	if (event->button == 3) {
	    if (event->y >= 0 && event->x >= 0 && 
		event->y <= widget->requisition.height &&
		event->x <= arrow_size + 2 * widget->style->xthickness) {
		if (click_child == GTK_ARROW_UP &&
		    event->y <= widget->requisition.height / 2) {
		    gdouble diff;

		    diff = spin->adjustment->upper - spin->adjustment->value;
		    if (diff > EPSILON)
			obs_button_real_spin (spin, diff);
		}
		else if (click_child == GTK_ARROW_DOWN &&
			 event->y > widget->requisition.height / 2) {
		    gdouble diff;

		    diff = spin->adjustment->value - spin->adjustment->lower;
		    if (diff > EPSILON)
			obs_button_real_spin (spin, -diff);
		}
	    }
	}		  
	obs_button_redraw (spin);

	return TRUE;
    }
    else
	return GTK_WIDGET_CLASS (parent_class)->button_release_event (widget, event);
}

static gint
obs_button_motion_notify (GtkWidget      *widget,
			  GdkEventMotion *event)
{
    ObsButton *spin = OBS_BUTTON (widget);

    if (spin->button)
	return FALSE;

    if (event->window == spin->panel) {
	gint y;
      
	gdk_window_get_pointer (spin->panel, NULL, &y, NULL);
  
	if (y <= widget->requisition.height / 2 && 
	    spin->in_child == GTK_ARROW_DOWN) {
	    spin->in_child = GTK_ARROW_UP;
	    obs_button_redraw (spin);
	}
	else if (y > widget->requisition.height / 2 && 
		 spin->in_child == GTK_ARROW_UP) {
	    spin->in_child = GTK_ARROW_DOWN;
	    obs_button_redraw (spin);
	}
      
	return FALSE;
    }
	  
    return GTK_WIDGET_CLASS (parent_class)->motion_notify_event (widget, event);
}

static gboolean
obs_button_timer (ObsButton *obs_button)
{
    gboolean retval = FALSE;
  
    GDK_THREADS_ENTER ();

    if (obs_button->timer) {
	if (obs_button->click_child == GTK_ARROW_UP)
	    obs_button_real_spin (obs_button, obs_button->timer_step);
	else
	    obs_button_real_spin (obs_button, -obs_button->timer_step);

	if (obs_button->need_timer) {
	    obs_button->need_timer = FALSE;
	    obs_button->timer = g_timeout_add 
		(OBS_BUTTON_TIMER_DELAY, (GSourceFunc) obs_button_timer, 
		 (gpointer) obs_button);
	} else {
	    retval = TRUE;
	}
    }

    GDK_THREADS_LEAVE ();

    return retval;
}

static void
obs_button_value_changed (GtkAdjustment *adjustment,
			  ObsButton *obs_button)
{
    g_return_if_fail(GTK_IS_ADJUSTMENT(adjustment));

    obs_button_default_output(obs_button);

    obs_button_redraw(obs_button);
  
    g_object_notify(G_OBJECT(obs_button), "value");
}

static void
obs_button_real_change_value (ObsButton *spin,
			      GtkScrollType  scroll)
{
    /* We don't test whether the entry is editable, since
     * this key binding conceptually corresponds to changing
     * the value with the buttons using the mouse, which
     * we allow for non-editable spin buttons.
     */
    switch (scroll) {
    case GTK_SCROLL_STEP_BACKWARD:
    case GTK_SCROLL_STEP_DOWN:
    case GTK_SCROLL_STEP_LEFT:
	obs_button_real_spin (spin, -spin->timer_step);
	break;
      
    case GTK_SCROLL_STEP_FORWARD:
    case GTK_SCROLL_STEP_UP:
    case GTK_SCROLL_STEP_RIGHT:
	obs_button_real_spin (spin, spin->timer_step);
	break;
      
    case GTK_SCROLL_PAGE_BACKWARD:
    case GTK_SCROLL_PAGE_DOWN:
    case GTK_SCROLL_PAGE_LEFT:
	obs_button_real_spin (spin, -spin->adjustment->page_increment);
	break;
      
    case GTK_SCROLL_PAGE_FORWARD:
    case GTK_SCROLL_PAGE_UP:
    case GTK_SCROLL_PAGE_RIGHT:
	obs_button_real_spin (spin, spin->adjustment->page_increment);
	break;
      
    case GTK_SCROLL_START:
	{
	    gdouble diff = spin->adjustment->value - spin->adjustment->lower;
	    if (diff > EPSILON)
		obs_button_real_spin (spin, -diff);
	    break;
	}
      
    case GTK_SCROLL_END:
	{
	    gdouble diff = spin->adjustment->upper - spin->adjustment->value;
	    if (diff > EPSILON)
		obs_button_real_spin (spin, diff);
	    break;
	}
      
    default:
	g_warning ("Invalid scroll type %d for ObsButton::change-value", scroll);
	break;
    }
  
    obs_button_update (spin);
}

static gint
obs_button_key_release (GtkWidget   *widget,
			GdkEventKey *event)
{
    ObsButton *spin = OBS_BUTTON(widget);

    /* We only get a release at the end of a key repeat run, so reset the timer_step */
    spin->timer_step = spin->adjustment->step_increment;
    spin->timer_calls = 0;
  
    return TRUE;
}

static void
obs_button_activate (GtkEntry *entry)
{
    if (entry->editable)
	obs_button_update(OBS_BUTTON(entry));

    /* Chain up so that entry->activates_default is honored */
    parent_class->activate(entry);
}

static void
obs_button_insert_text (GtkEditable *editable,
			const gchar *new_text,
			gint         new_text_length,
			gint        *position)
{
    GtkEditableClass *parent_editable_iface = g_type_interface_peek(parent_class, GTK_TYPE_EDITABLE);
 
    parent_editable_iface->insert_text(editable, new_text,
				       new_text_length, position);
}

static void
obs_button_real_spin (ObsButton *obs_button,
		      gdouble        increment)
{
    GtkAdjustment *adj;
    gdouble new_value = 0.0;
  
    adj = obs_button->adjustment;

    new_value = adj->value + increment;

    if (increment > 0) {
	new_value = MIN(new_value, adj->upper);
    } else if (increment < 0) {
	new_value = MAX(new_value, adj->lower);
    }

    if (fabs (new_value - adj->value) > EPSILON)
	gtk_adjustment_set_value (adj, new_value);

    obs_button_redraw (obs_button);
}

static gint
obs_button_default_input (ObsButton *obs_button,
			  gdouble *new_val)
{
    *new_val = dateton(gtk_entry_get_text(GTK_ENTRY (obs_button)), 
		       obs_button->pdinfo);

    if (*new_val < 0)
	return GTK_INPUT_ERROR;
    else
	return FALSE;
}

extern gboolean update_obs_label (GtkEditable *entry, gpointer data);

static gint
obs_button_default_output (ObsButton *obs_button)
{
    int ival = obs_button->adjustment->value;
    gchar buf[OBSLEN];
    gpointer data;

    ntodate_full(buf, ival, obs_button->pdinfo);

    if (strcmp(buf, gtk_entry_get_text(GTK_ENTRY(obs_button))))
	gtk_entry_set_text(GTK_ENTRY(obs_button), buf);

    data = g_object_get_data(G_OBJECT(obs_button), "rset");
    if (data != NULL) {
	update_obs_label(NULL, data);
    }

    return FALSE;
}

/***********************************************************
 ***********************************************************
 ***                  Public interface                   ***
 ***********************************************************
 ***********************************************************/

GtkWidget *
obs_button_new (GtkAdjustment *adjustment, const DATAINFO *pdinfo)
{
    ObsButton *spin;

    if (adjustment)
	g_return_val_if_fail(GTK_IS_ADJUSTMENT(adjustment), NULL);

    spin = g_object_new(GTK_TYPE_OBS_BUTTON, NULL);

    spin->pdinfo = pdinfo;

    obs_button_set_adjustment(spin, adjustment);
    gtk_adjustment_value_changed(adjustment);

    return GTK_WIDGET(spin);
}

/* Callback used when the spin button's adjustment changes.  We need to redraw
 * the arrows when the adjustment's range changes, and reevaluate our size request.
 */
static void
adjustment_changed_cb (GtkAdjustment *adjustment, gpointer data)
{
    ObsButton *obs_button = OBS_BUTTON(data);

    obs_button->timer_step = obs_button->adjustment->step_increment;
    gtk_widget_queue_resize(GTK_WIDGET(obs_button));
}

/**
 * obs_button_set_adjustment:
 * @obs_button: a #ObsButton
 * @adjustment: a #GtkAdjustment to replace the existing adjustment
 * 
 * Replaces the #GtkAdjustment associated with @obs_button.
 **/
static void
obs_button_set_adjustment (ObsButton *obs_button,
			   GtkAdjustment *adjustment)
{
    g_return_if_fail (GTK_IS_OBS_BUTTON (obs_button));

    if (obs_button->adjustment != adjustment) {
	if (obs_button->adjustment) {
	    g_signal_handlers_disconnect_by_func (obs_button->adjustment,
						  obs_button_value_changed,
						  obs_button);
	    g_signal_handlers_disconnect_by_func (obs_button->adjustment,
						  adjustment_changed_cb,
						  obs_button);
	    g_object_unref (obs_button->adjustment);
	}
	obs_button->adjustment = adjustment;
	if (adjustment) {
	    g_object_ref (adjustment);
#if (GTK_MAJOR_VERSION == 2 && GTK_MINOR_VERSION < 10)
	    gtk_object_sink (GTK_OBJECT (adjustment));
#else
	    g_object_ref_sink (GTK_OBJECT (adjustment));
#endif
	    g_signal_connect (adjustment, "value_changed",
			      G_CALLBACK (obs_button_value_changed),
			      obs_button);
	    g_signal_connect (adjustment, "changed",
			      G_CALLBACK (adjustment_changed_cb),
			      obs_button);
	    obs_button->timer_step = obs_button->adjustment->step_increment;
	}

	gtk_widget_queue_resize (GTK_WIDGET (obs_button));
    }
}

/**
 * obs_button_get_value:
 * @obs_button: a #ObsButton
 * 
 * Get the value in the @obs_button.
 * 
 * Return value: the value of @obs_button
 **/
gdouble
obs_button_get_value (ObsButton *obs_button)
{
    g_return_val_if_fail (GTK_IS_OBS_BUTTON (obs_button), 0.0);

    return obs_button->adjustment->value;
}

/**
 * obs_button_set_value:
 * @obs_button: a #ObsButton
 * @value: the new value
 * 
 * Set the value of @obs_button.
 **/
void 
obs_button_set_value (ObsButton *obs_button, 
		      gdouble        value)
{
    g_return_if_fail (GTK_IS_OBS_BUTTON (obs_button));

    if (fabs(value - obs_button->adjustment->value) > EPSILON) {
	gtk_adjustment_set_value (obs_button->adjustment, value);
    } else {
	obs_button_default_output (obs_button);
    }
}

static gint
obs_button_get_arrow_size (ObsButton *obs_button)
{
    gint size = pango_font_description_get_size (GTK_WIDGET (obs_button)->style->font_desc);
    gint arrow_size;

    arrow_size = MAX (PANGO_PIXELS (size), MIN_ARROW_WIDTH);

    return arrow_size - arrow_size % 2; /* force even */
}

static gint
obs_button_get_shadow_type (ObsButton *spin_button)
{
    GtkShadowType rc_shadow_type;

    gtk_widget_style_get(GTK_WIDGET(spin_button), "shadow_type", &rc_shadow_type, NULL);

    return rc_shadow_type;
}

/**
 * obs_button_update:
 * @obs_button: a #ObsButton 
 * 
 * Manually force an update of the spin button.
 **/
void 
obs_button_update (ObsButton *obs_button)
{
    gdouble val;
    gint error = 0;
    gint return_val;

    g_return_if_fail(GTK_IS_OBS_BUTTON(obs_button));

    return_val = obs_button_default_input(obs_button, &val);
    error = (return_val == GTK_INPUT_ERROR);
    
    obs_button_redraw(obs_button);

    if (val < obs_button->adjustment->lower)
	val = obs_button->adjustment->lower;
    else if (val > obs_button->adjustment->upper)
	val = obs_button->adjustment->upper;

    if (fabs(val - obs_button->adjustment->value) > EPSILON) {
	gtk_adjustment_set_value(obs_button->adjustment, val);
    } else {
	obs_button_default_output(obs_button);
    }
}

static void
obs_button_redraw (ObsButton *obs_button)
{
    GtkWidget *widget;

    widget = GTK_WIDGET(obs_button);

    if (GTK_WIDGET_DRAWABLE(widget)) {
	gtk_widget_queue_draw(widget);

	/* We must invalidate the panel window ourselves, because it
	 * is not a child of widget->window
	 */
	gdk_window_invalidate_rect(obs_button->panel, NULL, TRUE);
    }
}        
