/* GTK - The GIMP Toolkit
 *
 * Copyright (C) 2007 John Stowers, Neil Jagdish Patel.
 * Copyright (C) 2009 Bastien Nocera, David Zeuthen
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA  02111-1307, USA.
 *
 * Code adapted from egg-spinner
 * by Christian Hergert <christian.hergert@gmail.com>
 */

/*
 * Modified by the GTK+ Team and others 2007.  See the AUTHORS
 * file for a list of people on the GTK+ Team.  See the ChangeLog
 * files for a list of changes.  These files are distributed with
 * GTK+ at ftp://ftp.gtk.org/pub/gtk/.
 */

/*
 * Back-ported and simplified for use with earlier versions of GTK+
 * (prior to 2.20) by Allin Cottrell, 2010-06-30.
 */

#include <math.h>
#include <gtk/gtk.h>
#include "spinner.h"

#define GTK_SPINNER_GET_PRIVATE(obj) (G_TYPE_INSTANCE_GET_PRIVATE ((obj), GTK_TYPE_SPINNER, GtkSpinnerPrivate))

G_DEFINE_TYPE (GtkSpinner, gtk_spinner, GTK_TYPE_DRAWING_AREA);

enum {
    PROP_0,
    PROP_ACTIVE
};

struct _GtkSpinnerPrivate
{
    guint current;
    guint num_steps;
    guint cycle_duration;
    gboolean active;
    guint timeout;
};

static void gtk_spinner_class_init     (GtkSpinnerClass *klass);
static void gtk_spinner_init           (GtkSpinner      *spinner);
static void gtk_spinner_dispose        (GObject         *gobject);
static void gtk_spinner_realize        (GtkWidget       *widget);
static void gtk_spinner_unrealize      (GtkWidget       *widget);
static gboolean gtk_spinner_expose     (GtkWidget       *widget,
                                        GdkEventExpose  *event);
static void gtk_spinner_screen_changed (GtkWidget       *widget,
                                        GdkScreen       *old_screen);
static void gtk_spinner_style_set      (GtkWidget       *widget,
                                        GtkStyle        *prev_style);
static void gtk_spinner_get_property   (GObject         *object,
                                        guint            param_id,
                                        GValue          *value,
                                        GParamSpec      *pspec);
static void gtk_spinner_set_property   (GObject         *object,
                                        guint            param_id,
                                        const GValue    *value,
                                        GParamSpec      *pspec);
static void gtk_spinner_set_active     (GtkSpinner      *spinner,
                                        gboolean         active);

static void
gtk_spinner_class_init (GtkSpinnerClass *klass)
{
    GObjectClass *gobject_class;
    GtkWidgetClass *widget_class;

    gobject_class = G_OBJECT_CLASS(klass);
    g_type_class_add_private (gobject_class, sizeof (GtkSpinnerPrivate));
    gobject_class->dispose = gtk_spinner_dispose;
    gobject_class->get_property = gtk_spinner_get_property;
    gobject_class->set_property = gtk_spinner_set_property;

    widget_class = GTK_WIDGET_CLASS(klass);
    widget_class->expose_event = gtk_spinner_expose;
    widget_class->realize = gtk_spinner_realize;
    widget_class->unrealize = gtk_spinner_unrealize;
    widget_class->screen_changed = gtk_spinner_screen_changed;
    widget_class->style_set = gtk_spinner_style_set;
    widget_class->get_accessible = NULL;

    g_object_class_install_property (gobject_class,
				     PROP_ACTIVE,
				     g_param_spec_boolean ("active",
							   "Active",
							   "Whether the spinner is active",
							   FALSE,
							   G_PARAM_READWRITE));

  gtk_widget_class_install_style_property (widget_class,
                                           g_param_spec_uint ("num-steps",
							      "Number of steps",
							      "Steps to complete a full loop.",
							      1,
							      G_MAXUINT,
							      12,
							      G_PARAM_READABLE));

  gtk_widget_class_install_style_property (widget_class,
                                           g_param_spec_uint ("cycle-duration",
							      "Animation duration",
							      "Milliseconds to complete a full loop",
							      500,
							      G_MAXUINT,
							      1000,
							      G_PARAM_READABLE));
}

static void
gtk_spinner_get_property (GObject    *object,
                          guint       param_id,
                          GValue     *value,
                          GParamSpec *pspec)
{
    GtkSpinnerPrivate *priv;

    priv = GTK_SPINNER (object)->priv;

    switch (param_id) {
    case PROP_ACTIVE:
        g_value_set_boolean (value, priv->active);
        break;
    default:
        G_OBJECT_WARN_INVALID_PROPERTY_ID (object, param_id, pspec);
    }
}

static void
gtk_spinner_set_property (GObject      *object,
                          guint         param_id,
                          const GValue *value,
                          GParamSpec   *pspec)
{
    switch (param_id) {
    case PROP_ACTIVE:
        gtk_spinner_set_active (GTK_SPINNER (object), g_value_get_boolean (value));
        break;
    default:
        G_OBJECT_WARN_INVALID_PROPERTY_ID (object, param_id, pspec);
    }
}

static void gtk_spinner_init (GtkSpinner *spinner)
{
    GtkSpinnerPrivate *priv;

    priv = GTK_SPINNER_GET_PRIVATE (spinner);
    priv->current = 0;
    priv->timeout = 0;

    spinner->priv = priv;

#if GTK_MAJOR_VERSION == 2 && GTK_MINOR_VERSION < 18
    GTK_WIDGET_SET_FLAGS (GTK_WIDGET (spinner), GTK_NO_WINDOW);
#else
    gtk_widget_set_has_window (GTK_WIDGET (spinner), FALSE);
#endif
}

static void
gtk_paint_spinner (GtkWidget    *widget,
		   GdkWindow    *window,
		   GtkStateType  state_type,
		   GdkRectangle *area,
		   const gchar  *detail,
		   guint         step,
		   gint          x,
		   gint          y,
		   gint          width,
		   gint          height)
{
    GtkStyle *style = widget->style;
    GdkColor *color;
    cairo_t *cr;
    guint num_steps = 12; /* the GTK default, hardwired */
    gdouble dx, dy;
    gdouble radius;
    gdouble half;
    gint i;
    guint real_step;

    real_step = step % num_steps;

    /* get cairo context */
    cr = gdk_cairo_create (window);

    /* set a clip region for the expose event */
    cairo_rectangle (cr, x, y, width, height);
    cairo_clip (cr);

    cairo_translate (cr, x, y);

    /* draw clip region */
    cairo_set_operator (cr, CAIRO_OPERATOR_OVER);

    color = &style->fg[state_type];
    dx = width / 2;
    dy = height / 2;
    radius = MIN (width / 2, height / 2);
    half = num_steps / 2;

    for (i = 0; i < num_steps; i++)
	{
	    gint inset = 0.7 * radius;

	    /* transparency is a function of time and intial value */
	    gdouble t = (gdouble) ((i + num_steps - real_step)
				   % num_steps) / num_steps;

	    cairo_save (cr);

	    cairo_set_source_rgba (cr,
				   color->red / 65535.,
				   color->green / 65535.,
				   color->blue / 65535.,
				   t);

	    cairo_set_line_width (cr, 2.0);
	    cairo_move_to (cr,
			   dx + (radius - inset) * cos (i * G_PI / half),
			   dy + (radius - inset) * sin (i * G_PI / half));
	    cairo_line_to (cr,
			   dx + radius * cos (i * G_PI / half),
			   dy + radius * sin (i * G_PI / half));
	    cairo_stroke (cr);

	    cairo_restore (cr);
	}

    /* free memory */
    cairo_destroy (cr);
}

static gboolean gtk_spinner_expose (GtkWidget *widget,
				    GdkEventExpose *event)
{
    GtkStateType state_type;
    GtkSpinnerPrivate *priv;
    int width, height;

    priv = GTK_SPINNER (widget)->priv;

    width = widget->allocation.width;
    height = widget->allocation.height;

    if ((width < 12) || (height <12))
	gtk_widget_set_size_request (widget, 12, 12);

    state_type = GTK_STATE_NORMAL;

#if GTK_MAJOR_VERSION == 2 && GTK_MINOR_VERSION < 18
    if (!GTK_WIDGET_IS_SENSITIVE (widget))
	state_type = GTK_STATE_INSENSITIVE;
#else
    if (!gtk_widget_is_sensitive (widget))
	state_type = GTK_STATE_INSENSITIVE;
#endif

    gtk_paint_spinner (widget,
		       widget->window,
		       state_type,
		       &event->area,
		       "spinner",
		       priv->current,
		       event->area.x, event->area.y,
		       event->area.width, event->area.height);

    return FALSE;
}

static gboolean gtk_spinner_timeout (gpointer data)
{
    GtkSpinnerPrivate *priv;

    priv = GTK_SPINNER (data)->priv;

    if (priv->current + 1 >= priv->num_steps) {
	priv->current = 0;
    } else {
	priv->current++;
    }

    gtk_widget_queue_draw (GTK_WIDGET (data));

    return TRUE;
}

static void gtk_spinner_add_timeout (GtkSpinner *spinner)
{
    GtkSpinnerPrivate *priv;

    priv = spinner->priv;

    priv->timeout = gdk_threads_add_timeout ((guint) priv->cycle_duration / priv->num_steps, 
					     gtk_spinner_timeout, spinner);
}

static void gtk_spinner_remove_timeout (GtkSpinner *spinner)
{
    GtkSpinnerPrivate *priv;

    priv = spinner->priv;

    g_source_remove (priv->timeout);
    priv->timeout = 0;
}

static void gtk_spinner_realize (GtkWidget *widget)
{
    GtkSpinnerPrivate *priv;

    priv = GTK_SPINNER (widget)->priv;

    GTK_WIDGET_CLASS (gtk_spinner_parent_class)->realize (widget);

    if (priv->active)
	gtk_spinner_add_timeout (GTK_SPINNER (widget));
}

static void gtk_spinner_unrealize (GtkWidget *widget)
{
    GtkSpinnerPrivate *priv;

    priv = GTK_SPINNER (widget)->priv;

    if (priv->timeout != 0) {
	gtk_spinner_remove_timeout (GTK_SPINNER (widget));
    }

    GTK_WIDGET_CLASS (gtk_spinner_parent_class)->unrealize (widget);
}

static void
gtk_spinner_screen_changed (GtkWidget* widget, GdkScreen* old_screen)
{
    GdkScreen* new_screen;
    GdkColormap* colormap;

    new_screen = gtk_widget_get_screen (widget);
    colormap = gdk_screen_get_rgba_colormap (new_screen);

    if (!colormap) {
	colormap = gdk_screen_get_rgb_colormap (new_screen);
    }

    gtk_widget_set_colormap (widget, colormap);
}

static void gtk_spinner_style_set (GtkWidget *widget,
				   GtkStyle  *prev_style)
{
    GtkSpinnerPrivate *priv;

    priv = GTK_SPINNER (widget)->priv;

    gtk_widget_style_get (GTK_WIDGET (widget),
			  "num-steps", &(priv->num_steps),
			  "cycle-duration", &(priv->cycle_duration),
			  NULL);

    if (priv->current > priv->num_steps)
	priv->current = 0;
}

static void gtk_spinner_dispose (GObject *gobject)
{
    GtkSpinnerPrivate *priv;

    priv = GTK_SPINNER (gobject)->priv;

    if (priv->timeout != 0) {
	gtk_spinner_remove_timeout (GTK_SPINNER (gobject));
    }

    G_OBJECT_CLASS (gtk_spinner_parent_class)->dispose (gobject);
}

static void gtk_spinner_set_active (GtkSpinner *spinner, gboolean active)
{
    GtkSpinnerPrivate *priv;

    active = active != FALSE;

    priv = GTK_SPINNER (spinner)->priv;

    if (priv->active != active) {
	gboolean spinner_realized;

	priv->active = active;
	g_object_notify (G_OBJECT (spinner), "active");

#if GTK_MAJOR_VERSION == 2 && GTK_MINOR_VERSION < 20
	spinner_realized = GTK_WIDGET_REALIZED (GTK_WIDGET (spinner));
#else
	spinner_realized = gtk_widget_get_realized (GTK_WIDGET (spinner));
#endif

	if (active && spinner_realized && priv->timeout == 0) {
	    gtk_spinner_add_timeout (spinner);
	} else if (!active && priv->timeout != 0) {
	    gtk_spinner_remove_timeout (spinner);
        }
    }
}

GtkWidget *gtk_spinner_new (void)
{
    return g_object_new (GTK_TYPE_SPINNER, NULL);
}

void gtk_spinner_start (GtkSpinner *spinner)
{
    g_return_if_fail (GTK_IS_SPINNER (spinner));

    gtk_spinner_set_active (spinner, TRUE);
}

void gtk_spinner_stop (GtkSpinner *spinner)
{
    g_return_if_fail (GTK_IS_SPINNER (spinner));

    gtk_spinner_set_active (spinner, FALSE);
}
