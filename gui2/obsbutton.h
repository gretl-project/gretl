/* GTK - The GIMP Toolkit
 * Copyright (C) 1995-1997 Peter Mattis, Spencer Kimball and Josh MacDonald
 *
 * GtkSpinButton widget for GTK+
 * Copyright (C) 1998 Lars Hamann and Stefan Jeske
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
 * Boston, MA 02111-1307, USA.
 */

/*
 * Modified by the GTK+ Team and others 1997-2000.  See the AUTHORS
 * file for a list of people on the GTK+ Team.  See the ChangeLog
 * files for a list of changes.  These files are distributed with
 * GTK+ at ftp://ftp.gtk.org/pub/gtk/. 
 */

#ifndef OBS_BUTTON_H__
#define OBS_BUTTON_H__

#include <gdk/gdk.h>
#include <gtk/gtkentry.h>
#include <gtk/gtkadjustment.h>

#define GTK_TYPE_OBS_BUTTON              (obs_button_get_type ())

#define OBS_BUTTON(obj)                  (G_TYPE_CHECK_INSTANCE_CAST ((obj), GTK_TYPE_OBS_BUTTON, ObsButton))
#define OBS_BUTTON_CLASS(klass)          (G_TYPE_CHECK_CLASS_CAST ((klass), GTK_TYPE_OBS_BUTTON, ObsButtonClass))
#define GTK_IS_OBS_BUTTON(obj)           (G_TYPE_CHECK_INSTANCE_TYPE ((obj), GTK_TYPE_OBS_BUTTON))
#define GTK_IS_OBS_BUTTON_CLASS(klass)   (G_TYPE_CHECK_CLASS_TYPE ((klass), GTK_TYPE_OBS_BUTTON))
#define OBS_BUTTON_GET_CLASS(obj)        (G_TYPE_INSTANCE_GET_CLASS ((obj), GTK_TYPE_OBS_BUTTON, ObsButtonClass))

#define GTK_INPUT_ERROR -1

typedef struct _ObsButton	    ObsButton;
typedef struct _ObsButtonClass  ObsButtonClass;


struct _ObsButton
{
    GtkEntry entry;
  
    GtkAdjustment *adjustment;
  
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
			  GtkScrollType  scroll);
};

GType		obs_button_get_type	   (void) G_GNUC_CONST;

GtkWidget*	obs_button_new		   (GtkAdjustment  *adjustment, const DATAINFO *pdinfo);

gdouble		obs_button_get_value       (ObsButton  *obs_button);

void		obs_button_set_value	   (ObsButton  *obs_button, 
					    gdouble	    value);

void            obs_button_update          (ObsButton  *obs_button);



#endif /* OBS_BUTTON_H__ */
