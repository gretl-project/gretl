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
#define OBS_BUTTON(obj)                  (GTK_CHECK_CAST ((obj), GTK_TYPE_OBS_BUTTON, ObsButton))
#define OBS_BUTTON_CLASS(klass)          (GTK_CHECK_CLASS_CAST ((klass), GTK_TYPE_OBS_BUTTON, ObsButtonClass))
#define GTK_IS_OBS_BUTTON(obj)           (GTK_CHECK_TYPE ((obj), GTK_TYPE_OBS_BUTTON))
#define GTK_IS_OBS_BUTTON_CLASS(klass)   (GTK_CHECK_CLASS_TYPE ((klass), GTK_TYPE_OBS_BUTTON))

typedef struct _ObsButton	    ObsButton;
typedef struct _ObsButtonClass  ObsButtonClass;

struct _ObsButton
{
    GtkEntry entry;
  
    GtkAdjustment *adjustment;
  
    GdkWindow *panel;
    GtkShadowType shadow_type;
  
    guint32 timer;
    guint32 ev_time;

    gfloat timer_step;

    guint in_child : 2;
    guint click_child : 2; 
    guint button : 2;
    guint need_timer : 1;
    guint timer_calls : 3;
};

struct _ObsButtonClass
{
    GtkEntryClass parent_class;
};

GtkType		obs_button_get_type	   (void);

GtkWidget*	obs_button_new		   (GtkAdjustment  *adjustment);

gfloat  	obs_button_get_value       (ObsButton  *obs_button);

void		obs_button_set_value	   (ObsButton  *obs_button, 
					    gfloat	    value);

void            obs_button_update          (ObsButton  *obs_button);


#endif /* OBS_BUTTON_H__ */
