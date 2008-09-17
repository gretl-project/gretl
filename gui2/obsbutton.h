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
 * Modified by the GTK+ Team and others 1997-2000.  See the AUTHORS
 * file for a list of people on the GTK+ Team.  See the ChangeLog
 * files for a list of changes.  These files are distributed with
 * GTK+ at ftp://ftp.gtk.org/pub/gtk/. 
 */

#ifndef OBS_BUTTON_H__
#define OBS_BUTTON_H__

#define GTK_TYPE_OBS_BUTTON            (obs_button_get_type())

#define OBS_BUTTON(obj)                (G_TYPE_CHECK_INSTANCE_CAST ((obj), GTK_TYPE_OBS_BUTTON, ObsButton))
#define OBS_BUTTON_CLASS(klass)        (G_TYPE_CHECK_CLASS_CAST ((klass), GTK_TYPE_OBS_BUTTON, ObsButtonClass))
#define GTK_IS_OBS_BUTTON(obj)         (G_TYPE_CHECK_INSTANCE_TYPE ((obj), GTK_TYPE_OBS_BUTTON))
#define GTK_IS_OBS_BUTTON_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE ((klass), GTK_TYPE_OBS_BUTTON))
#define OBS_BUTTON_GET_CLASS(obj)      (G_TYPE_INSTANCE_GET_CLASS ((obj), GTK_TYPE_OBS_BUTTON, ObsButtonClass))

#define GTK_INPUT_ERROR -1

typedef struct _ObsButton ObsButton;
typedef struct _ObsButtonClass ObsButtonClass;

GType		obs_button_get_type	   (void) G_GNUC_CONST;
GtkWidget*	obs_button_new		   (GtkAdjustment *adjustment, const DATAINFO *pdinfo);
gdouble		obs_button_get_value       (ObsButton *obs_button);
void		obs_button_set_value	   (ObsButton *obs_button, gdouble value);
void            obs_button_update          (ObsButton *obs_button);

#endif /* OBS_BUTTON_H__ */
