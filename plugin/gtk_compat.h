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

#ifndef GTK_COMPAT_H
#define GTK_COMPAT_H

#if (GTK_MAJOR_VERSION == 2 && GTK_MINOR_VERSION < 14)
# define gtk_dialog_get_content_area(d) (d->vbox)
# define gtk_dialog_get_action_area(d) (d->action_area)
#endif

#if (GTK_MAJOR_VERSION == 2 && GTK_MINOR_VERSION < 18)
# define gtk_widget_is_sensitive(w) GTK_WIDGET_IS_SENSITIVE(w)
# define gtk_widget_has_focus(w) GTK_WIDGET_HAS_FOCUS(w)
# define gtk_widget_set_can_default(w,s) GTK_WIDGET_SET_FLAGS(w, GTK_CAN_DEFAULT)
# define gtk_widget_set_can_focus(w,s) GTK_WIDGET_SET_FLAGS(w, GTK_CAN_FOCUS)
#endif

#endif /* GTK_COMPAT_H */
