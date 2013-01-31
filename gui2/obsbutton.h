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

#ifndef OBS_BUTTON_H__
#define OBS_BUTTON_H__

typedef enum {
    OBS_BUTTON_NONE,
    OBS_BUTTON_T1,
    OBS_BUTTON_T2
} ObsButtonRole;

GtkWidget *obs_button_new (GtkAdjustment *adjustment, DATASET *dset,
			   ObsButtonRole role);

GtkWidget *data_start_button (GtkAdjustment *adj, DATASET *dset);

int obs_button_get_value (GtkWidget *button);

const gchar *obs_button_get_string (GtkWidget *button);

void obs_button_set_partner (GtkWidget *button, GtkWidget *partner);

#endif /* OBS_BUTTON_H__ */
