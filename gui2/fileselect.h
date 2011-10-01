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

typedef enum {
    FSEL_DATA_NONE = 0,
    FSEL_DATA_PRN,       /* a text save action: data in the form of a PRN */
    FSEL_DATA_VWIN,      /* action pertaining to a specific viewer window */
    FSEL_DATA_MISC,      /* the file-selector call carries some misc data */
    FSEL_DATA_STATUS     /* provides a means of returning action status */
} FselDataSrc;

void file_selector (int action, FselDataSrc src, gpointer data);

void file_selector_with_parent (int action, FselDataSrc src, 
				gpointer data, GtkWidget *w);
