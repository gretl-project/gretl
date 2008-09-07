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

/*  guiprint.h for gretl */ 

#ifndef GUIPRINT_H
#define GUIPRINT_H

#ifdef G_OS_WIN32
int winprint_graph (char *emfname);
#endif

#ifdef NATIVE_PRINTING
void winprint (char *fullbuf, char *selbuf);
# ifndef G_OS_WIN32
void gtk_print_graph (const char *fname);
# endif
#endif

void special_print_summary (const Summary *summ,
			    const DATAINFO *pdinfo,
			    PRN *prn);

void special_print_vmatrix (const VMatrix *corr,
			    const DATAINFO *pdinfo, 
			    PRN *prn);

void special_print_fit_resid (const FITRESID *fr, 
			      const DATAINFO *pdinfo, 
			      PRN *prn);

void special_print_forecast (const FITRESID *fr, 
			     const DATAINFO *pdinfo, 
			     PRN *prn);

void special_print_confints (const CoeffIntervals *cf, 
			     PRN *prn);

int csv_to_clipboard (void);

int csv_selected_to_clipboard (void);

int csv_copy_listed_vars (windata_t *vwin, int fmt, int action);

int font_has_minus (PangoFontDescription *desc);

#endif /* GUIPRINT_H */
