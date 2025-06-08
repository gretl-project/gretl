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

#ifndef GFN_TRANSLATIONS_H
#define GFN_TRANSLATIONS_H

typedef struct Translation_ Translation;

void destroy_translation (Translation *T);

const char *get_gfn_translation (Translation *T,
                                 const char *id);

Translation *read_translations_file (const char *fname,
                                     int *err);

Translation *read_translation_element (xmlNodePtr root,
                                       xmlDocPtr doc);

Translation *update_translation (Translation *T0,
                                 const char *trbuf);

void write_translation (Translation *T, PRN *prn);

gchar *get_translatable_content (const char **ps);

#endif /* GFN_TRANSLATIONS_H */
