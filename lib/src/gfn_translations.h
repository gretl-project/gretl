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

typedef struct Translations_ Translations;

void destroy_translations (Translations *T);

const char *get_gfn_translation (Translations *T,
                                 const char *id,
                                 const char *lang);

Translations *read_translations_file (const char *fname,
                                      int *err);

Translations *read_translations_element (xmlNodePtr root,
                                         xmlDocPtr doc);

void write_translations (Translations *T, PRN *prn);

#endif /* GFN_TRANSLATIONS_H */
