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

/* struct definition: basically private but needs to be shared between
   lib/src/dataset.c and plugin/purebin.c
*/

#ifndef VARINFO_PRIV_H
#define VARINFO_PRIV_H

struct VARINFO_ {
    char *label;
    char display_name[MAXDISP];
    char parent[VNAMELEN];
    VarFlags flags;
    char compact_method;
    gint64 mtime;
    short transform;    /* note: command index of transform */
    short lag;
    short stack_level;
    short midas_period;
    char midas_freq;
    short orig_pd;
    series_table *st;
    const char *list_parent;
};

#endif /* VARINFO_PRIV_H */
