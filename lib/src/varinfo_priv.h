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

/* Warning: as things stand (2021-06-21) this struct must not be
   modified, on pain of breaking the currently favoured means of
   saving gretl data in binary form -- for which see purebin.c
   in the plugin directory of the source tree. The "purebin"
   reader works on the assumption that the size and membership
   of this struct are known; that will be subverted if the struct
   has been changed since data were saved as gdtb, and opening
   a file containing an out-of-date VARINFO will almost surely
   result in a segfault.
*/

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
};

# if defined(G_OS_WIN32) && !defined(_WIN64)

/* In the writing of gdtb files the VARINFO struct must be a fixed
   size, invariant with respect to word length. But because of the
   @label and @st pointers -- which are always NULL in a gtdb file,
   but which differ in size -- the respective "native" structs are not
   the same size. We let the burden of adjustment fall on 32-bit
   builds: when writing to gdtb we transcribe from VARINFO to
   VARINFO64, which includes emulated 64-bit pointers, and when
   reading we transcribe back to the small-pointer variant. See
   plugin/purebin/c for details.
*/

struct VARINFO64 {
    guint64 p1; /* fake 64-bit pointer */
    char display_name[MAXDISP];
    char parent[VNAMELEN];
    VarFlags flags;
    char compact_method;
    gint64 mtime;
    short transform;
    short lag;
    short stack_level;
    short midas_period;
    char midas_freq;
    short orig_pd;
    guint64 p2; /* fake 64-bit pointer */
};

# endif

#endif /* VARINFO_PRIV_H */
