/* gtkplot - 2d scientific plots widget for gtk+
 * Copyright 1999-2001  Adrian E. Feiguin <feiguin@ifir.edu.ar>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef __GTK_PLOT_H__
#define __GTK_PLOT_H__

#include "gtkplotpc.h"
#include "gtkplotps.h"

gboolean psinit (GtkPlotPC *pc);
void psleave (GtkPlotPC *pc);

typedef enum {
    GTK_PLOT_BORDER_NONE,
    GTK_PLOT_BORDER_LINE,
    GTK_PLOT_BORDER_SHADOW
} GtkPlotBorderStyle;

#endif /* __GTK_PLOT_H__ */
