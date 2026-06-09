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

#include "gretl.h"

/* Apparatus to support multi-level "Go back" references in a
   GtkTextView object
*/

typedef struct backrefs_ {
    GtkTextMark **marks;
    int n_marks;
    int n_slots;
} backrefs;

/* Allocate a new backrefs struct. */

backrefs *backrefs_new (void)
{
    backrefs *refs = malloc(sizeof *backrefs);

    fprintf(stderr, "backrefs_new() was called\n");

    refs->marks = NULL;
    refs->n_marks = 0;
    refs->n_slots = 0;

    return refs;
}

/* Free a backrefs struct when its parent widget is destroyed. */

void backrefs_destroy (backrefs *refs)
{
    if (refs != NULL) {
	fprintf(stderr, "backrefs_destroy() was called\n");
	free(refs->marks);
	free(refs);
    }
}

/* Push @mark onto a stack of backward references for @w. */

void push_backref (GtkWidget *w, GtkTextMark *mark)
{
    backrefs *refs = g_object_get_data(G_OBJECT(w), "backrefs");

    fprintf(stderr, "push_backref() was called (refs %p)\n",
	    (void *) refs);

    if (refs == NULL) {
	refs = backrefs_new();
    }

    if (refs->n_marks == refs->n_slots) {
	/* @refs is full already */
	int ns = refs->n_slots + 1;

	refs->marks = realloc(refs->marks, ns * sizeof *refs->marks);
	refs->n_slots = ns;
    }

    refs->marks[refs->n_marks] = mark;
    refs->n_marks += 1;
    fprintf(stderr, "  refs->nmarks is now %d\n", refs->n_marks);
}

/* Pop a GtkTextMark off the stack of backward references for @w, if
   such a stack exists and is not empty.
*/

GtkTextMark *pop_backref (GtkWidget *w)
{
    backrefs *refs = g_object_get_data(G_OBJECT(w), "backrefs");
    GtkTextMark *ret = NULL;

    fprintf(stderr, "pop_backref() was called (refs %p)\n",
	    (void *) refs);

    if (refs != NULL && refs->n_marks > 0) {
	int n = refs->n_marks - 1;
	
	ret = ref->marks[n];
	refs->n_marks = n;
    }

    fprintf(stderr, "  refs->nmarks is now %d\n", refs->n_marks);

    return ret;
}
