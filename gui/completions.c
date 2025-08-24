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
#include "gen_public.h"
#include "completions.h"

#if GTKSOURCEVIEW_VERSION > 2
# define GTK_TYPE_SOURCE_COMPLETION_PROVIDER GTK_SOURCE_TYPE_COMPLETION_PROVIDER
#endif

#define AC_DEBUG 0

#include <gtksourceview/gtksourcecompletionprovider.h>
#include <gtksourceview/gtksourcecompletionitem.h>
#include <gtksourceview/completion-providers/words/gtksourcecompletionwords.h>

/* global, referenced in settings.c and elsewhere */
int hansl_completion;
int console_completion;

enum {
    PROV_WORDS,
    PROV_CMDS,
    PROV_FUNCS,
    PROV_SNIPPETS,
    PROV_SERIES,
    N_PROV
};

static const char *prov_names[] = {
    "words", "commands", "functions", "snippets", "series"
};

typedef struct prov_info_ prov_info;

struct prov_info_ {
    void *ptr;          /* pointer to the completer struct */
    GtkTextBuffer *buf; /* static buffer used by PROV_CMDS */
};

/* allocate the required number of "provider" structs */

static prov_info *prov_info_new (void)
{
    prov_info *pi = malloc(N_PROV * sizeof *pi);
    int i;

    for (i=0; i<N_PROV; i++) {
	pi[i].ptr = NULL;
	pi[i].buf = NULL;
    }

    return pi;
}

static void prov_info_destroy (prov_info *pi)
{
    int i;

    for (i=0; i<N_PROV; i++) {
	if (pi[i].buf != NULL) {
	    g_object_unref(pi[i].buf);
	}
    }
    free(pi);
}

static void destroy_completion_providers (GtkWidget *w, gpointer p)
{
    prov_info *pi = g_object_get_data(G_OBJECT(w), "prov_info");

    if (pi != NULL) {
	prov_info_destroy(pi);
	g_object_set_data(G_OBJECT(w), "prov_info", NULL);
    }
}

static gchar const **func_names;
static int n_func_names;

/* Create an array of function names for use with the functions
   completion provider */

static void make_function_names_array (void)
{
    const char *s;
    int nf = gen_func_count();
    int i, j, nn = 0;

    for (i=0; i<nf; i++) {
	s = gen_func_name(i);
	if (*s != '_') {
            /* exclude "hidden" functions */
            nn++;
        }
    }

    func_names = calloc(nn, sizeof *func_names);
    n_func_names = nn;

    for (i=0, j=0; i<nf; i++) {
	s = gen_func_name(i);
	if (*s != '_') {
            func_names[j++] = s;
        }
    }

#if AC_DEBUG
    fprintf(stderr, "make_function_names_array: got %d names\n", nn);
#endif
}

/* Create a GtkTextBuffer holding the names of gretl commands to serve
   as a completion provider.
*/

static GtkTextBuffer *command_names_buffer (void)
{
    GtkTextBuffer *tbuf;
    GString *str;
    gchar *cnames;
    int i;

    tbuf = gtk_text_buffer_new(NULL);
    str = g_string_sized_new(NC * 8);

    for (i=1; i<=NC; i++) {
	g_string_append(str, gretl_command_word(i));
	g_string_append_c(str, ' ');
    }

    cnames = g_string_free(str, FALSE);
    gtk_text_buffer_set_text(tbuf, cnames, -1);
    g_free(cnames);

    return tbuf;
}

/* Apparatus for providing "snippets" which may consist of several
   "words", and other gretl-specific material.
*/

typedef struct _GretlProvider GretlProvider;
typedef struct _GretlProviderClass GretlProviderClass;

struct _GretlProvider {
    GObject parent;
    GList *proposals;
    gint priority;
    const gchar *name;
    GdkPixbuf *icon;
    GtkSourceCompletionActivation activation;
    gint id;
};

struct _GretlProviderClass {
    GObjectClass parent_class;
};

static void gretl_provider_iface_init (GtkSourceCompletionProviderIface *iface);
GType gretl_provider_get_type (void);

G_DEFINE_TYPE_WITH_CODE (GretlProvider,
			 gretl_provider,
			 G_TYPE_OBJECT,
			 G_IMPLEMENT_INTERFACE(GTK_TYPE_SOURCE_COMPLETION_PROVIDER,
					       gretl_provider_iface_init))

static gchar *
gretl_provider_get_name (GtkSourceCompletionProvider *provider)
{
    return g_strdup(((GretlProvider *) provider)->name);
}

static gint
gretl_provider_get_priority (GtkSourceCompletionProvider *provider)
{
    return ((GretlProvider *) provider)->priority;
}

static GtkSourceCompletionActivation
gretl_provider_get_activation (GtkSourceCompletionProvider *provider)
{
    return ((GretlProvider *) provider)->activation;
}

#define valid(c) (c == '_' || isalnum(c))

static gboolean backward_word_start (GtkTextIter *iter)
{
    GtkTextIter prev = *iter;

    while (TRUE) {
	/* starting a line is OK */
	if (gtk_text_iter_starts_line(&prev)) {
	    break;
	}
	gtk_text_iter_backward_char(&prev);
	/* check if the previous character is a valid word character */
	if (!valid(gtk_text_iter_get_char(&prev))) {
	    break;
	}
	*iter = prev;
    }

    if (!valid(gtk_text_iter_get_char(iter))) {
	return FALSE;
    }

    return isalpha(gtk_text_iter_get_char(iter));
}

static gboolean forward_word_end (GtkTextIter *iter)
{
    while (TRUE) {
	/* ending a line is OK */
	if (gtk_text_iter_ends_line(iter)) {
	    break;
	}
	/* check if the next character is a valid word character */
	if (!valid(gtk_text_iter_get_char(iter))) {
	    break;
	}
	gtk_text_iter_forward_char(iter);
    }

    return TRUE;
}

static gchar *get_word_at_iter (GtkTextIter *iter)
{
    GtkTextIter end = *iter;

    if (!forward_word_end(iter) || !gtk_text_iter_equal(iter, &end)) {
	return NULL;
    } else if (!backward_word_start(iter)) {
	return NULL;
    } else if (gtk_text_iter_equal(iter, &end)) {
	return NULL;
    } else {
	return gtk_text_iter_get_text(iter, &end);
    }
}

static int snippet_proposal_get_cursor_offset (const gchar *s)
{
    if (!strncmp(s, "if", 2)) {
	return 3;
    } else if (!strncmp(s, "loop", 4)) {
	return 5;
    } else if (!strncmp(s, "function", 8)) {
	return 9;
    } else if (!strncmp(s, "outfile", 7)) {
	return 8;
    } else if (!strncmp(s, "plot", 4)) {
	return 5;
    } else if (!strncmp(s, "mpi", 3)) {
	return 5;
    } else {
	return 0;
    }
}

/* Back up over the trigger for completion; insert the replacement text;
   then back the cursor up to the point where the user will first have
   to add something to the boilerplate.
*/

static gboolean
snippet_activate_proposal (GtkSourceCompletionProvider *provider,
			   GtkSourceCompletionProposal *proposal,
			   GtkTextIter *iter)
{
    gchar *s = gtk_source_completion_proposal_get_text(proposal);
    int n = snippet_proposal_get_cursor_offset(s);

    if (n > 0) {
	GtkTextBuffer *buf = gtk_text_iter_get_buffer(iter);
	GtkTextIter start = *iter;

	backward_word_start(&start);
	gtk_text_buffer_delete(buf, &start, iter);
	gtk_text_buffer_insert(buf, iter, s, -1);
	gtk_text_iter_backward_chars(iter, strlen(s) - n);
	gtk_text_buffer_place_cursor(buf, iter);
	return TRUE;
    } else {
	return FALSE;
    }
}

static gboolean
func_activate_proposal (GtkSourceCompletionProvider *provider,
                        GtkSourceCompletionProposal *proposal,
                        GtkTextIter *iter)
{
    gchar *s = gtk_source_completion_proposal_get_text(proposal);
    GtkTextBuffer *buf = gtk_text_iter_get_buffer(iter);
    GtkTextIter start = *iter;

    backward_word_start(&start);
    gtk_text_buffer_delete(buf, &start, iter);
    gtk_text_buffer_insert(buf, iter, s, -1);
    gtk_text_buffer_insert(buf, iter, "()", -1);
    gtk_text_iter_backward_chars(iter, 1);
    gtk_text_buffer_place_cursor(buf, iter);

    return TRUE;
}

static gboolean
series_activate_proposal (GtkSourceCompletionProvider *provider,
			  GtkSourceCompletionProposal *proposal,
			  GtkTextIter *iter)
{
    /* just accept the gtksourceview default */
    return FALSE;
}

static gboolean
gretl_activate_proposal (GtkSourceCompletionProvider *provider,
			 GtkSourceCompletionProposal *proposal,
			 GtkTextIter *iter)
{
    gint id = ((GretlProvider *) provider)->id;
    gboolean ret = FALSE;

#if AC_DEBUG
    fprintf(stderr, "gretl_activate_proposal (%s)\n", prov_names[id]);
#endif

    if (id == PROV_SNIPPETS) {
	ret = snippet_activate_proposal(provider, proposal, iter);
    } else if (id == PROV_FUNCS) {
        ret = func_activate_proposal(provider, proposal, iter);
    } else if (id == PROV_SERIES) {
	ret = series_activate_proposal(provider, proposal, iter);
    }

    return ret;
}

static gboolean
gretl_provider_match (GtkSourceCompletionProvider *provider,
		      GtkSourceCompletionContext *context)
{
    return TRUE;
}

static void
snippet_provider_populate (GtkSourceCompletionProvider *provider,
			   GtkSourceCompletionContext *context)
{
    GList *L = ((GretlProvider *) provider)->proposals;
    GList *ret = NULL;
    GtkTextIter iter;
    gchar *word;
    int n;

    gtk_source_completion_context_get_iter(context, &iter);
    word = get_word_at_iter(&iter);

    if (word != NULL && (n = strlen(word)) >= 2) {
	GtkSourceCompletionItem *item;
	gchar *label;

	while (L != NULL) {
	    item = L->data;
	    g_object_get(item, "label", &label, NULL);
	    if (!strncmp(label, word, n)) {
		ret = g_list_prepend(ret, item);
	    }
	    g_free(label);
	    L = L->next;
	}
	g_free(word);
    }

    if (ret != NULL) {
	ret = g_list_reverse(ret);
    }

    gtk_source_completion_context_add_proposals(context, provider, ret, TRUE);
    g_list_free(ret);
}

static GtkSourceCompletionItem *comp_item_new (const gchar *label,
					       const gchar *text)
{
    GtkSourceCompletionItem *item;

#if GTKSOURCEVIEW_VERSION == 4
    item = gtk_source_completion_item_new();
    gtk_source_completion_item_set_label(item, label);
    gtk_source_completion_item_set_text(item, text);
#else
    item = gtk_source_completion_item_new(label, text, NULL, NULL);
#endif

    return item;
}

static void
func_provider_populate (GtkSourceCompletionProvider *provider,
                        GtkSourceCompletionContext *context)
{
    GList *ret = NULL;
    GtkTextIter iter;
    gchar *word;
    int n;

    gtk_source_completion_context_get_iter(context, &iter);
    word = get_word_at_iter(&iter);

    if (word != NULL && (n = strlen(word)) > 1) {
	GtkSourceCompletionItem *item;
	const char *fname;
	int i;

	for (i=0; i<n_func_names; i++) {
	    fname = func_names[i];
	    if (!strncmp(fname, word, n)) {
		item = comp_item_new(fname, fname);
		ret = g_list_prepend(ret, item);
	    }
	}
    }

    if (ret != NULL) {
	ret = g_list_reverse(ret);
    }

    gtk_source_completion_context_add_proposals(context, provider, ret, TRUE);
    g_list_free(ret);
}

#ifndef GRETL_EDIT

static void
series_provider_populate (GtkSourceCompletionProvider *provider,
			  GtkSourceCompletionContext *context)
{
    GList *ret = NULL;
    GtkTextIter iter;
    gchar *word;
    int n;

    gtk_source_completion_context_get_iter(context, &iter);
    word = get_word_at_iter(&iter);

    if (word != NULL && (n = strlen(word)) > 1) {
	GtkSourceCompletionItem *item;
	const char *vname;
	int i;

	for (i=0; i<dataset->v; i++) {
	    vname = dataset->varname[i];
	    if (!strncmp(vname, word, n)) {
		item = comp_item_new(vname, vname);
		ret = g_list_prepend(ret, item);
	    }
	}
    }

    if (ret != NULL) {
	ret = g_list_reverse(ret);
    }

    gtk_source_completion_context_add_proposals(context, provider, ret, TRUE);
    g_list_free(ret);
}

#endif

static void
gretl_provider_populate (GtkSourceCompletionProvider *provider,
			 GtkSourceCompletionContext *context)
{
    gint id = ((GretlProvider *) provider)->id;

    if (id == PROV_SNIPPETS) {
	snippet_provider_populate(provider, context);
    } else if (id == PROV_FUNCS) {
        func_provider_populate(provider, context);
    }
#ifndef GRETL_EDIT
    else if (id == PROV_SERIES) {
	series_provider_populate(provider, context);
    }
#endif
}

static void
gretl_provider_iface_init (GtkSourceCompletionProviderIface *iface)
{
    iface->get_name = gretl_provider_get_name;
    iface->populate = gretl_provider_populate;
    iface->match = gretl_provider_match;
    iface->get_priority = gretl_provider_get_priority;
    iface->get_activation = gretl_provider_get_activation;
    iface->activate_proposal = gretl_activate_proposal;
}

static void
gretl_provider_set_activation (GretlProvider *gp,
			       GtkSourceCompletionActivation A)
{
    gp->activation = A;
}

static void
words_provider_set_activation (GObject *obj,
			       GtkSourceCompletionActivation A)
{
    g_object_set(obj, "activation", A, NULL);
}

static void providers_set_activation (prov_info *pi,
				      int userval)
{
    GtkSourceCompletionActivation A =
	GTK_SOURCE_COMPLETION_ACTIVATION_NONE;
    int i;

    if (userval == COMPLETE_USER) {
	A = GTK_SOURCE_COMPLETION_ACTIVATION_USER_REQUESTED;
    } else if (userval == COMPLETE_AUTO) {
	A = GTK_SOURCE_COMPLETION_ACTIVATION_INTERACTIVE;
    }

    for (i=0; i<N_PROV; i++) {
	if (pi[i].ptr != NULL) {
#if AC_DEBUG
	    fprintf(stderr, "set activation %d on %s\n", A, prov_names[i]);
#endif
	    if (i >= PROV_FUNCS) {
		gretl_provider_set_activation(pi[i].ptr, A);
	    } else {
		words_provider_set_activation(pi[i].ptr, A);
	    }
	}
    }
}

static void gretl_provider_class_init (GretlProviderClass *klass)
{
    return;
}

struct snippet {
    const char *label;
    const char *text;
};

/* a few simple example snippets for now */

struct snippet snippets[] = {
    { "loop..endloop", "loop \n\t\nendloop\n" },
    { "if..endif", "if \n\t\nendif\n" },
    { "function...", "function type name ()\n\t\nend function\n\n" },
    { "outfile...", "outfile \n\t\nend outfile\n" },
    { "plot...", "plot \n\t\nend plot\n" },
    { "mpi...", "mpi\n\t\nend mpi\n" }
};

static void snippet_provider_init (GretlProvider *self)
{
    GtkSourceCompletionItem *item;
    GList *proposals = NULL;
    int i, n = G_N_ELEMENTS(snippets);

    for (i=0; i<n; i++) {
	item = comp_item_new(snippets[i].label, snippets[i].text);
	proposals = g_list_prepend(proposals, item);
    }
    self->proposals = proposals;
}

static void func_provider_init (GretlProvider *self)
{
    GtkSourceCompletionItem *item;
    GList *proposals = NULL;
    int i;

    if (n_func_names == 0) {
        /* create static array */
        make_function_names_array();
    }

    for (i=0; i<n_func_names; i++) {
	item = comp_item_new(func_names[i], func_names[i]);
	proposals = g_list_prepend(proposals, item);
    }
    self->proposals = proposals;
}

static int gretl_prov_id;

static void gretl_provider_init (GretlProvider *self)
{
    if (gretl_prov_id == PROV_SNIPPETS) {
	snippet_provider_init(self);
    } else if (gretl_prov_id == PROV_FUNCS) {
        func_provider_init(self);
    } else {
	self->proposals = NULL;
    }
}

#if AC_DEBUG

static void notify_activated (GtkSourceCompletion *comp,
			      gpointer p)
{
    fprintf(stderr, "+++ got activate-proposal signal +++\n");
}

static void notify_hidden (GtkSourceCompletion *comp,
			   gpointer p)
{
    fprintf(stderr, "+++ got hide signal +++\n");
}

static void notify_populate (GtkSourceCompletion *comp,
			     GtkSourceCompletionContext *context,
			     gpointer p)
{
    fprintf(stderr, "+++ got populate-context signal +++\n");
}

#endif

static void add_gretl_provider (GtkSourceCompletion *comp,
				gint id, gint priority,
				prov_info *pi)
{
    GretlProvider *gp;

    gretl_prov_id = id; /* hack! */
    gp = g_object_new(gretl_provider_get_type(), NULL);
    pi[id].ptr = gp;
    gp->id = id;
    gp->priority = priority;
    gp->name = prov_names[id];
    gtk_source_completion_add_provider(comp,
				       GTK_SOURCE_COMPLETION_PROVIDER(gp),
				       NULL);
    g_object_unref(gp);
}

/* end gretl-specific provider apparatus */

static void add_words_provider (GtkSourceCompletion *comp,
				gint8 id, gint priority,
				GtkWidget *w, prov_info *pi)
{
    const char *name = prov_names[id];
    GtkSourceCompletionWords *cw;
    GtkTextBuffer *buf;

    cw = gtk_source_completion_words_new(name, NULL);
    pi[id].ptr = cw;
    g_object_set(cw, "priority", priority, NULL);

    if (id == PROV_CMDS) {
	buf = pi[id].buf = command_names_buffer();
    } else {
	/* plain PROV_WORDS */
	buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(w));
    }
    gtk_source_completion_words_register(cw, buf);
    gtk_source_completion_add_provider(comp,
				       GTK_SOURCE_COMPLETION_PROVIDER(cw),
				       NULL);
    g_object_set_data(G_OBJECT(w), name, cw);
    g_object_unref(cw);
}

void set_sv_completion (GtkWidget *text, int role)
{
    GtkSourceCompletion *comp;
    prov_info *pi = NULL;
    int compval;

    comp = gtk_source_view_get_completion(GTK_SOURCE_VIEW(text));
    pi = g_object_get_data(G_OBJECT(text), "prov_info");
    compval = (role == CONSOLE)? console_completion :
	hansl_completion;

#if AC_DEBUG
    fprintf(stderr, "set_sv_completion: comp %s, prov_info %s, completion %d\n",
	    comp==NULL ? "null" : "present", pi==NULL? "null" : "present",
	    compval);
#endif

    if (compval && pi == NULL) {
	/* set up and activate */
	g_object_set(G_OBJECT(comp), "accelerators", 10,
		     "remember-info-visibility", TRUE, NULL);
	pi = prov_info_new();
	g_object_set_data(G_OBJECT(text), "prov_info", pi);
	if (role == CONSOLE) {
	    add_words_provider(comp, PROV_CMDS,   4, text, pi);
	    add_gretl_provider(comp, PROV_SERIES, 3, pi);
	    add_gretl_provider(comp, PROV_FUNCS,  2, pi);
	    add_words_provider(comp, PROV_WORDS,  1, text, pi);
	} else {
	    /* context is script editor */
#ifdef GRETL_EDIT
	    add_gretl_provider(comp, PROV_SNIPPETS, 4, pi);
	    add_words_provider(comp, PROV_CMDS,     3, text, pi);
	    add_gretl_provider(comp, PROV_FUNCS,    2, pi);
	    add_words_provider(comp, PROV_WORDS,    1, text, pi);
#else
	    add_gretl_provider(comp, PROV_SNIPPETS, 5, pi);
	    add_words_provider(comp, PROV_CMDS,     4, text, pi);
	    add_gretl_provider(comp, PROV_FUNCS,    3, pi);
	    add_gretl_provider(comp, PROV_SERIES,   2, pi);
	    add_words_provider(comp, PROV_WORDS,    1, text, pi);
#endif
	}
	g_signal_connect(G_OBJECT(text), "destroy",
			 G_CALLBACK(destroy_completion_providers), NULL);
#if AC_DEBUG
	g_signal_connect(G_OBJECT(comp), "activate-proposal",
			 G_CALLBACK(notify_activated), NULL);
	g_signal_connect(G_OBJECT(comp), "hide",
			 G_CALLBACK(notify_hidden), NULL);
	g_signal_connect(G_OBJECT(comp), "populate-context",
			 G_CALLBACK(notify_populate), NULL);
	fprintf(stderr, "providers set up\n");
#endif
    }

    if (pi != NULL) {
	providers_set_activation(pi, compval);
    }
}

void call_user_completion (GtkWidget *w)
{
#if AC_DEBUG
    fprintf(stderr, "call_user_completion...\n");
#endif
    g_signal_emit_by_name(GTK_SOURCE_VIEW(w), "show-completion", NULL);
}
