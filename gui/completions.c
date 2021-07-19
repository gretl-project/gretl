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
#include "genparse.h"

#ifdef USE_GTKSOURCEVIEW_3
# define GTK_TYPE_SOURCE_COMPLETION_PROVIDER GTK_SOURCE_TYPE_COMPLETION_PROVIDER
#endif

#define AC_DEBUG 0

#include <gtksourceview/gtksourcecompletionprovider.h>
#include <gtksourceview/gtksourcecompletionitem.h>
#include <gtksourceview/completion-providers/words/gtksourcecompletionwords.h>

int script_auto_complete;

enum {
    PROV_WORDS,
    PROV_FUNCS,
    PROV_CMDS,
    NPROV
};

typedef struct prov_info_ prov_info;

struct prov_info_ {
    gint8 id[NPROV];
    gint8 priority[NPROV];
    const char *name[NPROV];
    GtkTextBuffer *buf[NPROV];
};

static prov_info pi_defaults = {
    { PROV_WORDS, PROV_FUNCS, PROV_CMDS},
    { 1, 2, 3 },
    { "words", "functions", "commands" },
    { NULL, NULL, NULL }
};

static prov_info *prov_info_new (void)
{
    prov_info *pi = malloc(sizeof *pi);

    *pi = pi_defaults;
    return pi;
}

static void prov_info_destroy (prov_info *pi)
{
    int i;

    for (i=0; i<NPROV; i++) {
	if (pi->buf[i] != NULL) {
	    g_object_unref(pi->buf[i]);
	}
    }
    free(pi);
}

/* Create a GtkTextBuffer holding the names of built-in
   gretl functions, to serve as a completion provider.
*/

static GtkTextBuffer *function_names_buffer (void)
{
    GtkTextBuffer *tbuf;
    GString *str;
    gchar *fnames;
    const char *s;
    int i, nf;

    nf = gen_func_count();
    tbuf = gtk_text_buffer_new(NULL);
    str = g_string_sized_new(nf * 8);

    for (i=0; i<nf; i++) {
	s = gen_func_name(i);
	if (*s != '_') {
	    g_string_append(str, s);
	    g_string_append_c(str, ' ');
	}
    }

    fnames = g_string_free(str, FALSE);
    gtk_text_buffer_set_text(tbuf, fnames, -1);
    g_free(fnames);

    return tbuf;
}

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

#ifdef HAVE_ALT_COMPLETE

static void maybe_set_user_activation (GObject *obj)
{
    if (script_auto_complete == COMPLETE_TAB) {
	g_object_set(obj, "activation",
		     GTK_SOURCE_COMPLETION_ACTIVATION_USER_REQUESTED, NULL);
    }
}

#endif

static void add_words_provider (GtkSourceCompletion *comp,
				gint8 id, windata_t *vwin,
				prov_info *pi)
{
    const char *name = pi->name[id];
    GtkSourceCompletionWords *cw;
    GtkTextBuffer *buf;

    cw = gtk_source_completion_words_new(name, NULL);

    g_object_set(cw, "priority", pi->priority[id], NULL);
#ifdef HAVE_ALT_COMPLETE
    maybe_set_user_activation(G_OBJECT(cw));
#endif

    if (id == PROV_WORDS) {
	buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text));
    } else if (id == PROV_FUNCS) {
	buf = pi->buf[id] = function_names_buffer();
    } else {
	buf = pi->buf[id] = command_names_buffer();
    }
    gtk_source_completion_words_register(cw, buf);
    gtk_source_completion_add_provider(comp,
				       GTK_SOURCE_COMPLETION_PROVIDER(cw),
				       NULL);
    g_object_set_data(G_OBJECT(vwin->text), name, cw);
}

static void delete_words_provider (GtkSourceCompletion *comp,
				   gint8 id, windata_t *vwin)
{
    const char *name;
    GtkSourceCompletionWords *cw;
    GtkTextBuffer *buf;
    prov_info *pi;

    pi = g_object_get_data(G_OBJECT(vwin->text), "prov_info");
    if (pi == NULL) {
	return;
    }

    name = pi->name[id];
    cw = g_object_get_data(G_OBJECT(vwin->text), name);
    if (cw == NULL) {
	return;
    }

    if (id == PROV_WORDS) {
	buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text));
    } else {
	buf = pi->buf[id];
    }
    gtk_source_completion_words_unregister(cw, buf);
    gtk_source_completion_remove_provider(comp,
					  GTK_SOURCE_COMPLETION_PROVIDER(cw),
					  NULL);
    g_object_unref(G_OBJECT(cw));
    g_object_set_data(G_OBJECT(vwin->text), name, NULL);

    if (id != PROV_WORDS) {
	g_object_unref(G_OBJECT(pi->buf[id]));
	pi->buf[id] = NULL;
    }
}

/* apparatus for providing "snippets" which may consist of
   several "words"
*/

typedef struct _SnippetProvider SnippetProvider;
typedef struct _SnippetProviderClass SnippetProviderClass;

struct _SnippetProvider
{
    GObject parent;
    GList *proposals;
    gint priority;
    gchar *name;
    GdkPixbuf *icon;
    GtkSourceCompletionActivation activation;
};

struct _SnippetProviderClass
{
    GObjectClass parent_class;
};

static void snippet_provider_iface_init (GtkSourceCompletionProviderIface *iface);
GType snippet_provider_get_type (void);

G_DEFINE_TYPE_WITH_CODE (SnippetProvider,
			 snippet_provider,
			 G_TYPE_OBJECT,
			 G_IMPLEMENT_INTERFACE(GTK_TYPE_SOURCE_COMPLETION_PROVIDER,
					       snippet_provider_iface_init))

static gchar *
snippet_provider_get_name (GtkSourceCompletionProvider *provider)
{
    return g_strdup(((SnippetProvider *) provider)->name);
}

static gint
snippet_provider_get_priority (GtkSourceCompletionProvider *provider)
{
    return ((SnippetProvider *) provider)->priority;
}

static GtkSourceCompletionActivation
snippet_provider_get_activation (GtkSourceCompletionProvider *provider)
{
    return ((SnippetProvider *) provider)->activation;
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

static int proposal_get_cursor_offset (const gchar *s)
{
    if (!strncmp(s, "if", 2)) {
	return 3;
    } else if (!strncmp(s, "loop", 4)) {
	return 5;
    } else if (!strncmp(s, "function", 8)) {
	return 9;
    } else {
	return 0;
    }
}

/* Back up over the trigger for completion; insert the replacement
   text; then back the cursor up to the point where the user will
   first have to add something to the boilerplate.
*/

static gboolean
snippet_activate_proposal (GtkSourceCompletionProvider *provider,
			   GtkSourceCompletionProposal *proposal,
			   GtkTextIter *iter)
{
    gchar *s = gtk_source_completion_proposal_get_text(proposal);
    int n = proposal_get_cursor_offset(s);

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
snippet_provider_match (GtkSourceCompletionProvider *provider,
			GtkSourceCompletionContext *context)
{
    return TRUE;
}

static void
snippet_provider_populate (GtkSourceCompletionProvider *provider,
			   GtkSourceCompletionContext *context)
{
    GList *L = ((SnippetProvider *) provider)->proposals;
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

static void
snippet_provider_iface_init (GtkSourceCompletionProviderIface *iface)
{
    iface->get_name = snippet_provider_get_name;
    iface->populate = snippet_provider_populate;
    iface->match = snippet_provider_match;
    iface->get_priority = snippet_provider_get_priority;
    iface->get_activation = snippet_provider_get_activation;
    iface->activate_proposal = snippet_activate_proposal;
}

static void snippet_provider_class_init (SnippetProviderClass *klass)
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
    { "function...", "function type name ()\n\t\nend function\n\n" }
};

static void snippet_provider_init (SnippetProvider *self)
{
    GtkSourceCompletionItem *item;
    GList *proposals = NULL;
    int i, n = G_N_ELEMENTS(snippets);

    for (i=0; i<n; i++) {
	item = gtk_source_completion_item_new(snippets[i].label,
					      snippets[i].text,
					      NULL, NULL);
	proposals = g_list_prepend(proposals, item);
    }
    self->proposals = proposals;
}

static void add_snippets_provider (GtkSourceCompletion *comp,
				   windata_t *vwin)
{
    SnippetProvider *sp;

    sp = g_object_new(snippet_provider_get_type(), NULL);
    sp->priority = 1;
#ifdef HAVE_ALT_COMPLETE
    if (script_auto_complete == COMPLETE_TAB) {
	sp->activation = GTK_SOURCE_COMPLETION_ACTIVATION_USER_REQUESTED;
    } else {
	sp->activation = GTK_SOURCE_COMPLETION_ACTIVATION_INTERACTIVE;
    }
#else
    sp->activation = GTK_SOURCE_COMPLETION_ACTIVATION_INTERACTIVE;
#endif
    sp->name = "snippets";
    gtk_source_completion_add_provider(comp,
				       GTK_SOURCE_COMPLETION_PROVIDER(sp),
				       NULL);
    g_object_set_data(G_OBJECT(vwin->text), sp->name, sp);
}

static void delete_snippets_provider (GtkSourceCompletion *comp,
				      windata_t *vwin)
{
    SnippetProvider *sp = g_object_get_data(G_OBJECT(vwin->text),
					    "snippets");

    if (sp == NULL) {
	return;
    }

    gtk_source_completion_remove_provider(comp,
					  GTK_SOURCE_COMPLETION_PROVIDER(sp),
					  NULL);
    g_object_unref(G_OBJECT(sp));
}

/* end snippets apparatus */

static void destroy_providers (GtkWidget *w, gpointer p)
{
    prov_info *pi = g_object_get_data(G_OBJECT(w), "prov_info");

    if (pi != NULL) {
	prov_info_destroy(pi);
    }
}

void set_sv_auto_completion (windata_t *vwin)
{
    GtkSourceCompletion *comp;
    GList *L;

    comp = gtk_source_view_get_completion(GTK_SOURCE_VIEW(vwin->text));
    L = gtk_source_completion_get_providers(comp);

#if 0
    fprintf(stderr, "set_sv_auto_completion: comp %p, L %p\n",
	    (void *) comp, (void *) L);
#endif

    if (script_auto_complete && L == NULL) {
	/* set up and activate */
	prov_info *pi = prov_info_new();

	g_object_set_data(G_OBJECT(vwin->text), "prov_info", pi);
	add_words_provider(comp, PROV_WORDS, vwin, pi);
	add_words_provider(comp, PROV_FUNCS, vwin, pi);
	if (vwin->role == CONSOLE) {
	    add_words_provider(comp, PROV_CMDS, vwin, pi);
	} else {
	    add_snippets_provider(comp, vwin);
	}
	g_signal_connect(G_OBJECT(vwin->text), "destroy",
			 G_CALLBACK(destroy_providers), NULL);
    } else if (!script_auto_complete && L != NULL) {
	/* de-activate and clean up */
	delete_words_provider(comp, PROV_WORDS, vwin);
	delete_words_provider(comp, PROV_FUNCS, vwin);
	delete_words_provider(comp, PROV_CMDS, vwin);
	delete_snippets_provider(comp, vwin);
    }
}

/* Note: @w must be a vwin->text instance, and @priorities
   must contain 4 small integers, high values representing
   high priority. The order of the values is PROV_WORDS,
   PROV_FUNCS, PROV_CMDS, PROV_SNIPPETS.
*/

void set_auto_completion_priority (GtkWidget *w, gint8 *order)
{
    GtkSourceCompletion *comp;
    GList *L;
    int i = 0;

    comp = gtk_source_view_get_completion(GTK_SOURCE_VIEW(w));
    L = gtk_source_completion_get_providers(comp);

    while (L != NULL) {
	g_object_set(G_OBJECT(L->data), "priority", order[i++], NULL);
	L = L->next;
    }
}

void tab_auto_complete (GtkWidget *w)
{
    g_signal_emit_by_name(GTK_SOURCE_VIEW(w), "show-completion", NULL);
}
