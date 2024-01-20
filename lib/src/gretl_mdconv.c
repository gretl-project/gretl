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

/* Conversion from simple markdown to gretl's custom help-file
   markup. The regular expressions used below are adapted from
   Uwe Jugel's md2pango (https://github.com/ubunatic/md2pango),
   which is under the MIT license.
*/

#include "libgretl.h"
#include "gretl_mdconv.h"

struct m2g_data_ {
    const char *re;
    const char *sub;
    GRegexCompileFlags cf;
    GRegex *rx;
};

typedef struct m2g_data_ m2g_data;

static m2g_data lr_quotes = {"(\")([^\"]*)(\")", "“\\2”", G_REGEX_UNGREEDY, NULL};

/* m2g_sections define how to detect special markdown sections.
   These expressions scan the full line to detect headings, lists,
   and code.
*/
static m2g_data m2g_sections[] = {
    { "^(#\\s+)(.*)(\\s*)$",  "<@bld=\"\\2\">\n", 0, NULL },
    { "^(##\\s+)(.*)(\\s*)$", "<@itl=\"\\2\">\n", 0, NULL },
    { "^(\\s*[\\*\\-]\\s)(.*)(\\s*)$",  "• \\2", 0, NULL },
    { "^(\\s*[0-9]+\\.\\s)(.*)(\\s*)$", "\\1\\2\n", 0, NULL },
    { "^```[a-z_]*$", NULL, 0, NULL }
};

/* symbols for the m2g_sections above, plus paragraph */
enum {HD1, HD2, UL, OL, CODE, PARA};

static int n_m2g_sections = G_N_ELEMENTS(m2g_sections);

/* m2g_styles handle replacement of inline styles, URIs */
static m2g_data m2g_styles[] = {
    { "(^|[^\\*])(\\*\\*)(.*)(\\*\\*)", "\\1<@bld=\"\\3\">", G_REGEX_UNGREEDY, NULL },
    { "(\\*\\*)(.*)(\\*\\*)([^\\*]|$)", "<@bld=\"\\3\">\\4", G_REGEX_UNGREEDY, NULL },
    { "(^|[^\\*])(\\*)(.*)(\\*)", "\\1<@itl=\"\\3\">", G_REGEX_UNGREEDY, NULL },
    { "(\\*)(.*)(\\*)([^\\*]|$)", "<@itl=\"\\3\">\\4", G_REGEX_UNGREEDY, NULL },
    { "(`)([^`]*)(`)", "<@lit=\"\\2\">", G_REGEX_UNGREEDY, NULL },
#if 0 /* not currently applicable */
    { "(!)?(\\[)(.*)(\\]\\()(.+)(\\))", "<@url=\"\\5\">", 0, NULL },
    { "(!)?(\\[)(.*)(\\]\\(\\))", "<@url=\"\\3\">", 0, NULL },
#endif
    { "(http[s]?:\\/\\/[^\\s',]*)", "<@url=\"\\1\">", 0, NULL },
#if 0 /* messes up '_' in URIs */
    { "(^|[^_])(_)(.*)(_)", "\\1<@itl=\"\\3\">", NULL },
#endif
    { "(^|\\s)(_)(.*)(_)([^_]|$)", "\\1<@itl=\"\\3\">\\5", G_REGEX_UNGREEDY, NULL }
};

static int n_m2g_styles = G_N_ELEMENTS(m2g_styles);

static GError *conv_err;

static int handle_conv_error (void)
{
    if (conv_err != NULL) {
	gretl_errmsg_set(conv_err->message);
	g_error_free(conv_err);
	conv_err = NULL;
	return E_DATA;
    }

    return 0;
}

static gchar **strv_append (gchar **S, const char *s, int trim)
{
    guint n = S == NULL ? 2 : g_strv_length(S) + 2;
    gchar *scpy;

    S = g_realloc(S, n * sizeof *S);
    if (trim) {
	s += strspn(s, " ");
    }
    scpy = S[n-2] = g_strdup(s);
    S[n-1] = NULL;
    if (trim) {
	n = strlen(scpy) - 1;
	if (scpy[n] == '\n') {
	    scpy[n] = '\0';
	}
    }

    return S;
}

void mdconv_cleanup (void)
{
    int i;

    for (i=0; i<n_m2g_sections; i++) {
        if (m2g_sections[i].rx != NULL) {
            g_regex_unref(m2g_sections[i].rx);
            m2g_sections[i].rx = NULL;
        }
    }
    for (i=0; i<n_m2g_styles; i++) {
        if (m2g_styles[i].rx != NULL) {
            g_regex_unref(m2g_styles[i].rx);
            m2g_styles[i].rx = NULL;
        }
    }
}

static gboolean rx_match (const char *s, m2g_data *d)
{
    if (d->rx == NULL) {
	d->rx = g_regex_new(d->re, d->cf, 0, &conv_err);
	if (d->rx == NULL) {
	    return 0;
	}
    }
    return g_regex_match(d->rx, s, 0, NULL);
}

static gchar *rx_replace (const char *s, m2g_data *d)
{
    if (d->rx == NULL) {
	d->rx = g_regex_new(d->re, d->cf, 0, &conv_err);
	if (d->rx == NULL) {
	    return NULL;
	}
    }
    return g_regex_replace(d->rx, s, -1, 0, d->sub, 0, &conv_err);
}

static void output_para (gchar **S, PRN *prn)
{
    guint i, n = g_strv_length(S);

    for (i=0; i<n; i++) {
	if (i > 0) {
	    pputc(prn, ' ');
	}
	pputs(prn, S[i]);
    }
    pputc(prn, '\n');
}

static void output_list_item (gchar **S, PRN *prn)
{
    guint i, n = g_strv_length(S);

    pputs(prn, "<indent>\n");
    for (i=0; i<n; i++) {
	if (i > 0) {
	    pputc(prn, ' ');
	}
	pputs(prn, S[i]);
    }
    pputs(prn, "\n</indent>\n");
}

static void output_code_block (gchar **S, PRN *prn)
{
    guint i, n = g_strv_length(S);

    pputs(prn, "<code>\n");
    for (i=0; i<n; i++) {
	pputs(prn, S[i]);
    }
    pputs(prn, "</code>\n");
}

static void finish_chunk (char ***lines, int mode, PRN *prn)
{
    gchar **S = *lines;

    if (S != NULL) {
        if (mode == UL || mode == OL) {
            output_list_item(S, prn);
        } else if (mode == CODE) {
            output_code_block(S, prn);
        } else {
            output_para(S, prn);
        }
	g_strfreev(S);
	*lines = NULL;
    }
}

/**
 * md_to_gretl:
 * @buf: input buffer containing gretl markdown.
 * @prn: printing struct to which output should be written.
 *
 * The input must conform to gretl's markdown 'flavor', for
 * the details of which see the Gretl Function Package Guide.
 * On success, output is in the form of the internal mark-up
 * used by gretl for display via GtkTextView.
 *
 * Returns: 0 on success or error code on error.
 **/

int md_to_gretl (const char *buf, PRN *prn)
{
    char line[4096];
    gchar **code_lines = NULL;
    gchar **item_lines = NULL;
    gchar **para_lines = NULL;
    int code_start;
    int is_code = 0;
    int list_item = 0;
    int i;

    if (buf == NULL) {
        return 0;
    }

    bufgets_init(buf);

    while (bufgets(line, sizeof line, buf) && !conv_err) {
	gchar *result = line;
	gchar *tmp = NULL;
	int free_result = 0;
	int is_blank;

	is_blank = string_is_blank(line);
#if 0
	printf("@ blank %d, input '%s'\n", is_blank, line);
#endif

        if (is_blank && !is_code) {
            if (list_item) {
                list_item = 0;
                finish_chunk(&item_lines, UL, prn);
            } else if (para_lines != NULL) {
                finish_chunk(&para_lines, PARA, prn);
            }
            gretl_print_ensure_vspace(prn);
            continue;
        }

	if (!is_code && strchr(line, '"')) {
            result = rx_replace(line, &lr_quotes);
	    strcpy(line, result);
	    g_free(result);
	    result = line;
	}

	code_start = 0; /* start of a code block? */

	for (i=0; i<n_m2g_sections && !conv_err; i++) {
	    m2g_data *sect = &m2g_sections[i];

            if (rx_match(result, sect)) {
                if (i == CODE) {
		    /* start/end of code-block */
                    if (!is_code) {
                        /* starts a block */
                        code_start = 1;
                        is_code = 1;
			result = NULL;
                    } else {
                        /* ends a block */
                        is_code = 0;
                        finish_chunk(&code_lines, CODE, prn);
                        result = NULL;
                    }
                } else if (is_code) {
		    /* note: no style replacement */
		    result = line;
		} else if (i == UL || i == OL) {
		    if (list_item) {
			finish_chunk(&item_lines, UL, prn);
		    }
		    list_item = 1;
		    result = rx_replace(line, sect);
		    free_result = 1;
		} else {
		    result = rx_replace(line, sect);
		    free_result = 1;
		}
            }
        }

	if (result == NULL) {
	    continue;
	}

        if (is_code && !code_start && !conv_err) {
            code_lines = strv_append(code_lines, result, 0);
            continue;
        }

        /* all text other than code can be styled */
	for (i=0; i<n_m2g_styles && !conv_err; i++) {
            tmp = rx_replace(result, &m2g_styles[i]);
	    if (free_result) {
		g_free(result);
	    }
	    result = tmp;
	    free_result = 1;
        }

        if (!conv_err) {
            if (list_item) {
                item_lines = strv_append(item_lines, result, 1);
                continue;
            } else {
                para_lines = strv_append(para_lines, result, 1);
                continue;
            }
        }

	if (result != NULL && free_result) {
	    g_free(result);
	}
    }

    /* catch any chunks not flushed by blank line */
    if (para_lines != NULL) {
	finish_chunk(&para_lines, PARA, prn);
    } else if (item_lines != NULL) {
	finish_chunk(&item_lines, UL, prn);
    }

    bufgets_finalize(buf);
    pputc(prn, '\n');

    return handle_conv_error();
}

/* heuristic to determine whether the help text for a gretl
   function package (probably) employs markdown
*/

int help_text_is_markdown (const char *buf)
{
    if (buf == NULL) {
	return 0;
    } else if (strchr(buf, '#') || strchr(buf, '`')) {
	/* these bytes are unlikely to occur in non-markdown text */
	return 1;
    } else {
	/* look for "*word*" */
	const char *s = buf;
	int i;

	while (*s) {
	    if (*s == '*') {
		for (i=1; s[i]; i++) {
		    if (!isalpha(s[i])) {
			break;
		    }
		}
		if (s[i] == '*') {
		    return 1;
		}
		s += i;
	    }
	    s++;
	}
    }

    return 0;
}
