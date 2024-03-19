/* Auxiliary code to output a gretl language-spec file for gtksourceview,
   or (depending on the option supplied) a block of material for syntax
   highlighting in emacs.
*/

#include <stdio.h>

#include "libgretl.h"
#include "gen_public.h"

const char *special_keyword[] = {
    "for",
    "foreach",
    "funcerr",
    "return",
    "while",
    "elif",
    "eval",
    "3sls",
    "liml",
    "fiml",
    "sur",
    "params",
    "deriv",
    "orthog",
    "weights",
    "hessian",
    "option",
    "options",
    "literal",
    "printf",
    "identity",
    "endog",
    "instr",
    "equations",
    "tsplots",
    NULL
};

const char *gretl_data_types[] = {
    "bool",
    "int",
    "scalar",
    "series",
    "matrix",
    "list",
    "string",
    "bundle",
    "strings",
    "matrices",
    "bundles",
    "lists",
    "arrays",
    "null",
    "void",
    "const",
    "numeric",
    "empty",
    NULL
};

const char *builtin_strings[] = {
    "gretldir",
    "dotdir",
    "workdir",
    "gnuplot",
    "x12a",
    "x12adir",
    "tramo",
    "tramodir",
    "seats",
    "pkgdir",
    "lang",
    NULL
};

/* Similar to strcmp, except that in case two strings match
   to the length of the shorter one, we order the longer
   one first. This is needed to avoid the shorter string
   "masking" the longer one when doing syntax highlighting.
*/

static int mklang_strcmp (const void *a, const void *b)
{
    const char **sa = (const char **) a;
    const char **sb = (const char **) b;
    int na = strlen(*sa);
    int nb = strlen(*sb);
    int n = (na > nb) ? nb : na;
    int cmp = strncmp(*sa, *sb, n);

    if (cmp == 0) {
	cmp = nb - na;
    }

    return cmp;
}

char **make_var_name_list (int *pn)
{
    char **S;
    const char *s;
    int n1, n2, n3, n4;
    int n, m;
    int i;

    n1 = model_var_count();
    n2 = data_var_count();
    n3 = bundle_var_count();
    n4 = gretl_const_count();
    n = n1 + n2 + n3 + n4;

    S = strings_array_new(n);
    if (S == NULL) {
	return NULL;
    }

    m = 0;

    for (i=0; i<n1; i++) {
	s = model_var_name(i);
	if (s == NULL) {
	    continue;
	}
	if (*s == '$') s++;
	S[m++] = gretl_strdup(s);
    }

    for (i=0; i<n2; i++) {
	s = data_var_name(i);
	if (s == NULL) {
	    continue;
	}
	if (*s == '$') s++;
	S[m++] = gretl_strdup(s);
    }

    for (i=0; i<n3; i++) {
	s = bundle_var_name(i);
	if (s == NULL) {
	    continue;
	}
	if (*s == '$') s++;
	S[m++] = gretl_strdup(s);
    }

    for (i=0; i<n4; i++) {
	s = gretl_const_name(i);
	if (s == NULL) {
	    continue;
	}
	if (*s == '$') s++;
	S[m++] = gretl_strdup(s);
    }

    qsort(S, m, sizeof *S, mklang_strcmp);

    *pn = n;

    return S;
}

void output_emacs_block (void)
{
    char **strs;
    int nopts;
    int i, n;

    n = 1;
    fputs("(defvar gretl-command-words\n '(", stdout);
    for (i=1; i<NC; i++) {
	printf("\"%s\"", gretl_command_word(i));
	if (n % 8 == 0) {
	    fputs("\n   ", stdout);
	} else {
	    putchar(' ');
	}
	n++;
    }
    for (i=0; special_keyword[i] != NULL; i++) {
	printf("\"%s\"", special_keyword[i]);
	if (special_keyword[i+1] != NULL) {
	    if (n % 8 == 0) {
		fputs("\n   ", stdout);
	    } else {
		putchar(' ');
	    }
	}
	n++;
    }
    puts(")\n  \"Commands in gretl.\")\n");

    fputs("(defvar gretl-genr-functions\n '(", stdout);
    n = gen_func_count();
    for (i=0; i<n; i++) {
	printf("\"%s\"", gen_func_name(i));
	if (i < n-1) {
	    if ((i+1) % 8 == 0) {
		fputs("\n   ", stdout);
	    } else {
		putchar(' ');
	    }
	}
    }

    puts(")\n  \"Built-in functions in gretl.\")\n");

    strs = get_all_option_strings(&nopts);
    if (strs != NULL) {
	n = 1;
	fputs("(defvar gretl-option-flags\n '(", stdout);
	for (i=1; i<nopts; i++) {
	    printf("\"%s\"", strs[i]);
	    if (i < nopts-1) {
		if (n % 4 == 0) {
		    fputs("\n   ", stdout);
		} else {
		    putchar(' ');
		}
	    }
	    n++;
	}
	puts(")\n  \"Gretl option flags.\")\n");
	strings_array_free(strs, nopts);
    }

    /* "dollar" variables */
    fputs("(defvar gretl-internal-vars\n '(", stdout);

    strs = make_var_name_list(&n);
    if (strs != NULL) {
	for (i=0; i<n; i++) {
	    printf("\"%s\"", strs[i]);
	    if ((i+1) % 7 == 0) {
		fputs("\n   ", stdout);
	    } else {
		putchar(' ');
	    }
	}
	strings_array_free(strs, n);
    }

    puts(")\n  \"Model- and dataset-related variables.\")\n");
}

static int compare_options (const void *a, const void *b)
{
    const char *sa = *(const char **) a;
    const char *sb = *(const char **) b;
    int ret = strcmp(sa, sb);

    return ret == 0 ? ret : -ret;
}

#define DO_FOREIGN 1

#if DO_FOREIGN

static void output_octave_specials (void)
{
    puts("    <context id=\"octave-keyword\" style-ref=\"command\">");
    puts("      <keyword>assert</keyword>");
    puts("      <keyword>break</keyword>");
    puts("      <keyword>case</keyword>");
    puts("      <keyword>catch</keyword>");
    puts("      <keyword>continue</keyword>");
    puts("      <keyword>do</keyword>");
    puts("      <keyword>elseif</keyword>");
    puts("      <keyword>else</keyword>");
    puts("      <keyword>endfor</keyword>");
    puts("      <keyword>endfunction</keyword>");
    puts("      <keyword>endif</keyword>");
    /* here's the required modification */
    puts("      <keyword>end(?! foreign)</keyword>");
    puts("      <keyword>endswitch</keyword>");
    puts("      <keyword>end_try_catch</keyword>");
    puts("      <keyword>end_unwind_protect</keyword>");
    puts("      <keyword>endwhile</keyword>");
    puts("      <keyword>for</keyword>");
    puts("      <keyword>function</keyword>");
    puts("      <keyword>global</keyword>");
    puts("      <keyword>if</keyword>");
    puts("      <keyword>nargin</keyword>");
    puts("      <keyword>nargout</keyword>");
    puts("      <keyword>otherwise</keyword>");
    puts("      <keyword>return</keyword>");
    puts("      <keyword>switch</keyword>");
    puts("      <keyword>try</keyword>");
    puts("      <keyword>until</keyword>");
    puts("      <keyword>unwind_protect_cleanup</keyword>");
    puts("      <keyword>unwind_protect</keyword>");
    puts("      <keyword>while</keyword>");
    puts("    </context>\n");

    puts("    <context id=\"gretl-octave\">");
    puts("      <include>");
    puts("        <context ref=\"octave:line-comment\"/>");
    puts("        <context ref=\"c:string\"/>");
    puts("        <context ref=\"octave:boolean\"/>");
    puts("        <context ref=\"octave:reserved-constant\"/>");
    /* note that here we use our modified keywords list */
    puts("        <context ref=\"octave-keyword\"/>");
    puts("        <context ref=\"def:decimal\"/>");
    puts("        <context ref=\"def:float\"/>");
    puts("        <context ref=\"def:hexadecimal\"/>");
    puts("      </include>");
    puts("    </context>\n");
}

static void output_foreign_context (const char *id, const char *ctxt)
{
    printf("    <context id=\"foreign-%s\" style-inside=\"true\">\n", id);
    printf("      <start>^\\s*(foreign)\\s+language=\\%%{%s}\\%%{foreign-opt}</start>\n", id);
    puts("      <end>^\\s*end\\s+foreign</end>");
    puts("      <include>");
    puts("        <context sub-pattern=\"1\" where=\"start\" style-ref=\"command\"/>");
    puts("        <context sub-pattern=\"2\" where=\"start\" style-ref=\"option\"/>");
    puts("        <context sub-pattern=\"3\" where=\"start\" style-ref=\"option\"/>");
    puts("        <context sub-pattern=\"0\" where=\"end\" style-ref=\"command\"/>");
    printf("        <context ref=\"%s\"/>\n", ctxt);
    puts("      </include>");
    puts("    </context>\n");
}

#endif

void output_lang2_file (int all_foreign)
{
    char **strs;
    int nopts;
    int i, n;

    puts("<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
    puts("<language id=\"gretl\" name=\"gretl\" version=\"2.0\" _section=\"Script\">");
    puts("  <metadata>");
    puts("    <property name=\"mimetypes\">application/x-gretlscript</property>");
    puts("    <property name=\"globs\">*.inp</property>");
    puts("    <property name=\"line-comment-start\">#</property>");
    puts("  </metadata>\n");

    puts("  <styles>");
    puts("    <style id=\"comment\" _name=\"Comment\" map-to=\"def:comment\"/>");
    puts("    <style id=\"function\" _name=\"Function\" map-to=\"def:function\"/>");
    puts("    <style id=\"data-type\" _name=\"Data Type\" map-to=\"def:type\"/>");
    puts("    <style id=\"command\" _name=\"Command\" map-to=\"def:keyword\"/>");
    puts("    <style id=\"option\" _name=\"Option\" map-to=\"def:type\"/>");
    puts("  </styles>\n");

    puts("  <definitions>\n");

#if DO_FOREIGN
    puts("    <define-regex id=\"R\">R|r</define-regex>");
    puts("    <define-regex id=\"python\">Python|python</define-regex>");
    puts("    <define-regex id=\"octave\">Octave|octave</define-regex>");
    if (all_foreign) {
	puts("    <define-regex id=\"stata\">Stata|stata|STATA</define-regex>");
	puts("    <define-regex id=\"ox\">Ox|ox</define-regex>");
    }
    puts("    <define-regex id=\"foreign-opt\">(\\s+--(?:quiet|send-data))?"
	 "(\\s+--(?:send-data|quiet))?</define-regex>\n");

    output_foreign_context("R", "r:r");
    output_foreign_context("python", "python:python");
    output_foreign_context("octave", "gretl-octave");
    if (all_foreign) {
	output_foreign_context("stata", "stata:stata");
	output_foreign_context("ox", "cpp:cpp");
    }
#endif

    puts("    <context id=\"block-comment\" style-ref=\"comment\">");
    puts("      <start>/\\*</start>");
    puts("      <end>\\*/</end>");
    puts("      <include>");
    puts("        <context ref=\"def:escape\"/>");
    puts("        <context ref=\"def:line-continue\"/>");
    puts("      </include>");
    puts("    </context>\n");

    /* gretl data types */
    puts("    <context id=\"gretl-types\" style-ref=\"data-type\">");
    for (i=0; gretl_data_types[i] != NULL; i++) {
	printf("      <keyword>%s</keyword>\n", gretl_data_types[i]);
    }
    puts("    </context>\n");

    /* gretl functions */
    puts("    <context id=\"genr-functions\" style-ref=\"function\">");
    n = gen_func_count();
    for (i=0; i<n; i++) {
	printf("      <keyword>%s</keyword>\n", gen_func_name(i));
    }
    printf("      <keyword>catch</keyword>\n");
    puts("    </context>\n");

    /* gretl commands */
    puts("    <context id=\"commands\" style-ref=\"command\">");
    puts("      <prefix>(^|\\040|\\011)</prefix>");
    puts("      <suffix>(?![\\w\\-\\.\\(])</suffix>");
    for (i=1; i<NC; i++) {
	printf("      <keyword>%s</keyword>\n", gretl_command_word(i));
    }
    /* plus a few specials */
    for (i=0; special_keyword[i] != NULL; i++) {
	printf("      <keyword>%s</keyword>\n", special_keyword[i]);
    }
    puts("    </context>\n");

    /* command option strings */
    strs = get_all_option_strings(&nopts);
    qsort(strs, nopts, sizeof *strs, compare_options);
    if (strs != NULL) {
	puts("    <context id=\"options\" style-ref=\"option\">");
	puts("      <prefix>--</prefix>");
	for (i=1; i<nopts; i++) {
	    printf("      <keyword>%s</keyword>\n", strs[i]);
	}
	puts("    </context>\n");
	strings_array_free(strs, nopts);
    }

    /* dollar variables */
    puts("    <context id=\"internalvars\" style-ref=\"data-type\">");
    puts("      <prefix>\\$</prefix>");
    puts("      <suffix></suffix>");
    strs = make_var_name_list(&n);
    if (strs != NULL) {
	for (i=0; i<n; i++) {
            printf("      <keyword>%s</keyword>\n", strs[i]);
	}
	strings_array_free(strs, n);
    }
    strs = (char **) builtin_strings;
    for (i=0; strs[i] != NULL; i++) {
	printf("      <keyword>%s</keyword>\n", strs[i]);
    }
    puts("    </context>\n");

#if DO_FOREIGN
    /* octave special: prevent octave from eating "end foreign" */
    output_octave_specials();
#endif

    puts("    <context id=\"gretl\">");
    puts("      <include>");
    puts("        <context ref=\"def:shell-like-comment\"/>");
    puts("        <context ref=\"def:string\"/>");
    puts("        <context ref=\"block-comment\"/>");
#if DO_FOREIGN
    puts("        <context ref=\"foreign-R\"/>");
    puts("        <context ref=\"foreign-python\"/>");
    puts("        <context ref=\"foreign-octave\"/>");
    if (all_foreign) {
	puts("        <context ref=\"foreign-stata\"/>");
	puts("        <context ref=\"foreign-ox\"/>");
    }
#endif
    puts("        <context ref=\"gretl-types\"/>");
    puts("        <context ref=\"commands\"/>");
    puts("        <context ref=\"genr-functions\"/>");
    puts("        <context ref=\"options\"/>");
    puts("        <context ref=\"internalvars\"/>");
    puts("      </include>");
    puts("    </context>\n");
    puts("  </definitions>");

    puts("</language>");
}

int main (int argc, char **argv)
{
    char *opt = NULL;

    if (argc == 2) {
	opt = argv[1];
    }

    if (opt != NULL) {
	if (!strcmp(opt, "--emacs")) {
	    output_emacs_block();
	} else if (!strcmp(opt, "--gtksv")) {
	    output_lang2_file(0);
	}
    } else {
	output_lang2_file(1);
    }

    return 0;
}
