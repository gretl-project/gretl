#include <stdio.h>

#include "libgretl.h"
#include "genparse.h"

const char *special_keyword[] = {
    "for",
    "foreach",
    "funcerr",
    "return", 
    "while",
    "elif",
    "eval",
    "const",
    "3sls",
    "liml",
    "fiml",
    "sur",
    "params",    
    "deriv",
    "orthog",
    "weights",
    "hessian",
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
    "null",
    "void",
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
    int n1, n2, n, m;
    int i;

    n1 = model_var_count();
    n2 = data_var_count();
    n = n1 + n2;

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

    qsort(S, m, sizeof *S, mklang_strcmp);

    *pn = n;

    return S;
}

void output_emacs_block (void)
{
    char **strs;
    int nopts;
    int i, n;

    /* gretl commands */
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
    puts(")\n  \"Commands in Gretl (these names are also reserved).\")\n");

    /* functions in "genr" command */
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

    puts(")\n  \"Builtin functions for Gretl's genr command.\")\n");

    /* option strings */
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

    /* internal "dollar" variables */
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

void output_lang2_file (void)
{
    char **strs;
    int nopts;
    int i, n;

    puts("<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
    puts("<language id=\"gretl\" _name=\"gretl\" version=\"2.0\" _section=\"Scripts\">");
    puts("<metadata>");
    puts(" <property name=\"mimetypes\">application/x-gretlscript</property>");
    puts(" <property name=\"globs\">*.inp</property>");
    puts("</metadata>");

    puts("<styles>");
    puts(" <style id=\"comment\" _name=\"Comment\" map-to=\"def:comment\"/>");
    puts(" <style id=\"function\" _name=\"Function\" map-to=\"def:function\"/>");
    puts(" <style id=\"data-type\" _name=\"Data Type\" map-to=\"def:type\"/>");
    puts(" <style id=\"string\" _name=\"String\" map-to=\"def:string\"/>");
    puts(" <style id=\"keyword\" _name=\"Keyword\" map-to=\"def:keyword\"/>");
    puts("</styles>");

    puts("<definitions>");
    puts(" <context id=\"line-comment\" style-ref=\"comment\" end-at-line-end=\"true\">");
    puts("  <start>#</start>");
    puts("  <include>");
    puts("   <context ref=\"def:escape\"/>");
    puts("   <context ref=\"def:line-continue\"/>");
    puts("  </include>");
    puts(" </context>");
    puts(" <context id=\"block-comment\" style-ref=\"comment\">");
    puts("  <start>/\\*</start>");
    puts("  <end>\\*/</end>");
    puts("  <include>");
    puts("   <context ref=\"def:escape\"/>");
    puts("   <context ref=\"def:line-continue\"/>");
    puts(" </include>");
    puts(" </context>");
    puts(" <context id=\"string\" style-ref=\"string\" end-at-line-end=\"true\">");
    puts("  <start>\"</start>");
    puts("  <end>\"</end>");
    puts("  <include>");
    puts("   <context ref=\"def:escape\"/>");
    puts("   <context ref=\"def:line-continue\"/>");
    puts("  </include>");
    puts(" </context>");
#if 0 /* not yet */
    puts(" <context id=\"foreign\" style-inside=\"true\" style-ref=\"comment\">");
    puts("  <start>(?&lt;=foreign language)</start>");
    puts("  <end>(?=end foreign)</end>");
    puts(" </context>");
#endif

    /* gretl data types */
    puts(" <context id=\"gretl-types\" style-ref=\"data-type\">");
    for (i=0; gretl_data_types[i] != NULL; i++) {
	printf("  <keyword>%s</keyword>\n", gretl_data_types[i]);  
    }
    puts(" </context>");

    /* gretl functions */
    puts(" <context id=\"genr-functions\" style-ref=\"function\">");
    n = gen_func_count();
    for (i=0; i<n; i++) {
	printf("  <keyword>%s</keyword>\n", gen_func_name(i));
    }
    printf("  <keyword>catch</keyword>\n");
    puts(" </context>");

    /* gretl commands */
    puts(" <context id=\"commands\" style-ref=\"keyword\">");
    puts("  <prefix>(^|\\040|\\011)</prefix>");
    puts("  <suffix>(?![\\w\\-\\.\\(])</suffix>");
    for (i=1; i<NC; i++) {
	printf("  <keyword>%s</keyword>\n", gretl_command_word(i));
    }
    /* plus a few specials */
    for (i=0; special_keyword[i] != NULL; i++) {
	printf("  <keyword>%s</keyword>\n", special_keyword[i]);
    }
    puts(" </context>");

    /* command option strings */
    strs = get_all_option_strings(&nopts);
    qsort(strs, nopts, sizeof *strs, compare_options);
    if (strs != NULL) {
	puts(" <context id=\"options\" style-ref=\"data-type\">");
	puts(" <prefix>--</prefix>");
	for (i=1; i<nopts; i++) {
	    printf("  <keyword>%s</keyword>\n", strs[i]);
	}    
	puts(" </context>");
	strings_array_free(strs, nopts);
    }

    /* dollar variables */
    puts(" <context id=\"internalvars\" style-ref=\"data-type\">");
    puts("  <prefix>\\$</prefix>");
    puts("  <suffix></suffix>");

    strs = make_var_name_list(&n);
    if (strs != NULL) {
	for (i=0; i<n; i++) {
	   printf("  <keyword>%s</keyword>\n", strs[i]); 
	}
	strings_array_free(strs, n);
    }

    puts(" </context>");	
    
    puts(" <context id=\"gretl\">");
    puts("  <include>");
    puts("   <context ref=\"line-comment\"/>");
    puts("   <context ref=\"block-comment\"/>");
    puts("   <context ref=\"string\"/>");
    puts("   <context ref=\"gretl-types\"/>");
    puts("   <context ref=\"commands\"/>");
    puts("   <context ref=\"genr-functions\"/>");
    puts("   <context ref=\"options\"/>");
    puts("   <context ref=\"internalvars\"/>");
    puts("  </include>");
    puts(" </context>");
    puts("</definitions>");

    puts("</language>");
}

int main (int argc, char **argv)
{
    if (argc == 2 && !strcmp(argv[1], "--emacs")) {
	output_emacs_block();
    } else {
	output_lang2_file();
    }

    return 0;
}
