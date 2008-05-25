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
    "const",
    "3sls",
    "liml",
    "fiml",
    "sur",
    NULL
};

const char *gretl_data_types[] = {
    "bool",
    "int",
    "scalar",
    "series",
    "matrix",
    "list",
    "null",
    "void",
    NULL
};

void output_emacs_block (void)
{
    const char *s;
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
	free_strings_array(strs, nopts);
    }

    /* internal "dollar" variables */
    fputs("(defvar gretl-internal-vars\n '(", stdout);
    /* model variables */
    n = model_var_count();
    for (i=0; i<n; i++) {
	s = model_var_name(i);
	if (s == NULL) {
	    continue;
	}
	if (*s == '$') s++;
	printf("\"%s\"", s + 1);
	if ((i+1) % 8 == 0) {
	    fputs("\n   ", stdout);
	} else {
	    putchar(' ');
	}
    }
    /* dataset variables */
    n = data_var_count();
    for (i=0; i<n; i++) {
	s = data_var_name(i);
	if (s == NULL) {
	    continue;
	}
	if (*s == '$') s++;
	printf("\"%s\"", s);
	if ((i+1) % 8 == 0) {
	    fputs("\n   ", stdout);
	} else {
	    putchar(' ');
	}
    }
    puts(")\n  \"Model- and dataset-related variables.\")\n");    
}

void output_lang2_file (void)
{
    const char *s;
    char **strs;
    int nopts;
    int i, n;

    puts("<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
    puts("<language id=\"gretl\" _name=\"gretl\" version=\"2.0\" _section=\"Sources\">");
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

    /* gretl data types */
    puts(" <context id=\"gretl-types\" style-ref=\"data-type\">");
    for (i=0; gretl_data_types[i] != NULL; i++) {
	printf("  <keyword>%s</keyword>\n", gretl_data_types[i]);  
    }
    puts(" </context>");

    /* gretl commands (N.B. old case-sensitive="TRUE"?) */
    puts(" <context id=\"commands\" style-ref=\"keyword\">");
    for (i=1; i<NC; i++) {
	if (strcmp(gretl_command_word(i), "matrix")) {
	    printf("  <keyword>%s</keyword>\n", gretl_command_word(i));
	}
    }
    /* plus a few specials */
    for (i=0; special_keyword[i] != NULL; i++) {
	printf("  <keyword>%s</keyword>\n", special_keyword[i]);
    }
    puts(" </context>");

    /* functions in "genr" command */
    puts(" <context id=\"genr-functions\" style-ref=\"function\">");
    n = gen_func_count();
    for (i=0; i<n; i++) {
	printf("  <keyword>%s</keyword>\n", gen_func_name(i));
    }    
    puts(" </context>");

    /* command option strings */
    strs = get_all_option_strings(&nopts);
    if (strs != NULL) {
	puts(" <context id=\"options\" style-ref=\"data-type\">");
	puts(" <prefix>--</prefix>");
	for (i=1; i<nopts; i++) {
	    printf("  <keyword>%s</keyword>\n", strs[i]);
	}    
	puts(" </context>");
	free_strings_array(strs, nopts);
    }

    /* dollar variables */
    puts(" <context id=\"internalvars\" style-ref=\"data-type\">");
    puts("  <prefix>\\$</prefix>");
    puts("  <suffix></suffix>");
    n = model_var_count();
    for (i=0; i<n; i++) {
	s = model_var_name(i);
	if (s == NULL) {
	    continue;
	}
	if (*s == '$') s++;
	printf("  <keyword>%s</keyword>\n", s);
    }
    n = data_var_count();
    for (i=0; i<n; i++) {
	s = data_var_name(i);
	if (s == NULL) {
	    continue;
	}
	if (*s == '$') s++;	
	printf("  <keyword>%s</keyword>\n", s);
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

void output_lang1_file (void)
{
    const char *s;
    char **strs;
    int nopts;
    int i, n;

    puts("<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
    puts("<!DOCTYPE language SYSTEM \"language.dtd\">");
    puts("<language _name=\"gretl\" version=\"1.0\" _section=\"Sources\" "
	 "mimetypes=\"application/x-gretlscript\">\n");

    puts("<escape-char>\\</escape-char>\n");

    puts("<line-comment _name = \"Line Comment\" style= \"Comment\">");   
    puts(" <start-regex>#</start-regex>");   
    puts("</line-comment>\n");   

    puts("<block-comment _name = \"Block Comment\" style = \"Comment\">");
    puts(" <start-regex>/\\*</start-regex>");
    puts(" <end-regex>\\*/</end-regex>");
    puts("</block-comment>\n");

    puts("<string _name = \"String\" style = \"String\" end-at-line-end = \"TRUE\">");
    puts(" <start-regex>&quot;</start-regex>");
    puts(" <end-regex>&quot;</end-regex>");
    puts("</string>\n");

    /* gretl data types */
    puts("<keyword-list _name = \"Gretl-types\" style = \"Data Type\" case-sensitive=\"TRUE\">");
    for (i=0; gretl_data_types[i] != NULL; i++) {
	printf(" <keyword>%s</keyword>\n", gretl_data_types[i]);  
    }
    puts("</keyword-list>\n");

    /* gretl commands */
    puts("<keyword-list _name = \"Commands\" style = \"Keyword\" case-sensitive=\"TRUE\">");
    for (i=1; i<NC; i++) {
	printf(" <keyword>%s</keyword>\n", gretl_command_word(i));
    }
    /* plus a few specials */
    for (i=0; special_keyword[i] != NULL; i++) {
	printf(" <keyword>%s</keyword>\n", special_keyword[i]);
    }
    puts("</keyword-list>\n");

    /* functions in "genr" command */
    puts("<keyword-list _name = \"Genr-functions\" style = \"Function\" case-sensitive=\"TRUE\">");
    n = gen_func_count();
    for (i=0; i<n; i++) {
	printf(" <keyword>%s</keyword>\n", gen_func_name(i));
    }    
    puts("</keyword-list>\n");

    /* command option strings */
    strs = get_all_option_strings(&nopts);
    if (strs != NULL) {
	puts("<keyword-list _name = \"Options\" style = \"Data Type\" case-sensitive=\"TRUE\"");
	puts(" match-empty-string-at-beginning = \"FALSE\" beginning-regex=\"--\">");
	for (i=1; i<nopts; i++) {
	    printf(" <keyword>%s</keyword>\n", strs[i]);
	}    
	puts("</keyword-list>\n");
	free_strings_array(strs, nopts);
    }

    /* dollar variables */
    puts("<keyword-list _name = \"InternalVars\" style = \"Data Type\" case-sensitive=\"TRUE\"");
    puts(" match-empty-string-at-beginning = \"FALSE\" match-empty-string-at-end = \"FALSE\"");
    puts(" beginning-regex=\"\\$\">");

    n = model_var_count();
    for (i=0; i<n; i++) {
	s = model_var_name(i);
	if (s == NULL) {
	    continue;
	}
	if (*s == '$') s++;
	printf(" <keyword>%s</keyword>\n", s);
    }

    n = data_var_count();
    for (i=0; i<n; i++) {
	s = data_var_name(i);
	if (s == NULL) {
	    continue;
	}
	if (*s == '$') s++;	
	printf(" <keyword>%s</keyword>\n", s);
    }
	
    puts("</keyword-list>\n");

#if 0
    puts("<string _name = \"Character Constant\" style = \"String\" end-at-line-end = \"TRUE\">");
    puts(" <start-regex>&apos;</start-regex>");
    puts(" <end-regex>&apos;</end-regex>");
    puts("</string>\n");
#endif

    puts("</language>");
}

int main (int argc, char **argv)
{
    if (argc == 2 && !strcmp(argv[1], "--emacs")) {
	output_emacs_block();
    } else if (argc == 2 && !strcmp(argv[1], "--lang2")) {
	output_lang2_file();
    } else {
	output_lang1_file();
    }

    return 0;
}
