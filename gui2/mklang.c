#include <stdio.h>

#include "libgretl.h"
#include "genstack.h"

void output_emacs_block (void)
{
    int i, n;

    /* gretl commands */
    n = 1;
    fputs("(defvar gretl-command-words\n '(", stdout);
    for (i=1; i<NC; i++) {
	if (strcmp(gretl_command_word(i), "matrix")) {
	    printf("\"%s\"", gretl_command_word(i));
	}
	if (i < NC-1) {
	    if (n % 8 == 0) {
		fputs("\n   ", stdout);
	    } else {
		putchar(' ');
	    }
	}
	n++;
    } 
    puts(")\n  \"Commands in Gretl (these names are also reserved).\")");

    /* functions in "genr" command */
    n = 1;
    fputs("(defvar gretl-genr-functions\n '(", stdout);
    for (i=1; i<T_IDENTITY; i++) {
	printf("\"%s\"", get_genr_func_word(i));
	if (i < T_IDENTITY-1) {
	    if (n % 8 == 0) {
		fputs("\n   ", stdout);
	    } else {
		putchar(' ');
	    }
	}
	n++;
    }
    puts(")\n  \"Builtin functions for Gretl's genr command.\")\n");
}

void output_lang_file (void)
{
    char **strs;
    int nopts;
    int i;

    puts("<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
    puts("<!DOCTYPE language SYSTEM \"language.dtd\">");
    puts("<language _name=\"gretl\" version=\"1.0\" _section=\"Sources\" "
	 "mimetypes=\"application/x-gretlsession\">\n");

    puts("<escape-char>\\</escape-char>");

    puts("<line-comment _name = \"Line Comment\" style= \"Comment\">");   
    puts(" <start-regex>#</start-regex>");   
    puts("</line-comment>\n");   

    puts("<block-comment _name = \"Block Comment\" style = \"Comment\">");
    puts(" <start-regex>\\(\\*</start-regex>");
    puts(" <end-regex>\\*\\)</end-regex>");
    puts("</block-comment>\n");

    puts("<block-comment _name = \"C-style Block Comment\" style = \"Comment\">");
    puts(" <start-regex>/\\*</start-regex>");
    puts(" <end-regex>\\*/</end-regex>");
    puts("</block-comment>\n");

    puts("<string _name = \"String\" style = \"String\" end-at-line-end = \"TRUE\">");
    puts(" <start-regex>&quot;</start-regex>");
    puts(" <end-regex>&quot;</end-regex>");
    puts("</string>\n");

    /* gretl data types */
    puts("<keyword-list _name = \"Grel-types\" style = \"Data Type\"");
    puts("  case-sensitive=\"TRUE\">");
    puts(" <keyword>scalar</keyword>");
    puts(" <keyword>series</keyword>");
    puts(" <keyword>matrix</keyword>");
    puts(" <keyword>list</keyword>");
    puts("</keyword-list>\n");

    /* gretl commands */
    puts("<keyword-list _name = \"Commands\" style = \"Keyword\" case-sensitive=\"TRUE\"");
    puts("  beginning-regex=\"\\b\" end-regex=\"\\b\">");
    for (i=1; i<NC; i++) {
	if (strcmp(gretl_command_word(i), "matrix")) {
	    printf(" <keyword>%s</keyword>\n", gretl_command_word(i));
	}
    }
    puts("</keyword-list>\n");

    /* functions in "genr" command */
    puts("<keyword-list _name = \"Genr-functions\" style = \"Function\"");
    puts("  case-sensitive=\"TRUE\">");
    for (i=1; i<T_IDENTITY; i++) {
	printf(" <keyword>%s</keyword>\n", get_genr_func_word(i));
    }    
    puts("</keyword-list>\n");

    /* command option strings */
    strs = get_all_option_strings(&nopts);
    if (strs != NULL) {
	puts("<keyword-list _name = \"Options\" style = \"Data Type\" "
	     "case-sensitive=\"TRUE\">");
	for (i=1; i<nopts; i++) {
	    printf(" <keyword>%s</keyword>\n", strs[i]);
	}    
	puts("</keyword-list>\n");
	free_strings_array(strs, nopts);
    }

    puts("<pattern-item _name = \"Internal-Variables\" style = \"Data Type\">");
    puts(" <regex>[$][$]?[a-zA-Z_][a-zA-Z0-9_]*</regex>");
    puts("</pattern-item>\n");

    puts("<string _name = \"Character Constant\" style = \"String\" end-at-line-end = \"TRUE\">");
    puts(" <start-regex>&apos;</start-regex>");
    puts(" <end-regex>&apos;</end-regex>");
    puts("</string>\n");

    puts("</language>");
}

int main (int argc, char **argv)
{
    if (argc == 2 && !strcmp(argv[1], "--emacs")) {
	output_emacs_block();
    } else {
	output_lang_file();
    }

    return 0;
}
