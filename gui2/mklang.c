#include <stdio.h>

#include "libgretl.h"
#include "genstack.h"

void output_lang_file (void)
{
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

    puts("<string _name = \"String\" style = \"String\" end-at-line-end = \"TRUE\">");
    puts(" <start-regex>&quot;</start-regex>");
    puts(" <end-regex>&quot;</end-regex>");
    puts("</string>\n");

    /* gretl commands */
    puts("<keyword-list _name = \"Commands\" style = \"Keyword\" case-sensitive=\"TRUE\"");
    puts("  beginning-regex=\"\\b\" end-regex=\"\\b\">");
    for (i=1; i<NC; i++) {
	printf(" <keyword>%s</keyword>\n", gretl_command_word(i));
    }
    puts("</keyword-list>\n");

    /* functions in "genr" command */
    puts("<keyword-list _name = \"Genr-functions\" style = \"Function\"");
    puts("  case-sensitive=\"TRUE\">");
    for (i=1; i<T_IDENTITY; i++) {
	printf(" <keyword>%s</keyword>\n", get_genr_func_word(i));
    }    
    puts("</keyword-list>\n");

    puts("<pattern-item _name = \"Internal-Variables\" style = \"Data Type\">");
    puts(" <regex>[$][$]?[a-zA-Z_][a-zA-Z0-9_]*</regex>");
    puts("</pattern-item>\n");

    puts("<pattern-item _name = \"Options\" style = \"Data Type\">");
    puts(" <regex>--[a-z-]*[ \\n]*</regex>");
    puts("</pattern-item>\n");

    puts("<string _name = \"Character Constant\" style = \"String\" end-at-line-end = \"TRUE\">");
    puts(" <start-regex>&apos;</start-regex>");
    puts(" <end-regex>&apos;</end-regex>");
    puts("</string>\n");

    puts("</language>");
}

int main (void)
{
    output_lang_file();

    return 0;
}
