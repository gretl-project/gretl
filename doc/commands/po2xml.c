/* This program can be used to retrieve translations from
   gretl *.po files and write them out in the form wanted for
   XSLT.

   The example below works on the names of "funclist" elements
   in translations of gretl_functions_en.xml, with the aim of
   making them come out correctly in the TeX/PDF version of the
   Gretl Reference. The LANG setting when the program is run
   determines the language of the translated strings.

   Allin Cottrell, August 2021
*/

#include "libgretl.h"
#include <locale.h>

void nls_init (void)
{
    setlocale(LC_ALL, "");
    bindtextdomain(PACKAGE, LOCALEDIR);
    textdomain(PACKAGE);
    bind_textdomain_codeset(PACKAGE, "UTF-8");
}

int main (void)
{
    /* translatable strings from gui/textbuf.c */
    const char *names_en[] = {
	"Accessors",
	"Built-in strings",
	"Functions proper",
	NULL
    };
    int i;

    nls_init();

    /* output lines that can be inserted into hlpstrs_<lang>.xml */
    for (i=0; names_en[i] != NULL; i++) {
	printf(" <gentext key=\"%s\" text=\"%s\" />\n", names_en[i],
	       _(names_en[i]));
    }

    return 0;
}
