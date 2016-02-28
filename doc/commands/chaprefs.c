#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* 
   Example aux file line:

   \newlabel{chap:timeseries}{{25}{216}{Univariate time series models}{chapter.25}{}}

   We want to extract the chapter title (or maybe the chapter number?)
   given the label (e.g. "chap:timeseries").
*/

void get_auxfile_info (const char *src)
{
    char targ[64];
    char tmp[512];
    char title[80];
    FILE *fp;
    char *s;
    int ch, pg, ch2;
    int n;

    /* try opening @src.aux */
    sprintf(tmp, "%s.aux", src);
    fp = fopen(tmp, "r");
    if (fp == NULL) {
	return;
    }

    /* the string we're looking for */
    sprintf(targ, "newlabel{chap:%s}", src);

    while (fgets(tmp, sizeof tmp, fp)) {
	if ((s = strstr(tmp, targ)) != NULL) {
	    s = strchr(s, '}') + 1;
	    n = sscanf(s, "{{%d}{%d}{%79[^}]}{chapter.%d}", &ch, &pg, title, &ch2);
	    if (n == 4) {
		printf(" <chapref key=\"chap:%s\" title=\"%s\" />\n", src, title);
	    }
	    break;
	}
    }

    fclose(fp);
}

int main (void)
{
    FILE *fp;
    char line[512];
    char chap[32];
    char *s;

    fp = fopen("gretl-guide.tex", "r");

    if (fp == NULL) {
	exit(EXIT_FAILURE);
    }

    puts("<?xml version=\"1.0\" encoding=\"utf-8\" ?>\n");
    puts("<localization language=\"en\">");

    while (fgets(line, sizeof line, fp)) {
	s = line;
	s += strspn(s, " ");
	if (!strncmp(s, "\\include{", 9)) {
	    /* get info for an included chapter file */
	    sscanf(s + 9, "%31[^}]", chap);
	    get_auxfile_info(chap);
	}
    }

    puts("</localization>");

    fclose(fp);

    return 0;
}
