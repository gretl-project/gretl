#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define VERBOSE 0

void write_chapter_info (const char *path, const char *src,
			 int *chapnum)
{
    const char *targ = "\\label{chap:";
    char label[32];
    char tmp[512];
    FILE *fp;
    char *s;

    /* try opening @src.tex */
    sprintf(tmp, "%s%s.tex", path, src);
    fp = fopen(tmp, "r");
    if (fp == NULL) {
	fprintf(stderr, "no tex file for %s (?)\n", src);
	return;
    }

    while (fgets(tmp, sizeof tmp, fp)) {
	if ((s = strstr(tmp, targ)) != NULL) {
	    sscanf(s + 12, "%31[^}]", label);
	    printf(" <ref key=\"chap:%s\" chapter=\"%d\"/>\n", label, *chapnum);
	    *chapnum += 1;
	    break;
	}
    }

    fclose(fp);
}

int main (int argc, char **argv)
{
    FILE *fp;
    char path[FILENAME_MAX];
    char line[512];
    char chap[32];
    char *guide;
    char *s;
    int chapnum = 1;

    if (argc < 2) {
	fprintf(stderr, "Give the path to gretl-guide.tex\n");
	exit(EXIT_FAILURE);
    }

    guide = argv[1];

    *path = '\0';
    s = strrchr(guide, '/');
    if (s != NULL) {
	strncat(path, guide, s - guide + 1);
    }

    fp = fopen(guide, "r");
    if (fp == NULL) {
	fprintf(stderr, "Couldn't open %s\n", guide);
	exit(EXIT_FAILURE);
    }    

    puts("<?xml version=\"1.0\" encoding=\"utf-8\" ?>");
    puts("<refsets>");
    puts("<refset id=\"guide-chapters\">");

    /* We search in gretl-guide.tex for \include{} lines
       that give the names of chapter files: having found
       one, we look for the corresponding .tex file and
       check that it has a chap:foo label; if so, we write
       the chapter number (could also grab the title)
       and write the info into an XML file for use with
       the XSL transformation of the "online" help files.
       
       In this way we're able to give a chapter number
       instead of just saying "see the Gretl User's Guide".
       This does depend on having well-regimented labels
       in the TeX chapter sources.
    */

    while (fgets(line, sizeof line, fp)) {
	s = line;
	s += strspn(s, " ");
	if (!strncmp(s, "\\guidechap{", 11)) {
	    /* get info for an included chapter file */
	    sscanf(s + 11, "%31[^}]", chap);
#if VERBOSE
	    fprintf(stderr, "got include of '%s'\n", chap);
#endif
	    write_chapter_info(path, chap, &chapnum);
	}
    }

    puts("</refset>");
    puts("</refsets>");

    fclose(fp);

    return 0;
}
