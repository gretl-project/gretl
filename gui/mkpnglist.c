/* Little program to generate the XML listing of gretl PNG
   icons, which is required if we're to use glib-compile-resources
   to compile the icons into a GResource for quick loading.

   Allin Cottrell, 2021-07-09
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <dirent.h>

int main (int argc, char **argv)
{
    struct dirent *de;
    char dirname[512];
    DIR *icond[2];
    int px[2] = {16, 24};
    FILE *fp;
    int i;

    if (argc < 2) {
	fprintf(stderr, "%s: need path to share dir\n", argv[0]);
	exit(EXIT_FAILURE);
    }

    sprintf(dirname, "%s/icons/16x16", argv[1]);
    icond[0] = opendir(dirname);
    if (icond[0] == NULL) {
	fprintf(stderr, "%s: couldn't open %s\n", argv[0], dirname);
	exit(EXIT_FAILURE);
    }

    sprintf(dirname, "%s/icons/24x24", argv[1]);
    icond[1] = opendir(dirname);
    if (icond[1] == NULL) {
	fprintf(stderr, "%s: couldn't open %s\n", argv[0], dirname);
	exit(EXIT_FAILURE);
    }    

    fp = fopen("gretl-icons.xml", "w");
    if (fp == NULL) {
	fprintf(stderr, "%s: couldn't open %s\n", argv[0], "gretl-icons.xml");
	fclose(fp);
	exit(EXIT_FAILURE);
    }

    fputs("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n", fp);
    fputs("<gresources>\n", fp);
    fputs("  <gresource prefix=\"/gretl\">\n", fp);

    for (i=0; i<2; i++) {
	while ((de = readdir(icond[i])) != NULL) {
	    if (strstr(de->d_name, ".png")) {
		fprintf(fp, "    <file>icons/%dx%d/%s</file>\n",
			px[i], px[i], de->d_name);
	    }
	}
    }

    fputs("  </gresource>\n", fp);
    fputs("</gresources>\n", fp);

    closedir(icond[0]);
    closedir(icond[1]);
    fclose(fp);

    return 0;
}
