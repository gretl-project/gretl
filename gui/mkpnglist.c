/* A little program to generate the XML listing of gretl PNG icons,
   which is required if we're to use glib-compile-resources to compile
   the icons into a GResource for quick loading.

   Allin Cottrell, 2021-07-09. Modified by Dirk Eddelbeuttel 2026-10-28
   to alphabetize the PNG filenames via scandir(), so as to ensure a
   reproducible build.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <dirent.h>

static int cmp_png (const struct dirent **a,
		    const struct dirent **b)
{
    return strcmp((*a)->d_name, (*b)->d_name);
}

#ifdef _WIN32

static int win32_cmp_png (const void *a, const void *p)
{
    const struct dirent *da = *(const struct dirent **) a;
    const struct dirent *db = *(const struct dirent **) b;
    
    return strcmp(da->d_name, db->d_name);
}

int win32_scandir (const char *dirname, struct dirent ***list)
{
    DIR *dirp;
    struct dirent *d;
    struct dirent **mylist = NULL;
    int i, n = 0;

    dirp = opendir(dirname);
    if (dirp == NULL) {
	return -1;
    }

    while ((d = readdir(dirp)) != NULL) {
        struct dirent *copy = malloc(sizeof *copy);
	struct dirent **new_list;

        if (copy == NULL) goto fail;
        memcpy(copy, d, sizeof *copy);
        new_list = realloc(mylist, (n + 1) * sizeof *d);
        if (new_list == NULL) {
            free(copy);
            goto fail;
        }
        mylist = new_list;
        mylist[n++] = copy;
    }
    closedir(dirp);

    if (n > 0) {
        qsort(mylist, n, sizeof *d, win32_cmp_png);
    }

    *list = mylist;
    return n;

 fail:
    for (i=0; i<n; i++) {
	free(mylist[i]);
    }
    free(mylist);
    closedir(dirp);
    return -1;
}

#else

static int cmp_png (const struct dirent **a,
		    const struct dirent **b)
{
    return strcmp((*a)->d_name, (*b)->d_name);
}

#endif

int main (int argc, char **argv)
{
    char dirname[1024];
    DIR *icond[2];
    int px[2] = {16, 24};
    FILE *fp;
    int i;

    if (argc < 2) {
	fprintf(stderr, "%s: need path to share dir\n", argv[0]);
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
	struct dirent **list = NULL;
	int n, j;

	sprintf(dirname, "%s/icons/%dx%d", argv[1], px[i], px[i]);
#ifdef _WIN32
	n = win32_scandir(dirname, &list);
#else
	n = scandir(dirname, &list, NULL, cmp_png);
#endif
	if (n < 0) {
	    fprintf(stderr, "%s: couldn't open %s\n", argv[0], dirname);
	    exit(EXIT_FAILURE);
	}
	for (j=0; j<n; j++) {
	    if (strstr(list[j]->d_name, ".png")) {
		fprintf(fp, "    <file>icons/%dx%d/%s</file>\n",
			px[i], px[i], list[j]->d_name);
	    }
	    free(list[j]);
	}
	free(list);
    }

    fputs("  </gresource>\n", fp);
    fputs("</gresources>\n", fp);

    fclose(fp);

    return 0;
}
