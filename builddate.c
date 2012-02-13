#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

void current_ymd (int *y, int *m, int *d)
{
    time_t t = time(NULL);
    struct tm *lt = localtime(&t);

    *y = lt->tm_year + 1900;
    *m = lt->tm_mon + 1;
    *d = lt->tm_mday;
}

/* See if build.h exists and is up to date, and if not, 
   create/update it. */

int main (void)
{
    int yb = 0, mb = 0, db = 0;
    int y = 0, m = 0, d = 0;
    int n, update = 1;
    char *s, line[128];
    FILE *fp;

    current_ymd(&y, &m, &d);

    fp = fopen("build.h", "r");
    if (fp != NULL) {
	if (fgets(line, sizeof line, fp) != NULL) {
	    s = strstr(line, "20");
	    if (s != NULL) {
		n = sscanf(s, "%d-%d-%d", &yb, &mb, &db);
		if (n == 3 && y == yb && m == mb && d == db) {
		    /* dates agree */
		    update = 0;
		}
	    }
	}
	fclose(fp);
    }

    if (update) {
	fp = fopen("build.h", "w");
	if (fp == NULL) {
	    fprintf(stderr, "Can't write to build.h\n");
	    exit(EXIT_FAILURE);
	} else {
	    printf("Updating build.h\n");
	    fprintf(fp, "#define BUILD_DATE \"%d-%02d-%02d\"\n", 
		    y, m, d);
	    fclose(fp);
	}
    } else {
	printf("build.h is current\n");
    }

    return 0;
}
