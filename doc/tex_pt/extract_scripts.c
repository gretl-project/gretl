/* extract code fragments from the gretl manual .tex files
   so they can be checked for up-to-dateness.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <dirent.h>
#include <sys/stat.h>

char line[1024];
const char *START = "begin{code}";
const char *END = "end{code}";

int first_char (const char *s)
{
    return *(s + strspn(s, " \t"));
}

int process_tex_file (const char *fname)
{
    FILE *fp;

    fp = fopen(fname, "r");
    if (fp == NULL) {
	fprintf(stderr, "couldn't open %s\n", fname);
	return 1;
    }

    while (fgets(line, sizeof line, fp)) {
	if (strstr(line, START)) {
	    printf("\n\n# %s", line);
	    while (fgets(line, sizeof line, fp)) {
		if (strstr(line, END)) {
		    printf("# %s", line);
		    break;
		} else if (first_char(line) == '\\') {
		    printf("# %s", line);
		} else {
		    fputs(line, stdout);
		}
	    }
	}
    }

    fclose(fp);
	
    return 0;
}

int main (int argc, char **argv)
{
    const char *dirname = ".";
    const char *fname;
    DIR *dir = NULL;
    struct dirent *dirent;
    
    if ((dir = opendir(dirname)) == NULL) {
        printf("%s: can't open directory\n", argv[0]);
	exit(EXIT_FAILURE);
    }

    while ((dirent = readdir(dir)) != NULL) {
	fname = dirent->d_name;
	if (strstr(fname, ".tex")) {
	    fprintf(stderr, "processing %s\n", fname);
	    process_tex_file(fname);
	}
    }

    closedir(dir);

    return 0;
}
