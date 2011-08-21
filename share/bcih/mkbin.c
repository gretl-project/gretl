#include <stdio.h>
#include <stdlib.h>
#include <string.h>

float reverse_float (float x)
{
    float y;
    char *src = (char *) &x;
    char *targ = (char *) &y;

    /* swap the bytes */
    targ[0] = src[3];
    targ[1] = src[2];
    targ[2] = src[1];
    targ[3] = src[0];

    return y;
}

int main (int argc, char **argv)
{
    const char *binfile = "fedstl.bin";
    const char *datfile = "fedstl.dat";
    char datpath[512];
    int swap_ends = 0;
    FILE *fdat, *fbin;
    double xx;
    float x;

    if (argc == 2 && !strcmp(argv[1], "--swap-ends")) {
	swap_ends = 1;
	argc--;
    } else if (argc == 3 && !strcmp(argv[2], "--swap-ends")) {
	swap_ends = 1;
	argc--;
    }

    if (swap_ends) {
	puts("*** making data file with swapped endianness");
    }

    if (argc > 1) {
	sprintf(datpath, "%s/%s", argv[1], datfile);
    } else {
	strcpy(datpath, datfile);
    }

    fdat = fopen(datpath, "r");
    if (fdat == NULL) {
	fprintf(stderr, "Couldn't open %s\n", datfile);
	exit(EXIT_FAILURE);
    }

    fbin = fopen(binfile, "wb");
    if (fbin == NULL) {
	fprintf(stderr, "Couldn't open %s\n", binfile);
	exit(EXIT_FAILURE);
    }

    while (fscanf(fdat, "%lf", &xx) == 1) {
	x = xx;
	if (swap_ends) {
	    reverse_float(x);
	}
	fwrite(&x, sizeof x, 1, fbin);
    }

    fclose(fbin);
    fclose(fdat);

    return 0;
}
