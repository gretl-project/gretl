#include <stdio.h>
#include <stdlib.h>

int main (void)
{
    const char *binfile = "bcih.bin";
    const char *datfile = "bcih.dat";
    FILE *fdat, *fbin;
    double xx;
    float x;

    fdat = fopen(datfile, "r");
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
	fwrite(&x, sizeof x, 1, fbin);
    }

    fclose(fbin);
    fclose(fdat);

    return 0;
}
