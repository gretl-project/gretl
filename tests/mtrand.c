#include <stdio.h>
#include <stdlib.h>

#include "libgretl.h"

#define NINTS 2867200

int main (int argc, char *argv[])
{
    int i, seed = 0;
    unsigned r;
    FILE *fp;

    if (argc > 1) {
	seed = atoi(argv[1]);
    } 

    gretl_rand_init();

    if (seed) {
	gretl_rand_set_seed((unsigned) seed);
    }

    fp = fopen("gretl_output.32", "wb");

    for (i=0; i<NINTS; i++) {
	r = gretl_rand_int();
	fwrite(&r, sizeof r, 1, fp);
    }
    fclose(fp);

    gretl_rand_free();

    return 0;
}
