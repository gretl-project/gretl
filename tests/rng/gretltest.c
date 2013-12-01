#include "unif01.h"
#include "bbattery.h"

#include <gretl/libgretl.h>

static double u01_from_normal (void)
{
    double x = gretl_one_snormal();

    return normal_cdf(x);
}

static double u01_direct (void)
{
    return gretl_rand_01();
}

static void show_help (const char *prog)
{
    printf("options for %s:\n\n", prog);

    fputs("--size={0,1,2} : for SmallCrush, Crush, BigCrush\n", stdout);
    fputs("                 (default is Small)\n\n", stdout);
    
    fputs("--normal : base test on gretl's normal RNG\n", stdout);
    fputs("           (default, uniform RNG)\n\n", stdout);

    fputs("--box-muller : when using normal RNG, use Box-Muller\n", stdout);
    fputs("               instead of Ziggurat (Ziggurat is default)\n\n", stdout);

    exit(EXIT_SUCCESS);
}

enum {
    NORMAL_RNG = 1 << 0,
    BOX_MULLER = 1 << 1
};

static int read_opts (int argc, char **argv, int *flags,
		      int *testsize)
{
    char *s;
    int i, sz = -1;
    int err = 0;

    for (i=1; i<argc && !err; i++) {
	s = argv[i];
	if (!strcmp(s, "--help") || !strcmp(s, "-h")) {
	    show_help(argv[0]);
	} else if (!strncmp(s, "--size=", 7)) {
	    sz = atoi(s + 7);
	    if (sz < 0 || sz > 2) {
		fprintf(stderr, "%s: invalid test size %d\n", argv[0], sz);
		err = 1;
	    }
	} else if (!strcmp(s, "--uniform")) {
	    ; /* the default, np-op */
	} else if (!strcmp(s, "--normal")) {
	    *flags |= NORMAL_RNG;
	} else if (!strcmp(s, "--box-muller")) {
	    *flags |= BOX_MULLER;
	} else {
	    fprintf(stderr, "%s: invalid option '%s'\n", argv[0], s);
	}
    }

    if (!err) {
	const char *size_string[] = {
	    "SmallCrush",
	    "Crush",
	    "BigCrush"
	};
	    
	if (sz > 0) {
	    *testsize = sz;
	}
	printf("*** %s: running %s using %s generator\n", argv[0],
	       size_string[*testsize], (*flags & NORMAL_RNG) ?
	       "normal" : "uniform");
	if (*flags & NORMAL_RNG) {
	    printf(" using %s for normals\n", (*flags & BOX_MULLER) ?
		   "Box-Muller" : "Ziggurat");
	}
	putchar('\n');
    }

    return err;
}

/* arguments: 

   --size={0,1,2} : for SmallCrush, Crush, BigCrush 
                    (default is Small)

   --normal : base test on gretl's normal RNG (default, use
              uniform RNG)

   --box-muller : when using normal RNG, use Box-Muller
                  instead of Ziggurat (Ziggurat is default)
*/

int main (int argc, char **argv)
{
    unif01_Gen *gen;
    int testsize = 0;
    int flags = 0;
    int err;

    err = read_opts(argc, argv, &flags, &testsize);

    if (err) {
	exit(EXIT_FAILURE);
    } 

    gretl_rand_init();

    if (flags & BOX_MULLER) {
	gretl_rand_set_box_muller(1);
    }    

    if (flags & NORMAL_RNG) {
	gen = unif01_CreateExternGen01("libgretl", u01_from_normal);
    } else {
	gen = unif01_CreateExternGen01("libgretl", u01_direct);
    }	

    /* select the test battery */
    if (testsize == 0) {
	bbattery_SmallCrush(gen);
    } else if (testsize == 1) {
	bbattery_Crush(gen);
    } else if (testsize == 2) {
	bbattery_BigCrush(gen);
    }

    unif01_DeleteExternGen01(gen); 
    gretl_rand_free();

    return 0;
}
