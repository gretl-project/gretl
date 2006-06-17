/* Grab some sorts of command info from libgretl, and write it out in
   tabular form.  Designed to keep certain parts of the manual up to
   date automatically.

   Allin Cottrell, June 2006.
*/

#include "generate.c"

#define COLS 8

void print_tabsep (int *n)
{
    *n += 1;

    if (*n == COLS) {
	fputs(" \\\\\n", stdout);
	*n = 0;
    } else {
	fputs(" & ", stdout);
    }
}

void print_tabtop (void)
{
    int i;

    fputs("\\begin{tabular}{", stdout);
    for (i=0; i<COLS; i++) {
	putchar('l');
    }
    fputs("}\n", stdout);
}

void print_tabfoot (int n)
{
    if (n < COLS) {
	fputs("\\\\\n", stdout);
    }

    fputs("\\end{tabular}\n\n", stdout);
}

void print_constants (void)
{
    int n1 = sizeof res1 / sizeof res1[0];
    int i, n = 0;

    print_tabtop();

    for (i=0; i<n1; i++) {
	printf("%s", res1[i]);
	print_tabsep(&n);
    }	

    print_tabfoot(n);
}

void print_internals (void)
{
    int n2 = sizeof res2 / sizeof res2[0];
    int n3 = sizeof res3 / sizeof res3[0];
    int i, n = 0;

    print_tabtop();

    for (i=0; i<n2; i++) {
	printf("%s", res2[i]);
	print_tabsep(&n);
    }

    for (i=0; i<n3; i++) {
	printf("%s", res3[i]);
	print_tabsep(&n);
    }

    print_tabfoot(n);
}

void print_func_words (void)
{
    int i, n = 0;

    print_tabtop();

    for (i=0; funcs[i].fnum != 0; i++) {
	printf("%s", funcs[i].fword);
	print_tabsep(&n);
    }

    print_tabfoot(n);
}

void print_loop_commands (void)
{
    int i, n = 0;

    print_tabtop();

    for (i=0; i<NC; i++) {
	if (ok_in_loop(i)) {
	    printf("%s", gretl_command_word(i));
	    print_tabsep(&n);
	}
    }

    print_tabfoot(n);
}

enum {
    NOTHING,
    CONSTANTS,
    INTERNALS,
    FUNCTIONS,
    LOOPCMDS
};    

int ok_opt (const char *str)
{
    const char *opts[] = {
	"--constants",
	"--internals",
	"--functions",
	"--loopcmds",
	NULL
    };
    int i;

    for (i=0; opts[i] != NULL; i++) {
	if (!strcmp(str, opts[i])) {
	    return i+1;
	}
    }
    
    return 0;
}

static void usage (const char *prog)
{
    fprintf(stderr, "%s: needs one valid option\n", prog);
    exit(EXIT_FAILURE);

}

int main (int argc, char **argv)
{
    int opt = 0;

    if (argc != 2) {
	usage(argv[0]);
    }

    opt = ok_opt(argv[1]);
    if (opt == 0) {
	usage(argv[0]);
    }

    if (opt == CONSTANTS) {
	print_constants();
    } else if (opt == INTERNALS) {
	print_internals();
    } else if (opt == FUNCTIONS) {
	print_func_words();
    } else if (opt == LOOPCMDS) {
	print_loop_commands();
    } else {
	/* impossible */
	usage(argv[0]);
    }

    return 0;
}
