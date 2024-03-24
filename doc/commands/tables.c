/* Grab some sorts of command info from libgretl, and write it out in
   tabular form.  Designed to keep certain parts of the manual up to
   date automatically.

   Allin Cottrell, June 2006.
*/

#include "genmain.c" /* need access to some private stuff */
#include "gen_public.h"
#include "monte_carlo.h"

int sort_strings (const void *a, const void *b)
{
    const char **sa = (const char **) a;
    const char **sb = (const char **) b;

    return strcmp(*sa, *sb);
}

int push_string_on_array (char ***arr, const char *s, int i)
{
    char **S = realloc(*arr, (i + 1) * sizeof *S);

    if (S == NULL) {
        return 1;
    }

    S[i] = gretl_strdup(s);
    if (S[i] == NULL) {
        return 1;
    }

    *arr = S;

    return 0;
}

void print_tabsep (int cols, int *n)
{
    *n += 1;

    if (*n == cols) {
        fputs(" \\\\\n", stdout);
        *n = 0;
    } else {
        fputs(" & ", stdout);
    }
}

void print_tabtop (int cols, int small)
{
    int i;

    if (small) {
        fputs("{\\small\n", stdout);
    }
    fputs("\\begin{tabular}{", stdout);
    for (i=0; i<cols; i++) {
        putchar('l');
    }
    fputs("}\n", stdout);
}

void print_tabfoot (int cols, int n, int small)
{
    if (n < cols) {
        fputs("\\\\\n", stdout);
    }
    if (small) {
        fputs("\\end{tabular}\n}\n\n", stdout);
    } else {
        fputs("\\end{tabular}\n\n", stdout);
    }
}

void sort_and_print_text (char **S, int n)
{
    int i;

    qsort(S, n, sizeof *S, sort_strings);

    for (i=0; i<n; i++) {
        printf("\\texttt{%s}", S[i]);
        if (i < n - 1) {
            fputs(", ", stdout);
        } else {
            fputs(".\n\n", stdout);
        }
    }

    strings_array_free(S, n);
}

void print_text_unsorted (char **S, int n)
{
    int i;

    for (i=0; i<n; i++) {
        printf("\\texttt{%s}", S[i]);
        if (i < n - 1) {
            fputs(", ", stdout);
        } else {
            fputs(".\n\n", stdout);
        }
    }

    strings_array_free(S, n);
}

/* Allow for the possibility that a function name
   contains an underscore */

char *tex_escape (char *targ, const char *src)
{
    char *ret = targ;

    while (*src) {
        if (*src == '_') {
            *targ++ = '\\';
            *targ++ = '_';
        } else {
            *targ++ = *src;
        }
        src++;
    }

    *targ = '\0';

    return ret;
}

void sort_and_print_tabular (char **S, int n, int cols,
                             int small)
{
    char tmp[32];
    int i, t = 0;

    qsort(S, n, sizeof *S, sort_strings);

    if (n < cols) {
        cols = n;
    }

    print_tabtop(cols, small);
    for (i=0; i<n; i++) {
        printf("\\texttt{%s}", tex_escape(tmp, S[i]));
        print_tabsep(cols, &t);
    }
    print_tabfoot(cols, t, small);

    strings_array_free(S, n);
}

void print_internals (void)
{
    int nr = sizeof reswords / sizeof reswords[0];
    char **S = NULL;
    int i, n = 0;
    int err = 0;

    for (i=0; i<nr && !err; i++) {
        err = push_string_on_array(&S, reswords[i], n++);
    }

    if (!err) {
        print_text_unsorted(S, n);
    }
}

void print_func_words (void)
{
    char **S = NULL;
    const char *s;
    int i, n, nf;
    int err = 0;

    nf = gen_func_count();

    for (i=0, n=0; i<nf && !err; i++) {
        s = gen_func_name(i);
        if (*s != '_') {
            err = push_string_on_array(&S, s, n++);
        }
    }

    if (!err) {
        sort_and_print_tabular(S, n, 8, 1);
    }
}

int print_loop_commands (void)
{
    char **S = NULL;
    int i, n = 0;
    int err = 0;

    for (i=0; i<NC && !err; i++) {
        if (ok_in_loop(i) && !(HIDDEN_COMMAND(i))) {
            err = push_string_on_array(&S, gretl_command_word(i), n++);
        }
    }

    if (!err) {
        sort_and_print_tabular(S, n, 8, 0);
    }

    return err;
}

int print_non_loop_commands (void)
{
    char **S = NULL;
    int i, n = 0;
    int err = 0;

    for (i=1; i<NC && !err; i++) {
        if (!ok_in_loop(i) && !(HIDDEN_COMMAND(i))) {
            err = push_string_on_array(&S, gretl_command_word(i), n++);
        }
    }

    if (!err) {
        sort_and_print_tabular(S, n, 8, 0);
    }

    return err;
}

enum {
    NOTHING,
    INTERNALS,
    FUNCTIONS,
    LOOPCMDS,
    NONLOOPCMDS,
};

int ok_opt (const char *str)
{
    const char *opts[] = {
        "--internals",
        "--functions",
        "--loopcmds",
        "--nonloopcmds",
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

    if (opt == INTERNALS) {
        print_internals();
    } else if (opt == FUNCTIONS) {
        print_func_words();
    } else if (opt == LOOPCMDS) {
        print_loop_commands();
    } else if (opt == NONLOOPCMDS) {
        print_non_loop_commands();
    } else {
        /* impossible */
        usage(argv[0]);
    }

    return 0;
}
