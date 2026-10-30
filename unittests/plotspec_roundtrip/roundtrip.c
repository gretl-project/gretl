/* Round-trip test for the gnuplot-script reader in libgretl
   (plotspec_read_file() in lib/src/plotspec_read.c) and the
   corresponding writer, plotspec_print() in lib/src/plotspec.c.

   For each gnuplot script named on the command line (as written by
   libgretl) we do the following:

   1. Read the script into a GPT_SPEC, S1, and write S1 out as script
      P1.
   2. Read P1 into S2 and write S2 out as P2.
   3. Check that S1 and S2 have the same basic structure, and that
      P1 and P2 are identical.
   4. If gnuplot is available, render the original script and P1 on
      gnuplot's text-only "dumb" terminal and check that the results
      are identical.

   Step 3 does not demand that P1 reproduce the original script,
   which it generally will not (comments, ordering of "set" lines and
   so on will differ), but it does require that the reader and writer
   are consistent: whatever the writer emits, the reader must be able
   to take back in without loss.

   Step 3 alone cannot catch information that is dropped by the
   reader (or the writer) on the first pass, since the loss is then
   the same on the second pass. That is what step 4 is for: it
   compares P1 with the original in terms of what gnuplot draws.

   Usage: roundtrip [-v] [-k] [-x name] file1.plt [file2.plt ...]

   -v  report successes as well as failures
   -x  expect the named file (basename) to fail; may be repeated
   -k  keep the intermediate files P1 and P2 (in a scratch directory
       whose name is printed) for inspection

   The exit code is non-zero if any file fails.
*/

#include "libgretl.h"
#include "plotspec.h"
#include "graphing.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static int verbose;
static int keep;
static int have_gnuplot;
static const char *xfail[32];
static int n_xfail;

static int read_spec (GPT_SPEC *spec, const char *fname)
{
    PlotReadInfo info = {0};

    if (strlen(fname) >= sizeof spec->fname) {
        return E_DATA;
    }
    strcpy(spec->fname, fname);
    gretl_error_clear();

    return plotspec_read_file(spec, NULL, &info);
}

static int write_spec (GPT_SPEC *spec, const char *fname)
{
    FILE *fp = fopen(fname, "wb");
    int err;

    if (fp == NULL) {
        return E_FOPEN;
    }
    err = plotspec_print(spec, fp);
    if (fclose(fp) != 0 && !err) {
        err = E_FOPEN;
    }

    return err;
}

static char *slurp (const char *fname, size_t *len)
{
    FILE *fp = fopen(fname, "rb");
    char *buf = NULL;
    long n;

    if (fp == NULL) {
        return NULL;
    }
    fseek(fp, 0, SEEK_END);
    n = ftell(fp);
    rewind(fp);
    buf = malloc(n + 1);
    if (buf != NULL) {
        *len = fread(buf, 1, n, fp);
        buf[*len] = '\0';
    }
    fclose(fp);

    return buf;
}

static int same_double (double x, double y)
{
    return (na(x) && na(y)) || x == y;
}

/* the writer omits blank strings, so NULL and "" count as equal */

static int same_string (const char *s, const char *t)
{
    if (string_is_blank(s) || string_is_blank(t)) {
        return string_is_blank(s) && string_is_blank(t);
    }
    return strcmp(s, t) == 0;
}

/* compare the "shape" of two specs: report the first mismatch */

static const char *spec_mismatch (const GPT_SPEC *a, const GPT_SPEC *b)
{
    int i;

    if (a->code != b->code) return "plot type";
    if (a->nobs != b->nobs) return "nobs";
    if (a->n_datacols != b->n_datacols) return "n_datacols";
    if (a->n_lines != b->n_lines) return "n_lines";
    if (a->n_labels != b->n_labels) return "n_labels";
    if (a->n_arrows != b->n_arrows) return "n_arrows";
    if (a->n_literal != b->n_literal) return "n_literal";
    if (a->flags != b->flags) return "flags";
    if (a->fit != b->fit) return "fit";
    if (a->pd != b->pd) return "pd";

    for (i=0; i<5; i++) {
        if (!same_string(a->titles[i], b->titles[i])) return "titles";
    }
    for (i=0; i<6; i++) {
        if (!same_double(a->range[i][0], b->range[i][0]) ||
            !same_double(a->range[i][1], b->range[i][1])) return "range";
    }
    for (i=0; i<3; i++) {
        if (a->logbase[i] != b->logbase[i]) return "logbase";
    }
    for (i=0; i<a->n_lines; i++) {
        if (a->lines[i].ncols != b->lines[i].ncols) return "line ncols";
        if (a->lines[i].style != b->lines[i].style) return "line style";
        if (a->lines[i].type != b->lines[i].type) return "line type";
        if (a->lines[i].yaxis != b->lines[i].yaxis) return "line yaxis";
        if (!same_string(a->lines[i].title, b->lines[i].title)) return "line title";
    }

    return NULL;
}

static int is_xfail (const char *fname)
{
    const char *base = strrchr(fname, '/');
    int i;

    base = (base != NULL) ? base + 1 : fname;
    for (i=0; i<n_xfail; i++) {
        if (!strcmp(base, xfail[i])) {
            return 1;
        }
    }

    return 0;
}

/* Run gnuplot on @script, using the "dumb" (text) terminal, with the
   output going to @outname. Any "set term" and "set output" lines in
   the script itself are dropped. The script is assumed to contain
   its own inline data.
*/

static int render_dumb (const char *script, const char *outname,
                        const char *tmpdir)
{
    char wname[FILENAME_MAX], cmd[3 * FILENAME_MAX];
    char line[1024];
    FILE *fi, *fo;

    fi = fopen(script, "r");
    if (fi == NULL) {
        return E_FOPEN;
    }
    sprintf(wname, "%s/render.gp", tmpdir);
    fo = fopen(wname, "w");
    if (fo == NULL) {
        fclose(fi);
        return E_FOPEN;
    }

    fprintf(fo, "set terminal dumb size 110,40\n");
    fprintf(fo, "set output '%s'\n", outname);
    while (fgets(line, sizeof line, fi)) {
        if (!strncmp(line, "set term ", 9) ||
            !strncmp(line, "set terminal ", 13) ||
            !strncmp(line, "set output", 10)) {
            continue;
        }
        fputs(line, fo);
    }
    fputs("unset output\n", fo);
    fclose(fi);
    fclose(fo);

    sprintf(cmd, "GNUTERM=dumb gnuplot '%s' > '%s.err' 2>&1", wname, outname);

    return system(cmd) != 0;
}

/* returns 0 if gnuplot renders the two scripts identically */

static int compare_rendering (const char *orig, const char *p1,
                              const char *tmpdir, const char **why)
{
    char o1[FILENAME_MAX], o2[FILENAME_MAX];
    char *t1, *t2;
    size_t n1 = 0, n2 = 0;
    int err = 0;

    sprintf(o1, "%s/orig.txt", tmpdir);
    sprintf(o2, "%s/p1.txt", tmpdir);

    if (render_dumb(orig, o1, tmpdir)) {
        *why = "gnuplot failed on original";
        return 1;
    }
    if (render_dumb(p1, o2, tmpdir)) {
        *why = "gnuplot failed on P1";
        return 1;
    }
    t1 = slurp(o1, &n1);
    t2 = slurp(o2, &n2);
    if (t1 == NULL || t2 == NULL || n1 == 0) {
        *why = "no gnuplot output";
        err = 1;
    } else if (n1 != n2 || memcmp(t1, t2, n1) != 0) {
        *why = "gnuplot renders original and P1 differently";
        err = 1;
    }
    free(t1);
    free(t2);

    return err;
}

static int check_file (const char *fname, const char *tmpdir)
{
    char p1[FILENAME_MAX], p2[FILENAME_MAX];
    GPT_SPEC *s1 = NULL, *s2 = NULL;
    char *t1 = NULL, *t2 = NULL;
    size_t n1 = 0, n2 = 0;
    const char *why = NULL;
    int skipped = 0;
    int err = 0;

    sprintf(p1, "%s/p1.plt", tmpdir);
    sprintf(p2, "%s/p2.plt", tmpdir);

    s1 = plotspec_new();
    s2 = plotspec_new();
    if (s1 == NULL || s2 == NULL) {
        err = E_ALLOC;
        why = "allocation";
        goto done;
    }

    err = read_spec(s1, fname);
    if (err) {
        why = "reading original";
        goto done;
    }
    if (cant_edit(s1->code)) {
        /* the reader deliberately stops short for these (the GUI
           plot editor does not support them), so there is nothing
           to round-trip */
        if (verbose) {
            printf("skip %s (not editable)\n", fname);
        }
        skipped = 1;
        goto done;
    }
    err = write_spec(s1, p1);
    if (err) {
        why = "writing P1";
        goto done;
    }
    err = read_spec(s2, p1);
    if (err) {
        why = "reading P1";
        goto done;
    }
    err = write_spec(s2, p2);
    if (err) {
        why = "writing P2";
        goto done;
    }

    why = spec_mismatch(s1, s2);
    if (why != NULL) {
        err = 1;
        goto done;
    }

    t1 = slurp(p1, &n1);
    t2 = slurp(p2, &n2);
    if (t1 == NULL || t2 == NULL) {
        err = E_FOPEN;
        why = "slurping output";
    } else if (n1 != n2 || memcmp(t1, t2, n1) != 0) {
        err = 1;
        why = "P1 and P2 differ";
    }

    if (!err && have_gnuplot) {
        err = compare_rendering(fname, p1, tmpdir, &why);
    }

 done:

    if (is_xfail(fname)) {
        /* a known problem: failure is expected, success is news */
        if (err) {
            printf("XFAIL %s: %s\n", fname, why);
            err = 0;
        } else {
            printf("XPASS %s: expected to fail\n", fname);
            err = 1;
        }
    } else if (err) {
        const char *msg = gretl_errmsg_get();

        printf("FAIL %s: %s (err=%d%s%s)\n", fname, why, err,
               *msg ? ": " : "", msg);
    } else if (verbose && !skipped) {
        printf("ok   %s (%d lines, %d obs)\n", fname, s1->n_lines, s1->nobs);
    }

    plotspec_destroy(s1);
    plotspec_destroy(s2);
    free(t1);
    free(t2);

    return err;
}

int main (int argc, char **argv)
{
    char tmpdir[] = "/tmp/plotspec_rt_XXXXXX";
    char p[FILENAME_MAX];
    int nfail = 0, nfiles = 0;
    int i;

    while (argc > 1 && argv[1][0] == '-') {
        if (!strcmp(argv[1], "-v")) {
            verbose = 1;
        } else if (!strcmp(argv[1], "-k")) {
            keep = 1;
        } else if (!strcmp(argv[1], "-x") && argc > 2 && n_xfail < 32) {
            xfail[n_xfail++] = argv[2];
            argc--;
            argv++;
        } else {
            break;
        }
        argc--;
        argv++;
    }
    if (argc < 2) {
        fprintf(stderr, "usage: roundtrip [-v] [-k] [-x name] file.plt ...\n");
        return 2;
    }
    have_gnuplot = system("command -v gnuplot >/dev/null 2>&1") == 0;
    if (!have_gnuplot) {
        printf("gnuplot not found: skipping the rendering comparison\n");
    }
    if (mkdtemp(tmpdir) == NULL) {
        perror("mkdtemp");
        return 2;
    }

    libgretl_init();

    for (i=1; i<argc; i++) {
        nfail += check_file(argv[i], tmpdir) != 0;
        nfiles++;
    }

    printf("%d files, %d failed\n", nfiles, nfail);

    libgretl_cleanup();

    if (keep) {
        printf("intermediate files left in %s (those of the last file)\n", tmpdir);
    } else {
        sprintf(p, "%s/p1.plt", tmpdir);
        remove(p);
        sprintf(p, "%s/p2.plt", tmpdir);
        remove(p);
        sprintf(p, "%s/render.gp", tmpdir);
        remove(p);
        sprintf(p, "%s/orig.txt", tmpdir);
        remove(p);
        sprintf(p, "%s/orig.txt.err", tmpdir);
        remove(p);
        sprintf(p, "%s/p1.txt", tmpdir);
        remove(p);
        sprintf(p, "%s/p1.txt.err", tmpdir);
        remove(p);
        rmdir(tmpdir);
    }

    return nfail > 0;
}
