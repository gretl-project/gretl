/*
 *  gretl -- Gnu Regression, Econometrics and Time-series Library
 *  Copyright (C) 2001 Allin Cottrell and Riccardo "Jack" Lucchetti
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "libgretl.h"
#include "libset.h"
#include "gretl_mt.h"
#include "gretl_func.h"
#include "gretl_foreign.h"
#include "matrix_extra.h"
#include "gretl_typemap.h"
#include "uservar.h"

#include <unistd.h>

#ifdef HAVE_MPI
# include "gretl_mpi.h"
# include "gretl_xml.h"
# ifndef G_OS_WIN32
#  define MPI_REALTIME 1 /* won't work on Windows */
# endif
#endif

#ifdef USE_RLIB
# include <Rinternals.h> /* for SEXP and friends */
#endif

#ifdef G_OS_WIN32
# include "gretl_win32.h"
#else
# include <sys/types.h>
# include <sys/stat.h>
# include <signal.h>
#endif

#define FDEBUG 0

static char **foreign_lines;
static int foreign_started;
static int foreign_n_lines;
static int foreign_lang;

/* foreign_opt may include OPT_D to send data, OPT_Q to operate
   quietly (don't display output from foreign program)
*/
static gretlopt foreign_opt;

/* dotdir filenames for R */
static gchar *gretl_Rprofile;
static gchar *gretl_Rsetup;
static gchar *gretl_Rsrc;
static gchar *gretl_Rout;
static gchar *gretl_Rmsg;

#ifdef HAVE_MPI
static gchar *gretl_mpi_script;
#endif

struct fmap {
    int lang;
    const char *scriptname;
    const char *iofile;
    gchar *scriptpath;
};

static void write_R_io_file (FILE *fp, const char *ddir);
static void do_stata_printout (PRN *prn);

#ifdef MPI_REALTIME

struct iodata {
    int fd;
    char buf[4096];
    int len;
    int finished;
    int got_all;
    int *err;
    PRN *prn;
};

static int run_mpi_with_pipes (char **argv, struct iodata *io,
                               gretlopt opt, PRN *prn);
#endif /* MPI_REALTIME */

static struct fmap foreign_map[] = {
     { LANG_OX,     "gretltmp.ox", "gretl_io.ox", NULL },
     { LANG_OCTAVE, "gretltmp.m",  "gretl_io.m", NULL },
     { LANG_STATA,  "gretltmp.do", "gretl_export.ado", NULL },
     { LANG_PYTHON, "gretltmp.py", "gretl_io.py", NULL },
     { LANG_JULIA,  "gretltmp.jl", "gretl_io.jl", NULL }
};

static const char *get_io_filename (int lang)
{
    int i, n = G_N_ELEMENTS(foreign_map);

    for (i=0; i<n; i++) {
        if (foreign_map[i].lang == lang) {
            return foreign_map[i].iofile;
        }
    }

    return NULL;
}

static gchar *get_foreign_scriptpath (int lang)
{
    int i, n = G_N_ELEMENTS(foreign_map);

    for (i=0; i<n; i++) {
        if (foreign_map[i].lang == lang) {
            if (foreign_map[i].scriptpath == NULL) {
                foreign_map[i].scriptpath =
                    gretl_make_dotpath(foreign_map[i].scriptname);
            }
            return foreign_map[i].scriptpath;
        }
    }

    return NULL;
}

static void delete_foreign_script (int lang)
{
    int i, n = G_N_ELEMENTS(foreign_map);

    for (i=0; i<n; i++) {
        if (foreign_map[i].lang == lang) {
            if (foreign_map[i].scriptpath != NULL) {
                gretl_remove(foreign_map[i].scriptpath);
            }
            break;
        }
    }
}

void foreign_destroy (void)
{
    if (foreign_lines != NULL) {
        strings_array_free(foreign_lines, foreign_n_lines);
        foreign_lines = NULL;
    }
    foreign_started = 0;
    foreign_n_lines = 0;
    foreign_opt = OPT_NONE;
}

/* Get the user's "dotdir" in a form suitable for writing
   into files to be read by third-party programs. On
   Windows, such programs will presumably expect the
   directory name to be in the locale encoding. In most
   cases (but Stata?) they will also want forward slashes.

   On platforms other than Windows we just pass a copy of
   dotdir as is.
*/

static const gchar *get_export_dotdir (void)
{
    static gchar *fdot;

    if (fdot == NULL) {
        fdot = g_strdup(gretl_dotdir());
#ifdef G_OS_WIN32
        /* recode to locale if necessary */
        if (utf8_encoded(fdot)) {
            gsize bytes;
            gchar *locdot;

            locdot = g_locale_from_utf8(fdot, -1, NULL, &bytes, NULL);
            if (locdot != NULL) {
                g_free(fdot);
                fdot = locdot;
            }
        }

        /* ensure forward slashes? is stata OK with this? */
        if (1) {
            char *s = fdot;

            while (*s) {
                if (*s == '\\') {
                    *s = '/';
                }
                s++;
            }
        }
#endif
    }

    return fdot;
}

static int set_foreign_lang (const char *lang, PRN *prn)
{
    int err = 0;

    if (g_ascii_strcasecmp(lang, "R") == 0) {
        foreign_lang = LANG_R;
    } else if (g_ascii_strcasecmp(lang, "ox") == 0) {
        foreign_lang = LANG_OX;
    } else if (g_ascii_strcasecmp(lang, "octave") == 0) {
        foreign_lang = LANG_OCTAVE;
    } else if (g_ascii_strcasecmp(lang, "stata") == 0) {
        foreign_lang = LANG_STATA;
    } else if (g_ascii_strcasecmp(lang, "python") == 0) {
        foreign_lang = LANG_PYTHON;
    } else if (g_ascii_strcasecmp(lang, "julia") == 0) {
        foreign_lang = LANG_JULIA;
    } else if (g_ascii_strcasecmp(lang, "mpi") == 0) {
#ifdef HAVE_MPI
        if (gretl_mpi_initialized()) {
            gretl_errmsg_set(_("MPI is already initialized"));
            err = 1;
        } else {
            foreign_lang = LANG_MPI;
        }
#else
        gretl_errmsg_set(_("MPI is not supported in this gretl build"));
        err = E_NOTIMP;
#endif
    } else {
        pprintf(prn, _("%s: unknown language\n"), lang);
        err = E_DATA;
    }

    return err;
}

static int lib_run_prog_sync (char **argv, gretlopt opt,
                              int lang, PRN *prn)
{
    gchar *sout = NULL;
    gchar *errout = NULL;
    gint status = 0;
    GError *gerr = NULL;
    int err = 0;

    g_spawn_sync(gretl_workdir(), argv,
                 NULL, G_SPAWN_SEARCH_PATH,
                 NULL, NULL, &sout, &errout,
                 &status, &gerr);

    if (gerr != NULL) {
        pprintf(prn, "g_spawn_sync: %s\n", gerr->message);
        g_error_free(gerr);
        err = 1;
    } else if (status != 0) {
        pprintf(prn, _("%s exited with status %d\n"), argv[0], status);
        if (sout != NULL && *sout != '\0') {
            pputs(prn, "stdout:\n");
            pputs(prn, sout);
            pputc(prn, '\n');
        }
        if (errout != NULL && *errout != '\0') {
            pprintf(prn, "\nstderr from %s:\n", argv[0]);
            pputs(prn, errout);
            pputc(prn, '\n');
        }
        err = 1;
    } else if (sout != NULL) {
        if (lang == LANG_MPI || !(opt & OPT_Q)) {
            /* with OPT_Q, don't print non-error output,
               unless we're running MPI
            */
            if (foreign_lang == LANG_STATA) {
                do_stata_printout(prn);
            } else if (*sout != '\0') {
                pputs(prn, sout);
                if (sout[strlen(sout) -1] != '\n') {
                    pputc(prn, '\n');
                }
            }
        }
        if (opt & OPT_V) {
            /* also print stderr output, if any */
            if (errout != NULL && *errout != '\0') {
                pprintf(prn, "\nstderr from %s:\n", argv[0]);
                pputs(prn, errout);
                pputc(prn, '\n');
            }
        }
    } else {
        pprintf(prn, "%s: %s\n", argv[0], _("Got no output"));
        err = 1;
    }

    g_free(sout);
    g_free(errout);

    return err;
}

#ifdef HAVE_MPI /* start common MPI-driver code block */

enum {
    MPI_OPENMPI,
    MPI_MPICH,
    MPI_MSMPI
};

/* OS-specific default values */
# ifdef G_OS_WIN32
static int mpi_variant = MPI_MSMPI;
# else
static int mpi_variant = MPI_OPENMPI;
# endif

void set_mpi_variant (const char *pref)
{
    if (!strcmp(pref, "OpenMPI")) {
        mpi_variant = MPI_OPENMPI;
    } else if (!strcmp(pref, "MPICH")) {
        mpi_variant = MPI_MPICH;
    } else if (strcmp(pref, "MS-MPI")) {
        mpi_variant = MPI_MSMPI;
    }
}

int gretl_max_mpi_processes (void)
{
    const char *hostfile = gretl_mpi_hosts();
    int procmax = gretl_n_processors();

    if (hostfile != NULL && *hostfile != '\0') {
        FILE *fp = gretl_fopen(hostfile, "r");

        if (fp != NULL) {
            const char *fmt;
            char line[256], host[128];
            int nf, slots, allslots = 0;
            int err = 0;

            if (mpi_variant == MPI_MSMPI) {
                fmt = "%127s %d";
            } else if (mpi_variant == MPI_MPICH) {
                fmt = "%127[^:]:%d";
            } else {
                fmt = "%127s slots=%d";
            }

            while (fgets(line, sizeof line, fp) && !err) {
                if (*line != '#' && !string_is_blank(line)) {
                    nf = sscanf(line, fmt, host, &slots);
                    if (nf == 2) {
                        allslots += slots;
                    } else {
                        err = E_DATA;
                    }
                }
            }

            if (!err && allslots > 0) {
                procmax = allslots;
            }

            fclose(fp);
        }
    }

    return procmax;
}

int check_for_mpiexec (void)
{
    const char *prog = gretl_mpiexec();

    return check_for_program(prog);
}

static const gchar *get_mpi_scriptname (void)
{
    if (gretl_mpi_script == NULL) {
        gretl_mpi_script = gretl_make_dotpath("gretltmp-mpi.inp");
    }

    return gretl_mpi_script;
}

static const char *get_mpi_hostfile (void)
{
    const char *hfile = gretl_mpi_hosts();

    if (*hfile == '\0') {
        hfile = getenv("GRETL_MPI_HOSTS");
    }

    return hfile;
}

/* handle the number-of-processes MPI option */

static int get_user_mpi_procs (int *err)
{
    int np = get_optval_int(MPI, OPT_N, err);

    if (!*err && (np <= 0 || np > 9999999)) {
        *err = E_DATA;
    }

    return np;
}

# ifdef G_OS_WIN32

static int win32_run_mpi_sync (gretlopt opt, PRN *prn)
{
    const char *hostfile = get_mpi_hostfile();
    char np = 0;
    int err = 0;

    if (opt & OPT_N) {
        np = get_user_mpi_procs(&err);
    }

    if (!err) {
        GString *cmd = g_string_new(NULL);

        g_string_append_printf(cmd, "\"%s\"", gretl_mpiexec());

        if (!(opt & OPT_L) && hostfile != NULL && *hostfile != '\0') {
            /* note: OPT_L corresponds to --local, meaning that we
               should not use a hosts file even if one is present
            */
            g_string_append(cmd, " /machinefile ");
            g_string_append(cmd, hostfile);
        } else if (np == 0) {
            /* no user spec, so supply a default np value */
            if (libset_get_bool(MPI_USE_SMT)) {
                /* use max number of processes */
                np = gretl_n_processors();
            } else {
                /* don't use hyper-threads */
                np = gretl_n_physical_cores();
            }
        }
        if (np > 0) {
            g_string_append_printf(cmd, " /np %d", np);
        }
        g_string_append_printf(cmd, " \"%sgretlmpi\"", gretl_bindir());
        if (opt & OPT_S) {
            g_string_append(cmd, " --single-rng");
        }
        if (opt & OPT_Q) {
            g_string_append(cmd, " --quiet");
        }
        if (opt & OPT_K) {
	    const char *s = optval_string(MPI, OPT_K);

            g_string_append_printf(cmd, " --key=\"%s\"", s);
        }
        g_string_append_printf(cmd, " \"%s\"", get_mpi_scriptname());

        if (opt & OPT_V) {
            pprintf(prn, "gretl mpi command:\n %s\n", cmd->str);
        }

        err = gretl_win32_pipe_output(cmd->str, gretl_workdir(), prn);
        g_string_free(cmd, TRUE);
    }

    return err;
}

# else /* non-Windows MPI-calling variant */

/* The following should probably be redundant on Linux but
   may be needed for the macOS package, where the gretl bin
   directory may not be in PATH.
*/

static gchar *gretl_mpi_binary (void)
{
    gchar *tmp = g_strdup(gretl_home());
    gchar *p = strstr(tmp, "/share/gretl");
    gchar *ret = NULL;

    if (p != NULL) {
        *p = '\0';
        ret = g_strdup_printf("%s/bin/gretlmpi", tmp);
    } else {
        ret = g_strdup("gretlmpi");
    }

    g_free(tmp);

    return ret;
}

/* called when gretlmpi is run with the --verbose flag */

static void print_mpi_command (char **argv, PRN *prn)
{
    int i;

    pputs(prn, "gretl mpi command:\n ");
    for (i=0; argv[i] != NULL; i++) {
        pprintf(prn, "%s ", argv[i]);
    }
    pputc(prn, '\n');
}

static int lib_run_mpi_sync (gretlopt opt, void *ptr, PRN *prn)
{
    const char *hostfile = get_mpi_hostfile();
    int use_valgrind = 0;
    char np = 0;
    int err = 0;

    if (opt & OPT_N) {
        np = get_user_mpi_procs(&err);
    }

    if (getenv("MPI_VGRIND") != NULL) {
        use_valgrind = 1;
    }

    if (!err) {
        const char *mpiexec = gretl_mpiexec();
        gchar *mpiprog = gretl_mpi_binary();
	gchar *keystr = NULL;
        const char *hostsopt = NULL;
        char *argv[14] = {0};
        char npnum[8] = {0};
        int nproc, i = 0;

        nproc = gretl_n_processors();

        if (!(opt & OPT_L) && hostfile != NULL && *hostfile != '\0') {
            /* note: OPT_L corresponds to --local, meaning that we
               should not use a hosts file even if one is present
            */
            if (mpi_variant == MPI_MPICH) {
                hostsopt = "-machinefile";
            } else {
                hostsopt = "--hostfile";
            }
        } else if (np == 0) {
            /* no user spec, so supply a default np value */
            if (libset_get_bool(MPI_USE_SMT)) {
                /* use max number of processes */
                np = nproc;
            } else {
                /* don't use hyper-threads */
                np = gretl_n_physical_cores();
            }
        }

        argv[i++] = (char *) mpiexec;
        if (hostsopt != NULL) {
            argv[i++] = (char *) hostsopt;
            argv[i++] = (char *) hostfile;
        }
        if (np > 0) {
            sprintf(npnum, "%d", np);
            argv[i++] = "-np";
            argv[i++] = npnum;
            if (mpi_variant == MPI_OPENMPI && np > nproc/2) {
                argv[i++] = "--oversubscribe";
            }
        }
        if (use_valgrind) {
            argv[i++] = "valgrind";
            argv[i++] = "--leak-check=full";
        }
        argv[i++] = mpiprog;
        if (opt & OPT_S) {
            argv[i++] = "--single-rng";
        }
        if (opt & OPT_Q) {
            argv[i++] = "--quiet";
        }
        if (opt & OPT_K) {
	    const char *s = get_optval_string(MPI, OPT_K);

	    keystr = g_strdup_printf("--key=%s", s);
            argv[i++] = keystr;
        }
        argv[i++] = (char *) get_mpi_scriptname();
        argv[i] = NULL;

        if (opt & OPT_V) {
            print_mpi_command(argv, prn);
        }

#ifdef MPI_REALTIME
        err = run_mpi_with_pipes(argv, ptr, opt, prn);
#else
        err = lib_run_prog_sync(argv, opt, LANG_MPI, prn);
#endif
        g_free(mpiprog);
	g_free(keystr);
    }

    return err;
}

# endif /* Windows vs not */

# ifdef MPI_REALTIME

/* This approach works on Linux */

static void mpi_childwatch (GPid pid, gint status, gpointer p)
{
    struct iodata *io = p;
    GError *gerr = NULL;

    if (!g_spawn_check_exit_status(status, &gerr)) {
        pprintf(io->prn, "gretlmpi: %s\n", gerr->message);
        *(io->err) = 1;
        io->got_all = 1;
    }
    g_spawn_close_pid(pid);
    io->finished = 1;
}

static void relay_mpi_output (struct iodata *io)
{
    int got;

    memset(io->buf, 0, io->len);
    got = read(io->fd, io->buf, io->len - 1);

    if (got > 0) {
        char *s = strstr(io->buf, "__GRETLMPI_EXIT__");

        if (s != NULL) {
            io->got_all = 1;
            *s = '\0';
        }
        pputs(io->prn, io->buf);
        gretl_flush(io->prn);
    }
}

static int run_mpi_with_pipes (char **argv, struct iodata *io,
                               gretlopt opt, PRN *prn)
{
    gint sout;
    GError *gerr = NULL;
    GPid child_pid;
    int err = 0;

    g_spawn_async_with_pipes(gretl_workdir(),
                             argv,
                             NULL, /* envp */
                             G_SPAWN_SEARCH_PATH |
                             G_SPAWN_DO_NOT_REAP_CHILD,
                             NULL, /* child_setup */
                             NULL, /* data for child_setup */
                             &child_pid,
                             NULL, /* stdin */
                             &sout,
                             NULL,
                             &gerr);

    if (gerr != NULL) {
        pprintf(prn, "%s\n", gerr->message);
        g_error_free(gerr);
        err = 1;
    } else {
        io->fd = sout;
        io->len = sizeof io->buf;
        io->buf[0] = '\0';
        io->err = &err;
        io->prn = prn;
        io->finished = 0;
        io->got_all = 0;

        g_child_watch_add(child_pid, mpi_childwatch, io);

        while (!io->finished && !io->got_all) {
            relay_mpi_output(io);
            if (!io->got_all) {
                g_usleep(100000); /* 0.10 seconds */
            }
            g_main_context_iteration(NULL, FALSE);
        }
        while (!err && !io->got_all) {
            relay_mpi_output(io);
        }
        close(sout);
    }

    return err;
}

# endif /* MPI_REALTIME */

#else /* no MPI support */

int gretl_max_mpi_processes (void)
{
    return 0;
}

#endif /* HAVE_MPI or not */

/* special: print to @prn Stata's batch logfile */

static void do_stata_printout (PRN *prn)
{
    gchar *buf = NULL;

    /* we need to be located in the directory in which
       gretltmp.log is written at this point */
    gretl_chdir(gretl_workdir());

    if (g_file_get_contents("gretltmp.log", &buf, NULL, NULL)) {
        pputs(prn, buf);
        g_free(buf);
        pputc(prn, '\n');
    }

    gretl_remove("gretltmp.log");
}

static void make_gretl_R_names (void)
{
    static int done;

    if (!done) {
        const char *ddir = get_export_dotdir();

        gretl_Rprofile = g_strdup_printf("%sgretl.Rprofile", ddir);
        gretl_Rsetup = g_strdup_printf("%sgretl.Rsetup", ddir);
        gretl_Rsrc = g_strdup_printf("%sRsrc", ddir);
        gretl_Rout = g_strdup_printf("%sR.out", ddir);
        gretl_Rmsg = g_strdup_printf("%sR.msg", ddir);
        done = 1;
    }
}

#ifdef G_OS_WIN32 /* Windows specific */

static char *win32_get_rscript_path (void)
{
    const char *rbin = gretl_rbin_path();
    char *p, *rscript;
    int err = 0;

    rscript = calloc(strlen(rbin) + 16, 1);
    strcpy(rscript, rbin);
    p = strrchr(rscript, 'R');

    if (p != NULL) {
        *p = '\0';
        strcat(p, "Rscript.exe");
        err = gretl_stat(rscript, NULL);
    } else {
        err = 1;
    }

    if (err) {
        free(rscript);
        rscript = NULL;
    }

    return rscript;
}

/* FIXME Windows console? */

static void win32_put_R_output_line (const char *line, PRN *prn)
{
    if (gretl_in_gui_mode() && !g_utf8_validate(line, -1, NULL)) {
        gsize bytes;

        gchar *lconv = g_locale_to_utf8(line, -1, NULL,
                                        &bytes, NULL);
        if (lconv != NULL) {
            pputs(prn, lconv);
            g_free(lconv);
        } else {
            pputs(prn, _("line could not be converted to UTF-8\n"));
        }
    } else {
        pputs(prn, line);
    }
}

static int win32_lib_run_R_sync (gretlopt opt, PRN *prn)
{
    char *rscript = win32_get_rscript_path();
    gchar *cmd;
    int err = 0;

    /* ensure that we don't get stale output */
    gretl_remove(gretl_Rout);
    gretl_remove(gretl_Rmsg);

    /* Note that here we're calling R with gretl_Rprofile
       as an argument, as opposed to getting R to source
       it via the environment, since the latter seemed not
       to be working on Windows.
    */

    if (rscript != NULL) {
        cmd = g_strdup_printf("\"%s\" --vanilla \"%s\"", rscript,
                              gretl_Rprofile);
        free(rscript);
    } else {
        cmd = g_strdup_printf("\"%s\" CMD BATCH --no-save --no-init-file "
                              "--no-restore-data --slave \"%s\"",
                              gretl_rbin_path(), gretl_Rprofile);
    }

    err = win_run_sync(cmd, NULL);

#if FDEBUG
    fprintf(stderr, "win32_lib_run_R_sync: err = %d\n cmd='%s'\n", err, cmd);
#endif

    if (!(opt & OPT_Q)) {
        const gchar *outname;
        FILE *fp;

        outname = err ? gretl_Rmsg : gretl_Rout;
        fp = gretl_fopen(outname, "r");

        if (fp != NULL) {
            char line[1024];

            while (fgets(line, sizeof line, fp)) {
                win32_put_R_output_line(line, prn);
            }
            fclose(fp);
            gretl_remove(outname);
        }
    }

    g_free(cmd);

    return err;
}

static int win32_lib_run_other_sync (gretlopt opt, PRN *prn)
{
    const char *exe;
    const char *fname;
    gchar *cmd = NULL;
    int err;

    fname = get_foreign_scriptpath(foreign_lang);

    if (foreign_lang == LANG_OX) {
        exe = gretl_oxl_path();
        cmd = g_strdup_printf("\"%s\" \"%s\"", exe, fname);
    } else if (foreign_lang == LANG_OCTAVE) {
        exe = gretl_octave_path();
        cmd = g_strdup_printf("\"%s\" --silent -H \"%s\"", exe, fname);
    } else if (foreign_lang == LANG_STATA) {
        exe = gretl_stata_path();
        cmd = g_strdup_printf("\"%s\" /q /e do \"%s\"", exe, fname);
    } else if (foreign_lang == LANG_PYTHON) {
        exe = gretl_python_path();
        cmd = g_strdup_printf("\"%s\" \"%s\"", exe, fname);
    } else if (foreign_lang == LANG_JULIA) {
        exe = gretl_julia_path();
        if (opt & OPT_N) {
            cmd = g_strdup_printf("\"%s\" --compile=no \"%s\"", exe, fname);
        } else {
            cmd = g_strdup_printf("\"%s\" \"%s\"", exe, fname);
        }
    } else {
        return 1;
    }

    if (foreign_lang == LANG_OCTAVE) {
        /* try to suppress octave history mechanism */
        gretl_setenv("OCTAVE_HISTFILE", "c:\\nul");
    }

    err = gretl_win32_pipe_output(cmd, gretl_workdir(), prn);

    if (!err && foreign_lang == LANG_STATA && !(opt & OPT_Q)) {
        /* output will be in log file, not stdout */
        do_stata_printout(prn);
    }

    g_free(cmd);

    return err;
}

#else /* non-Windows code follows */

#ifdef OS_OSX

static int lib_run_R_sync (gretlopt opt, PRN *prn)
{
    gchar *argv[] = {
        NULL,
        "--no-save",
        "--no-init-file",
        "--no-restore-data",
        "--slave",
        NULL
    };
    int ret;

    argv[0] = g_find_program_in_path("R");
    if (argv[0] == NULL) {
        argv[0] = g_strdup("/Library/Frameworks/R.framework/Resources/bin/R");
    }
#if FDEBUG
    fprintf(stderr, "lib_run_R_sync: argv[0] = '%s'\n", argv[0]);
#endif
    ret = lib_run_prog_sync(argv, opt, LANG_R, prn);
    g_free(argv[0]);

    return ret;
}

#else

static int lib_run_R_sync (gretlopt opt, PRN *prn)
{
    char *argv[] = {
        "R",
        "--no-save",
        "--no-init-file",
        "--no-restore-data",
        "--slave",
        NULL
    };

    return lib_run_prog_sync(argv, opt, LANG_R, prn);
}

#endif /* macOS or not */

static int lib_run_other_sync (gretlopt opt, PRN *prn)
{
    char *scriptpath;
    char *argv[6];
    int err;

    scriptpath = (char *) get_foreign_scriptpath(foreign_lang);

    if (foreign_lang == LANG_OX) {
        argv[0] = (char *) gretl_oxl_path();
        argv[1] = scriptpath;
        argv[2] = NULL;
    } else if (foreign_lang == LANG_OCTAVE) {
        argv[0] = (char *) gretl_octave_path();
        argv[1] = "--silent";
        argv[2] = "-H";
        argv[3] = scriptpath;
        argv[4] = NULL;
    } else if (foreign_lang == LANG_PYTHON) {
        argv[0] = (char *) gretl_python_path();
        argv[1] = scriptpath;
        argv[2] = NULL;
    } else if (foreign_lang == LANG_JULIA) {
        argv[0] = (char *) gretl_julia_path();
        if (opt & OPT_N) {
            argv[1] = "--compile=no";
            argv[2] = scriptpath;
            argv[3] = NULL;
        } else {
            argv[1] = scriptpath;
            argv[2] = NULL;
        }
    } else if (foreign_lang == LANG_STATA) {
        argv[0] = (char *) gretl_stata_path();
        argv[1] = "-q";
        argv[2] = "-b";
        argv[3] = "do";
        argv[4] = scriptpath;
        argv[5] = NULL;
    }

    if (foreign_lang == LANG_OCTAVE) {
        /* suppress history mechanism */
        gretl_setenv("OCTAVE_HISTFILE", "/dev/null");
    }

    err = lib_run_prog_sync(argv, opt, foreign_lang, prn);

    return err;
}

#endif /* end non-Windows block */

static FILE *write_open_dotfile (const char *fname)
{
    gchar *path = gretl_make_dotpath(fname);
    FILE *fp;

    fp = gretl_fopen(path, "w");
    g_free(path);

    return fp;
}

static int dotfile_exists (const char *fname)
{
    gchar *path = gretl_make_dotpath(fname);
    struct stat buf = {0};
    int ret = 0;

    if (gretl_stat(path, &buf) == 0 && buf.st_size > 32) {
        ret = 1;
    }

    g_free(path);

    return ret;
}

static void write_ox_io_file (FILE *fp, const char *ddir)
{
    fputs("gretl_dotdir ()\n{\n", fp);
    fprintf(fp, "  return \"%s\";\n", ddir);
    fputs("}\n\n", fp);

    fputs("gretl_export_nodot (const X, const str)\n{\n", fp);
    fputs("  decl fp = fopen(str, \"w\");\n", fp);
    fputs("  fprint(fp, \"%d \", rows(X), \"%d\", columns(X));\n", fp);
    fputs("  fprint(fp, \"%.15g\", X);\n", fp);
    fputs("  fclose(fp);\n}\n\n", fp);

    fputs("gretl_export (const X, const str)\n{\n", fp);
    fputs("  decl dname = gretl_dotdir();\n", fp);
    fputs("  decl fp = fopen(dname ~ str, \"w\");\n", fp);
    fputs("  fprint(fp, \"%d \", rows(X), \"%d\", columns(X));\n", fp);
    fputs("  fprint(fp, \"%.15g\", X);\n", fp);
    fputs("  fclose(fp);\n}\n\n", fp);

    fputs("gretl_loadmat (const str)\n{\n", fp);
    fputs("  decl dname = gretl_dotdir();\n", fp);
    fputs("  decl X = loadmat(dname ~ str);\n", fp);
    fputs("  return X;\n}\n", fp);
}

static void write_octave_io_file (FILE *fp, const char *ddir)
{
    fputs("# not a 'function file' as such\n1;\n\n", fp);
    fputs("function dotdir = gretl_dotdir()\n", fp);
    fprintf(fp, "  dotdir = \"%s\";\n", ddir);
    fputs("endfunction\n\n", fp);

    fputs("function gretl_export(X, fname, autodot=1)\n", fp);
    fputs("  if (autodot && !is_absolute_filename(fname))\n", fp);
    fputs("    dname = gretl_dotdir();\n", fp);
    fputs("    fd = fopen(strcat(dname, fname), \"w\");\n", fp);
    fputs("  else\n", fp);
    fputs("    fd = fopen(fname, \"w\");\n", fp);
    fputs("  endif\n", fp);
    fputs("  ext = fname(end-3:end);\n", fp);
    fputs("  if (strcmp(ext, '.bin'))\n", fp);
    fputs("    if (iscomplex(X))\n", fp);
    fputs("      fwrite(fd, 'gretl_binar_cmatrix', \"uchar\");\n", fp);
    fputs("      fwrite(fd, 2*rows(X), \"int32\", 0, \"l\");\n", fp);
    fputs("      fwrite(fd, columns(X), \"int32\", 0, \"l\");\n", fp);
    fputs("      A = real(X);\n", fp);
    fputs("      B = imag(X);\n", fp);
    fputs("      N = rows(X) * columns(X);\n", fp);
    fputs("      for i = 1:N\n", fp);
    fputs("        fwrite(fd, A(i), \"double\", 0, \"l\");\n", fp);
    fputs("        fwrite(fd, B(i), \"double\", 0, \"l\");\n", fp);
    fputs("      endfor\n", fp);
    fputs("    else\n", fp);
    fputs("      fwrite(fd, 'gretl_binary_matrix', \"uchar\");\n", fp);
    fputs("      fwrite(fd, rows(X), \"int32\", 0, \"l\");\n", fp);
    fputs("      fwrite(fd, columns(X), \"int32\", 0, \"l\");\n", fp);
    fputs("      fwrite(fd, X, \"double\", 0, \"l\");\n", fp);
    fputs("    endif\n", fp);
    fputs("  else\n", fp);
    fputs("    fprintf(fd, \"%d %d\\n\", size(X));\n", fp);
    fputs("    c = columns(X);\n", fp);
    fputs("    fs = strcat(strrep(sprintf(\"%d \", ones(1, c)), "
          "\"1\", \"%.15g\"), \"\\n\");", fp);
    fputc('\n', fp);
    fputs("    fprintf(fd, fs, X');\n", fp);
    fputs("  endif\n", fp);
    fputs("  fclose(fd);\n", fp);
    fputs("endfunction\n\n", fp);

    fputs("function A = gretl_loadmat(fname, autodot=1)\n", fp);
    fputs("  if (autodot && !is_absolute_filename(fname))\n", fp);
    fputs("    dname = gretl_dotdir();\n", fp);
    fputs("    fd = fopen(strcat(dname, fname), \"r\");\n", fp);
    fputs("  else\n", fp);
    fputs("    fd = fopen(fname, \"r\");\n", fp);
    fputs("  endif\n", fp);
    fputs("  ext = fname(end-3:end);\n", fp);
    fputs("  if (strcmp(ext, '.bin'))\n", fp);
    fputs("    hdr = fread(fd, 19, \"uchar\")';\n", fp);
    fputs("    if sum(hdr == 'gretl_binary_matrix') == 19\n", fp);
    fputs("      cmplx = 0;\n", fp);
    fputs("    elseif sum(hdr == 'gretl_binar_cmatrix') == 19\n", fp);
    fputs("      cmplx = 1;\n", fp);
    fputs("    else\n", fp);
    fputs("      error(\"not a valid gretl binary matrix\");\n", fp);
    fputs("    endif\n", fp);
    fputs("    r = fread(fd, 1, \"int32\", 0, \"l\");\n", fp);
    fputs("    c = fread(fd, 1, \"int32\", 0, \"l\");\n", fp);
    fputs("    A = fread(fd, [r, c], \"double\", 0, \"l\");\n", fp);
    fputs("    if cmplx == 1;\n", fp);
    fputs("      Re = A(1:2:end,:);\n", fp);
    fputs("      Im = A(2:2:end,:);\n", fp);
    fputs("      A = complex(Re, Im);\n", fp);
    fputs("    endif\n", fp);
    fputs("  else\n", fp);
    fputs("    [r,c] = fscanf(fd, \"%d %d\", \"C\");\n", fp);
    fputs("    A = reshape(fscanf(fd, \"%g\", r*c),c,r)';\n", fp);
    fputs("  endif\n", fp);
    fputs("  fclose(fd);\n", fp);
    fputs("endfunction\n\n", fp);
}

static void write_python_io_file (FILE *fp, const char *ddir)
{
    fprintf(fp, "gretl_dotdir = \"%s\"\n\n", ddir);
    fputs("import os\n", fp);

    /* export matrix for reading by gretl */
    fputs("\ndef gretl_export(X, fname, autodot=1):\n", fp);
    fputs("  binwrite = 0\n", fp);
    fputs("  if fname[-4:] == '.bin':\n", fp);
    fputs("    binwrite = 1\n", fp);
    fputs("    from numpy import asmatrix, asarray, iscomplexobj, real, imag\n", fp);
    fputs("    from struct import pack\n", fp);
    fputs("  else:\n", fp);
    fputs("    from numpy import asmatrix, savetxt\n", fp);
    fputs("  M = asmatrix(X)\n", fp);
    fputs("  r, c = M.shape\n", fp);
    fputs("  if autodot and not os.path.isabs(fname):\n", fp);
    fputs("    fname = gretl_dotdir + fname\n", fp);
    fputs("  if binwrite:\n", fp);
    fputs("    from sys import byteorder\n", fp);
    fputs("    cmplx = iscomplexobj(X)\n", fp);
    fputs("    f = open(fname, 'wb')\n", fp);
    fputs("    if cmplx:\n", fp);
    fputs("      f.write(b'gretl_binar_cmatrix')\n", fp);
    fputs("      f.write(pack('<i', 2*r))\n", fp);
    fputs("    else:\n", fp);
    fputs("      f.write(b'gretl_binary_matrix')\n", fp);
    fputs("      f.write(pack('<i', r))\n", fp);
    fputs("    f.write(pack('<i', c))\n", fp);
    fputs("    if cmplx:\n", fp);
    fputs("      for j in range(0, c):\n", fp);
    fputs("        for i in range(0, r):\n", fp);
    fputs("          f.write(pack('<d', real(M[i,j])))\n", fp);
    fputs("          f.write(pack('<d', imag(M[i,j])))\n", fp);
    fputs("    elif byteorder == 'big':\n", fp);
    fputs("      for j in range(0, c):\n", fp);
    fputs("        for i in range(0, r):\n", fp);
    fputs("          f.write(pack('<d', M[i,j]))\n", fp);
    fputs("    else:\n", fp);
    fputs("      A = asarray(X, dtype=float)\n", fp);
    fputs("      f.write(A.tobytes('F'))\n", fp);
    fputs("    f.close()\n", fp);
    fputs("  else:\n", fp);
    fputs("    ghead = repr(r) + ' ' + repr(c)\n", fp);
    fputs("    savetxt(fname, M, header=ghead, comments='')\n", fp);

    /* import matrix from gretl */
    fputs("\ndef gretl_loadmat(fname, autodot=1):\n", fp);
    fputs("  if autodot and not os.path.isabs(fname):\n", fp);
    fputs("    fname = gretl_dotdir + fname\n", fp);
    fputs("  if fname[-4:] == '.bin':\n", fp);
    fputs("    from numpy import fromfile\n", fp);
    fputs("    from struct import unpack\n", fp);
    fputs("    f = open(fname, 'rb')\n", fp);
    fputs("    buf = f.read(19)\n", fp);
    fputs("    if buf == b'gretl_binary_matrix':\n", fp);
    fputs("      cmplx = 0\n", fp);
    fputs("    elif buf == b'gretl_binar_cmatrix':\n", fp);
    fputs("      cmplx = 1\n", fp);
    fputs("    else:\n", fp);
    fputs("      raise ValueError('Not a gretl binary matrix')\n", fp);
    fputs("    r = unpack('<i', f.read(4))[0]\n", fp);
    fputs("    c = unpack('<i', f.read(4))[0]\n", fp);
    fputs("    if cmplx == 1:\n", fp);
    fputs("      r = r // 2\n", fp);
    fputs("      M = fromfile(f, dtype=complex, count=-1)\n", fp);
    fputs("    else:\n", fp);
    fputs("      M = fromfile(f, dtype=float, count=-1)\n", fp);
    fputs("    f.close()\n", fp);
    fputs("    M = M.reshape(r, c, order='F')\n", fp);
    fputs("  else:\n", fp);
    fputs("    from numpy import loadtxt\n", fp);
    fputs("    M = loadtxt(fname, skiprows=1)\n", fp);
    fputs("  return M\n\n", fp);
}

static void write_julia_io_file (FILE *fp, const char *ddir)
{
    fprintf(fp, "gretl_dotdir = \"%s\"\n\n", ddir);
    fputs("using Printf\n", fp);
    fputs("using DelimitedFiles\n\n", fp);
    fputs("function gretl_export(M, fname, autodot=1)\n", fp);
    fputs("  if VERSION < v\"1.0\"\n", fp);
    fputs("    error(\"gretl_export requires Julia >= 1.0\")\n", fp);
    fputs("  end\n", fp);
    fputs("  r,c = size(M)\n", fp);
    fputs("  if autodot != 0 && !isabspath(fname)\n", fp);
    fputs("    fname = gretl_dotdir * fname\n", fp);
    fputs("  end\n", fp);
    fputs("  f = open(fname, \"w\")\n", fp);
    fputs("  n = lastindex(fname)\n", fp);
    fputs("  if fname[n-3:n] == \".bin\"\n", fp);
    fputs("    # binary mode\n", fp);
    fputs("    local cmplx::Bool\n", fp);
    fputs("    if typeof(M) == Array{Complex{Float64},2}\n", fp);
    fputs("      cmplx = 1\n", fp);
    fputs("      rr = 2*r\n", fp);
    fputs("      write(f, b\"gretl_binar_cmatrix\")\n", fp);
    fputs("    else\n", fp);
    fputs("      cmplx = 0\n", fp);
    fputs("      rr = r\n", fp);
    fputs("      write(f, b\"gretl_binary_matrix\")\n", fp);
    fputs("    end\n", fp);
    fputs("    write(f, htol(Int32(rr)))\n", fp);
    fputs("    write(f, htol(Int32(c)))\n", fp);
    fputs("    for j=1:c\n", fp);
    fputs("      for i=1:r\n", fp);
    fputs("        if cmplx\n", fp);
    fputs("          write(f, htol(Float64(real(M[i,j]))))\n", fp);
    fputs("          write(f, htol(Float64(imag(M[i,j]))))\n", fp);
    fputs("        else\n", fp);
    fputs("          write(f, htol(Float64(M[i,j])))\n", fp);
    fputs("        end\n", fp);
    fputs("      end\n", fp);
    fputs("    end\n", fp);
    fputs("  else\n", fp);
    fputs("    # text mode\n", fp);
    fputs("    @printf(f, \"%d\\t%d\\n\", r, c)\n", fp);
    fputs("    for i = 1:r\n", fp);
    fputs("      for j = 1:c\n", fp);
    fputs("        @printf(f, \"%.18e \", M[i,j])\n", fp);
    fputs("      end\n", fp);
    fputs("      @printf(f, \"\\n\")\n", fp);
    fputs("    end\n", fp);
    fputs("  end\n", fp);
    fputs("  close(f)\n", fp);
    fputs("end\n\n", fp);

    fputs("function gretl_loadmat(fname, autodot=1)\n", fp);
    fputs("  if VERSION < v\"1.0\"\n", fp);
    fputs("    error(\"gretl_loadmat requires Julia >= 1.0\")\n", fp);
    fputs("  end\n", fp);
    fputs("  if autodot != 0 && !isabspath(fname)\n", fp);
    fputs("    fname = gretl_dotdir * fname\n", fp);
    fputs("  end\n", fp);
    fputs("  n = lastindex(fname)\n", fp);
    fputs("  if fname[n-3:n] == \".bin\"\n", fp);
    fputs("    # binary mode\n", fp);
    fputs("    local cmplx::Bool\n", fp);
    fputs("    f = open(fname, \"r\")\n", fp);
    fputs("    hdr = read(f, 19)\n", fp);
    fputs("    if hdr == b\"gretl_binary_matrix\"\n", fp);
    fputs("      cmplx = 0\n", fp);
    fputs("    elseif hdr == b\"gretl_binar_cmatrix\"\n", fp);
    fputs("      cmplx = 1\n", fp);
    fputs("    else\n", fp);
    fputs("      error(\"Not a gretl binary matrix\")\n", fp);
    fputs("    end\n", fp);
    fputs("    r = ltoh(read(f, Int32))\n", fp);
    fputs("    c = ltoh(read(f, Int32))\n", fp);
    fputs("    if cmplx\n", fp);
    fputs("      r::Int64 = floor(r/2)\n", fp);
    fputs("      M = Array{Complex{Float64}, 2}(undef, r, c)\n", fp);
    fputs("      for j=1:c\n", fp);
    fputs("        for i=1:r\n", fp);
    fputs("          x = ltoh(read(f, Float64))\n", fp);
    fputs("          y = ltoh(read(f, Float64))\n", fp);
    fputs("          M[i,j] = x + y*im\n", fp);
    fputs("        end\n", fp);
    fputs("      end\n", fp);
    fputs("    else\n", fp);
    fputs("      M = Array{Float64, 2}(undef, r, c)\n", fp);
    fputs("      for j=1:c\n", fp);
    fputs("        for i=1:r\n", fp);
    fputs("          M[i,j] = ltoh(read(f, Float64))\n", fp);
    fputs("        end\n", fp);
    fputs("      end\n", fp);
    fputs("    end\n", fp);
    fputs("    close(f)\n", fp);
    fputs("  else\n", fp);
    fputs("    # text mode\n", fp);
    fputs("    M = readdlm(fname, skipstart=1)\n", fp);
    fputs("  end\n", fp);
    fputs("  M\n", fp);
    fputs("end\n\n", fp);
}

static void write_stata_io_file (FILE *fp, const char *ddir)
{
    fputs("program define gretl_export\n", fp);
    /* not sure about req'd version, but see mat2txt.ado */
    fputs("version 8.2\n", fp);
    fputs("local matrix `1'\n", fp);
    fputs("local fname `2'\n", fp);
    fputs("tempname myfile\n", fp);
    fprintf(fp, "file open `myfile' using \"%s`fname'\", "
            "write text replace\n", ddir);
    fputs("local nrows = rowsof(`matrix')\n", fp);
    fputs("local ncols = colsof(`matrix')\n", fp);
    fputs("file write `myfile' %8.0g (`nrows') %8.0g (`ncols') _n\n", fp);
    fputs("forvalues r=1/`nrows' {\n", fp);
    fputs("  forvalues c=1/`ncols' {\n", fp);
    fputs("    file write `myfile' %15.0e (`matrix'[`r',`c']) _n\n", fp);
    fputs("  }\n", fp);
    fputs("}\n", fp);
    fputs("file close `myfile'\n", fp);
    fputs("end\n", fp);
}

static int ensure_foreign_io_file (int lang, const char *fname)
{
    FILE *fp = NULL;
    int err = 0;

    if (fname != NULL) {
        fp = gretl_fopen(fname, "wb");
    } else {
        const char *iofile = get_io_filename(lang);

        if (iofile == NULL) {
            return E_DATA;
        }
        if (dotfile_exists(iofile)) {
            /* already present */
            return 0;
        } else {
            fp = write_open_dotfile(iofile);
        }
    }

    if (fp == NULL) {
        err = E_FOPEN;
    } else {
        const char *ddir = get_export_dotdir();

        if (lang == LANG_PYTHON) {
            write_python_io_file(fp, ddir);
        } else if (lang == LANG_OCTAVE) {
            write_octave_io_file(fp, ddir);
        } else if (lang == LANG_JULIA) {
            write_julia_io_file(fp, ddir);
        } else if (lang == LANG_OX) {
            write_ox_io_file(fp, ddir);
        } else if (lang == LANG_STATA) {
            write_stata_io_file(fp, ddir);
        } else if (lang == LANG_R) {
            /* used only under --io-funcs */
            write_R_io_file(fp, ddir);
        }

        fclose(fp);
    }

    return err;
}

static void add_gretl_include (int lang, gretlopt opt, FILE *fp)
{
    const char *ddir = get_export_dotdir();

    if (lang == LANG_PYTHON) {
        fputs("from gretl_io import gretl_dotdir, gretl_loadmat, "
              "gretl_export\n", fp);
        return;
    }

#ifdef G_OS_WIN32
    if (lang == LANG_STATA) {
        /* or leave path with backslashes? */
        fprintf(fp, "quietly adopath + \"%s\"\n", ddir);
        return;
    }
#endif

    if (lang == LANG_OX) {
        if (strchr(ddir, ' ')) {
            fprintf(fp, "#include \"%sgretl_io.ox\"\n", ddir);
        } else {
            fprintf(fp, "#include <%sgretl_io.ox>\n", ddir);
        }
    } else if (lang == LANG_OCTAVE) {
        fprintf(fp, "source(\"%sgretl_io.m\")\n", ddir);
    } else if (lang == LANG_JULIA) {
        fprintf(fp, "include(\"%sgretl_io.jl\")\n", ddir);
    } else if (lang == LANG_STATA) {
        if (opt & OPT_Q) {
            fputs("set output error\n", fp);
        }
        fprintf(fp, "quietly adopath + \"%s\"\n", ddir);
    }
}

static int get_foreign_indent (void)
{
    const char *s;
    int i, n, ret = 100;

    for (i=0; i<foreign_n_lines; i++) {
        n = 0;
        s = foreign_lines[i];
        while (*s == ' ' || *s == '\t') {
            n++;
            s++;
        }
        if (n < ret) {
            ret = n;
        }
    }

    return ret;
}

static void put_foreign_lines (FILE *fp)
{
    int i, n = get_foreign_indent();

    for (i=0; i<foreign_n_lines; i++) {
        fprintf(fp, "%s\n", foreign_lines[i] + n);
        if (foreign_lang == LANG_OX) {
            if (strstr(foreign_lines[i], "oxstd.h")) {
                add_gretl_include(LANG_OX, 0, fp);
            }
        } else if (foreign_lang == LANG_OCTAVE) {
            if (strstr(foreign_lines[i], "dynare ") &&
                !strstr(foreign_lines[i], "noclearall")) {
                add_gretl_include(LANG_OCTAVE, 0, fp);
            }
        }
    }
}

static void put_foreign_buffer (const char *buf, FILE *fp)
{
    char line[1024];

    bufgets_init(buf);

    while (bufgets(line, sizeof line, buf)) {
        fputs(line, fp);
        if (foreign_lang == LANG_OX) {
            if (strstr(line, "oxstd.h")) {
                add_gretl_include(LANG_OX, 0, fp);
            }
        } else if (foreign_lang == LANG_OCTAVE) {
            if (strstr(line, "dynare ") &&
                !strstr(line, "noclearall")) {
                add_gretl_include(LANG_OCTAVE, 0, fp);
            }
        }
    }

    bufgets_finalize(buf);
}

static int no_data_check (const DATASET *dset, int ci)
{
    if (dset == NULL || dset->v == 0) {
	return E_NODATA;
    } else if (ci != MPI && dset->n == 0) {
        return E_NODATA;
    } else {
        return 0;
    }
}

static int *get_send_data_list (int ci, const DATASET *dset,
                                int *err)
{
    const char *dname = get_optval_string(ci, OPT_D);
    static int list1[2] = {1, 0};
    int *list = NULL;

    if (dname != NULL) {
        list = get_list_by_name(dname);
        if (list != NULL) {
            int i;

            for (i=1; i<=list[0] && !*err; i++) {
                if (list[i] < 0 || list[i] >= dset->v) {
                    *err = E_DATA;
                }
            }
        } else {
            int vi = current_series_index(dset, dname);

            if (vi >= 0) {
                list1[1] = vi;
                list = list1;
            } else {
                gretl_errmsg_sprintf(_("'%s': no such list"), dname);
                *err = E_DATA;
            }
        }
    }

    return list;
}

#ifdef HAVE_MPI

static int mpi_send_data_setup (const DATASET *dset,
				gretlopt opt,
				FILE *fp)
{
    gretlopt gdt_opt = OPT_NONE;
    int zlist[2] = {1, 0};
    int *list = NULL;
    size_t datasize;
    int nvars;
    gchar *fname;
    int err;

    err = no_data_check(dset, MPI);
    if (err) {
        return err;
    }

    if (opt & OPT_R) {
	/* --full-range */
	gdt_opt = OPT_F;
    }

    if (opt & OPT_M) {
	/* --send-metadata (only) */
	nvars = 0; /* nvars is net of const */
	list = zlist;
    } else {
	list = get_send_data_list(MPI, dset, &err);
	if (err) {
	    return err;
	} else if (list != NULL) {
	    nvars = list[0];
	} else {
	    nvars = dset->v;
	}
    }

    datasize = dset->n * nvars;

    if (!(opt & OPT_R) && datasize > 10000) {
        /* write "big" data as binary? */
        fname = gretl_make_dotpath("mpi-data.gdtb");
    } else {
        fname = gretl_make_dotpath("mpi-data.gdt");
    }

    err = gretl_write_gdt(fname, list, dset, gdt_opt, 0);

    if (!err) {
        /* here we're writing into the file to be run
           by gretlmpi */
        fprintf(fp, "open \"%s\" --quiet\n", fname);
    }

    g_free(fname);

    return err;
}

static int mpi_send_funcs_setup (FILE *fp)
{
    gchar *fname;
    int err;

    fname = gretl_make_dotpath("mpi-funcs-tmp.xml");
    err = write_loaded_functions_file(fname, 1);

    if (!err) {
        fprintf(fp, "include \"%s\"\n", fname);
    }

    g_free(fname);

    return err;
}

static int write_gretl_mpi_script (gretlopt opt, DATASET *dset)
{
    const gchar *fname = get_mpi_scriptname();
    FILE *fp = NULL;
    int err;

    err = incompatible_options(opt, OPT_D | OPT_M);
    if (err) {
	return err;
    }

    fp = gretl_fopen(fname, "w");
    if (fp == NULL) {
        return E_FOPEN;
    }

    if (opt & (OPT_D | OPT_M)) {
        /* honor the --send-data option */
        err = mpi_send_data_setup(dset, opt, fp);
    }

    if (opt & OPT_F) {
        /* honor the --send-functions option */
        if (n_user_functions() > 0) {
            err = mpi_send_funcs_setup(fp);
        }
    }

#ifdef _OPENMP
    if (!err) {
        if (opt & OPT_T) {
            /* respect the --omp-threads option */
            int nt = get_optval_int(MPI, OPT_T, &err);

            if (nt == -1) {
                ; /* unlimited/auto */
            } else {
                if (!err && (nt <= 0 || nt > 9999999)) {
                    err = E_DATA;
                }
                if (!err) {
                    fprintf(fp, "set omp_num_threads %d\n", nt);
                }
            }
        } else {
            /* by default, don't use OMP threading */
            fputs("set omp_num_threads 1\n", fp);
        }
    }
#endif

    if (!err) {
        /* put out the stored 'foreign' lines */
        put_foreign_lines(fp);
#ifdef MPI_REALTIME
        /* plus an easily recognized trailer */
        fputs("flush\n", fp);
        fputs("mpibarrier()\n", fp);
        fputs("if $mpirank == 0\n", fp);
        fputs("  print \"__GRETLMPI_EXIT__\"\n", fp);
        fputs("endif\n", fp);
#endif
    }

    fclose(fp);

    return err;
}

#endif /* HAVE_MPI */

static int write_data_for_stata (const DATASET *dset,
                                 FILE *fp)
{
    int *list = NULL;
    char save_na[8];
    int err;

    err = no_data_check(dset, FOREIGN);
    if (err) {
        return err;
    }

    list = get_send_data_list(FOREIGN, dset, &err);

    if (!err) {
        gchar *sdata;

        *save_na = '\0';
        strncat(save_na, get_csv_na_write_string(), 7);
        set_csv_na_write_string(".");
        sdata = gretl_make_dotpath("stata.csv");
        err = write_data(sdata, list, dset, OPT_C, NULL);
        set_csv_na_write_string(save_na);
        g_free(sdata);
    }

    if (err) {
        gretl_errmsg_sprintf("write_data_for_stata: failed with err = %d\n", err);
    } else {
        fputs("* load data from gretl\n", fp);
        fprintf(fp, "insheet using \"%sstata.csv\"\n", get_export_dotdir());
    }

    return err;
}

/* write out current dataset as an octave matrix, and, if this succeeds,
   write appropriate octave commands to @fp to source the data
*/

static int write_data_for_octave (const DATASET *dset,
                                  FILE *fp)
{
    int *list = NULL;
    int err;

    err = no_data_check(dset, FOREIGN);
    if (err) {
        return err;
    }

    list = get_send_data_list(FOREIGN, dset, &err);

    if (!err) {
        gchar *mdata = gretl_make_dotpath("mdata.tmp");

        err = write_data(mdata, list, dset, OPT_M, NULL);
        g_free(mdata);
    }

    if (err) {
        gretl_errmsg_sprintf("write_data_for_octave: failed with err = %d\n", err);
    } else {
        fputs("% load data from gretl\n", fp);
        fprintf(fp, "load '%smdata.tmp'\n", get_export_dotdir());
    }

    return err;
}

static int put_dynare_script (const char *buf, FILE *fp)
{
    const char *wdir = gretl_workdir();
    gchar *modpath;
    FILE *fd;
    int err = 0;

    modpath = g_build_filename(wdir, "gretltmp.mod", NULL);
    fd = gretl_fopen(modpath, "w");

    if (fd == NULL) {
        err = E_FOPEN;
    } else {
        fputs(buf, fd);
        fclose(fd);
	fprintf(fp, "cd \"%s\"\n", wdir);
        fputs("dynare gretltmp.mod noclearall\n", fp);
    }

    g_free(modpath);

    return err;
}

/**
 * write_gretl_foreign_script:
 * @buf: text buffer containing foreign code: Ox, Octave, Stata,
 * Python or Julia; R is handled separately.
 * @lang: language identifier.
 * @opt: should contain %OPT_G for use from GUI.
 * @dset: pointer to dataset (or %NULL if no data are to be passed).
 * @pfname: location to receive name of file written, or %NULL.
 *
 * Writes the content of @buf into a file in the gretl user's
 * "dotdir".
 *
 * Returns: 0 on success, non-zero on error.
 */

int write_gretl_foreign_script (const char *buf, int lang,
                                gretlopt opt, const DATASET *dset,
                                const char **pfname)
{
    const gchar *fname = get_foreign_scriptpath(lang);
    FILE *fp = gretl_fopen(fname, "w");
    int err = 0;

    ensure_foreign_io_file(lang, NULL);

    if (fp == NULL) {
        err = E_FOPEN;
    } else {
        /* source the I-O functions? */
        if (lang != LANG_OX) {
            /* this comes later for Ox */
            add_gretl_include(lang, opt, fp);
        }
        if (dset != NULL && (opt & OPT_D)) {
            /* --send-data */
            if (lang == LANG_OCTAVE) {
                err = write_data_for_octave(dset, fp);
            } else if (lang == LANG_STATA) {
                err = write_data_for_stata(dset, fp);
            }
        }
        if (!err && buf != NULL) {
            /* pass on the material supplied in the @buf argument */
            if (lang == LANG_OCTAVE && (opt & OPT_Y)) {
                /* handle a dynare .mod file */
                err = put_dynare_script(buf, fp);
            } else {
                /* regular script */
                put_foreign_buffer(buf, fp);
            }
        } else if (!err) {
            /* put out the stored 'foreign' lines */
            put_foreign_lines(fp);
        }
        if (!err && pfname != NULL) {
            *pfname = fname;
        }
    }

    if (fp != NULL) {
        fclose(fp);
    }

    return err;
}

static gretl_matrix *make_coded_vec (int *list,
                                     const DATASET *dset)
{
    gretl_matrix *coded = NULL;
    int free_list = 0;
    int i, nc = 0;

    if (list == NULL) {
        list = full_var_list(dset, NULL);
        free_list = 1;
    }

    if (list != NULL) {
        for (i=1; i<=list[0]; i++) {
            if (series_is_coded(dset, list[i])) {
                nc++;
            }
        }
    }

    if (nc > 0) {
        coded = gretl_matrix_alloc(1, nc);
        if (coded != NULL) {
            int j = 0;

            for (i=1; i<=list[0]; i++) {
                if (series_is_coded(dset, list[i])) {
                    coded->val[j++] = i;
                }
            }
        }
    }

    if (free_list) {
        free(list);
    }

    return coded;
}

/* write out current dataset in R format, and, if this succeeds,
   write appropriate R commands to @fp to source the data
*/

static int write_data_for_R (const DATASET *dset,
                             gretlopt opt,
                             FILE *fp)
{
    gretl_matrix *coded = NULL;
    int *list = NULL;
    int ts, err;

    err = no_data_check(dset, FOREIGN);
    if (err) {
        return err;
    }

    /* FIXME: can R's "ts" handle daily data, weekly data, etc.? */
    ts = annual_data(dset) || quarterly_or_monthly(dset);

    list = get_send_data_list(FOREIGN, dset, &err);

    if (!err) {
        gchar *Rdata = gretl_make_dotpath("Rdata.tmp");

        coded = make_coded_vec(list, dset);
        err = write_data(Rdata, list, dset, OPT_R, NULL);
        g_free(Rdata);
    }

    if (err) {
        gretl_errmsg_sprintf("write_data_for_R: failed with err = %d\n", err);
        gretl_matrix_free(coded);
        return err;
    }

    if (coded != NULL) {
        gchar *tmp = gretl_make_dotpath("Rcoded.mat");
        int write_err;

        /* ensure we don't load a stale file */
        gretl_remove(tmp);
        g_free(tmp);
        write_err = gretl_matrix_write_to_file(coded, "Rcoded.mat", 1);
        if (write_err) {
           gretl_matrix_free(coded);
           coded = NULL;
        }
    }

    fputs("# load data from gretl\n", fp);
    fprintf(fp, "gretldata <- read.table(\"%sRdata.tmp\", header=TRUE)\n",
            get_export_dotdir());

    if (ts) {
        char *p, datestr[OBSLEN];
        int subper = 1;

        ntolabel(datestr, dset->t1, dset);
        p = strchr(datestr, ':');
        if (p != NULL) {
            subper = atoi(p + 1);
        }

        if (opt & OPT_F) {
            /* treat as data frame (but set columns as "ts") */
            fputs("if (length(class(gretldata)) > 1) {m <- ncol(x)} else {m <- 1}\n", fp);
            fputs("for (i in 1:m) {\n", fp);
            fprintf(fp, "  gretldata[,i] <- ts(gretldata[,i], start=c(%d, %d), "
                  "frequency=%d)\n", atoi(datestr), subper, dset->pd);
            fputs("}\n", fp);
            fputs("attach(gretldata)\n", fp);
        } else {
            /* convert to "mts" (multiple time series object) */
            fprintf(fp, "gretldata <- ts(gretldata, start=c(%d, %d), frequency = %d)\n",
                    atoi(datestr), subper, dset->pd);
        }
    } else {
        fputs("attach(gretldata)\n", fp);
    }

    if (coded != NULL) {
        fputs("Coded <- gretl.loadmat(\"Rcoded.mat\")\n", fp);
        fputs("for (i in Coded) {gretldata[,i] <- as.factor(gretldata[,i])}\n", fp);
    }

    gretl_matrix_free(coded);

    if (opt & OPT_I) {
        /* let the (interactive) user see that this worked */
        if (ts) {
            fputs("gretlmsg <- \"current data loaded as ts object \\\"gretldata\\\"\\n\"\n", fp);
        } else {
            fputs("gretlmsg <- \"current data loaded as data frame \\\"gretldata\\\"\\n\"\n", fp);
        }
        fputs("cat(gretlmsg)\n", fp);
    }

    return err;
}

/* define an R function for passing data back to gretl */

static void write_R_io_file (FILE *fp, const char *ddir)
{
    const char *export_body =
        "  objname <- as.character(substitute(x))\n"
        "  if (missing(sx)) {\n"
        "    sx <- objname\n"
        "  }\n"
        "  if (is.character(x) && length(x) == 1) {\n"
        "    fname <- paste(prefix, sx, \".txt\", sep=\"\")\n"
        "    cat(x, file=fname, sep='\n')\n"
        "    gretlmsg <- paste(\"wrote txt data\", fname, \"\\n\")\n"
        "  } else if (is.ts(x)) {\n"
        "    fname <- paste(prefix, sx, \".csv\", sep=\"\")\n"
        "    dfx <- data.frame(x)\n"
        "    if (ncol(dfx) == 1) {\n"
        "      colnames(dfx) <- sx;\n"
        "    }\n"
        "    write.csv(dfx, file=fname, row.names=F)\n"
        "    gretlmsg <- paste(\"wrote CSV data\", fname, \"\\n\")\n"
        "  } else if (is.data.frame(x)) {\n"
        "    fname <- paste(prefix, sx, \".csv\", sep=\"\")\n"
        "    write.csv(x, file=fname, row.names=F)\n"
        "    gretlmsg <- paste(\"wrote CSV data\", fname, \"\\n\")\n"
        "  } else if (is.matrix(x)) {\n"
        "    len = nchar(sx)\n"
        "    if (substr(sx, len-3, len) == \".bin\") {\n"
        "      fname <- paste(prefix, sx, sep=\"\")\n"
        "      mout = file(fname, \"wb\")\n"
        "      b <- charToRaw(\"gretl_binary_matrix\")\n"
        "      writeBin(as.integer(b), mout, size=1)\n"
        "      writeBin(nrow(x), mout, size=4)\n"
        "      writeBin(ncol(x), mout, size=4)\n"
        "      writeBin(c(x), mout)\n"
        "      close(mout)\n"
        "    } else {\n"
        "      fname <- paste(prefix, sx, \".mat\", sep=\"\")\n"
        "      write(dim(x), fname)\n"
        "      write(format(t(x), digits=15), file=fname, ncolumns=ncol(x), append=TRUE)\n"
        "    }\n"
        "    gretlmsg <- paste(\"wrote matrix\", fname, \"\\n\")\n"
        "  } else {\n"
        "    gretlmsg <- paste(\"gretl.export: don't know how to write object\", objname, "
        " \"(try as.matrix?)\\n\")\n"
        "  }\n"
        "  if (!quiet) {\n"
        "    cat(gretlmsg)\n"
        "  }\n"
        "}\n";
    const char *loadmat_body =
        "  len <- nchar(mname)\n"
        "  binfile <- substr(mname, len-3, len) == \".bin\"\n"
        "  if (binfile) {\n"
        "    mfile = file(fname, \"rb\")\n"
        "    hdr = readBin(mfile, integer(), size=1, n=19)\n"
        "    s <- rawToChar(as.raw(hdr))\n"
        "    if (s != \"gretl_binary_matrix\") {\n"
        "      stop('Not a gretl binary matrix')\n"
        "    } else {\n"
        "      nr <- readBin(mfile, integer(), size=4, endian = \"little\")\n"
        "      nc <- readBin(mfile, integer(), size=4, endian = \"little\")\n"
        "      matvals <- readBin(mfile, double(), nr*nc, endian = \"little\")\n"
        "      m <- matrix(matvals, nr, nc)\n"
        "    }\n"
        "  } else {\n"
        "    m <- as.matrix(read.table(fname, skip=1))\n"
        "  }\n"
        "  return(m)\n"
        "}\n";

    fprintf(fp, "gretl.dotdir <- \"%s\"\n", ddir);

#if 0
    fputs("is_abspath <- function(path) {\n", fp);
    fputs("grepl(\"^(/|[A-Za-z]:|\\\\\\\\|~)\", path)\n");
    fputs("}\n", fp);
#endif

    fputs("gretl.export <- function(x, sx, quiet=0) {\n", fp);
    fprintf(fp, "  prefix <- \"%s\"\n", ddir);
    fputs(export_body, fp);

    fputs("gretl.loadmat <- function(mname, dotpath=1) {\n", fp);
    fputs("  if (dotpath) {\n", fp);
    fprintf(fp, "    prefix <- \"%s\"\n", ddir);
    fputs("    fname <- paste(prefix, mname, sep=\"\")\n", fp);
    fputs("  } else {\n", fp);
    fputs("    fname <- mname\n", fp);
    fputs("  }\n", fp);
    fputs(loadmat_body, fp);
}

static gchar *get_checkdeps_buffer (const char *deps)
{
    const char *checkdeps_body =
        "checkdeps <- function(deps) {\n"
        "  terms <- unlist(strsplit(deps, \" \"))\n"
        "  n <- length(terms) / 2\n"
        "  n_ok <- n\n"
        "  j <- 1\n"
        "  for (i in 1:n) {\n"
        "    obj <- terms[j]\n"
        "    minver <- terms[j+1]\n"
        "    cat(sprintf(\"check for %s >= %s: \", obj, minver))\n"
        "    if (obj == \"R\") {\n"
        "      rv <- getRversion()\n"
        "        if (getRversion() >= minver) {\n"
        "          cat(sprintf(\"%s OK\\n\", rv))\n"
        "        } else {\n"
        "          cat(sprintf(\"%s too old\\n\", rv))\n"
        "          n_ok <- n_ok - 1\n"
        "        }\n"
        "    } else {\n"
        "      if (!require(obj, character.only=TRUE, quietly=TRUE)) {\n"
        "        cat(\"not found\\n\")\n"
        "        n_ok <- n_ok - 1\n"
        "      } else {\n"
        "        pv <- packageVersion(obj)\n"
        "        if (packageVersion(obj) >= minver) {\n"
        "          cat(sprintf(\"%s OK\\n\", pv))\n"
        "        } else {\n"
        "          cat(sprintf(\"%s too old\\n\", pv))\n"
        "          n_ok <- n_ok - 1\n"
        "        }\n"
        "      }\n"
        "    }\n"
        "    j <- j + 2\n"
        "  }\n"
        "  if (n_ok < n) {\n"
        "    cat(\"R requirements are not met\\n\")\n"
        "  }\n"
        "  err <- n_ok < n\n"
        "}\n";
    const char *suppress = "suppressPackageStartupMessages";
    GString *gs = g_string_new(checkdeps_body);

    g_string_append_printf(gs, "%s(checkdeps(\"%s\"))\n", suppress, deps);

    return g_string_free(gs, FALSE);
}

/* basic content which can either go into gretl.Rsetup or into
   Rsrc for sourcing */

static void put_R_setup_content (FILE *fp)
{
    fputs("vnum <- as.double(R.version$major) + (as.double(R.version$minor) / 10.0)\n",
          fp);
    fputs("if (vnum > 2.41) library(utils)\n", fp);
    fputs("library(stats)\n", fp);
    fputs("if (vnum <= 1.89) library(ts)\n", fp);
    write_R_io_file(fp, get_export_dotdir());
}

/* Set up a gretl-specific R profile, and put notice of its existence
   into the environment. Used when exec'ing the R binary (only) */

static int write_gretl_R_profile (gretlopt opt)
{
    FILE *f1, *f2;
    int err = 0;

#if FDEBUG
    fprintf(stderr, "writing R profile: interactive = %d\n",
            (opt & OPT_I)? 1 : 0);
#endif

    /* On Windows we'll not use the environment-variable
       mechanism unless we're in interactive (async) mode
    */

#ifdef G_OS_WIN32
    if (opt & OPT_I) {
        err = gretl_setenv("R_PROFILE", gretl_Rprofile);
        if (err) {
            return err;
        }
    }
#else
    err = gretl_setenv("R_PROFILE", gretl_Rprofile);
    if (err) {
        return err;
    }
#endif

    f1 = gretl_fopen(gretl_Rprofile, "w");
    f2 = gretl_fopen(gretl_Rsetup, "w");

    if (f1 == NULL || f2 == NULL) {
        err = E_FOPEN;
    } else {
        fprintf(f1, "source(\"%s\")\n", gretl_Rsetup);
        fprintf(f1, "source(\"%s\", %s = TRUE)\n",
                gretl_Rsrc, (opt & OPT_V)? "echo" : "print.eval");
        fclose(f1);
        put_R_setup_content(f2);
        fclose(f2);
    }

#if FDEBUG
    fprintf(stderr, "writing R profile: returning %d\n", err);
#endif

    return err;
}

/* Write an R command file to be sourced by R.  @buf may contain R
   commands assembled via the GUI; if it is NULL the current "foreign"
   block (if any) is used as input.

   @opt may contain the following:

   OPT_I: indicates that we're in the context of an interactive R
   session.

   OPT_D: indicates that the current gretl dataset should be sent
   to R.

   OPT_G: we're being called via the gretl GUI.

   OPT_L: indicates that the source file is intended for use
   via the R shared library.
*/

static int write_R_source_file (const char *buf,
                                const DATASET *dset,
                                gretlopt opt)
{
    FILE *fp = gretl_fopen(gretl_Rsrc, "w");
    int err = 0;

#if FDEBUG
    fprintf(stderr, "write R source file: interactive %d, library %d\n",
            (opt & OPT_I)? 1 : 0, (opt & OPT_L)? 1 : 0);
#endif

    if (fp == NULL) {
        err = E_FOPEN;
    } else {
        int sunk = 0;

#ifdef G_OS_WIN32
        if (!(opt & OPT_I) && !(opt & OPT_L)) {
            /* Windows, non-interactive, not using Rlib */
            fprintf(fp, "sink(\"%s\", type=\"output\")\n", gretl_Rout);
            fputs("if (sink.number(\"message\") == 2) {\n", fp);
            fprintf(fp, "  errout <- file(\"%s\", open=\"wt\")\n", gretl_Rmsg);
            fputs("  sink(errout, type=\"message\")\n", fp);
            fputs("}\n", fp);
            sunk = 1;
        }
#endif

        if (opt & OPT_L) {
            /* we're using the R shared library */
            static int setup_done;

            if (!setup_done) {
#if FDEBUG
                fprintf(stderr, "Rlib: writing 'setup' material\n");
#endif
                put_R_setup_content(fp);
                setup_done = 1;
            }
            fprintf(fp, "sink(\"%s\", type=\"output\")\n", gretl_Rout);
            if (!(opt & OPT_I)) {
                fputs("if (sink.number(\"message\") == 2) {\n", fp);
                fprintf(fp, "  errout <- file(\"%s\", open=\"wt\")\n", gretl_Rmsg);
                fputs("  sink(errout, type=\"message\")\n", fp);
                fputs("}\n", fp);
            }
            sunk = 1;
        }

        if (opt & OPT_D) {
            /* --send-data */
            err = write_data_for_R(dset, opt, fp);
            if (err) {
                fclose(fp);
                return err;
            }
        }

        if (buf != NULL) {
            /* pass on the script supplied in @buf */
            fputs("# load script from gretl\n", fp);
            fputs(buf, fp);
        } else if (!(opt & OPT_G)) {
            /* non-GUI */
            put_foreign_lines(fp);
        }

        if (sunk) {
            fputs("sink()\n", fp);
        }

        fclose(fp);
    }

#if FDEBUG
    fprintf(stderr, "write R source file: returning %d\n", err);
#endif

    return err;
}

/* Write files to be read by R: profile to be read on startup, and
   command source file.  This is called when we're exec'ing the R
   binary.  OPT_G in @opt indicates that this function is being called
   from the GUI program; @buf may contain R commands taken from a GUI
   window, or may be NULL.
*/

int write_gretl_R_files (const char *buf,
                         const DATASET *dset,
                         gretlopt opt)
{
    int err = 0;

#if FDEBUG
    fprintf(stderr, "write_gretl_R_files: starting\n");
#endif

    make_gretl_R_names();

    /* write a temporary R profile so R knows what to do */
    err = write_gretl_R_profile(opt);
    if (err) {
        fprintf(stderr, "error writing gretl.Rprofile\n");
    }

    if (!err) {
        /* write commands and/or data to file, to be sourced in R */
        err = write_R_source_file(buf, dset, opt);
        if (err) {
            fprintf(stderr, "error writing gretl's Rsrc\n");
        }
    }

#if FDEBUG
    fprintf(stderr, "write_gretl_R_files: returning %d\n", err);
#endif

    return err;
}

void delete_gretl_R_files (void)
{
#if FDEBUG
    fprintf(stderr, "deleting gretl R files...\n");
#endif

    if (gretl_Rprofile != NULL) {
        gretl_remove(gretl_Rprofile);
    }
    if (gretl_Rsrc != NULL) {
        gretl_remove(gretl_Rsrc);
    }
}

#ifdef HAVE_MPI

static void delete_mpi_script (void)
{
    if (gretl_mpi_script != NULL) {
        gretl_remove(gretl_mpi_script);
    }
}

#endif

/* The following code block is used if we're implementing
   gretl's R support by dlopening the R shared library
   (as opposed to executing the R binary).
*/

#ifdef USE_RLIB

static void *Rhandle;  /* handle to the R library */
static int Rlib_err;   /* initialization error record */
static int Rinit;      /* are we initialized or not? */

static SEXP current_arg;
static SEXP current_call;

/* pointers to, and renamed versions of, the R global variables
   we'll need */

SEXP *PR_GlobalEnv;
SEXP *PR_NilValue;
SEXP *PR_UnboundValue;
SEXP *PR_NamesSymbol;
SEXP *PR_DimNamesSymbol;

SEXP VR_GlobalEnv;
SEXP VR_NilValue;
SEXP VR_UnboundValue;
SEXP VR_NamesSymbol;
SEXP VR_DimNamesSymbol;

/* renamed, pointerized versions of the R functions we need */

static double *(*R_REAL) (SEXP);
static int *(*R_INT) (SEXP);
static const char *(*R_STRING) (SEXP);
static SEXP *(*R_STRING_PTR) (SEXP);
static SEXP (*R_STRING_ELT) (SEXP, int);

static SEXP (*R_CDR) (SEXP);
static SEXP (*R_allocList) (int);
static SEXP (*R_allocMatrix) (SEXPTYPE, int, int);
static SEXP (*R_allocVector) (SEXPTYPE, R_len_t);
static SEXP (*R_findFun) (SEXP, SEXP);
static SEXP (*R_findVar) (SEXP, SEXP);
static SEXP (*R_SETCAR) (SEXP, SEXP);
static SEXP (*R_protect) (SEXP);
static SEXP (*R_ScalarReal) (double);
static SEXP (*R_catch) (SEXP, SEXP, int *);
static SEXP (*R_install) (const char *);
static SEXP (*R_mkString) (const char *);
static SEXP (*R_mkChar) (const char *);
static SEXP (*R_mkNamed) (SEXPTYPE, const char **);
static SEXP (*R_getAttrib) (SEXP, SEXP);
static SEXP (*R_setAttrib) (SEXP, SEXP, SEXP);

static Rboolean (*R_isMatrix) (SEXP);
static Rboolean (*R_isVector) (SEXP);
static Rboolean (*R_isLogical) (SEXP);
static Rboolean (*R_isInteger) (SEXP);
static Rboolean (*R_isReal) (SEXP);
static Rboolean (*R_isString) (SEXP);
static Rboolean (*R_isList) (SEXP);

static int (*R_initEmbeddedR) (int, char **);
static int (*R_ncols) (SEXP);
static int (*R_nrows) (SEXP);
static int (*R_length) (SEXP);
static int (*R_TYPEOF) (SEXP);
static int *(*R_LOGICAL) (SEXP);

static void (*R_endEmbeddedR) (int);
static void (*R_unprotect) (int);
static void (*R_PrintValue) (SEXP);
static void (*R_SET_TYPEOF) (SEXP, int);
static void (*R_SET_TAG) (SEXP, SEXP);

static SEXP (*R_VECTOR_ELT) (SEXP, int);
static void (*R_SET_VECTOR_ELT) (SEXP, int, SEXP);

static const char *(*R_errmsg) (void);

#ifdef WIN32
static char *(*R_get_HOME) (void);
#endif

/* utility function to cumulate errors from dlsym */

static void *dlget (void *handle, const char *name, int *err)
{
    void *p = gretl_dlsym(handle, name);

    if (p == NULL) {
        fprintf(stderr, "dlget: couldn't find '%s'\n", name);
        *err += 1;
    }

    return p;
}

/* dlopen the R library and grab all the symbols we need:
   several function pointers and a few global variables
*/

static int load_R_symbols (void)
{
    const char *libpath = gretl_rlib_path();
    int err = 0;

#if FDEBUG
    fprintf(stderr, "Loading libR symbols from '%s'\n", libpath);
#endif

    Rhandle = gretl_dlopen(libpath, 1);
    if (Rhandle == NULL) {
        err = E_EXTERNAL;
        goto bailout;
    }

    R_CDR           = dlget(Rhandle, "CDR", &err);
    R_REAL          = dlget(Rhandle, "REAL", &err);
    R_INT           = dlget(Rhandle, "INTEGER", &err);
    R_STRING        = dlget(Rhandle, "R_CHAR", &err);
    R_STRING_PTR    = dlget(Rhandle, "STRING_PTR", &err);
    R_STRING_ELT    = dlget(Rhandle, "STRING_ELT", &err);
    R_allocList     = dlget(Rhandle, "Rf_allocList", &err);
    R_allocMatrix   = dlget(Rhandle, "Rf_allocMatrix", &err);
    R_allocVector   = dlget(Rhandle, "Rf_allocVector", &err);
    R_endEmbeddedR  = dlget(Rhandle, "Rf_endEmbeddedR", &err);
    R_findFun       = dlget(Rhandle, "Rf_findFun", &err);
    R_findVar       = dlget(Rhandle, "Rf_findVar", &err);
    R_initEmbeddedR = dlget(Rhandle, "Rf_initEmbeddedR", &err);
    R_install       = dlget(Rhandle, "Rf_install", &err);
    R_isMatrix      = dlget(Rhandle, "Rf_isMatrix", &err);
    R_isVector      = dlget(Rhandle, "Rf_isVector", &err);
    R_isLogical     = dlget(Rhandle, "Rf_isLogical", &err);
    R_isInteger     = dlget(Rhandle, "Rf_isInteger", &err);
    R_isReal        = dlget(Rhandle, "Rf_isReal", &err);
    R_isString      = dlget(Rhandle, "Rf_isString", &err);
    R_isList        = dlget(Rhandle, "Rf_isList", &err);
    R_mkString      = dlget(Rhandle, "Rf_mkString", &err);
    R_mkChar        = dlget(Rhandle, "Rf_mkChar", &err);
    R_mkNamed       = dlget(Rhandle, "Rf_mkNamed", &err);
    R_ncols         = dlget(Rhandle, "Rf_ncols", &err);
    R_nrows         = dlget(Rhandle, "Rf_nrows", &err);
    R_length        = dlget(Rhandle, "Rf_length", &err);
    R_getAttrib     = dlget(Rhandle, "Rf_getAttrib", &err);
    R_setAttrib     = dlget(Rhandle, "Rf_setAttrib", &err);
    R_PrintValue    = dlget(Rhandle, "Rf_PrintValue", &err);
    R_protect       = dlget(Rhandle, "Rf_protect", &err);
    R_ScalarReal    = dlget(Rhandle, "Rf_ScalarReal", &err);
    R_unprotect     = dlget(Rhandle, "Rf_unprotect", &err);
    R_catch         = dlget(Rhandle, "R_tryEval", &err);
    R_errmsg        = dlget(Rhandle, "R_curErrorBuf", &err);
    R_SETCAR        = dlget(Rhandle, "SETCAR", &err);
    R_SET_TYPEOF    = dlget(Rhandle, "SET_TYPEOF", &err);
    R_TYPEOF        = dlget(Rhandle, "TYPEOF", &err);
    R_SET_TAG       = dlget(Rhandle, "SET_TAG", &err);
    R_LOGICAL       = dlget(Rhandle, "LOGICAL", &err);
    R_VECTOR_ELT    = dlget(Rhandle, "VECTOR_ELT", &err);
    R_SET_VECTOR_ELT = dlget(Rhandle, "SET_VECTOR_ELT", &err);

#ifdef WIN32
    R_get_HOME = dlget(Rhandle, "get_R_HOME", &err);
#endif

    if (!err) {
        PR_GlobalEnv    = (SEXP *) dlget(Rhandle, "R_GlobalEnv", &err);
        PR_NilValue     = (SEXP *) dlget(Rhandle, "R_NilValue", &err);
        PR_UnboundValue = (SEXP *) dlget(Rhandle, "R_UnboundValue", &err);
        PR_NamesSymbol  = (SEXP *) dlget(Rhandle, "R_NamesSymbol", &err);
        PR_DimNamesSymbol = (SEXP *) dlget(Rhandle, "R_DimNamesSymbol", &err);
    }

    if (err) {
        fprintf(stderr, "ERROR loading R symbols\n");
        close_plugin(Rhandle);
        Rhandle = NULL;
        err = E_EXTERNAL;
    }

 bailout:

#if FDEBUG || defined(WIN32)
    fprintf(stderr, "load_R_symbols: returning %d\n", err);
#endif

    return err;
}

void gretl_R_cleanup (void)
{
#if FDEBUG
    fprintf(stderr, "gretl_R_cleanup: Rinit = %d\n", Rinit);
#endif

    if (Rinit) {
        R_endEmbeddedR(0);
        close_plugin(Rhandle);
        Rhandle = NULL;
    }
}

/* called from gretl_paths.c on revising the Rlib path:
   allow for the possibility that the path was wrong but is
   now OK
*/

void gretl_R_reset_error (void)
{
    Rlib_err = 0;
}

/* Code for trying to ensure that there's a correct
   value for R_HOME in the environment.
*/

static char RHOME[MAXLEN];

#ifdef G_OS_WIN32

static void set_path_for_Rlib (const char *libpath)
{
    char *oldpath = getenv("PATH");
    gchar *Rpath = g_strdup(libpath);
    gchar *s;

    /* chop off the last libpath element */
    s = strrchr(Rpath, '\\');
    *s = '\0';

    if (oldpath != NULL && strstr(oldpath, Rpath) != NULL) {
        ; /* nothing to be done */
    } else if (oldpath == NULL) {
        /* very unlikely, but... */
        gretl_setenv("PATH", Rpath);
    } else {
        gchar *modpath;

        fprintf(stderr, "Adding '%s' to PATH\n", Rpath);
        modpath = g_strdup_printf("%s;%s", oldpath, Rpath);
        gretl_setenv("PATH", modpath);
        g_free(modpath);
    }

    g_free(Rpath);
}

int check_set_R_home (void)
{
    const char *path = gretl_rlib_path();
    const char *orig = path;
    char tmp[MAX_PATH];
    int err = 0;

    if (RHOME[0] != '\0') {
	/* we've already done this */
	return 0;
    }

    err = gretl_test_fopen(path, "r");
    /* an error here is not necessarily fatal, if we can get
       the R path from the registry.
    */
    if (err) {
	path = NULL;
	err = 0;
    }

    err = R_home_from_registry(tmp);
    if (!err) {
	strcpy(RHOME, tmp);
	err = win32_R_path(tmp, RLIB);
	path = tmp;
    }

    if (path == NULL) {
	gretl_errmsg_set(_("Can't find the R shared library"));
	return E_EXTERNAL;
    }

    if (!err) {
	if (strcmp(path, orig)) {
	    gretl_set_path_by_name("rlibpath", path);
	}
	if (RHOME[0] == '\0') {
	    char *s;

	    strcpy(RHOME, path);
	    s = strstr(RHOME, "\\bin");
	    if (s != NULL) {
		*s = '\0';
	    }
	}
	/* put R_HOME into env */
	fprintf(stderr, "setting R_HOME='%s'\n", RHOME);
	gretl_setenv("R_HOME", RHOME);
        set_path_for_Rlib(path);
    }

    return err;
}

#else /* non-Windows variant */

static int absolutize_R_path (char *targ, const char *path)
{
    gchar *s = g_strdup(path);
    gchar *p = strrchr(s, '/');
    int err = 0;

    if (p != NULL) {
        gchar *abspath;

        *p = '\0';
        abspath = g_strdup_printf("%s/%s", s, targ);
        strcpy(targ, abspath);
        g_free(abspath);
        fprintf(stderr, " absolute path '%s'\n", targ);
    } else {
        err = E_EXTERNAL;
    }

    g_free(s);

    return err;
}

int check_set_R_home (void)
{
    const char *Rhome = NULL;
    const char *path = gretl_rlib_path();
    const char *orig = path;
    struct stat buf = {0};
    gchar *lp = NULL;
    char tmp[2048] = {0};
    int do_setenv = 0;
    int err = 0;

    if (RHOME[0] != '\0') {
	/* we've already done this */
	return 0;
    }

    err = (lstat(path, &buf) != 0);
    if (err) {
	fprintf(stderr, "R_home: user path '%s', err %d\n", path, err);
    }
    /* an error here is not necessarily fatal, if R_HOME
       is defined and works OK.
    */
    if (err) {
	path = NULL;
	err = 0;
    }

    Rhome = getenv("R_HOME");

    if (Rhome == NULL || *Rhome == '\0') {
	fprintf(stderr, "No R_HOME\n");
	do_setenv = 1;
    } else {
	/* verify R_HOME */
#ifdef OS_OSX
        const char *libname = "libR.dylib";
#else
        const char *libname = "libR.so";
#endif
	lp = g_strdup_printf("%s/lib/%s", Rhome, libname);
	err = (lstat(lp, &buf) != 0);
	if (err) {
	    /* needs correcting */
	    do_setenv = 1;
	} else {
	    strcpy(RHOME, Rhome);
	    path = lp;
	}
    }

    if (path == NULL) {
	g_free(lp);
	gretl_errmsg_set(_("Can't find the R shared library"));
	return E_EXTERNAL;
    }

    if (S_ISLNK(buf.st_mode)) {
	/* the path is a symlink */
	ssize_t n = readlink(path, tmp, sizeof(tmp) - 1);

	fprintf(stderr, "The R library path '%s' seems to be a symlink\n", path);
	if (n < 0 || n >= sizeof tmp) {
	    fprintf(stderr, " couldn't resolve it (n=%d, tmp='%s'\n", (int) n, tmp);
	    err = E_EXTERNAL;
	} else {
	    fprintf(stderr, " resolved to '%s'\n", tmp);
	    if (!g_path_is_absolute(tmp)) {
		err = absolutize_R_path(tmp, path);
	    }
	}
    } else if (S_ISREG(buf.st_mode)) {
	/* the given path is a regular file */
	strcpy(tmp, path);
    }

    if (!err && strcmp(tmp, orig)) {
	gretl_set_path_by_name("rlibpath", tmp);
    }

    if (!err && do_setenv) {
	char *s = strstr(tmp, "/lib/libR");

	if (s != NULL) {
	    *s = '\0';
	    fprintf(stderr, "setting R_HOME='%s'\n", tmp);
	    gretl_setenv("R_HOME", tmp);
	    strcpy(RHOME, tmp);
	} else {
	    gretl_errmsg_set(_("Can't determine R_HOME"));
	    err = E_EXTERNAL;
	}
    }

    g_free(lp);

    return err;
}

#endif /* Windows vs non-Windows variants */

/* Initialize the R library for use with gretl.  Note that we only
   need do this once per gretl session.  We need to check that the
   environment is set to R's liking first, otherwise initialization
   will fail -- and will abort gretl too!
*/

static int gretl_Rlib_init (void)
{
    int err = 0;

#if FDEBUG
    fprintf(stderr, "gretl_Rlib_init: starting\n");
#endif

    err = check_set_R_home();

    if (!err) {
	err = load_R_symbols();
    }

    if (err) {
        fprintf(stderr, "gretl_Rlib_init: failed to load R functions\n");
        return E_EXTERNAL;
    }

    /* ensure common filenames are in place */
    make_gretl_R_names();

    /* and ensure that gretl.Rprofile doesn't get in the way */
    gretl_remove(gretl_Rprofile);

    if (!err) {
        char *argv[] = {
            "gretl",
            "--no-save",
            "--silent",
        };
        int ok, argc = 3;

#if FDEBUG
        fprintf(stderr, "calling R_initEmbeddedR\n");
#endif
        ok = R_initEmbeddedR(argc, argv);
        if (ok) {
            VR_GlobalEnv = *PR_GlobalEnv;
            VR_NilValue = *PR_NilValue;
            VR_UnboundValue = *PR_UnboundValue;
            VR_NamesSymbol = *PR_NamesSymbol;
            VR_DimNamesSymbol = *PR_DimNamesSymbol;
            Rinit = 1;
        } else {
            close_plugin(Rhandle);
            Rhandle = NULL;
            err = Rlib_err = E_EXTERNAL;
        }
    }

#if FDEBUG
    fprintf(stderr, "gretl_Rlib_init: returning %d\n", err);
#endif

    return err;
}

/* run R's source() function on an R command file written by
   gretl, shared library version */

static int lib_run_Rlib_sync (gretlopt opt, PRN *prn)
{
    int err = 0;

#if FDEBUG
    fprintf(stderr, "lib_run_Rlib_sync: starting\n");
#endif

    if (!Rinit) {
        err = gretl_Rlib_init();
    }

    if (!err) {
        SEXP expr, p;

        /* make echo/print.eval argument */
        R_protect(p = R_allocVector(LGLSXP, 1));
        R_LOGICAL(p)[0] = TRUE;

        /* expression source(f, print.eval=p) */
        R_protect(expr = R_allocVector(LANGSXP, 3));
        R_SETCAR(expr, R_install("source"));
        R_SETCAR(R_CDR(expr), R_mkString(gretl_Rsrc));
        R_SETCAR(R_CDR(R_CDR(expr)), p);
        R_SET_TAG(R_CDR(R_CDR(expr)),
                  R_install((opt & OPT_V)? "echo" : "print.eval"));

        R_catch(expr, NULL, &err);
        R_unprotect(2);
    }

    if (prn != NULL) {
        const gchar *outname;
        FILE *fp;

        outname = (err)? gretl_Rmsg : gretl_Rout;
        fp = gretl_fopen(outname, "r");

        if (fp != NULL) {
            char line[512];

            while (fgets(line, sizeof line, fp)) {
#ifdef G_OS_WIN32
                win32_put_R_output_line(line, prn);
#else
                pputs(prn, line);
#endif
            }
            fclose(fp);
            gretl_remove(outname);
        }
    }

#if FDEBUG
    fprintf(stderr, "lib_run_Rlib_sync: returning %d\n", err);
#endif

    return (err)? E_EXTERNAL : 0;
}

static SEXP find_R_function (const char *name)
{
    SEXP fun;
    SEXPTYPE t;

    fun = R_findVar(R_install(name), VR_GlobalEnv);
    t = R_TYPEOF(fun);

    if (t == PROMSXP) {
        /* eval promise if need be */
        int err = 1;

        fun = R_catch(fun, VR_GlobalEnv, &err);
        if (!err) {
            t = R_TYPEOF(fun);
        }
    }

    if (t != CLOSXP && t != BUILTINSXP && t != SPECIALSXP) {
        return VR_UnboundValue;
    }

    return fun;
}

/* Check if we should be using the R shared library for executing the
   code in a "foreign" block.  This is disabled if the user has
   done "set R_lib off", and can be prohibited by the environment
   variable GRETL_NO_RLIB.  It may also be blocked if we already tried
   and failed to initialize the library for gretl's use.  (The
   fallback will be to call the R binary.)
*/

static int gretl_use_Rlib (void)
{
    int ret = 0;

#if FDEBUG
    fprintf(stderr, "gretl_use_Rlib: starting\n");
#endif

    if (!Rlib_err && libset_get_bool(R_LIB) && !getenv("GRETL_NO_RLIB")) {
        /* use of library is not blocked */
        if (Rinit) {
            /* already opened, fine */
            ret = 1;
        } else {
            /* try opening library */
            Rlib_err = gretl_Rlib_init();
            ret = !Rlib_err;
        }
    }

#if FDEBUG
    fprintf(stderr, "gretl_use_Rlib: using %s\n", (ret)? "library" : "executable");
#endif

    return ret;
}

/* Used in "genr", to see if @name denotes an R function,
   either built-in or possibly user-defined.  The lookup
   is conditional on the user's doing "set R_functions on".
*/

int get_R_function_by_name (const char *name)
{
    int ret = 0;

    if (libset_get_bool(R_FUNCTIONS) && gretl_use_Rlib()) {
        SEXP fun = find_R_function(name);

        ret = (fun == VR_UnboundValue)? 0 : 1;
    }

    return ret;
}

/* R_function_add... : these functions are used to convert
   from gretl types to R constructs for passing to R
   functions. Theye are called (indirectly) from geneval.c
*/

static int R_function_add_scalar (double x)
{
    R_SETCAR(current_arg, R_ScalarReal(x));
    return 0;
}

static int R_function_add_string (const char *s)
{
    R_SETCAR(current_arg, R_mkString(s));
    return 0;
}

static int R_function_add_vector (const double *x, int t1, int t2)
{
    SEXP res = R_allocVector(REALSXP, t2 - t1 + 1);
    int i;

    if (res == NULL) {
        return E_ALLOC;
    }

    for (i=t1; i<=t2; i++) {
        R_REAL(res)[i-t1] = x[i];
    }

    R_SETCAR(current_arg, res);

    return 0;
}

static int R_function_add_factor (const DATASET *dset, int v)
{
    SEXP res = R_allocVector(STRSXP, dset->t2 - dset->t1 + 1);
    const char *si;
    int i;

    if (res == NULL) {
        return E_ALLOC;
    }

    for (i=dset->t1; i<=dset->t2; i++) {
        si = series_get_string_for_obs(dset, v, i);
        if (si == NULL) {
            R_STRING_PTR(res)[i-dset->t1] = R_mkChar("");
        } else {
            R_STRING_PTR(res)[i-dset->t1] = R_mkChar(si);
        }
    }

    R_SETCAR(current_arg, res);

    return 0;
}

static SEXP make_R_array (gretl_array *a, int *err)
{
    GretlType t = gretl_array_get_type(a);
    int n = gretl_array_get_length(a);
    const char *si;
    SEXP as;
    int i;

    if (t != GRETL_TYPE_STRINGS) {
        gretl_errmsg_set(_("Only strings arrays are accepted as R-function arguments"));
        *err = E_TYPES;
        return NULL;
    }

    as = R_allocVector(STRSXP, n);
    if (as == NULL) {
        *err = E_ALLOC;
    } else {
        for (i=0; i<n; i++) {
            si = gretl_array_get_data(a, i);
            R_STRING_PTR(as)[i] = R_mkChar(si);
        }
    }

    return as;
}

static SEXP make_R_strings (const char **S, int ns)
{
    SEXP as = R_allocVector(STRSXP, ns);
    int i;

    if (as != NULL) {
        for (i=0; i<ns; i++) {
            R_STRING_PTR(as)[i] = R_mkChar(S[i]);
        }
    }

    return as;
}

static void maybe_add_export_dimnames (SEXP ms, const gretl_matrix *m)
{
    const char **S;
    SEXP rn = NULL;
    SEXP cn = NULL;

    S = gretl_matrix_get_rownames(m);
    if (S != NULL) {
	rn = make_R_strings(S, m->rows);
    }
    S = gretl_matrix_get_colnames(m);
    if (S != NULL) {
	cn = make_R_strings(S, m->cols);
    }
    if (rn != NULL || cn != NULL) {
	/* FIXME PROTECT? */
	SEXP dimnames = R_allocVector(VECSXP, 2);

	if (rn != NULL) {
	    R_SET_VECTOR_ELT(dimnames, 0, rn);
	}
	if (cn != NULL) {
	    R_SET_VECTOR_ELT(dimnames, 1, cn);
	}
	R_setAttrib(ms, VR_DimNamesSymbol, dimnames);
    }
}

static SEXP make_R_matrix (const gretl_matrix *m, int *err)
{
    SEXP ms = R_allocMatrix(REALSXP, m->rows, m->cols);

    if (ms == NULL) {
        *err = E_ALLOC;
    } else {
	size_t sz = m->rows * m->cols * sizeof *m->val;

	memcpy(R_REAL(ms), m->val, sz);
	maybe_add_export_dimnames(ms, m);
    }

    return ms;
}

static SEXP make_R_bundle (gretl_bundle *b, int *err)
{
    GretlType type;
    int size;
    void *ptr;
    char **keys;
    SEXP res;
    int i, n;

    if (gretl_bundle_get_n_members(b) == 0) {
        const char *nokeys[] = { "" };

        return R_mkNamed(VECSXP, nokeys);
    }

    keys = gretl_bundle_get_keys_raw(b, &n);
    if (keys == NULL) {
        *err = E_DATA;
        return NULL;
    }

    keys[n] = ""; /* the termination wanted by R */
    res = R_mkNamed(VECSXP, (const char **) keys);

    for (i=0; i<n && !*err; i++) {
        ptr = gretl_bundle_get_data(b, keys[i], &type, &size, err);
        if (*err) {
            break;
        }
        if (type == GRETL_TYPE_DOUBLE) {
            double x = *(double *) ptr;

            R_SET_VECTOR_ELT(res, i, R_ScalarReal(x));
        } else if (type == GRETL_TYPE_STRING) {
            R_SET_VECTOR_ELT(res, i, R_mkString(ptr));
        } else if (type == GRETL_TYPE_MATRIX) {
            SEXP ms = make_R_matrix(ptr, err);

            if (!*err) {
                R_SET_VECTOR_ELT(res, i, ms);
            }
        } else if (type == GRETL_TYPE_ARRAY) {
            SEXP as = make_R_array(ptr, err);

            if (!*err) {
                R_SET_VECTOR_ELT(res, i, as);
            }
        } else if (type == GRETL_TYPE_BUNDLE) {
            SEXP bs = make_R_bundle(ptr, err);

            if (!*err) {
                R_SET_VECTOR_ELT(res, i, bs);
            }
        } else {
            gretl_errmsg_sprintf(_("%s: not handled\n"),
                                 gretl_type_get_name(type));
            *err = E_TYPES;
        }
    }

    keys[n] = NULL;
    strings_array_free(keys, n);

    return res;
}

static int R_function_add_object (void *ptr, GretlType t)
{
    SEXP res = NULL;
    int err = 0;

    if (t == GRETL_TYPE_MATRIX) {
        res = make_R_matrix(ptr, &err);
    } else if (t == GRETL_TYPE_ARRAY) {
        res = make_R_array(ptr, &err);
    } else if (t == GRETL_TYPE_BUNDLE) {
        res = make_R_bundle(ptr, &err);
    } else {
        gretl_errmsg_sprintf(_("%s: not handled as R-function argument\n"),
                             gretl_type_get_name(t));
        err = E_TYPES;
    }

    if (!err) {
        R_SETCAR(current_arg, res);
    }

    return err;
}

/* public because called from geneval.c */

int gretl_R_function_add_arg (void *ptr, GretlType type)
{
    int err;

    current_arg = R_CDR(current_arg);

    if (type == GRETL_TYPE_DOUBLE) {
        double x = *(double *) ptr;

        err = R_function_add_scalar(x);
    } else if (type == GRETL_TYPE_STRING) {
        const char *s = ptr;

        err = R_function_add_string(s);
    } else {
        /* matrix, array, bundle */
        err = R_function_add_object(ptr, type);
    }

    return err;
}

/* public because called from geneval.c */

int gretl_R_function_add_series (double *x, const DATASET *dset, int v)
{
    int err;

    current_arg = R_CDR(current_arg);

    if (v > 0 && is_string_valued(dset, v)) {
        err = R_function_add_factor(dset, v);
    } else {
        err = R_function_add_vector(x, dset->t1, dset->t2);
    }

    return err;
}

/* called from geneval.c only, and should be pre-checked */

int gretl_R_get_call (const char *name, int argc)
{
    SEXP call, e;

    call = R_findFun(R_install(name), VR_GlobalEnv);

    if (call == VR_NilValue) {
        fprintf(stderr, "gretl_R_get_call: no definition for function %s\n",
                name);
        R_unprotect(1); /* is this OK? */
        return E_EXTERNAL;
    }

    R_protect(e = R_allocList(argc + 1));
    R_SET_TYPEOF(e, LANGSXP);
    R_SETCAR(e, R_install(name));
    current_call = current_arg = e;

    return 0;
}

static int numeric_ok (SEXP s)
{
    return R_isReal(s) || R_isInteger(s) || R_isLogical(s);
}

#define R_NULL(p) (p == NULL || (SEXP) p == VR_NilValue)

static GretlType R_type_to_gretl_type (SEXP s, const char *name, int *err);
static void *object_from_R (SEXP res, const char *name, GretlType type,
                            int *err);

static void maybe_add_import_dimnames (SEXP s, gretl_matrix *m)
{
    gretl_array *a;
    char **S = NULL;
    SEXP dimnames;
    SEXP d[2];
    GretlType t;
    int i, ns;
    int err = 0;

    dimnames = R_getAttrib(s, VR_DimNamesSymbol);

    if (R_NULL(dimnames)) {
        return;
    }

    d[0] = R_VECTOR_ELT(dimnames, 0);
    d[1] = R_VECTOR_ELT(dimnames, 1);

    for (i=0; i<2; i++) {
	if (R_NULL(d[i])) {
	    continue;
	}
	t = R_type_to_gretl_type(d[i], "dimnames", &err);
	if (t == GRETL_TYPE_ARRAY) {
	    a = object_from_R(d[i], "dimnames", t, &err);
	    ns = i == 0 ? m->rows : m->cols;
	    if (a != NULL && gretl_array_get_length(a) == ns) {
		S = gretl_array_steal_strings(a, &ns);
		if (i == 0) {
		    gretl_matrix_set_rownames(m, S);
		} else {
		    gretl_matrix_set_colnames(m, S);
		}
		gretl_array_destroy(a);
	    }
	}
    }
}

static gretl_matrix *matrix_from_R (SEXP s, const char *name,
                                    int *err)
{
    gretl_matrix *m = NULL;
    int nr = R_nrows(s);
    int nc = R_ncols(s);

#if FDEBUG
    fprintf(stderr, "matrix_from_R: name='%s', %d x %d\n", name, nr, nc);
#endif

    if (nr >= 0 && nc >= 0) {
        m = gretl_matrix_alloc(nr, nc);
        if (m == NULL) {
            *err = E_ALLOC;
        }
    } else {
        gretl_errmsg_sprintf(_("%s: invalid matrix dimensions, %d x %d"),
                             name, nr, nc);
        *err = E_DATA;
    }

    if (!*err && R_isReal(s)) {
	size_t sz = nr * nc * sizeof *m->val;

	memcpy(m->val, R_REAL(s), sz);
    } else if (!*err) {
        int i, j;

        for (j=0; j<nc; j++) {
            for (i=0; i<nr; i++) {
		gretl_matrix_set(m, i, j, R_INT(s)[i + j * nr]);
            }
        }
    }

    if (!*err) {
	maybe_add_import_dimnames(s, m);
    }

    return m;
}

/* at present only arrays of strings are handled */

static gretl_array *array_from_R (SEXP res, const char *name,
                                  int *err)
{
    gretl_array *a = NULL;
    int nr = R_nrows(res);

    if (nr >= 0) {
        a = gretl_array_new(GRETL_TYPE_STRINGS, nr, err);
    } else {
        gretl_errmsg_sprintf(_("%s: invalid array length %d"),
                             name, nr);
        *err = E_DATA;
    }
    if (a != NULL && nr > 0) {
        const char *s;
        int i;

        for (i=0; i<nr && !*err; i++) {
            s = R_STRING(R_STRING_PTR(res)[i]);
            *err = gretl_array_set_string(a, i, (char *) s, 1);
        }
    }

    return a;
}

static int vector_can_be_bundle (SEXP s, const char *name)
{
    SEXP si, names;
    const char *key;
    GretlType gtype;
    int i, k;
    int err = 0;

    if (R_TYPEOF(s) != VECSXP) {
        gretl_errmsg_sprintf(_("%s: not generic vector, can't make into bundle"), name);
        return 0;
    }

    k = R_nrows(s);
    if (k == 0) {
        /* empty bundle, should be OK */
        return 1;
    }

    names = R_getAttrib(s, VR_NamesSymbol);
    if (R_NULL(names)) {
        gretl_errmsg_sprintf(_("%s: R list has no tags, can't make into bundle"), name);
        return 0;
    }

    for (i=0; i<k && !err; i++) {
        si = R_VECTOR_ELT(s, i);
        key = R_STRING(R_STRING_ELT(names, i));
        if (R_NULL(key)) {
            err = E_TYPES;
            break;
        }
#if 0
        int rt = R_TYPEOF(si);
        fprintf(stderr, "element %d, R-type %d, key '%s'\n", i, rt, key);
#endif
        gtype = R_type_to_gretl_type(si, key, &err);
        if (!err && gtype == GRETL_TYPE_NONE) {
            err = E_TYPES;
        }
    }

    return err == 0;
}

static gretl_bundle *bundle_from_R (SEXP s, int *err)
{
    int i, k = R_nrows(s);
    SEXP si, names = R_getAttrib(s, VR_NamesSymbol);
    const char *key;
    GretlType type;
    void *ptr;
    gretl_bundle *b;

    b = gretl_bundle_new();
    if (b == NULL) {
        *err = E_ALLOC;
        return NULL;
    }

    for (i=0; i<k && !*err; i++) {
        si = R_VECTOR_ELT(s, i);
        key = R_STRING(R_STRING_ELT(names, i));
        type = R_type_to_gretl_type(si, key, err);
        if (!*err) {
            ptr = object_from_R(si, key, type, err);
        }
        if (!*err) {
            *err = gretl_bundle_set_data(b, key, ptr, type, 0);
        }
    }

    if (*err) {
        gretl_bundle_destroy(b);
        b = NULL;
    }

    return b;
}

static void *object_from_R (SEXP res, const char *name,
                            GretlType type, int *err)
{
    if (type == GRETL_TYPE_MATRIX) {
        return matrix_from_R(res, name, err);
    } else if (gretl_scalar_type(type)) {
        return R_REAL(res);
    } else if (type == GRETL_TYPE_STRING) {
        SEXP *rsp = R_STRING_PTR(res);

        return (void *) R_STRING(*rsp);
    } else if (type == GRETL_TYPE_ARRAY) {
        return array_from_R(res, name, err);
    } else if (type == GRETL_TYPE_BUNDLE) {
        return bundle_from_R(res, err);
    } else if (type != GRETL_TYPE_NONE) {
        *err = E_TYPES;
    }

    return NULL;
}

static GretlType R_type_to_gretl_type (SEXP s, const char *name, int *err)
{
    GretlType t = GRETL_TYPE_NONE;

    if (R_isList(s)) {
        fprintf(stderr, " %s: list, length %d\n", name, R_length(s));
    } else if (R_isMatrix(s)) {
        if (numeric_ok(s)) {
            t = GRETL_TYPE_MATRIX;
        } else {
            gretl_errmsg_sprintf(_("%s: got 'matrix' result, but not numeric"), name);
            *err = E_TYPES;
        }
    } else if (R_isVector(s)) {
        if (numeric_ok(s)) {
            t = GRETL_TYPE_MATRIX;
        } else if (R_isString(s)) {
            if (R_nrows(s) == 1) {
                t = GRETL_TYPE_STRING;
            } else {
                t = GRETL_TYPE_ARRAY;
            }
        } else if (vector_can_be_bundle(s, name)) {
            t = GRETL_TYPE_BUNDLE;
        } else {
            *err = E_TYPES;
        }
    } else if (R_isLogical(s)) {
        t = GRETL_TYPE_BOOL;
    } else if (R_isInteger(s)) {
        t = GRETL_TYPE_INT;
    } else if (R_isReal(s)) {
        t = GRETL_TYPE_DOUBLE;
    } else if (R_isString(s)) {
        t = GRETL_TYPE_STRING;
    }

    return t;
}

/* execute an R function and try to convert the value returned
   into a gretl type
*/

int gretl_R_function_exec (const char *name, int *rtype, void **ret)
{
    void *data = NULL;
    SEXP res;
    int err = 0;

#if FDEBUG
    printf("gretl_R_function_exec...\n");
    R_PrintValue(current_call);
#endif

    if (0 /* gretl_messages_on() */) {
        R_PrintValue(current_call);
    }

    res = R_catch(current_call, VR_GlobalEnv, &err);
    if (err) {
        const char *Rmsg = R_errmsg();

        fprintf(stderr, "gretl_R_function_exec: R_catch failed on %s\n", name);
        if (Rmsg != NULL) {
            gretl_errmsg_set(Rmsg);
        }
        return E_EXTERNAL;
    }

    *rtype = R_type_to_gretl_type(res, name, &err);

#if FDEBUG
    printf("R return value: got type %d (%s)\n", *rtype,
           gretl_type_get_name(*rtype));
    printf("Calling R_PrintValue() on @res\n");
    R_PrintValue(res);
#endif

    if (!err) {
        data = object_from_R(res, name, *rtype, &err);
    }

    if (!err) {
        if (gretl_scalar_type(*rtype)) {
            double *dret = *ret;

            *dret = *(double *) data;
        } else if (*rtype == GRETL_TYPE_STRING) {
            *ret = gretl_strdup(data);
        } else {
            /* matrix, array, bundle */
            *ret = data;
        }
        R_unprotect(1);
    }

    return err;
}

static int run_R_lib (const char *buf,
                      const DATASET *dset,
                      gretlopt opt,
                      PRN *prn)
{
    int err;

#if FDEBUG
    fprintf(stderr, "run_R_lib\n");
#endif

    /* we don't want gretl.Rprofile in the way */
    gretl_remove(gretl_Rprofile);

    gretl_push_c_numeric_locale();

    /* by passing OPT_L below we indicate that we're
       using the library */
    err = write_R_source_file(buf, dset, opt | OPT_L);
    if (!err) {
        err = lib_run_Rlib_sync(opt, prn);
    }

    gretl_pop_c_numeric_locale();

    return err;
}

#endif /* USE_RLIB */

/**
 * foreign_start:
 * @ci: either FOREIGN or MPI.
 * @param: string specifying language, for FOREIGN.
 * @opt: may include %OPT_V for verbose operation.
 * @prn: struct for printing output.
 *
 * Starts a new "foreign" block if no such block is
 * currently defined.
 *
 * Returns: 0 on success, non-zero on error.
 */

int foreign_start (int ci, const char *param, gretlopt opt,
                   PRN *prn)
{
    int err = 0;

    if (foreign_started) {
        gretl_errmsg_sprintf(_("%s: a block is already started"),
                             gretl_command_word(ci));
        return E_DATA;
    }

    foreign_opt = OPT_NONE;

    if (ci == FOREIGN) {
        if (param == NULL || *param == '\0') {
            err = E_ARGS;
        } else {
            char lang[16];

            if (sscanf(param, "language=%15s", lang) == 1) {
                err = set_foreign_lang(lang, prn);
            } else {
                err = E_PARSE;
            }
        }
    } else if (ci == MPI) {
        err = set_foreign_lang("mpi", prn);
    }

    if (!err) {
        foreign_started = 1;
        foreign_opt = opt;
    }

    return err;
}

/**
 * foreign_append:
 * @line: line to append.
 * @context: either FOREIGN or MPI.
 *
 * Appends @line to an internally stored block of "foreign"
 * or MPI commands, if such a block is currently defined.
 *
 * Returns: 0 on success, non-zero on error.
 */

int foreign_append (const char *line, int context)
{
    int err = 0;

#if 0
    fprintf(stderr, "foreign_append: '%s'\n", line);
#endif

    if (!foreign_started) {
        gretl_errmsg_sprintf(_("%s: no block is in progress"),
                             gretl_command_word(context));
        err = E_DATA;
    } else if (!string_is_blank(line)) {
        err = strings_array_add(&foreign_lines, &foreign_n_lines, line);
        if (err) {
            foreign_destroy();
        }
    }

    return err;
}

/* write profile (perhaps) and Rsrc files */

static int run_R_binary (const char *buf,
                         const DATASET *dset,
                         gretlopt opt,
                         PRN *prn)
{
    int err;

    gretl_push_c_numeric_locale();
    err = write_gretl_R_files(buf, dset, opt);

    if (err) {
        delete_gretl_R_files();
    } else {
#ifdef G_OS_WIN32
        err = win32_lib_run_R_sync(opt, prn);
#else
        err = lib_run_R_sync(opt, prn);
#endif
    }

    gretl_pop_c_numeric_locale();

    return err;
}

static int write_foreign_io_file (int lang, PRN *prn)
{
    const char *fname;
    int err = 0;

    fname = get_optval_string(FOREIGN, OPT_I);
    if (fname == NULL) {
        err = E_ARGS;
    } else {
        fname = gretl_maybe_switch_dir(fname);
        err = ensure_foreign_io_file(lang, fname);
        if (!err && gretl_messages_on()) {
            pprintf(prn, _("Wrote '%s'\n"), fname);
        }
    }

    return err;
}

#ifdef HAVE_MPI

static void try_for_mpi_errmsg (PRN *prn)
{
    gchar *tmp = gretl_make_dotpath("mpi.fail");
    gchar *msg = NULL;

    g_file_get_contents(tmp, &msg, NULL, NULL);
    if (msg != NULL) {
        pputs(prn, msg);
        g_free(msg);
    }
    gretl_remove(tmp);
    g_free(tmp);
}

#endif

/**
 * foreign_execute:
 * @dset: dataset struct.
 * @opt: may include %OPT_V for verbose operation
 * @prn: struct for printing output.
 *
 * Executes a block of commands previously established via
 * calls to foreign_append_line().
 *
 * Returns: 0 on success, non-zero on error.
 */

int foreign_execute (DATASET *dset, gretlopt opt, PRN *prn)
{
#if defined(MPI_REALTIME)
    static struct iodata io;
    void *ptr = &io;
#elif defined(HAVE_MPI) && !defined(G_OS_WIN32)
    void *ptr = NULL;
#endif
    int i, err = 0;

    if (opt & OPT_I) {
        /* just write IO file */
        return write_foreign_io_file(foreign_lang, prn);
    }

    if (foreign_lang == LANG_R) {
        make_gretl_R_names();
    }

    if (opt & OPT_V) {
        /* verbose: echo the stored commands */
        for (i=0; i<foreign_n_lines; i++) {
            pprintf(prn, "> %s\n", foreign_lines[i]);
        }
    }

    foreign_opt |= opt;

#ifdef HAVE_MPI
    if (foreign_lang == LANG_MPI) {
        err = write_gretl_mpi_script(foreign_opt, dset);
        if (err) {
            delete_mpi_script();
        } else {
# ifdef G_OS_WIN32
            err = win32_run_mpi_sync(foreign_opt, prn);
# else
            err = lib_run_mpi_sync(foreign_opt, ptr, prn);
# endif
            if (err) {
                try_for_mpi_errmsg(prn);
            }
        }
        foreign_destroy();
        return err; /* handled */
    }
#endif /* HAVE_MPI */

    if (foreign_lang == LANG_R) {
#ifdef USE_RLIB
        if (gretl_use_Rlib()) {
# if FDEBUG
            fprintf(stderr, "HERE 1 LANG_R, using Rlib\n");
# endif
            err = run_R_lib(NULL, dset, foreign_opt, prn);
        } else {
# if FDEBUG
            fprintf(stderr, "HERE 2 LANG_R, using R binary\n");
# endif
            err = run_R_binary(NULL, dset, foreign_opt, prn);
        }
#else
# if FDEBUG
        fprintf(stderr, "HERE 3 LANG_R, using R binary\n");
# endif
        err = run_R_binary(NULL, dset, foreign_opt, prn);
#endif
        foreign_destroy();
        return err; /* handled */
    }

    err = write_gretl_foreign_script(NULL, foreign_lang,
                                     foreign_opt, dset, NULL);

    if (err) {
        delete_foreign_script(foreign_lang);
    } else {
#ifdef G_OS_WIN32
        err = win32_lib_run_other_sync(foreign_opt, prn);
#else
        err = lib_run_other_sync(foreign_opt, prn);
#endif
    }

    foreign_destroy();

    return err;
}

/**
 * execute_R_buffer:
 * @buf: buffer containing commands.
 * @dset: dataset struct.
 * @opt: may include %OPT_D to send data from gretl, %OPT_T
 * when called from gretl_edit.
 * @prn: struct for printing output.
 *
 * Returns: 0 on success, non-zero on error.
 */

int execute_R_buffer (const char *buf,
                      const DATASET *dset,
                      gretlopt opt,
                      PRN *prn)
{
    int err = 0;

    make_gretl_R_names();

#ifdef USE_RLIB
    if (!(opt & OPT_T) && gretl_use_Rlib()) {
        err = run_R_lib(buf, dset, opt, prn);
    } else {
        err = run_R_binary(buf, dset, opt, prn);
    }
#else
    err = run_R_binary(buf, dset, opt, prn);
#endif

    return err;
}

int check_R_depends (const char *pkgname, const char *deps,
                     PRN *prn)
{
    gchar *rbuf;
    int err;

    rbuf = get_checkdeps_buffer(deps);
    err = execute_R_buffer(rbuf, NULL, OPT_NONE, prn);

    if (err) {
        pprintf(prn, "%s: check for R dependencies failed\n", pkgname);
    } else {
        const char *pbuf = gretl_print_get_buffer(prn);

        //fprintf(stderr, "prn %p\n", (void *) prn);
        //fprintf(stderr, "prn buf: '%s'\n", pbuf);
        if (pbuf != NULL && strstr(pbuf, "not met")) {
            err = E_DEPENDS;
        }
    }

    g_free(rbuf);

    return err;
}
