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
#include "gretl_func.h"
#include "gretl_foreign.h"
#include "matrix_extra.h"
#include "gretl_typemap.h"

#ifdef HAVE_MPI
# include "gretl_mpi.h"
# include "gretl_xml.h"
#endif

#ifdef USE_RLIB
# include <Rinternals.h> /* for SEXP and friends */
#endif

#ifdef G_OS_WIN32
# include "gretl_win32.h"
#else
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
static gchar *gretl_Rsrc;
static gchar *gretl_Rout;
static gchar *gretl_Rmsg;

/* names of dotdir scripts to be passed to foreign programs */
static gchar *gretl_ox_script;
static gchar *gretl_octave_script;
static gchar *gretl_stata_script;
static gchar *gretl_python_script;
static gchar *gretl_julia_script;
#ifdef HAVE_MPI
static gchar *gretl_mpi_script;
#endif

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
	    err = E_EXTERNAL;
	} else {
	    foreign_lang = LANG_MPI;
	}
#else
	gretl_errmsg_set(_("MPI is not supported in this gretl build"));
	err = E_NOTIMP;
#endif
    } else {
	pprintf(prn, "%s: unknown language\n", lang);
	err = E_DATA;
    }

    return err;
}

#ifdef HAVE_MPI

enum {
    MPI_OPENMPI,
    MPI_MPICH,
    MPI_MSMPI
};

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

#else

int gretl_max_mpi_processes (void)
{
    return 0;
}

#endif /* HAVE_MPI */

static const gchar *get_ox_scriptname (void)
{
    if (gretl_ox_script == NULL) {
	const char *dotdir = gretl_dotdir();

	gretl_ox_script = g_strdup_printf("%sgretltmp.ox", dotdir);
    }

    return gretl_ox_script;
}

static const gchar *get_octave_scriptname (void)
{
    if (gretl_octave_script == NULL) {
	const char *dotdir = gretl_dotdir();

	gretl_octave_script = g_strdup_printf("%sgretltmp.m", dotdir);
    }

    return gretl_octave_script;
}

static const gchar *get_stata_scriptname (void)
{
    if (gretl_stata_script == NULL) {
	const char *dotdir = gretl_dotdir();

	gretl_stata_script = g_strdup_printf("%sgretltmp.do", dotdir);
    }

    return gretl_stata_script;
}

static const gchar *get_python_scriptname (void)
{
    if (gretl_python_script == NULL) {
	const char *dotdir = gretl_dotdir();

	gretl_python_script = g_strdup_printf("%sgretltmp.py", dotdir);
    }

    return gretl_python_script;
}

static const gchar *gretl_julia_scriptname (void)
{
    if (gretl_julia_script == NULL) {
	const char *dotdir = gretl_dotdir();

	gretl_julia_script = g_strdup_printf("%sgretltmp.jl", dotdir);
    }

    return gretl_julia_script;
}

#ifdef HAVE_MPI

static const gchar *get_mpi_scriptname (void)
{
    if (gretl_mpi_script == NULL) {
	const char *dotdir = gretl_dotdir();

	gretl_mpi_script = g_strdup_printf("%sgretltmp-mpi.inp", dotdir);
    }

    return gretl_mpi_script;
}

#endif

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
	gretl_Rsrc = g_strdup_printf("%sRsrc", ddir);
	gretl_Rout = g_strdup_printf("%sR.out", ddir);
	gretl_Rmsg = g_strdup_printf("%sR.msg", ddir);
	done = 1;
    }
}

#ifdef G_OS_WIN32

static char *get_rscript_path (void)
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

/* Windows variant */

static int lib_run_R_sync (gretlopt opt, PRN *prn)
{
    char *rscript = get_rscript_path();
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
    fprintf(stderr, "lib_run_R_sync: err = %d\n cmd='%s'\n", err, cmd);
#endif

    if (!(opt & OPT_Q)) {
	const gchar *outname;
	FILE *fp;

	outname = err ? gretl_Rmsg : gretl_Rout;
	fp = gretl_fopen(outname, "r");

	if (fp != NULL) {
	    char line[1024];

	    while (fgets(line, sizeof line, fp)) {
		pputs(prn, line);
	    }
	    fclose(fp);
	    gretl_remove(outname);
	}
    }

    g_free(cmd);

    return err;
}

/* Windows variant */

static int lib_run_other_sync (gretlopt opt, PRN *prn)
{
    const char *exe;
    const char *fname;
    gchar *cmd, *sout = NULL;
    int err;

    if (foreign_lang == LANG_OX) {
	exe = gretl_oxl_path();
	fname = get_ox_scriptname();
	cmd = g_strdup_printf("\"%s\" \"%s\"", exe, fname);
    } else if (foreign_lang == LANG_OCTAVE) {
	exe = gretl_octave_path();
	fname = get_octave_scriptname();
	cmd = g_strdup_printf("\"%s\" --silent \"%s\"", exe, fname);
    } else if (foreign_lang == LANG_STATA) {
	exe = gretl_stata_path();
	fname = get_stata_scriptname();
	cmd = g_strdup_printf("\"%s\" /q /e do \"%s\"", exe, fname);
    } else if (foreign_lang == LANG_PYTHON) {
	exe = gretl_python_path();
	fname = get_python_scriptname();
	cmd = g_strdup_printf("\"%s\" \"%s\"", exe, fname);
    } else if (foreign_lang == LANG_JULIA) {
	exe = gretl_julia_path();
	fname = gretl_julia_scriptname();
	cmd = g_strdup_printf("\"%s\" \"%s\"", exe, fname);
    } else {
	return 1;
    }

    err = gretl_win32_grab_output(cmd, gretl_workdir(), &sout);

    if (sout != NULL && *sout != '\0') {
	pputs(prn, sout);
    }

    if (!err && foreign_lang == LANG_STATA && !(opt & OPT_Q)) {
	/* output will be in log file, not stdout */
	do_stata_printout(prn);
    }

    g_free(sout);
    g_free(cmd);

    return err;
}

# ifdef HAVE_MPI

/* Windows: for now we'll not attempt to support anything
   other than "native" MS-MPI
*/

static int lib_run_mpi_sync (gretlopt opt, PRN *prn)
{
    const char *mpiexec = gretl_mpiexec();
    const char *hostfile = gretl_mpi_hosts();
    int np = 0;
    int err = 0;

    if (*hostfile == '\0') {
	hostfile = getenv("GRETL_MPI_HOSTS");
    }

    if (opt & OPT_N) {
	/* handle the number-of-processes option */
	np = get_optval_int(MPI, OPT_N, &err);
	if (!err && (np <= 0 || np > 9999999)) {
	    err = E_DATA;
	}
    }

    if (!err) {
	gchar *hostbit, *npbit, *rngbit, *qopt;
	gchar *cmd, *sout = NULL;

	if (!(opt & OPT_L) && hostfile != NULL && *hostfile != '\0') {
	    /* note: OPT_L corresponds to --local, meaning that we
	       should not use a hosts file even if one is present
	    */
	    hostbit = g_strdup_printf(" /machinefile \"%s\"", hostfile);
	} else {
	    hostbit = g_strdup("");
	    if (np == 0) {
		/* no hosts file: supply a default np value */
		np = gretl_n_processors();
	    }
	}

	if (np > 0) {
	    npbit = g_strdup_printf(" /np %d", np);
	} else {
	    npbit = g_strdup("");
	}

	if (opt & OPT_S) {
	    rngbit = g_strdup(" --single-rng");
	} else {
	    rngbit = g_strdup("");
	}

	if (opt & OPT_Q) {
	    qopt = g_strdup(" --quiet");
	} else {
	    qopt = g_strdup("");
	}

	cmd = g_strdup_printf("%s%s%s \"%sgretlmpi\"%s%s \"%s\"",
			      mpiexec, hostbit, npbit, gretl_bindir(), rngbit,
			      qopt, get_mpi_scriptname());

	if (opt & OPT_V) {
	    pputs(prn, "gretl mpi command:\n ");
	    pputs(prn, cmd);
	    pputc(prn, '\n');
	}

	err = gretl_win32_grab_output(cmd, gretl_workdir(), &sout);
	if (sout != NULL && *sout != '\0') {
	    pputs(prn, sout);
	}

	g_free(hostbit);
	g_free(npbit);
	g_free(rngbit);
	g_free(sout);
	g_free(cmd);
    }

    return err;
}

# endif /* HAVE_MPI (&& G_OS_WIN32) */

#else /* !G_OS_WIN32 */

static int lib_run_prog_sync (char **argv, gretlopt opt,
			      PRN *prn)
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
	pprintf(prn, "%s\n", gerr->message);
	g_error_free(gerr);
	err = 1;
    } else if (status != 0) {
	pprintf(prn, "%s exited with status %d", argv[0], status);
	if (sout != NULL && *sout != '\0') {
	    pputs(prn, "stdout:\n");
	    pputs(prn, sout);
	    pputc(prn, '\n');
	}
	if (errout != NULL && *errout != '\0') {
	    pputs(prn, "\nstderr:\n");
	    pputs(prn, errout);
	    pputc(prn, '\n');
	}
	err = 1;
    } else if (sout != NULL) {
	if (!(opt & OPT_Q)) {
	    /* with OPT_Q, don't print non-error output */
	    if (foreign_lang == LANG_STATA) {
		do_stata_printout(prn);
	    } else {
		pputs(prn, sout);
		pputc(prn, '\n');
	    }
	}
	if (opt & OPT_V) {
	    /* also print stderr output, if any */
	    if (errout != NULL && *errout != '\0') {
		pputs(prn, "\nstderr:\n");
		pputs(prn, errout);
		pputc(prn, '\n');
	    }
	}
    } else {
	pprintf(prn, "%s: %s\n", argv[0], "Got no output");
	err = 1;
    }

    g_free(sout);
    g_free(errout);

    return err;
}

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

    return lib_run_prog_sync(argv, opt, prn);
}

static int lib_run_other_sync (gretlopt opt, PRN *prn)
{
    char *argv[6];
    int err;

    if (foreign_lang == LANG_OX) {
	argv[0] = (char *) gretl_oxl_path();
	argv[1] = (char *) get_ox_scriptname();
	argv[2] = NULL;
    } else if (foreign_lang == LANG_OCTAVE) {
	argv[0] = (char *) gretl_octave_path();
	argv[1] = "--silent";
	argv[2] = (char *) get_octave_scriptname();
	argv[3] = NULL;
    } else if (foreign_lang == LANG_PYTHON) {
	argv[0] = (char *) gretl_python_path();
	argv[1] = (char *) get_python_scriptname();
	argv[2] = NULL;
    } else if (foreign_lang == LANG_JULIA) {
	argv[0] = (char *) gretl_julia_path();
	argv[1] = (char *) gretl_julia_scriptname();
	argv[2] = NULL;
    } else if (foreign_lang == LANG_STATA) {
	argv[0] = (char *) gretl_stata_path();
	argv[1] = "-q";
	argv[2] = "-b";
	argv[3] = "do";
	argv[4] = (char *) get_stata_scriptname();
	argv[5] = NULL;
    }

    err = lib_run_prog_sync(argv, opt, prn);

    return err;
}

#ifdef HAVE_MPI

static void print_mpi_command (char **argv, PRN *prn)
{
    int i;

    pputs(prn, "gretl mpi command:\n ");
    for (i=0; argv[i] != NULL; i++) {
	pprintf(prn, "%s ", argv[i]);
    }
    pputc(prn, '\n');
}

/* The following should probably be redundant on Linux but may
   be needed for the OS X package, where the gretl bin
   directory may not be in PATH.
*/

static gchar *gretl_mpi_binary (void)
{
    gchar *tmp = g_strdup(gretl_home());
    gchar *p = strstr(tmp, "/share/gretl");
    gchar *ret;

    if (p != NULL) {
	*p = '\0';
	ret = g_strdup_printf("%s/bin/gretlmpi", tmp);
    } else {
	ret = g_strdup("gretlmpi");
    }

    g_free(tmp);

    return ret;
}

static int lib_run_mpi_sync (gretlopt opt, PRN *prn)
{
    const char *hostfile = gretl_mpi_hosts();
    char npnum[8] = {0};
    int err = 0;

    if (*hostfile == '\0') {
	hostfile = getenv("GRETL_MPI_HOSTS");
    }

    if (opt & OPT_N) {
	/* handle the number-of-processes option */
	int np = get_optval_int(MPI, OPT_N, &err);

	if (!err && (np <= 0 || np > 9999999)) {
	    err = E_DATA;
	}
	if (!err) {
	    sprintf(npnum, "%d", np);
	}
    }

    if (!err) {
	const char *mpiexec = gretl_mpiexec();
	gchar *mpiprog = gretl_mpi_binary();
	const char *hostsopt = NULL;
	char *argv[11];
	int npmax, i = 0;

	npmax = gretl_n_processors();

	if (!(opt & OPT_L) && hostfile != NULL && *hostfile != '\0') {
	    hostsopt = (mpi_variant == MPI_MPICH)? "-machinefile" :
		"--hostfile";
	} else if (*npnum == '\0') {
	    /* no hosts file: supply a default np value */
	    sprintf(npnum, "%d", npmax);
	}

	argv[i++] = (char *) mpiexec;
	if (hostsopt != NULL) {
	    argv[i++] = (char *) hostsopt;
	    argv[i++] = (char *) hostfile;
	}
	if (*npnum != '\0') {
	    argv[i++] = "-np";
	    argv[i++] = npnum;
	    if (mpi_variant == MPI_OPENMPI && atoi(npnum) > npmax/2) {
		argv[i++] = "--oversubscribe";
	    }
	}
	argv[i++] = mpiprog;
	if (opt & OPT_S) {
	    argv[i++] = "--single-rng";
	}
	if (opt & OPT_Q) {
	    argv[i++] = "--quiet";
	}
	argv[i++] = (char *) get_mpi_scriptname();
	argv[i] = NULL;

	if (opt & OPT_V) {
	    print_mpi_command(argv, prn);
	}

	err = lib_run_prog_sync(argv, opt, prn);
	g_free(mpiprog);
    }

    return err;
}

#endif /* HAVE_MPI */

#endif /* switch on MS Windows or not */

static FILE *write_open_dotfile (const char *fname)
{
    FILE *fp;
    gchar *path;

    path = g_strdup_printf("%s%s", gretl_dotdir(), fname);
    fp = gretl_fopen(path, "w");
    g_free(path);

    return fp;
}

static int write_ox_io_file (void)
{
    static int written;

    if (!written) {
	FILE *fp = write_open_dotfile("gretl_io.ox");

	if (fp == NULL) {
	    return E_FOPEN;
	} else {
	    const char *ddir = get_export_dotdir();

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

	    fclose(fp);
	    written = 1;
	}
    }

    return 0;
}

static int real_write_octave_io_file (void)
{
    FILE *fp = write_open_dotfile("gretl_io.m");

    if (fp == NULL) {
	return E_FOPEN;
    } else {
	const char *ddir = get_export_dotdir();

	fputs("# not a 'function file' as such\n1;\n", fp);
	fputs("function dotdir = gretl_dotdir()\n", fp);
	fprintf(fp, "  dotdir = \"%s\";\n", ddir);
	fputs("endfunction\n\n", fp);

	fputs("function gretl_export(X, str, autodot=1)\n", fp);
	fputs("  if (autodot)\n", fp);
	fputs("    dname = gretl_dotdir();\n", fp);
	fputs("    fd = fopen(strcat(dname, str), \"w\");\n", fp);
	fputs("  else\n", fp);
	fputs("    fd = fopen(str, \"w\");\n", fp);
	fputs("  endif\n", fp);
	fputs("  fprintf(fd, \"%d %d\\n\", size(X));\n", fp);
	fputs("  c = columns(X);\n", fp);
	fputs("  fs = strcat(strrep(sprintf(\"%d \", ones(1, c)), \"1\", \"%.15g\"), \"\\n\");",
	      fp);
	fputc('\n', fp);
	fputs("  fprintf(fd, fs, X');\n", fp);
	fputs("  fclose(fd);\n", fp);
	fputs("endfunction\n\n", fp);

	fputs("function A = gretl_loadmat(str, autodot=1)\n", fp);
	fputs("  if (autodot)\n", fp);
	fputs("    dname = gretl_dotdir();\n", fp);
	fputs("    fd = fopen(strcat(dname, str), \"r\");\n", fp);
	fputs("  else\n", fp);
	fputs("    fd = fopen(str, \"r\");\n", fp);
	fputs("  endif\n", fp);
	fputs("  [r,c] = fscanf(fd, \"%d %d\", \"C\");\n", fp);
	fputs("  A = reshape(fscanf(fd, \"%g\", r*c),c,r)';\n", fp);
	fputs("  fclose(fd);\n", fp);
	fputs("endfunction\n\n", fp);

	fclose(fp);
    }

    return 0;
}

static int write_octave_io_file (void)
{
    static int written;
    int err = 0;

    if (!written) {
	err = real_write_octave_io_file();
	if (!err) {
	    written = 1;
	}
    }

    return err;
}

static int write_python_io_file (void)
{
    static int written;

    if (!written) {
	FILE *fp = write_open_dotfile("gretl_io.py");

	if (fp == NULL) {
	    return E_FOPEN;
	} else {
	    const char *ddir = get_export_dotdir();

	    fprintf(fp, "gretl_dotdir = \"%s\"\n\n", ddir);
	    /* export matrix for reading by gretl */
	    fputs("def gretl_export(X, fname, autodot=1):\n", fp);
	    fputs("  binwrite = 0\n", fp);
	    fputs("  if fname[-4:] == '.bin':\n", fp);
	    fputs("    binwrite = 1\n", fp);
	    fputs("    from numpy import asmatrix, asarray\n", fp);
	    fputs("    from struct import pack\n", fp);
	    fputs("  else:\n", fp);
	    fputs("    from numpy import asmatrix, savetxt\n", fp);
	    fputs("  M = asmatrix(X)\n", fp);
	    fputs("  r, c = M.shape\n", fp);
	    fputs("  if autodot:\n", fp);
            fputs("    fname = gretl_dotdir + fname\n", fp);
	    fputs("  if binwrite:\n", fp);
	    fputs("    from sys import byteorder\n", fp);
	    fputs("    f = open(fname, 'wb')\n", fp);
	    fputs("    f.write(b'gretl_binary_matrix')\n", fp);
	    fputs("    f.write(pack('<i', r))\n", fp);
	    fputs("    f.write(pack('<i', c))\n", fp);
	    fputs("    if byteorder == 'big':\n", fp);
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
	    fputs("def gretl_loadmat(fname, autodot=1):\n", fp);
	    fputs("  if autodot:\n", fp);
	    fputs("    fname = gretl_dotdir + fname\n", fp);
	    fputs("  if fname[-4:] == '.bin':\n", fp);
	    fputs("    from numpy import ndarray, asmatrix\n", fp);
	    fputs("    from struct import unpack\n", fp);
	    fputs("    f = open(fname, 'rb')\n", fp);
	    fputs("    buf = f.read(19)\n", fp);
	    fputs("    if buf != b'gretl_binary_matrix':\n", fp);
	    fputs("      raise ValueError('Not a gretl binary matrix')\n", fp);
	    fputs("    r = unpack('<i', f.read(4))[0]\n", fp);
	    fputs("    c = unpack('<i', f.read(4))[0]\n", fp);
	    fputs("    M = ndarray(shape=(r,c), dtype=float, order='F')\n", fp);
	    fputs("    for j in range(0, c):\n", fp);
	    fputs("      for i in range(0, r):\n", fp);
	    fputs("        M[i,j] = unpack('<d', f.read(8))[0]\n", fp);
	    fputs("    f.close()\n", fp);
	    fputs("    M = asmatrix(M)\n", fp);
	    fputs("  else:\n", fp);
	    fputs("    from numpy import loadtxt\n", fp);
	    fputs("    M = loadtxt(fname, skiprows=1)\n", fp);
	    fputs("  return M\n\n", fp);

	    fclose(fp);
	    written = 1;
	}
    }

    return 0;
}

static int write_julia_io_file (void)
{
    static int written;

    if (!written) {
	FILE *fp = write_open_dotfile("gretl_io.jl");

	if (fp == NULL) {
	    return E_FOPEN;
	} else {
	    const char *ddir = get_export_dotdir();

	    fprintf(fp, "gretl_dotdir = \"%s\"\n\n", ddir);
	    /* Julia 1.0 requires more library-loading */
	    fputs("v1 = VERSION > v\"0.6.9\"\n", fp);
	    fputs("if v1\n", fp);
	    fputs("  using Printf\n", fp);
	    fputs("  using DelimitedFiles\n", fp);
	    fputs("end\n\n", fp);
	    fputs("function gretl_export(M, fname, autodot=1)\n", fp);
	    fputs("  r,c = size(M)\n", fp);
	    fputs("  if autodot != 0\n", fp);
	    fputs("    fname = gretl_dotdir * fname\n", fp);
	    fputs("  end\n", fp);
	    fputs("  f = open(fname, \"w\")\n", fp);
	    fputs("  if v1\n", fp);
	    fputs("    n = lastindex(fname)\n", fp);
	    fputs("  else\n", fp);
	    fputs("    n = endof(fname)\n", fp);
	    fputs("  end\n", fp);
	    fputs("  if fname[n-3:n] == \".bin\"\n", fp);
	    fputs("    # binary mode\n", fp);
	    fputs("    write(f, b\"gretl_binary_matrix\")\n", fp);
	    fputs("    if ENDIAN_BOM == 0x01020304\n", fp);
	    fputs("      # host is big-endian\n", fp);
	    fputs("      write(f, htol(Int32(r)))\n", fp);
	    fputs("      write(f, htol(Int32(c)))\n", fp);
	    fputs("      for j=1:c\n", fp);
	    fputs("        for i=1:r\n", fp);
	    fputs("          write(f, htol(Float64(M[i,j])))\n", fp);
	    fputs("        end\n", fp);
	    fputs("      end\n", fp);
	    fputs("    else\n", fp);
	    fputs("      write(f, Int32(r))\n", fp);
	    fputs("      write(f, Int32(c))\n", fp);
	    fputs("      write(f, M)\n", fp);
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
	    fputs("  if autodot != 0\n", fp);
	    fputs("    fname = gretl_dotdir * fname\n", fp);
	    fputs("  end\n", fp);
	    fputs("  if v1\n", fp);
	    fputs("    n = lastindex(fname)\n", fp);
	    fputs("  else\n", fp);
	    fputs("    n = endof(fname)\n", fp);
	    fputs("  end\n", fp);
	    fputs("  if fname[n-3:n] == \".bin\"\n", fp);
	    fputs("    # binary mode\n", fp);
	    fputs("    f = open(fname, \"r\")\n", fp);
	    fputs("    hdr = read(f, UInt8, 19)\n", fp);
	    fputs("    if hdr != b\"gretl_binary_matrix\"\n", fp);
	    fputs("      error(\"Not a gretl binary matrix\")\n", fp);
	    fputs("    end\n", fp);
	    fputs("    if ENDIAN_BOM == 0x01020304\n", fp);
	    fputs("      # host is big-endian\n", fp);
	    fputs("      r = ltoh(read(f, Int32))\n", fp);
	    fputs("      c = ltoh(read(f, Int32))\n", fp);
	    fputs("      M = Array{Float64, 2}(r, c)\n", fp);
	    fputs("      for j=1:c\n", fp);
	    fputs("        for i=1:r\n", fp);
	    fputs("          M[i,j] = ltoh(read(f, Float64))\n", fp);
	    fputs("        end\n", fp);
	    fputs("      end\n", fp);
	    fputs("    else\n", fp);
	    fputs("      r = read(f, Int32)\n", fp);
	    fputs("      c = read(f, Int32)\n", fp);
	    fputs("      M = Array{Float64, 2}(r, c)\n", fp);
	    fputs("      for j=1:c\n", fp);
	    fputs("        for i=1:r\n", fp);
	    fputs("          M[i,j] = read(f, Float64)\n", fp);
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

	    fclose(fp);
	    written = 1;
	}
    }

    return 0;
}

static int write_stata_io_file (void)
{
    static int written;

    if (!written) {
	FILE *fp = write_open_dotfile("gretl_export.ado");

	if (fp == NULL) {
	    return E_FOPEN;
	} else {
	    const char *ddir = get_export_dotdir();

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

	    fclose(fp);
	    written = 1;
	}
    }

    return 0;
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
	}
    }

    bufgets_finalize(buf);
}

/**
 * write_gretl_ox_script:
 * @buf: text buffer containing Ox code.
 * @opt: should contain %OPT_G for use from GUI.
 * @pfname: location to receive name of file written, or %NULL.
 *
 * Writes the content of @buf into a file in the gretl user's
 * "dotdir".
 *
 * Returns: 0 on success, non-zero on error.
 */

int write_gretl_ox_script (const char *buf, gretlopt opt,
			   const char **pfname)
{
    const gchar *fname = get_ox_scriptname();
    FILE *fp = gretl_fopen(fname, "w");

    write_ox_io_file();

    if (fp == NULL) {
	return E_FOPEN;
    } else {
	if (buf != NULL) {
	    /* pass on the material supplied in the 'buf' argument */
	    put_foreign_buffer(buf, fp);
	} else {
	    /* put out the stored 'foreign' lines */
	    put_foreign_lines(fp);
	}
	fclose(fp);
	if (pfname != NULL) {
	    *pfname = fname;
	}
    }

    return 0;
}

/**
 * write_gretl_python_script:
 * @buf: text buffer containing Python code.
 * @opt: should contain %OPT_G for use from GUI.
 * @pfname: location to receive name of file written, or %NULL.
 *
 * Writes the content of @buf into a file in the gretl user's
 * "dotdir".
 *
 * Returns: 0 on success, non-zero on error.
 */

int write_gretl_python_script (const char *buf, gretlopt opt,
			       const char **pfname)
{
    const gchar *fname = get_python_scriptname();
    FILE *fp = gretl_fopen(fname, "w");

    write_python_io_file();

    if (fp == NULL) {
	return E_FOPEN;
    } else {
	add_gretl_include(LANG_PYTHON, opt, fp);
	if (buf != NULL) {
	    /* pass on the material supplied in the 'buf' argument */
	    put_foreign_buffer(buf, fp);
	} else {
	    /* put out the stored 'foreign' lines */
	    put_foreign_lines(fp);
	}
	fclose(fp);
	if (pfname != NULL) {
	    *pfname = fname;
	}
    }

    return 0;
}

/**
 * write_gretl_julia_script:
 * @buf: text buffer containing Julia code.
 * @opt: should contain %OPT_G for use from GUI.
 * @pfname: location to receive name of file written, or %NULL.
 *
 * Writes the content of @buf into a file in the gretl user's
 * "dotdir".
 *
 * Returns: 0 on success, non-zero on error.
 */

int write_gretl_julia_script (const char *buf, gretlopt opt,
			      const char **pfname)
{
    /* FIXME more to be done here! */
    const gchar *fname = gretl_julia_scriptname();
    FILE *fp = gretl_fopen(fname, "w");

    write_julia_io_file();

    if (fp == NULL) {
	return E_FOPEN;
    } else {
	add_gretl_include(LANG_JULIA, opt, fp);
	if (buf != NULL) {
	    /* pass on the material supplied in the 'buf' argument */
	    put_foreign_buffer(buf, fp);
	} else {
	    /* put out the stored 'foreign' lines */
	    put_foreign_lines(fp);
	}
	fclose(fp);
	if (pfname != NULL) {
	    *pfname = fname;
	}
    }

    return 0;
}

static int no_data_check (const DATASET *dset)
{
    if (dset == NULL || dset->n == 0 || dset->v == 0) {
	return E_NODATA;
    } else {
	return 0;
    }
}

static int *get_send_data_list (const DATASET *dset, int ci, int *err)
{
    const char *lname = get_optval_string(FOREIGN, OPT_D);
    int *list = NULL;

    if (lname != NULL) {
	list = get_list_by_name(lname);
	if (list == NULL) {
	    *err = E_DATA;
	} else {
	    int i;

	    for (i=1; i<=list[0] && !*err; i++) {
		if (list[i] < 0 || list[i] >= dset->v) {
		    *err = E_DATA;
		}
	    }
	}
    }

    return list;
}

#ifdef HAVE_MPI

static int mpi_send_data_setup (const DATASET *dset, FILE *fp)
{
    const char *dotdir = gretl_dotdir();
    int *list = NULL;
    size_t datasize;
    int nvars;
    gchar *fname;
    int err;

    err = no_data_check(dset);
    if (err) {
	return err;
    }

    list = get_send_data_list(dset, MPI, &err);
    if (list != NULL) {
	nvars = list[0];
    } else {
	nvars = dset->v;
    }

    datasize = dset->n * nvars;

    if (datasize > 10000) {
	/* write "big" data as binary? */
	fname = g_strdup_printf("%smpi-data.gdtb", dotdir);
    } else {
	fname = g_strdup_printf("%smpi-data.gdt", dotdir);
    }

    err = gretl_write_gdt(fname, list, dset, OPT_NONE, 0);

    if (!err) {
	/* here we're writing into the file to be run
	   by gretlmpi */
	fprintf(fp, "open \"%s\"\n", fname);
    }

    g_free(fname);

    return err;
}

static int mpi_send_funcs_setup (FILE *fp)
{
    gchar *fname;
    int err;

    fname = g_strdup_printf("%smpi-funcs-tmp.xml", gretl_dotdir());
    err = write_loaded_functions_file(fname, 1);

    if (!err) {
	fprintf(fp, "include \"%s\"\n", fname);
    }

    g_free(fname);

    return err;
}

static int write_gretl_mpi_script (gretlopt opt, const DATASET *dset)
{
    const gchar *fname = get_mpi_scriptname();
    FILE *fp = gretl_fopen(fname, "w");
    int err = 0;

    if (fp == NULL) {
	return E_FOPEN;
    }

    if (opt & OPT_D) {
	/* honor the --send-data option */
	err = mpi_send_data_setup(dset, fp);
    }

    if (opt & OPT_F) {
	/* honor the --send-functions option */
	if (n_user_functions() > 0) {
	    err = mpi_send_funcs_setup(fp);
	}
    }

#if defined(_OPENMP)
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
	fclose(fp);
    }

    return err;
}

#endif /* HAVE_MPI */

static int write_data_for_stata (const DATASET *dset,
				 FILE *fp)
{
    int *list = NULL;
    gchar *sdata = NULL;
    char save_na[8];
    int err;

    err = no_data_check(dset);
    if (err) {
	return err;
    }

    list = get_send_data_list(dset, FOREIGN, &err);

    if (!err) {
	*save_na = '\0';
	strncat(save_na, get_csv_na_write_string(), 7);
	set_csv_na_write_string(".");
	sdata = g_strdup_printf("%sstata.csv", gretl_dotdir());
	err = write_data(sdata, list, dset, OPT_C, NULL);
	set_csv_na_write_string(save_na);
	g_free(sdata);
    }

    if (err) {
	gretl_errmsg_sprintf("write_data_for_stata: failed with err = %d\n", err);
    } else {
	fputs("* load data from gretl\n", fp);
	sdata = g_strdup_printf("%sstata.csv", get_export_dotdir());
	fprintf(fp, "insheet using \"%s\"\n", sdata);
	g_free(sdata);
    }

    return err;
}

int write_gretl_stata_script (const char *buf, gretlopt opt,
			      const DATASET *dset,
			      const char **pfname)
{
    const gchar *fname = get_stata_scriptname();
    FILE *fp = gretl_fopen(fname, "w");

    write_stata_io_file();

    if (fp == NULL) {
	return E_FOPEN;
    } else {
	int err;

	/* source the I-O functions */
	add_gretl_include(LANG_STATA, opt, fp);
	if (opt & OPT_D) {
	    /* --send-data */
	    err = write_data_for_stata(dset, fp);
	    if (err) {
		fclose(fp);
		return err;
	    }
	}
	if (buf != NULL) {
	    /* pass on the material supplied in the 'buf' argument */
	    put_foreign_buffer(buf, fp);
	} else {
	    /* put out the stored 'foreign' lines */
	    put_foreign_lines(fp);
	}
	fclose(fp);
	if (pfname != NULL) {
	    *pfname = fname;
	}
    }

    return 0;
}

/* write out current dataset as an octave matrix, and, if this succeeds,
   write appropriate octave commands to @fp to source the data
*/

static int write_data_for_octave (const DATASET *dset,
				  FILE *fp)
{
    int *list = NULL;
    gchar *mdata;
    int err;

    err = no_data_check(dset);
    if (err) {
	return err;
    }

    list = get_send_data_list(dset, FOREIGN, &err);

    if (!err) {
	mdata = g_strdup_printf("%smdata.tmp", gretl_dotdir());
	err = write_data(mdata, list, dset, OPT_M, NULL);
	g_free(mdata);
    }

    if (err) {
	gretl_errmsg_sprintf("write_data_for_octave: failed with err = %d\n", err);
    } else {
	fputs("% load data from gretl\n", fp);
	fprintf(fp, "load '%smdata.tmp'\n", get_export_dotdir());
    }

    g_free(mdata);

    return err;
}

int write_gretl_octave_script (const char *buf, gretlopt opt,
			       const DATASET *dset,
			       const char **pfname)
{
    const gchar *fname = get_octave_scriptname();
    FILE *fp = gretl_fopen(fname, "w");

    write_octave_io_file();

    if (fp == NULL) {
	return E_FOPEN;
    } else {
	int err;

	/* source the I-O functions */
	add_gretl_include(LANG_OCTAVE, opt, fp);
	if (opt & OPT_D) {
	    /* --send-data */
	    err = write_data_for_octave(dset, fp);
	    if (err) {
		fclose(fp);
		return err;
	    }
	}
	if (buf != NULL) {
	    /* pass on the material supplied in the 'buf' argument */
	    put_foreign_buffer(buf, fp);
	} else {
	    /* put out the stored 'foreign' lines */
	    put_foreign_lines(fp);
	}
	fclose(fp);
	if (pfname != NULL) {
	    *pfname = fname;
	}
    }

    return 0;
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
    gchar *Rdata;
    int ts, err;

    err = no_data_check(dset);
    if (err) {
	return err;
    }

    /* FIXME: can R's "ts" handle daily data, weekly data, etc.? */
    ts = annual_data(dset) || quarterly_or_monthly(dset);

    list = get_send_data_list(dset, FOREIGN, &err);

    if (!err) {
	Rdata = g_strdup_printf("%sRdata.tmp", gretl_dotdir());
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
	gchar *tmp = g_strdup_printf("%sRcoded.mat", gretl_dotdir());
	int mwerr;

	/* ensure we don't load a stale file */
	gretl_remove(tmp);
	g_free(tmp);
	mwerr = gretl_matrix_write_to_file(coded, "Rcoded.mat", 1);
	if (mwerr) {
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

	ntodate(datestr, dset->t1, dset);
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

static void write_R_io_funcs (FILE *fp)
{
    const char *ddir = get_export_dotdir();

    fprintf(fp, "gretl.dotdir <- \"%s\"\n", ddir);

    fputs("gretl.export <- function(x, sx) {\n", fp);
    fprintf(fp, "  prefix <- \"%s\"\n", ddir);
    fputs("  objname <- as.character(substitute(x))\n", fp);
    fputs("  if (missing(sx)) {\n", fp);
    fputs("    sx <- objname\n", fp);
    fputs("  }\n", fp);
    fputs("  if (is.ts(x)) {\n", fp);
    fputs("    fname <- paste(prefix, sx, \".csv\", sep=\"\")\n", fp);
    fputs("    dfx <- data.frame(x)\n", fp);
    fputs("    if (ncol(dfx) == 1) {\n", fp);
    fputs("      colnames(dfx) <- sx;\n", fp);
    fputs("    }\n", fp);
    fputs("    write.csv(dfx, file=fname, row.names=F)\n", fp);
    fputs("    gretlmsg <- paste(\"wrote CSV data\", fname, \"\\n\")\n", fp);
    fputs("  } else if (is.data.frame(x)) {\n", fp);
    fputs("    fname <- paste(prefix, sx, \".csv\", sep=\"\")\n", fp);
    fputs("    write.csv(x, file=fname, row.names=F)\n", fp);
    fputs("    gretlmsg <- paste(\"wrote CSV data\", fname, \"\\n\")\n", fp);
    fputs("  } else if (is.matrix(x)) {\n", fp);
    fputs("    fname <- paste(prefix, sx, \".mat\", sep=\"\")\n", fp);
    fputs("    write(dim(x), fname)\n", fp);
    fputs("    write(format(t(x), digits=15), file=fname, ncolumns=ncol(x), append=TRUE)\n", fp);
    fputs("    gretlmsg <- paste(\"wrote matrix\", fname, \"\\n\")\n", fp);
    fputs("  } else {\n", fp);
    fputs("    gretlmsg <- paste(\"gretl.export: don't know how to write object\", objname, "
	  " \"(try as.matrix?)\\n\")\n", fp);
    fputs("  }\n", fp);
    fputs("  cat(gretlmsg)\n", fp);
    fputs("}\n", fp);

    fputs("gretl.loadmat <- function(mname) {\n", fp);
    fprintf(fp, "  prefix <- \"%s\"\n", ddir);
    fputs("  fname <- paste(prefix, mname, sep=\"\")\n", fp);
    fputs("  m <- as.matrix(read.table(fname, skip=1))\n", fp);
    fputs("  return(m)\n", fp);
    fputs("}\n", fp);
}

/* basic content which can either go into gretl.Rprofile or into
   Rsrc for sourcing */

static void put_R_startup_content (FILE *fp)
{
    fputs("vnum <- as.double(R.version$major) + (as.double(R.version$minor) / 10.0)\n",
	  fp);
    fputs("if (vnum > 2.41) library(utils)\n", fp);
    fputs("library(stats)\n", fp);
    fputs("if (vnum <= 1.89) library(ts)\n", fp);
    write_R_io_funcs(fp);
}

/* Set up a gretl-specific R profile, and put notice of its existence
   into the environment. Used when exec'ing the R binary (only) */

static int write_gretl_R_profile (gretlopt opt)
{
    FILE *fp;
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

    fp = gretl_fopen(gretl_Rprofile, "w");

    if (fp == NULL) {
	err = E_FOPEN;
    } else {
	put_R_startup_content(fp);
	fprintf(fp, "source(\"%s\", %s = TRUE)\n",
		gretl_Rsrc, (opt & OPT_V)? "echo" : "print.eval");
	fclose(fp);
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
	    fprintf(fp, "errout <- file(\"%s\", open=\"wt\")\n", gretl_Rmsg);
	    fputs("sink(errout, type=\"message\")\n", fp);
	    sunk = 1;
	}
#endif

	if (opt & OPT_L) {
	    /* we're using the R shared library */
	    static int startup_done;

	    if (!startup_done) {
#if FDEBUG
		fprintf(stderr, "Rlib: writing 'startup' material\n");
#endif
		put_R_startup_content(fp);
		startup_done = 1;
	    }
	    fprintf(fp, "sink(\"%s\", type=\"output\")\n", gretl_Rout);
	    if (!(opt & OPT_I)) {
		fprintf(fp, "errout <- file(\"%s\", open=\"wt\")\n", gretl_Rmsg);
		fputs("sink(errout, type=\"message\")\n", fp);
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

static void delete_ox_script (void)
{
    if (gretl_ox_script != NULL) {
	gretl_remove(gretl_ox_script);
    }
}

static void delete_octave_script (void)
{
    if (gretl_octave_script != NULL) {
	gretl_remove(gretl_octave_script);
    }
}

static void delete_stata_script (void)
{
    if (gretl_stata_script != NULL) {
	gretl_remove(gretl_stata_script);
    }
}

static void delete_python_script (void)
{
    if (gretl_python_script != NULL) {
	gretl_remove(gretl_python_script);
    }
}

static void delete_julia_script (void)
{
    if (gretl_julia_script != NULL) {
	gretl_remove(gretl_julia_script);
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

SEXP VR_GlobalEnv;
SEXP VR_NilValue;
SEXP VR_UnboundValue;

/* renamed, pointerized versions of the R functions we need */

static double *(*R_REAL) (SEXP);
static const char *(*R_STRING) (SEXP);
static SEXP *(*R_STRING_PTR) (SEXP);

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

static Rboolean (*R_isMatrix) (SEXP);
static Rboolean (*R_isLogical) (SEXP);
static Rboolean (*R_isInteger) (SEXP);
static Rboolean (*R_isReal) (SEXP);
static Rboolean (*R_isString) (SEXP);

static int (*R_initEmbeddedR) (int, char **);
static int (*R_ncols) (SEXP);
static int (*R_nrows) (SEXP);
static int (*R_TYPEOF) (SEXP);

static void (*R_endEmbeddedR) (int);
static void (*R_unprotect) (int);
static void (*R_PrintValue) (SEXP);
static void (*R_SET_TYPEOF) (SEXP, int);
static void (*R_SET_TAG) (SEXP, SEXP);

static int *(*R_LOGICAL) (SEXP);

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
    R_STRING        = dlget(Rhandle, "R_CHAR", &err);
    R_STRING_PTR    = dlget(Rhandle, "STRING_PTR", &err);
    R_allocList     = dlget(Rhandle, "Rf_allocList", &err);
    R_allocMatrix   = dlget(Rhandle, "Rf_allocMatrix", &err);
    R_allocVector   = dlget(Rhandle, "Rf_allocVector", &err);
    R_endEmbeddedR  = dlget(Rhandle, "Rf_endEmbeddedR", &err);
    R_findFun       = dlget(Rhandle, "Rf_findFun", &err);
    R_findVar       = dlget(Rhandle, "Rf_findVar", &err);
    R_initEmbeddedR = dlget(Rhandle, "Rf_initEmbeddedR", &err);
    R_install       = dlget(Rhandle, "Rf_install", &err);
    R_isMatrix      = dlget(Rhandle, "Rf_isMatrix", &err);
    R_isLogical     = dlget(Rhandle, "Rf_isLogical", &err);
    R_isInteger     = dlget(Rhandle, "Rf_isInteger", &err);
    R_isReal        = dlget(Rhandle, "Rf_isReal", &err);
    R_isString      = dlget(Rhandle, "Rf_isString", &err);
    R_mkString      = dlget(Rhandle, "Rf_mkString", &err);
    R_ncols         = dlget(Rhandle, "Rf_ncols", &err);
    R_nrows         = dlget(Rhandle, "Rf_nrows", &err);
    R_PrintValue    = dlget(Rhandle, "Rf_PrintValue", &err);
    R_protect       = dlget(Rhandle, "Rf_protect", &err);
    R_ScalarReal    = dlget(Rhandle, "Rf_ScalarReal", &err);
    R_unprotect     = dlget(Rhandle, "Rf_unprotect", &err);
    R_catch         = dlget(Rhandle, "R_tryEval", &err);
    R_SETCAR        = dlget(Rhandle, "SETCAR", &err);
    R_SET_TYPEOF    = dlget(Rhandle, "SET_TYPEOF", &err);
    R_TYPEOF        = dlget(Rhandle, "TYPEOF", &err);
    R_SET_TAG       = dlget(Rhandle, "SET_TAG", &err);
    R_LOGICAL       = dlget(Rhandle, "LOGICAL", &err);

#ifdef WIN32
    R_get_HOME = dlget(Rhandle, "get_R_HOME", &err);
#endif

    if (!err) {
	PR_GlobalEnv    = (SEXP *) dlget(Rhandle, "R_GlobalEnv", &err);
	PR_NilValue     = (SEXP *) dlget(Rhandle, "R_NilValue", &err);
	PR_UnboundValue = (SEXP *) dlget(Rhandle, "R_UnboundValue", &err);
    }

    if (err) {
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

#ifdef WIN32

/* try to ensure that the directory holding the R
   DLL is in PATH
*/

static void set_path_for_Rlib (const char *Rhome)
{
#ifdef _WIN64
    const char *arch = "x64";
#else
    const char *arch = "i386";
#endif
    char *oldpath = getenv("PATH");
    gchar *Rpath;

    Rpath = g_strdup_printf("%s\\bin\\%s", Rhome, arch);

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

#else /* !WIN32 */

/* non-Windows: attempt to remedy the absence of the
   R_HOME environment variable. We try to infer the
   required directory from take the path to libR.so and
   push it into the environment.
*/

static void try_set_R_home (void)
{
    const char *libpath = gretl_rlib_path();
    char *s, *tmp;

    tmp = gretl_strdup(libpath);
    s = strstr(tmp, "/lib/libR");
    if (s != NULL) {
	*s = '\0';
	gretl_setenv("R_HOME", tmp);
    }
    free(tmp);
}

#endif /* WIN32 or not */

/* Initialize the R library for use with gretl.  Note that we only
   need do this once per gretl session.  We need to check that the
   environment is set to R's liking first, otherwise initialization
   will fail -- and will abort gretl too!
*/

static int gretl_Rlib_init (void)
{
    char *Rhome;
    int err = 0;

#if FDEBUG
    fprintf(stderr, "gretl_Rlib_init: starting\n");
#endif

#ifndef WIN32
    Rhome = getenv("R_HOME");
    if (Rhome == NULL) {
	try_set_R_home();
    }
#endif

    err = load_R_symbols();
    if (err) {
	fprintf(stderr, "gretl_Rlib_init: failed to load R functions\n");
	goto bailout;
    }

#ifdef WIN32
    Rhome = R_get_HOME();
    fprintf(stderr, "R_get_HOME() gave '%s'\n", Rhome);
    if (Rhome == NULL) {
	fprintf(stderr, "To use Rlib, the variable R_HOME must be set\n");
	err = E_EXTERNAL;
	goto bailout;
    } else {
	set_path_for_Rlib(Rhome);
    }
#endif

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
	    Rinit = 1;
	} else {
	    close_plugin(Rhandle);
	    Rhandle = NULL;
	    err = Rlib_err = E_EXTERNAL;
	}
    }

 bailout:

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
		pputs(prn, line);
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

    if (t != CLOSXP && t != BUILTINSXP &&
	t != BUILTINSXP && t != SPECIALSXP) {
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

/* gretl_R_function_add... : these functions are used in geneval.c
   to convert from gretl types to R constructs for passing to R
   functions
*/

int gretl_R_function_add_scalar (double x)
{
    current_arg = R_CDR(current_arg);
    R_SETCAR(current_arg, R_ScalarReal(x));

    return 0;
}

int gretl_R_function_add_string (const char *s)
{
    current_arg = R_CDR(current_arg);
    R_SETCAR(current_arg, R_mkString(s));

    return 0;
}

int gretl_R_function_add_vector (const double *x, int t1, int t2)
{
    SEXP res = R_allocVector(REALSXP, t2 - t1 + 1);
    int i;

    if (res == NULL) {
	return E_ALLOC;
    }

    current_arg = R_CDR(current_arg);

    for (i=t1; i<=t2; i++) {
    	R_REAL(res)[i-t1] = x[i];
    }

    R_SETCAR(current_arg, res);

    return 0;
}

int gretl_R_function_add_matrix (const gretl_matrix *m)
{
    int nr = gretl_matrix_rows(m);
    int nc = gretl_matrix_cols(m);
    SEXP res;
    int i, j;

    current_arg = R_CDR(current_arg);
    res = R_allocMatrix(REALSXP, nr, nc);
    if (res == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<nr; i++) {
	for (j=0; j<nc; j++) {
	    R_REAL(res)[i + j * nr] = gretl_matrix_get(m, i, j);
    	}
    }

    R_SETCAR(current_arg, res);

    return 0;
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

static int R_type_to_gretl_type (SEXP s)
{
    if (R_isMatrix(s)) {
	return GRETL_TYPE_MATRIX;
    } else if (R_isLogical(s)) {
	return GRETL_TYPE_BOOL;
    } else if (R_isInteger(s)) {
	return GRETL_TYPE_INT;
    } else if (R_isReal(s)) {
	return GRETL_TYPE_DOUBLE;
    } else if (R_isString(s)) {
	return GRETL_TYPE_STRING;
    } else {
	return GRETL_TYPE_NONE;
    }
}

/* execute an R function and try to convert the value returned
   into a gretl type
*/

int gretl_R_function_exec (const char *name, int *rtype, void **ret)
{
    SEXP res;
    int err = 0;

    if (gretl_messages_on()) {
	R_PrintValue(current_call);
    }

    res = R_catch(current_call, VR_GlobalEnv, &err);
    if (err) {
	fprintf(stderr, "gretl_R_function_exec: R_catch failed on %s\n", name);
	return E_EXTERNAL;
    }

    *rtype = R_type_to_gretl_type(res);

#if FDEBUG
    printf("R return value: got type %d (%s)\n", *rtype,
	   gretl_type_get_name(*rtype));
    printf("Calling R_PrintValue() on @res\n");
    R_PrintValue(res);
#endif

    if (*rtype == GRETL_TYPE_MATRIX) {
	gretl_matrix *m = NULL;
	int nr = 0, nc = 0;

	if (!R_isReal(res)) {
	    gretl_errmsg_sprintf("%s: got 'matrix' result, but not of type real", name);
	    err = E_TYPES;
	} else {
	    nr = R_nrows(res);
	    nc = R_ncols(res);

	    if (nr > 0 && nc > 0) {
		m = gretl_matrix_alloc(nr, nc);
		if (m == NULL) {
		    err = E_ALLOC;
		}
	    } else if (nr == 0 && nc == 0) {
		m = gretl_null_matrix_new();
	    } else {
		gretl_errmsg_sprintf("%s: invalid matrix dimensions, %d x %d",
				     name, nr, nc);
		err = E_DATA;
	    }
	}

	if (m != NULL && nr > 0) {
	    int i, j;

	    for (i=0; i<nr; i++) {
		for (j=0; j<nc; j++) {
		    gretl_matrix_set(m, i, j, R_REAL(res)[i + j * nr]);
		}
	    }
	}
	R_unprotect(1);
	*ret = m;
    } else if (gretl_scalar_type(*rtype)) {
	double *realres = R_REAL(res);
	double *dret = *ret;

	*dret = *realres;
    	R_unprotect(1);
    } else if (*rtype == GRETL_TYPE_STRING) {
	SEXP *rsp = R_STRING_PTR(res);
	const char *s = R_STRING(*rsp);

	*ret = gretl_strdup(s);
	R_unprotect(1);
    } else {
	err = E_TYPES;
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

    /* by passing OPT_L below we indicate that we're
       using the library */
    err = write_R_source_file(buf, dset, opt | OPT_L);
    if (!err) {
	err = lib_run_Rlib_sync(opt, prn);
    }

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
	gretl_errmsg_sprintf("%s: a block is already started",
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
	gretl_errmsg_sprintf("%s: no block is in progress",
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
    int err = write_gretl_R_files(buf, dset, opt);

    if (err) {
	delete_gretl_R_files();
    } else {
	err = lib_run_R_sync(opt, prn);
    }

    return err;
}

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

int foreign_execute (const DATASET *dset,
		     gretlopt opt, PRN *prn)
{
    int i, err = 0;

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
	    err = lib_run_mpi_sync(foreign_opt, prn);
	}
	foreign_destroy();
	return err; /* handled */
    }
#endif

    if (foreign_lang == LANG_R) {
#ifdef USE_RLIB
	if (gretl_use_Rlib()) {
	    err = run_R_lib(NULL, dset, foreign_opt, prn);
	} else {
	    err = run_R_binary(NULL, dset, foreign_opt, prn);
	}
#else
	err = run_R_binary(NULL, dset, foreign_opt, prn);
#endif
	foreign_destroy();
	return err; /* handled */
    }

    if (foreign_lang == LANG_OX) {
	err = write_gretl_ox_script(NULL, foreign_opt, NULL);
	if (err) {
	    delete_ox_script();
	}
    } else if (foreign_lang == LANG_OCTAVE) {
	err = write_gretl_octave_script(NULL, foreign_opt,
					dset, NULL);
	if (err) {
	    delete_octave_script();
	}
    } else if (foreign_lang == LANG_STATA) {
	err = write_gretl_stata_script(NULL, foreign_opt,
				       dset, NULL);
	if (err) {
	    delete_stata_script();
	}
    } else if (foreign_lang == LANG_PYTHON) {
	err = write_gretl_python_script(NULL, foreign_opt,
					NULL);
	if (err) {
	    delete_python_script();
	}
    } else if (foreign_lang == LANG_JULIA) {
	err = write_gretl_julia_script(NULL, foreign_opt,
				       NULL);
	if (err) {
	    delete_julia_script();
	}
    } else {
	/* "can't happen" */
	err = E_DATA;
    }

    if (!err) {
	err = lib_run_other_sync(foreign_opt, prn);
    }

    foreign_destroy();

    return err;
}

/**
 * execute_R_buffer:
 * @buf: buffer containing commands.
 * @dset: dataset struct.
 * @opt: may include %OPT_D to send data from gretl.
 * @prn: struct for printing output.
 *
 * This is used only for MS Windows, working around
 * breakage in previously coded non-interactive calls to R
 * executable(s).
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
    if (gretl_use_Rlib()) {
	err = run_R_lib(buf, dset, opt, prn);
    } else {
	err = run_R_binary(buf, dset, opt, prn);
    }
#else
    err = run_R_binary(buf, dset, opt, prn);
#endif


    return err;
}
