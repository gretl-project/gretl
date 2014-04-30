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

#ifdef HAVE_MPI
# include "gretl_mpi.h"
#endif

#include <glib.h>

#ifdef USE_RLIB
# include <Rinternals.h> /* for SEXP and friends */
#endif

#ifdef G_OS_WIN32
# include <windows.h>
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

/* "dotdir" filenames */
static gchar *gretl_dot_dir;
static gchar *gretl_Rprofile;
static gchar *gretl_Rsrc;
static gchar *gretl_Rout;
static gchar *gretl_Rmsg;
static gchar *gretl_ox_prog;
static gchar *gretl_octave_prog;
static gchar *gretl_stata_prog;
static gchar *gretl_python_prog;
#ifdef HAVE_MPI
static gchar *gretl_mpi_prog;
#endif

static void destroy_foreign (void)
{
    strings_array_free(foreign_lines, foreign_n_lines);
    foreign_lines = NULL;
    foreign_started = 0;
    foreign_n_lines = 0;
    foreign_opt = OPT_NONE;
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

static const gchar *gretl_ox_filename (void)
{
    if (gretl_ox_prog == NULL) {
	const char *dotdir = gretl_dotdir();

	gretl_ox_prog = g_strdup_printf("%sgretltmp.ox", dotdir);
    }

    return gretl_ox_prog;
}

static const gchar *gretl_octave_filename (void)
{
    if (gretl_octave_prog == NULL) {
	const char *dotdir = gretl_dotdir();

	gretl_octave_prog = g_strdup_printf("%sgretltmp.m", dotdir);
    }

    return gretl_octave_prog;
}

static const gchar *gretl_stata_filename (void)
{
    if (gretl_stata_prog == NULL) {
	const char *dotdir = gretl_dotdir();

	gretl_stata_prog = g_strdup_printf("%sgretltmp.do", dotdir);
    }

    return gretl_stata_prog;
}

static const gchar *gretl_python_filename (void)
{
    if (gretl_python_prog == NULL) {
	const char *dotdir = gretl_dotdir();

	gretl_python_prog = g_strdup_printf("%sgretltmp.py", dotdir);
    }

    return gretl_python_prog;
}

#ifdef HAVE_MPI

static const gchar *gretl_mpi_filename (void)
{
    if (gretl_mpi_prog == NULL) {
	const char *dotdir = gretl_dotdir();

	gretl_mpi_prog = g_strdup_printf("%sgretltmp-mpi.inp", dotdir);
    }

    return gretl_mpi_prog;
}

#endif

static void make_gretl_R_names (void)
{
    static int done;

    if (!done) {
	gretl_dot_dir = g_strdup(gretl_dotdir());
#ifdef G_OS_WIN32
	slash_convert(gretl_dot_dir, FROM_BACKSLASH);
#endif
	gretl_Rprofile = g_strdup_printf("%sgretl.Rprofile", gretl_dot_dir);
	gretl_Rsrc = g_strdup_printf("%sRsrc", gretl_dot_dir);
	gretl_Rout = g_strdup_printf("%sR.out", gretl_dot_dir);
	gretl_Rmsg = g_strdup_printf("%sR.msg", gretl_dot_dir);
	done = 1;
    }
}

#ifdef G_OS_WIN32

static int lib_run_R_sync (gretlopt opt, PRN *prn)
{
    gchar *cmd;
    int err = 0;

    cmd = g_strdup_printf("\"%s\" --no-save --no-init-file --no-restore-data "
			  "--slave", gretl_rbin_path());

#if FDEBUG
    fprintf(stderr, "Running R binary, '%s'\n", gretl_rbin_path());
#endif

    err = win_run_sync(cmd, NULL);

    if (!(opt & OPT_Q)) {
	const gchar *outname;
	FILE *fp;

	outname = (err)? gretl_Rmsg : gretl_Rout;
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

#if FDEBUG
    fprintf(stderr, "win_run_sync: err = %d\n", err);
#endif

    g_free(cmd);

    return err;
}

static int lib_run_other_sync (gretlopt opt, PRN *prn)
{
    const char *path;
    const char *fname;
    gchar *cmd, *sout = NULL;
    int err;

    if (foreign_lang == LANG_OX) {
	path = gretl_oxl_path();
	fname = gretl_ox_filename();
	cmd = g_strdup_printf("\"%s\" \"%s\"", path, fname);
    } else if (foreign_lang == LANG_OCTAVE) {
	path = gretl_octave_path();
	fname = gretl_octave_filename();
	cmd = g_strdup_printf("\"%s\" --silent \"%s\"", path, fname);
    } else if (foreign_lang == LANG_STATA) {
	path = gretl_stata_path();
	cmd = g_strdup_printf("%s /q /e gretltmp.do", path);
    } else if (foreign_lang == LANG_PYTHON) {
	path = gretl_python_path();
	fname = gretl_python_filename();
	cmd = g_strdup_printf("%s \"%s\"", path, fname);
    } else {
	return 1;
    }

    err = gretl_win32_grab_output(cmd, &sout);

    if (sout != NULL && *sout != '\0') {
	pputs(prn, sout);
    }

    g_free(sout);
    g_free(cmd);

    return err;
}

#ifdef HAVE_MPI

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
	
	cmd = g_strdup_printf("%s%s%s \"%sgretlcli-mpi\"%s%s \"%s\"",
			      mpiexec, hostbit, npbit, gretl_home(), rngbit,
			      qopt, gretl_mpi_filename());

	err = gretl_win32_grab_output(cmd, &sout);
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

#endif

static char *win32_dotpath (void)
{
    gchar *dotpath = g_strdup(gretl_dotdir());

    gretl_charsub(dotpath, '\\', '/');
    return dotpath;
}

#else /* !G_OS_WIN32 */

/* print from Stata's batch logfile */

static void do_stata_printout (PRN *prn)
{
    gchar *buf = NULL;

    /* we're located in dotdir at this point */

    if (g_file_get_contents("gretltmp.log", &buf, NULL, NULL)) {
	pputs(prn, buf);
	g_free(buf);
	pputc(prn, '\n');
    }
}

static int lib_run_prog_sync (char **argv, gretlopt opt, PRN *prn)
{
    gchar *sout = NULL;
    gchar *errout = NULL;
    gint status = 0;
    GError *gerr = NULL;
    int err = 0;

    g_spawn_sync(NULL, argv, NULL, G_SPAWN_SEARCH_PATH,
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
	argv[1] = (char *) gretl_ox_filename();
	argv[2] = NULL;
    } else if (foreign_lang == LANG_OCTAVE) {
	argv[0] = (char *) gretl_octave_path();
	argv[1] = "--silent";
	argv[2] = (char *) gretl_octave_filename();
	argv[3] = NULL;
    } else if (foreign_lang == LANG_PYTHON) {
	argv[0] = (char *) gretl_python_path();
	argv[1] = (char *) gretl_python_filename();
	argv[2] = NULL;
    } else if (foreign_lang == LANG_STATA) {
	argv[0] = (char *) gretl_stata_path();
	argv[1] = "-q";
	argv[2] = "-b";
	argv[3] = "do";
	argv[4] = "gretltmp.do";
	argv[5] = NULL;
	/* otherwise there's no way to control the location
	   of the stata output (gretltmp.log)
	*/
	gretl_chdir(gretl_dotdir());
    }

    err = lib_run_prog_sync(argv, opt, prn);

    return err;
}

#ifdef HAVE_MPI

#define MPI_DEBUG 0

#if MPI_DEBUG
static void print_mpi_command (char **argv)
{
    int i;

    fputs("gretl/MPI argv array:\n ", stderr);
    for (i=0; argv[i] != NULL; i++) {
	fprintf(stderr, "%s ", argv[i]);
    }
    fputc('\n', stderr);
}
#endif

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
	ret = g_strdup_printf("%s/bin/gretlcli-mpi", tmp);
    } else {
	ret = g_strdup("gretlcli-mpi");
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
	char *argv[10];
	int i = 0;

	if (!(opt & OPT_L) && hostfile != NULL && *hostfile != '\0') {
	    hostsopt = (mpi_variant == MPI_MPICH)? "-machinefile" :
		"--hostfile";
	} else if (*npnum == '\0') {
	    /* no hosts file: supply a default np value */
	    sprintf(npnum, "%d", gretl_n_processors());
	}

	argv[i++] = (char *) mpiexec;
	if (hostsopt != NULL) {
	    argv[i++] = (char *) hostsopt;
	    argv[i++] = (char *) hostfile;
	}
	if (*npnum != '\0') {
	    argv[i++] = "-np";
	    argv[i++] = npnum;
	}
	argv[i++] = mpiprog;
	if (opt & OPT_S) {
	    argv[i++] = "--single-rng";
	}
	if (opt & OPT_Q) {
	    argv[i++] = "--quiet";
	}
	argv[i++] = (char *) gretl_mpi_filename();
	argv[i] = NULL;

#if MPI_DEBUG
	print_mpi_command(argv);
#endif
	err = lib_run_prog_sync(argv, opt, prn);
	g_free(mpiprog);
    }

    return err;
}

#endif /* HAVE_MPI */

#endif /* switch on MS Windows or not */

static int write_ox_io_file (void)
{
    static int written;

    if (!written) {
	const char *dotdir = gretl_dotdir();
	gchar *fname;
	FILE *fp;

	fname = g_strdup_printf("%sgretl_io.ox", dotdir);
	fp = gretl_fopen(fname, "w");
	g_free(fname);

	if (fp == NULL) {
	    return E_FOPEN;
	} else {
#ifdef G_OS_WIN32
            gchar *dotcpy = win32_dotpath();

	    fputs("gretl_dotdir ()\n{\n", fp);
	    fprintf(fp, "  return \"%s\";\n", dotcpy);
	    g_free(dotcpy);
#else
	    fputs("gretl_dotdir ()\n{\n", fp);
	    fprintf(fp, "  return \"%s\";\n", dotdir);
#endif	
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

static int write_octave_io_file (void)
{
    static int written;

    if (!written) {
	const char *dotdir = gretl_dotdir();
	gchar *fname;
	FILE *fp;

	fname = g_strdup_printf("%sgretl_io.m", dotdir);
	fp = gretl_fopen(fname, "w");
	g_free(fname);

	if (fp == NULL) {
	    return E_FOPEN;
	} else {
#ifdef G_OS_WIN32
            gchar *dotcpy = win32_dotpath();

	    fputs("function dotdir = gretl_dotdir()\n", fp);
	    fprintf(fp, "  dotdir = \"%s\";\n", dotcpy);
	    g_free(dotcpy);
#else
	    fputs("function dotdir = gretl_dotdir()\n", fp);
	    fprintf(fp, "  dotdir = \"%s\";\n", dotdir);
#endif	
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
	    written = 1;
	}
    }

    return 0;
}

static int write_python_io_file (void)
{
    static int written;

    if (!written) {
	const char *dotdir = gretl_dotdir();
	gchar *fname;
	FILE *fp;

	fname = g_strdup_printf("%sgretl_io.py", dotdir);
	fp = gretl_fopen(fname, "w");
	g_free(fname);

	if (fp == NULL) {
	    return E_FOPEN;
	} else {
#ifdef G_OS_WIN32
            gchar *dotcpy = win32_dotpath();

	    fprintf(fp, "gretl_dotdir = \"%s\"\n\n", dotcpy);
	    g_free(dotcpy);
#else
	    fprintf(fp, "gretl_dotdir = \"%s\"\n\n", dotdir);
#endif	

	    fputs("def gretl_export(X, fname, autodot=1):\n", fp);
	    fputs("  from numpy import asmatrix\n", fp);
	    fputs("  M = asmatrix(X)\n", fp);
	    fputs("  r, c = M.shape\n", fp);
	    fputs("  if autodot:\n", fp);
            fputs("    fname = gretl_dotdir + fname\n", fp);
	    fputs("  f = open(fname, 'w')\n", fp);
	    fputs("  f.write(repr(r) + '\\t' + repr(c) + '\\n')\n", fp);
	    fputs("  for i in range(0, r):\n", fp);
	    fputs("    for j in range(0, c):\n", fp);
	    fputs("      f.write('%.18e ' % M[i,j])\n", fp);
            fputs("    f.write('\\n')\n", fp);
	    fputs("  f.close()\n\n", fp);

	    fputs("def gretl_loadmat(fname, autodot=1):\n", fp);
	    fputs("  from numpy import loadtxt\n", fp);
	    fputs("  if autodot:\n", fp);
	    fputs("    fname = gretl_dotdir + fname\n", fp);
	    fputs("  M = loadtxt(fname, skiprows=1)\n", fp);
	    fputs("  return M\n\n", fp);

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
	const char *dotdir = gretl_dotdir();
	gchar *fname;
	FILE *fp;

	fname = g_strdup_printf("%sgretl_export.ado", dotdir);
	fp = gretl_fopen(fname, "w");
	g_free(fname);

	if (fp == NULL) {
	    return E_FOPEN;
	} else {
	    fputs("program define gretl_export\n", fp);
	    /* not sure about req'd version, but see mat2txt.ado */
	    fputs("version 8.2\n", fp);
	    fputs("local matrix `1'\n", fp);
	    fputs("local fname `2'\n", fp);
	    fputs("tempname myfile\n", fp);
	    fputs("file open `myfile' using \"`fname'\", "
		  "write text replace\n", fp);
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

static void add_gretl_include (int lang, FILE *fp)
{
    if (lang == LANG_PYTHON) {
	fputs("from gretl_io import gretl_dotdir, gretl_loadmat, "
	      "gretl_export\n", fp);
	return;
    }

#ifdef G_OS_WIN32
    gchar *dotcpy;

    if (lang == LANG_OX) {
	dotcpy = win32_dotpath();
	if (strchr(dotcpy, ' ')) {
	    fprintf(fp, "#include \"%sgretl_io.ox\"\n", dotcpy);
	} else {
	    fprintf(fp, "#include <%sgretl_io.ox>\n", dotcpy);
	}
	g_free(dotcpy);
    } else if (lang == LANG_OCTAVE) {
	dotcpy = win32_dotpath();
	fprintf(fp, "source(\"%sgretl_io.m\")\n", dotcpy);
	g_free(dotcpy);
    } else if (lang == LANG_STATA) {
	fprintf(fp, "quietly adopath + \"%s\"\n", gretl_dotdir());
    }
#else
    const char *dotdir = gretl_dotdir();

    if (lang == LANG_OX) {
	if (strchr(dotdir, ' ')) {
	    fprintf(fp, "#include \"%sgretl_io.ox\"\n", dotdir);
	} else {
	    fprintf(fp, "#include <%sgretl_io.ox>\n", dotdir);
	}
    } else if (lang == LANG_OCTAVE) {
	fprintf(fp, "source(\"%sgretl_io.m\")\n", dotdir);
    } else if (lang == LANG_STATA) {
	fprintf(fp, "quietly adopath + \"%s\"\n", dotdir);
    }
#endif
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
		add_gretl_include(LANG_OX, fp);
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
		add_gretl_include(LANG_OX, fp);
	    }
	}
    }

    bufgets_finalize(buf);
}

/**
 * write_gretl_ox_file:
 * @buf: text buffer containing Ox code.
 * @opt: should contain %OPT_G for use from GUI.
 * @pfname: location to receive name of file written, or %NULL.
 *
 * Writes the content of @buf into a file in the gretl user's
 * "dotdir".
 *
 * Returns: 0 on success, non-zero on error.
 */

int write_gretl_ox_file (const char *buf, gretlopt opt, const char **pfname)
{
    const gchar *fname = gretl_ox_filename();
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
 * write_gretl_python_file:
 * @buf: text buffer containing Python code.
 * @opt: should contain %OPT_G for use from GUI.
 * @pfname: location to receive name of file written, or %NULL.
 *
 * Writes the content of @buf into a file in the gretl user's
 * "dotdir".
 *
 * Returns: 0 on success, non-zero on error.
 */

int write_gretl_python_file (const char *buf, gretlopt opt, const char **pfname)
{
    const gchar *fname = gretl_python_filename();
    FILE *fp = gretl_fopen(fname, "w");

    write_python_io_file();

    if (fp == NULL) {
	return E_FOPEN;
    } else {
	add_gretl_include(LANG_PYTHON, fp);
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

#ifdef HAVE_MPI

static int mpi_send_funcs_setup (FILE *fp)
{
    const char *dotdir = gretl_dotdir();
    gchar *fname;
    int err;

    fname = g_strdup_printf("%smpi-funcs-tmp.xml", dotdir);
    err = write_session_functions_file(fname);

    if (!err) {
	fprintf(fp, "include \"%s\"\n", fname);
    }

    g_free(fname);

    return err;
}

static int write_gretl_mpi_file (gretlopt opt)
{
    const gchar *fname = gretl_mpi_filename();
    FILE *fp = gretl_fopen(fname, "w");
    int err = 0;

    if (fp == NULL) {
	return E_FOPEN;
    }

    if (opt & OPT_F) {
	/* honor the --send-functions option */
	if (n_user_functions() > 0) {
	    err = mpi_send_funcs_setup(fp);
	}
    }

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
    char save_na[8];
    gchar *sdata;
    int err;

    *save_na = '\0';
    strncat(save_na, get_csv_na_write_string(), 7);
    set_csv_na_write_string(".");
    sdata = g_strdup_printf("%sstata.csv", gretl_dotdir());
    err = write_data(sdata, NULL, dset, OPT_C, NULL);
    set_csv_na_write_string(save_na);
 
    if (err) {
	gretl_errmsg_sprintf("write_data_for_stata: failed with err = %d\n", err);
    } else {
	fputs("* load data from gretl\n", fp);
	fputs("insheet using \"stata.csv\"\n", fp);
    }

    g_free(sdata);

    return err;
}

int write_gretl_stata_file (const char *buf, gretlopt opt, 
			    const DATASET *dset,
			    const char **pfname)
{
    const gchar *fname = gretl_stata_filename();
    FILE *fp = gretl_fopen(fname, "w");

    write_stata_io_file();

    if (fp == NULL) {
	return E_FOPEN;
    } else {
	int err;

	/* source the I-O functions */
	add_gretl_include(LANG_STATA, fp);
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
    gchar *mdata;
    int err;

    mdata = g_strdup_printf("%smdata.tmp", gretl_dotdir());
    err = write_data(mdata, NULL, dset, OPT_M, NULL);
 
    if (err) {
	gretl_errmsg_sprintf("write_data_for_octave: failed with err = %d\n", err);
    } else {
	fputs("% load data from gretl\n", fp);
	fprintf(fp, "load '%s'\n", mdata);
    }

    g_free(mdata);

    return err;
}

int write_gretl_octave_file (const char *buf, gretlopt opt, 
			     const DATASET *dset,
			     const char **pfname)
{
    const gchar *fname = gretl_octave_filename();
    FILE *fp = gretl_fopen(fname, "w");

    write_octave_io_file();

    if (fp == NULL) {
	return E_FOPEN;
    } else {
	int err;

	/* source the I-O functions */
	add_gretl_include(LANG_OCTAVE, fp);
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

/* write out current dataset in R format, and, if this succeeds,
   write appropriate R commands to @fp to source the data
*/

static int write_data_for_R (const DATASET *dset,
			     gretlopt opt,
			     FILE *fp)
{
    int ts = dataset_is_time_series(dset);
    gchar *Rdata;
    int err;

    Rdata = g_strdup_printf("%sRdata.tmp", gretl_dot_dir);

    err = write_data(Rdata, NULL, dset, OPT_R, NULL);
    if (err) {
	gretl_errmsg_sprintf("write_data_for_R: failed with err = %d\n", err);
	g_free(Rdata);
	return err;
    }

    fputs("# load data from gretl\n", fp);

    if (ts) {
	char *p, datestr[OBSLEN];
	int subper = 1;
	    
	ntodate(datestr, dset->t1, dset);
	p = strchr(datestr, ':');
	if (p != NULL) {
	    subper = atoi(p + 1);
	}

	fprintf(fp, "gretldata <- read.table(\"%s\", header=TRUE)\n", Rdata);
	fprintf(fp, "gretldata <- ts(gretldata, start=c(%d, %d), frequency = %d)\n", 
		atoi(datestr), subper, dset->pd);
    } else {	
	fprintf(fp, "gretldata <- read.table(\"%s\", header=TRUE)\n", Rdata);
	fputs("attach(gretldata)\n", fp);
    }

    g_free(Rdata);

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
    fprintf(fp, "gretl.dotdir <- \"%s\"\n", gretl_dot_dir);

    fputs("gretl.export <- function(x) {\n", fp);
    fprintf(fp, "  prefix <- \"%s\"\n", gretl_dot_dir);
    fputs("  sx <- as.character(substitute(x))\n", fp);
    fputs("  if (is.ts(x)) {\n", fp);
    fputs("    fname <- paste(prefix, sx, \".csv\", sep=\"\")\n", fp);
    fputs("    dfx <- data.frame(x)\n", fp);
    fputs("    if (ncol(dfx) == 1) {\n", fp);
    fputs("      colnames(dfx) <- sx;\n", fp);
    fputs("    }\n", fp);
    fputs("    write.csv(dfx, file=fname, row.names=F)\n", fp);
    fputs("  } else if (is.data.frame(x)) {\n", fp);
    fputs("    fname <- paste(prefix, sx, \".csv\", sep=\"\")\n", fp);
    fputs("    write.csv(x, file=fname, row.names=F)\n", fp);
    fputs("  } else if (is.matrix(x)) {\n", fp);
    fputs("    fname <- paste(prefix, sx, \".mat\", sep=\"\")\n", fp);
    fputs("    write(dim(x), fname)\n", fp);
    fputs("    write(t(x), file=fname, ncolumns=ncol(x), append=TRUE)\n", fp);
    fputs("  }\n", fp);
    fputs("  gretlmsg <- paste(\"wrote\", fname, \"\\n\")\n", fp);
    fputs("  cat(gretlmsg)\n", fp);
    fputs("}\n", fp);

    fputs("gretl.loadmat <- function(mname) {\n", fp);
    fprintf(fp, "  prefix <- \"%s\"\n", gretl_dot_dir);
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
    int err;

#if FDEBUG
    printf("writing R profile: starting\n");
#endif

    err = gretl_setenv("R_PROFILE", gretl_Rprofile);
    if (err) {
	return err;
    }     

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
    printf("writing R profile: returning %d\n", err);
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
    printf("write R source file: starting\n");
#endif

    if (fp == NULL) {
	err = E_FOPEN;
    } else {
	int sunk = 0;

#ifdef G_OS_WIN32
	if (!(opt & (OPT_I | OPT_L))) {
	    /* non-interactive, but not using Rlib */
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
		printf("Rlib: writing 'startup' material\n");
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
#if 0 /* ifdef G_OS_WIN32 */
	    /* not working, 2012-12-27 */
	    maybe_print_R_path_addition(fp);
#endif
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
	    /* pass on the script supplied in the 'buf' argument */
	    fputs("# load script from gretl\n", fp);
	    fputs(buf, fp);
	} else if (!(opt & OPT_G)) {
	    /* non-GUI */
	    put_foreign_lines(fp);
	}

	if (sunk) {
	    fputs("sink()\n", fp);
	}

#ifdef G_OS_WIN32
	if (!(opt & (OPT_I | OPT_L))) {
	    /* Rterm on Windows won't exit without this? */
	    fputs("q()\n", fp);
	}
#endif

	fclose(fp);
    }

#if FDEBUG
    printf("write R source file: returning %d\n", err);
#endif

    return err;
}

/* Write files to be read by R: profile to be read on startup and
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
    printf("write_gretl_R_files: starting\n");
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
    printf("write_gretl_R_files: returning %d\n", err);
#endif

    return err;
}

void delete_gretl_R_files (void)
{
#if FDEBUG
    printf("deleting gretl R files...\n");
#endif

    if (gretl_Rprofile != NULL) {
	gretl_remove(gretl_Rprofile);
    }
    if (gretl_Rsrc != NULL) {
	gretl_remove(gretl_Rsrc);
    }
}

static void delete_gretl_ox_file (void)
{
    if (gretl_ox_prog != NULL) {
	gretl_remove(gretl_ox_prog);
    }
}

static void delete_gretl_octave_file (void)
{
    if (gretl_octave_prog != NULL) {
	gretl_remove(gretl_octave_prog);
    }
}

static void delete_gretl_stata_file (void)
{
    if (gretl_stata_prog != NULL) {
	gretl_remove(gretl_stata_prog);
    }
}

static void delete_gretl_python_file (void)
{
    if (gretl_python_prog != NULL) {
	gretl_remove(gretl_python_prog);
    }
}

#ifdef HAVE_MPI

static void delete_gretl_mpi_file (void)
{
    if (gretl_mpi_prog != NULL) {
	gretl_remove(gretl_mpi_prog);
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
    printf("Loading libR symbols from '%s'\n", libpath);
#endif

    Rhandle = gretl_dlopen(libpath, 1);
    if (Rhandle == NULL) {
	err = E_EXTERNAL;
	goto bailout;
    } 

    R_CDR           = dlget(Rhandle, "CDR", &err);
    R_REAL          = dlget(Rhandle, "REAL", &err);
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

#if FDEBUG
    printf("load_R_symbols: returning %d\n", err);
#endif

    return err;
}

void gretl_R_cleanup (void)
{
#if FDEBUG
    printf("gretl_R_cleanup: Rinit = %d\n", Rinit);
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

static void set_path_for_Rlib (const char *Rhome)
{
#ifdef _WIN64
    const char *arch = "x64";
#else
    const char *arch = "i386";
#endif
    char *path = getenv("PATH");
    gchar *Rpath;

    Rpath = g_strdup_printf("%s\\bin\\%s", Rhome, arch);
    fprintf(stderr, "Rpath = '%s'\n", Rpath);

    if (path != NULL && strstr(path, Rpath) != NULL) {
	; /* nothing to be done */
    } else {
	g_free(Rpath);
	Rpath = g_strdup_printf("%s;%s\\bin\\%s", path, Rhome, arch);
	gretl_setenv("PATH", Rpath);
	g_free(Rpath);
	Rpath = NULL;
    }

    if (Rpath != NULL) {
	g_free(Rpath);
    }
}

#endif /* WIN32 */

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
    printf("gretl_Rlib_init: starting\n");
#endif

#ifndef WIN32
    Rhome = getenv("R_HOME");
    if (Rhome == NULL) {
	fprintf(stderr, "To use Rlib, the variable R_HOME must be set\n");
	err = E_EXTERNAL;
	goto bailout;
    }
#endif

    err = load_R_symbols();
    if (err) {
	fprintf(stderr, "gretl_Rlib_init: failed to load R functions\n");
	goto bailout;
    }

#ifdef WIN32
    Rhome = R_get_HOME();
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
    printf("lib_run_Rlib_sync: starting\n");
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
    printf("lib_run_Rlib_sync: returning %d\n", err);
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
    printf("gretl_use_Rlib: starting\n");
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
    printf("gretl_use_Rlib: using %s\n", (ret)? "library" : "executable");
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

/* gretl_R_function_add... : these functions are used to convert from
   gretl types to R constructs for passing to R functions
*/

int gretl_R_function_add_scalar (double x) 
{
    current_arg = R_CDR(current_arg);
    R_SETCAR(current_arg, R_ScalarReal(x));

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
	return E_EXTERNAL;
    }

    *rtype = R_type_to_gretl_type(res);

    if (*rtype == GRETL_TYPE_MATRIX) {
	gretl_matrix *pm;
	int nc = R_ncols(res);
	int nr = R_nrows(res);
	int i, j;

	pm = gretl_matrix_alloc(nr, nc);
	if (pm == NULL) {
	    return E_ALLOC;
	}

	for (i=0; i<nr; i++) {
	    for (j=0; j<nc; j++) {
		gretl_matrix_set(pm, i, j, R_REAL(res)[i + j * nr]);
	    }
	}
	R_unprotect(1);
	*ret = pm;
    } else if (gretl_scalar_type(*rtype)) {
	double *realres = R_REAL(res);
	double *dret = *ret;

	*dret = *realres;

    	R_unprotect(1);
    } else {
	err = E_EXTERNAL;
    }

    return err;
}

static int run_R_lib (const DATASET *dset, 
		      gretlopt opt, PRN *prn)
{
    int err;

    /* we don't want gretl.Rprofile in the way */
    gretl_remove(gretl_Rprofile);

    /* by passing OPT_L below we indicate that we're
       using the library */
    err = write_R_source_file(NULL, dset, opt | OPT_L);
    if (!err) {
	err = lib_run_Rlib_sync(opt, prn);
    }

    return err;
}

#endif /* USE_RLIB */

static int foreign_block_init (const char *line, gretlopt opt, PRN *prn)
{
    int err = 0;

    foreign_opt = OPT_NONE;

    if (!strncmp(line, "foreign ", 8)) {
	char lang[16];

	line += 8;
	line += strspn(line, " ");
	if (!strncmp(line, "language", 8)) {
	    line += 8;
	    line += strspn(line, " =");
	    if (sscanf(line, "%15s", lang) == 1) {
		err = set_foreign_lang(lang, prn);
	    } else {
		err = E_PARSE;
	    }
	} else {
	    err = E_PARSE;
	}
    } else if (!strncmp(line, "mpi", 3)) {
	err = set_foreign_lang("mpi", prn);
    } else {
	/* we'll default to R for now */
	foreign_lang = LANG_R;
    }

    if (!err) {
	foreign_opt = opt;
    }

    return err;
}

/**
 * foreign_append_line:
 * @line: command line.
 * @opt: may include %OPT_V for verbose operation
 * @prn: struct for printing output.
 *
 * Appends @line to an internally stored block of "foreign"
 * commands, or starts a new block if no such block is
 * currently defined.
 * 
 * Returns: 0 on success, non-zero on error.
 */

int foreign_append_line (const char *line, gretlopt opt, PRN *prn)
{
    int err = 0;

    if (string_is_blank(line)) {
	return 0;
    }

#if 0
    fprintf(stderr, "foreign_append_line: '%s'\n", line);
#endif

    if (!foreign_started) {
	/* starting from scratch */
	err = foreign_block_init(line, opt, prn);
	if (!err) {
	    foreign_started = 1;
	}
    } else {
	/* appending */
	err = strings_array_add(&foreign_lines, &foreign_n_lines, line);
	if (err) {
	    destroy_foreign();
	}
    }	

    return err;
}

static int run_R_binary (const DATASET *dset, 
			 gretlopt opt, PRN *prn)
{
    int err;

    /* write both profile and Rsrc files */

    err = write_gretl_R_files(NULL, dset, opt);
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
	err = write_gretl_mpi_file(foreign_opt);
	if (err) {
	    delete_gretl_mpi_file();
	} else {
	    err = lib_run_mpi_sync(foreign_opt, prn);
	}
	destroy_foreign();
	return err;
    }
#endif

    if (foreign_lang == LANG_R) {
#ifdef USE_RLIB
	if (gretl_use_Rlib()) {
	    err = run_R_lib(dset, foreign_opt, prn);
	} else {
	    err = run_R_binary(dset, foreign_opt, prn);
	}
#else
	err = run_R_binary(dset, foreign_opt, prn);
#endif
    } else if (foreign_lang == LANG_OX) {
	err = write_gretl_ox_file(NULL, foreign_opt, NULL);
	if (err) {
	    delete_gretl_ox_file();
	} else {
	    err = lib_run_other_sync(foreign_opt, prn);
	}
    } else if (foreign_lang == LANG_OCTAVE) {
	err = write_gretl_octave_file(NULL, foreign_opt, 
				      dset, NULL);
	if (err) {
	    delete_gretl_octave_file();
	} else {
	    err = lib_run_other_sync(foreign_opt, prn);
	}
    } else if (foreign_lang == LANG_STATA) {
	err = write_gretl_stata_file(NULL, foreign_opt, 
				     dset, NULL);
	if (err) {
	    delete_gretl_stata_file();
	} else {
	    err = lib_run_other_sync(foreign_opt, prn);
	}
    } else if (foreign_lang == LANG_PYTHON) {
	err = write_gretl_python_file(NULL, foreign_opt, 
				      NULL);
	if (err) {
	    delete_gretl_python_file();
	} else {
	    err = lib_run_other_sync(foreign_opt, prn);
	}	
    } else {
	/* "can't happen" */
	err = E_DATA;
    }
    
    destroy_foreign();

    return err;
}
