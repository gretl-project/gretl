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
#include <glib.h>

#ifdef G_OS_WIN32
# include <windows.h>
#else
# include <signal.h>
#endif

static char **foreign_lines;
static int foreign_n_lines; 
static int foreign_lang;
static gretlopt foreign_opt;

/* foreign_opt may include OPT_D to send data, OPT_Q to operate
   quietly (don't display output from foreign program)
*/

enum {
    LANG_R = 1,
    LANG_OX
};

static void destroy_foreign (void)
{
    free_strings_array(foreign_lines, foreign_n_lines);
    foreign_lines = NULL;
    foreign_n_lines = 0;
    foreign_opt = OPT_NONE;
}

static int set_foreign_lang (const char *lang, PRN *prn)
{
    int err = 0;

    if (!strcmp(lang, "R")) {
	foreign_lang = LANG_R;
    } else if (!strcmp(lang, "ox")) {
	foreign_lang = LANG_OX;
    } else {
	pprintf(prn, "%s: unknown language\n");
	err = E_DATA;
    }

    return err;
}

#ifdef G_OS_WIN32

static int lib_run_R_sync (PRN *prn)
{
    gchar *Rterm, *cmd;
    int err = 0;

    Rterm = R_path_from_registry();
    if (Rterm == NULL) {
	return E_EXTERNAL;
    }

    cmd = g_strdup_printf("\"%s\" --no-save --no-init-file --no-restore-data "
			  "--slave", Rterm);

    err = winfork(cmd, NULL, SW_SHOWMINIMIZED, CREATE_NEW_CONSOLE);

    if (!err) {
	gchar *Rout = g_strdup_printf("%sR.out", gretl_dot_dir());
	FILE *fp;

	fp = gretl_fopen(Rout, "r");
	if (fp == NULL) {
	    err = E_FOPEN;
	} else {
	    char line[1024];

	    while (fgets(line, sizeof line, fp)) {
		pputs(prn, line);
	    }
	    fclose(fp);
	}
	g_free(Rout);
    }

    g_free(cmd);
    g_free(Rterm);

    return err;
}

#else

static int lib_run_R_sync (gretlopt opt, PRN *prn)
{
    gchar *argv[6];
    gchar *out = NULL;
    gchar *errout = NULL;
    gint status = 0;
    GError *gerr = NULL;
    int err = 0;

    argv[0] = "R";
    argv[1] = "--no-save";
    argv[2] = "--no-init-file";
    argv[3] = "--no-restore-data";
    argv[4] = "--slave";
    argv[5] = NULL;

    signal(SIGCHLD, SIG_DFL);

    g_spawn_sync(NULL, argv, NULL, G_SPAWN_SEARCH_PATH,
		 NULL, NULL, &out, &errout,
		 &status, &gerr);

    if (gerr != NULL) {
	pprintf(prn, "%s\n", gerr->message); 
	g_error_free(gerr);
	err = 1;
    } else if (status != 0) {
	if (errout != NULL) {
	    if (*errout == '\0') {
		pprintf(prn, "R exited with status %d", status);
	    } else {
		pprintf(prn, "%s\n", errout); 
	    } 
	}
	err = 1;
    } else if (out != NULL) {
	if (!(opt & OPT_Q)) {
	    /* with OPT_Q, don't print non-error output */
	    pprintf(prn, "%s\n", out);
	}
    } else {
	pprintf(prn, "%s\n", "Got no output");
	err = 1;
    }

    g_free(out);
    g_free(errout);

    return err;
}

#endif

/* write out current dataset in R format, and, if this succeeds,
   write appropriate R commands to @fp to source the data
*/

static int write_data_for_R (const double **Z, 
			     const DATAINFO *pdinfo,
			     const gchar *dotdir,
			     FILE *fp)
{
    gchar *Rline = NULL, *Rdata = NULL;
    int err;

    Rdata = g_strdup_printf("%sRdata.tmp", dotdir);
    Rline = g_strdup_printf("store \"%s\" -r", Rdata);
    g_free(Rline);

    err = write_data(Rdata, NULL, Z, pdinfo, OPT_R, NULL);
    if (err) {
	g_free(Rdata);
	return err;
    }

    fputs("# load data from gretl\n", fp);
    fprintf(fp, "gretldata <- read.table(\"%s\", header=TRUE)\n", Rdata);
    g_free(Rdata);

    if (dataset_is_time_series(pdinfo)) {
	char *p, datestr[OBSLEN];
	int subper = 1;
	    
	ntodate_full(datestr, pdinfo->t1, pdinfo);
	p = strchr(datestr, ':');
	if (p != NULL) {
	    subper = atoi(p + 1);
	}
	    
	fprintf(fp, "gretldata <- ts(gretldata, start=c(%d, %d), frequency = %d)\n", 
		atoi(datestr), subper, pdinfo->pd);
    } else {
	fputs("attach(gretldata)\n", fp);
    }

    return err;
}

/* define an R function for passing data or matrices back to gretl */

static void write_R_export_func (const gchar *dotdir, FILE *fp) 
{
    fprintf(fp, "gretl.dotdir <- \"%s\"\n", dotdir);
    fputs("gretl.export <- function(x) {\n", fp);
    fprintf(fp, "  prefix <- \"%s\"\n", dotdir);
    fputs("  sx <- as.character(substitute(x))\n", fp);
    fputs("  if (is.ts(x)) {\n", fp);
    fputs("    fname <- paste(prefix, sx, \".csv\", sep=\"\")\n", fp);
    fputs("    dfx <- data.frame(x)\n", fp);
    fputs("    if (dim(dfx)[2] == 1) {\n", fp);
    fputs("      dfx <- as.data.frame(x);\n", fp);
    fputs("      colnames(dfx) <- sx;\n", fp);
    fputs("    }\n", fp);
    fputs("    write.csv(dfx, file=fname, row.names=F)\n", fp);
    fputs("  } else if (is.data.frame(x)) {\n", fp);
    fputs("    fname <- paste(prefix, sx, \".csv\", sep=\"\")\n", fp);
    fputs("    write.csv(x, file=fname, row.names=F)\n", fp);
    fputs("  } else if (is.matrix(x)) {\n", fp);
    fputs("    fname <- paste(prefix, sx, \".mat\", sep=\"\")\n", fp);
    fputs("    write(dim(x), fname)\n", fp);
    fputs("    write(t(x), file=fname, ncolumns=dim(x)[2], append=TRUE)\n", fp);
    fputs("  }\n", fp);
    fputs("}\n", fp);
}

/* set up for a temporary R profile, so R knows what to do */

static int write_gretl_R_profile (const char *Rprofile, const char *Rsrc,
				  const gchar *dotdir)
{
    FILE *fp;
    int err;

#ifdef G_OS_WIN32
    err = !SetEnvironmentVariable("R_PROFILE", Rprofile);
#else
    err = setenv("R_PROFILE", Rprofile, 1);
#endif
    if (err) {
	return err;
    }     

    fp = gretl_fopen(Rprofile, "w");
    if (fp == NULL) {
	return E_FOPEN;
    }

    /* profile preamble */
    fputs("vnum <- as.double(R.version$major) + (as.double(R.version$minor) / 10.0)\n", 
	  fp);
    fputs("if (vnum > 2.41) library(utils)\n", fp);
    fputs("library(stats)\n", fp);
    fputs("if (vnum <= 1.89) library(ts)\n", fp);

#if 1
    write_R_export_func(dotdir, fp);
#endif

    /* source the commands and/or data from gretl */
    fprintf(fp, "source(\"%s\", echo=TRUE)\n", Rsrc);

    fclose(fp);

    return 0;
}

static int write_R_source_file (const char *Rsrc, const char *buf,
				const double **Z, const DATAINFO *pdinfo,
				const gchar *dotdir, gretlopt opt)
{
    FILE *fp = gretl_fopen(Rsrc, "w");
    int i, err = 0;

    if (fp == NULL) {
	err = E_FOPEN;
    } else {
#ifdef G_OS_WIN32
	if (!(opt & OPT_I)) {
	    /* non-interactive */
	    gchar *Rout = g_strdup_printf("%sR.out", dotdir);

	    fprintf(fp, "sink(\"%s\")\n", Rout);
	    g_free(Rout);
	}
#endif

	if (opt & OPT_D) {
	    /* send data */
	    err = write_data_for_R(Z, pdinfo, dotdir, fp);
	}

	if (buf != NULL) {
	    /* pass on the script supplied in the 'buf' argument */
	    fputs("# load script from gretl\n", fp);
	    fputs(buf, fp);
	} else if (!(opt & OPT_G)) {
	    /* non-GUI */
	    for (i=0; i<foreign_n_lines; i++) { 
		fprintf(fp, "%s\n", foreign_lines[i]);
	    }
	}

#ifdef G_OS_WIN32
	if (!(opt & OPT_I)) {
	    /* Rterm on Windows won't exit without this? */
	    fputs("q()\n", fp);
	}
#endif

	fclose(fp);
    }

    return err;
}

int write_gretl_R_files (const char *buf,
			 const double **Z, const DATAINFO *pdinfo,
			 gretlopt opt)
{ 
    const char *dotdir = gretl_dot_dir();
    gchar *Rprofile, *Rsrc;
    gchar *dotcpy = g_strdup(dotdir);
    int err = 0;

#ifdef G_OS_WIN32
    slash_convert(dotcpy, FROM_BACKSLASH);
#endif

    /* organize filenames */
    Rprofile = g_strdup_printf("%sgretl.Rprofile", dotcpy);
    Rsrc = g_strdup_printf("%sRsrc", dotcpy);
    
    /* first write a temporary R profile, so R knows what to do */
    err = write_gretl_R_profile(Rprofile, Rsrc, dotcpy);
    if (err) {
	return err;
    }

    /* then write commands and/or data to file, to be sourced
       by the R profile */
    err = write_R_source_file(Rsrc, buf, Z, pdinfo, dotcpy, opt);

    g_free(dotcpy);
    g_free(Rprofile);
    g_free(Rsrc);

    return err;
}

void delete_gretl_R_files (void)
{
    const char *dotdir = gretl_dot_dir();
    gchar *Rprofile, *Rsrc;
    
    Rprofile = g_strdup_printf("%sgretl.Rprofile", dotdir);
    Rsrc = g_strdup_printf("%sRsrc", dotdir);

    remove(Rprofile);
    remove(Rsrc);

    g_free(Rprofile);
    g_free(Rsrc);
}

int foreign_append_line (const char *line, gretlopt opt, PRN *prn)
{
    int err = 0;

    if (string_is_blank(line)) {
	return 0;
    }

    if (foreign_n_lines == 0 && !strncmp(line, "foreign ", 8)) {
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
	if (!err) {
	    foreign_opt |= opt;
	}
    } else {
	err = strings_array_add(&foreign_lines, &foreign_n_lines, line);
	if (err) {
	    destroy_foreign();
	}
    }	

    return err;
}

/* respond to "end foreign" [opts] */

int foreign_execute (const double **Z, const DATAINFO *pdinfo, 
		     gretlopt opt, PRN *prn)
{
    int err = 0;

#if FDEBUG
    int i;

    for (i=0; i<foreign_n_lines; i++) {
	pprintf(prn, "> %s\n", foreign_lines[i]);
    }
#endif

    foreign_opt |= opt;

    if (foreign_lang == LANG_R) {
	err = write_gretl_R_files(NULL, Z, pdinfo, foreign_opt);
	if (err) {
	    delete_gretl_R_files();
	} else {
	    lib_run_R_sync(foreign_opt, prn);
	}
    } else if (foreign_lang == LANG_OX) {
	pputs(prn, "language=ox: not supported yet\n");
	err = E_NOTIMP;
    }
    
    destroy_foreign();

    return err;
}
