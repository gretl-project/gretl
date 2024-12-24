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
#include "gretl_help.h"

/* for PDF opening */
#if defined(G_OS_WIN32)
# include <windows.h>
#elif defined(__APPLE__)
# if defined(USE_CARBON)
#  include <Carbon/Carbon.h>
# else
#  include <CoreFoundation/CoreFoundation.h>
#  include <CoreServices/CoreServices.h>
# endif
#endif

static int maybe_need_recode (void)
{
    const gchar *cset = NULL;
    int utf = g_get_charset(&cset);

    return !utf;
}

static int recode_print_line (const char *s, PRN *prn)
{
    gchar *trs;
    gsize bytes;
    GError *err = NULL;

    trs = g_locale_from_utf8(s, -1, NULL, &bytes, &err);

    if (err != NULL) {
	pprintf(prn, "%s\n", err->message);
	g_error_free(err);
    } else {
	pputs(prn, trs);
    }

    if (trs != NULL) {
	g_free(trs);
    }

    return 0;
}

/* list the topics available in the functions help file */

static int func_help_topics (const char *helpfile, PRN *prn)
{
    char line[128], word[12];
    FILE *fp;
    int j, k;

    if ((fp = gretl_fopen(helpfile, "r")) == NULL) {
	printf(_("Unable to access the file %s.\n"), helpfile);
	return E_FOPEN;
    }

    j = 1;
    k = 0;
    while (fgets(line, sizeof line, fp) != NULL) {
	if (!strncmp(line, "## ", 3)) {
	    /* sub-heading */
	    tailstrip(line);
	    if (k++ > 0) {
		pputc(prn, '\n');
	    }
	    pprintf(prn, "\n%s:\n", line + 3);
	    j = 1;
	} else if (*line == '#') {
	    /* actual entry */
	    sscanf(line + 2, "%10s", word);
	    pprintf(prn, "%-10s", word);
	    if (j % 7 == 0) {
		pputc(prn, '\n');
	    } else {
		pputc(prn, ' ');
	    }
	    j++;
	}
    }

    pputs(prn, _("\n\nFor help on a specific function, type: help funname"));
    pputs(prn, _(" (e.g. help qrdecomp)\n"));

    fclose(fp);

    return 0;
}

static void output_help_line (const char *line, PRN *prn, int recode)
{
    if (recode > 0) {
	recode_print_line(line, prn);
    } else {
	pputs(prn, line);
    }
}

/* check in the CLI helpfile for a line of the form "  @s: ...",
   where @s has been recognized as a libset variable */

static int got_setvar_line (const char *s, int n, char *line)
{
    if (!strncmp(line, "  ", 2) &&
	!strncmp(line + 2, s, n) &&
	line[2 + n] == ':') {
	return 1;
    } else {
	return 0;
    }
}

#define HELPLEN 128

/* special case: the user has done "help set @setvar" */

static int do_set_help (const char *setvar, FILE *fp,
			int recode, PRN *prn)
{
    char s[9], line[HELPLEN];
    int n = strlen(setvar);
    int ok = 0, count = 0;

    while (count < 2 && fgets(line, HELPLEN, fp) != NULL) {
	if (*line != '#') {
	    continue;
	}
	sscanf(line + 2, "%8s", s);
	if (!strcmp(s, "set")) {
	    while (count < 2 && fgets(line, HELPLEN, fp)) {
		if (got_setvar_line(setvar, n, line)) {
		    pputc(prn, '\n');
		    output_help_line(line + 2, prn, recode);
		    count++;
		} else if (count > 0) {
		    if (string_is_blank(line)) {
			/* reached the end of the item */
			pputc(prn, '\n');
			count++;
		    } else {
			output_help_line(line + 2, prn, recode);
		    }
		} else if (*line == '#') {
		    /* overshot: not found */
		    count = 999;
		}
	    }
	}
    }

    ok = count > 0 && count < 999;

    if (!ok) {
	pprintf(prn, _("%s: sorry, no help available.\n"), setvar);
    }

    return ok;
}

static int output_help_text (const char *cmd, FILE *fp,
			     int recode, PRN *prn)
{
    char word[16], line[HELPLEN];
    int ok = 0;

    while (fgets(line, sizeof line, fp) != NULL) {
	if (*line != '#') {
	    continue;
	}
	sscanf(line + 2, "%15s", word);
	if (!strcmp(cmd, word)) {
	    ok = 1;
	    pprintf(prn, "\n%s\n", word);
	    while (fgets(line, sizeof line, fp)) {
		if (*line == '#') {
		    break;
		}
		output_help_line(line, prn, recode);
	    }
	    break;
	}
    }

    if (!ok) {
	pprintf(prn, _("%s: sorry, no help available.\n"), cmd);
    }

    return ok;
}

static void do_help_on_help (PRN *prn)
{
    int i, j;

    pputs(prn, _("\nValid gretl commands are:\n"));
    j = 1;
    for (i=1; i<NC; i++) {
	if (HIDDEN_COMMAND(i)) {
	    continue;
	}
	pprintf(prn, "%-9s", gretl_command_word(i));
	if (j % 8 == 0) {
	    pputc(prn, '\n');
	} else {
	    pputc(prn, ' ');
	}
	j++;
    }

    pputs(prn, _("\n\nFor help on a specific command, type: help cmdname"));
    pputs(prn, _(" (e.g. help smpl)\n"));
    pputs(prn, _("You can also do 'help functions' for a list of functions\n"));
}

static int is_functions_dir (const char *path)
{
    int n = strlen(path) - 9;

    return n > 0 && !strcmp(path + n, "functions");
}

static int find_pkg_in_dir (const char *targ,
			    GDir *dir, const char *path,
			    char **gfn, char **pdf)
{
    const gchar *basename;
    char fullname[MAXLEN];
    char test[128];
    int found = 0;

    *fullname = '\0';

    /* first look for pdf or gfn in own subdir */

    if (is_functions_dir(path)) {
	while ((basename = g_dir_read_name(dir)) != NULL && !found) {
	    if (!strcmp(basename, ".") ||
		!strcmp(basename, "..")) {
		continue;
	    }
	    if (!strcmp(basename, targ)) {
		gretl_build_path(fullname, path, basename, NULL);
		if (gretl_isdir(fullname)) {
		    /* construct functions/foo/foo.pdf */
		    strcat(fullname, SLASHSTR);
		    strcat(fullname, basename);
		    strcat(fullname, ".pdf");
		    if (gretl_file_exists(fullname)) {
			*pdf = gretl_strdup(fullname);
			found = 1;
		    } else {
			switch_ext(fullname, fullname, "gfn");
			if (gretl_file_exists(fullname)) {
			    *gfn = gretl_strdup(fullname);
			    found = 1;
			}
		    }
		}
		gretl_error_clear();
	    }
	}
	g_dir_rewind(dir);
    }

    /* then look for "plain gfn" files */

    while (!found && (basename = g_dir_read_name(dir)) != NULL) {
	if (has_suffix(basename, ".gfn")) {
	    sprintf(test, "%s.gfn", targ);
	    if (!strcmp(basename, test)) {
		gretl_build_path(fullname, path, basename, NULL);
		*gfn = gretl_strdup(fullname);
		found = 1;
	    }
	}
    }

    return found;
}

static int show_pkg_pdf (const char *fname)
{
    int err = 0;

#if defined(G_OS_WIN32)
    if ((int) ShellExecute(NULL, "open", fname, NULL, NULL, SW_SHOW) <= 32) {
	err = E_FOPEN;
    }
#elif defined(__APPLE__)
    CFURLRef u;

    u = CFURLCreateFromFileSystemRepresentation(NULL,
                                                (const UInt8 *) fname,
                                                strlen(fname),
                                                false);
    if (u != NULL) {
        err = LSOpenCFURLRef(u, NULL);
        CFRelease(u);
    } else {
        err = 1;
    }
#else
    char *syscmd = g_strdup_printf("xdg-open \"%s\"", fname);

    err = system(syscmd);
    g_free(syscmd);
#endif

    return err;
}

static int find_function_package_help (const char *targ,
				       char **pbuf, PRN *prn,
                                       int *err)
{
    char *gfnname = NULL;
    char *pdfname = NULL;
    char **dnames = NULL;
    int i, n_dirs = 0;
    int found = 0;

    /* get names of directories to search */
    dnames = get_plausible_search_dirs(FUNCS_SEARCH, &n_dirs);

    for (i=0; i<n_dirs && !found; i++) {
	GDir *dir = gretl_opendir(dnames[i]);

	if (dir != NULL) {
	    found = find_pkg_in_dir(targ, dir, dnames[i],
                                    &gfnname, &pdfname);
	    g_dir_close(dir);
	}
    }

    if (pdfname != NULL) {
	*err = show_pkg_pdf(pdfname);
	free(pdfname);
    } else if (gfnname != NULL) {
	*err = print_function_package_help(gfnname, pbuf, prn);
	free(gfnname);
    }

    strings_array_free(dnames, n_dirs);

    return found;
}

static int maybe_do_markup (const char *hlpword, char **pbuf,
			    PRN *prn)
{
    PRN *hprn;
    int err = 0;

    hprn = gretl_print_new(GRETL_PRINT_BUFFER, &err);

    if (err) {
	return user_function_help(hlpword, OPT_NONE, prn);
    } else {
	err = user_function_help(hlpword, OPT_M, hprn);
	if (!err) {
	    *pbuf = gretl_print_steal_buffer(hprn);
	    gretl_print_destroy(hprn);
	}
	return err;
    }
}

/**
 * cli_help:
 * @hlpword: the word (usually a command) on which help is wanted.
 * @param: parameter usable when @hlpword is (literally) "set".
 * @opt: may include %OPT_F to give priority to functions
 * rather than commands.
 * @pbuf: location to receive text buffer with GUI-friendly markup,
 * if applicable.
 * @prn: pointer to gretl printing struct.
 *
 * Searches in the gretl help files for help on @hlpword and,
 * if help is found, prints it to @prn.  If @hlpword is %NULL,
 * lists the valid commands.
 *
 * Returns: 0 on success, non-zero if no help could be found.
 */

int cli_help (const char *hlpword, const char *param,
	      gretlopt opt, char **pbuf, PRN *prn)
{
    static int recode = -1;
    char helpfile[FILENAME_MAX];
    FILE *fp;
    int noword, funhelp = (opt & OPT_F);
    int ok = 0;
    int err = 0;

    noword = (hlpword == NULL || *hlpword == '\0');

    /* @param should be NULL unless hlpword is "set" */
    if ((noword || strcmp(hlpword, "set")) && param != NULL) {
	return 1;
    }

    if (noword && !funhelp) {
	do_help_on_help(prn);
	return 0; /* handled */
    }

    if ((noword && funhelp) || !strcmp(hlpword, "functions")) {
	strcpy(helpfile, helpfile_path(GRETL_FUNCREF, 1, 0));
	return func_help_topics(helpfile, prn); /* handled */
    }

    if (!funhelp) {
	int ci = gretl_help_index(hlpword);

	if (ci == SET && param != NULL) {
            /* respond to "help set <setvar>" */
	    ok = libset_help_available(param);
	    if (!ok) {
		pprintf(prn, _("%s: sorry, no help available.\n"), param);
		return 1;
	    }
	} else {
	    ok = ci > 0;
	}
    }

    if (ok) {
	strcpy(helpfile, helpfile_path(GRETL_CMDREF, 1, 0));
    } else if (genr_function_word(hlpword)) {
	strcpy(helpfile, helpfile_path(GRETL_FUNCREF, 1, 0));
    } else if (gretl_is_public_user_function(hlpword)) {
	if (pbuf != NULL) {
	    return maybe_do_markup(hlpword, pbuf, prn);
	} else {
	    return user_function_help(hlpword, OPT_NONE, prn);
	}
    } else if (find_function_package_help(hlpword, pbuf, prn, &err)) {
	return err; /* handled */
    } else {
	pprintf(prn, _("\"%s\" is not a gretl command.\n"), hlpword);
	return 1;
    }

    fp = gretl_fopen(helpfile, "r");

    if (fp == NULL) {
	printf(_("Unable to access the file %s.\n"), helpfile);
	err = E_FOPEN;
    } else {
	if (!gretl_in_gui_mode() && recode < 0) {
	    recode = maybe_need_recode();
	}
	/* actually output the relevant help text */
	if (param != NULL) {
	    ok = do_set_help(param, fp, recode, prn);
	} else {
	    ok = output_help_text(hlpword, fp, recode, prn);
	}

	fclose(fp);
    }

    return err;
}
