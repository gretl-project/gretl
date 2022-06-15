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

/* R-related functions that are not used by gretl_edit */

#ifdef G_OS_WIN32 /* specific to Windows */

static void win32_run_R_sync (const char *buf, gretlopt opt)
{
    PRN *prn = NULL;
    int err;

    if (bufopen(&prn)) {
	return;
    }

    err = execute_R_buffer(buf, dataset, opt, prn);

    if (err) {
	gui_errmsg(err);
    } else {
	view_buffer(prn, 78, 350, _("gretl: script output"),
		    PRINT, NULL);
    }
}

#else /* non-Windows functions */

/* Start an R session in asynchronous (interactive) mode.
   Note that there's a separate win32 function for this
   in gretlwin32.c. We don't do interactive when in editor
   mode.
*/

static void start_R_async (void)
{
    char *s0 = NULL, *s1 = NULL, *s2 = NULL;
    int n = -1;

    s0 = mymalloc(64);
    s1 = mymalloc(32);
    s2 = mymalloc(32);

    if (s0 != NULL && s1 != NULL && s2 != NULL) {
	*s0 = *s1 = *s2 = '\0';
	/* probably "xterm -e R" or similar */
	n = sscanf(Rcommand, "%63s %31s %31s", s0, s1, s2);
    }

    if (n == 0) {
	errbox(_("No command was supplied to start R"));
    } else if (n > 0) {
	char *supp1 = "--no-init-file";
	char *supp2 = "--no-restore-data";
	gchar *argv[6];
	GError *error = NULL;
	gboolean ok;
	int i = 0;

	argv[i++] = s0;
	if (n > 1) {
	    argv[i++] = s1;
	}
	if (n > 2) {
	    argv[i++] = s2;
	}
	argv[i++] = supp1;
	argv[i++] = supp2;
	argv[i++] = NULL;

	ok = g_spawn_async(NULL,
			   argv,
			   NULL,
			   G_SPAWN_SEARCH_PATH,
			   NULL,
			   NULL,
			   NULL,
			   &error);

	if (error != NULL) {
	    errbox(error->message);
	    g_error_free(error);
	} else if (!ok) {
	    gui_errmsg(E_EXTERNAL);
	    g_error_free(error);
	}
    }

    free(s0);
    free(s1);
    free(s2);
}

static void run_R_sync (void)
{
    gchar *argv[] = {
	"R",
	"--no-save",
	"--no-init-file",
	"--no-restore-data",
	"--slave",
	NULL
    };

    run_prog_sync(argv, LANG_R);
}

#endif /* Windows vs non-Windows, common code follows */

void start_R (const char *buf, int send_data, int interactive)
{
    gretlopt Ropt = OPT_G;
    int err = 0;

    if (send_data && !data_status) {
	warnbox(_("Please open a data file first"));
	return;
    }

    if (interactive) {
	Ropt |= OPT_I;
    }

    if (send_data) {
	Ropt |= OPT_D;
	if (annual_data(dataset) || quarterly_or_monthly(dataset)) {
	    const char *opts[] = {
		N_("multiple time series object"),
		N_("data frame")
	    };
	    int resp;

	    resp = radio_dialog(NULL, _("Send data as"), opts, 2, 0, 0, NULL);
	    if (resp < 0) {
		return;
	    } else if (resp == 1) {
		Ropt |= OPT_F;
	    }
	}
    }

    /* On Windows in non-interactive mode, don't write
       these files here; that will be handled later
    */
#ifdef G_OS_WIN32
    if (interactive) {
	err = write_gretl_R_files(buf, dataset, Ropt);
    }
#else
    err = write_gretl_R_files(buf, dataset, Ropt);
#endif

    if (err) {
	gui_errmsg(err);
	delete_gretl_R_files();
    } else if (interactive) {
#ifdef G_OS_WIN32
	win32_start_R_async();
#else
	start_R_async();
#endif
    } else {
	/* non-interactive */
#ifdef G_OS_WIN32
	win32_run_R_sync(buf, Ropt);
#else
	run_R_sync();
#endif
    }
}

