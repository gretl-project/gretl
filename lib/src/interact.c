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

/* interact.c for gretl */

#include "libgretl.h"
#include "monte_carlo.h"
#include "var.h"
#include "johansen.h"
#include "gretl_func.h"
#include "compat.h"
#include "system.h"
#include "forecast.h"
#include "cmd_private.h"
#include "libset.h"
#include "uservar.h"
#include "gretl_panel.h"
#include "texprint.h"
#include "gretl_join.h"
#include "gretl_xml.h"
#include "gretl_string_table.h"
#include "gretl_typemap.h"
#include "gretl_midas.h"
#include "dbread.h"
#include "gretl_foreign.h"
#include "boxplots.h"
#include "gretl_plot.h"
#include "flow_control.h"
#include "gretl_drivers.h"
#include "csvdata.h"
#include "gretl_zip.h"
#include "matrix_extra.h"
#include "addons_utils.h"
#include "gretl_gridplot.h"
#include "gretl_sampler.h"
#include "gretl_untar.h"
#ifdef USE_CURL
# include "gretl_www.h"
#endif
#ifdef HAVE_MPI
# include "gretl_mpi.h"
#endif

#include <errno.h>

/* for the "shell" command */
#ifdef WIN32
# include "gretl_win32.h"
#else
# ifdef HAVE_PATHS_H
#  include <paths.h>
# endif
#endif

#define CMD_DEBUG 0
#define ECHO_DEBUG 0

#include "tokenize.c"

#define bare_quote(p,s) (*p == '"' && (p-s==0 || *(p-1) != '\\'))
#define starts_ccmt(p)  (*p == '/' && *(p+1) == '*')
#define ends_ccmt(p)    (*p == '*' && *(p+1) == '/')

static int install_package (const char *pkgname,
			    gretlopt opt,
			    ExecState *s,
			    PRN *prn);

static int run_script (const char *fname, ExecState *s,
                       DATASET *dset, gretlopt opt,
                       PRN *prn);

static int strip_inline_comments (char *s)
{
    char *p = strchr(s, '#');

    /* We return 1 if and only if the entire line is a #-comment;
       otherwise we try to ensure that if '#' occurs in its role
       as comment-opener, the comment is trimmed off. The only
       cases where '#' is not a comment-opener are when it
       occurs in a string literal, or inside "{...}" in the
       supplement to a plotting command.
    */

    if (p - s == 0) {
        /* the entire line is a comment */
        return 1;
    } else if (p == NULL) {
        /* no '#' in line */
        return 0;
    }

    if (strchr(p+1, '"') == NULL &&
        strchr(p+1, '}') == NULL) {
        /* '#' must start a comment */
        *p = '\0';
    } else {
        int quoted = 0;
        int braced = 0;

        p = s;
        while (*p) {
            if (bare_quote(p, s)) {
                quoted = !quoted;
            } else if (!quoted) {
                if (*p == '{') {
                    braced++;
                } else if (*p == '}') {
                    braced--;
                }
            }
            if (!quoted && !braced) {
                if (*p == '#') {
                    *p = '\0';
                    break;
                }
            }
            p++;
        }
    }

    return 0;
}

static void toggle_ccmt (CMD *cmd, int state)
{
    if (state) {
        cmd->flags |= CMD_CCMT;
    } else {
        cmd->flags &= ~CMD_CCMT;
    }
}

static int tail_is_blank (const char *s)
{
    while (*s) {
	if (!isspace(*s)) {
	    return 0;
	}
	s++;
    }

    return 1;
}

/* filter_comments: strip comments out of line; return non-zero if
   the whole line is a comment
*/

static int filter_comments (char *s, CMD *cmd, int preserve)
{
    char tmp[MAXLINE];
    char *p = s;
    int ccmt, quoted = 0;
    int j = 0, filt = 0;

    if (strlen(s) >= MAXLINE) {
        cmd->err = E_TOOLONG;
        return 0;
    }

    ccmt = (cmd->flags & CMD_CCMT);

    while (*p) {
        if (!quoted && !ccmt && *p == '#') {
            break;
        }
        if (!ccmt && bare_quote(p, s)) {
            quoted = !quoted;
        }
        if (!quoted) {
            if (starts_ccmt(p)) {
                ccmt = 1;
                p += 2;
            } else if (ends_ccmt(p)) {
                if (!ccmt) {
                    cmd->err = E_PARSE;
                    fprintf(stderr, "unbalanced close comment\n");
                    return 0;
                }
                ccmt = 0;
                p += 2;
                p += strspn(p, " ");
            }
        }
        if (!ccmt && *p != '\r') {
            tmp[j++] = *p;
        }
        if (*p) {
            p++;
        }
    }

    tmp[j] = '\0';
    if (preserve && tail_is_blank(tmp)) {
        cmd->ci = CMD_COMMENT;
        toggle_ccmt(cmd, ccmt);
        return 1;
    }
    strcpy(s, tmp);
    tailstrip(s);

    if (*s == '\0') {
        filt = 1;
    } else if (!ccmt) {
        /* '#' comments */
        filt = strip_inline_comments(s);
        tailstrip(s);
    }

    if (filt) {
        /* the whole line is a comment */
        cmd->ci = CMD_COMMENT;
    }

    toggle_ccmt(cmd, ccmt);

    return filt;
}

#define MODIFIES_LIST(c) (c == DIFF ||          \
                          c == DUMMIFY ||       \
                          c == LDIFF ||         \
                          c == SDIFF ||         \
                          c == LAGS ||          \
                          c == LOGS ||          \
                          c == SQUARE ||        \
                          c == ORTHDEV ||       \
                          c == STDIZE)

static int has_param (const CMD *cmd)
{
    return cmd->param != NULL && *cmd->param != '\0';
}

/* Look for a line with an "implicit genr", such as
   y = 3*x, x += 10, etc. This is used in nls.c to
   assess auxiliary genrs in nls, mle, gmm.
*/

int plausible_genr_start (const char *s, const DATASET *dset)
{
    int ret = 0;

    if (strchr(s, '=') || strstr(s, "++") || strstr(s, "--")) {
        const char *ok = ".+-*/%^~|=[";
        char word[VNAMELEN] = {0};
        char fmt[20];

        sprintf(fmt, "%%%d[^[ .+*/%%^~|=-]", VNAMELEN - 1);

        if (sscanf(s, fmt, word)) {
            s += strlen(word);
            while (*s == ' ') s++;
            if (strspn(s, ok) > 0 && check_identifier(word) == 0) {
                ret = 1;
            }
        }
    } else if (gretl_type_from_name(s, dset) != 0) {
        ret = 1;
    }

    return ret;
}

static int ends_foreign_block (const char *s)
{
    s += strspn(s, " \t");

    if (!strncmp(s, "end ", 4)) {
        s += 3;
        s += strspn(s, " \t");
        if (!strncmp(s, "foreign", 7)) {
            return 1;
        } else if (!strncmp(s, "mpi", 3)) {
            return 1;
        }
    }

    return 0;
}

/**
 * parse_command_line:
 * @s: pointer to execution-state struct.
 * @dset: dataset struct.
 * @ptr: pointer for use with "compilation" of
 * conditionals.
 *
 * Parses @line and fills out @cmd accordingly.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int parse_command_line (ExecState *s, DATASET *dset, void *ptr)
{
    char *line = s->line;
    CMD *cmd = s->cmd;

    gretl_cmd_clear(cmd);
    gretl_error_clear();

#if CMD_DEBUG
    fprintf(stderr, "\nparse_command_line: '%s' (nosub=%d)\n",
            line, cmd_nosub(cmd) ? 1 : 0);
#endif

    if (cmd_nosub(cmd)) {
        cmd->flags &= ~CMD_SUBST;
    } else {
        int subst = 0;

        cmd->err = substitute_named_strings(line, &subst);
        if (cmd->err) {
            return cmd->err;
        } else if (subst) {
            /* record the fact that substitution has been done */
            cmd->flags |= CMD_SUBST;
        } else {
            cmd->flags &= ~CMD_SUBST;
        }
    }

#if CMD_DEBUG
    if (cmd->flags & CMD_SUBST) {
        fprintf(stderr, "after substitution: '%s'\n", line);
    }
#endif

    if ((cmd->context == FOREIGN || cmd->context == MPI) &&
        !ends_foreign_block(line)) {
        cmd->opt = OPT_NONE;
        cmd->ci = cmd->context;
        return 0;
    }

    if ((cmd->flags & CMD_SUBST) || !gretl_looping_currently()) {
        /* normalize line spaces */
        compress_spaces(line);

        /* trap lines that are nothing but comments */
        if (filter_comments(line, cmd, 0)) {
            return 0;
        }

        /* catch errors associated with comment syntax */
        if (cmd->err) {
            return cmd->err;
        }
    }

    cmd->err = real_parse_command(s, dset, 0, ptr);

    if (cmd->err) {
        gretl_cmd_destroy_context(cmd);
    }

    return cmd->err;
}

#ifndef WIN32

static int gretl_shell_async (const char *cmdline, PRN *prn)
{
    GError *gerr = NULL;
    int err = 0;

    g_spawn_command_line_async(cmdline, &gerr);

    if (gerr != NULL) {
        pprintf(prn, "%s\n", gerr->message);
        g_error_free(gerr);
        err = 1;
    }

    return err;
}

#define SHELL_KISS 1 /* Keep It Simple, Stupid */

static int gretl_shell_sync (const char *arg, gchar **psout,
                             PRN *prn)
{
    gchar *sout = NULL;
    gchar *serr = NULL;
    GError *gerr = NULL;
    gchar *argv[5];
    int i, status;
    int err = 0;

#if SHELL_KISS
    argv[0] = g_strdup("/bin/sh");
    argv[1] = g_strdup("-sh");
    argv[2] = g_strdup("-c");
    argv[3] = g_strdup(arg);
    argv[4] = NULL;
#else
    /* the following fails with zsh on macOS */
    const char *theshell = getenv("SHELL");
    const char *namep;
    char shellnam[40];

    if (theshell == NULL) {
# ifdef HAVE_PATHS_H
        theshell =_PATH_BSHELL;
# else
        theshell = "/bin/sh";
# endif
    }

    namep = strrchr(theshell, '/');
    if (namep == NULL) {
        namep = theshell;
    }

    strcpy(shellnam, "-");
    strcat(shellnam, ++namep);
    if (strcmp(namep, "sh") != 0) {
        shellnam[0] = '+';
    }

    argv[0] = g_strdup(theshell);
    argv[1] = g_strdup(shellnam);
    argv[2] = g_strdup("-c");
    argv[3] = g_strdup(arg);
    argv[4] = NULL;
#endif /* SHELL_KISS or not */

    g_spawn_sync(gretl_workdir(), argv, NULL, 0, NULL, NULL,
                 &sout, &serr, &status, &gerr);

    for (i=0; i<4; i++) {
	g_free(argv[i]);
    }

    if (gerr != NULL) {
        if (prn != NULL) {
            pprintf(prn, "%s\n", gerr->message);
        } else {
            gretl_errmsg_set(gerr->message);
        }
        g_error_free(gerr);
        err = 1;
    }

    if (psout != NULL) {
        *psout = sout;
    } else if (sout != NULL) {
        pputs(prn, sout);
        g_free(sout);
    }

    if (serr != NULL) {
        pputs(prn, serr);
        g_free(serr);
    }

    return err;
}

/**
 * gretl_shell_grab:
 * @arg: command line to be executed.
 * @sout: location to receive output from command.
 *
 * Calls the shell to execute @arg syncronously and captures the
 * standard output, if any, in @sout.
 *
 * Returns: 0 on successful completion, non-zero on error.
 */

int gretl_shell_grab (const char *arg, char **sout)
{
    /* note: the win32 implementation of gretl_shell_grab()
       is defined in gretl_win32.c
    */
    return gretl_shell_sync(arg, sout, NULL);
}

static int gretl_shell (const char *arg, gretlopt opt, PRN *prn)
{
    int err = 0;

    if (arg == NULL || *arg == '\0') {
        return 0;
    }

    if (!libset_get_bool(SHELL_OK)) {
        gretl_errmsg_set(_("The shell command is not activated."));
        return 1;
    }

    arg += strspn(arg, " \t");

    if (opt & OPT_A) {
        /* "launch" */
        err = gretl_shell_async(arg, prn);
    } else {
        err = gretl_shell_sync(arg, NULL, prn);
    }

    return err;
}

#endif /* ! WIN32 */

#define SAFELEN 78 /* ? */

static void trim_to_length (char *s)
{
    int i, n = strlen(s);

    if (n < SAFELEN - 1) return;

    for (i=n-1; i>0; i--) {
        if (s[i] == ' ') {
            s[i] = '\0';
            break;
        }
    }
}

void safe_print_line (const char *line, int *plen, PRN *prn)
{
    char tmp[SAFELEN];
    const char *q, *p = line;
    int n, m, rem, out = 0;
    int len0 = *plen;

    rem = n = strlen(line);

    while (out < n) {
        *tmp = 0;
        q = p;
        strncat(tmp, p, SAFELEN - 1);
        len0 = 0;
        trim_to_length(tmp - len0);
        len0 = 0;
        m = strlen(tmp);
        out += m;
        rem = n - out;
        p = q + m;
        if (rem > 0) {
            pprintf(prn, "%s \\\n ", tmp);
            *plen = 1;
        } else {
            pprintf(prn, "%s", tmp);
            *plen += m;
        }
    }
}

static void new_trim_to_length (char *s, int len)
{
    int n = strlen(s);

    if (n > len) {
        int i, quoted = 0;
        int bp0 = 0, bp1 = 0;

        for (i=1; i<n-1; i++) {
            if (s[i] == '"' && s[i-1] != '\\') {
                quoted = !quoted;
            }
            if (!quoted && s[i] == ' ') {
                if (i < len) {
                    bp0 = i;
                } else {
                    bp1 = i;
                    break;
                }
            }
        }
        if (bp0 > 0) {
            s[bp0] = '\0';
        } else if (bp1 > 0) {
            s[bp1] = '\0';
        }
    }
}

static void basic_trim_to_length (char *s, int len)
{
    int n = strlen(s);

    if (n > len) {
        int i;
        int bp0 = 0, bp1 = 0;

        for (i=1; i<n-1; i++) {
            if (s[i] == ' ') {
                if (i < len) {
                    bp0 = i;
                } else {
                    bp1 = i;
                    break;
                }
            }
        }
        if (bp0 > 0) {
            s[bp0] = '\0';
        } else if (bp1 > 0) {
            s[bp1] = '\0';
        }
    }
}

#define TESTLEN 256
#define LINELEN 70

static void reflow_line (const char *line, const CMD *cmd,
                         const char *leader, PRN *prn)
{
    int maxline = LINELEN;

    if (leader != NULL) {
        maxline -= 2;
        pputs(prn, leader);
    }

    if (cmd != NULL && (cmd->ciflags & CI_EXPR)) {
        /* "genr"-type lines: be more generous? */
        maxline += 10;
    } else if (gretl_in_gui_mode()) {
        /* we can handle a little more width */
        maxline += 4;
    }

    if (strlen(line) < maxline) {
        pputs(prn, line);
    } else {
        const char *p = line;
        char buf[TESTLEN];
        int linenum = 0;

        while (*p) {
            *buf = '\0';
            strncat(buf, p, TESTLEN - 1);
            if (linenum > 0 && leader == NULL) {
                new_trim_to_length(buf, maxline - 2);
            } else {
                new_trim_to_length(buf, maxline);
            }
            p += strlen(buf);
            if (!string_is_blank(buf)) {
                if (linenum > 0) {
                    pputs(prn, "  ");
                }
                pputs(prn, (*buf == ' ')? buf + 1 : buf);
                if (*p) {
                    pputs(prn, " \\\n");
                }
            }
            linenum++;
        }
    }
}

static int command_is_silent (const CMD *cmd, const char *line)
{
    if (cmd == NULL) {
        return 0;
    }

    if (cmd->ci == FUNCERR || cmd->ci == PRINTF ||
        (cmd->ci == PRINT && strchr(line, '"'))) {
        return 1;
    }

    if (!strcmp(line, "set echo off") ||
        !strcmp(line, "set verbose off") ||
        !strcmp(line, "flush")) {
        return 1;
    }

    if (!strncmp(line, "quit", 4) && string_is_blank(line + 4)) {
        return 1;
    }

    if (cmd->ci == SET && cmd->param != NULL &&
        !strcmp(cmd->param, "echo") &&
        gretl_function_depth() > 0) {
        return 1;
    }

    if (cmd->ci == END && cmd->param != NULL &&
	!strcmp(cmd->param, "outfile")) {
        return 1;
    }

    if (*line == '!') {
        return 1;
    }

    return 0;
}

/*
 * real_echo_command:
 * @cmd: pointer to #CMD struct.
 * @line: "raw" command line associated with @cmd.
 * @recording: echo is going to command log (0/1).
 * @prn: pointer to gretl printing struct.
 *
 * Echoes the user command represented by @cmd and @line to
 * @prn.  This is used for two distinct purposes: to give
 * visual feedback on the command supplied, and (in some
 * contexts) to record a command that was executed interactively.
 */

static void real_echo_command (CMD *cmd, const char *line,
                               int recording, PRN *prn)
{
    const char *leader = NULL;
    int commented_store = 0;
    int compiling = 0;

    if (line == NULL || *line == '\0' || prn == NULL) {
        return;
    }

    if (cmd != NULL && cmd->ci >= NC) {
        return;
    }

    if (gretl_compiling_function() || gretl_compiling_loop()) {
        compiling = 1;
    }

#if ECHO_DEBUG
    if (cmd != NULL) {
        fprintf(stderr, "echo_cmd:\n*** line='%s'\n param='%s' parm2='%s'\n",
                line, cmd->param, cmd->parm2);
        fprintf(stderr, " cmd->opt=%d, recording=%d, compiling=%d\n",
                cmd->opt, recording, compiling);
        fprintf(stderr, " cmd->ci = %d (%s), context = %d\n", cmd->ci,
                gretl_command_word(cmd->ci), cmd->context);
        fprintf(stderr, " cmd->savename = '%s'\n", cmd->savename);
        if (cmd->list != NULL) {
            printlist(cmd->list, "cmd->list");
        }
    }
#endif

    /* certain things don't get echoed at all, if not recording or
       compiling a function or loop */
    if (!recording && !compiling && command_is_silent(cmd, line)) {
#if ECHO_DEBUG
        fprintf(stderr, " silent: no echo\n");
#endif
        return;
    }

    /* print leading string before echo? */
    if (recording) {
        if (cmd != NULL && cmd->ci == STORE) {
            commented_store = 1;
        }
    } else if (compiling) {
        leader = "> ";
    } else {
        leader = "? ";
    }

    if (commented_store) {
        pputs(prn, "# ");
        pputs(prn, line);
    } else if (cmd != NULL && (cmd->context == FOREIGN || cmd->context == MPI)) {
        if (leader != NULL) {
            pputs(prn, leader);
        }
        pputs(prn, line);
    } else {
        reflow_line(line, cmd, leader, prn);
    }

    pputc(prn, '\n');

    gretl_print_flush_stream(prn);
}

void gretl_echo_command (CMD *cmd, const char *line, PRN *prn)
{
    real_echo_command(cmd, line, 0, prn);
}

void gretl_record_command (CMD *cmd, const char *line, PRN *prn)
{
    real_echo_command(cmd, line, 1, prn);
}

static int set_var_info (const int *list,
                         const char *parm1,
                         const char *parm2,
                         gretlopt opt,
                         DATASET *dset)
{
    int vi, v = list[1];
    int i, err = 0;

    if (dset == NULL || dset->varinfo == NULL) {
        return E_NODATA;
    } else if (v <= 0 || v >= dset->v) {
        return E_DATA;
    }

    if (opt & OPT_M) {
        err = gretl_list_set_midas(list, dset);
        if (err) {
            return err;
        }
    }

    for (i=1; i<=list[0]; i++) {
        vi = list[i];
        if (opt & OPT_D) {
            series_set_discrete(dset, vi, 1);
        } else if (opt & OPT_C) {
            series_set_discrete(dset, vi, 0);
        }
        if (opt & OPT_F) {
            /* --coded */
            int ivals = series_is_integer_valued(dset, vi);
            int isdum = gretl_isdummy(0, dset->n - 1, dset->Z[vi]);

            if (ivals && !isdum) {
                series_set_flag(dset, vi, VAR_CODED);
            } else {
                gretl_errmsg_sprintf(_("%s cannot be set as 'coded' (%s)"),
                                     dset->varname[vi], ivals ?
                                     _("is a 0/1 variable") :
                                     _("not integer-valued"));
                err = E_TYPES;
            }
        } else if (opt & OPT_N) {
            /* -- numeric (i.e. not-coded) */
            series_unset_flag(dset, vi, VAR_CODED);
        }
    }

    if (err) {
        return err;
    }

    /* below: we'll accept multi-series lists, but the
       string-setting facility will apply to just the
       first member, as "representative" of the list
    */

    if (opt & OPT_I) {
        const char *s = get_optval_string(SETINFO, OPT_I);

        if (s == NULL) {
            err = E_ARGS;
        } else {
            series_record_label(dset, v, s);
        }
    } else if (parm1 != NULL) {
        /* backward compatibility */
        series_record_label(dset, v, parm1);
    }

    if (opt & OPT_G) {
        const char *s = get_optval_string(SETINFO, OPT_G);

        if (s == NULL) {
            err = E_ARGS;
        } else {
            series_record_display_name(dset, v, s);
        }
    } else if (parm2 != NULL) {
        /* backward compatibility */
        series_record_display_name(dset, v, parm2);
    }

    return err;
}

static void reflow_label (const char *line, PRN *prn)
{
    int maxline = 72;

    if (strlen(line) < maxline) {
        pputc(prn, ' ');
        pputs(prn, line);
        pputc(prn, '\n');
    } else {
        const char *p = line;
        char buf[TESTLEN];
        int lnum = 0;

        while (*p) {
            *buf = '\0';
            strncat(buf, p, TESTLEN - 1);
            if (lnum == 1) {
                maxline -= 2;
            }
            basic_trim_to_length(buf, maxline);
            p += strlen(buf);
            if (!string_is_blank(buf)) {
                if (lnum == 0) {
                    pputc(prn, ' ');
                } else {
                    pputs(prn, "   ");
                }
                pputs(prn, (*buf == ' ')? buf + 1 : buf);
                pputc(prn, '\n');
            } else {
                pputc(prn, '\n');
            }
            lnum++;
        }
    }
}

static void showlabels (const int *list, gretlopt opt,
                        const DATASET *dset, PRN *prn)
{
    const char *label;
    gchar *tmp;
    int i, v, vmax, nl = 0;

    if (dset == NULL || dset->v == 0) {
        pprintf(prn, _("No series are defined\n"));
        return;
    }

    vmax = list == NULL ? dset->v - 1 : list[0];

    for (i=1; i<=vmax; i++) {
        v = list == NULL ? i : list[i];
        if (v >= 0 && v < dset->v) {
            label = series_get_label(dset, v);
            if (label != NULL && *label != '\0') {
                nl++;
            }
        }
    }

    if (nl == 0) {
        pprintf(prn, _("No series labels are defined\n"));
        return;
    }

    pputc(prn, '\n');
    for (i=1; i<=vmax; i++) {
        v = list == NULL ? i : list[i];
        if (v >= 0 && v < dset->v) {
            label = series_get_label(dset, v);
            if (label != NULL && *label != '\0') {
                if (opt & OPT_Q) {
                    pprintf(prn, "%s: %s\n", dset->varname[v], label);
                } else {
                    tmp = g_strdup_printf("%s: %s", dset->varname[v], label);
                    reflow_label(tmp, prn);
                    g_free(tmp);
                }
            }
        }
    }
    pputc(prn, '\n');
}

static int cwd_is_workdir (void)
{
    gchar *thisdir = g_get_current_dir();
    int ret = 0;

    if (thisdir != NULL) {
        int n = strlen(thisdir);

        ret = (strncmp(thisdir, gretl_workdir(), n) == 0);
        g_free(thisdir);
    }

    return ret;
}

/* elements of "outfile" state */
static char of_name[MAXLEN];
/* allow for 7 levels of redirection */
#define OF_MAX 7
static guint8 of_parms[OF_MAX+1];

enum {
    OF_ECHO = 1 << 0, /* echo was on before redirection */
    OF_MSGS = 1 << 1, /* messages were on before redirection */
    OF_NODP = 1 << 2  /* decimal comma was on before redirection */
};

static int outfile_redirect (PRN *prn, FILE *fp, const char *strvar,
                             const char *fname, gretlopt opt)
{
    int r = print_redirection_level(prn);
    int err;

    err = print_start_redirection(prn, fp, fname, strvar);

    if (!err) {
	of_parms[r] = 0;
	if (opt & OPT_Q) {
	    guint8 eo = gretl_echo_on();
	    guint8 mo = gretl_messages_on();

	    if (eo || mo) {
		if (eo) {
		    of_parms[r] |= OF_ECHO;
		}
		if (mo) {
		    of_parms[r] |= OF_MSGS;
		}
		set_gretl_echo(0);
		set_gretl_messages(0);
	    }
	}
	if (opt & OPT_D) {
	    int fd = libset_get_bool(FORCE_DECPOINT);
	    int dp = get_local_decpoint();

	    if (!fd && dp == ',') {
		of_parms[r] |= OF_NODP;
		libset_set_bool(FORCE_DECPOINT, 1);
	    }
	}
    }

    return err;
}

static void maybe_restore_outfile_parms (PRN *prn)
{
    int r = print_redirection_level(prn);

    if (of_parms[r] & OF_ECHO) {
        set_gretl_echo(1);
    }
    if (of_parms[r] & OF_MSGS) {
        set_gretl_messages(1);
    }
    if (of_parms[r] & OF_NODP) {
	libset_set_bool(FORCE_DECPOINT, 0);
    }
    of_parms[r] = 0;
}

static int redirection_ok (PRN *prn)
{
    int r = print_redirection_level(prn);
    int fd = gretl_function_depth();

    if ((fd == 0 && r > 0) || print_redirected_at_level(prn, fd)) {
	gretl_errmsg_set(_("Output has already been redirected"));
        return 0;
    } else if (r == OF_MAX) {
	gretl_errmsg_set(_("Output redirection: maximum depth reached"));
	return 0;
    } else {
        return 1;
    }
}

int check_stringvar_name (const char *name, int allow_new,
			  const DATASET *dset)
{
    GretlType t = gretl_type_from_name(name, dset);
    int err = 0;

    if (t != GRETL_TYPE_NONE && t != GRETL_TYPE_STRING) {
        gretl_errmsg_sprintf(_("'%s' is of type %s"), name,
                             gretl_type_get_name(t));
        err = E_TYPES;
    } else if (t == GRETL_TYPE_NONE) {
        if (allow_new) {
            /* create the variable if possible */
            err = check_identifier(name);
            if (!err) {
                err = user_var_add(name, GRETL_TYPE_STRING, NULL);
            }
        } else {
            /* otherwise require that the variable already exists */
            gretl_errmsg_sprintf(_("'%s' : not a string variable"), name);
            err = E_DATA;
        }
    }

    return err;
}

#define TMPFILE_DEBUG 0

/* We come here in the --tempfile and --buffer cases of "outfile". The
   @strvar argument, which names a string variable, plays a different
   role in each case:

   * With --tempfile (OPT_T) @strvar should hold the name of the
   temporary file to be created.

   * With --buffer (OPT_B) @strvar should hold the name of the string
   variable that gets the _content_ of the tempfile when the
   restriction ends.

   We accomplish the tempfile effect here, but arrange for the buffer
   effect by passing @strvar to outfile_redirect().
*/

static int redirect_to_tempfile (const char *strvar, PRN *prn,
                                 gretlopt opt)
{
    gchar *tempname = NULL;
    FILE *fp = NULL;
    int err = 0;

#if TMPFILE_DEBUG
    fprintf(stderr, "redirect_to_tempfile, strvar '%s'\n", strvar);
#endif

    if (opt & OPT_T) {
	const char *s = get_string_by_name(strvar);

	if (s != NULL && strstr(s, "XXXXXX")) {
	    tempname = gretl_make_dotpath(s);
	}
    }
    if (tempname == NULL) {
	tempname = gretl_make_dotpath("outfile.XXXXXX");
    }

    if (opt & OPT_B) {
        fp = gretl_mktemp(tempname, "wb+");
    } else {
        fp = gretl_mktemp(tempname, "wb");
    }

#if TMPFILE_DEBUG
    fprintf(stderr, " tempname = '%s', fp %p\n", tempname, (void *) fp);
#endif

    if (fp == NULL) {
        err = E_FOPEN;
    } else if (opt & OPT_B) {
	/* the buffer variant */
        err = outfile_redirect(prn, fp, strvar, tempname, opt);
    } else {
	/* the explicit tempfile variant */
        err = outfile_redirect(prn, fp, NULL, tempname, opt);
    }
    if (!err && (opt & OPT_T)) {
        /* write @tempname into @strvar */
        user_string_reset(strvar, tempname, &err);
        if (err) {
            fclose(fp);
        }
    }

    g_free(tempname);

    return err;
}

/* Note, 2024-03-22: only one element of the old outfile syntax
   (in place prior to gretl 2018d) is still supported, namely

   outfile <strname> --buffer

   where the specified string is created automatically if it's
   not already present. This is found in some older function
   packages.
*/

static int
do_outfile_command (gretlopt opt, const char *fname,
                    const DATASET *dset, PRN *prn)
{
    const char *strvar = NULL;
    int err = 0;

    if (prn == NULL) {
	gretl_errmsg_set(_("output cannot be redirected in this context"));
        return E_DATA;
    }

    /* options: allow at most one of --append, --buffer, --tempfile */
    err = incompatible_options(opt, (OPT_A | OPT_B | OPT_T));
    if (err) {
        return err;
    }

    /* check for invalid nesting of "outfile", or hitting max depth */
    if (!redirection_ok(prn)) {
        return 1;
    }

    if (opt & (OPT_B | OPT_T)) {
	/* handle the --buffer and --tempfile cases */
	if (opt & OPT_B) {
	    strvar = get_optval_string(OUTFILE, OPT_B);
	    if (strvar == NULL) {
		/* backward compatibility */
		strvar = fname;
	    }
	} else {
	    strvar = get_optval_string(OUTFILE, OPT_T);
	}
	if (strvar == NULL) {
	    return E_ARGS;
	}
	err = check_stringvar_name(strvar, 1, dset);
        if (!err) {
            err = redirect_to_tempfile(strvar, prn, opt);
        }
        *of_name = '\0';
        return err; /* we're done */
    }

    /* Handle the remaining cases, all of which need @fname, either
       a genuine filename or a dummy constant such as "stderr".
    */

    if (fname == NULL || *fname == '\0') {
	return E_ARGS;
    }

    if (!strcmp(fname, "null")) {
        if (gretl_messages_on()) {
            pputs(prn, _("Now discarding output\n"));
        }
        err = outfile_redirect(prn, NULL, NULL, fname, opt);
        *of_name = '\0';
    } else if (!strcmp(fname, "stderr")) {
        err = outfile_redirect(prn, stderr, NULL, fname, opt);
        *of_name = '\0';
    } else if (!strcmp(fname, "stdout")) {
        err = outfile_redirect(prn, stdout, NULL, fname, opt);
        *of_name = '\0';
    } else {
        char outname[FILENAME_MAX];
        const char *targ;
        FILE *fp;

        /* switch to workdir if needed */
        strcpy(outname, fname);
        gretl_maybe_prepend_dir(outname);
        if (opt & OPT_A) {
            /* appending to a file */
            fp = gretl_fopen(outname, "ab");
        } else {
            /* (over-)writing a file */
            fp = gretl_fopen(outname, "wb");
        }
        if (fp == NULL) {
            pprintf(prn, _("Couldn't open %s for writing\n"), outname);
            return E_FOPEN;
        }

        /* string to identify the output stream for display */
        targ = cwd_is_workdir() ? fname : outname;

        if (gretl_messages_on()) {
            /* print message before actual redirection! */
            if (opt & OPT_A) {
                pprintf(prn, _("Now appending output to '%s'\n"), targ);
            } else {
                pprintf(prn, _("Now writing output to '%s'\n"), targ);
            }
        }

        err = outfile_redirect(prn, fp, NULL, targ, opt);
        if (err) {
            fclose(fp);
            remove(outname);
        } else {
            strcpy(of_name, targ);
        }
    }

    return err;
}

static int close_outfile (PRN *prn)
{
    int rlevel = print_redirection_level(prn);
    int err = 0;

    if (rlevel == 0) {
	pputs(prn, _("Output is not currently diverted to file\n"));
	err = 1;
    } else {
	print_end_redirection(prn);
	maybe_restore_outfile_parms(prn);
	if (gretl_messages_on() && *of_name != '\0') {
	    pprintf(prn, _("Closed output file '%s'\n"), of_name);
	}
    }

    return err;
}

int call_pca_plugin (VMatrix *cmat, DATASET *dset,
                     gretlopt opt, PRN *prn)
{
    int (*pca_from_cmatrix) (VMatrix *, DATASET *,
                             gretlopt, PRN *);

    gretl_error_clear();

    pca_from_cmatrix = get_plugin_function("pca_from_cmatrix");
    if (pca_from_cmatrix == NULL) {
        return 1;
    }

    return (*pca_from_cmatrix) (cmat, dset, opt, prn);
}

static int do_pca (int *list, DATASET *dset, gretlopt opt, PRN *prn)
{
    int freelist = 0;
    int err = 0;

    if (list != NULL && list[0] == 0) {
        return 0;
    }

    if (list == NULL) {
        list = full_var_list(dset, NULL);
        freelist = 1;
    }

    if (list != NULL) {
        VMatrix *cmat = NULL;

        /* adding OPT_N ensures a uniform sample for the correlation
           or covariance matrix
	*/
        cmat = corrlist(PCA, list, dset, opt | OPT_N, &err);
        if (!err) {
            err = call_pca_plugin(cmat, dset, opt, prn);
            if (!err && (opt & (OPT_O | OPT_A))) {
                /* results saved as series */
                if (gretl_messages_on()) {
                    pputs(prn, _("Generated principal component series\n"));
                }
            }
            free_vmatrix(cmat);
        }
        if (freelist) {
            free(list);
        }
    }

    return err;
}

static void query_package (const char *pkgname,
                           gretlopt opt, PRN *prn)
{
    char path[MAXLEN];
    fnpkg *pkg = NULL;
    int err = 0;

    pkg = get_function_package_by_name(pkgname);

    if (pkg != NULL) {
        const char *p = function_package_get_string(pkg, "fname");

        strcpy(path, p);
    } else if (has_suffix(pkgname, ".gfn")) {
        err = get_full_filename(pkgname, path, OPT_I);
    } else {
        gchar *gfn = g_strdup_printf("%s.gfn", pkgname);

        err = get_full_filename(gfn, path, OPT_I);
        g_free(gfn);
    }

    if (opt & OPT_Q) {
        /* --quiet */
        gretl_bundle *b = gretl_bundle_new();

	if (b != NULL) {
	    if (!err) {
		bundle_function_package_info(path, b);
	    }
            set_last_result_data(b, GRETL_TYPE_BUNDLE);
        }
    } else {
        if (err) {
            pprintf(prn, _("%s: not found\n\n"), pkgname);
        } else {
            print_function_package_info(path, 0, prn);
        }
    }
}

static int lib_run_pkg_sample (const char *pkgname,
			       char *runfile,
			       DATASET *dset,
			       ExecState *s,
			       PRN *prn)
{
    char *buf = NULL;
    int err = 0;

    err = grab_package_sample(pkgname, &buf);

    if (!err) {
	gchar *base = g_strdup_printf("%s_sample.inp", pkgname);
	gchar *fname = gretl_make_dotpath(base);
	FILE *fp = gretl_fopen(fname, "w");

	if (fp == NULL) {
	    err = 1;
	} else {
	    fputs(buf, fp);
	    fclose(fp);
	    strcpy(runfile, fname);
	    err = run_script(runfile, s, dset, OPT_NONE, prn);
	    gretl_remove(runfile);
	}
	g_free(base);
	g_free(fname);
	free(buf);
    }

    return err;
}

static int do_pkg_command (char *readfile,
			   DATASET *dset,
                           ExecState *s,
                           PRN *prn)
{
    const char *action = s->cmd->param;
    const char *pkgname = s->cmd->parm2;
    gretlopt opt = s->cmd->opt;
    int err = 0;

    if (!strcmp(action, "install")) {
        err = install_package(pkgname, opt, s, prn);
    } else if (!strcmp(action, "unload")) {
        err = uninstall_function_package(pkgname, OPT_NONE, prn);
    } else if (!strcmp(action, "remove")) {
        err = uninstall_function_package(pkgname, OPT_P, prn);
    } else if (!strcmp(action, "query")) {
        query_package(pkgname, opt, prn);
    } else if (!strcmp(action, "run-sample")) {
	err = lib_run_pkg_sample(pkgname, readfile, dset, s, prn);
    } else if (!strcmp(action, "index") && !strcmp(pkgname, "addons")) {
        update_addons_index((opt & OPT_V)? prn : NULL);
    } else {
        gretl_errmsg_sprintf(_("pkg: unknown action '%s'"), action);
        err = E_PARSE;
    }

    return err;
}

static void print_info (gretlopt opt, DATASET *dset, PRN *prn)
{
    if (dset != NULL && dset->descrip != NULL) {
        pprintf(prn, "%s\n", dset->descrip);
    } else {
        pputs(prn, _("No data information is available.\n"));
    }
}

/* After estimating a model, check its errcode member to see
   if anything went wrong, and reset gretl_errno to zero.

   If we're looping (that is, if a loop is in progress at the
   current level of function execution) and @loop_force is 0,
   that's all, but if not then:

   (a) print the model (this may require special handling inside
   loops);

   (b) if the user has employed the "name <- command" mechanism,
   attach the supplied name to the model;

   (c) conditionally add the model to the stack in objstack.c,
   and if this is done, signal the fact by setting the 'pmod'
   member of @ExecState.

   (d) if we're called by the GUI program and the model has
   been assigned a name, activate the callback that adds the
   model to the GUI session.
*/

static int print_save_model (MODEL *pmod, DATASET *dset,
                             gretlopt opt, int loop_force,
                             PRN *prn, ExecState *s)
{
    int err = pmod->errcode;

    if (!err) {
        set_gretl_errno(0);
        if (!gretl_looping_currently() || loop_force) {
            int havename = *s->cmd->savename != '\0';
            int window = (opt & OPT_W) != 0;
            gretlopt popt;

            if (havename) {
                gretl_model_set_name(pmod, s->cmd->savename);
            }

            popt = get_printmodel_opt(pmod, opt);
            printmodel(pmod, dset, popt, prn);
            attach_subsample_to_model(pmod, dset);
            s->pmod = maybe_stack_model(pmod, s->cmd, prn, &err);
            if (!err && gretl_in_gui_mode() && s->callback != NULL &&
                (havename || window)) {
                if (s->pmod != NULL && (opt & OPT_Q) && !window) {
                    /* With OPT_Q (--quiet) and without the --window
                       flag, this model will not have a unique ID;
                       but that will be needed if it's going to be
                       a fully fledged gui model, as requested by
                       its having been given a "savename".
                    */
                    set_model_id(s->pmod, OPT_NONE);
                }
                s->callback(s, s->pmod, GRETL_OBJ_EQN);
            }
        }
    }

    return err;
}

static void save_var_vecm (ExecState *s)
{
    maybe_stack_var(s->var, s->cmd);

    if (gretl_in_gui_mode() && s->callback != NULL) {
        int havename = *s->cmd->savename != '\0';
        int window = (s->cmd->opt & OPT_W) != 0;

        if (havename || window) {
            s->callback(s, s->var, GRETL_OBJ_VAR);
        }
    }
}

static void gui_save_system (ExecState *s)
{
    equation_system *sys = s->sys;

    if (!gretl_in_gui_mode() || s->callback == NULL) {
        return;
    }
    /* note: with GRETL_OBJ_SYS, the business of calling
       "maybe_stack" is handled within system.c, so here
       all we have to do is invoke the GUI callback, if
       appropriate
    */
    if (sys == NULL && s->cmd->ci == ESTIMATE && s->cmd->param != NULL) {
        sys = get_equation_system_by_name(s->cmd->param);
    }
    if (sys != NULL && (*s->cmd->savename != '\0' || (s->cmd->opt & OPT_W))) {
        s->callback(s, sys, GRETL_OBJ_SYS);
    }
}

static int model_test_check (CMD *cmd, DATASET *dset, PRN *prn)
{
    int err = last_model_test_ok(cmd->ci, cmd->opt, dset, prn);

    if (err == E_DATA && cmd->ci == RESTRICT && has_param(cmd)) {
        /* try for a not-yet estimated anonymous system */
        if (get_anonymous_equation_system() != NULL) {
            gretl_error_clear();
            err = 0;
        }
    }

    return err;
}

static int get_line_continuation (char *line, FILE *fp, PRN *prn)
{
    char tmp[MAXLINE];
    int err = 0;

    if (!strncmp(line, "quit", 4)) {
        return 0;
    }

    while (top_n_tail(line, MAXLINE, &err)) {
        if (err) {
            break;
        }
        *tmp = '\0';
        if (fgets(tmp, sizeof tmp, fp) && *tmp != '\0') {
            if (strlen(line) + strlen(tmp) > MAXLINE - 1) {
                pprintf(prn, _("Maximum length of command line "
                               "(%d bytes) exceeded"), MAXLINE);
                pputc(prn, '\n');
                err = E_TOOLONG;
                break;
            } else {
                strcat(line, tmp);
                compress_spaces(line);
            }
        }
    }

    return err;
}

static int run_script (const char *fname, ExecState *s,
                       DATASET *dset, gretlopt opt,
                       PRN *prn)
{
    int indent = gretl_if_state_record();
    int echo = gretl_echo_on();
    int messages = gretl_messages_on();
    FILE *fp;
    int iferr, err = 0;

    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
        gretl_errmsg_sprintf(_("Couldn't open %s"), fname);
        return E_FOPEN;
    }

    strcpy(s->runfile, fname);

    if (opt & OPT_Q) {
        set_gretl_echo(0);
        set_gretl_messages(0);
    }

    if (gretl_echo_on()) {
        pprintf(prn, "run \"%s\"\n", fname);
    }

    while (fgets(s->line, MAXLINE - 1, fp) && !err) {
        err = get_line_continuation(s->line, fp, prn);
        if (!err) {
            err = maybe_exec_line(s, dset, NULL);
        }
    }

    fclose(fp);

    if (opt & OPT_Q) {
        set_gretl_echo(echo);
        set_gretl_messages(messages);
    }

    iferr = gretl_if_state_check(indent);
    if (iferr && !err) {
        err = iferr;
    }

    return err;
}

static int lib_clear_data (ExecState *s, DATASET *dset)
{
    int err = 0;

    if (dset->Z != NULL) {
        err = restore_full_sample(dset, NULL);
        free_Z(dset);
    }

    clear_model(s->model);
    clear_datainfo(dset, CLEAR_FULL);
    libgretl_session_cleanup(SESSION_CLEAR_DATASET);
    set_model_count(0);
    gretl_cmd_destroy_context(s->cmd);

    return err;
}

static int join_aggregation_method (const char *s, int *seqval,
                                    char **auxname, int *err)
{
    int ret = -1;

    if (!strncmp(s, "seq:", 4)) {
        char *endptr;

        *seqval = (int) strtol(s + 4, &endptr, 10);
        if (*endptr == '\0' && *seqval != 0) {
            ret = AGGR_SEQ;
        } else {
            gretl_errmsg_sprintf(_("%s: invalid input '%s'\n"), "--seq", s + 4);
            *err = E_DATA;
        }
    } else if (!strcmp(s, "count")) {
        ret = AGGR_COUNT;
    } else if (!strcmp(s, "avg")) {
        ret = AGGR_AVG;
    } else if (!strcmp(s, "sum")) {
        ret = AGGR_SUM;
    } else if (!strcmp(s, "min")) {
        ret = AGGR_MIN;
    } else if (!strcmp(s, "max")) {
        ret = AGGR_MAX;
    } else if (!strcmp(s, "none")) {
        ret = AGGR_NONE;
    } else if (!strcmp(s, "spread")) {
        ret = AGGR_MIDAS;
    } else if (!strncmp(s, "min(", 4) ||
               !strncmp(s, "max(", 4)) {
        const char *p = strchr(s + 4, ')');

        if (p != NULL && strlen(p) == 1) {
            int len = p - (s + 4);

            if (len > 0) {
                *auxname = gretl_strndup(s + 4, len);
                if (*auxname == NULL) {
                    *err = E_ALLOC;
                } else {
                    ret = (s[1] == 'a')? AGGR_MAX : AGGR_MIN;
                }
            }
        } else {
            *err = E_PARSE;
        }
    } else {
        *err = E_PARSE;
    }

    return ret;
}

static int get_inner_key_id (const char *s, int n,
                             const DATASET *dset,
                             int *err)
{
    char vname[VNAMELEN];
    int id = -1;

    if (n == 0 || n >= VNAMELEN) {
        *err = E_PARSE;
    } else {
        *vname = '\0';
        strncat(vname, s, n);
        if (gretl_namechar_spn(vname) != n) {
            gretl_errmsg_sprintf(_("field '%s' in command is invalid"), vname);
            *err = E_PARSE;
        } else {
            id = current_series_index(dset, vname);
            if (id < 0) {
                gretl_errmsg_sprintf(_("'%s': no such series"), vname);
                *err = E_UNKVAR;
            }
        }
    }

    return id;
}

static int *get_inner_keys (const char *s, DATASET *dset,
                            int *err)
{
    int *klist = NULL;
    int ikey1 = -1, ikey2 = -1;
    int nkeys = 0;

    if (strchr(s, ',') == NULL) {
        /* just one key, fine */
        ikey1 = current_series_index(dset, s);
        if (ikey1 < 0) {
            gretl_errmsg_sprintf(_("'%s': no such series"), s);
            *err = E_UNKVAR;
        } else {
            nkeys = 1;
        }
    } else {
        /* we should have a double key */
        int n = strcspn(s, ",");

        ikey1 = get_inner_key_id(s, n, dset, err);

        if (!*err) {
            s += n + 1;
            n = strlen(s);
            ikey2 = get_inner_key_id(s, n, dset, err);
        }

        if (!*err) {
            nkeys = 2;
        }
    }

    if (!*err) {
        klist = gretl_list_new(nkeys);
        if (klist == NULL) {
            *err = E_ALLOC;
        } else {
            klist[1] = ikey1;
            if (nkeys == 2) {
                klist[2] = ikey2;
            }
        }
    }

    return klist;
}

static int check_import_names (char **S, int ns,
			       int ci, DATASET *dset)
{
    int i, err = 0;

    for (i=0; i<ns && !err; i++) {
        if (S[i] == NULL || S[i][0] == '\0') {
            err = E_DATA;
        } else if (ci == JOIN && (strchr(S[i], '*') || strchr(S[i], '?'))) {
            ; /* wildcards: may be OK? */
        } else if (ci == JOIN && current_series_index(dset, S[i]) < 0) {
            err = check_varname(S[i]);
            if (!err && gretl_type_from_name(S[i], NULL)) {
                err = E_TYPES;
            }
        }
    }

    return err;
}

static char **names_from_array_arg (gretl_array *A,
                                    int *ns,
                                    int *err)
{
    char **S = NULL;

    if (gretl_array_get_type(A) != GRETL_TYPE_STRINGS) {
        *err = E_TYPES;
    } else {
        S = gretl_array_get_strings(A, ns);
        if (S == NULL) {
            *err = E_DATA;
        }
    }

    return S;
}

static char **strings_array_singleton (const char *s,
                                       int *err)
{
    int len = strlen(s) + 1;
    char **S = strings_array_new_with_length(1, len);

    if (S == NULL) {
        *err = E_ALLOC;
    } else {
        strcat(S[0], s);
    }

    return S;
}

static int get_selected_import_names (const char *s,
				      int ci,
				      DATASET *dset,
				      char ***pvnames,
				      int *pnvars)
{
    char **S = NULL;
    gretl_array *A = NULL;
    int ns = 0;
    int err = 0;

    if (s == NULL) {
        return E_PARSE;
    }

    if (strchr(s, ' ') != NULL) {
        /* @s should hold two or more names */
        S = gretl_string_split(s, &ns, NULL);
        if (S == NULL) {
            err = E_DATA;
        }
    } else if ((A = get_array_by_name(s)) != NULL) {
        /* @s should be the name of an array of strings */
        S = names_from_array_arg(A, &ns, &err);
    } else {
        /* @s should be a legit series name */
        S = strings_array_singleton(s, &err);
        ns = 1;
    }

    if (S != NULL) {
        err = check_import_names(S, ns, ci, dset);
    }

    if (!err) {
        if (A != NULL) {
            /* copy strings "borrowed" from array */
            *pvnames = strings_array_dup(S, ns);
            if (*pvnames == NULL) {
                err = E_ALLOC;
            }
        } else {
            /* grab strings allocated here */
            *pvnames = S;
        }
        *pnvars = ns;
    }

    return err;
}

int lib_join_data (const char *param,
                   const char *filename,
                   DATASET *dset,
                   gretlopt opt,
                   PRN *prn)
{
    gretlopt opts[] = {
        OPT_I, /* ikey: inner key(s) */
        OPT_O, /* okey: outer key(s) */
        OPT_F, /* filter: filter expression */
        OPT_A, /* aggr: aggregation method */
        OPT_D, /* data: "payload" spec */
        OPT_K, /* tkey: outer time-key name,format */
        OPT_X, /* tconvert: date columns for conversion */
        OPT_T, /* tconv-fmt: format for "tconvert" */
        OPT_P, /* pd: outer data frequency */
        0
    };
    char *okey = NULL, *filter = NULL;
    char **vnames = NULL;
    char *dataname = NULL;
    char *auxname = NULL;
    char *tconvstr = NULL;
    char *tconvfmt = NULL;
    int *ikeyvars = NULL;
    int aggr = 0, seqval = 0;
    int tseries, nvars = 1;
    int midas_pd = 0;
    int i, err = 0;

    tseries = dataset_is_time_series(dset);

    if (opt & OPT_K) {
        /* --tkey implies special handling of keys */
        if (opt & (OPT_I | OPT_O)) {
            return E_BADOPT;
        } else if (!tseries) {
            return E_PDWRONG;
        }
    }

    err = get_selected_import_names(param, JOIN, dset, &vnames, &nvars);

    if (!err && nvars > 1) {
        /* multiple series: we can't handle the --data option */
        if (opt & OPT_D) {
            err = E_BADOPT;
        }
    }

    for (i=0; opts[i] && !err; i++) {
        gretlopt jopt = opts[i];
        const char *optparm;

        if (opt & jopt) {
            optparm = get_optval_string(JOIN, jopt);
            if (optparm == NULL) {
                gretl_errmsg_set(_("Missing option parameter"));
                err = E_DATA;
            } else if (jopt == OPT_I) {
                /* --ikey: the inner key(s) string */
                ikeyvars = get_inner_keys(optparm, dset, &err);
            } else if (jopt == OPT_O) {
                /* --okey: the outer key(s) string */
                okey = gretl_strdup(optparm);
            } else if (jopt == OPT_F) {
                /* --filter: string specifying a row filter */
                filter = gretl_strdup(optparm);
            } else if (jopt == OPT_A) {
                /* --aggr: aggregation */
                aggr = join_aggregation_method(optparm, &seqval,
                                               &auxname, &err);
            } else if (jopt == OPT_D) {
                /* --data: string specifying the outer data series */
                dataname = gretl_strdup(optparm);
            } else if (jopt == OPT_K) {
                /* --tkey: string specifying outer time key */
                okey = gretl_strdup(optparm);
            } else if (jopt == OPT_X) {
                /* --tconvert: list of time/date cols */
                tconvstr = gretl_strdup(optparm);
            } else if (jopt == OPT_T) {
                /* --tconv-fmt: format for tconvert columns */
                tconvfmt = gretl_strdup(optparm);
            } else if (jopt == OPT_P) {
                midas_pd = atoi(optparm);
            }
        }
    }

    if (!err && okey != NULL && ikeyvars == NULL && !(opt & OPT_K)) {
        /* We can't have an outer key but no inner one, unless
           we're matching by the time-series structure of the
           left-hand dataset (implied by OPT_K)
        */
        gretl_errmsg_set(_("Inner key is missing"));
        err = E_PARSE;
    }

    if (!err && aggr != 0 && ikeyvars == NULL && !tseries) {
        /* aggregation requires ikeyvars, unless there's
           an implicit time-series inner key
        */
        gretl_errmsg_set(_("Inner key is missing"));
        err = E_ARGS;
    }

    if (!err) {
        err = gretl_join_data(filename,
                              (const char **) vnames,
                              nvars, dset,
                              ikeyvars, okey, filter,
                              dataname, aggr, seqval,
                              auxname, tconvstr,
                              tconvfmt, midas_pd,
                              opt, prn);
    }

    strings_array_free(vnames, nvars);
    free(ikeyvars);
    free(okey);
    free(filter);
    free(dataname);
    free(auxname);
    free(tconvstr);
    free(tconvfmt);

    return err;
}

static int is_http (const char *s)
{
    return (strncmp(s, "http://", 7) == 0 ||
            strncmp(s, "https://", 8) == 0);
}

static void open_op_init (CMD *cmd, OpenOp *op)
{
    op->fname[0] = '\0';
    op->quiet = (cmd->opt & OPT_Q)? 1 : 0;
    op->http = 0;
    op->dbdata = 0;
    op->ftype = -1;
}

/* selection of series by name: applicable only for "open"
   for native gretl data files
*/

static int check_import_subsetting (CMD *cmd, OpenOp *op)
{
    if (cmd->ci != OPEN || (op->ftype != GRETL_XML_DATA &&
			    op->ftype != GRETL_BINARY_DATA)) {
	return E_BADOPT;
    } else {
	const char *s = get_optval_string(OPEN, OPT_E);

	if (get_array_by_name(s)) {
	    /* protect array from deletion */
	    cmd->opt |= OPT_P;
	}
	return 0;
    }
}

static int open_append_stage_1 (CMD *cmd,
                                DATASET *dset,
                                OpenOp *op,
                                PRN *prn)
{
    gretlopt opt = cmd->opt;
    int pkgdata = 0;
    int err = 0;

    open_op_init(cmd, op);

    /* initial detection of some open/append variants */

    if (opt & OPT_W) {
        /* --www: remote databases only, no other options
           applicable
        */
        if (opt != OPT_W) {
            return E_BADOPT;
        }
        op->ftype = GRETL_NATIVE_DB_WWW;
        op->dbdata = 1;
        strncat(op->fname, cmd->param, MAXLEN - 1);
    } else if (cmd->ci != JOIN && (opt & OPT_O)) {
        op->ftype = GRETL_ODBC;
        op->dbdata = 1;
    } else if ((cmd->ci == OPEN && (opt & OPT_K)) ||
	       (cmd->ci == APPEND && (opt & OPT_K)) ||
	       (cmd->ci == JOIN && (opt & OPT_R))) {
        /* --frompkg=whatever */
        pkgdata = 1;
    } else if (!strcmp(cmd->param, "dbnomics")) {
	strcpy(op->fname, "dbnomics");
        op->ftype = GRETL_DBNOMICS;
        op->dbdata = 1;
    } else if (is_http(cmd->param)) {
        op->http = 1;
    }

    if (!op->dbdata) {
        if (op->http) {
#ifdef USE_CURL
            err = try_http(cmd->param, op->fname, NULL);
#else
            gretl_errmsg_set("http resource: cURL is not available");
            err = E_DATA;
#endif
        } else if (pkgdata) {
            err = get_package_data_path(cmd->ci, cmd->param, op->fname);
        } else {
            err = get_full_filename(cmd->param, op->fname, OPT_NONE);
        }
    }

    if (!err) {
	if (op->ftype < 0) {
	    op->ftype = detect_filetype(op->fname, OPT_P);
	}
	if (opt & OPT_E) {
	    err = check_import_subsetting(cmd, op);
	}
    }

    if (op->ftype == GRETL_SESSION) {
	if (gretl_in_gui_mode() && (cmd->opt & OPT_U)) {
	    ; /* experimental, could be OK */
	} else {
	    gretl_errmsg_set(_("gretl session files can only be opened via the GUI program"));
	    err = E_DATA;
	}
    }

    if (err) {
        errmsg(err, prn);
        return err;
    }

    if (cmd->ci == JOIN) {
        if (op->ftype == GRETL_CSV || op->ftype == GRETL_XML_DATA ||
            op->ftype == GRETL_BINARY_DATA) {
            set_dataset_is_changed(dset, 0);
            err = lib_join_data(cmd->parm2, op->fname, dset, opt, prn);
        } else {
            gretl_errmsg_set(_("join: only CSV and gdt[b] files are supported"));
            err = E_DATA;
        }
        if (err) {
            errmsg(err, prn);
        }
    } else {
        if (!op->dbdata) {
            op->dbdata = (op->ftype == GRETL_NATIVE_DB ||
                          op->ftype == GRETL_RATS_DB ||
                          op->ftype == GRETL_PCGIVE_DB);
        }
        if (cmd->ci == OPEN && !op->dbdata && gretl_function_depth() > 0) {
            gretl_errmsg_sprintf(_("The \"%s\" command cannot be used in this context"),
                                 gretl_command_word(cmd->ci));
            err = E_DATA;
        }
    }

    return err;
}

/* respond to --select (select specific series) on OPEN
   for native gdt or gdtb data files
*/

static int handle_gdt_selection (const char *fname,
				 DATASET *dset,
				 gretlopt opt,
				 PRN *prn)
{
    const char *s = get_optval_string(OPEN, OPT_E);
    char **S_sel = NULL;
    char **S_ok = NULL;
    int n_ok, n_sel = 0;
    int err = 0;

    if (s == NULL || *s == '\0') {
	return E_BADOPT;
    }

    err = get_selected_import_names(s, OPEN, dset, &S_sel, &n_sel);
    if (!err) {
	err = gretl_read_gdt_varnames(fname, &S_ok, &n_ok);
    }

    if (!err) {
	int *list = gretl_list_new(n_sel);
	int i, j, k = 0;

	for (j=0; j<n_sel; j++) {
	    for (i=1; i<n_ok; i++) {
		if (!strcmp(S_sel[j], S_ok[i])) {
		    list[++k] = i;
		    break;
		}
	    }
	}
	if (k != n_sel) {
	    pputs(prn, _("Invalid selection"));
	    pputc(prn, '\n');
	    err = E_DATA;
	} else {
	    err = gretl_read_gdt_subset(fname, dset, list, opt);
	}
	free(list);
    }

    strings_array_free(S_ok, n_ok);
    strings_array_free(S_sel, n_sel);

    return err;
}

static int lib_open_append (ExecState *s,
                            DATASET *dset,
                            char *newfile,
                            PRN *prn)
{
    CMD *cmd = s->cmd;
    gretlopt opt = cmd->opt;
    int catch = cmd->flags & CMD_CATCH;
    OpenOp op = {0};
    PRN *vprn = prn;
    int err;

    *newfile = '\0';

    err = open_append_stage_1(cmd, dset, &op, prn);

    if (cmd->ci == JOIN) {
        /* handled by stage 1 */
        return err;
    } else if (err) {
        /* error at stage 1, don't proceed */
        return err;
    }

    if (cmd->ci == OPEN && op.ftype == GRETL_SESSION) {
	/* open session as bundle in GUI */
	s->callback(s, NULL, GRETL_OBJ_SESSION);
	return 0;
    }

    if (cmd->ci == OPEN && !op.dbdata) {
        if (s->callback != NULL) {
            /* allow aborting a destructive "open" */
            int stop = s->callback(s, dset, GRETL_OBJ_DSET);

            if (stop) return 0;
        } else {
            lib_clear_data(s, dset);
        }
    }

    if (op.quiet) {
	if (catch) {
	    vprn = NULL;
	} else {
	    /* in case we hit any problems below... */
	    vprn = gretl_print_new(GRETL_PRINT_BUFFER, NULL);
	}
    }

    if (op.ftype == GRETL_XML_DATA || op.ftype == GRETL_BINARY_DATA) {
	if (opt & OPT_E) {
	    err = handle_gdt_selection(op.fname, dset, opt, vprn);
	} else {
	    err = gretl_read_gdt(op.fname, dset, opt, vprn);
	}
    } else if (op.ftype == GRETL_CSV) {
        err = import_csv(op.fname, dset, opt, vprn);
    } else if (SPREADSHEET_IMPORT(op.ftype)) {
        err = import_spreadsheet(op.fname, op.ftype, cmd->list, cmd->parm2,
                                 dset, opt, vprn);
    } else if (OTHER_IMPORT(op.ftype)) {
        err = import_other(op.fname, op.ftype, dset, opt, vprn);
    } else if (op.ftype == GRETL_ODBC) {
        err = set_odbc_dsn(cmd->param, vprn);
    } else if (op.dbdata) {
        err = set_db_name(op.fname, op.ftype, vprn);
    } else {
        err = gretl_seek_data(op.fname, dset, opt, vprn);
    }

    if (vprn != NULL && vprn != prn) {
        if (err) {
            /* The user asked for --quiet operation, but something
               went wrong so let's print any info we got on
               @vprn.
            */
            const char *buf = gretl_print_get_buffer(vprn);

            if (buf != NULL && *buf != '\0') {
                pputs(prn, buf);
            }
        } else if (gretl_messages_on() &&
                   op.ftype != GRETL_NATIVE_DB_WWW &&
                   op.ftype != GRETL_ODBC) {
#ifdef HAVE_MPI
	    /* print minimal success message? */
	    if (!gretl_mpi_initialized()) {
		pprintf(prn, _("Read datafile %s\n"), op.fname);
	    }
#else
	    pprintf(prn, _("Read datafile %s\n"), op.fname);
#endif
        }
        gretl_print_destroy(vprn);
    }

    if (err) {
	if (op.quiet && catch) {
	    ; /* keep quiet */
	} else {
	    errmsg(err, prn);
	}
    } else {
        if (dset->v > 0 && !op.dbdata && !op.quiet) {
            list_series(dset, OPT_NONE, prn);
        }
        if (op.http) {
            remove(op.fname);
        }
        if (!op.dbdata && !op.http) {
            /* redundant? */
            strcpy(newfile, op.fname);
        }
        if (cmd->ci == OPEN && !op.dbdata) {
            if (s->callback != NULL) {
                if (op.http) {
                    op.fname[0] = '\0';
                }
                s->callback(s, &op, GRETL_OBJ_ANY);
            }
        }
    }

    return err;
}

static int check_clear (gretlopt opt)
{
    int err = 0;

    if (gretl_function_depth() > 0) {
        gretl_errmsg_sprintf(_("The \"%s\" command cannot be used in this context"),
                             gretl_command_word(CLEAR));
        err = E_DATA;
    } else {
	err = incompatible_options(opt, OPT_A | OPT_D | OPT_F);
    }

    return err;
}

static EXEC_CALLBACK gui_callback;

static void schedule_callback (ExecState *s)
{
    if (s->callback != NULL) {
        s->flags |= CALLBACK_EXEC;
    }
}

static void maybe_schedule_set_callback (ExecState *s)
{
    if (s->callback != NULL && s->cmd->param != NULL) {
	if (!strcmp(s->cmd->param, "plot_collection")) {
	    s->flags |= CALLBACK_EXEC;
	}
    }
}

static int callback_scheduled (ExecState *s)
{
    return (s->flags & CALLBACK_EXEC) ? 1 : 0;
}

static void callback_exec (ExecState *s, char *fname, int err)
{
    if (s->cmd->ci == MODELTAB) {
	if (gui_callback != NULL) {
	    gui_callback(s, NULL, 0);
	}
    } else if (!err && s->callback != NULL) {
        if (s->cmd->ci == OPEN) {
            s->callback(s, fname, 0);
        } else {
            s->callback(s, NULL, 0);
        }
    }

    s->flags &= ~CALLBACK_EXEC;
    *s->cmd->savename = '\0';
}

/* for use in contexts where we don't have an
   ExecState (or CMD) to hand
*/

void manufacture_gui_callback (int ci)
{
    if (gui_callback != NULL) {
        ExecState s = {0};
        CMD cmd = {0};

        cmd.ci = ci;
        if (ci == FLUSH) {
            cmd.opt = OPT_Q;
        }
        s.cmd = &cmd;
        gui_callback(&s, NULL, 0);
    }
}

/* whether we're in gui or cli mode, try to ensure
   that printed output gets displayed
*/

void gretl_flush (PRN *prn)
{
    if (gui_callback != NULL) {
        ExecState s = {0};
        CMD cmd = {0};

        cmd.ci = FLUSH;
        cmd.opt = OPT_Q;
        s.cmd = &cmd;
        gui_callback(&s, NULL, 0);
    } else {
        gretl_print_flush_stream(prn);
    }
}

static int do_end_restrict (ExecState *s, DATASET *dset)
{
    GretlObjType otype = gretl_restriction_get_type(s->rset);
    gretlopt ropt = gretl_restriction_get_options(s->rset);
    gretlopt opt = s->cmd->opt | ropt;
    int err = 0;

    if (opt & OPT_F) {
        /* restrict --full */
        if (otype == GRETL_OBJ_VAR) {
            s->var = gretl_restricted_vecm(s->rset, dset,
                                           opt, s->prn, &err);
            if (!err && s->var != NULL) {
                save_var_vecm(s);
            }
        } else if (otype == GRETL_OBJ_EQN) {
            err = gretl_restriction_finalize_full(s, s->rset, dset,
                                                  opt, s->prn);
            if (!err) {
                gretlopt printopt = OPT_NONE;

                if (opt & (OPT_Q | OPT_S)) {
                    printopt = OPT_Q;
                }
                print_save_model(s->pmod, dset, printopt, 1,
                                 s->prn, s);
            }
        }
    } else {
        err = gretl_restriction_finalize(s->rset, dset,
                                         opt, s->prn);
    }

    s->rset = NULL;

    return err;
}

static void exec_state_prep (ExecState *s)
{
    s->flags &= ~CALLBACK_EXEC;
    s->pmod = NULL;
}

int gretl_delete_variables (int *list,
                            const char *param,
                            gretlopt opt,
                            DATASET *dset,
                            int *renumber,
                            PRN *prn)
{
    int err;

    err = incompatible_options(opt, OPT_T | OPT_D | OPT_F);

    if (!err) {
        if (opt & OPT_T) {
            /* delete all vars of given type */
            if (list != NULL || param != NULL) {
                err = E_BADOPT;
            }
        } else if (opt & OPT_D) {
            /* delete named vars from database */
            if (list != NULL || param == NULL) {
                err = E_BADOPT;
            }
        }
    }

    if (err) {
        return err;
    }

    if (opt & OPT_T) {
        const char *s = get_optval_string(DELEET, OPT_T);

        if (s == NULL) {
            err = E_ARGS;
        } else {
            GretlType type = gretl_type_from_string(s);

            err = delete_user_vars_of_type(type, prn);
        }
    } else if (opt & OPT_D) {
        err = db_delete_series_by_name(param, prn);
    } else if (param != NULL) {
        err = gretl_delete_var_by_name(param, prn);
    } else if (list != NULL) {
        /* delete listed series from dataset */
        if (renumber == NULL && !(opt & OPT_F)) {
            /* lacking the --force option */
            pputs(prn, _("You cannot delete series in this context\n"));
            err = E_DATA;
        } else {
            err = dataset_drop_listed_variables(list, dset,
                                                renumber, prn);
        }
    } else {
        err = E_DATA;
    }

    return err;
}

/* OMIT and ADD: if we're estimating a revised model, should
   we be saving it as the "last model", or are we just treating
   the command as a stand-alone test?
*/

static int add_omit_save (CMD *cmd, MODEL *pmod)
{
    int ret = 1;

    if (cmd->ci == ADD) {
        if (cmd->opt & OPT_L) {
            /* not saving if given the --lm option */
            ret = 0;
        }
    } else if (cmd->ci == OMIT) {
        if (cmd->opt & OPT_W) {
            /* not saving under the --test-only (Wald) option */
            ret = 0;
        } else if (!(cmd->opt & OPT_A)) {
            /* not in --auto mode */
            if (cmd->list != NULL && pmod->list != NULL &&
                cmd->list[0] == pmod->list[0] - 1) {
                /* omitting everything, can't save */
                ret = 0;
            }
        }
    }

    return ret;
}

static int VAR_omit_driver (CMD *cmd, DATASET *dset, PRN *prn)
{
    GRETL_VAR *var = get_last_model(NULL);
    int err = 0;

    if (cmd->opt & OPT_W) {
        /* Wald test using VCV */
        err = gretl_VAR_wald_omit_test(var, cmd->list, dset,
                                       cmd->opt, prn);
    } else {
        /* the full deal: estimate reduced system */
        GRETL_VAR *vnew;

        vnew = gretl_VAR_omit_test(var, cmd->list, dset, cmd->opt,
                                   prn, &err);
        if (!err) {
            err = maybe_stack_var(vnew, cmd);
        }
    }

    return err;
}

static int model_print_driver (MODEL *pmod, DATASET *dset,
                               int ci, const char *param,
                               gretlopt opt, PRN *prn)
{
    int err = incompatible_options(opt, OPT_R | OPT_C);

    if (!err) {
        char fname[FILENAME_MAX];

        *fname = '\0';

        if (param != NULL) {
            /* the legacy mechanism */
            strcpy(fname, param);
        } else if (opt & OPT_U) {
            /* try for --output=filename, and if found let
               the suffix determine the output type
            */
            const char *s = get_optval_string(ci, OPT_U);

            if (s != NULL && *s != '\0') {
                strcpy(fname, s);
                if (has_suffix(fname, ".rtf")) {
                    opt |= OPT_R;
                } else if (has_suffix(fname, ".csv")) {
                    opt |= OPT_C;
                }
            }
        }

        if (*fname == '\0') {
            /* fallback */
            const char *sfx = (opt & OPT_R)? "rtf" :
                (opt & OPT_C)? "csv" : "tex";

            if (pmod->ID > 0) {
                sprintf(fname, "model_%d.%s", pmod->ID, sfx);
            } else {
                /* FIXME: this needs to be checked! */
                sprintf(fname, "model_%" G_GINT64_FORMAT ".%s",
                        pmod->esttime, sfx);
            }
        }

        if (opt & OPT_R) {
            err = rtfprint(pmod, dset, fname, opt);
        } else if (opt & OPT_C) {
            err = csvprint(pmod, dset, fname, opt);
        } else {
            gretlopt texopt = opt;

            if (ci == EQNPRINT) {
                texopt |= OPT_E;
            }
            err = texprint(pmod, dset, fname, texopt);
        }
        if (!err) {
            pprintf(prn, _("Model printed to %s\n"), fname);
        }
    }

    return err;
}

#if USE_CURL

static int package_check_dependencies (const char *fname,
                                       ExecState *s,
                                       PRN *prn)
{
    char **depends;
    int ndeps;
    int err = 0;

    if (has_suffix(fname, ".zip")) {
        /* we need to map from @fname to the gfn name */
        gchar *tmp, *gfnname;
        const char *p;

        tmp = g_strndup(fname, strlen(fname) - 4);
        p = path_last_element(tmp);
        gfnname = g_strdup_printf("%s%c%s.gfn", tmp, SLASH, p);
        depends = package_peek_dependencies(gfnname, &ndeps);
        g_free(tmp);
        g_free(gfnname);
    } else {
        depends = package_peek_dependencies(fname, &ndeps);
    }

    if (depends != NULL) {
        char *pkgpath;
        int i;

        for (i=0; i<ndeps && !err; i++) {
            pkgpath = gretl_function_package_get_path(depends[i], PKG_ALL);
            if (pkgpath == NULL) {
                err = install_package(depends[i], OPT_D, s, prn);
		fprintf(stderr, "pkg install dependency %s: err = %d\n",
			depends[i], err);
            }
            free(pkgpath);
        }
        strings_array_free(depends, ndeps);
    }

    return err;
}

/* For now DO_BINPGK (that is, handle binary packages -- x13as and
   tramo-seats -- via the "pkg" command) is specific to macOS. It
   could in principle be extended to other OSes.
*/

#ifdef __APPLE__
# define DO_BINPKG 1
#else
# define DO_BINPKG 0
#endif

#if DO_BINPKG

static int handle_binary_package (const char *fname,
				  ExecState *es,
				  PRN *prn)
{
    const char *sf = "https://sourceforge.net/projects/gretl/files";
    const char *sfdir, *path_id = NULL;
    gchar *uri, *fullname;
    int err;

    fullname = gretl_make_dotpath(fname);
    if (strstr(fname, "tramo") != NULL) {
	sfdir = "tramo";
    } else {
	sfdir = "x13as";
    }
    uri = g_strdup_printf("%s/%s/%s", sf, sfdir, fname);
    err = retrieve_public_file(uri, fullname);

    if (!err) {
        gchar *s, *topdir = NULL;
        int try = 0;

    try_again:

        if (try == 0) {
            topdir = g_strdup(gretl_bindir());
            s = strstr(topdir, "/bin/");
        } else {
            /* in case "/bin/" didn't work */
            topdir = g_strdup(gretl_package_install_path("functions"));
            s = strstr(topdir, "/functions");
        }

        if (s == NULL) {
            err = E_DATA;
        } else {
            *s = '\0';
            fprintf(stderr, "topdir %d: '%s'\n", try, topdir);
            err = gretl_chdir(topdir);
            if (err) {
                fprintf(stderr, "chdir %d error\n", try);
            } else {
                err = gretl_untar(fullname);
                fprintf(stderr, "untar %d: err = %d\n", try, err);
            }
        }

        if (err && try == 0) {
            g_free(topdir);
            err = 0;
            try = 1;
            goto try_again;
        }
        if (!err) {
            if (strstr(fname, "tramo") != NULL) {
                path_id = "tramo";
                s = g_strdup_printf("%s/bin/tramo", topdir);
                gretl_set_path_by_name(path_id, s);
                g_free(s);
            } else if (strstr(fname, "x13as") != NULL) {
                path_id = "x12a";
                s = g_strdup_printf("%s/bin/x13as", topdir);
                gretl_set_path_by_name(path_id, s);
                g_free(s);
            }
        }
        gretl_remove(fullname);
        g_free(topdir);
    }

    if (!err && gretl_messages_on()) {
        pprintf(prn, _("Installed %s\n"), fname);
    }

    if (!err && es != NULL && path_id != NULL && gui_callback != NULL) {
        gretl_bundle *b = gretl_bundle_new();

        gretl_bundle_set_string(b, "path_id", path_id);
        gretl_bundle_set_int(b, "binpkg", 1);
        gui_callback(es, b, GRETL_OBJ_BUNDLE);
    }

    g_free(fullname);
    g_free(uri);

    return err;
}

#endif /* DO_BINPKG */

static int really_local (const char *pkgname)
{
#ifdef WIN32
    if (isalpha(pkgname[0]) && pkgname[1] == ':') {
        return 1;
    } else if (strchr(pkgname, '\\')) {
        return 1;
    }
#else
    if (strchr(pkgname, ':') == NULL && strchr(pkgname, '/')) {
        return 1;
    }
#endif
    return 0;
}

static void pkg_install_invoke_callback (ExecState *s,
					 const char *basename,
					 const char *fullname,
					 int filetype)
{
    gretl_bundle *b = gretl_bundle_new();

    if (strstr(basename, ".tar.gz")) {
	/* a data-file collection */
	gretl_bundle_set_string(b, "pkgname", basename);
    } else if (strstr(basename, "scripts")) {
        /* a script-file collection */
        gretl_bundle_set_string(b, "pkgname", basename);
        gretl_bundle_set_int(b, "zipfile", filetype == 2);
        gretl_bundle_set_int(b, "scripts", 1);
    } else {
	/* a function package */
	char *nosfx = g_strdup(basename);
	char *p = strrchr(nosfx, '.');

	if (p != NULL) {
	    *p = '\0';
	}
	gretl_bundle_set_string(b, "pkgname", nosfx);
	gretl_bundle_set_int(b, "zipfile", filetype == 2);
	g_free(nosfx);
    }

    gretl_bundle_set_string(b, "filename", fullname);
    gui_callback(s, b, GRETL_OBJ_BUNDLE);
}

static int install_package (const char *pkgname,
			    gretlopt opt,
			    ExecState *s,
			    PRN *prn)
{
    char *fname = NULL;
    int filetype = 0;
    int local = 0;
    int staging = 0;
    int scripts = 0;
    int addons = 0;
    int http = 0;
    int err;

    err = incompatible_options(opt, OPT_L | OPT_S);
    if (err) {
	return err;
    } else if (opt & OPT_L) {
	local = 1;
    } else if (opt & OPT_S) {
	staging = 1;
    }

#if DO_BINPKG
    if (!local && has_suffix(pkgname, ".tar.gz") &&
	(strstr(pkgname, "tramo") || strstr(pkgname, "x13as"))) {
	/* specific to macOS at present */
        return handle_binary_package(pkgname, s, prn);
    }
#endif

    if (!local && is_gretl_addon(pkgname)) {
	gretl_errmsg_set("Downloading individual addons is not supported");
	return E_DATA;
    }

    if (!strncmp(pkgname, "http://", 7) ||
        !strncmp(pkgname, "https://", 8)) {
        http = 1;
    }

    if (!local && !http && really_local(pkgname)) {
        local = 1;
    }

    if (has_suffix(pkgname, ".gfn")) {
        filetype = 1;
    } else if (has_suffix(pkgname, ".zip")) {
        filetype = 2;
        if (strstr(pkgname, "scripts")) {
            scripts = 1;
        }
    } else if (has_suffix(pkgname, ".tar.gz")) {
	filetype = 3;
	if (local) {
	    /* maybe FIXME ? */
	    gretl_errmsg_set("Installation of local tar.gz files is not supported");
	    return E_BADOPT;
	}
	if (strstr(pkgname, "addons")) {
	    addons = 1;
	}
    } else if (local || http) {
        /* must have suitable suffix */
        err = E_DATA;
    } else {
        /* from gretl server: determine the correct suffix */
        fname = retrieve_remote_pkg_filename(pkgname, &err);
        if (!err) {
            filetype = strstr(fname, ".zip") ? 2 : 1;
        }
    }

    if (!err) {
	if (filetype == 3) {
	    /* a tar.gz file */
	    gchar *dlpath = get_download_path(pkgname, &err);

	    if (!err && addons) {
		err = retrieve_addons_package(dlpath);
	    } else if (!err) {
		err = retrieve_remote_files_package(pkgname, dlpath);
	    }
	    if (!err) {
		err = unpack_files_collection(dlpath);
	    }
	    if (!err) {
		pprintf(prn, _("Installed %s\n"), pkgname);
	    }
	    gretl_remove(dlpath);
	    if (!err && addons) {
		update_addons_index(NULL);
	    } else if (!err && s != NULL && gui_callback != NULL) {
		pkg_install_invoke_callback(s, pkgname, dlpath, filetype);
	    }
	    g_free(dlpath);
	    return err;
	} else if (http) {
            /* get @fname as last portion of URL */
            const char *p = strrchr(pkgname, '/');

            if (p == NULL) {
                err = E_DATA;
            } else {
                fname = gretl_strdup(p + 1);
            }
        } else if (local) {
            const char *p;

            gretl_maybe_switch_dir(pkgname);
            p = strrslash(pkgname);
            if (p != NULL) {
                fname = gretl_strdup(p + 1);
            }
        }
    }

    if (!err && filetype) {
        const char *basename = fname != NULL ? fname : pkgname;
        const char *instpath;
        gchar *fullname;

        instpath = gretl_package_install_path(scripts ? "scripts" : "functions");
        fullname = g_strdup_printf("%s%s", instpath, basename);

        if (local) {
            err = gretl_copy_file(pkgname, fullname);
        } else if (http) {
            /* get file from a specified server */
            err = retrieve_public_file(pkgname, fullname);
        } else if (scripts) {
            /* get scripts collection from sourceforge */
            err = retrieve_remote_files_package(basename, fullname);
        } else {
            /* get function package from sourceforge */
            err = retrieve_remote_function_package(basename,
						   fullname,
						   staging);
        }

        if (!err && filetype == 2) {
            err = gretl_unzip_into(fullname, instpath);
            if (!err) {
                /* delete the zipfile */
                gretl_remove(fullname);
            }
        }

        if (!err && !scripts) {
	    package_check_dependencies(fullname, s, prn);
        }

        if (!err && gretl_messages_on()) {
            if (opt & OPT_D) {
                pprintf(prn, _("Installed dependency %s\n"), basename);
            } else {
                pprintf(prn, _("Installed %s\n"), basename);
            }
        }

        if (!err && s != NULL && gui_callback != NULL && !(opt & OPT_D)) {
            /* FIXME: handling of OPT_D (dependency) here? */
	    pkg_install_invoke_callback(s, basename, fullname, filetype);
        }

        g_free(fullname);
    }

    free(fname);

    return err;
}

#else /* !USE_CURL */

/* in this case we can only install from a local file */

static int install_package (const char *pkgname,
			    gretlopt opt,
			    ExecState *s,
			    PRN *prn)
{
    char *fname = NULL;
    int filetype = 0;
    int scripts = 0;
    int err = 0;

    if (!strncmp(pkgname, "http://", 7) ||
        !strncmp(pkgname, "https://", 8)) {
        gretl_errmsg_set(_("Internet access not supported"));
        return E_DATA;
    }

    if (strstr(pkgname, ".gfn")) {
        filetype = 1;
    } else if (strstr(pkgname, ".zip")) {
        filetype = 2;
        if (strstr(pkgname, "scripts")) {
            scripts = 1;
        }
    } else {
        /* must have suitable suffix */
        err = E_DATA;
    }

    if (!err) {
        /* get last portion of local filename */
        const char *p;

        gretl_maybe_switch_dir(pkgname);
        p = strrslash(pkgname);
        if (p != NULL) {
            fname = gretl_strdup(p + 1);
        }
    }

    if (!err && filetype) {
        const char *basename = fname != NULL ? fname : pkgname;
        const char *instpath;
        gchar *fullname;

        instpath = gretl_package_install_path(scripts ? "scripts" : "functions");
        fullname = g_strdup_printf("%s%s", instpath, basename);

        /* copy file into place */
        err = gretl_copy_file(pkgname, fullname);

        if (!err && filetype == 2) {
            err = gretl_unzip_into(fullname, instpath);
            if (!err) {
                /* delete the zipfile */
                gretl_remove(fullname);
            }
        }

        g_free(fullname);

        if (!err && gretl_messages_on()) {
            pprintf(prn, _("Installed %s\n"), basename);
        }
    }

    free(fname);

    return err;
}

#endif

static int abort_execution (ExecState *s)
{
    *s->cmd->savename = '\0';
    gretl_cmd_destroy_context(s->cmd);
    return E_STOP;
}

static int plot_ok;

void set_plot_produced (void)
{
    plot_ok = 1;
}

int is_plotting_command (CMD *cmd)
{
    if (GRAPHING_COMMAND(cmd->ci)) {
        return cmd->ci;
    } else if (cmd->ci == END && cmd->param != NULL) {
	int ci = gretl_command_number(cmd->param);

	if (ci == PLOT) {
	    return ci;
	}
    }

    return 0;
}

static void maybe_schedule_graph_callback (ExecState *s)
{
    int gui_mode = gretl_in_gui_mode();

    if (graph_written_to_file() || plot_output_to_buffer()) {
        if (gui_mode && *s->cmd->savename != '\0') {
            pprintf(s->prn, _("Warning: ignoring \"%s <-\"\n"), s->cmd->savename);
        }
	if (!plot_output_to_buffer()) {
	    report_plot_written(s->prn);
	}
    } else if (gui_mode && !graph_displayed()) {
        schedule_callback(s);
    }
}

static int execute_plot_call (CMD *cmd, DATASET *dset,
                              char *line, PRN *prn)
{
    gretlopt opt = cmd->opt;
    int err = 0;

    if (gretl_in_gui_mode() && *cmd->savename != '\0') {
        /* saving plot "as icon": add internal option to
           override production of a "gpttmp" file
        */
        opt |= OPT_G;
    }

    if (gretl_no_plots() && cmd->ci != END) {
        return 0;
    }

    if (cmd->ci == END) {
        /* end of a "plot" block */
        err = gretl_plot_finalize(line, dset, opt);
    } else if (opt & OPT_X) {
        err = matrix_command_driver(cmd->ci, cmd->list, cmd->param,
                                    dset, opt, prn);
    } else if (cmd->ci == GNUPLOT) {
        if (opt & (OPT_I | OPT_i)) {
            err = gnuplot_process_input(cmd->param, opt, prn);
        } else if (opt & OPT_C) {
            err = xy_plot_with_control(cmd->list, cmd->param,
                                       dset, opt);
        } else {
            err = gnuplot(cmd->list, cmd->param, dset, opt);
        }
    } else if (cmd->ci == SCATTERS) {
        err = multi_plots(cmd->list, dset, opt);
    } else if (cmd->ci == BXPLOT) {
        err = boxplots(cmd->list, cmd->param, dset, opt);
    } else if (cmd->ci == HFPLOT) {
        err = hf_plot(cmd->list, cmd->param, dset, opt);
    } else if (cmd->ci == PANPLOT) {
        err = cli_panel_plot(cmd->list, cmd->param, dset, opt);
    } else if (cmd->ci == QQPLOT) {
        err = qq_plot(cmd->list, dset, opt);
    } else if (cmd->ci == KDPLOT) {
	err = kd_plot(cmd->list, dset, opt);
    }

    return err;
}

static int execute_gridplot_call (CMD *cmd, PRN *prn)
{
    gretlopt opt = cmd->opt;

#if 0 /* maybe try to enable this? */
    if (gretl_in_gui_mode() && *cmd->savename != '\0') {
        /* saving multiplot "as icon": add internal option to
           override production of a "gpttmp" file
        */
        opt |= OPT_G;
    }
#endif

    if (gretl_no_plots()) {
        return 0;
    }

    return gretl_gridplot_finalize(opt);
}

static int smpl_restrict (gretlopt opt)
{
    opt &= ~OPT_Q;
    opt &= ~OPT_T;
    opt &= ~OPT_D;
    return opt != OPT_NONE;
}

static int check_smpl_full (gretlopt opt)
{
    opt ^= OPT_F;
    if (opt & OPT_Q) {
        opt ^= OPT_Q;
    }
    return opt == OPT_NONE ? 0 : E_BADOPT;
}

static int panel_smpl_special (gretlopt opt)
{
    /* the option --unit or --time is given */
    return opt & (OPT_U | OPT_X);
}

static int cant_do_smpl (DATASET *dset)
{
    if (dset != NULL && dset->n > 0) {
        return 0; /* OK */
    } else {
        const char *db = get_db_name();

        return *db == '\0';
    }
}

static void maybe_print_error_message (CMD *cmd, int err, PRN *prn)
{
    if (gretl_function_depth() > 0) {
        ; /* defer printing */
    } else if (cmd->flags & CMD_CATCH) {
        /* print only if messages on */
        if (gretl_messages_on()) {
            errmsg(err, prn);
        }
    } else {
        /* otherwise go ahead and print */
        errmsg(err, prn);
    }
}

int gretl_cmd_exec (ExecState *s, DATASET *dset)
{
    CMD *cmd = s->cmd;
    char *line = s->line;
    MODEL *model = s->model;
    PRN *prn = s->prn;
    char readfile[MAXLEN];
    int *listcpy = NULL;
    int err = 0;

    exec_state_prep(s);
    plot_ok = 0;

    if (get_user_stop()) {
        /* the user called a halt to execution */
        return abort_execution(s);
    }

    if (NEEDS_MODEL_CHECK(cmd->ci)) {
        err = model_test_check(cmd, dset, prn);
    } else if (MODIFIES_LIST(cmd->ci)) {
        if (cmd->list[0] == 0) {
            /* no-op */
            return 0;
        } else {
            /* list is potentially modified -> make a copy */
            listcpy = gretl_list_copy(cmd->list);
            if (listcpy == NULL) {
                err = E_ALLOC;
            }
        }
    }

    if (!err) {
	/* the --buffer option (OPT_b) is always incompatible
	   with --output (OPT_U)
	*/
	err = incompatible_options(cmd->opt, OPT_b | OPT_U);
    }

    if (err) {
        goto bailout;
    }

    *readfile = '\0';

    if (cmd->ci == OLS && dataset_is_panel(dset)) {
        cmd->ci = PANEL;
        cmd->opt |= (OPT_P | OPT_L); /* panel pooled OLS flag */
    }

#if 0
    fprintf(stderr, "gretl_cmd_exec: '%s' (ci %d, %s) \n", line,
	    cmd->ci, gretl_command_word(cmd->ci));
#endif

    switch (cmd->ci) {

    case APPEND:
    case JOIN:
    case OPEN:
	if (cmd->ci == JOIN && (dset == NULL || dset->v == 0)) {
	    /* "join" is not applicable in the absence of
	       a dataset */
	    err = E_NODATA;
	} else {
	    err = lib_open_append(s, dset, readfile, prn);
	}
        if (!err && cmd->ci != OPEN) {
            schedule_callback(s);
        }
        break;

    case CLEAR:
        err = check_clear(cmd->opt);
        if (!err) {
	    if (gretl_in_gui_mode()) {
		schedule_callback(s);
	    } else if (cmd->opt & OPT_A) {
		gretl_functions_cleanup();
		lib_clear_data(s, dset);
	    } else if (cmd->opt & OPT_F) {
		gretl_functions_cleanup();
	    } else {
                lib_clear_data(s, dset);
            }
        }
        break;

    case FLUSH:
        if (gretl_in_gui_mode()) {
            schedule_callback(s);
        } else {
            gretl_print_flush_stream(prn);
        }
        break;

    case ANOVA:
        err = anova(cmd->list, dset, cmd->opt, prn);
        break;

    case ADF:
        err = adf_test(cmd->order, cmd->list, dset, cmd->opt, prn);
        break;

    case KPSS:
        err = kpss_test(cmd->order, cmd->list, dset, cmd->opt, prn);
        break;

    case LEVINLIN:
        err = llc_test_driver(cmd->param, cmd->list, dset,
                              cmd->opt, prn);
        break;

    case COINT:
        err = engle_granger_test(cmd->order, cmd->list, dset,
                                 cmd->opt, prn);
        break;

    case COINT2:
        err = johansen_test_simple(cmd->order, cmd->list,
                                   dset, cmd->opt, prn);
        break;

    case CORR:
        err = incompatible_options(cmd->opt, OPT_N | OPT_S | OPT_K);
        if (!err) {
            err = incompatible_options(cmd->opt, OPT_X | OPT_S | OPT_K);
        }
        if (err) {
            break;
        }
        if (cmd->opt & OPT_K) {
            err = kendall_tau(cmd->list, dset, cmd->opt, prn);
        } else if (cmd->opt & OPT_S) {
            err = spearman_rho(cmd->list, dset, cmd->opt, prn);
        } else if (cmd->opt & OPT_X) {
            err = matrix_command_driver(CORR, cmd->list, cmd->param,
                                        dset, cmd->opt, prn);
        } else {
            err = gretl_corrmx(cmd->list, dset, cmd->opt, prn);
        }
        break;

    case CORRGM:
        err = corrgram(cmd->list[1], cmd->order, 0, dset,
                       cmd->opt, prn);
        break;

    case XCORRGM:
        err = xcorrgram(cmd->list, cmd->order, dset,
                        cmd->opt, prn);
        break;

    case PERGM:
        err = periodogram(cmd->list[1], cmd->order, dset,
                          cmd->opt, prn);
        break;

    case FRACTINT:
        err = fractint(cmd->list[1], cmd->order, dset,
                       cmd->opt, prn);
        break;

    case BREAK:
    case CONTINUE:
    case ENDLOOP:
        pprintf(prn, _("You can't end a loop here, "
                       "you haven't started one\n"));
        err = 1;
        break;

    case FCAST:
        err = do_forecast(cmd->param, dset, cmd->opt, prn);
        break;

    case FREQ:
        if (cmd->opt & OPT_X) {
            err = matrix_freq_driver(cmd->list, cmd->opt, prn);
        } else {
            err = freqdist(cmd->list[1], dset, cmd->opt, prn);
        }
        break;

    case DISCRETE:
        err = list_makediscrete(cmd->list, dset, cmd->opt);
        break;

    case ESTIMATE:
        err = estimate_named_system(cmd->param, cmd->parm2, dset,
                                    cmd->opt, prn);
        if (!err && (cmd->opt & OPT_W) && gretl_in_gui_mode()) {
            gui_save_system(s);
        }
        break;

    case FUNC:
        err = gretl_start_compiling_function(cmd->param, dset, prn);
        break;

    case GENR:
    case EVAL:
        if (cmd->flags & CMD_CATCH) {
            err = generate(cmd->vstart, dset, cmd->gtype,
                           cmd->opt | OPT_C, prn);
            if (err == E_BADCATCH) {
                cmd->flags ^= CMD_CATCH;
            }
        } else {
            err = generate(cmd->vstart, dset, cmd->gtype,
                           cmd->opt, prn);
        }
        break;

    case PCA:
        err = do_pca(cmd->list, dset, cmd->opt, prn);
        break;

    case DATA:
        err = db_get_series(cmd->param, dset, cmd->opt, prn);
        break;

    case DATAMOD:
        err = modify_dataset(dset, cmd->auxint, cmd->list,
                             cmd->parm2, cmd->opt, prn);
        if (!err) {
            schedule_callback(s);
        }
        break;

    case DIFF:
    case LDIFF:
    case SDIFF:
        err = list_diffgenr(listcpy, cmd->ci, dset);
        if (!err) {
            maybe_list_series(dset, prn);
            set_dataset_is_changed(dset, 1);
        }
        break;

    case ORTHDEV:
        err = list_orthdev(listcpy, dset);
        if (!err) {
            maybe_list_series(dset, prn);
            set_dataset_is_changed(dset, 1);
        }
        break;

    case DUMMIFY:
        err = list_dumgenr(&listcpy, dset, cmd->opt);
        if (!err) {
            maybe_list_series(dset, prn);
            set_dataset_is_changed(dset, 1);
        }
        break;

    case LAGS:
        err = list_laggenr(&listcpy, 1, cmd->order, NULL,
                           dset, 0, cmd->opt);
        if (!err) {
            maybe_list_series(dset, prn);
            set_dataset_is_changed(dset, 1);
        }
        break;

    case LOGS:
        err = list_loggenr(listcpy, dset);
        if (!err) {
            maybe_list_series(dset, prn);
            set_dataset_is_changed(dset, 1);
        }
        break;

    case STDIZE:
        err = list_stdgenr(listcpy, dset, cmd->opt);
        if (!err) {
            maybe_list_series(dset, prn);
            set_dataset_is_changed(dset, 1);
        }
        break;

    case SQUARE:
        err = list_xpxgenr(&listcpy, dset, cmd->opt);
        if (!err) {
            maybe_list_series(dset, prn);
            set_dataset_is_changed(dset, 1);
        }
        break;

    case TEXTPLOT:
        err = textplot(cmd->list, dset, cmd->opt, prn);
        break;

    case RMPLOT:
        err = rmplot(cmd->list, dset, cmd->opt, prn);
        break;

    case HURST:
        err = hurstplot(cmd->list, dset, cmd->opt, prn);
        break;

    case INFO:
        if (cmd->opt & (OPT_F | OPT_T)) {
            err = read_or_write_dset_description(cmd->opt, dset, prn);
        } else {
            print_info(cmd->opt, dset, prn);
        }
        break;

    case RENAME:
        err = rename_series(dset, cmd->auxint, cmd->parm2, cmd->opt);
        if (!err && !(cmd->opt & OPT_Q)) {
            maybe_list_series(dset, prn);
        }
        break;

    case SET:
        err = execute_set(cmd->param, cmd->parm2, dset, cmd->opt, prn);
	if (!err && cmd->parm2 != NULL) {
	    maybe_schedule_set_callback(s);
	}
        break;

    case SETINFO:
        err = set_var_info(cmd->list, cmd->param, cmd->parm2,
                           cmd->opt, dset);
        break;

    case SETMISS:
        err = set_miss(cmd->list, cmd->param, dset, prn);
        break;

    case LABELS:
        if (cmd->opt & (OPT_D | OPT_F | OPT_T | OPT_A | OPT_R)) {
            err = read_or_write_var_labels(cmd->opt, dset, prn);
            if (!err && (cmd->opt & (OPT_D | OPT_F))) {
                schedule_callback(s);
            }
        } else {
            showlabels(cmd->list, cmd->opt, dset, prn);
        }
        break;

    case MARKERS:
        err = read_or_write_obs_markers(cmd->opt, dset, prn);
        if (!err && (cmd->opt & (OPT_D | OPT_F))) {
            schedule_callback(s);
        }
        break;

    case VARLIST:
        if (cmd->opt & OPT_T) {
            list_user_vars_of_type(dset, prn);
        } else if (cmd->opt & OPT_S) {
            print_scalars(prn);
        } else if (cmd->opt & OPT_A) {
            list_ok_dollar_vars(dset, prn);
        } else {
            list_series(dset, cmd->opt, prn);
        }
        break;

    case PRINT:
        if (cmd->opt & OPT_L) {
            err = printdata(NULL, cmd->param, dset, OPT_NONE, prn);
        } else if (cmd->param != NULL) {
            /* directly printing a string literal */
            pputs(prn, cmd->param);
            pputc(prn, '\n');
        } else {
            err = printdata(cmd->list, cmd->parm2, dset, cmd->opt, prn);
        }
        break;

    case PRINTF:
        err = do_printf_command(cmd->param, cmd->vstart, dset,
				prn, cmd_arg1_quoted(cmd));
        break;

    case PVAL:
        err = batch_pvalue(cmd->param, dset, prn);
        break;

    case SUMMARY:
        err = incompatible_options(cmd->opt, OPT_B | OPT_W);
        if (err) {
            break;
        }
        if (cmd->opt & OPT_B) {
            err = summary_statistics_by(cmd->list, dset, cmd->opt, prn);
        } else if (cmd->opt & OPT_X) {
            err = matrix_command_driver(cmd->ci, cmd->list, cmd->param,
                                        dset, cmd->opt, prn);
        } else {
            err = list_summary_driver(cmd->list, dset, cmd->opt, prn);
        }
        break;

    case XTAB:
        if (cmd->opt & OPT_X) {
            err = crosstab_from_matrix(cmd->opt, prn);
        } else {
            err = crosstab(cmd->list, dset, cmd->opt, prn);
        }
        break;

    case MAHAL:
        err = mahalanobis_distance(cmd->list, dset, cmd->opt, prn);
        break;

    case MEANTEST:
        err = means_test(cmd->list, dset, cmd->opt, prn);
        break;

    case VARTEST:
        err = vars_test(cmd->list, dset, cmd->opt, prn);
        break;

    case RUNS:
        err = runs_test(cmd->list[1], dset, cmd->opt, prn);
        break;

    case SPEARMAN:
        err = spearman_rho(cmd->list, dset, cmd->opt, prn);
        break;

    case DIFFTEST:
        err = diff_test(cmd->list, dset, cmd->opt, prn);
        break;

    case OUTFILE:
        err = do_outfile_command(cmd->opt, cmd->param, dset, prn);
        break;

    case SETOBS:
        err = set_obs(cmd->param, cmd->parm2, dset, cmd->opt);
        if (!err) {
            if (dset->n > 0) {
                if (!(cmd->opt & (OPT_I | OPT_G))) {
                    print_smpl(dset, 0, OPT_NONE, prn);
                }
                schedule_callback(s);
            } else {
                pprintf(prn, _("data frequency = %d\n"), dset->pd);
            }
        }
        break;

    case SETOPT:
        err = set_options_for_command(cmd->param, cmd->parm2, cmd->opt);
        if (!err && gretl_messages_on()) {
            pprintf(prn, _("Set option(s) for command \"%s\"\n"), cmd->param);
        }
        break;

    case SMPL:
        if (cant_do_smpl(dset)) {
            err = E_NODATA;
        } else if (cmd->opt & OPT_F) {
            err = check_smpl_full(cmd->opt);
            if (!err) {
                err = restore_full_sample(dset, s);
            }
        } else if (cmd->opt == OPT_T && cmd->param == NULL) {
            /* --permanent, by itself */
            err = perma_sample(dset, cmd->opt, prn, NULL);
        } else if (panel_smpl_special(cmd->opt)) {
            /* panel data: --unit or --time option */
            err = set_panel_sample(cmd->param, cmd->parm2, cmd->opt, dset,
				   s, prn);
        } else if (smpl_restrict(cmd->opt)) {
            /* --restrict, --dummy, etc. */
            err = restrict_sample(cmd->param, cmd->list, dset,
                                  s, cmd->opt, prn, NULL);
        } else if (cmd->param == NULL && cmd->parm2 == NULL) {
            /* no args given: give a report */
            print_smpl(dset, get_full_length_n(), OPT_F, prn);
            break;
        } else {
            /* simple setting of t1, t2 business */
            err = set_sample(cmd->param, cmd->parm2, dset, cmd->opt);
        }
        if (!err && gretl_messages_on() && !(cmd->opt & OPT_Q)) {
            print_smpl(dset, get_full_length_n(), OPT_NONE, prn);
        }
        break;

    case PKG:
        err = do_pkg_command(readfile, dset, s, prn);
        break;

    case MAKEPKG:
        err = create_and_write_function_package(cmd->param, cmd->opt, prn);
        break;

    case STORE:
        if (!has_param(cmd)) {
            pputs(prn, _("store: no filename given\n"));
            err = E_PARSE;
        } else if (cmd->opt & OPT_A) {
            err = write_matrix_as_dataset(cmd->param, cmd->opt, prn);
        } else if (cmd->opt & OPT_U) {
            err = write_loaded_functions_file(cmd->param, cmd->opt);
        } else if (dset == NULL || dset->Z == NULL) {
            err = E_NODATA;
        } else {
            err = write_data(cmd->param, cmd->list, dset, cmd->opt, prn);
        }
        break;

    case SHELL:
        err = gretl_shell(cmd->vstart, cmd->opt, prn);
        break;

    case OLS:
    case WLS:
        clear_model(model);
        *model = lsq(cmd->list, dset, cmd->ci, cmd->opt);
        err = print_save_model(model, dset, cmd->opt, 0, prn, s);
        break;

    case MPOLS:
        clear_model(model);
        *model = mp_ols(cmd->list, dset, cmd->opt);
        err = print_save_model(model, dset, cmd->opt, 0, prn, s);
        break;

    case AR:
    case AR1:
    case ARMA:
    case ARCH:
        clear_model(model);
        if (cmd->ci == AR) {
            *model = ar_model(cmd->list, dset, cmd->opt, prn);
        } else if (cmd->ci == AR1) {
            *model = ar1_model(cmd->list, dset, cmd->opt, prn);
        } else if (cmd->ci == ARMA) {
            *model = arma(cmd->list, cmd->auxlist, dset, cmd->opt, prn);
        } else {
            *model = arch_model(cmd->list, cmd->order, dset, cmd->opt);
        }
        if (cmd->ci == ARMA && (cmd->opt & OPT_Z)) {
            ; /* no real MODEL is returned */
        } else {
            err = print_save_model(model, dset, cmd->opt, 0, prn, s);
        }
        break;

    case PANEL:
    case DPANEL:
        if (!dataset_is_panel(dset)) {
            gretl_errmsg_set(_("This estimator requires panel data"));
            err = E_DATA;
            break;
        }
        /* Falls through. */
    case GARCH:
    case HECKIT:
    case HSK:
    case INTREG:
    case IVREG:
    case LAD:
    case LOGISTIC:
    case LOGIT:
    case POISSON:
    case NEGBIN:
    case PROBIT:
    case QUANTREG:
    case TOBIT:
    case DURATION:
    case BIPROBIT:
    case MIDASREG:
        clear_model(model);
        if (cmd->ci == LOGIT || cmd->ci == PROBIT) {
            *model = logit_probit(cmd->list, dset, cmd->ci, cmd->opt, prn);
        } else if (cmd->ci == HSK) {
            *model = hsk_model(cmd->list, dset, cmd->opt);
        } else if (cmd->ci == LOGISTIC) {
            *model = logistic_driver(cmd->list, dset, cmd->opt);
        } else if (cmd->ci == TOBIT) {
            *model = tobit_driver(cmd->list, dset, cmd->opt, prn);
        } else if (cmd->ci == POISSON || cmd->ci == NEGBIN) {
            *model = count_model(cmd->list, cmd->ci, dset, cmd->opt, prn);
        } else if (cmd->ci == HECKIT) {
            *model = heckit_model(cmd->list, dset, cmd->opt, prn);
        } else if (cmd->ci == IVREG) {
            *model = ivreg(cmd->list, dset, cmd->opt);
        } else if (cmd->ci == LAD) {
            *model = lad_model(cmd->list, dset, cmd->opt);
        } else if (cmd->ci == QUANTREG) {
            *model = quantreg_driver(cmd->param, cmd->list, dset, cmd->opt, prn);
        } else if (cmd->ci == DURATION) {
            *model = duration_model(cmd->list, dset, cmd->opt, prn);
        } else if (cmd->ci == GARCH) {
            *model = garch(cmd->list, dset, cmd->opt, prn);
        } else if (cmd->ci == PANEL) {
            *model = panel_model(cmd->list, dset, cmd->opt, prn);
        } else if (cmd->ci == DPANEL) {
            *model = dpd_model(cmd->list, cmd->auxlist, cmd->param, dset, cmd->opt, prn);
        } else if (cmd->ci == INTREG) {
            *model = interval_model(cmd->list, dset, cmd->opt, prn);
        } else if (cmd->ci == BIPROBIT) {
            *model = biprobit_model(cmd->list, dset, cmd->opt, prn);
        } else if (cmd->ci == MIDASREG) {
            *model = midas_model(cmd->list, cmd->param, dset, cmd->opt, prn);
        } else {
            /* can't happen */
            err = 1;
            break;
        }
        err = print_save_model(model, dset, cmd->opt, 0, prn, s);
        break;

    case GMM:
    case MLE:
    case NLS:
        err = nl_parse_line(cmd->ci, cmd->vstart, dset, prn);
        if (!err) {
            gretl_cmd_set_context(cmd, cmd->ci);
        }
        break;

    case FOREIGN:
    case MPI:
        if (cmd->context == FOREIGN || cmd->context == MPI) {
            err = foreign_append(line, cmd->context);
        } else {
            err = foreign_start(cmd->ci, cmd->param, cmd->opt, prn);
            if (!err) {
                gretl_cmd_set_context(cmd, cmd->ci);
            }
        }
        break;

    case GIBBS:
        if (cmd->context == GIBBS) {
            err = gibbs_block_append(line);
        } else {
            err = gibbs_block_start(cmd->vstart, prn);
            if (!err) {
                gretl_cmd_set_context(cmd, cmd->ci);
            }
        }
        break;

    case PLOT:
        if (!cmd->context) {
            err = gretl_plot_start(cmd->param, dset);
        } else {
            err = gretl_plot_append_line(line, dset);
        }
        if (!err && !cmd->context) {
            gretl_cmd_set_context(cmd, cmd->ci);
        }
        break;

    case GPBUILD:
    case GRIDPLOT:
	if (gretl_gridplot_collecting()) {
	    gretl_errmsg_set(_("gpbuild/gridplot: cannot be nested"));
	    err = E_DATA;
	} else if (cmd->ci == GRIDPLOT) {
	    err = check_gridplot_options(cmd->opt);
	}
	if (!err) {
	    if (cmd->ci == GPBUILD) {
		/* the block-start case */
		err = gretl_gridplot_start(cmd->param, cmd->opt, dset);
	    } else {
		/* the gridplot case */
		err = gretl_gridplot_from_array(cmd->param, cmd->opt);
	    }
	}
        break;

    case ADD:
    case OMIT:
        if (get_last_model_type() == GRETL_OBJ_VAR) {
            err = VAR_omit_driver(cmd, dset, prn);
        } else if (add_omit_save(cmd, model)) {
            MODEL mymod;

            gretl_model_init(&mymod, dset);
            if (cmd->ci == ADD) {
                err = add_test_full(model, &mymod, cmd->list,
                                    dset, cmd->opt, prn);
            } else {
                err = omit_test_full(model, &mymod, cmd->list,
                                     dset, cmd->opt, prn);
            }
            if (!err && mymod.ncoeff > 0) {
                gretlopt popt = OPT_NONE;

                if (cmd->opt & (OPT_I | OPT_Q)) {
                    popt = OPT_Q;
                } else if (cmd->opt & OPT_O) {
                    popt = OPT_O; /* --vcv printing option */
                }
                clear_model(model);
                *model = mymod;
                print_save_model(model, dset, popt, 1, prn, s);
            }
        } else if (cmd->ci == ADD) {
            err = add_test(model, cmd->list, dset, cmd->opt, prn);
        } else {
            err = omit_test(model, cmd->list, dset, cmd->opt, prn);
        }
        if (err == E_NOOMIT) {
            /* auto-omit was a no-op */
            err = 0;
        }
        break;

    case BKW:
    case COEFFSUM:
    case CUSUM:
    case RESET:
    case CHOW:
    case QLRTEST:
    case VIF:
        if (cmd->ci == COEFFSUM) {
            err = gretl_sum_test(cmd->list, model, dset, cmd->opt, prn);
        } else if (cmd->ci == CUSUM) {
            err = cusum_test(model, dset, cmd->opt, prn);
        } else if (cmd->ci == RESET) {
            err = reset_test(model, dset, cmd->opt, prn);
        } else if (cmd->ci == CHOW) {
            err = chow_test_driver(cmd->param, model, dset, cmd->opt, prn);
        } else if (cmd->ci == QLRTEST) {
            err = QLR_test(model, dset, cmd->opt, prn);
        } else if (cmd->ci == VIF) {
            err = vif_test(model, dset, cmd->opt, prn);
        } else if (cmd->ci == BKW) {
            err = bkw_test(model, dset, cmd->opt, prn);
        }
        break;

    case BDS:
        err = bds_test_driver(cmd->order, cmd->list, dset,
			      cmd->opt, prn);
        break;

    case NORMTEST:
        err = gretl_normality_test(cmd->list[1], dset, cmd->opt, prn);
        break;

    case PANSPEC:
        err = panel_specification_test(model, dset, cmd->opt, prn);
        break;

    case MODTEST:
        err = model_test_driver(cmd->order, dset, cmd->opt, prn);
        break;

    case LEVERAGE:
        err = leverage_test(model, dset, cmd->opt, prn);
        if (!err && (cmd->opt & OPT_S) && !(cmd->opt & OPT_Q)) {
            maybe_list_series(dset, prn);
        }
        break;

    case EQNPRINT:
    case TABPRINT:
        if (model->errcode == E_NAN) {
            pprintf(prn, _("Couldn't format model\n"));
        } else {
            err = model_print_driver(model, dset, cmd->ci,
                                     cmd->param, cmd->opt,
                                     prn);
        }
        break;

    case RESTRICT:
        /* joint hypothesis test on model */
        if (s->rset == NULL) {
            if (!has_param(cmd)) {
                /* if param is non-blank, we're restricting a named system */
                err = model_test_check(cmd, dset, prn);
                if (err) break;
            }
            s->rset = restriction_set_start(cmd->param, cmd->opt, &err);
            if (!err) {
                /* FIXME redundant? */
                gretl_cmd_set_context(cmd, RESTRICT);
            }
        } else {
            err = restriction_set_parse_line(s->rset, line, dset);
            if (err) {
                s->rset = NULL;
            }
        }
        break;

    case SYSTEM:
        if (s->sys == NULL) {
            /* no equation system is defined currently */
            s->sys = equation_system_start(cmd->param, cmd->savename,
                                           cmd->opt, &err);
            if (!err) {
                gretl_cmd_set_context(cmd, SYSTEM);
            }
        } else {
            /* tokenize: use of @line OK here? */
            err = system_parse_line(s->sys, line, dset);
            if (err) {
                s->sys = NULL;
            }
        }
        break;

    case EQUATION:
        if (cmd->opt & OPT_M) {
            err = equation_system_append_multi(s->sys, cmd->param,
                                               cmd->parm2, dset);
        } else {
            err = equation_system_append(s->sys, cmd->list);
        }
        if (err) {
            s->sys = NULL;
        }
        break;

    case END:
        if (!strcmp(cmd->param, "system")) {
            err = equation_system_finalize(s->sys, dset, cmd->opt, prn);
            if (!err) {
                gui_save_system(s);
            }
            /* clear for next use */
            s->sys = NULL;
        } else if (!strcmp(cmd->param, "mle") ||
                   !strcmp(cmd->param, "nls") ||
                   !strcmp(cmd->param, "gmm")) {
            clear_model(model);
            *model = nl_model(dset, cmd->opt, prn);
            err = print_save_model(model, dset, cmd->opt, 0, prn, s);
        } else if (!strcmp(cmd->param, "restrict")) {
            err = do_end_restrict(s, dset);
        } else if (!strcmp(cmd->param, "foreign")) {
            err = foreign_execute(dset, cmd->opt, prn);
        } else if (!strcmp(cmd->param, "mpi")) {
            err = foreign_execute(dset, cmd->opt, prn);
        } else if (!strcmp(cmd->param, "plot")) {
            err = execute_plot_call(cmd, dset, line, prn);
        } else if (!strcmp(cmd->param, "outfile")) {
            err = close_outfile(prn);
        } else if (!strcmp(cmd->param, "gpbuild")) {
            err = execute_gridplot_call(cmd, prn);
        } else if (!strcmp(cmd->param, "gibbs")) {
            err = gibbs_execute(cmd->opt, prn);
        } else {
            err = E_PARSE;
        }
        break;

    case VAR:
    case VECM:
        if (cmd->ci == VAR) {
            s->var = gretl_VAR(cmd->order, cmd->auxlist, cmd->list,
                               dset, cmd->opt, prn, &err);
        } else {
            s->var = gretl_VECM(cmd->order, cmd->auxint, cmd->list,
                                dset, cmd->opt, prn, &err);
        }
        if (!err && s->var != NULL) {
            save_var_vecm(s);
        }
        break;

    case RUN:
    case INCLUDE:
        if (cmd->ci == RUN) {
            err = get_full_filename(cmd->param, readfile, OPT_S);
        } else {
            err = get_full_filename(cmd->param, readfile, OPT_I);
            cmd->opt |= OPT_Q;
        }
        if (err) {
            break;
        }
        if (gretl_messages_on()) {
            pprintf(prn, " %s\n", readfile);
        }
        if (cmd->ci == INCLUDE && gretl_is_xml_file(readfile)) {
            err = load_XML_functions_file(readfile, cmd->opt, prn);
            break;
        } else if (cmd->ci == INCLUDE && gfn_is_loaded(readfile)) {
            break;
        }
        if (!strcmp(readfile, s->runfile)) {
            pprintf(prn, _("Infinite loop detected in script\n"));
            err = 1;
            break;
        }
        err = run_script(readfile, s, dset, cmd->opt, prn);
        break;

    case FUNCERR:
    case FUNCRET:
        if (gretl_function_depth() == 0) {
            gretl_errmsg_sprintf(_("'%s': can only be used within a function"),
                                 gretl_command_word(cmd->ci));
            err = 1;
        } else if (cmd->ci == FUNCERR) {
            err = E_FUNCERR;
        }
        break;

    case DELEET:
        err = gretl_delete_variables(cmd->list, cmd->param, cmd->opt,
                                     dset, NULL, prn);
        break;

    case MODPRINT:
        err = do_modprint(cmd->param, cmd->parm2, cmd->opt, prn);
        break;

    case GNUPLOT:
    case BXPLOT:
    case SCATTERS:
    case HFPLOT:
    case PANPLOT:
    case QQPLOT:
    case KDPLOT:
        err = execute_plot_call(cmd, dset, NULL, prn);
        break;

    case MODELTAB:
    case GRAPHPG:
        if (gretl_in_gui_mode()) {
            schedule_callback(s);
        } else {
            pprintf(prn, _("%s: command not available\n"),
                    gretl_command_word(cmd->ci));
        }
        break;

    default:
        {
            const char *word = gretl_command_word(cmd->ci);

            if (*word != '\0') {
                pprintf(prn, _("The \"%s\" command cannot be used in this context"),
                        word);
                pputc(prn, '\n');
            } else {
                pprintf(prn, "What??\n");
            }
        }
        err = 1;
        break;
    }

    if (listcpy != NULL) {
        free(listcpy);
    }

    if (err == E_OK) {
        err = 0;
    }

    if (!err && plot_ok) {
        maybe_schedule_graph_callback(s);
        plot_ok = 0;
    }

    if (callback_scheduled(s)) {
        callback_exec(s, readfile, err);
    }

 bailout:

    if (err) {
        maybe_print_error_message(cmd, err, prn);
        err = process_command_error(s, err);
    }

    if (err) {
        gretl_cmd_destroy_context(cmd);
    } else {
        /* this is a no-op if there's no warning */
        warnmsg(prn);
    }

    return err;
}

/**
 * get_command_index:
 * @s: pointer to execution state
 * @cmode: compilation mode: LOOP or FUNC
  *
 * Parse @line and assign to the %ci field of @s->cmd the index number of
 * the command embedded in @s->line.  Note: this is a "lite" version of
 * parse_command_line().  It is used when commands are being stacked
 * for execution within a loop. Command options are not parsed out.
 *
 * Returns: 1 on error, otherwise 0.
 */

int get_command_index (ExecState *s, int cmode, int preserve)
{
    CMD *cmd = s->cmd;
    char *line = s->line;
    int err = 0;

    gretl_cmd_clear(cmd);

#if CMD_DEBUG
    fprintf(stderr, "get_command_index: line='%s'\n", line);
#endif

    if ((cmd->context == FOREIGN || cmd->context == MPI) &&
        !ends_foreign_block(line)) {
        cmd->opt = OPT_NONE;
        cmd->ci = cmd->context;
        return 0;
    }

    if (filter_comments(line, cmd, preserve)) {
        return 0;
    }

    err = real_parse_command(s, NULL, cmode, NULL);

    if (!err && cmd->ci == 0) {
        /* maybe genr? */
        const char *s = cmd->toks[0].s;

        if (s != NULL) {
            if (*s == '$' || *s == '@') s++;
            if (strlen(s) == gretl_namechar_spn(s)) {
                cmd->ci = GENR;
            }
        }
    }

    if (!err && cmd->ci == 0) {
        /* FIXME watch out for fallout! (2012-03-01) */
        cmd->ci = CMD_NULL;
        err = E_PARSE;
    }

    if (err) {
        return err;
    }

    if (cmd->ci == END) {
        cmd->context = 0;
    } else if (cmd->context) {
        cmd->ci = cmd->context;
    }

    if (cmd->ci == NLS || cmd->ci == MLE ||
        cmd->ci == GMM || cmd->ci == FOREIGN ||
        cmd->ci == PLOT || cmd->ci == MPI) {
        cmd->context = cmd->ci;
    }

#if CMD_DEBUG
    fprintf(stderr, " get_command_index: cmd->ci set to %d (%s)\n",
	    cmd->ci, gretl_command_word(cmd->ci));
#endif

    return 0;
}

void gretl_cmd_set_context (CMD *cmd, int ci)
{
    cmd->context = ci;
}

void gretl_cmd_destroy_context (CMD *cmd)
{
    if (cmd->context == FOREIGN || cmd->context == MPI) {
        foreign_destroy();
    }
    cmd->context = 0;
    *cmd->savename = '\0';
}

gretlopt gretl_cmd_get_opt (const CMD *cmd)
{
    return cmd->opt;
}

void gretl_cmd_set_opt (CMD *cmd, gretlopt opt)
{
    cmd->opt = opt;
}

const char *gretl_cmd_get_savename (CMD *cmd)
{
    return cmd->savename;
}

void function_state_init (CMD *cmd, ExecState *state, int *indent0)
{
    cmd->list = NULL;
    cmd->auxlist = NULL;
    cmd->param = NULL;
    cmd->parm2 = NULL;
    /* FIXME tokenize more needed? */

    state->cmd = NULL;
    state->model = NULL;
    state->submask = NULL;

    state->padded = 0;

    *indent0 = gretl_if_state_record();
}

void gretl_exec_state_set_callback (ExecState *s,
                                    EXEC_CALLBACK callback,
                                    gretlopt opt)
{
    s->callback = callback;
    s->pmod = NULL;

    if (opt & OPT_G) {
        gui_callback = callback;
    }
}

EXEC_CALLBACK get_gui_callback (void)
{
    return gui_callback;
}

void set_gui_callback (EXEC_CALLBACK callback)
{
    gui_callback = callback;
}

void gretl_exec_state_clear (ExecState *s)
{
    gretl_cmd_free(s->cmd);

    if (s->free_line) {
	free(s->line);
	s->line = NULL;
    }

    if (s->flags & FUNCTION_EXEC) {
        /* Restore whatever was the 'last model' before
           function execution. Note that this includes
           the case where there was no 'last model', in
           which case we restore the null state. Drop
           the extra refcount for the model we put into
           last model position (if any), so we don't end
           up leaking memory.
        */
        set_as_last_model(s->prev_model, s->prev_type);
        if (s->prev_model != NULL) {
            gretl_object_unref(s->prev_model, s->prev_type);
        }
#if PMDEBUG
        fprintf(stderr, "ExecState %p, set prev_model %p as last_model\n",
                (void *) s, (void *) s->prev_model);
#endif
        /* restore the previous model count */
        if (s->prev_model_count >= 0) {
            set_model_count(s->prev_model_count);
        }
    }

    destroy_working_model(s->model);

    s->prev_model = NULL;
    s->prev_type = GRETL_OBJ_NULL;
    s->prev_model_count = -1;

    free_subsample_mask(s->submask);
}

void gretl_exec_state_destroy (ExecState *s)
{
    free_subsample_mask(s->submask);
    gretl_abort_compiling_loop();
    free(s);
}

void gretl_exec_state_uncomment (ExecState *s)
{
    s->in_comment = 0;
    s->cmd->flags &= ~CMD_CCMT;
}

void gretl_exec_state_transcribe_flags (ExecState *s, CMD *cmd)
{
    s->in_comment = (cmd_ccmt(cmd))? 1 : 0;
}

void gretl_exec_state_set_model (ExecState *s, MODEL *pmod)
{
    s->pmod = pmod;
}

int process_command_error (ExecState *s, int err)
{
    int ret = err;

    if (err) {
        if (gretl_compiling_function() ||
            gretl_compiling_loop()) {
            ; /* pass the error through */
        } else if (s->cmd->flags & CMD_CATCH) {
            /* local "continue on error" */
            set_gretl_errno(err);
            s->cmd->flags ^= CMD_CATCH;
            ret = 0;
        }
    }

    if (ret) {
	if (gretl_gridplot_collecting()) {
	    gretl_gridplot_destroy();
	}
	if (print_redirection_level(s->prn) > 0) {
	    print_end_redirection(s->prn);
	    pputs(s->prn, _("An error occurred when 'outfile' was active\n"));
	}
    }

    return (err == E_STOP)? err : ret;
}
