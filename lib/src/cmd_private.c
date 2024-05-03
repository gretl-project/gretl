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

/* cmd_private.c for libgretl: functions used in a small set of
   libgretl C files concerned with command execution.
*/

#include "libgretl.h"
#include "gretl_func.h"
#include "monte_carlo.h"
#include "cmd_private.h"

#define PMDEBUG 0
#define COMP_DEBUG 0

void gretl_exec_state_init (ExecState *s,
                            ExecFlags flags,
                            char *line,
                            CMD *cmd,
                            MODEL *model,
                            PRN *prn)
{
    s->flags = flags;

    s->line = line;
    if (s->line != NULL) {
        s->line[0] = '\0';
    }

    s->cmd = cmd;
    if (s->cmd != NULL) {
        s->cmd->ci = 0;
    }

    *s->runfile = '\0';

    s->model = model;
    s->prn = prn;

    s->pmod = NULL;
    s->sys = NULL;
    s->rset = NULL;
    s->var = NULL;
    s->in_comment = 0;
    s->padded = 0;

    if (flags == FUNCTION_EXEC) {
        /* On entry to function execution we check if there's
           a 'last model' in place. If so, we want to make
           this invisible within the function, but set things
           up so that we can restore it as last model on
           exit from the function -- the idea being that
           executing a function should not change the 'last
           model' state at caller level. To achieve this we
           need to take out a 'private' reference to the
           model, stored in the ExecState, and then remove
           it from last model position for the present.
        */
        s->prev_model = get_last_model(&s->prev_type);
        if (s->prev_model != NULL) {
#if PMDEBUG
            fprintf(stderr, "ExecState %p: set prev_model %p\n",
                    (void *) s, s->prev_model);
#endif
            gretl_object_ref(s->prev_model, s->prev_type);
            set_as_last_model(NULL, GRETL_OBJ_NULL);
        }
        s->prev_model_count = get_model_count();
    } else {
        s->prev_model = NULL;
        s->prev_type = GRETL_OBJ_NULL;
        s->prev_model_count = -1;
    }

    s->loop = NULL;
    s->submask = NULL;
    s->callback = NULL;
}

#define GENCOMP_DEBUG 0

static int try_compile_func_genr (ExecState *s,
                                  DATASET *dset,
                                  void *ptr,
                                  int *done)
{
    GENERATOR **pgen = ptr;
    GretlType gtype = s->cmd->gtype;
    const char *line = s->cmd->vstart;
    gretlopt gopt = OPT_NONE;
    int err = 0;

    if (s->cmd->opt & OPT_O) {
	gopt |= OPT_O;
    }

#if !COMPILE_FEVAL
    /* see comment on COMPILE_FEVAL in cmd_private.h */
    if (strstr(line, "feval")) {
	return 0;
    }
#endif

#if GENCOMP_DEBUG
    fprintf(stderr, "cmd_private: calling genr_compile on '%s'\n", line);
#endif
    *pgen = genr_compile(line, dset, gtype, gopt, s->prn, &err);
    if (!err && *pgen != NULL) {
#if GENCOMP_DEBUG
	const char *funname = NULL;

	current_function_info(&funname, NULL);
        fprintf(stderr, "cmd_private: attached genr %p to function %s, depth %d\n",
                (void *) *pgen, funname, gretl_function_depth());
        fprintf(stderr, "  '%s'\n", line);
#endif
        *done = 1;
    } else if (err == E_EQN) {
	/* may be a non-compilable special such as "genr time",
           or perhaps a bare declaration */
        gretl_error_clear();
	err = 0;
    }

    return err;
}

/* Called by functions, and by scripts executed from within
   functions. Augmented 2022-08-11 to support a 3rd argument,
   for use in a function that is being called from a loop or
   internal iteration: we'll attempt to "compile" GENR
   statements in this context, for reuse when the function is
   next called.
*/

int maybe_exec_line (ExecState *s, DATASET *dset, void *ptr)
{
    int done = 0;
    int err = 0;

    if (string_is_blank(s->line)) {
        return 0;
    }

    if (gretl_compiling_loop()) {
        err = get_command_index(s, LOOP, 0);
    } else {
        err = parse_command_line(s, dset, ptr);
        if (s->cmd->ci == GENR && !err && ptr != NULL) {
            if (!(s->cmd->flags & (CMD_SUBST | CMD_CATCH))) {
                err = try_compile_func_genr(s, dset, ptr, &done);
            }
        }
    }

    if (err) {
	errmsg(err, s->prn);
	if (s->cmd->flags & CMD_CATCH) {
	    set_gretl_errno(err);
	    return 0;
	} else {
	    return err;
	}
    }

    gretl_exec_state_transcribe_flags(s, s->cmd);

    if (s->cmd->ci < 0) {
        return 0; /* nothing there, or a comment */
    }

    if (s->cmd->ci == LOOP || gretl_compiling_loop()) {
        /* accumulating loop commands */
        err = gretl_loop_append_line_full(s, dset, ptr);
        if (err) {
            errmsg(err, s->prn);
        }
        return err;
    }

    s->pmod = NULL; /* be on the safe side */

    if (s->cmd->ci == FUNCERR) {
        err = E_FUNCERR;
    } else if (!done) {
        /* note: error messages may be printed to s->prn */
        err = gretl_cmd_exec(s, dset);
    }

    return err;
}

/* Try to determine if a function calls itself, in case we want to
   block "compilation" of such functions. We have to be careful to
   avoid false positives.
*/

static int check_for_recursion (const char *line, const char *targ)
{
    const char *p = strstr(line, targ);
    int ret = 0;

    if (p != NULL) {
	if (p - line > 0) {
	    char c = *(p-1);

	    ret = !isalnum(c) && c != '_';
	} else {
	    ret = 1;
	}
    }

    return ret;
}

/* Takes an array of "stmt" structs (see cmd_private.h) from a
   function or loop and examines them for structure, in terms of
   conditionality (if/elif/else/endif). The object of the exercise is
   to figure out which line we should skip to in case of a false
   if-condition, and on reaching the end of a block enabled by a true
   if-condition. We thereby construct a set of efficiency-enhancing
   "gotos", which are recorded in the "next" member of each relevant
   line.

   In the case of functions only we also determine and record embedded
   loop structure -- and we ignore conditionality within such loops,
   which will be handled by the specific loop-compilation code. (Note
   the asymmetry: a function may contain one or more loops but a loop
   cannot contain a function definition.)

   The @context argument should be FUNC or LOOP to distinguish the two
   cases. The @name argument is used only in debugging, to identify
   the function we're working on; in the LOOP case we just pass
   "loop".
*/

int statements_get_structure (stmt *lines, int n_lines,
			      const char *name,
			      int *recurses)
{
    gchar *self_call = NULL;
    stmt *line;
    int *match_start = NULL;
    int *match_end = NULL;
    int *next_from_loop = NULL;
    int if_depth = 0;
    int max_if_depth = 0;
    int loop_depth = 0;
    int max_loop_depth = 0;
    int loop_imax = 0;
    int loop_imin = -1;
    int target, src;
    int context;
    int i, j, ci;
    int err = 0;

    if (!strcmp(name, "loop")) {
	context = LOOP;
    } else {
	context = FUNC;
	if (recurses != NULL) {
	    self_call = g_strdup_printf("%s(", name);
	}
    }

    /* first pass: determine the maximum depth of conditionality (and
       looping, if in a function). Also record the indices of the
       first and last lines "of interest" (loop_imin and loop_imax).
    */
    for (i=0; i<n_lines; i++) {
	ci = lines[i].ci;
#if COMP_DEBUG > 1
        if (context == LOOP) {
            fprintf(stderr, "LOOP L%d: '%s'\n", i, lines[i].s);
        }
#endif
        if (context == FUNC) {
	    if (recurses != NULL && *recurses == 0 && ci == GENR) {
		*recurses = check_for_recursion(lines[i].s, self_call);
	    }
            if (ci == LOOP) {
                if (loop_imin < 0) {
                    loop_imin = i;
                }
                if (++loop_depth > max_loop_depth) {
                    max_loop_depth = loop_depth;
                }
            } else if (ci == ENDLOOP) {
                if (i > loop_imax) {
                    loop_imax = i;
                }
                loop_depth--;
            } else if (loop_depth > 0) {
                continue;
            }
        }
	if (ci == IF) {
            if (loop_imin < 0) {
                loop_imin = i;
            }
            if (++if_depth > max_if_depth) {
                max_if_depth = if_depth;
            }
	} else if (ci == ENDIF) {
            if (i > loop_imax) {
                loop_imax = i;
            }
            if_depth--;
	}
    }

#if COMP_DEBUG
    fprintf(stderr, "\n%s: max if-depth %d, max loop-depth %d\n",
            name, max_if_depth, max_loop_depth);
#endif

    if (if_depth != 0 || loop_depth != 0) {
	gretl_errmsg_sprintf(_("Broken syntax in %s: unmatched %s"), name,
			     if_depth ? "if/endif" : "loop/endloop");
	return E_PARSE;
    }

    if (max_if_depth > 0) {
        match_start = gretl_list_new(max_if_depth);
        match_end = gretl_list_new(max_if_depth);
    }
    if (max_loop_depth > 0) {
        next_from_loop = gretl_list_new(max_loop_depth);
    }

    /* second pass: analysis */
    for (i=loop_imin; i<=loop_imax && !err; i++) {
	line = &lines[i];
	if (context == FUNC) {
            if (line->ci == LOOP) {
                loop_depth++;
                next_from_loop[loop_depth] = i;
            } else if (line->ci == ENDLOOP) {
                if (loop_depth == 0) {
                    err = E_PARSE;
                } else {
                    j = next_from_loop[loop_depth];
                    lines[j].next = i;
                    loop_depth--;
                }
            } else if (loop_depth > 0) {
                continue;
            }
        }
	if (line->ci == IF) {
	    match_start[++if_depth] = i;
	} else if (line->ci == ENDIF) {
	    if (if_depth == 0) {
		err = E_PARSE;
	    } else {
                line->next = -if_depth;
		j = match_start[if_depth];
		lines[j].next = i;
                if_depth--;
	    }
	} else if (line->ci == ELIF || line->ci == ELSE) {
	    if (if_depth == 0) {
		err = E_PARSE;
	    } else {
                if (lines[i-1].ci != ENDIF) {
                    lines[i-1].next = -if_depth;
                }
		j = match_start[if_depth];
		lines[j].next = i;
		match_start[if_depth] = i;
	    }
	}
    }

    if (!err) {
        /* third pass: fill in goto's for true-block terminators */
        loop_depth = if_depth = target = src = 0;
        for (i=loop_imax; i>=loop_imin; i--) {
            line = &lines[i];
            if (context == FUNC) {
                if (line->ci == ENDLOOP) {
                    loop_depth++;
                } else if (line->ci == LOOP) {
                    loop_depth--;
                } else if (loop_depth > 0) {
                    continue;
                }
            }
            if (line->ci == ENDIF) {
                if_depth++;
                target = match_start[if_depth] = line->next;
                src = match_end[if_depth] = i;
            } else if (line->ci == IF) {
                target = (if_depth == 0)? 0 : match_start[if_depth];
                src = (if_depth == 0)? 0 : match_end[if_depth];
                if_depth--;
            } else if (target < 0 && line->next == target) {
                line->next = src;
            }
        }
    }

    free(match_start);
    free(match_end);
    free(next_from_loop);

#if COMP_DEBUG
    /* display what we figured out */
    fputc('\n', stderr);
    for (i=0; i<n_lines; i++) {
	line = &lines[i];
        j = line->next;
	if (j <= 0) {
	    continue;
	}
	if (line->ci == IF || line->ci == ELIF || line->ci == ELSE) {
	    fprintf(stderr, "L%d ('%s'): next-on-false = %d ('%s')\n",
		    i, line->s, j, lines[j].s);
        } else if (context == FUNC && line->ci == LOOP) {
	    fprintf(stderr, "L%d ('%s'): end-of-loop = %d\n", i, line->s, j);
	} else {
	    fprintf(stderr, "L%d ('%s'): next-on-true = %d ('%s')\n",
		    i, line->s, j, lines[j].s);
	}
    }
#endif

    return err;
}
